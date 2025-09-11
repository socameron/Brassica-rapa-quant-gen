#### PROJECT: Brassica rapa Va/W Study (Data collected by Cameron So 2019 at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Construct within heated environment models using MCMCglmm to extract Va, Vd, and Vm
#### AUTHOR: Cameron So

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Necessary Packages
library(tidyverse)
library(MCMCglmm)
library(MuMIn)
library(logitnorm)
library(nadiv)
library(QGglmm)

#Importing Data using readr from tidyverse (Base R confuses the class for certain vectors). 
# d for double, f for factor, i for integer
col_types_list2 <- cols_only(posID = "d", individual = "d", plot = "f",
                             animal = "f",
                             patID = "f", matID = "f", famID = "f",
                             treatment = col_factor(levels=c("A", "H")),
                             germ = "d", flower = "d", 
                             seed = "d", leaf = "d", flwr_clstr = "d", seed_pods = "d", seed_number = "d", 
                             height = "d", stem_diam = "d",germ_census = "d", flwr_census = "d")

MCMC.2019 <- read_csv("Rdata/MCMC_2019_cleaned.csv", col_names = TRUE, na = "NA",
                      col_types = col_types_list2)
spec(MCMC.2019)

ped <- read.csv("Rdata/heatarrays_animal_pedigree.csv") #note don't use readr to import the csv file.
for (x in 1:3) ped[, x] <- as.factor(ped[, x])

#Changing tibbles into dataframes because MCMCglmm cannot read tibbles
MCMC.2019 <- as.data.frame(MCMC.2019)
ped <- as.data.frame(ped)

MCMC.heated <- MCMC.2019 %>%
  filter(treatment=="H")

MCMC.heated$animal <- as.factor(MCMC.heated$animal) #double check
lapply(ped, class)
lapply(MCMC.heated, class)

#Function to stack MCMC chains, borrowed from Paul Maurizio 
#Note, requires a mcmc coda object, not a list. So any parallel runs will need to be converted. 
#See : https://github.com/mauriziopaul/treatmentResponseDiallel
mcmc.stack <- function (coda.object, ...){
  ## This function is from Will; also part of BayesDiallel
  if (inherits(coda.object, "mcmc")) {
    return(coda.object)
  }
  if (!inherits(coda.object, "mcmc.list")) {
    stop("Non-mcmc object passed to function\n")
  }
  chain <- coda.object[[1]]
  for (i in 2:nchain(coda.object)) {
    chain <- rbind(chain, coda.object[[i]])
  }
  as.mcmc(chain)
}

#'####################################################################'#
##############      MCMCglmm ANALYSIS FOR HEATED       ###############
#'####################################################################'#

#Some terminology: 
#R argument list: prior for the residual variance
#G argument list: random effects variance
#For each random effect, there must be a prior (e.g G1, G2, G3)
#V = variance
#nu = degree of belief
#fix = number of fixed effects
# if using more than one fixed response variable (e.g seed_pods ~ sex ~ colour), need to use rcov which binds the two variables.
# Thus, you get rcov = ~us (trait):units ... similar to the cbind function
# Note, ~us can be interchanged with another function idh.
# us = genetic correlation b/w the processes to be estimated
# idh = genetic correlation is 0

#For binary data, need to add alpha mu and V 
#alpha.mu = prior mean
#alpha.V = prior covaraince matrix

#Creating Dominance Matrix for Vd estimates
Ainv <- inverseA(ped[, 1:3])$Ainv
Dinv <- makeD(ped[, 1:3])$Dinv
MCMC.heated$animalDom <- MCMC.heated$animal















#### No. 1 Lifetime Fitness - "Univariate" | Zero-inflated Poisson ####
#We set two different priors for the zero inflation process, and the Poisson process. However, according to Hadfield (https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/018802.html), the genetic correlation between the two processes must be set to 0 and cannot be estimated. 

MCMC.heated %>%
  summarise(mean=mean(seed_pods), sd=sd(seed_pods), var=(sd)^2) #var = 519.2583

#Model 1.0 with additive + dominance + maternal random effects and plot fixed effects
#Zero-Inflated Poisson distribution + 2 different priors
#Level 1 = Poisson process, Level 2 = Zero inflation
heated.total <- MCMC.heated
prior1.0 <- list(R = list(V = diag(2), nu = 1, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G5 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G6 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) # Chi Square Prior with SD 1

HW_model1.0 <- MCMCglmm(seed_pods ~ trait - 1 + plot,
                     random = ~idh(at.level(trait,1)):animal + idh(at.level(trait,2)):animal +
                       idh(at.level(trait,1)):animalDom + idh(at.level(trait,2)):animalDom + 
                       idh(at.level(trait,1)):matID + idh(at.level(trait,2)):matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv), 
                     rcov = ~idh(trait):units,
                     family = "zipoisson", data = heated.total, prior=prior1.0, 
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model1.0.RData")
plot(HW_model1.0$VCV) # good trace plots
summary(HW_model1.0) # good effective sample size
autocorr.diag(HW_model1.0$VCV) # autocorrelation for dominance at Poisson process
heidel.diag(HW_model1.0$VCV) # convergence success

#Heritability - Zero inflation process
HW_herit1.0_ZI <- HW_model1.0$VCV[, "at.level(trait, 2).animal"]/(HW_model1.0$VCV[, "at.level(trait, 2).animal"] + HW_model1.0$VCV[, "at.level(trait, 2).animalDom"] + HW_model1.0$VCV[, "at.level(trait, 2).matID"] + HW_model1.0$VCV[, "traitzi_seed_pods.units"])
mean(HW_herit1.0_ZI)
posterior.mode(HW_herit1.0_ZI)
HPDinterval(HW_herit1.0_ZI)

#Heritability - Poisson Process
HW_herit1.0_PP <- HW_model1.0$VCV[, "at.level(trait, 1).animal"]/(HW_model1.0$VCV[, "at.level(trait, 1).animal"] + HW_model1.0$VCV[, "at.level(trait, 1).animalDom"] + HW_model1.0$VCV[, "at.level(trait, 1).matID"] + HW_model1.0$VCV[, "traitseed_pods.units"])
mean(HW_herit1.0_PP)
posterior.mode(HW_herit1.0_PP)
HPDinterval(HW_herit1.0_PP)

#Transforming from latent to data scale using package "QGglmm"
#Zero-inflation process
yhat_H1.0_ZI <- predict(HW_model1.0, type = "terms")
vp_H1.0_ZI <- mean(HW_model1.0$VCV[, "at.level(trait, 2).animal"] + HW_model1.0$VCV[, "at.level(trait, 2).animalDom"] + HW_model1.0$VCV[, "at.level(trait, 2).matID"] + HW_model1.0$VCV[, "traitzi_seed_pods.units"])
va_H1.0_ZI <- mean(HW_model1.0[["VCV"]][ , "at.level(trait, 2).animal"]) #ZERO INFLATION
vd_H1.0_ZI <- mean(HW_model1.0[["VCV"]][ , "at.level(trait, 2).animalDom"]) #ZERO INFLATION
vm_H1.0_ZI <- mean(HW_model1.0[["VCV"]][ , "at.level(trait, 2).matID"]) #ZERO INFLATION

#Additive Variance
QGA_H1.0_ZI <- QGparams(predict=yhat_H1.0_ZI, var.a=va_H1.0_ZI, var.p=vp_H1.0_ZI, model = "binom1.probit")
QGA_H1.0_ZI #Va = 0.003005507 .. 1.64% of Vp is Va .. Obs Mean = 0.759

#Dominance Variance
QGD_H1.0_ZI <- QGicc(predict=yhat_H1.0_ZI, var.comp=vd_H1.0_ZI, var.p=vp_H1.0_ZI, model = "binom1.probit")
QGD_H1.0_ZI #Vd = 0.008303287 .. 4.54% of Vp is Vd

#Maternal Variance
QGM_H1.0_ZI <- QGicc(predict=yhat_H1.0_ZI, var.comp=vm_H1.0_ZI, var.p=vp_H1.0_ZI, model = "binom1.probit")
QGM_H1.0_ZI #Vm = 0.00410977 .. 2.24% of Vp is Vm



#Poisson Process
yhat_H1.0_PP <- predict(HW_model1.0, type = "terms")
vp_H1.0_PP <- mean(HW_model1.0$VCV[, "at.level(trait, 1).animal"] + HW_model1.0$VCV[, "at.level(trait, 1).animalDom"] + HW_model1.0$VCV[, "at.level(trait, 1).matID"] + HW_model1.0$VCV[, "traitseed_pods.units"])
va_H1.0_PP <- mean(HW_model1.0[["VCV"]][ , "at.level(trait, 1).animal"]) #POISSON PROCESS
vd_H1.0_PP <- mean(HW_model1.0[["VCV"]][ , "at.level(trait, 1).animalDom"]) #POISSON PROCESS
vm_H1.0_PP <- mean(HW_model1.0[["VCV"]][ , "at.level(trait, 1).matID"]) #POISSON PROCESS

#Additive Variance
QGA_H1.0_PP <- QGparams(predict=yhat_H1.0_PP, var.a=va_H1.0_PP, var.p=vp_H1.0_PP, model = "Poisson.log")
QGA_H1.0_PP #Va = 20.97739 .. 1.66% of Vp is Va .. Obs Mean = 14.129

#Dominance Variance
QGD_H1.0_PP <- QGicc(predict=yhat_H1.0_PP, var.comp=vd_H1.0_PP, var.p=vp_H1.0_PP, model = "Poisson.log")
QGD_H1.0_PP #Vd = 46.79955 .. 3.71% of Vp is Vd

#Maternal Variance
QGM_H1.0_PP <- QGicc(predict=yhat_H1.0_PP, var.comp=vm_H1.0_PP, var.p=vp_H1.0_PP, model = "Poisson.log")
QGM_H1.0_PP #Vm = 13.86292 .. 1.10% of Vp is Vm


#Model 1.1
heated.total <- MCMC.heated
prior1.1 <- list(R = list(V = diag(2), nu = 1, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1), # Expanded Fisher prior for Poisson process
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G5 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1), # Expanded Fisher prior for Poisson process
                          G6 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) # Chi Square Prior with SD 1

HW_model1.1 <- MCMCglmm(seed_pods ~ trait - 1 + plot,
                        random = ~idh(at.level(trait,1)):animal + idh(at.level(trait,2)):animal +
                          idh(at.level(trait,1)):animalDom + idh(at.level(trait,2)):animalDom + 
                          idh(at.level(trait,1)):matID + idh(at.level(trait,2)):matID, 
                        ginverse = list(animal = Ainv, animalDom = Dinv), 
                        rcov = ~idh(trait):units,
                        family = "zipoisson", data = heated.total, prior=prior1.1, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model1.1.RData")
plot(HW_model1.1$VCV) # 
summary(HW_model1.1) # 
autocorr.diag(HW_model1.1$VCV) # 
heidel.diag(HW_model1.1$VCV) #

#Model 1.2 - No dominance, Poisson family
load(file="Routput/Heated/HW_model1.2.RData")
plot(HW_model1.2$VCV) # 
summary(HW_model1.2) #
autocorr.diag(HW_model1.2$VCV) #
heidel.diag(HW_model1.2$VCV) #

#Model 1.3 - No dominance, Zero inflated Poisson family
load(file="Routput/Heated/HW_model1.3.RData")
plot(HW_model1.3$VCV) # 
summary(HW_model1.3) #
autocorr.diag(HW_model1.3$VCV) #
heidel.diag(HW_model1.3$VCV) #












#### No. 2 Lifetime Fitness V2: Survival x Fecundity | Bivariate ####
#Note: For this trait, we include all plants from the experiment
#Survival to Flowering = Bernoulli Dist | Fecundity = Poisson Dist

#Calculating Variance for Fecundity for the Non-informative Prior (equal variance)
MCMC.heated %>%
  summarise(mean=mean(seed_pods), sd=sd(seed_pods), var=(sd)^2) #mean = 9.261649
log(9.261649) # = 2.225882

MCMC.heated %>%
  summarise(sum.f=sum(flower), n=n(), mean=sum.f/n) #mean prob = 0.4423536
logit(0.4423536*(1-0.4423536))

#Model 2.0 with additive, dominance, and maternal random effects with plot fixed effects
#Using prior similar Bemmels & Anderson 2019
#Note: we set rcov = ~units because covariance can't be estimated in residual variance for a binary trait, AND we set residual variance to fix = 2 (the number of multivariate dimensions)
heated.Survival <- MCMC.heated
prior2.0 <- list(R = list(V = diag(2)*0.02, nu = 3, fix = 2),
                 G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000))),
                          G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000))),
                          G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000)))
                 ))
prior2.0$R$V[2,2] <- 10 # fixing last element of residual variance to 10 to improve mixing, following Bemmels & Anderson

HW_model2.0 <- MCMCglmm(cbind(seed_pods, flower) ~ trait + plot - 1,
                        random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv), rcov = ~us(trait):units,
                        family = c("poisson", "threshold"), 
                        data = heated.Survival, prior = prior2.0, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Heated/HW_model2.0.RData")
plot(HW_model2.0$VCV) # poor trace plots
summary(HW_model2.0) # low effective sample size
autocorr.diag(HW_model2.0$VCV) # high autocorrelation
heidel.diag(HW_model2.0$VCV) # convergence fail

#Model 2.1 without dominance covariance
heated.Survival <- MCMC.heated
prior2.1 <- list(R = list(V = diag(2)*0.02, nu = 3, fix = 2),
                 G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000))),
                          G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000))),
                          G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000)))
                 ))
prior2.1$R$V[2,2] <- 10 # fixing last element of residual variance to 10 to improve mixing, following Bemmels & Anderson

HW_model2.1 <- MCMCglmm(cbind(seed_pods, flower) ~ trait + plot - 1,
                        random = ~us(trait):animal + idh(trait):animalDom + us(trait):matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv), rcov = ~us(trait):units,
                        family = c("poisson", "threshold"), 
                        data = heated.Survival, prior = prior2.1, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Heated/HW_model2.1.RData")
plot(HW_model2.1$VCV) # poor trace plots
summary(HW_model2.1) # low sample sizes
autocorr.diag(HW_model2.1$VCV) # bad autocorrelation
heidel.diag(HW_model2.1$VCV) # convergence fail

#Calculating Genetic Correlations = Covariance b/w trait 1 & 2 / sqrt (var [trait 1] * var [trait 2])
# (1) Survival & Fecundity
gen.corrH2.1 <-HW_model2.1$VCV[,'traitflower:traitseed_pods.animal']/
  sqrt(HW_model2.1$VCV[,'traitflower:traitflower.animal']*HW_model2.1$VCV[,'traitseed_pods:traitseed_pods.animal']) 
mean(gen.corrH2.1) #
posterior.mode(gen.corrH2.1) #
HPDinterval(gen.corrH2.1) #

phen.corrH2.1 <- cor.test(heated.Survival$flower, heated.Survival$seed_pods, test="pearson")
phen.corrH2.1 # Corr = 0.4564391 | Posterior 95% CI =  (0.4291918, 0.4828564)

#Model 2.2 - No dominance
load(file="Routput/Heated/HW_model2.2.RData")
plot(HW_model2.2$VCV) # 
summary(HW_model2.2) #
autocorr.diag(HW_model2.2$VCV) #
heidel.diag(HW_model2.2$VCV) #

#Model 2.3 - No dominance and maternal effects, nu = 10
load(file="Routput/Heated/HW_model2.3.RData")
plot(HW_model2.3$VCV) # 
summary(HW_model2.3) #
autocorr.diag(HW_model2.3$VCV) #
heidel.diag(HW_model2.3$VCV) #

#Model 2.4 - Full model agian but attempting new prior
heated.Survival <- MCMC.heated
prior2.4 <- list(R = list(V = diag(2), nu = 2, fix = 2),
                 G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G3 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))
                 ))

set.seed(1)
HW_model2.4 <- mclapply(1:20, function(i) {
	MCMCglmm(cbind(seed_pods, flower) ~ trait + plot - 1,
                        random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv), rcov = ~us(trait):units,
                        family = c("poisson", "threshold"), 
                        data = heated.Survival, prior = prior2.4, 
                        nitt = 200000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)
}, mc.cores = 20)

load(file="Routput/Heated/HW_model2.4.RData")
HW_2.4_sol <- lapply(HW_model2.4, function(m) m$Sol)
HW_2.4_sol <- do.call(mcmc.list, HW_2.4_sol)

HW_2.4_vcv <- lapply(HW_model2.4, function(m) m$VCV)
HW_2.4_vcv <- do.call(mcmc.list, HW_2.4_vcv)
HW_2.4_vcv <- mcmc.stack(HW_2.4_vcv)
plot(HW_2.4_vcv) #poor trace plots
autocorr.diag(HW_2.4_vcv, relative=FALSE) #high autocorrelation
heidel.diag(HW_2.4_vcv) #convergence fail
effectiveSize(HW_2.4_vcv) #low effective size

#Since the multiple models failed to converge, we estimate the genetic correlation by extracting breeding values from univariate models. 
#See 05_MCMC_full_analysis for the estimation of the genetic correlation between fecundity and survival.

load(file="Routput/heatarrays_data_explore.RData")
bivariate <- dat %>%
  group_by(treatment, famID) %>%
  summarise(sum_flower=sum(flower), n=n(), mean_pods=mean(seed_pods)) %>%
  mutate(mean_flower=(sum_flower/n))

bivariate %>%
  ggplot(aes(x=mean_flower, y=mean_pods, color=treatment)) + 
  geom_point(size=2, alpha=0.5) +
  labs(x="Survival Success", y="Fecundity") +
  ggtitle("Bivariate Model - Family Means") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated"))







#### No. 3 Survival to Flowering - Univariate | Bernoulli ####
#Note: For this trait, we include all plants from the experiment

#Probability of Survival to Flowering
heated.Survival <- MCMC.heated
heated.Survival %>% 
  summarise(tot.flwr=sum(flower), n=n(), prob=tot.flwr/n) #prob = 0.4426181
log(0.4426*(1-0.4426)) # -1.399

#Model 3.0 with additive, dominance and maternal random effects + plot fixed effects
#We use a prior suggested by Villemereuil et al. 2013 for binary traits
heated.Survival <- MCMC.heated
prior3.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model3.0 <- MCMCglmm(flower ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = heated.Survival, prior = prior3.0, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)

load(file="Routput/Heated/HW_model3.0.RData")
plot(HW_model3.0$VCV) # good trace plots
summary(HW_model3.0) # good effective sample size
autocorr.diag(HW_model3.0$VCV) # no autocorrelation
heidel.diag(HW_model3.0$VCV) # convergence success for all

#Latent Scale 
mean(HW_model3.0[["VCV"]][ , "animal"]) #Va
mean(HW_model3.0[["VCV"]][ , "animalDom"]) #Vd
mean(HW_model3.0[["VCV"]][ , "matID"]) #Vm
HW_herit3.0 <- HW_model3.0$VCV[, "animal"]/(HW_model3.0$VCV[, "animal"] + HW_model3.0$VCV[, "animalDom"] + HW_model3.0$VCV[, "matID"] + HW_model3.0$VCV[, "units"])
mean(HW_herit3.0)
posterior.mode(HW_herit3.0) #Heritability
HPDinterval(HW_herit3.0) #h2 HPD

mean(HW_model3.0$VCV[, "animal"] + HW_model3.0$VCV[, "animalDom"] + HW_model3.0$VCV[, "matID"]) #Vg
HPDinterval(HW_model3.0$VCV[, "animal"] + HW_model3.0$VCV[, "animalDom"] + HW_model3.0$VCV[, "matID"])

#Transforming from latent to data scale using package "QGglmm"
df_H3.0 <- data.frame(
  va = as.vector(HW_model3.0[["VCV"]][, "animal"]),
  vd = as.vector(HW_model3.0[["VCV"]][, "animalDom"]),
  vm = as.vector(HW_model3.0[["VCV"]][, "matID"]),
  vp = rowSums(HW_model3.0[["VCV"]]))
yhat_H3.0 <- predict(HW_model3.0, type = "terms")

#Va
Va_H3.0 <- do.call("rbind", apply(df_H3.0, 1, function(row){
  QGparams(predict = yhat_H3.0, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
  }))
mean(as.mcmc(Va_H3.0$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_H3.0$var.a.obs))

mean(as.mcmc(Va_H3.0$h2.obs)) #h2
HPDinterval(as.mcmc(Va_H3.0$h2.obs))

mu_H3.0 <- mean(as.mcmc(Va_H3.0$mean.obs))
mu_H3.0
HPDinterval(as.mcmc(Va_H3.0$mean.obs))

# Check for significant difference between environments. Requires ambient data
mu_diff_3.0 <- Va_H3.0$mean.obs - Va_A3.0$mean.obs
mean(mu_diff_3.0)
quantile(mu_diff_3.0, probs = c(0.025, 0.975))
mean(mu_diff_3.0 > 0)

mean(as.mcmc(Va_H3.0$h2.obs))/mu_H3.0 #Rw
HPDinterval(as.mcmc(Va_H3.0$h2.obs))/mu_H3.0

#Vd
Vd_H3.0 <- do.call("rbind", apply(df_H3.0, 1, function(row){
  QGicc(predict = yhat_H3.0, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vd_H3.0$var.comp.obs))
HPDinterval(as.mcmc(Vd_H3.0$var.comp.obs))

#Vm
Vm_H3.0 <- do.call("rbind", apply(df_H3.0, 1, function(row){
  QGicc(predict = yhat_H3.0, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_H3.0$var.comp.obs))
HPDinterval(as.mcmc(Vm_H3.0$var.comp.obs))

#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_H3.0 <- data.frame(
  Va = Va_H3.0$var.a.obs,
  Vd = Vd_H3.0$var.comp.obs,
  Vm = Vm_H3.0$var.comp.obs,
  mu = Va_H3.0$mean.obs,
  h2 = Va_H3.0$h2.obs,
  Rw = (Va_H3.0$var.a.obs / Va_H3.0$mean.obs)) %>%
  pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_H3.0 <- data.frame(
  Component = c("Va", "Vd", "Vm", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(Va_H3.0$var.a.obs)),
    mean(as.mcmc(Vd_H3.0$var.comp.obs)),
    mean(as.mcmc(Vm_H3.0$var.comp.obs)),
    mean(as.mcmc(Va_H3.0$mean.obs)),
    mean(as.mcmc(Va_H3.0$h2.obs)),
    mean(as.mcmc(Va_H3.0$var.a.obs / Va_H3.0$mean.obs))
  ),
  Lower = c(
    HPDinterval(as.mcmc(Va_H3.0$var.a.obs))[1],
    HPDinterval(as.mcmc(Vd_H3.0$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vm_H3.0$var.comp.obs))[1],
    HPDinterval(as.mcmc(Va_H3.0$mean.obs))[1],
    HPDinterval(as.mcmc(Va_H3.0$h2.obs))[1],
    HPDinterval(as.mcmc(Va_H3.0$var.a.obs / Va_H3.0$mean.obs))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(Va_H3.0$var.a.obs))[2],
    HPDinterval(as.mcmc(Vd_H3.0$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vm_H3.0$var.comp.obs))[2],
    HPDinterval(as.mcmc(Va_H3.0$mean.obs))[2],
    HPDinterval(as.mcmc(Va_H3.0$h2.obs))[2],
    HPDinterval(as.mcmc(Va_H3.0$var.a.obs / Va_H3.0$mean.obs))[2]
  )
)

posteriors_H3.0$Component <- factor(posteriors_H3.0$Component, levels = c("Vm", "Va", "Vd", "mu", "Rw", "h2"))
save(posteriors_H3.0, stats_H3.0, file="Routput/Heated/HW_stats3.0.RData")
load("Routput/Heated/HW_stats3.0.RData")

var_H3.0_fig <- posteriors_H3.0 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, 0.03) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown()) +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H3.0_fig

ggsave(plot = var_H3.0_fig, filename = "Routput/Figures/fig_4_H3.0_var.tif", width = 8, height = 8, units="cm")

var_H3.0_fig2 <- posteriors_H3.0 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  scale_color_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  labs(x = expression("Proportion of variance (" * h^2 * " or " * Δ[evol] * bar(W) * ")"),
       y = "Density") +
  xlim(0, 0.04) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H3.0_fig2

ggsave(plot = var_H3.0_fig2, filename = "Routput/Figures/fig_4_H3.0_h2.tif", width = 8, height = 8, units="cm")


#Plotting density plot of breeding values
HW_3.0_breed <- HW_model3.0$Sol[,7:7657]
HW_3.0_BV <- posterior.mode(as.mcmc(HW_breed)) #breeding values for each individual
HW_3.0_BV_CI <- HPDinterval(as.mcmc(HW_breed)) #confidence intervals of breeding values
HW_3.0_BV_plot <- as.data.frame(HW_3.0_BV)
HW_3.0_BV_plot %>% 
  ggplot(aes(x=HW_BV)) +
  geom_density(fill="tomato2", color="tomato2", alpha=0.3) +
  labs(x="Ambient Survival Breeding Values", y="Frequency")  + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

#Model 3.1 - No dominance
load(file="Routput/Heated/HW_model3.1.RData")
plot(HW_model3.1$VCV) # good trace plots 
summary(HW_model3.1) # good effective sample size
autocorr.diag(HW_model3.1$VCV) # no autocorrelation
heidel.diag(HW_model3.1$VCV) # convergence success

mean(HW_model3.1[["VCV"]][ , "animal"]) #Va
mean(HW_model3.1[["VCV"]][ , "matID"]) #Vm
mean(HW_model3.1[["VCV"]][ , "units"]) #Vr


#Transforming from latent to data scale using package "QGglmm"
df_H3.1 <- data.frame(
  va = as.vector(HW_model3.1[["VCV"]][, "animal"]),
  vm = as.vector(HW_model3.1[["VCV"]][, "matID"]),
  vr = as.vector(HW_model3.1[["VCV"]][, "units"]),
  vp = rowSums(HW_model3.1[["VCV"]]))
yhat_H3.1 <- predict(HW_model3.1, type = "terms")

#Va
Va_H3.1 <- do.call("rbind", apply(df_H3.1, 1, function(row){
  QGparams(predict = yhat_H3.1, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Va_H3.1$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_H3.1$var.a.obs))

#Vm
Vm_H3.1 <- do.call("rbind", apply(df_H3.1, 1, function(row){
  QGicc(predict = yhat_H3.1, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_H3.1$var.comp.obs))
HPDinterval(as.mcmc(Vm_H3.1$var.comp.obs))


#Vr
Vr_H3.1 <- do.call("rbind", apply(df_H3.1, 1, function(row){
  QGicc(predict = yhat_H3.1, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vd_H3.1$var.comp.obs)) #shouldn't work because residual variance is fixed at 1
HPDinterval(as.mcmc(Vd_H3.1$var.comp.obs))






#### No. 4 Fecundity of Flowering Plants - Univariate | Poisson ####
#Note: For Fecundity, we only include plants that survived to reach flowering

#Calculating Variance for Fecundity for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1 & flower==1) %>%
  summarise(sd=sd(seed_pods), mean=mean(seed_pods), var=(sd)^2) #var = mean = 20.9372 in Poisson, but Vp = 929
#Variance doesn't matter with variance in QGglmm... maybe because variance is calculated of a different distribution?

#Model 4.0 with additive, dominance, and maternal random effects and plot fixed effects
#Prior with inverse gamma (0.001, 0.001)
heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior4.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model4.0 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "poisson", data = heated.Fecundity, prior = prior4.0,
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model4.0.RData")
plot(HW_model4.0$VCV) # poor trace plots
summary(HW_model4.0) # low effective sample size
autocorr.diag(HW_model4.0$VCV) # high autocorrelation
heidel.diag(HW_model4.0$VCV) # convergence fail

#Model 4.1 with Fisher Parameter Expanded Prior
heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior4.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model4.1 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "poisson", data = heated.Fecundity, prior = prior4.1,
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model4.1.RData")
plot(HW_model4.1$VCV) # good trace plots
summary(HW_model4.1) # good effective sample size
autocorr.diag(HW_model4.1$VCV) # no autocorrelation
heidel.diag(HW_model4.1$VCV) # convergence success

#Latent Scale
mean(HW_model4.1[["VCV"]][ , "animal"]) #Va
mean(HW_model4.1[["VCV"]][ , "animalDom"]) #Vd
mean(HW_model4.1[["VCV"]][ , "matID"]) #Vm
mean(HW_model4.1[["VCV"]][ , "units"]) #Vr
HW_herit4.1 <- HW_model4.1$VCV[, "animal"]/(HW_model4.1$VCV[, "animal"] + HW_model4.1$VCV[, "animalDom"] + HW_model4.1$VCV[, "matID"] + HW_model4.1$VCV[, "units"]) #
mean(HW_herit4.1)
posterior.mode(HW_herit4.1) #Heritability
HPDinterval(HW_herit4.1) #H2 HPD

mean(HW_model4.1$VCV[, "animal"] + HW_model4.1$VCV[, "animalDom"] + HW_model4.1$VCV[, "matID"]) #Vg
HPDinterval(HW_model4.1$VCV[, "animal"] + HW_model4.1$VCV[, "animalDom"] + HW_model4.1$VCV[, "matID"])

#Transforming from latent to data scale using package "QGglmm"
df_H4.1 <- data.frame(
  va = as.vector(HW_model4.1[["VCV"]][, "animal"]),
  vd = as.vector(HW_model4.1[["VCV"]][, "animalDom"]),
  vm = as.vector(HW_model4.1[["VCV"]][, "matID"]),
  vr = as.vector(HW_model4.1[["VCV"]][, "units"]),
  vp = rowSums(HW_model4.1[["VCV"]]))
yhat_H4.1 <- predict(HW_model4.1, type = "terms")

#Va
Va_H4.1 <- do.call("rbind", apply(df_H4.1, 1, function(row){
  QGparams(predict = yhat_H4.1, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Va_H4.1$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_H4.1$var.a.obs))

mean(as.mcmc(Va_H4.1$h2.obs)) #h2
HPDinterval(as.mcmc(Va_H4.1$h2.obs))

mu_H4.1 <- mean(as.mcmc(Va_H4.1$mean.obs))
mu_H4.1
HPDinterval(as.mcmc(Va_H4.1$mean.obs))

# Check for significant difference between environments. Requires ambient data
mu_diff_4.1 <- Va_H4.1$mean.obs - Va_A4.1$mean.obs
mean(mu_diff_4.1)
quantile(mu_diff_4.1, probs = c(0.025, 0.975))
mean(mu_diff_4.1 > 0)

mean(as.mcmc(Va_H4.1$var.a.obs))/mu_H4.1 #Rw
(HPDinterval(as.mcmc(Va_H4.1$var.a.obs)))/mu_H4.1

#Vd
Vd_H4.1 <- do.call("rbind", apply(df_H4.1, 1, function(row){
  QGicc(predict = yhat_H4.1, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vd_H4.1$var.comp.obs))
HPDinterval(as.mcmc(Vd_H4.1$var.comp.obs))

#Vm
Vm_H4.1 <- do.call("rbind", apply(df_H4.1, 1, function(row){
  QGicc(predict = yhat_H4.1, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vm_H4.1$var.comp.obs))
HPDinterval(as.mcmc(Vm_H4.1$var.comp.obs))

#Vr
Vr_H4.1 <- do.call("rbind", apply(df_H4.1, 1, function(row){
  QGicc(predict = yhat_H4.1, 
        var.comp = row[["vr"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vr_H4.1$var.comp.obs))
HPDinterval(as.mcmc(Vr_H4.1$var.comp.obs))

#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_H4.1 <- data.frame(
  Va = Va_H4.1$var.a.obs,
  Vd = Vd_H4.1$var.comp.obs,
  Vm = Vm_H4.1$var.comp.obs,
  Vr = Vr_H4.1$var.comp.obs,
  mu = Va_H4.1$mean.obs,
  h2 = Va_H4.1$h2.obs,
  Rw = (Va_H4.1$var.a.obs / Va_H4.1$mean.obs)) %>%
  pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_H4.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(Va_H4.1$var.a.obs)),
    mean(as.mcmc(Vd_H4.1$var.comp.obs)),
    mean(as.mcmc(Vm_H4.1$var.comp.obs)),
    mean(as.mcmc(Vr_H4.1$var.comp.obs)),
    mean(as.mcmc(Va_H4.1$mean.obs)),
    mean(as.mcmc(Va_H4.1$h2.obs)),
    mean(as.mcmc(Va_H4.1$var.a.obs / Va_H4.1$mean.obs))
  ),
  Lower = c(
    HPDinterval(as.mcmc(Va_H4.1$var.a.obs))[1],
    HPDinterval(as.mcmc(Vd_H4.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vm_H4.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vr_H4.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Va_H4.1$mean.obs))[1],
    HPDinterval(as.mcmc(Va_H4.1$h2.obs))[1],
    HPDinterval(as.mcmc(Va_H4.1$var.a.obs / Va_H4.1$mean.obs))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(Va_H4.1$var.a.obs))[2],
    HPDinterval(as.mcmc(Vd_H4.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vm_H4.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vr_H4.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Va_H4.1$mean.obs))[2],
    HPDinterval(as.mcmc(Va_H4.1$h2.obs))[2],
    HPDinterval(as.mcmc(Va_H4.1$var.a.obs / Va_H4.1$mean.obs))[2]
  )
)

posteriors_H4.1$Component <- factor(posteriors_H4.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "h2", "Rw"))
save(posteriors_H4.1, stats_H4.1, file="Routput/Heated/HW_stats4.1.RData")
load("Routput/Heated/HW_stats4.1.RData")

var_H4.1_fig <- posteriors_H4.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, 250) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown()) +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H4.1_fig

ggsave(plot = var_H4.1_fig, filename = "Routput/Figures/fig_4_H4.1_var.tif", width = 8, height = 8, units="cm")

var_H4.1_fig2 <- posteriors_H4.1 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  scale_color_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  labs(x = expression("Proportion of variance (" * h^2 * " or " * Δ[evol] * bar(W) * ")"),
       y = "Density") +
  xlim(0, 1) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H4.1_fig2

ggsave(plot = var_H4.1_fig2, filename = "Routput/Figures/fig_4_H4.1_h2.tif", width = 8, height = 8, units="cm")

#Plotting density plot of breeding values
HW_4.1_breed <- HW_model4.1$Sol[,7:7657]
HW_4.1_BV <- posterior.mode(as.mcmc(HW_4.1_breed)) #breeding values for each individual
HW_4.1_BV_CI <- HPDinterval(as.mcmc(HW_4.1_breed)) #confidence intervals of breeding values
HW_4.1_BV_plot <- as.data.frame(HW_4.1_BV)
HW_4.1_BV_plot %>% 
  ggplot(aes(x=HW_4.1_BV)) +
  geom_density(fill="tomato2", color="tomato2", alpha=0.3) +
  labs(x="Ambient Fecundity Breeding Values", y="Frequency")  + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

#Model 4.2 with Fisher Parameter Expanded Prior (nu = 0.002 for prior sensitivity analysis)
heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior4.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

HW_model4.2 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "poisson", data = heated.Fecundity, prior = prior4.2,
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model4.2.RData")
plot(HW_model4.2$VCV) # good trace plots
summary(HW_model4.2) # good effective sample size
autocorr.diag(HW_model4.2$VCV) # no autocorrelation
heidel.diag(HW_model4.2$VCV) # convergence success
#No difference. Insensitive to prior change

#Model 4.3 - No dominance
load(file="Routput/Heated/HW_model4.3.RData")
plot(HW_model4.3$VCV) # good trace plots
summary(HW_model4.3) # good effective sample sizes
autocorr.diag(HW_model4.3$VCV) # no autocorrelation
heidel.diag(HW_model4.3$VCV) # convergence success

mean(HW_model4.3[["VCV"]][ , "animal"]) #Va
mean(HW_model4.3[["VCV"]][ , "matID"]) #Vm
mean(HW_model4.3[["VCV"]][ , "units"]) #Vr

#Transforming from latent to data scale using package "QGglmm"
df_H4.3 <- data.frame(
  va = as.vector(HW_model4.3[["VCV"]][, "animal"]),
  vm = as.vector(HW_model4.3[["VCV"]][, "matID"]),
  vr = as.vector(HW_model4.3[["VCV"]][, "units"]),
  vp = rowSums(HW_model4.3[["VCV"]]))
yhat_H4.3 <- predict(HW_model4.3, type = "terms")

#Va
Va_H4.3 <- do.call("rbind", apply(df_H4.3, 1, function(row){
  QGparams(predict = yhat_H4.3, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Va_H4.3$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_H4.3$var.a.obs))

#Vm
Vm_H4.3 <- do.call("rbind", apply(df_H4.3, 1, function(row){
  QGicc(predict = yhat_H4.3, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vm_H4.3$var.comp.obs))
HPDinterval(as.mcmc(Vm_H4.3$var.comp.obs))

#Vr
Vr_H4.3 <- do.call("rbind", apply(df_H4.3, 1, function(row){
  QGicc(predict = yhat_H4.3, 
        var.comp = row[["vr"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vr_H4.3$var.comp.obs))
HPDinterval(as.mcmc(Vr_H4.3$var.comp.obs))


#Model 4.4 - Full with different smaller V prior
heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior4.4 <- list(R = list(V = 0.01, nu = 0.002),
                 G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

HW_model4.4 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "poisson", data = heated.Fecundity, prior = prior4.4,
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model4.4.RData")
plot(HW_model4.4$VCV) # poor trace plots
summary(HW_model4.4) # low effective sample size
autocorr.diag(HW_model4.4$VCV) # high autocorrelation
heidel.diag(HW_model4.4$VCV) # convergence fail









#### No. 5 Overwintering Survival - Univariate | Bernoulli ####
#Note: For this trait, we include all plants from the study

#Model 5.0 with additive, dominance and maternal random effects + plot fixed effects
#We use a prior suggested by Villemereuil et al. 2013 for binary traits
heated.Germination <- MCMC.heated
prior5.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model5.0 <- MCMCglmm(germ ~ plot, random = ~animal + animalDom + matID, 
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "threshold", data = heated.Germination, prior = prior5.0, #Bernoulli distribution
                        nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Heated/HW_model5.0.RData")
plot(HW_model5.0$VCV) # good trace plots
summary(HW_model5.0) # good effective sample size
autocorr.diag(HW_model5.0$VCV) # no autocorrelation
heidel.diag(HW_model5.0$VCV) # convergence success

#Heritability
mean(HW_model5.0[["VCV"]][ , "animal"])
mean(HW_model5.0[["VCV"]][ , "animalDom"])
mean(HW_model5.0[["VCV"]][ , "matID"])
HW_herit5.0 <- HW_model5.0$VCV[, "animal"]/(HW_model5.0$VCV[, "animal"] + HW_model5.0$VCV[, "animalDom"] + HW_model5.0$VCV[, "matID"] + 1)
mean(HW_herit5.0)
posterior.mode(HW_herit5.0)
HPDinterval(HW_herit5.0)

#Transforming from latent to data scale using package "QGglmm"
df_H5.0 <- data.frame(
  va = as.vector(HW_model5.0[["VCV"]][, "animal"]),
  vd = as.vector(HW_model5.0[["VCV"]][, "animalDom"]),
  vm = as.vector(HW_model5.0[["VCV"]][, "matID"]),
  vp = rowSums(HW_model5.0[["VCV"]]))
yhat_H5.0 <- predict(HW_model5.0, type = "terms")

#Va
Va_H5.0 <- do.call("rbind", apply(df_H5.0, 1, function(row){
  QGparams(predict = yhat_H5.0, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Va_H5.0$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_H5.0$var.a.obs))

mean(as.mcmc(Va_H5.0$h2.obs)) #h2
HPDinterval(as.mcmc(Va_H5.0$h2.obs))

mu_H5.0 <- mean(as.mcmc(Va_H5.0$mean.obs))
mean(as.mcmc(Va_H5.0$h2.obs))/mu_H5.0 #Rw
HPDinterval(as.mcmc(Va_H5.0$h2.obs))/mu_H5.0

#Vd
Vd_H5.0 <- do.call("rbind", apply(df_H5.0, 1, function(row){
  QGicc(predict = yhat_H5.0, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vd_H5.0$var.comp.obs))
HPDinterval(as.mcmc(Vd_H5.0$var.comp.obs))

#Vm
Vm_H5.0 <- do.call("rbind", apply(df_H5.0, 1, function(row){
  QGicc(predict = yhat_H5.0, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_H5.0$var.comp.obs))
HPDinterval(as.mcmc(Vm_H5.0$var.comp.obs))


#Model 5.1 - No dominance
load(file="Routput/Heated/HW_model5.1.RData")
plot(HW_model5.1$VCV) # good trace plots
summary(HW_model5.1) # good effective sample sizes
autocorr.diag(HW_model5.1$VCV) # no autocorrelation
heidel.diag(HW_model5.1$VCV) # convergence success













#### No. 6 Spring-Summer Survival - Univariate | Bernoulli ####
#Note: For this trait, we include only include plants that germinated

#Model with additive, dominance and maternal random effects + plot fixed effects
#We use a prior suggested by Villemereuil et al. 2013 for binary traits
heated.Flowering <- MCMC.heated %>% filter(germ==1)
prior6.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model6.0 <- MCMCglmm(flower ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = heated.Flowering, prior = prior6.0, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Heated/HW_model6.0.RData")
plot(HW_model6.0$VCV) # good trace plots
summary(HW_model6.0) # good effective sample size
autocorr.diag(HW_model6.0$VCV) # no autocorrelation
heidel.diag(HW_model6.0$VCV) # convergence success

#Heritability
mean(HW_model6.0[["VCV"]][ , "animal"])
mean(HW_model6.0[["VCV"]][ , "animalDom"])
mean(HW_model6.0[["VCV"]][ , "matID"])
HW_herit6.0 <- HW_model6.0$VCV[, "animal"]/(HW_model6.0$VCV[, "animal"] + HW_model6.0$VCV[, "animalDom"] + HW_model6.0$VCV[, "matID"] + 1)
mean(HW_herit6.0)
posterior.mode(HW_herit6.0)
HPDinterval(HW_herit6.0)

#Transforming from latent to data scale using package "QGglmm"
df_H6.0 <- data.frame(
  va = as.vector(HW_model6.0[["VCV"]][, "animal"]),
  vd = as.vector(HW_model6.0[["VCV"]][, "animalDom"]),
  vm = as.vector(HW_model6.0[["VCV"]][, "matID"]),
  vp = rowSums(HW_model6.0[["VCV"]]))
yhat_H6.0 <- predict(HW_model6.0, type = "terms")

#Va
Va_H6.0 <- do.call("rbind", apply(df_H6.0, 1, function(row){
  QGparams(predict = yhat_H6.0, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Va_H6.0$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_H6.0$var.a.obs))

mean(as.mcmc(Va_H6.0$h2.obs)) #h2
HPDinterval(as.mcmc(Va_H6.0$h2.obs))

mu_H6.0 <- mean(as.mcmc(Va_H6.0$mean.obs))
mean(as.mcmc(Va_H6.0$h2.obs))/mu_H6.0 #Rw
HPDinterval(as.mcmc(Va_H6.0$h2.obs))/mu_H6.0

#Vd
Vd_H6.0 <- do.call("rbind", apply(df_H6.0, 1, function(row){
  QGicc(predict = yhat_H6.0, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vd_H6.0$var.comp.obs))
HPDinterval(as.mcmc(Vd_H6.0$var.comp.obs))

#Vm
Vm_H6.0 <- do.call("rbind", apply(df_H6.0, 1, function(row){
  QGicc(predict = yhat_H6.0, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_H6.0$var.comp.obs))
HPDinterval(as.mcmc(Vm_H6.0$var.comp.obs))


#Model 6.1 - No dominance
load(file="Routput/Heated/HW_model6.1.RData")
plot(HW_model6.1$VCV) # good trace plots
summary(HW_model6.1) # good effective sample sizes
autocorr.diag(HW_model6.1$VCV) # no autocorrelation
heidel.diag(HW_model6.1$VCV) # convergence success













#### No. 7 Fruiting Success - Univariate | Bernoulli ####
#Note: For this trait, we include only include plants that germinated and flowered

#Model 7.0 with additive, dominance and maternal random effects + plot fixed effects
#We use a prior suggested by de Villemereuil et al. 2013 for binary traits
heated.SeedMaturation <- MCMC.heated %>% filter(germ==1, flower==1)
prior7.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model7.0 <- MCMCglmm(seed ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = heated.SeedMaturation, prior = prior7.0, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Heated/HW_model7.0.RData")
plot(HW_model7.0$VCV) # good trace plots
summary(HW_model7.0) #good effective sample size
autocorr.diag(HW_model7.0$VCV) # no autocorrelation
heidel.diag(HW_model7.0$VCV) # convergence success

#Heriability
mean(HW_model7.0[["VCV"]][ , "animal"])
mean(HW_model7.0[["VCV"]][ , "animalDom"])
mean(HW_model7.0[["VCV"]][ , "matID"])
HW_herit7.0 <- HW_model7.0$VCV[, "animal"]/(HW_model7.0$VCV[, "animal"] + HW_model7.0$VCV[, "animalDom"] + HW_model7.0$VCV[, "matID"] + 1)
mean(HW_herit7.0)
posterior.mode(HW_herit7.0)
HPDinterval(HW_herit7.0)

#Transforming from latent to data scale using package "QGglmm"
df_H7.0 <- data.frame(
  va = as.vector(HW_model7.0[["VCV"]][, "animal"]),
  vd = as.vector(HW_model7.0[["VCV"]][, "animalDom"]),
  vm = as.vector(HW_model7.0[["VCV"]][, "matID"]),
  vp = rowSums(HW_model7.0[["VCV"]]))
yhat_H7.0 <- predict(HW_model7.0, type = "terms")

#Va
Va_H7.0 <- do.call("rbind", apply(df_H7.0, 1, function(row){
  QGparams(predict = yhat_H7.0, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Va_H7.0$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_H7.0$var.a.obs))

mean(as.mcmc(Va_H7.0$h2.obs)) #h2
HPDinterval(as.mcmc(Va_H7.0$h2.obs))

mu_H7.0 <- mean(as.mcmc(Va_H7.0$mean.obs))
mean(as.mcmc(Va_H7.0$h2.obs))/mu_H7.0 #Rw
HPDinterval(as.mcmc(Va_H7.0$h2.obs))/mu_H7.0

#Vd
Vd_H7.0 <- do.call("rbind", apply(df_H7.0, 1, function(row){
  QGicc(predict = yhat_H7.0, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vd_H7.0$var.comp.obs))
HPDinterval(as.mcmc(Vd_H7.0$var.comp.obs))

#Vm
Vm_H7.0 <- do.call("rbind", apply(df_H7.0, 1, function(row){
  QGicc(predict = yhat_H7.0, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_H7.0$var.comp.obs))
HPDinterval(as.mcmc(Vm_H7.0$var.comp.obs))

#Model 7.1 - No dominance
load(file="Routput/Heated/HW_model7.1.RData")
plot(HW_model7.1$VCV) # good trace plots
summary(HW_model7.1) # good effective sample size
autocorr.diag(HW_model7.1$VCV) # no autocorrelation
heidel.diag(HW_model7.1$VCV) # convergence success













#### No. 8 Leaf Number - Univariate | Gaussian ####
#Note: For leaf number, we only include plants that germinated

#Calculating Variance for leaf number for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1) %>%
  summarise(sd=sd(leaf), var=(sd)^2) #var = 6.928249 *******

#Model 8.0 with additive, dominance, and maternal random effects + plot fixed effects
#Prior using inverse gamma (0.001, 0.001)
heated.leaf <- MCMC.heated %>% filter(germ==1)
prior8.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model8.0 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = heated.leaf, prior = prior8.0, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model8.0.RData")
plot(HW_model8.0$VCV) # poor trace plots
summary(HW_model8.0) # bad effective sample size
autocorr.diag(HW_model8.0$VCV) # autocorrelation
heidel.diag(HW_model8.0$VCV) # convergence fail

#Model 8.1 with a Fisher Expanded Parameter Prior
heated.leaf <- MCMC.heated %>% filter(germ==1)
prior8.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model8.1 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = heated.leaf, prior = prior8.1, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model8.1.RData")
plot(HW_model8.1$VCV) # good trace plots
summary(HW_model8.1) # good effective sample size
autocorr.diag(HW_model8.1$VCV) # no autocorrelation
heidel.diag(HW_model8.1$VCV) # convergence success


#Latent/Data Scale
HW_herit8.1 <- HW_model8.1$VCV[, "animal"]/(HW_model8.1$VCV[, "animal"] + HW_model8.1$VCV[, "matID"] + HW_model8.1$VCV[, "animalDom"] + HW_model8.1$VCV[, "units"])
mean(HW_herit8.1)
posterio.mode(HW_herit8.1) #h2
HPDinterval(HW_herit8.1)

mean(HW_model8.1[["VCV"]][ , "animal"]) #Va
HPDinterval(HW_model8.1[["VCV"]][ , "animal"])

mean(HW_model8.1[["VCV"]][ , "animalDom"]) #Vd
HPDinterval(HW_model8.1[["VCV"]][ , "animalDom"]) 

mean(HW_model8.1[["VCV"]][ , "matID"]) #Vm
HPDinterval(HW_model8.1[["VCV"]][ , "matID"]) 

mean(HW_model8.1[["VCV"]][ , "units"]) #Vr
HPDinterval(HW_model8.1[["VCV"]][ , "units"]) 

mean(HW_model8.1$VCV[, "animal"] + HW_model8.1$VCV[, "matID"] + HW_model8.1$VCV[, "animalDom"]) #Vg
HPDinterval(HW_model8.1$VCV[, "animal"] + HW_model8.1$VCV[, "matID"] + HW_model8.1$VCV[, "animalDom"])

# Estimate trait mean as average across all plots (fixed) effects using posterior distributions
X_H8.1 <- HW_model8.1$X   # [n_obs x n_fixed]
fixed_effect_names_H8.1 <- colnames(X_H8.1) # fixed effect columns (intercept, plot2, plot3, etc)
Sol_fixed_H8.1 <- HW_model8.1$Sol[, fixed_effect_names_H8.1]  # [n_iter x n_fixed]; posterior subset to fixed effects only
fitted_vals_H8.1 <- X_H8.1 %*% t(Sol_fixed_H8.1) #design matrix multiplied by transpose of fixed effects
marginal_means_H8.1 <- colMeans(fitted_vals_H8.1) # average across individuals for each posterior sample
mean(marginal_means_H8.1) # summarize posterior mean and CI
HPDinterval(as.mcmc(marginal_means_H8.1))

mu_H8.1 <- mean(HW_model8.1[["Sol"]][,"(Intercept)"]) # mean based on intercept
mean(HW_model8.1[["VCV"]][ , "animal"]/marginal_means_H8.1) #Rw based on marginal means
HPDinterval(HW_model8.1[["VCV"]][ , "animal"]/marginal_means_H8.1)

# Check for significant difference between environments. Requires ambient data
mu_diff_8.1 <- marginal_means_H8.1 - marginal_means_A8.1
mean(mu_diff_8.1)
quantile(mu_diff_8.1, probs = c(0.025, 0.975))
mean(mu_diff_8.1 > 0)


#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_H8.1 <- data.frame(
  Va = as.vector(HW_model8.1[["VCV"]][ , "animal"]),
  Vd = as.vector(HW_model8.1[["VCV"]][ , "animalDom"]),
  Vm = as.vector(HW_model8.1[["VCV"]][ , "matID"]),
  Vr = as.vector(HW_model8.1[["VCV"]][ , "units"]),
  mu = as.vector(marginal_means_H8.1),
  h2 = as.vector(HW_model8.1$VCV[, "animal"]/(HW_model8.1$VCV[, "animal"] +
                                    HW_model8.1$VCV[, "matID"] + 
                                    HW_model8.1$VCV[, "animalDom"] + 
                                    HW_model8.1$VCV[, "units"])),
  Rw = as.vector(HW_model8.1[["VCV"]][ , "animal"]/marginal_means_H8.1))
posteriors_H8.1 <- posteriors_H8.1 %>% pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_H8.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(HW_model8.1[["VCV"]][ , "animal"])),
    mean(as.mcmc(HW_model8.1[["VCV"]][ , "animalDom"])),
    mean(as.mcmc(HW_model8.1[["VCV"]][ , "matID"])),
    mean(as.mcmc(HW_model8.1[["VCV"]][ , "units"])),
    mean(as.mcmc(marginal_means_H8.1)),
    mean(as.mcmc(HW_model8.1$VCV[, "animal"]/(HW_model8.1$VCV[, "animal"] +
                                        HW_model8.1$VCV[, "matID"] + 
                                        HW_model8.1$VCV[, "animalDom"] + 
                                        HW_model8.1$VCV[, "units"]))),
    mean(as.mcmc(HW_model8.1[["VCV"]][ , "animal"]/marginal_means_H8.1))
  ),
  Lower = c(
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "animal"]))[1],
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "animalDom"]))[1],
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "matID"]))[1],
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "units"]))[1],
    HPDinterval(as.mcmc(marginal_means_H8.1))[1],
    HPDinterval(as.mcmc(HW_model8.1$VCV[, "animal"]/(HW_model8.1$VCV[, "animal"] +
                                                  HW_model8.1$VCV[, "matID"] + 
                                                  HW_model8.1$VCV[, "animalDom"] + 
                                                  HW_model8.1$VCV[, "units"])))[1],
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "animal"]/marginal_means_H8.1))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "animal"]))[2],
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "animalDom"]))[2],
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "matID"]))[2],
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "units"]))[2],
    HPDinterval(as.mcmc(marginal_means_H8.1))[2],
    HPDinterval(as.mcmc(HW_model8.1$VCV[, "animal"]/(HW_model8.1$VCV[, "animal"] +
                                                       HW_model8.1$VCV[, "matID"] + 
                                                       HW_model8.1$VCV[, "animalDom"] + 
                                                       HW_model8.1$VCV[, "units"])))[2],
    HPDinterval(as.mcmc(HW_model8.1[["VCV"]][ , "animal"]/marginal_means_H8.1))[2]
  )
)

posteriors_H8.1$Component <- factor(posteriors_H8.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "Rw", "h2"))
save(posteriors_H8.1, stats_H8.1, file="Routput/Heated/HW_stats8.1.RData")
load("Routput/Heated/HW_stats8.1.RData")

var_H8.1_fig <- posteriors_H8.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, 4) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown()) +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H8.1_fig

ggsave(plot = var_H8.1_fig, filename = "Routput/Figures/fig_4_H8.1_var.tif", width = 8, height = 8, units="cm")

var_H8.1_fig2 <- posteriors_H8.1 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  scale_color_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  labs(x = expression("Proportion of variance (" * h^2 * " or " * Δ[evol] * bar(W) * ")"),
       y = "Density") +
  xlim(0, 0.5) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H8.1_fig2

ggsave(plot = var_H8.1_fig2, filename = "Routput/Figures/fig_4_H8.1_h2.tif", width = 8, height = 8, units="cm")



#Model 8.2 with a Fisher Expanded Parameter Prior (nu = 0.002 for prior sensitivity analysis)
heated.leaf <- MCMC.heated %>% filter(germ==1)
prior8.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

HW_model8.2 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "gaussian", data = heated.leaf, prior = prior8.2, #gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model8.2.RData")
plot(HW_model8.2$VCV) # good trace plots
summary(HW_model8.2) # good (but lower) effective sample sizes
autocorr.diag(HW_model8.2$VCV) # no autocorrelation
heidel.diag(HW_model8.2$VCV) # convergence success
#No difference. Insensitive to prior change.

#Model 8.3 - No dominance
load(file="Routput/Heated/HW_model8.3.RData")
plot(HW_model8.3$VCV) # good trace plots
summary(HW_model8.3) # good (but lower than expected) effective sample sizes
autocorr.diag(HW_model8.3$VCV) #no autocorrelation
heidel.diag(HW_model8.3$VCV) # convergence success

mean(HW_model8.3[["VCV"]][ , "animal"]) #Va
HPDinterval(HW_model8.3[["VCV"]][ , "animal"])

mean(HW_model8.3[["VCV"]][ , "matID"]) #Vm
HPDinterval(HW_model8.3[["VCV"]][ , "matID"]) 

mean(HW_model8.3[["VCV"]][ , "units"]) #Vr
HPDinterval(HW_model8.3[["VCV"]][ , "units"]) 

#Model 8.4 - Full with different smaller V prior
heated.leaf <- MCMC.heated %>% filter(germ==1)
prior8.4 <- list(R = list(V = 0.01, nu = 0.002),
                 G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

HW_model8.4 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "gaussian", data = heated.leaf, prior = prior8.4, #gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model8.4.RData")
plot(HW_model8.4$VCV) # poor trace plots
summary(HW_model8.4) # low effective sample size
autocorr.diag(HW_model8.4$VCV) # high autocorrelation
heidel.diag(HW_model8.4$VCV) # convergence fail











#### No. 9 Height - Univariate | Gaussian ####
#Note: For height, we only include plants that germinated

#Calculating Variance for height for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1, !height==0) %>%
  summarise(sd=sd(height), var=(sd)^2) #var = 213.7212 *******

#Model 9.0 with additive, dominance, and maternal random effects + plot fixed effects
#Prior with inverse gamma (0.001, 0.001)
heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior9.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model9.0 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = heated.height, prior = prior9.0, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model9.0.RData")
plot(HW_model9.0$VCV) # poor trace plots
summary(HW_model9.0) # poor effective sample size
autocorr.diag(HW_model9.0$VCV) # high autocorrelation
heidel.diag(HW_model9.0$VCV) # convergence fail

#Model 9.1 with a Fisher Parameter Expanded Prior
heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior9.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model9.1 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = heated.height, prior = prior9.1, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model9.1.RData")
plot(HW_model9.1$VCV) # good trace plots
summary(HW_model9.1) # good effective sample size
autocorr.diag(HW_model9.1$VCV) # no autocorrelation
heidel.diag(HW_model9.1$VCV) # convergence success

#Latent/Data Scale
HW_herit9.1 <- HW_model9.1$VCV[, "animal"]/(HW_model9.1$VCV[, "animal"] + HW_model9.1$VCV[, "animalDom"] + 
                                              HW_model9.1$VCV[, "matID"] + HW_model9.1$VCV[, "units"])
mean(HW_herit9.1)
posterior.mode(HW_herit9.1) #h2
HPDinterval(HW_herit9.1)

mean(HW_model9.1[["VCV"]][ , "animal"]) #Va
HPDinterval(HW_model9.1[["VCV"]][ , "animal"]) 

mean(HW_model9.1[["VCV"]][ , "animalDom"]) #Vd
HPDinterval(HW_model9.1[["VCV"]][ , "animalDom"]) 

mean(HW_model9.1[["VCV"]][ , "matID"]) #Vm
HPDinterval(HW_model9.1[["VCV"]][ , "matID"]) 

mean(HW_model9.1[["VCV"]][ , "units"]) #Vr
HPDinterval(HW_model9.1[["VCV"]][ , "units"]) 

mean(HW_model9.1$VCV[, "animal"] + HW_model9.1$VCV[, "animalDom"] + HW_model9.1$VCV[, "matID"]) #Vg
HPDinterval(HW_model9.1$VCV[, "animal"] + HW_model9.1$VCV[, "animalDom"] + HW_model9.1$VCV[, "matID"])


# Estimate trait mean as average across all plots (fixed) effects using posterior distributions
X_H9.1 <- HW_model9.1$X    # [n_obs x n_fixed]
fixed_effect_names_H9.1 <- colnames(X_H9.1) # fixed effect columns (intercept, plot2, plot3, etc)
Sol_fixed_H9.1 <- HW_model9.1$Sol[, fixed_effect_names_H9.1]  # [n_iter x n_fixed]; posterior subset to fixed effects only
fitted_vals_H9.1 <- X_H9.1 %*% t(Sol_fixed_H9.1) #design matrix multiplied by transpose of fixed effects
marginal_means_H9.1 <- colMeans(fitted_vals_H9.1) # average across individuals for each posterior sample
mean(marginal_means_H9.1) # summarize posterior mean and CI
HPDinterval(as.mcmc(marginal_means_H9.1))

mu_H9.1 <- mean(HW_model9.1[["Sol"]][,"(Intercept)"]) # mu based on intercept
mean(HW_model9.1[["VCV"]][ , "animal"]/marginal_means_H9.1) #Rw based on marginal means
HPDinterval(HW_model9.1[["VCV"]][ , "animal"]/marginal_means_H9.1)

# Check for significant difference between environments. Requires ambient data
mu_diff_9.1 <- marginal_means_H9.1 - marginal_means_A9.1
mean(mu_diff_9.1)
quantile(mu_diff_9.1, probs = c(0.025, 0.975))
mean(mu_diff_9.1 > 0)

#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_H9.1 <- data.frame(
  Va = as.vector(HW_model9.1[["VCV"]][ , "animal"]),
  Vd = as.vector(HW_model9.1[["VCV"]][ , "animalDom"]),
  Vm = as.vector(HW_model9.1[["VCV"]][ , "matID"]),
  Vr = as.vector(HW_model9.1[["VCV"]][ , "units"]),
  mu = as.vector(marginal_means_H9.1),
  h2 = as.vector(HW_model9.1$VCV[, "animal"]/(HW_model9.1$VCV[, "animal"] +
                                                HW_model9.1$VCV[, "matID"] + 
                                                HW_model9.1$VCV[, "animalDom"] + 
                                                HW_model9.1$VCV[, "units"])),
  Rw = as.vector(HW_model9.1[["VCV"]][ , "animal"]/marginal_means_H9.1))
posteriors_H9.1 <- posteriors_H9.1 %>% pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_H9.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(HW_model9.1[["VCV"]][ , "animal"])),
    mean(as.mcmc(HW_model9.1[["VCV"]][ , "animalDom"])),
    mean(as.mcmc(HW_model9.1[["VCV"]][ , "matID"])),
    mean(as.mcmc(HW_model9.1[["VCV"]][ , "units"])),
    mean(as.mcmc(marginal_means_H9.1)),
    mean(as.mcmc(HW_model9.1$VCV[, "animal"]/(HW_model9.1$VCV[, "animal"] +
                                                HW_model9.1$VCV[, "matID"] + 
                                                HW_model9.1$VCV[, "animalDom"] + 
                                                HW_model9.1$VCV[, "units"]))),
    mean(as.mcmc(HW_model9.1[["VCV"]][ , "animal"]/marginal_means_H9.1))
  ),
  Lower = c(
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "animal"]))[1],
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "animalDom"]))[1],
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "matID"]))[1],
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "units"]))[1],
    HPDinterval(as.mcmc(marginal_means_H9.1))[1],
    HPDinterval(as.mcmc(HW_model9.1$VCV[, "animal"]/(HW_model9.1$VCV[, "animal"] +
                                                       HW_model9.1$VCV[, "matID"] + 
                                                       HW_model9.1$VCV[, "animalDom"] + 
                                                       HW_model9.1$VCV[, "units"])))[1],
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "animal"]/marginal_means_H9.1))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "animal"]))[2],
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "animalDom"]))[2],
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "matID"]))[2],
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "units"]))[2],
    HPDinterval(as.mcmc(marginal_means_H9.1))[2],
    HPDinterval(as.mcmc(HW_model9.1$VCV[, "animal"]/(HW_model9.1$VCV[, "animal"] +
                                                       HW_model9.1$VCV[, "matID"] + 
                                                       HW_model9.1$VCV[, "animalDom"] + 
                                                       HW_model9.1$VCV[, "units"])))[2],
    HPDinterval(as.mcmc(HW_model9.1[["VCV"]][ , "animal"]/marginal_means_H9.1))[2]
  )
)

posteriors_H9.1$Component <- factor(posteriors_H9.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "Rw", "h2"))
save(posteriors_H9.1, stats_H9.1, file="Routput/Heated/HW_stats9.1.RData")
load("Routput/Heated/HW_stats9.1.RData")

var_H9.1_fig <- posteriors_H9.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, 45) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown()) +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H9.1_fig

ggsave(plot = var_H9.1_fig, filename = "Routput/Figures/fig_4_H9.1_var.tif", width = 8, height = 8, units="cm")

var_H9.1_fig2 <- posteriors_H9.1 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  scale_color_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  labs(x = expression("Proportion of variance (" * h^2 * " or " * Δ[evol] * bar(W) * ")"),
       y = "Density") +
  xlim(0, 1) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H9.1_fig2

ggsave(plot = var_H9.1_fig2, filename = "Routput/Figures/fig_4_H9.1_h2.tif", width = 8, height = 8, units="cm")


#Model 9.2 with a Fisher Parameter Expanded Prior (nu = 0.002 for prior sensitivity analysis)
heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior9.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

HW_model9.2 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "gaussian", data = heated.height, prior = prior9.2, #gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model9.2.RData")
plot(HW_model9.2$VCV) # ok (but worse) trace plots
summary(HW_model9.2) # poor effective sample size
autocorr.diag(HW_model9.2$VCV) # autocorrelation
heidel.diag(HW_model9.2$VCV) # convergence fail
#Sensitive prior. Model 9.1 is better fit. 

#Model 9.3 - No dominance
load(file="Routput/Heated/HW_model9.3.RData")
plot(HW_model9.3$VCV) # good trace plots
summary(HW_model9.3) # good effective sample size
autocorr.diag(HW_model9.3$VCV) # no autocorrelation
heidel.diag(HW_model9.3$VCV) # convergence success

mean(HW_model9.3[["VCV"]][ , "animal"]) #Va
HPDinterval(HW_model9.3[["VCV"]][ , "animal"])

mean(HW_model9.3[["VCV"]][ , "matID"]) #Vm
HPDinterval(HW_model9.3[["VCV"]][ , "matID"]) 

mean(HW_model9.3[["VCV"]][ , "units"]) #Vr
HPDinterval(HW_model9.3[["VCV"]][ , "units"]) 


#Model 9.4 - Full with different smaller V prior
heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior9.4 <- list(R = list(V = 0.01, nu = 0.002),
                 G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

HW_model9.4 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "gaussian", data = heated.height, prior = prior9.4, #gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model9.4.RData")
plot(HW_model9.4$VCV) # poor trace plots
summary(HW_model9.4) # low effective sample size
autocorr.diag(HW_model9.4$VCV) # high autocorrelation
heidel.diag(HW_model9.4$VCV) # convergence fail












#### No. 10 Stem Diameter - Univariate | Gaussian ####
#Note: For stem diameter, we only include plants that germinated

#Calculating Variance for stem diameter for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1, !stem_diam==0) %>%
  summarise(mean=mean(stem_diam), sd=sd(stem_diam), var=(sd)^2) #var = 3.102643 *******

#Model 10.0 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior10.0 <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 0.002),
                            G2 = list(V = 1, nu = 0.002), 
                            G3 = list(V = 1, nu = 0.002)))

HW_model10.0 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                       ginverse = list(animal = Ainv, animalDom = Dinv),
                       family = "gaussian", data = heated.stem, prior = prior10.0, #gaussian distribution
                       nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model10.0.RData")
plot(HW_model10.0$VCV) # poor trace plots
summary(HW_model10.0) # poor effective sample size
autocorr.diag(HW_model10.0$VCV) # high autocorrelation
heidel.diag(HW_model10.0$VCV) # convergence fail

#Model 10.1 using a Fisher Parameter Expanded Prior
heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior10.1 <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                            G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                            G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model10.1 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                       ginverse = list(animal = Ainv, animalDom = Dinv),
                       family = "gaussian", data = heated.stem, prior = prior10.1, #gaussian distribution
                       nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model10.1.RData")
plot(HW_model10.1$VCV) # good trace plots
summary(HW_model10.1) # good effective sample size
autocorr.diag(HW_model10.1$VCV) # no autocorrelation
heidel.diag(HW_model10.1$VCV) # convergence success

#Latent/Data Scale
HW_herit10.1 <- HW_model10.1$VCV[, "animal"]/(HW_model10.1$VCV[, "animal"] + HW_model10.1$VCV[, "animalDom"] + HW_model10.1$VCV[, "matID"] + HW_model10.1$VCV[, "units"])
mean(HW_herit10.1)
posterior.mode(HW_herit10.1) #h2
HPDinterval(HW_herit10.1)

mean(HW_model10.1[["VCV"]][ , "animal"]) #Va
HPDinterval(HW_model10.1[["VCV"]][ , "animal"])

mean(HW_model10.1[["VCV"]][ , "animalDom"]) #Vd
HPDinterval(HW_model10.1[["VCV"]][ , "animalDom"])

mean(HW_model10.1[["VCV"]][ , "matID"]) #Vm
HPDinterval(HW_model10.1[["VCV"]][ , "matID"])

mean(HW_model10.1[["VCV"]][ , "units"]) #Vr
HPDinterval(HW_model10.1[["VCV"]][ , "units"]) 

mean(HW_model10.1$VCV[, "animal"] + HW_model10.1$VCV[, "animalDom"] + HW_model10.1$VCV[, "matID"]) #Vg
HPDinterval(HW_model10.1$VCV[, "animal"] + HW_model10.1$VCV[, "animalDom"] + HW_model10.1$VCV[, "matID"])


# Estimate trait mean as average across all plots (fixed) effects using posterior distributions
X_H10.1 <- HW_model10.1$X    # [n_obs x n_fixed]
fixed_effect_names_H10.1 <- colnames(X_H10.1) # fixed effect columns (intercept, plot2, plot3, etc)
Sol_fixed_H10.1 <- HW_model10.1$Sol[, fixed_effect_names_H10.1]  # [n_iter x n_fixed]; posterior subset to fixed effects only
fitted_vals_H10.1 <- X_H10.1 %*% t(Sol_fixed_H10.1) #design matrix multiplied by transpose of fixed effects
marginal_means_H10.1 <- colMeans(fitted_vals_H10.1) # average across individuals for each posterior sample
mean(marginal_means_H10.1) # summarize posterior mean and CI
HPDinterval(as.mcmc(marginal_means_H10.1))

mu_H10.1 <- mean(HW_model10.1[["Sol"]][,"(Intercept)"]) # mu based on intercept
mean(HW_model10.1[["VCV"]][ , "animal"]/marginal_means_H10.1) #Rw based on marginal means
HPDinterval(HW_model10.1[["VCV"]][ , "animal"]/marginal_means_H10.1)

# Check for significant difference between environments. Requires ambient data
mu_diff_10.1 <- marginal_means_H10.1 - marginal_means_A10.1
mean(mu_diff_10.1)
quantile(mu_diff_10.1, probs = c(0.025, 0.975))
mean(mu_diff_10.1 > 0)


#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_H10.1 <- data.frame(
  Va = as.vector(HW_model10.1[["VCV"]][ , "animal"]),
  Vd = as.vector(HW_model10.1[["VCV"]][ , "animalDom"]),
  Vm = as.vector(HW_model10.1[["VCV"]][ , "matID"]),
  Vr = as.vector(HW_model10.1[["VCV"]][ , "units"]),
  mu = as.vector(marginal_means_H10.1),
  h2 = as.vector(HW_model10.1$VCV[, "animal"]/(HW_model10.1$VCV[, "animal"] +
                                                HW_model10.1$VCV[, "matID"] + 
                                                HW_model10.1$VCV[, "animalDom"] + 
                                                HW_model10.1$VCV[, "units"])),
  Rw = as.vector(HW_model10.1[["VCV"]][ , "animal"]/marginal_means_H10.1))
posteriors_H10.1 <- posteriors_H10.1 %>% pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_H10.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(HW_model10.1[["VCV"]][ , "animal"])),
    mean(as.mcmc(HW_model10.1[["VCV"]][ , "animalDom"])),
    mean(as.mcmc(HW_model10.1[["VCV"]][ , "matID"])),
    mean(as.mcmc(HW_model10.1[["VCV"]][ , "units"])),
    mean(as.mcmc(marginal_means_H10.1)),
    mean(as.mcmc(HW_model10.1$VCV[, "animal"]/(HW_model10.1$VCV[, "animal"] +
                                                HW_model10.1$VCV[, "matID"] + 
                                                HW_model10.1$VCV[, "animalDom"] + 
                                                HW_model10.1$VCV[, "units"]))),
    mean(as.mcmc(HW_model10.1[["VCV"]][ , "animal"]/marginal_means_H10.1))
  ),
  Lower = c(
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "animal"]))[1],
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "animalDom"]))[1],
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "matID"]))[1],
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "units"]))[1],
    HPDinterval(as.mcmc(marginal_means_H10.1))[1],
    HPDinterval(as.mcmc(HW_model10.1$VCV[, "animal"]/(HW_model10.1$VCV[, "animal"] +
                                                       HW_model10.1$VCV[, "matID"] + 
                                                       HW_model10.1$VCV[, "animalDom"] + 
                                                       HW_model10.1$VCV[, "units"])))[1],
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "animal"]/marginal_means_H10.1))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "animal"]))[2],
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "animalDom"]))[2],
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "matID"]))[2],
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "units"]))[2],
    HPDinterval(as.mcmc(marginal_means_H10.1))[2],
    HPDinterval(as.mcmc(HW_model10.1$VCV[, "animal"]/(HW_model10.1$VCV[, "animal"] +
                                                       HW_model10.1$VCV[, "matID"] + 
                                                       HW_model10.1$VCV[, "animalDom"] + 
                                                       HW_model10.1$VCV[, "units"])))[2],
    HPDinterval(as.mcmc(HW_model10.1[["VCV"]][ , "animal"]/marginal_means_H10.1))[2]
  )
)

posteriors_H10.1$Component <- factor(posteriors_H10.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "Rw", "h2"))
save(posteriors_H10.1, stats_H10.1, file="Routput/Heated/HW_stats10.1.RData")
load("Routput/Heated/HW_stats10.1.RData")

var_H10.1_fig <- posteriors_H10.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, 1) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown()) +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H10.1_fig

ggsave(plot = var_H10.1_fig, filename = "Routput/Figures/fig_4_H10.1_var.tif", width = 8, height = 8, units="cm")

var_H10.1_fig2 <- posteriors_H10.1 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  scale_color_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  labs(x = expression("Proportion of variance (" * h^2 * " or " * Δ[evol] * bar(W) * ")"),
       y = "Density") +
  xlim(0, .5) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H10.1_fig2

ggsave(plot = var_H10.1_fig2, filename = "Routput/Figures/fig_4_H10.1_h2.tif", width = 8, height = 8, units="cm")


#Model 10.2 using a Fisher Parameter Expanded Prior (nu = 0.002 for prior sensitivity analysis)
heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior10.2 <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

HW_model10.2 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                         ginverse = list(animal = Ainv, animalDom = Dinv),
                         family = "gaussian", data = heated.stem, prior = prior10.2, #gaussian distribution
                         nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model10.2.RData")
plot(HW_model10.2$VCV) # good trace plots
summary(HW_model10.2) # good effective sample size
autocorr.diag(HW_model10.2$VCV) # no autocorrelation
heidel.diag(HW_model10.2$VCV) # convergence success
#No difference. Insensitive to prior change. 

#Model 10.3 - No dominance
load(file="Routput/Heated/HW_model10.3.RData")
plot(HW_model10.3$VCV) # good trace plots
summary(HW_model10.3) # good effective sample sizes
autocorr.diag(HW_model10.3$VCV) # no autocorrelation
heidel.diag(HW_model10.3$VCV) # convergence success

mean(HW_model10.3[["VCV"]][ , "animal"]) #Va
HPDinterval(HW_model10.3[["VCV"]][ , "animal"])

mean(HW_model10.3[["VCV"]][ , "matID"]) #Vm
HPDinterval(HW_model10.3[["VCV"]][ , "matID"])

mean(HW_model10.3[["VCV"]][ , "units"]) #Vr
HPDinterval(HW_model10.3[["VCV"]][ , "units"]) 


#Model 10.4 - Full with different smaller V prior
heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior10.4 <- list(R = list(V = 0.01, nu = 0.002),
                  G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

HW_model10.4 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                         ginverse = list(animal = Ainv, animalDom = Dinv),
                         family = "gaussian", data = heated.stem, prior = prior10.4, #gaussian distribution
                         nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model10.4.RData")
plot(HW_model10.4$VCV) # poor trace plots
summary(HW_model10.4) # low effective sample size
autocorr.diag(HW_model10.4$VCV) # high autocorrelation
heidel.diag(HW_model10.4$VCV) # convergence fail








#### No. 11 Number of Flowering Clusters - Univariate | Poisson ####
#Note: For Flowering Clusters, we only include plants that survived to reach flowering

#Calculating Variance for Flowering Clusters for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1 & flower==1) %>%
  summarise(mean=mean(flwr_clstr), sd=sd(flwr_clstr), var=(sd)^2) #var = 9.736969 

#Model 11.0 with additive, dominance, and maternal random effects + plot fixed effects
#Prior with inverse gamma (0.001, 0.001)
heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior11.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model11.0 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "poisson", data = heated.flwrclstr, prior = prior11.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model11.0.RData")
plot(HW_model11.0$VCV) # poor trace plot
summary(HW_model11.0) # poor effective sample size
autocorr.diag(HW_model11.0$VCV) # high autocorrelation
heidel.diag(HW_model11.0$VCV) # convergence fail

#Model 11.1 with Fisher Parameter Expanded Prior
heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
hist(heated.flwrclstr$flwr_clstr)

prior11.1 <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model11.1 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                      ginverse = list(animal = Ainv, animalDom = Dinv),
                      family = "poisson", data = heated.flwrclstr, prior = prior11.1,
                      nitt = 8500000, thin = 4000, burnin = 500000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model11.1-extend2.RData")
plot(HW_model11.1$VCV) # good trace plots
summary(HW_model11.1) # sufficient effective sample sizes
autocorr.diag(HW_model11.1$VCV) # moderate autocorrelation for dominance and residual variance
heidel.diag(HW_model11.1$VCV) # convergence success

#HW_model11.1 has a nitt of 3100000, thin = 1500
#HW_model11.1-extend has a nitt of 5100000, thin = 2500
#HW_model11.1-extend2 has nitt of 8500000, thin = 4000, burnin 500000
#HW_model11.1.priorAlpV1 has alpha.V = 1 

#Latent Scale
mean(HW_model11.1[["VCV"]][ , "animal"]) #Va
mean(HW_model11.1[["VCV"]][ , "animalDom"]) #Vd
mean(HW_model11.1[["VCV"]][ , "matID"]) #Vm
mean(HW_model11.1[["VCV"]][ , "units"]) #Vr
HW_herit11.1 <- HW_model11.1$VCV[, "animal"]/(HW_model11.1$VCV[, "animal"] + HW_model11.1$VCV[, "animalDom"] + HW_model11.1$VCV[, "matID"] + HW_model11.1$VCV[, "units"])
mean(HW_herit11.1) #h2
posterior.mode(HW_herit11.1)
HPDinterval(HW_herit11.1)

mean(HW_model11.1$VCV[, "animal"] + HW_model11.1$VCV[, "animalDom"] + HW_model11.1$VCV[, "matID"]) #Vg
HPDinterval(HW_model11.1$VCV[, "animal"] + HW_model11.1$VCV[, "animalDom"] + HW_model11.1$VCV[, "matID"])

#Transforming from latent to data scale using package "QGglmm"
df_H11.1 <- data.frame(
  va = as.vector(HW_model11.1[["VCV"]][, "animal"]),
  vd = as.vector(HW_model11.1[["VCV"]][, "animalDom"]),
  vm = as.vector(HW_model11.1[["VCV"]][, "matID"]),
  vr = as.vector(HW_model11.1[["VCV"]][, "units"]),
  vp = rowSums(HW_model11.1[["VCV"]]))
yhat_H11.1 <- predict(HW_model11.1, type = "terms")

#Va
Va_H11.1 <- do.call("rbind", apply(df_H11.1, 1, function(row){
  QGparams(predict = yhat_H11.1, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Va_H11.1$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_H11.1$var.a.obs))

mean(as.mcmc(Va_H11.1$h2.obs)) #h2
HPDinterval(as.mcmc(Va_H11.1$h2.obs))

mu_H11.1 <- mean(as.mcmc(Va_H11.1$mean.obs)) # mu
mu_H11.1
HPDinterval(as.mcmc(Va_H11.1$mean.obs))

# Check for significant difference between environments. Requires ambient data
mu_diff_11.1 <- Va_H11.1$mean.obs - Va_A11.1$mean.obs
mean(mu_diff_11.1)
quantile(mu_diff_11.1, probs = c(0.025, 0.975))
mean(mu_diff_11.1 > 0)


mean(as.mcmc(Va_H11.1$var.a.obs))/mu_H11.1 #Rw
(HPDinterval(as.mcmc(Va_H11.1$var.a.obs)))/mu_H11.1

#Vd
Vd_H11.1 <- do.call("rbind", apply(df_H11.1, 1, function(row){
  QGicc(predict = yhat_H11.1, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vd_H11.1$var.comp.obs))
HPDinterval(as.mcmc(Vd_H11.1$var.comp.obs))

#Vm
Vm_H11.1 <- do.call("rbind", apply(df_H11.1, 1, function(row){
  QGicc(predict = yhat_H11.1, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vm_H11.1$var.comp.obs))
HPDinterval(as.mcmc(Vm_H11.1$var.comp.obs))

#Vr
Vr_H11.1 <- do.call("rbind", apply(df_H11.1, 1, function(row){
  QGicc(predict = yhat_H11.1, 
        var.comp = row[["vr"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vr_H11.1$var.comp.obs))
HPDinterval(as.mcmc(Vr_H11.1$var.comp.obs))

#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_H11.1 <- data.frame(
  Va = Va_H11.1$var.a.obs,
  Vd = Vd_H11.1$var.comp.obs,
  Vm = Vm_H11.1$var.comp.obs,
  Vr = Vr_H11.1$var.comp.obs,
  mu = Va_H11.1$mean.obs,
  h2 = Va_H11.1$h2.obs,
  Rw = (Va_H11.1$var.a.obs / Va_H11.1$mean.obs)) %>%
  pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_H11.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(Va_H11.1$var.a.obs)),
    mean(as.mcmc(Vd_H11.1$var.comp.obs)),
    mean(as.mcmc(Vm_H11.1$var.comp.obs)),
    mean(as.mcmc(Vr_H11.1$var.comp.obs)),
    mean(as.mcmc(Va_H11.1$mean.obs)),
    mean(as.mcmc(Va_H11.1$h2.obs)),
    mean(as.mcmc(Va_H11.1$var.a.obs / Va_H11.1$mean.obs))
  ),
  Lower = c(
    HPDinterval(as.mcmc(Va_H11.1$var.a.obs))[1],
    HPDinterval(as.mcmc(Vd_H11.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vm_H11.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vr_H11.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Va_H11.1$mean.obs))[1],
    HPDinterval(as.mcmc(Va_H11.1$h2.obs))[1],
    HPDinterval(as.mcmc(Va_H11.1$var.a.obs / Va_H11.1$mean.obs))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(Va_H11.1$var.a.obs))[2],
    HPDinterval(as.mcmc(Vd_H11.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vm_H11.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vr_H11.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Va_H11.1$mean.obs))[2],
    HPDinterval(as.mcmc(Va_H11.1$h2.obs))[2],
    HPDinterval(as.mcmc(Va_H11.1$var.a.obs / Va_H11.1$mean.obs))[2]
  )
)

posteriors_H11.1$Component <- factor(posteriors_H11.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "h2", "Rw"))
save(posteriors_H11.1, stats_H11.1, file="Routput/Heated/HW_stats11.1.RData")
load("Routput/Heated/HW_stats11.1.RData")

var_H11.1_fig <- posteriors_H11.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "#ee5c42")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, 3) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown()) +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H11.1_fig

ggsave(plot = var_H11.1_fig, filename = "Routput/Figures/fig_4_H11.1_var.tif", width = 8, height = 8, units="cm")

var_H11.1_fig2 <- posteriors_H11.1 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  scale_color_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  labs(x = expression("Proportion of variance (" * h^2 * " or " * Δ[evol] * bar(W) * ")"),
       y = "Density") +
  xlim(0, 1) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) # remove x-y labels, y-ticks
var_H11.1_fig2

ggsave(plot = var_H11.1_fig2, filename = "Routput/Figures/fig_4_H11.1_h2.tif", width = 8, height = 8, units="cm")


#Model 11.2 with Fisher Parameter Expanded Prior (nu = 0.002 for prior sensitivity analysis + longer nitt)
heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior11.2 <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1),
                           G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1), 
                           G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1)))

HW_model11.2 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                         ginverse = list(animal = Ainv, animalDom = Dinv),
                         family = "poisson", data = heated.flwrclstr, prior = prior11.2,
                         nitt = 4250000, thin = 2000, burnin = 250000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model11.2.RData")
plot(HW_model11.2$VCV) # good trace plots
summary(HW_model11.2) # slightly higher effective sample sizes (without higher nitt)
autocorr.diag(HW_model11.2$VCV) # high autocorrelation for animal and dominance
heidel.diag(HW_model11.2$VCV) # convergence success
#No much difference. Low sensitivity to prior change. 

#Model 11.3 - No dominance
load(file="Routput/Heated/HW_model11.3.RData")
plot(HW_model11.3$VCV) # good trace plots
summary(HW_model11.3) # good effective sample sizes
autocorr.diag(HW_model11.3$VCV) # no autocorrelation
heidel.diag(HW_model11.3$VCV) # convergence success

mean(HW_model11.3[["VCV"]][ , "animal"]) #Va
mean(HW_model11.3[["VCV"]][ , "matID"]) #Vm
mean(HW_model11.3[["VCV"]][ , "units"]) #Vr


#Transforming from latent to data scale using package "QGglmm"
df_H11.3 <- data.frame(
  va = as.vector(HW_model11.3[["VCV"]][, "animal"]),
  vm = as.vector(HW_model11.3[["VCV"]][, "matID"]),
  vr = as.vector(HW_model11.3[["VCV"]][, "units"]),
  vp = rowSums(HW_model11.3[["VCV"]]))
yhat_H11.3 <- predict(HW_model11.3, type = "terms")

#Va
Va_H11.3 <- do.call("rbind", apply(df_H11.3, 1, function(row){
  QGparams(predict = yhat_H11.3, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Va_H11.3$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_H11.3$var.a.obs))

#Vm
Vm_H11.3 <- do.call("rbind", apply(df_H11.3, 1, function(row){
  QGicc(predict = yhat_H11.3, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vm_H11.3$var.comp.obs))
HPDinterval(as.mcmc(Vm_H11.3$var.comp.obs))

#Vr
Vr_H11.3 <- do.call("rbind", apply(df_H11.3, 1, function(row){
  QGicc(predict = yhat_H11.3, 
        var.comp = row[["vr"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vr_H11.3$var.comp.obs))
HPDinterval(as.mcmc(Vr_H11.3$var.comp.obs))


#Model 11.4 - Full with different smaller V prior
heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior11.4 <- list(R = list(V = 0.01, nu = 0.002),
                  G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1),
                           G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1), 
                           G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1)))

HW_model11.4 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                         ginverse = list(animal = Ainv, animalDom = Dinv),
                         family = "poisson", data = heated.flwrclstr, prior = prior11.4,
                         nitt = 4250000, thin = 2000, burnin = 250000, verbose = T, pr = TRUE)

load(file="Routput/Heated/HW_model11.4.RData")
plot(HW_model11.4$VCV) # poor trace plots
summary(HW_model11.4) # low effective sample size
autocorr.diag(HW_model11.4$VCV) # high autocorrelation
heidel.diag(HW_model11.4$VCV) # convergence fail
