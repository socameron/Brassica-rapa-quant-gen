#### PROJECT: Brassica rapa Va/W Study (Data collected by at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Construct within ambient environment models using MCMCglmm to extract Va, Vd, and Vm

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
                             height = "d", stem_diam = "d", germ_census = "d", flwr_census = "d")

MCMC.2019 <- read_csv("Rdata/MCMC_2019_cleaned.csv", col_names = TRUE, na = "NA", 
                      col_types=col_types_list2)
spec(MCMC.2019)

ped <- read.csv("Rdata/heatarrays_animal_pedigree.csv") #note don't use readr to import the csv file.
for (x in 1:3) ped[, x] <- as.factor(ped[, x])

#Changing tibbles into dataframes because MCMCglmm cannot read tibbles
MCMC.2019 <- as.data.frame(MCMC.2019)
ped <- as.data.frame(ped)

MCMC.ambient <- MCMC.2019 %>%
  filter(treatment=="A")

MCMC.ambient$animal <- as.factor(MCMC.ambient$animal) #double check
lapply(ped, class)
lapply(MCMC.ambient, class)

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
##############      MCMCglmm ANALYSIS FOR AMBIENT      ###############
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
MCMC.ambient$animalDom <- MCMC.ambient$animal















#### No. 1 Lifetime Fitness - "Univariate" | Zero-Inflated Poisson ####
#We set two different priors for the zero inflation process, and the Poisson process. However, according to Hadfield (https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/018802.html), the genetic correlation between the two processes must be set to 0 and cannot be estimated. 

MCMC.ambient %>%
  summarise(mean=mean(seed_pods), sd=sd(seed_pods), var=(sd)^2) #var = 136.3176

#Model 1.0 with additive + dominance + maternal random effects and plot fixed effects
#Zero-Inflated Poisson distribution + 2 different priors
#Level 1 = Poisson process, Level 2 = Zero inflation
#NOTE: the zero-inflation model is similar to breaking up the model into survival and fecundity. The ZI and PP are considered separate 'traits'. 
ambient.total <- MCMC.ambient
prior1.0 <- list(R = list(V = diag(2), nu = 1, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1 - Zero Inflation
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1 - Zero Inflation
                          G5 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G6 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) # Chi Square Prior with SD 1 - Zero Inflation

AM_model1.0 <- MCMCglmm(seed_pods ~ trait - 1 + plot,
                     random = ~idh(at.level(trait,1)):animal + idh(at.level(trait,2)):animal +
                       idh(at.level(trait,1)):animalDom + idh(at.level(trait,2)):animalDom + 
                       idh(at.level(trait,1)):matID + idh(at.level(trait,2)):matID, 
                     ginverse = list(animal=Ainv, animalDom =  Dinv), 
                     rcov = ~idh(trait):units,
                     family = "zipoisson", data = ambient.total, prior=prior1.0, 
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr=TRUE)

load(file="Routput/Ambient/AM_model1.0.RData")
plot(AM_model1.0$VCV) # good trace plots
summary(AM_model1.0) # good effective sample size
autocorr.diag(AM_model1.0$VCV) # no autocorrelation
heidel.diag(AM_model1.0$VCV) # convergence success

#Heritability - Zero inflation process
AM_herit1.0_ZI <- AM_model1.0$VCV[, "at.level(trait, 2).animal"]/(AM_model1.0$VCV[, "at.level(trait, 2).animal"] + AM_model1.0$VCV[, "at.level(trait, 2).animalDom"] + AM_model1.0$VCV[, "at.level(trait, 2).matID"] + AM_model1.0$VCV[, "traitzi_seed_pods.units"])
mean(AM_herit1.0_ZI)
posterior.mode(AM_herit1.0_ZI)
HPDinterval(AM_herit1.0_ZI)

#Heritability - Poisson Process
AM_herit1.0_PP <- AM_model1.0$VCV[, "at.level(trait, 1).animal"]/(AM_model1.0$VCV[, "at.level(trait, 1).animal"] + AM_model1.0$VCV[, "at.level(trait, 1).animalDom"] + AM_model1.0$VCV[, "at.level(trait, 1).matID"] + AM_model1.0$VCV[, "traitseed_pods.units"])
mean(AM_herit1.0_PP)
posterior.mode(AM_herit1.0_PP)
HPDinterval(AM_herit1.0_PP)

#Transforming from latent to data scale using package "QGglmm"
#Zero-inflation process
yhat_A1.0_ZI <- predict(AM_model1.0, type = "terms")
vp_A1.0_ZI <- mean(AM_model1.0$VCV[, "at.level(trait, 2).animal"] + AM_model1.0$VCV[, "at.level(trait, 2).animalDom"] + AM_model1.0$VCV[, "at.level(trait, 2).matID"] + AM_model1.0$VCV[, "traitzi_seed_pods.units"]) 
va_A1.0_ZI <- mean(AM_model1.0[["VCV"]][ , "at.level(trait, 2).animal"]) #ZERO INFLATION
vd_A1.0_ZI <- mean(AM_model1.0[["VCV"]][ , "at.level(trait, 2).animalDom"]) #ZERO INFLATION
vm_A1.0_ZI <- mean(AM_model1.0[["VCV"]][ , "at.level(trait, 2).matID"]) #ZERO INFLATION

#Additive Variance
QGA_A1.0_ZI <- QGparams(predict=yhat_A1.0_ZI, var.a=va_A1.0_ZI, var.p=vp_A1.0_ZI, model = "binom1.probit")
QGA_A1.0_ZI #Va = 0.003772693 .. 2.09% of Vp is Va .. Obs Mean = 0.764

#Dominance Variance
QGD_A1.0_ZI <- QGicc(predict=yhat_A1.0_ZI, var.comp=vd_A1.0_ZI, var.p=vp_A1.0_ZI, model = "binom1.probit")
QGD_A1.0_ZI #Vd = 0.006132166 .. 3.40% of Vp is Vd

#Maternal Variance
QGM_A1.0_ZI <- QGicc(predict=yhat_A1.0_ZI, var.comp=vm_A1.0_ZI, var.p=vp_A1.0_ZI, model = "binom1.probit")
QGM_A1.0_ZI #Vm = 0.002261614 .. 1.25% of Vp is Vm



#Poisson Process
yhat_A1.0_PP <- predict(AM_model1.0, type = "terms")
vp_A1.0_PP <- mean(AM_model1.0$VCV[, "at.level(trait, 1).animal"] + AM_model1.0$VCV[, "at.level(trait, 1).animalDom"] + AM_model1.0$VCV[, "at.level(trait, 1).matID"] + AM_model1.0$VCV[, "traitseed_pods.units"])
va_A1.0_PP <- mean(AM_model1.0[["VCV"]][ , "at.level(trait, 1).animal"]) #POISSON PROCESS
vd_A1.0_PP <- mean(AM_model1.0[["VCV"]][ , "at.level(trait, 1).animalDom"]) #POISSON PROCESS
vm_A1.0_PP <- mean(AM_model1.0[["VCV"]][ , "at.level(trait, 1).matID"]) #POISSON PROCESS

#Additive Variance
QGA_A1.0_PP <- QGparams(predict=yhat_A1.0_PP, var.a=va_A1.0_PP, var.p=vp_A1.0_PP, model = "Poisson.log")
QGA_A1.0_PP #Va = 5.197834 .. 2.88% of Vp is Va .. Obs Mean = 7.63

#Dominance Variance
QGD_A1.0_PP <- QGicc(predict=yhat_A1.0_PP, var.comp=vd_A1.0_PP, var.p=vp_A1.0_PP, model = "Poisson.log")
QGD_A1.0_PP #Vd = 8.63635 .. 4.79% of Vp is Vd

#Maternal Variance
QGM_A1.0_PP <- QGicc(predict=yhat_A1.0_PP, var.comp=vm_A1.0_PP, var.p=vp_A1.0_PP, model = "Poisson.log")
QGM_A1.0_PP #Vm = 0.9974059 .. 0.553% of Vp is Vm


#Model 1.1
#Note really sure what's different in Model 1.1 tbh but I'm keeping it here for now. 
prior1.1 <- list(R = list(V = diag(2), nu = 1, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1 - Zero Inflation
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1), # Expanded Fisher prior for Poisson process
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1 - Zero Inflation
                          G5 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1), # Expanded Fisher prior for Poisson process
                          G6 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) # Chi Square Prior with SD 1 - Zero Inflation

AM_model1.1 <- MCMCglmm(seed_pods ~ trait - 1 + plot,
                        random = ~idh(at.level(trait,1)):animal + idh(at.level(trait,2)):animal +
                          idh(at.level(trait,1)):animalDom + idh(at.level(trait,2)):animalDom + 
                          idh(at.level(trait,1)):matID + idh(at.level(trait,2)):matID, 
                        ginverse = list(animal=Ainv, animalDom =  Dinv), 
                        rcov = ~idh(trait):units,
                        family = "zipoisson", data = ambient.total, prior=prior1.1, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr=TRUE)

load(file="Routput/Ambient/AM_model1.1.RData")
plot(AM_model1.1$VCV) # 
summary(AM_model1.1) #
autocorr.diag(AM_model1.1$VCV) #
heidel.diag(AM_model1.1$VCV) #

#Model 1.2 - No dominance, Poisson family
load(file="Routput/Ambient/AM_model1.2.RData")
plot(AM_model1.2$VCV) # 
summary(AM_model1.2) #
autocorr.diag(AM_model1.2$VCV) #
heidel.diag(AM_model1.2$VCV) #

#Model 1.3 - No dominance, Zero inflated Poisson family
load(file="Routput/Ambient/AM_model1.3.RData")
plot(AM_model1.3$VCV) # 
summary(AM_model1.3) #
autocorr.diag(AM_model1.3$VCV) #
heidel.diag(AM_model1.3$VCV) #








#### No. 2 Lifetime Fitness V2: Survival x Fecundity | Bivariate ####
#Note: For this trait, we include all plants from the experiment
#Survival to Flowering = Bernoulli Dist | Fecundity = Poisson Dist

#Calculating Variance for Fecundity for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  summarise(mean(seed_pods), sd=sd(seed_pods), var=(sd)^2, n=n(), se=sd/sqrt(n)) #mean = 4.985826
log(4.985826) # = 1.606599

MCMC.ambient %>%
  summarise(sum.f=sum(flower), n=n(), mean=sum.f/n) # mean prob = 0.41105
logit(0.41105*(1-0.41105)) # -1.141267

#Model 2.0 with additive, dominance, and maternal random effects with plot fixed effects
#Using prior similar to plasticity models - See prior visualization tool in 09
#Note: we set rcov = ~units because covariance can't be estimated in residual variance for a binary trait, AND we set residual variance to fix = 2 (the number of multivariate dimensions)
ambient.Survival <- MCMC.ambient
prior2.0 <- list(R = list(V = diag(2)*0.02, nu = 3, fix = 2),
                 G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000))),
                          G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000))),
                          G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000)))
                 ))
prior2.0$R$V[2,2] <- 10 # fixing last element of residual variance to 10 to improve mixing, following Bemmels & Anderson

AM_model2.0 <- MCMCglmm(cbind(seed_pods, flower) ~ trait + plot - 1,
                        random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv), rcov = ~us(trait):units,
                        family = c("poisson", "threshold"), 
                        data = ambient.Survival, prior = prior2.0, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Ambient/AM_model2.0.RData")
plot(AM_model2.0$VCV) # poor trace plots
summary(AM_model2.0) # poor effective sample sizes
autocorr.diag(AM_model2.0$VCV) # high autocorrelation
heidel.diag(AM_model2.0$VCV) #convergence fail

#Model 2.1 without dominance variance
ambient.Survival <- MCMC.ambient
prior2.1 <- list(R = list(V = diag(2)*0.02, nu = 3, fix = 2),
                 G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000))),
                          G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000))),
                          G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1000,1000)))
                 ))
prior2.1$R$V[2,2] <- 10 # fixing last element of residual variance to 10 to improve mixing, following Bemmels & Anderson

AM_model2.1 <- MCMCglmm(cbind(seed_pods, flower) ~ trait + plot - 1,
                        random = ~us(trait):animal + idh(trait):animalDom + us(trait):matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv), rcov = ~us(trait):units,
                        family = c("poisson", "threshold"), 
                        data = ambient.Survival, prior = prior2.1, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Ambient/AM_model2.1.RData")
plot(AM_model2.1$VCV) # poor trace plots
summary(AM_model2.1) # poor sample size
autocorr.diag(AM_model2.1$VCV) # bad autocorrelation
heidel.diag(AM_model2.1$VCV) # convergence fail

#Calculating Genetic Correlations = Covariance b/w trait 1 & 2 / sqrt (var [trait 1] * var [trait 2])

# (1) Survival & Fecundity
gen.corrA2.1 <-AM_model2.1$VCV[,'traitflower:traitseed_pods.animal']/
  sqrt(AM_model2.1$VCV[,'traitflower:traitflower.animal']*AM_model2.1$VCV[,'traitseed_pods:traitseed_pods.animal']) 
mean(gen.corrA2.1) #
posterior.mode(gen.corrA2.1) #
HPDinterval(gen.corrA2.1) # 

phen.corrA2.1 <- cor.test(ambient.Survival$flower, ambient.Survival$seed_pods, method="spearman", exact=FALSE)
phen.corrA2.1 # Corr = 0.5112517 | Posterior 95% CI =  (0.4861976, 0.5354659)

#Model 2.2 - No dominance
load(file="Routput/Ambient/AM_model2.2.RData")
plot(AM_model2.2$VCV) # 
summary(AM_model2.2) #
autocorr.diag(AM_model2.2$VCV) #
heidel.diag(AM_model2.2$VCV) #

#Model 2.3 - No dominance and maternal effects
load(file="Routput/Ambient/AM_model2.3.RData")
plot(AM_model2.3$VCV) # 
summary(AM_model2.3) #
autocorr.diag(AM_model2.3$VCV) #
heidel.diag(AM_model2.3$VCV) #

#Model 2.4 - Attempting full model again with new prior
ambient.Survival <- MCMC.ambient
prior2.4 <- list(R = list(V = diag(2), nu = 2, fix = 2),
                 G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G3 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))
                 ))
set.seed(1)
AM_model2.4 <- mclapply(1:20, function(i) {
	MCMCglmm(cbind(seed_pods, flower) ~ trait + plot - 1,
                        random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv), rcov = ~us(trait):units, 
                        family = c("poisson", "threshold"), 
                        data = ambient.Survival, prior = prior2.4, 
                        nitt = 200000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)
}, mc.cores = 20)

load(file="Routput/Ambient/AM_model2.4.RData")
AM_2.4_sol <- lapply(AM_model2.4, function(m) m$Sol)
AM_2.4_sol <- do.call(mcmc.list, AM_2.4_sol)

AM_2.4_vcv <- lapply(AM_model2.4, function(m) m$VCV)
AM_2.4_vcv <- do.call(mcmc.list, AM_2.4_vcv)
AM_2.4_vcv <- mcmc.stack(AM_2.4_vcv)
plot(AM_2.4_vcv) #poor trace plots
autocorr.diag(AM_2.4_vcv, relative=FALSE) #high autocorrelation
heidel.diag(AM_2.4_vcv) #convergence fail
effectiveSize(AM_2.4vcv) #low effective size

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

# UPDATE DEC 2024: We do not report these in the revised manuscript. 





#### No. 3 Survival to Flowering - Univariate | Bernoulli ####
#Note: For this trait, we include all plants from the experiment

#Probability of Survival to Flowering
ambient.Survival <- MCMC.ambient
ambient.Survival %>% 
  summarise(tot.flwr=sum(flower), n=n(), prob=tot.flwr/n) #prob = 0.411169
log(0.411169*(1-0.411169)) # -1.418367

#Model 3.0 with additive + dominance + maternal random effects and plot fixed effects
#We use a prior suggested by de Villemereuil et al. 2013 for binary traits
ambient.Survival <- MCMC.ambient
prior3.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model3.0 <- MCMCglmm(flower ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal=Ainv, animalDom = Dinv),
                     family = "threshold", data = ambient.Survival, prior = prior3.0, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)

load(file="Routput/Ambient/AM_model3.0.RData")
plot(AM_model3.0$VCV) # good trace plots
summary(AM_model3.0) # good effective sample size
autocorr.diag(AM_model3.0$VCV) # no autocorrelation
heidel.diag(AM_model3.0$VCV) # convergence success

#Latent Scale 
mean(AM_model3.0[["VCV"]][ , "animal"]) #Va
mean(AM_model3.0[["VCV"]][ , "animalDom"]) #Vd
mean(AM_model3.0[["VCV"]][ , "matID"]) #Vm
mean(AM_model3.0[["VCV"]][ , "units"]) #Vr
AM_herit3.0 <- AM_model3.0$VCV[, "animal"]/(AM_model3.0$VCV[, "animal"] + AM_model3.0$VCV[, "animalDom"] + AM_model3.0$VCV[, "matID"] + AM_model3.0$VCV[, "units"]) # or since units is fixed, can just add + 1
mean(AM_herit3.0)
posterior.mode(AM_herit3.0)
HPDinterval(AM_herit3.0)

mean(AM_model3.0$VCV[, "animal"] + AM_model3.0$VCV[, "animalDom"] + AM_model3.0$VCV[, "matID"]) #Vg
HPDinterval(AM_model3.0$VCV[, "animal"] + AM_model3.0$VCV[, "animalDom"] + AM_model3.0$VCV[, "matID"])


#Transforming from latent to data scale using package "QGglmm"
df_A3.0 <- data.frame(
  va = as.vector(AM_model3.0[["VCV"]][, "animal"]),
  vd = as.vector(AM_model3.0[["VCV"]][, "animalDom"]),
  vm = as.vector(AM_model3.0[["VCV"]][, "matID"]),
  vp = rowSums(AM_model3.0[["VCV"]]))
yhat_A3.0 <- predict(AM_model3.0, type = "terms")

#Va
Va_A3.0 <- do.call("rbind", apply(df_A3.0, 1, function(row){
  QGparams(predict = yhat_A3.0, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Va_A3.0$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_A3.0$var.a.obs))

mean(as.mcmc(Va_A3.0$h2.obs)) #h2
HPDinterval(as.mcmc(Va_A3.0$h2.obs))

mu_A3.0 <- mean(as.mcmc(Va_A3.0$mean.obs)) #mean 
mu_A3.0
HPDinterval(as.mcmc(Va_A3.0$mean.obs))
            
mean(as.mcmc(Va_A3.0$h2.obs))/mu_A3.0 #Rw
HPDinterval(as.mcmc(Va_A3.0$h2.obs))/mu_A3.0

#Vd
Vd_A3.0 <- do.call("rbind", apply(df_A3.0, 1, function(row){
  QGicc(predict = yhat_A3.0, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vd_A3.0$var.comp.obs))
HPDinterval(as.mcmc(Vd_A3.0$var.comp.obs))

#Vm
Vm_A3.0 <- do.call("rbind", apply(df_A3.0, 1, function(row){
  QGicc(predict = yhat_A3.0, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_A3.0$var.comp.obs))
HPDinterval(as.mcmc(Vm_A3.0$var.comp.obs))

#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_A3.0 <- data.frame(
  Va = Va_A3.0$var.a.obs,
  Vd = Vd_A3.0$var.comp.obs,
  Vm = Vm_A3.0$var.comp.obs,
  mu = Va_A3.0$mean.obs,
  h2 = Va_A3.0$h2.obs,
  Rw = (Va_A3.0$var.a.obs / Va_A3.0$mean.obs)) %>%
  pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_A3.0 <- data.frame(
  Component = c("Va", "Vd", "Vm", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(Va_A3.0$var.a.obs)),
    mean(as.mcmc(Vd_A3.0$var.comp.obs)),
    mean(as.mcmc(Vm_A3.0$var.comp.obs)),
    mean(as.mcmc(Va_A3.0$mean.obs)),
    mean(as.mcmc(Va_A3.0$h2.obs)),
    mean(as.mcmc(Va_A3.0$var.a.obs / Va_A3.0$mean.obs))
  ),
  Lower = c(
    HPDinterval(as.mcmc(Va_A3.0$var.a.obs))[1],
    HPDinterval(as.mcmc(Vd_A3.0$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vm_A3.0$var.comp.obs))[1],
    HPDinterval(as.mcmc(Va_A3.0$mean.obs))[1],
    HPDinterval(as.mcmc(Va_A3.0$h2.obs))[1],
    HPDinterval(as.mcmc(Va_A3.0$var.a.obs / Va_A3.0$mean.obs))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(Va_A3.0$var.a.obs))[2],
    HPDinterval(as.mcmc(Vd_A3.0$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vm_A3.0$var.comp.obs))[2],
    HPDinterval(as.mcmc(Va_A3.0$mean.obs))[2],
    HPDinterval(as.mcmc(Va_A3.0$h2.obs))[2],
    HPDinterval(as.mcmc(Va_A3.0$var.a.obs / Va_A3.0$mean.obs))[2]
  )
)

posteriors_A3.0$Component <- factor(posteriors_A3.0$Component, levels = c("Vm", "Va", "Vd", "mu", "Rw", "h2"))
save(posteriors_A3.0, stats_A3.0, file="Routput/Ambient/AM_stats3.0.RData")
load("Routput/Ambient/AM_stats3.0.RData")

var_A3.0_fig <- posteriors_A3.0 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
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
  theme(axis.title.x = ggtext::element_markdown())
var_A3.0_fig

ggsave(plot = var_A3.0_fig, filename = "Routput/Figures/fig_3_A3.0_var.png", width = 8, height = 8, units="cm")

var_A3.0_fig2 <- posteriors_A3.0 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2, adjust = 1.5) +
  scale_fill_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  scale_color_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  labs(
    x = expression("Proportion of variance (" * h^2 * " or " * Δ[evol] * bar(W) * ")"),
    y = "Density") +
  xlim(0, 0.04) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none")
var_A3.0_fig2

ggsave(plot = var_A3.0_fig2, filename = "Routput/Figures/fig_3_A3.0_h2.png", width = 8, height = 8, units="cm")


#Plotting density plot of breeding values
AM_3.0_breed <- AM_model3.0$Sol[,7:7657]
AM_3.0_BV <- posterior.mode(as.mcmc(AM_3.0_breed)) #breeding values for each individual
AM_3.0_BV_CI <- HPDinterval(as.mcmc(AM_3.0_breed)) #confidence intervals of breeding values
AM_3.0_BV_plot <- as.data.frame(AM_3.0_BV)
AM_3.0_BV_plot %>% 
  ggplot(aes(x=AM_3.0_BV)) +
  geom_density(fill="dodgerblue2", color="dodgerblue2", alpha=0.3) +
  labs(x="Ambient Survival Breeding Values", y="Frequency")  + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

#Model 3.1 - No dominance
load(file="Routput/Ambient/AM_model3.1.RData")
plot(AM_model3.1$VCV) # trace plots
summary(AM_model3.1) # good effective sample size
autocorr.diag(AM_model3.1$VCV) # no autocorrelation
heidel.diag(AM_model3.1$VCV) # convergence success


mean(AM_model3.1[["VCV"]][ , "animal"]) #Va
mean(AM_model3.1[["VCV"]][ , "matID"]) #Vm
mean(AM_model3.1[["VCV"]][ , "units"]) #Vr


#Transforming from latent to data scale using package "QGglmm"
df_A3.1 <- data.frame(
  va = as.vector(AM_model3.1[["VCV"]][, "animal"]),
  vm = as.vector(AM_model3.1[["VCV"]][, "matID"]),
  vr = as.vector(AM_model3.1[["VCV"]][, "units"]),
  vp = rowSums(AM_model3.1[["VCV"]]))
yhat_A3.1 <- predict(AM_model3.1, type = "terms")

#Va
Va_A3.1 <- do.call("rbind", apply(df_A3.1, 1, function(row){
  QGparams(predict = yhat_A3.1, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Va_A3.1$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_A3.1$var.a.obs))

#Vm
Vm_A3.1 <- do.call("rbind", apply(df_A3.1, 1, function(row){
  QGicc(predict = yhat_A3.1, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_A3.1$var.comp.obs))
HPDinterval(as.mcmc(Vm_A3.1$var.comp.obs))

#Vr
Vm_A3.1 <- do.call("rbind", apply(df_A3.1, 1, function(row){
  QGicc(predict = yhat_A3.1, 
        var.comp = row[["vr"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_A3.1$var.comp.obs))
HPDinterval(as.mcmc(Vm_A3.1$var.comp.obs))














#### No. 4 Fecundity of Flowering Plants - Univariate | Poisson ####
#Note: For Fecundity, we only include plants that survived to reach flowering

#Calculating Variance for Fecundity for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1 & flower==1) %>%
  summarise(mean=mean(seed_pods), sd=sd(seed_pods), var=(sd)^2,n=n(), se=sd/sqrt(n)) #var = 244.9818 *******
log(12.12949)

#Model 4.0 with additive + dominance + maternal random effects and plot fixed effects.
#Prior with inverse gamma prior (0.001, 0.001)
ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior4.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model4.0 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "poisson", data = ambient.Fecundity, prior = prior4.0,
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model4.0.RData")
plot(AM_model4.0$VCV) # poor trace plots
summary(AM_model4.0) # low effective sample size
autocorr.diag(AM_model4.0$VCV) # autocorrelation
heidel.diag(AM_model4.0$VCV) # convergence failing

#Model 4.1 with Fisher Parameter Expanded Prior
ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior4.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model4.1 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal=Ainv, animalDom =  Dinv),
                        family = "poisson", data = ambient.Fecundity, prior = prior4.1,
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model4.1.RData")
plot(AM_model4.1$VCV) # good trace plots
summary(AM_model4.1) # good effective sample size
autocorr.diag(AM_model4.1$VCV) # no autocorrelation
heidel.diag(AM_model4.1$VCV) # convergence success

#Latent Scale
mean(AM_model4.1[["VCV"]][ , "animal"]) #Va
mean(AM_model4.1[["VCV"]][ , "animalDom"]) #Vd
mean(AM_model4.1[["VCV"]][ , "matID"]) #Vm
mean(AM_model4.1[["VCV"]][ , "units"]) #Vr
AM_herit4.1 <- AM_model4.1$VCV[, "animal"]/(AM_model4.1$VCV[, "animal"] + AM_model4.1$VCV[, "animalDom"] + AM_model4.1$VCV[, "matID"] + AM_model4.1$VCV[, "units"]) #
mean(AM_herit4.1)
posterior.mode(AM_herit4.1) #Heritability
HPDinterval(AM_herit4.1) #h2 HPD

posterior.mode(AM_model4.1$VCV[, "animal"] + AM_model4.1$VCV[, "animalDom"] + AM_model4.1$VCV[, "matID"]) #Vg
HPDinterval(AM_model4.1$VCV[, "animal"] + AM_model4.1$VCV[, "animalDom"] + AM_model4.1$VCV[, "matID"])


#Transforming from latent to data scale using package "QGglmm"
df_A4.1 <- data.frame(
  va = as.vector(AM_model4.1[["VCV"]][, "animal"]),
  vd = as.vector(AM_model4.1[["VCV"]][, "animalDom"]),
  vm = as.vector(AM_model4.1[["VCV"]][, "matID"]),
  vr = as.vector(AM_model4.1[["VCV"]][, "units"]),
  vp = rowSums(AM_model4.1[["VCV"]]))
yhat_A4.1 <- predict(AM_model4.1, type = "terms")

#Va
Va_A4.1 <- do.call("rbind", apply(df_A4.1, 1, function(row){
  QGparams(predict = yhat_A4.1, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Va_A4.1$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_A4.1$var.a.obs))

mean(as.mcmc(Va_A4.1$h2.obs)) #h2
HPDinterval(as.mcmc(Va_A4.1$h2.obs))

mu_A4.1 <- mean(as.mcmc(Va_A4.1$mean.obs))
mu_A4.1
HPDinterval(as.mcmc(Va_A4.1$mean.obs))

mean(as.mcmc(Va_A4.1$var.a.obs))/mu_A4.1 #Rw
(HPDinterval(as.mcmc(Va_A4.1$var.a.obs)))/mu_A4.1

#Vd
Vd_A4.1 <- do.call("rbind", apply(df_A4.1, 1, function(row){
  QGicc(predict = yhat_A4.1, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vd_A4.1$var.comp.obs))
HPDinterval(as.mcmc(Vd_A4.1$var.comp.obs))

#Vm
Vm_A4.1 <- do.call("rbind", apply(df_A4.1, 1, function(row){
  QGicc(predict = yhat_A4.1, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vm_A4.1$var.comp.obs))
HPDinterval(as.mcmc(Vm_A4.1$var.comp.obs))

#Vr
Vr_A4.1 <- do.call("rbind", apply(df_A4.1, 1, function(row){
  QGicc(predict = yhat_A4.1, 
        var.comp = row[["vr"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vr_A4.1$var.comp.obs))
HPDinterval(as.mcmc(Vr_A4.1$var.comp.obs))

#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_A4.1 <- data.frame(
  Va = Va_A4.1$var.a.obs,
  Vd = Vd_A4.1$var.comp.obs,
  Vm = Vm_A4.1$var.comp.obs,
  Vr = Vr_A4.1$var.comp.obs,
  mu = Va_A4.1$mean.obs,
  h2 = Va_A4.1$h2.obs,
  Rw = (Va_A4.1$var.a.obs / Va_A4.1$mean.obs)) %>%
  pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_A4.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(Va_A4.1$var.a.obs)),
    mean(as.mcmc(Vd_A4.1$var.comp.obs)),
    mean(as.mcmc(Vm_A4.1$var.comp.obs)),
    mean(as.mcmc(Vr_A4.1$var.comp.obs)),
    mean(as.mcmc(Va_A4.1$mean.obs)),
    mean(as.mcmc(Va_A4.1$h2.obs)),
    mean(as.mcmc(Va_A4.1$var.a.obs / Va_A4.1$mean.obs))
  ),
  Lower = c(
    HPDinterval(as.mcmc(Va_A4.1$var.a.obs))[1],
    HPDinterval(as.mcmc(Vd_A4.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vm_A4.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vr_A4.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Va_A4.1$mean.obs))[1],
    HPDinterval(as.mcmc(Va_A4.1$h2.obs))[1],
    HPDinterval(as.mcmc(Va_A4.1$var.a.obs / Va_A4.1$mean.obs))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(Va_A4.1$var.a.obs))[2],
    HPDinterval(as.mcmc(Vd_A4.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vm_A4.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vr_A4.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Va_A4.1$mean.obs))[2],
    HPDinterval(as.mcmc(Va_A4.1$h2.obs))[2],
    HPDinterval(as.mcmc(Va_A4.1$var.a.obs / Va_A4.1$mean.obs))[2]
  )
)

posteriors_A4.1$Component <- factor(posteriors_A4.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "h2", "Rw"))
save(posteriors_A4.1, stats_A4.1, file="Routput/Ambient/AM_stats4.1.RData")
load("Routput/Ambient/AM_stats4.1.RData")

var_A4.1_fig <- posteriors_A4.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, 75) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown())
var_A4.1_fig

ggsave(plot = var_A4.1_fig, filename = "Routput/Figures/fig_3_A4.1_var.png", width = 8, height = 8, units="cm")

var_A4.1_fig2 <- posteriors_A4.1 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2, adjust = 1.5) +
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
        legend.position  = "none")
var_A4.1_fig2

ggsave(plot = var_A4.1_fig2, filename = "Routput/Figures/fig_3_A4.1_h2.png", width = 8, height = 8, units="cm")


#Plotting density plot of breeding values
AM_4.1_breed <- AM_model4.1$Sol[,7:7657]
AM_4.1_BV <- posterior.mode(as.mcmc(AM_4.1_breed)) #breeding values for each individual
AM_4.1_BV_CI <- HPDinterval(as.mcmc(AM_4.1_breed)) #confidence intervals of breeding values
AM_4.1_BV_plot <- as.data.frame(AM_4.1_BV)
AM_4.1_BV_plot %>% 
  ggplot(aes(x=AM_4.1_BV)) +
  geom_density(fill="dodgerblue2", color="dodgerblue2", alpha=0.3) +
  labs(x="Ambient Fecundity Breeding Values", y="Frequency")  + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

#Model 4.2 with Fisher Parameter Expanded Prior - prior sensitivity check
ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior4.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

AM_model4.2 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "poisson", data = ambient.Fecundity, prior = prior4.2,
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model4.2.RData")
plot(AM_model4.2$VCV) # good trace plots
summary(AM_model4.2) # good effective sample size
autocorr.diag(AM_model4.2$VCV) # no autocorrelation
heidel.diag(AM_model4.2$VCV) # convergence success
#No difference. Insensitive to prior change

#Model 4.3 - No dominance
load(file="Routput/Ambient/AM_model4.3.RData")
plot(AM_model4.3$VCV) # good trace plots
summary(AM_model4.3) # good effective sample size
autocorr.diag(AM_model4.3$VCV) # no autocorrelation
heidel.diag(AM_model4.3$VCV) # convergence success

mean(AM_model4.3[["VCV"]][ , "animal"]) #Va
mean(AM_model4.3[["VCV"]][ , "matID"]) #Vm
mean(AM_model4.3[["VCV"]][ , "units"]) #Vr
AM_herit4.3 <- AM_model4.3$VCV[, "animal"]/(AM_model4.3$VCV[, "animal"] + AM_model4.3$VCV[, "matID"] + AM_model4.3$VCV[, "units"]) #
mean(AM_herit4.3)
posterior.mode(AM_herit4.3) #Heritability
HPDinterval(AM_herit4.3) #h2 HPD


#Transforming from latent to data scale using package "QGglmm"
df_A4.3 <- data.frame(
  va = as.vector(AM_model4.3[["VCV"]][, "animal"]),
  vm = as.vector(AM_model4.3[["VCV"]][, "matID"]),
  vr = as.vector(AM_model4.3[["VCV"]][, "units"]),
  vp = rowSums(AM_model4.3[["VCV"]]))
yhat_A4.3 <- predict(AM_model4.3, type = "terms")

#Va
Va_A4.3 <- do.call("rbind", apply(df_A4.3, 1, function(row){
  QGparams(predict = yhat_A4.3, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Va_A4.3$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_A4.3$var.a.obs))

#Vm
Vm_A4.3 <- do.call("rbind", apply(df_A4.3, 1, function(row){
  QGicc(predict = yhat_A4.3, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vm_A4.3$var.comp.obs))
HPDinterval(as.mcmc(Vm_A4.3$var.comp.obs))

#Vr
Vr_A4.3 <- do.call("rbind", apply(df_A4.3, 1, function(row){
  QGicc(predict = yhat_A4.3, 
        var.comp = row[["vr"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vr_A4.3$var.comp.obs))
HPDinterval(as.mcmc(Vr_A4.3$var.comp.obs))

#Model 4.4 - Full with different smaller V prior
ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior4.4 <- list(R = list(V = 0.01, nu = 0.002),
                 G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

AM_model4.4 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "poisson", data = ambient.Fecundity, prior = prior4.4,
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model4.4.RData")
plot(AM_model4.4$VCV) # bad trace plots
summary(AM_model4.4) # low effective sample size
autocorr.diag(AM_model4.4$VCV) # high autocorrelation
heidel.diag(AM_model4.4$VCV) # convergence fail







#### No. 5 Overwintering Survival - Univariate | Bernoulli ####
#Note: For this trait, we include all plants from the study

#Model 5.0 with additive, dominance and maternal random effects + plot fixed effects 
#We use a prior suggested by de Villemereuil et al. 2013 for binary traits
prior5.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model5.0 <- MCMCglmm(germ ~ plot, random = ~animal + animalDom + matID, 
                        ginverse = list(animal=Ainv, animalDom = Dinv),
                        family = "threshold", data = ambient.Germination, prior = prior5.0, #Bernoulli distribution
                        nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Ambient/AM_model5.0.RData")
plot(AM_model5.0$VCV) # trace plots look OK
summary(AM_model5.0) # good effective sample size
autocorr.diag(AM_model5.0$VCV) # no autocorrelation
heidel.diag(AM_model5.0$VCV) # convergence success

#Latent Scale
mean(AM_model5.0[["VCV"]][ , "animal"]) #Va
mean(AM_model5.0[["VCV"]][ , "animalDom"]) #Vd
mean(AM_model5.0[["VCV"]][ , "matID"]) #Vm
AM_herit5.0 <- AM_model5.0$VCV[, "animal"]/(AM_model5.0$VCV[, "animal"] + AM_model5.0$VCV[, "animalDom"] + AM_model5.0$VCV[, "matID"] + 1)
mean(AM_herit5.0)
posterior.mode(AM_herit5.0)
HPDinterval(AM_herit5.0)

#Transforming from latent to data scale using package "QGglmm"
df_A5.0 <- data.frame(
  va = as.vector(AM_model5.0[["VCV"]][, "animal"]),
  vd = as.vector(AM_model5.0[["VCV"]][, "animalDom"]),
  vm = as.vector(AM_model5.0[["VCV"]][, "matID"]),
  vp = rowSums(AM_model5.0[["VCV"]]))
yhat_A5.0 <- predict(AM_model5.0, type = "terms")

#Va
Va_A5.0 <- do.call("rbind", apply(df_A5.0, 1, function(row){
  QGparams(predict = yhat_A5.0, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Va_A5.0$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_A5.0$var.a.obs))

mean(as.mcmc(Va_A5.0$h2.obs)) #h2
HPDinterval(as.mcmc(Va_A5.0$h2.obs))

mu_A5.0 <- mean(as.mcmc(Va_A5.0$mean.obs))
mean(as.mcmc(Va_A5.0$h2.obs))/mu_A5.0 #Rw
HPDinterval(as.mcmc(Va_A5.0$h2.obs))/mu_A5.0

#Vd
Vd_A5.0 <- do.call("rbind", apply(df_A5.0, 1, function(row){
  QGicc(predict = yhat_A5.0, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vd_A5.0$var.comp.obs))
HPDinterval(as.mcmc(Vd_A5.0$var.comp.obs))

#Vm
Vm_A5.0 <- do.call("rbind", apply(df_A5.0, 1, function(row){
  QGicc(predict = yhat_A5.0, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_A5.0$var.comp.obs))
HPDinterval(as.mcmc(Vm_A5.0$var.comp.obs))

#Model 5.1 - No dominance
load(file="Routput/Ambient/AM_model5.1.RData")
plot(AM_model5.1$VCV) # good trace plots 
summary(AM_model5.1) # good effective sample size
autocorr.diag(AM_model5.1$VCV) # no autocorrelation
heidel.diag(AM_model5.1$VCV) # convergence success













#### No. 6 Spring-Summer Survival Success - Univariate | Bernoulli ####
#Note: For this trait, we include only include plants that germinated

#Model 6.0 with additive + dominance + maternal random effects and plot fixed effects
#We use a prior suggested by de Villemereuil et al. 2013 for binary traits
ambient.Flowering <- MCMC.ambient %>% filter(germ==1)
prior6.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model6.0 <- MCMCglmm(flower ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal=Ainv, animalDom = Dinv),
                     family = "threshold", data = ambient.Flowering, prior = prior6.0, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Ambient/AM_model6.0.RData")
plot(AM_model6.0$VCV) # good trace plots
summary(AM_model6.0) # good effective sample size
autocorr.diag(AM_model6.0$VCV) # no autocorrelation
heidel.diag(AM_model6.0$VCV) # convergence success

#Heritability
mean(AM_model6.0[["VCV"]][ , "animal"])
mean(AM_model6.0[["VCV"]][ , "animalDom"])
mean(AM_model6.0[["VCV"]][ , "matID"])
AM_herit6.0 <- AM_model6.0$VCV[, "animal"]/(AM_model6.0$VCV[, "animal"] + AM_model6.0$VCV[, "animalDom"] + AM_model6.0$VCV[, "matID"] + 1)
mean(AM_herit6.0)
posterior.mode(AM_herit6.0)
HPDinterval(AM_herit6.0)

#Transforming from latent to data scale using package "QGglmm"
df_A6.0 <- data.frame(
  va = as.vector(AM_model6.0[["VCV"]][, "animal"]),
  vd = as.vector(AM_model6.0[["VCV"]][, "animalDom"]),
  vm = as.vector(AM_model6.0[["VCV"]][, "matID"]),
  vp = rowSums(AM_model6.0[["VCV"]]))
yhat_A6.0 <- predict(AM_model6.0, type = "terms")

#Va
Va_A6.0 <- do.call("rbind", apply(df_A6.0, 1, function(row){
  QGparams(predict = yhat_A6.0, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Va_A6.0$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_A6.0$var.a.obs))

mean(as.mcmc(Va_A6.0$h2.obs)) #h2
HPDinterval(as.mcmc(Va_A6.0$h2.obs))

mu_A6.0 <- mean(as.mcmc(Va_A6.0$mean.obs))
mean(as.mcmc(Va_A6.0$h2.obs))/mu_A6.0 #Rw
HPDinterval(as.mcmc(Va_A6.0$h2.obs))/mu_A6.0

#Vd
Vd_A6.0 <- do.call("rbind", apply(df_A6.0, 1, function(row){
  QGicc(predict = yhat_A6.0, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vd_A6.0$var.comp.obs))
HPDinterval(as.mcmc(Vd_A6.0$var.comp.obs))

#Vm
Vm_A6.0 <- do.call("rbind", apply(df_A6.0, 1, function(row){
  QGicc(predict = yhat_A6.0, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_A6.0$var.comp.obs))
HPDinterval(as.mcmc(Vm_A6.0$var.comp.obs))

#Model 6.1 - No dominance
load(file="Routput/Ambient/AM_model6.1.RData")
plot(AM_model6.1$VCV) # good trace plots
summary(AM_model6.1) # good effective sample size
autocorr.diag(AM_model6.1$VCV) # no autocorrelation
heidel.diag(AM_model6.1$VCV) # convergence success













#### No. 7 Fruiting Success - Univariate | Bernoulli ####
#Note: For this trait, we include only include plants that germinated and flowered

#Producing filtered data set for Fruiting Success
ambient.SeedMaturation <- MCMC.ambient %>% filter(germ==1, flower==1)

#Model 7.0 with additive, dominance and maternal random effects + plot fixed effects
#We use a prior suggested by de Villemereuil et al. 2013 for binary traits
ambient.SeedMaturation <- MCMC.ambient %>% filter(germ==1, flower==1)
prior7.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model7.0 <- MCMCglmm(seed ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = ambient.SeedMaturation, prior = prior7.0, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Ambient/AM_model7.0.RData")
plot(AM_model7.0$VCV) # good trace plots
summary(AM_model7.0) # good effective sample size
autocorr.diag(AM_model7.0$VCV) # no autocorrelation
heidel.diag(AM_model7.0$VCV) # convergence success

#Heritability
mean(AM_model7.0[["VCV"]][, "animal"])
mean(AM_model7.0[["VCV"]][, "animalDom"])
mean(AM_model7.0[["VCV"]][, "matID"])
AM_herit7.0 <- AM_model7.0$VCV[, "animal"]/(AM_model7.0$VCV[, "animal"] + AM_model7.0$VCV[, "animalDom"] + AM_model7.0$VCV[, "matID"] + 1)
mean(AM_herit7.0)
posterior.mode(AM_herit7.0)
HPDinterval(AM_herit7.0)

#Transforming from latent to data scale using package "QGglmm"
df_A7.0 <- data.frame(
  va = as.vector(AM_model7.0[["VCV"]][, "animal"]),
  vd = as.vector(AM_model7.0[["VCV"]][, "animalDom"]),
  vm = as.vector(AM_model7.0[["VCV"]][, "matID"]),
  vp = rowSums(AM_model7.0[["VCV"]]))
yhat_A7.0 <- predict(AM_model7.0, type = "terms")

#Va
Va_A7.0 <- do.call("rbind", apply(df_A7.0, 1, function(row){
  QGparams(predict = yhat_A7.0, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Va_A7.0$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_A7.0$var.a.obs))

mean(as.mcmc(Va_A7.0$h2.obs)) #h2
HPDinterval(as.mcmc(Va_A7.0$h2.obs))

mu_A7.0 <- mean(as.mcmc(Va_A7.0$mean.obs))
mean(as.mcmc(Va_A7.0$h2.obs))/mu_A7.0 #Rw
HPDinterval(as.mcmc(Va_A7.0$h2.obs))/mu_A7.0

#Vd
Vd_A7.0 <- do.call("rbind", apply(df_A7.0, 1, function(row){
  QGicc(predict = yhat_A7.0, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vd_A7.0$var.comp.obs))
HPDinterval(as.mcmc(Vd_A7.0$var.comp.obs))

#Vm
Vm_A7.0 <- do.call("rbind", apply(df_A7.0, 1, function(row){
  QGicc(predict = yhat_A7.0, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "binom1.probit", verbose = FALSE)
}))
mean(as.mcmc(Vm_A7.0$var.comp.obs))
HPDinterval(as.mcmc(Vm_A7.0$var.comp.obs))


#Model 7.1 - No dominance
load(file="Routput/Ambient/AM_model7.1.RData")
plot(AM_model7.1$VCV) # good trace plots
summary(AM_model7.1) # good effective sample size
autocorr.diag(AM_model7.1$VCV) # no autocorrelation
heidel.diag(AM_model7.1$VCV) # convergence success













#### No. 8 Leaf Number - Univariate | Gaussian ####
#Note: For leaf number, we only include plants that germinated

#Calculating Variance for leaf number for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1) %>%
  summarise(sd=sd(leaf), var=(sd)^2) #var = 5.709704 *******

#Producing filtered data set for leaf number
ambient.leaf <- MCMC.ambient %>% filter(germ==1)

#Model 8.0 with additive, dominance, and maternal random effects and plot fixed effects
#Prior using inverse gamma (0.001, 0.001)
load(file="Routput/Ambient/AM_model8.0.RData")

ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior8.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model8.0 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.leaf, prior = prior8.0, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#AM_model8.0 used to have equal variance priors - redoing with Inverse Gamma Priors

plot(AM_model8.0$VCV) # poor trace plots
summary(AM_model8.0) # low effective sample size
autocorr.diag(AM_model8.0$VCV) # autocorrelation
heidel.diag(AM_model8.0$VCV) # convergence fail

#Model 8.1 with a Fisher Expanded Parameter Prior
ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior8.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model8.1 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.leaf, prior = prior8.1, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model8.1.RData")
plot(AM_model8.1$VCV) # good trace plots
summary(AM_model8.1) # good effective size
autocorr.diag(AM_model8.1$VCV) # no autocorrelation
heidel.diag(AM_model8.1$VCV) # convergence success

#Latent/Data Scale
AM_herit8.1 <- AM_model8.1$VCV[, "animal"]/(AM_model8.1$VCV[, "animal"] + AM_model8.1$VCV[, "animalDom"] + AM_model8.1$VCV[, "matID"] + AM_model8.1$VCV[, "units"])
mean(AM_herit8.1)
posterior.mode(AM_herit8.1) #Heritability
HPDinterval(AM_herit8.1) #h2 HPD

mean(AM_model8.1[["VCV"]][ , "animal"]) #Va
HPDinterval(AM_model8.1[["VCV"]][ , "animal"])

mean(AM_model8.1[["VCV"]][ , "animalDom"]) #Vd
HPDinterval(AM_model8.1[["VCV"]][ , "animalDom"])

mean(AM_model8.1[["VCV"]][ , "matID"]) #Vm
HPDinterval(AM_model8.1[["VCV"]][ , "matID"])

mean(AM_model8.1[["VCV"]][ , "units"]) #Vr
HPDinterval(AM_model8.1[["VCV"]][ , "units"])

mean(AM_model8.1$VCV[, "animal"] + AM_model8.1$VCV[, "animalDom"] + AM_model8.1$VCV[, "matID"]) # Vg
HPDinterval(AM_model8.1$VCV[, "animal"] + AM_model8.1$VCV[, "animalDom"] + AM_model8.1$VCV[, "matID"])


# Estimate trait mean as average across all plots (fixed) effects using posterior distributions
X_A8.1 <- AM_model8.1$X    # [n_obs x n_fixed]
fixed_effect_names_A8.1 <- colnames(X_A8.1) # fixed effect columns (intercept, plot2, plot3, etc)
Sol_fixed_A8.1 <- AM_model8.1$Sol[, fixed_effect_names_A8.1]  # [n_iter x n_fixed]; posterior subset to fixed effects only
fitted_vals_A8.1 <- X_A8.1 %*% t(Sol_fixed_A8.1) #design matrix multiplied by transpose of fixed effects
marginal_means_A8.1 <- colMeans(fitted_vals_A8.1) # average across individuals for each posterior sample
mean(marginal_means_A8.1) # summarize posterior mean and CI
HPDinterval(as.mcmc(marginal_means_A8.1))


mu_A8.1 <- mean(AM_model8.1[["Sol"]][,"(Intercept)"]) # mean based on intercept
mean(AM_model8.1[["VCV"]][ , "animal"]/marginal_means_A8.1) #Rw based on marginal means
HPDinterval(AM_model8.1[["VCV"]][ , "animal"]/marginal_means_A8.1)

# Check for significant difference between environments. Requires ambient data
mu_diff_8.1 <- marginal_means_A8.1 - marginal_means_A8.1
mean(mu_diff_8.1)
quantile(mu_diff_8.1, probs = c(0.025, 0.975))
mean(mu_diff_8.1 > 0)


#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_A8.1 <- data.frame(
  Va = as.vector(AM_model8.1[["VCV"]][ , "animal"]),
  Vd = as.vector(AM_model8.1[["VCV"]][ , "animalDom"]),
  Vm = as.vector(AM_model8.1[["VCV"]][ , "matID"]),
  Vr = as.vector(AM_model8.1[["VCV"]][ , "units"]),
  mu = as.vector(marginal_means_A8.1),
  h2 = as.vector(AM_model8.1$VCV[, "animal"]/(AM_model8.1$VCV[, "animal"] +
                                                AM_model8.1$VCV[, "matID"] + 
                                                AM_model8.1$VCV[, "animalDom"] + 
                                                AM_model8.1$VCV[, "units"])),
  Rw = as.vector(AM_model8.1[["VCV"]][ , "animal"]/marginal_means_A8.1))
posteriors_A8.1 <- posteriors_A8.1 %>% pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_A8.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(AM_model8.1[["VCV"]][ , "animal"])),
    mean(as.mcmc(AM_model8.1[["VCV"]][ , "animalDom"])),
    mean(as.mcmc(AM_model8.1[["VCV"]][ , "matID"])),
    mean(as.mcmc(AM_model8.1[["VCV"]][ , "units"])),
    mean(as.mcmc(marginal_means_A8.1)),
    mean(as.mcmc(AM_model8.1$VCV[, "animal"]/(AM_model8.1$VCV[, "animal"] +
                                                AM_model8.1$VCV[, "matID"] + 
                                                AM_model8.1$VCV[, "animalDom"] + 
                                                AM_model8.1$VCV[, "units"]))),
    mean(as.mcmc(AM_model8.1[["VCV"]][ , "animal"]/marginal_means_A8.1))
  ),
  Lower = c(
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "animal"]))[1],
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "animalDom"]))[1],
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "matID"]))[1],
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "units"]))[1],
    HPDinterval(as.mcmc(marginal_means_A8.1))[1],
    HPDinterval(as.mcmc(AM_model8.1$VCV[, "animal"]/(AM_model8.1$VCV[, "animal"] +
                                                       AM_model8.1$VCV[, "matID"] + 
                                                       AM_model8.1$VCV[, "animalDom"] + 
                                                       AM_model8.1$VCV[, "units"])))[1],
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "animal"]/marginal_means_A8.1))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "animal"]))[2],
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "animalDom"]))[2],
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "matID"]))[2],
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "units"]))[2],
    HPDinterval(as.mcmc(marginal_means_A8.1))[2],
    HPDinterval(as.mcmc(AM_model8.1$VCV[, "animal"]/(AM_model8.1$VCV[, "animal"] +
                                                       AM_model8.1$VCV[, "matID"] + 
                                                       AM_model8.1$VCV[, "animalDom"] + 
                                                       AM_model8.1$VCV[, "units"])))[2],
    HPDinterval(as.mcmc(AM_model8.1[["VCV"]][ , "animal"]/marginal_means_A8.1))[2]
  )
)

posteriors_A8.1$Component <- factor(posteriors_A8.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "Rw", "h2"))
save(posteriors_A8.1, stats_A8.1, file="Routput/Ambient/AM_stats8.1.RData")
load("Routput/Ambient/AM_stats8.1.RData")

var_A8.1_fig <- posteriors_A8.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
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
  theme(axis.title.x = ggtext::element_markdown())
var_A8.1_fig

ggsave(plot = var_A8.1_fig, filename = "Routput/Figures/fig_3_A8.1_var.png", width = 8, height = 8, units="cm")

var_A8.1_fig2 <- posteriors_A8.1 %>%
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
        legend.position  = "none")
var_A8.1_fig2

ggsave(plot = var_A8.1_fig2, filename = "Routput/Figures/fig_3_A8.1_h2.png", width = 8, height = 8, units="cm")

#Model 8.2 with a Fisher Expanded Parameter Prior (nu = 0.002 for prior sensitivity analysis)
ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior8.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

AM_model8.2 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "gaussian", data = ambient.leaf, prior = prior8.2, #gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model8.2.RData")
plot(AM_model8.2$VCV) # good trace plots
summary(AM_model8.2) # good (but lower) effective sample size
autocorr.diag(AM_model8.2$VCV) # no autocorrelation
heidel.diag(AM_model8.2$VCV) # convergence success
#No difference. Insensitive to prior change. 

#Model 8.3 - No dominance
load(file="Routput/Ambient/AM_model8.3.RData")
plot(AM_model8.3$VCV) # convergence success
summary(AM_model8.3) # good effective sample size
autocorr.diag(AM_model8.3$VCV) # no autocorrelation
heidel.diag(AM_model8.3$VCV) # convergence success

mean(AM_model8.3[["VCV"]][ , "animal"]) #Va
HPDinterval(AM_model8.3[["VCV"]][ , "animal"]) 

mean(AM_model8.3[["VCV"]][ , "matID"]) #Vm
HPDinterval(AM_model8.3[["VCV"]][ , "matID"])

mean(AM_model8.3[["VCV"]][ , "units"]) #Vr
HPDinterval(AM_model8.3[["VCV"]][ , "units"])

#Model 8.4 - Full with different smaller V prior
ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior8.4 <- list(R = list(V = 0.01, nu = 0.002),
                 G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

AM_model8.4 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "gaussian", data = ambient.leaf, prior = prior8.4, #gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model8.4.RData")
plot(AM_model8.4$VCV) # bad trace plots
summary(AM_model8.4) # poor effective sample size
autocorr.diag(AM_model8.4$VCV) # high autocorrelation
heidel.diag(AM_model8.4$VCV) # convergence fail






#### No. 9 Height - Univariate | Gaussian ####
#Note: For height, we only include plants that germinated

#Calculating Variance for height for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1, !height==0) %>%
  summarise(sd=sd(height), var=(sd)^2) #var = 147.7686 *******

#Model 9.0 with additive, dominance, and maternal random effects + plot fixed effects
#Prior with inverse gamma (0.001, 0.001)
ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior9.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model9.0 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.height, prior = prior9.0, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model9.0.RData")
plot(AM_model9.0$VCV) # poor trace plots
summary(AM_model9.0) # poor effective sample size
autocorr.diag(AM_model9.0$VCV) # high autocorrelation
heidel.diag(AM_model9.0$VCV) # convergence fail

#Model 9.1 with a Fisher Parameter Expanded Prior
ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior9.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model9.1 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.height, prior = prior9.1, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model9.1.RData")
plot(AM_model9.1$VCV) # good trace plots
summary(AM_model9.1) # good effective sample size
autocorr.diag(AM_model9.1$VCV) # no autocorrelation
heidel.diag(AM_model9.1$VCV) # convergence success

#Latent/Data Scale
AM_herit9.1 <- AM_model9.1$VCV[, "animal"]/(AM_model9.1$VCV[, "animal"] + AM_model9.1$VCV[, "animalDom"] + AM_model9.1$VCV[, "matID"] + AM_model9.1$VCV[, "units"])
mean(AM_herit9.1)
posterior.mode(AM_herit9.1) #h2
HPDinterval(AM_herit9.1)

mean(AM_model9.1[["VCV"]][ , "animal"]) #Va
HPDinterval(AM_model9.1[["VCV"]][ , "animal"])

mean(AM_model9.1[["VCV"]][ , "animalDom"]) #Vd
HPDinterval(AM_model9.1[["VCV"]][ , "animalDom"])

mean(AM_model9.1[["VCV"]][ , "matID"]) #Vm
HPDinterval(AM_model9.1[["VCV"]][ , "matID"])

mean(AM_model9.1[["VCV"]][ , "units"]) #Vr
HPDinterval(AM_model9.1[["VCV"]][ , "units"])

mean(AM_model9.1$VCV[, "animal"] + AM_model9.1$VCV[, "animalDom"] + AM_model9.1$VCV[, "matID"]) #Vg
HPDinterval(AM_model9.1$VCV[, "animal"] + AM_model9.1$VCV[, "animalDom"] + AM_model9.1$VCV[, "matID"])

# Estimate trait mean as average across all plots (fixed) effects using posterior distributions
X_A9.1 <- AM_model9.1$X    # [n_obs x n_fixed]
fixed_effect_names_A9.1 <- colnames(X_A9.1) # fixed effect columns (intercept, plot2, plot3, etc)
Sol_fixed_A9.1 <- AM_model9.1$Sol[, fixed_effect_names_A9.1]  # [n_iter x n_fixed]; posterior subset to fixed effects only
fitted_vals_A9.1 <- X_A9.1 %*% t(Sol_fixed_A9.1) #design matrix multiplied by transpose of fixed effects
marginal_means_A9.1 <- colMeans(fitted_vals_A9.1) # average across individuals for each posterior sample
mean(marginal_means_A9.1) # summarize posterior mean and CI
HPDinterval(as.mcmc(marginal_means_A9.1))

mu_A9.1 <- mean(AM_model9.1[["Sol"]][,"(Intercept)"]) # mean based on intercept
mean(AM_model9.1[["VCV"]][ , "animal"]/marginal_means_A9.1) #Rw based on marginal means
HPDinterval(AM_model9.1[["VCV"]][ , "animal"]/marginal_means_A9.1)

# Check for significant difference between environments. Requires ambient data
mu_diff_9.1 <- marginal_means_A9.1 - marginal_means_A9.1
mean(mu_diff_9.1)
quantile(mu_diff_9.1, probs = c(0.025, 0.975))
mean(mu_diff_9.1 > 0)


#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_A9.1 <- data.frame(
  Va = as.vector(AM_model9.1[["VCV"]][ , "animal"]),
  Vd = as.vector(AM_model9.1[["VCV"]][ , "animalDom"]),
  Vm = as.vector(AM_model9.1[["VCV"]][ , "matID"]),
  Vr = as.vector(AM_model9.1[["VCV"]][ , "units"]),
  mu = as.vector(marginal_means_A9.1),
  h2 = as.vector(AM_model9.1$VCV[, "animal"]/(AM_model9.1$VCV[, "animal"] +
                                                AM_model9.1$VCV[, "matID"] + 
                                                AM_model9.1$VCV[, "animalDom"] + 
                                                AM_model9.1$VCV[, "units"])),
  Rw = as.vector(AM_model9.1[["VCV"]][ , "animal"]/marginal_means_A9.1))
posteriors_A9.1 <- posteriors_A9.1 %>% pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_A9.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(AM_model9.1[["VCV"]][ , "animal"])),
    mean(as.mcmc(AM_model9.1[["VCV"]][ , "animalDom"])),
    mean(as.mcmc(AM_model9.1[["VCV"]][ , "matID"])),
    mean(as.mcmc(AM_model9.1[["VCV"]][ , "units"])),
    mean(as.mcmc(marginal_means_A9.1)),
    mean(as.mcmc(AM_model9.1$VCV[, "animal"]/(AM_model9.1$VCV[, "animal"] +
                                                AM_model9.1$VCV[, "matID"] + 
                                                AM_model9.1$VCV[, "animalDom"] + 
                                                AM_model9.1$VCV[, "units"]))),
    mean(as.mcmc(AM_model9.1[["VCV"]][ , "animal"]/marginal_means_A9.1))
  ),
  Lower = c(
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "animal"]))[1],
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "animalDom"]))[1],
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "matID"]))[1],
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "units"]))[1],
    HPDinterval(as.mcmc(marginal_means_A9.1))[1],
    HPDinterval(as.mcmc(AM_model9.1$VCV[, "animal"]/(AM_model9.1$VCV[, "animal"] +
                                                       AM_model9.1$VCV[, "matID"] + 
                                                       AM_model9.1$VCV[, "animalDom"] + 
                                                       AM_model9.1$VCV[, "units"])))[1],
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "animal"]/marginal_means_A9.1))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "animal"]))[2],
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "animalDom"]))[2],
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "matID"]))[2],
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "units"]))[2],
    HPDinterval(as.mcmc(marginal_means_A9.1))[2],
    HPDinterval(as.mcmc(AM_model9.1$VCV[, "animal"]/(AM_model9.1$VCV[, "animal"] +
                                                       AM_model9.1$VCV[, "matID"] + 
                                                       AM_model9.1$VCV[, "animalDom"] + 
                                                       AM_model9.1$VCV[, "units"])))[2],
    HPDinterval(as.mcmc(AM_model9.1[["VCV"]][ , "animal"]/marginal_means_A9.1))[2]
  )
)

posteriors_A9.1$Component <- factor(posteriors_A9.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "Rw", "h2"))
save(posteriors_A9.1, stats_A9.1, file="Routput/Ambient/AM_stats9.1.RData")
load("Routput/Ambient/AM_stats9.1.RData")

var_A9.1_fig <- posteriors_A9.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, 25) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown())
var_A9.1_fig

ggsave(plot = var_A9.1_fig, filename = "Routput/Figures/fig_3_A9.1_var.png", width = 8, height = 8, units="cm")

var_A9.1_fig2 <- posteriors_A9.1 %>%
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
        legend.position  = "none")
var_A9.1_fig2

ggsave(plot = var_A9.1_fig2, filename = "Routput/Figures/fig_3_A9.1_h2.png", width = 8, height = 8, units="cm")


#Model 9.2 with a Fisher Parameter Expanded Prior (nu = 0.002 for prior sensitivity analysis)
ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior9.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

AM_model9.2 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "gaussian", data = ambient.height, prior = prior9.2, #gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model9.2.RData")
plot(AM_model9.2$VCV) # ok trace plots. worse than before
summary(AM_model9.2) # low effective sample size
autocorr.diag(AM_model9.2$VCV) # high autocorrelation
heidel.diag(AM_model9.2$VCV) # convergence fail
#Sensitive to prior change. Model 9.1 is better. 

#Model 9.3 - No dominance
load(file="Routput/Ambient/AM_model9.3.RData")
plot(AM_model9.3$VCV) # good trace plots
summary(AM_model9.3) # good effective sample size
autocorr.diag(AM_model9.3$VCV) # no autocorrelation
heidel.diag(AM_model9.3$VCV) # convergence success

mean(AM_model9.3[["VCV"]][ , "animal"]) #Va
HPDinterval(AM_model9.3[["VCV"]][ , "animal"])

mean(AM_model9.3[["VCV"]][ , "matID"]) #Vm
HPDinterval(AM_model9.3[["VCV"]][ , "matID"])

mean(AM_model9.3[["VCV"]][ , "units"]) #Vr
HPDinterval(AM_model9.3[["VCV"]][ , "units"])

#Model 9.4 - Full with different smaller V prior
ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior9.4 <- list(R = list(V = 0.01, nu = 0.002),
                 G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

AM_model9.4 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "gaussian", data = ambient.height, prior = prior9.4, #gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model9.4.RData")
plot(AM_model9.4$VCV) # poor trace plots
summary(AM_model9.4) # low effective sample size
autocorr.diag(AM_model9.4$VCV) # high autocorrelation
heidel.diag(AM_model9.4$VCV) # convergence fail








#### No. 10 Stem Diameter - Univariate | Gaussian ####
#Note: For stem diameter, we only include plants that germinated

#Calculating Variance for stem diameter for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1, !stem_diam==0) %>%
  summarise(mean=mean(stem_diam), sd=sd(stem_diam), var=(sd)^2) #var = 1.623994 *******

#Model 10.0 additive, dominance, and maternal random effects + plot fixed effects
#Prior with inverse gamma (0.001, 0.001)
ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior10.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model10.0 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.stem, prior = prior10.0, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model10.0.RData")
plot(AM_model10.0$VCV) # poor trace plots 
summary(AM_model10.0) # poor effective sample size
autocorr.diag(AM_model10.0$VCV) # high autocorrelation
heidel.diag(AM_model10.0$VCV) # convergence fail

#Model 10.1 using a Fisher Parameter Expanded Prior
ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior10.1 <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                            G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                            G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model10.1 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                      ginverse = list(animal = Ainv, animalDom = Dinv),
                      family = "gaussian", data = ambient.stem, prior = prior10.1, #gaussian distribution
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model10.1.RData")
plot(AM_model10.1$VCV) # good trace plots
summary(AM_model10.1) # good effective sample size
autocorr.diag(AM_model10.1$VCV) # no autocorrelation
heidel.diag(AM_model10.1$VCV) # convergence success

#Latent/Data Scale
AM_herit10.1 <- AM_model10.1$VCV[, "animal"]/(AM_model10.1$VCV[, "animal"] + (AM_model10.1$VCV[, "animalDom"]) + AM_model10.1$VCV[, "matID"] + AM_model10.1$VCV[, "units"])
mean(AM_herit10.1)
posterior.mode(AM_herit10.1) #h2
HPDinterval(AM_herit10.1)

mean(AM_model10.1[["VCV"]][ , "animal"]) #Va
HPDinterval(AM_model10.1[["VCV"]][ , "animal"])

mean(AM_model10.1[["VCV"]][ , "animalDom"]) #Vd
HPDinterval(AM_model10.1[["VCV"]][ , "animalDom"])

mean(AM_model10.1[["VCV"]][ , "matID"]) #Vm
HPDinterval(AM_model10.1[["VCV"]][ , "matID"])

mean(AM_model10.1[["VCV"]][ , "units"]) #Vr
HPDinterval(AM_model10.1[["VCV"]][ , "units"])

mean(AM_model10.1$VCV[, "animal"] + (AM_model10.1$VCV[, "animalDom"]) + AM_model10.1$VCV[, "matID"]) #Vg
HPDinterval(AM_model10.1$VCV[, "animal"] + (AM_model10.1$VCV[, "animalDom"]) + AM_model10.1$VCV[, "matID"])

# Estimate trait mean as average across all plots (fixed) effects using posterior distributions
X_A10.1 <- AM_model10.1$X    # [n_obs x n_fixed]
fixed_effect_names_A10.1 <- colnames(X_A10.1) # fixed effect columns (intercept, plot2, plot3, etc)
Sol_fixed_A10.1 <- AM_model10.1$Sol[, fixed_effect_names_A10.1]  # [n_iter x n_fixed]; posterior subset to fixed effects only
fitted_vals_A10.1 <- X_A10.1 %*% t(Sol_fixed_A10.1) #design matrix multiplied by transpose of fixed effects
marginal_means_A10.1 <- colMeans(fitted_vals_A10.1) # average across individuals for each posterior sample
mean(marginal_means_A10.1) # summarize posterior mean and CI
HPDinterval(as.mcmc(marginal_means_A10.1))

mu_A10.1 <- mean(AM_model10.1[["Sol"]][,"(Intercept)"]) # mean based on intercept
mean(AM_model10.1[["VCV"]][ , "animal"]/marginal_means_A10.1) #Rw based on marginal means
HPDinterval(AM_model10.1[["VCV"]][ , "animal"]/marginal_means_A10.1)

# Check for significant difference between environments. Requires ambient data
mu_diff_10.1 <- marginal_means_A10.1 - marginal_means_A10.1
mean(mu_diff_10.1)
quantile(mu_diff_10.1, probs = c(0.025, 0.975))
mean(mu_diff_10.1 > 0)


#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_A10.1 <- data.frame(
  Va = as.vector(AM_model10.1[["VCV"]][ , "animal"]),
  Vd = as.vector(AM_model10.1[["VCV"]][ , "animalDom"]),
  Vm = as.vector(AM_model10.1[["VCV"]][ , "matID"]),
  Vr = as.vector(AM_model10.1[["VCV"]][ , "units"]),
  mu = as.vector(marginal_means_A10.1),
  h2 = as.vector(AM_model10.1$VCV[, "animal"]/(AM_model10.1$VCV[, "animal"] +
                                                AM_model10.1$VCV[, "matID"] + 
                                                AM_model10.1$VCV[, "animalDom"] + 
                                                AM_model10.1$VCV[, "units"])),
  Rw = as.vector(AM_model10.1[["VCV"]][ , "animal"]/marginal_means_A10.1))
posteriors_A10.1 <- posteriors_A10.1 %>% pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_A10.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(AM_model10.1[["VCV"]][ , "animal"])),
    mean(as.mcmc(AM_model10.1[["VCV"]][ , "animalDom"])),
    mean(as.mcmc(AM_model10.1[["VCV"]][ , "matID"])),
    mean(as.mcmc(AM_model10.1[["VCV"]][ , "units"])),
    mean(as.mcmc(marginal_means_A10.1)),
    mean(as.mcmc(AM_model10.1$VCV[, "animal"]/(AM_model10.1$VCV[, "animal"] +
                                                AM_model10.1$VCV[, "matID"] + 
                                                AM_model10.1$VCV[, "animalDom"] + 
                                                AM_model10.1$VCV[, "units"]))),
    mean(as.mcmc(AM_model10.1[["VCV"]][ , "animal"]/marginal_means_A10.1))
  ),
  Lower = c(
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "animal"]))[1],
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "animalDom"]))[1],
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "matID"]))[1],
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "units"]))[1],
    HPDinterval(as.mcmc(marginal_means_A10.1))[1],
    HPDinterval(as.mcmc(AM_model10.1$VCV[, "animal"]/(AM_model10.1$VCV[, "animal"] +
                                                       AM_model10.1$VCV[, "matID"] + 
                                                       AM_model10.1$VCV[, "animalDom"] + 
                                                       AM_model10.1$VCV[, "units"])))[1],
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "animal"]/marginal_means_A10.1))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "animal"]))[2],
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "animalDom"]))[2],
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "matID"]))[2],
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "units"]))[2],
    HPDinterval(as.mcmc(marginal_means_A10.1))[2],
    HPDinterval(as.mcmc(AM_model10.1$VCV[, "animal"]/(AM_model10.1$VCV[, "animal"] +
                                                       AM_model10.1$VCV[, "matID"] + 
                                                       AM_model10.1$VCV[, "animalDom"] + 
                                                       AM_model10.1$VCV[, "units"])))[2],
    HPDinterval(as.mcmc(AM_model10.1[["VCV"]][ , "animal"]/marginal_means_A10.1))[2]
  )
)

posteriors_A10.1$Component <- factor(posteriors_A10.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "Rw", "h2"))
save(posteriors_A10.1, stats_A10.1, file="Routput/Ambient/AM_stats10.1.RData")
load("Routput/Ambient/AM_stats10.1.RData")

var_A10.1_fig <- posteriors_A10.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, .4) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown())
var_A10.1_fig

ggsave(plot = var_A10.1_fig, filename = "Routput/Figures/fig_3_A10.1_var.png", width = 8, height = 8, units="cm")

var_A10.1_fig2 <- posteriors_A10.1 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  scale_color_manual(values = c("h2" = "#efbc82", "Rw" = "#675478")) +
  labs(x = expression("Proportion of variance (" * h^2 * " or " * Δ[evol] * bar(W) * ")"),
       y = "Density") +
  xlim(0, 0.2) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none")
var_A10.1_fig2

ggsave(plot = var_A10.1_fig2, filename = "Routput/Figures/fig_3_A10.1_h2.png", width = 8, height = 8, units="cm")

#Model 10.2 using a Fisher Parameter Expanded Prior
ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior10.2 <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

AM_model10.2 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                         ginverse = list(animal = Ainv, animalDom = Dinv),
                         family = "gaussian", data = ambient.stem, prior = prior10.2, #gaussian distribution
                         nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model10.2.RData")
plot(AM_model10.2$VCV) # good trace plots
summary(AM_model10.2) # good effective sample size
autocorr.diag(AM_model10.2$VCV) # no autocorrelation
heidel.diag(AM_model10.2$VCV) # convergence success
#No difference. Insensitive to prior change.

#Model 10.3 - No dominance
load(file="Routput/Ambient/AM_model10.3.RData")
plot(AM_model10.3$VCV) # good trace plots
summary(AM_model10.3) # good effective sample size
autocorr.diag(AM_model10.3$VCV) # no autocorrelation
heidel.diag(AM_model10.3$VCV) # convergence success

mean(AM_model10.3[["VCV"]][ , "animal"]) #Va
HPDinterval(AM_model10.3[["VCV"]][ , "animal"])

mean(AM_model10.3[["VCV"]][ , "matID"]) #Vm
HPDinterval(AM_model10.3[["VCV"]][ , "matID"])

mean(AM_model10.3[["VCV"]][ , "units"]) #Vr
HPDinterval(AM_model10.3[["VCV"]][ , "units"])

#Model 10.4 - Full with different smaller V prior
ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior10.4 <- list(R = list(V = 0.01, nu = 0.002),
                  G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

AM_model10.4 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                         ginverse = list(animal = Ainv, animalDom = Dinv),
                         family = "gaussian", data = ambient.stem, prior = prior10.4, #gaussian distribution
                         nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model10.4.RData")
plot(AM_model10.4$VCV) # bad trace plots
summary(AM_model10.4) # poor effective sample size
autocorr.diag(AM_model10.4$VCV) # high autocorrelation
heidel.diag(AM_model10.4$VCV) # convergence fail




#### No. 11 Number of Flowering Clusters - Univariate | Poisson ####
#Note: For Flowering Clusters, we only include plants that survived to reach flowering

#Calculating Variance for Flowering Clusters for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1 & flower==1) %>%
  summarise(mean=mean(flwr_clstr), sd=sd(flwr_clstr), var=(sd)^2) #var = 2.535085 *******

#Model 11.0 additive, dominance, and maternal random effects and plot fixed effects
#Prior with inverse gamma (0.001, 0.001)
ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior11.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model11.0 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "poisson", data = ambient.flwrclstr, prior = prior11.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model11.0.RData")
plot(AM_model11.0$VCV) # poor trace plot
summary(AM_model11.0) # poor effective sample size
autocorr.diag(AM_model11.0$VCV) # high autocorrelation
heidel.diag(AM_model11.0$VCV) # convergence fail

#Model 11.1 with Fisher Parameter Expanded Prior
ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior11.1 <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model11.1 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                      ginverse = list(animal = Ainv, animalDom = Dinv),
                      family = "poisson", data = ambient.flwrclstr, prior = prior11.1,
                      nitt = 8500000, thin = 4000, burnin = 500000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model11.1-extend2.RData")
plot(AM_model11.1$VCV) # moderate trace plots
summary(AM_model11.1) # good effective sample size
autocorr.diag(AM_model11.1$VCV) # no autocorrelation
heidel.diag(AM_model11.1$VCV) # convergence success

#AM_model11.1 has a nitt of 3100000, thin = 1500
#AM_model11.1-extend has a nitt of 5100000, thin = 2500
#AM_model11.1-extend2 has nitt of 8500000, thin = 4000, burnin 500000
#AM_model11.1.priorAlpV1 has alpha.V = 1 

#Latent Scale
mean(AM_model11.1[["VCV"]][ , "animal"]) #Va
mean(AM_model11.1[["VCV"]][ , "animalDom"]) #Vd
mean(AM_model11.1[["VCV"]][ , "matID"]) #Vm
mean(AM_model11.1[["VCV"]][ , "units"]) #Vr

HPDinterval(AM_model11.1$VCV[, "animal"])
AM_herit11.1 <- AM_model11.1$VCV[, "animal"]/(AM_model11.1$VCV[, "animal"] + AM_model11.1$VCV[, "animalDom"] + AM_model11.1$VCV[, "matID"] + AM_model11.1$VCV[, "units"])
mean(AM_herit11.1) #h2
posterior.mode(AM_herit11.1)
HPDinterval(AM_herit11.1)

mean(AM_model11.1$VCV[, "animal"] + AM_model11.1$VCV[, "animalDom"] + AM_model11.1$VCV[, "matID"]) #Vg
HPDinterval((AM_model11.1$VCV[, "animal"] + AM_model11.1$VCV[, "animalDom"] + AM_model11.1$VCV[, "matID"]))

#Transforming from latent to data scale using package "QGglmm"
df_A11.1 <- data.frame(
  va = as.vector(AM_model11.1[["VCV"]][, "animal"]),
  vd = as.vector(AM_model11.1[["VCV"]][, "animalDom"]),
  vm = as.vector(AM_model11.1[["VCV"]][, "matID"]),
  vr = as.vector(AM_model11.1[["VCV"]][, "units"]),
  vp = rowSums(AM_model11.1[["VCV"]]))
yhat_A11.1 <- predict(AM_model11.1, type = "terms")

#Va
Va_A11.1 <- do.call("rbind", apply(df_A11.1, 1, function(row){
  QGparams(predict = yhat_A11.1, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Va_A11.1$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_A11.1$var.a.obs))

mean(as.mcmc(Va_A11.1$h2.obs)) #h2
HPDinterval(as.mcmc(Va_A11.1$h2.obs))

mu_A11.1 <- mean(as.mcmc(Va_A11.1$mean.obs))
mu_A11.1
HPDinterval(as.mcmc(Va_A11.1$mean.obs))

mean(as.mcmc(Va_A11.1$var.a.obs))/mu_A11.1 #Rw
(HPDinterval(as.mcmc(Va_A11.1$var.a.obs)))/mu_A11.1

#Vd
Vd_A11.1 <- do.call("rbind", apply(df_A11.1, 1, function(row){
  QGicc(predict = yhat_A11.1, 
        var.comp = row[["vd"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vd_A11.1$var.comp.obs))
HPDinterval(as.mcmc(Vd_A11.1$var.comp.obs))

#Vm
Vm_A11.1 <- do.call("rbind", apply(df_A11.1, 1, function(row){
  QGicc(predict = yhat_A11.1, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vm_A11.1$var.comp.obs))
HPDinterval(as.mcmc(Vm_A11.1$var.comp.obs))

#Vr
Vr_A11.1 <- do.call("rbind", apply(df_A11.1, 1, function(row){
  QGicc(predict = yhat_A11.1, 
        var.comp = row[["vr"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vr_A11.1$var.comp.obs))
HPDinterval(as.mcmc(Vr_A11.1$var.comp.obs))


#Plotting density of variances on one figure
#First gather data into same dataframe
posteriors_A11.1 <- data.frame(
  Va = Va_A11.1$var.a.obs,
  Vd = Vd_A11.1$var.comp.obs,
  Vm = Vm_A11.1$var.comp.obs,
  Vr = Vr_A11.1$var.comp.obs,
  mu = Va_A11.1$mean.obs,
  h2 = Va_A11.1$h2.obs,
  Rw = (Va_A11.1$var.a.obs / Va_A11.1$mean.obs)) %>%
  pivot_longer(cols = everything(), names_to = "Component", values_to = "Values")

stats_A11.1 <- data.frame(
  Component = c("Va", "Vd", "Vm", "Vr", "mu", "h2", "Rw"),
  Mean = c(
    mean(as.mcmc(Va_A11.1$var.a.obs)),
    mean(as.mcmc(Vd_A11.1$var.comp.obs)),
    mean(as.mcmc(Vm_A11.1$var.comp.obs)),
    mean(as.mcmc(Vr_A11.1$var.comp.obs)),
    mean(as.mcmc(Va_A11.1$mean.obs)),
    mean(as.mcmc(Va_A11.1$h2.obs)),
    mean(as.mcmc(Va_A11.1$var.a.obs / Va_A11.1$mean.obs))
  ),
  Lower = c(
    HPDinterval(as.mcmc(Va_A11.1$var.a.obs))[1],
    HPDinterval(as.mcmc(Vd_A11.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vm_A11.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Vr_A11.1$var.comp.obs))[1],
    HPDinterval(as.mcmc(Va_A11.1$mean.obs))[1],
    HPDinterval(as.mcmc(Va_A11.1$h2.obs))[1],
    HPDinterval(as.mcmc(Va_A11.1$var.a.obs / Va_A11.1$mean.obs))[1]
  ),
  Upper = c(
    HPDinterval(as.mcmc(Va_A11.1$var.a.obs))[2],
    HPDinterval(as.mcmc(Vd_A11.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vm_A11.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Vr_A11.1$var.comp.obs))[2],
    HPDinterval(as.mcmc(Va_A11.1$mean.obs))[2],
    HPDinterval(as.mcmc(Va_A11.1$h2.obs))[2],
    HPDinterval(as.mcmc(Va_A11.1$var.a.obs / Va_A11.1$mean.obs))[2]
  )
)

posteriors_A11.1$Component <- factor(posteriors_A11.1$Component, levels = c("Vm", "Va", "Vd", "Vr", "mu", "Rw", "h2"))
save(posteriors_A11.1, stats_A11.1, file="Routput/Ambient/AM_stats11.1.RData")
load("Routput/Ambient/AM_stats11.1.RData")

var_A11.1_fig <- posteriors_A11.1 %>%
  filter(Component %in% c("Va", "Vd", "Vm")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  scale_fill_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  scale_color_manual(values = c("Vm" = "#297373", "Vd" = "#e2e260", "Va" = "dodgerblue2")) +
  labs(x = "Variance (Data Scale)", y = "Density") +
  xlim(0, 0.6) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = "none") +
  theme(axis.title.x = ggtext::element_markdown())
var_A11.1_fig

ggsave(plot = var_A11.1_fig, filename = "Routput/Figures/fig_3_A11.1_var.png", width = 8, height = 8, units="cm")

var_A11.1_fig2 <- posteriors_A11.1 %>%
  filter(Component %in% c("h2", "Rw")) %>%
  ggplot(aes(x = Values, fill = Component, color = Component)) +
  geom_density(alpha = 0.3, linewidth = 1.2, adjust = 1.5) +
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
        legend.position  = "none")
var_A11.1_fig2

ggsave(plot = var_A11.1_fig2, filename = "Routput/Figures/fig_3_A11.1_h2.png", width = 8, height = 8, units="cm")

#Model 11.2 with Fisher Parameter Expanded Prior (nu = 0.002 for prior sensitivity analysis + longer nitt)
ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior11.2 <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1),
                           G2 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1), 
                           G3 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1)))

AM_model11.2 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                         ginverse = list(animal = Ainv, animalDom = Dinv),
                         family = "poisson", data = ambient.flwrclstr, prior = prior11.2,
                         nitt = 4250000, thin = 2000, burnin = 250000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model11.2.RData")
plot(AM_model11.2$VCV) # OK trace plots
summary(AM_model11.2) # OK sample sizes (worse than 11.1)
autocorr.diag(AM_model11.2$VCV) # higher autocorrelation on additive and dominance
heidel.diag(AM_model11.2$VCV) # convergence success

#Model 11.3 - No dominance
load(file="Routput/Ambient/AM_model11.3.RData")
plot(AM_model11.3$VCV) # good trace plots
summary(AM_model11.3) # good effective sample size
autocorr.diag(AM_model11.3$VCV) # no autocorrelation
heidel.diag(AM_model11.3$VCV) # convergence success

mean(AM_model11.3[["VCV"]][ , "animal"]) #Va
mean(AM_model11.3[["VCV"]][ , "matID"]) #Vm
mean(AM_model11.3[["VCV"]][ , "units"]) #Vr


#Transforming from latent to data scale using package "QGglmm"
df_A11.3 <- data.frame(
  va = as.vector(AM_model11.3[["VCV"]][, "animal"]),
  vm = as.vector(AM_model11.3[["VCV"]][, "matID"]),
  vr = as.vector(AM_model11.3[["VCV"]][, "units"]),
  vp = rowSums(AM_model11.3[["VCV"]]))
yhat_A11.3 <- predict(AM_model11.3, type = "terms")

#Va
Va_A11.3 <- do.call("rbind", apply(df_A11.3, 1, function(row){
  QGparams(predict = yhat_A11.3, 
           var.a = row[["va"]],
           var.p = row[["vp"]],
           model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Va_A11.3$var.a.obs)) #Va
HPDinterval(as.mcmc(Va_A11.3$var.a.obs))

#Vm
Vm_A11.3 <- do.call("rbind", apply(df_A11.3, 1, function(row){
  QGicc(predict = yhat_A11.3, 
        var.comp = row[["vm"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vm_A11.3$var.comp.obs))
HPDinterval(as.mcmc(Vm_A11.3$var.comp.obs))

#Vr
Vr_A11.3 <- do.call("rbind", apply(df_A11.3, 1, function(row){
  QGicc(predict = yhat_A11.3, 
        var.comp = row[["vr"]],
        var.p = row[["vp"]],
        model = "Poisson.log", verbose = FALSE)
}))
mean(as.mcmc(Vr_A11.3$var.comp.obs))
HPDinterval(as.mcmc(Vr_A11.3$var.comp.obs))

#Model 11.4 - Full with different smaller V prior
ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior11.4 <- list(R = list(V = 0.01, nu = 0.002),
                  G = list(G1 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1),
                           G2 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1), 
                           G3 = list(V = 0.001, nu = 0.002, alpha.mu = 0, alpha.V = 1)))

AM_model11.4 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                         ginverse = list(animal = Ainv, animalDom = Dinv),
                         family = "poisson", data = ambient.flwrclstr, prior = prior11.4,
                         nitt = 4250000, thin = 2000, burnin = 250000, verbose = T, pr = TRUE)

load(file="Routput/Ambient/AM_model11.4.RData")
plot(AM_model11.4$VCV) # poor trace plots
summary(AM_model11.4) # low effective sample size
autocorr.diag(AM_model11.4$VCV) # high autocorrelation
heidel.diag(AM_model11.4$VCV) # convergence fail
