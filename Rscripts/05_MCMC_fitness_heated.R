#### PROJECT: Brassica rapa GxE Study (Data collected by Cameron So 2019 at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Calculate Va using MCMCglmm animal model
#### AUTHOR: Cameron So
#### DATE LAST MODIFIED: 2020/02/02

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
col_types_list2 <- cols_only(posID = "d", individual = "d", plot = col_factor(levels=c(1:12)),
                             animal = "f",
                             patID = "f", matID = "f", famID = col_factor(levels=c(1:62)),
                             treatment = col_factor(levels=c("A", "H")),
                             germ = "d", flower = "d", 
                             seed = "d", leaf = "d", flwr_clstr = "d", seed_pods = "d",
                             height = "d", stem_diam = "d", germ_census = "d", flwr_census = "d",
                             relative = "d")

MCMC.2019 <- read_csv("Rdata/MCMC_2019_cleaned.csv", col_names = TRUE, na = "NA", 
                      col_types=col_types_list2)
MCMC.2019$animal <- as.factor(MCMC.2019$animal)

ped <- read.csv("Rdata/heatarrays_animal_pedigree.csv") #note don't use readr to import the csv file.
for (x in 1:3) ped[, x] <- as.factor(ped[, x])

#Changing tibbles into dataframes because MCMCglmm cannot read tibbles
MCMC.2019 <- as.data.frame(MCMC.2019)
ped <- as.data.frame(ped)

MCMC.heated <- MCMC.2019 %>%
  filter(treatment=="H")

MCMC.heated$animal <- as.factor(MCMC.heated$animal)
lapply(ped, class)
lapply(MCMC.heated, class)

#'####################################################################'#
##############      MCMCglmm ANALYSIS FOR HEATED       ###############
#'####################################################################'#

#MCMCglmm uses an inverse-Gamma distribution, which is parametrized by two parameters nu and V^4. 

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

#### Creating Dominance Matrix for Vd estimates ####
Ainv <- inverseA(ped[, 1:3])$Ainv
Dinv <- makeD(ped[, 1:3])$Dinv
MCMC.heated$animalDom <- MCMC.heated$animal



















#### No. 1 Total Fitness - "Univariate" | Hurdle Poisson ####
#We set two different priors for the zero inflation process, and the Poisson process. However, according to Hadfield (https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/018802.html), the genetic correlation between the two processes must be set to 0 and cannot be estimated. 

MCMC.heated %>%
  summarise(mean=mean(seed_pods), sd=sd(seed_pods), var=(sd)^2) #var = 519.2583

#Model 1.2 with additive + dominance + maternal effects and plot fixed effects
#Zero-Inflated Poisson distribution + 2 different priors
#Level 1 = Poisson process, Level 2 = Zero inflation
load(file="Routput/MCMC_Heated_Model_1.2.RData")

heated.total <- MCMC.heated
prior1.2 <- list(R = list(V = diag(2), nu = 1, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G5 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G6 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) # Chi Square Prior with SD 1

HW_model1.2 <- MCMCglmm(seed_pods ~ trait - 1 + plot,
                     random = ~idh(at.level(trait,1)):animal + idh(at.level(trait,2)):animal +
                       idh(at.level(trait,1)):animalDom + idh(at.level(trait,2)):animalDom + 
                       idh(at.level(trait,1)):matID + idh(at.level(trait,2)):matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv), 
                     rcov = ~idh(trait):units,
                     family = "zipoisson", data = heated.total, prior=prior1.2, 
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(HW_model1.2, file="Routput/MCMC_Heated_Model_1.2.RData")

plot(HW_model1.2$VCV) # good trace plots
summary(HW_model1.2) # good effective sample size
autocorr.diag(HW_model1.2$VCV) # autocorrelation for dominance at Poisson process
heidel.diag(HW_model1.2$VCV) # convergence success

#Model 1.2b with Fisher Parameter Expanded Prior
#Note: we fit an overdispersed Poisson using the simple family='poisson' model. Apparently, Kruuk et al. 2014b had no problems running this model!
load(file="Routput/MCMC_Heated_Model_1.2b.RData")

heated.total <- MCMC.heated
prior1.2b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model1.2b <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID, 
                        ginverse = list(animal = Ainv, animalDom = Dinv),
                        family = "poisson", data = heated.total, prior=prior1.2b, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable have been removed. Use an informative prior!

plot(HW_model1.2b$VCV) # good trace plots
summary(HW_model1.2b) # good effective sample size
autocorr.diag(HW_model1.2b$VCV) # no autocorrelation
heidel.diag(HW_model1.2b$VCV) # convergence success for all

HW_herit1.2b <- HW_model1.2b$VCV[, "animal"]/(HW_model1.2b$VCV[, "animal"] + HW_model1.2b$VCV[, "animalDom"] +
                                              HW_model1.2b$VCV[, "matID"] + HW_model1.2b$VCV[, "units"] )
mean(HW_herit1.2b)
HPDinterval(HW_herit1.2b)

#Transformating data using "QGglmm" to get estimates of Va on data scale
yhat_H1.2b <- predict(HW_model1.2b, type = "terms")
mu_H1.2b <- mean(HW_model1.2b[["Sol"]][ , "(Intercept)"])
va_H1.2b <- mean(HW_model1.2b[["VCV"]][ , "animal"])
vp_H1.2b <- mean(rowSums(HW_model1.2b[["VCV"]]))

QG_H1.2b <- QGparams(predict=yhat_H1.2b, mu=mu_H1.2b, var.a=va_H1.2b, var.p=vp_H1.2b, model = "Poisson.log")
QG_H1.2b
va_H1.2b / vp_H1.2b


#Model 1.3 using Zero-Inflated Poisson distribution + 2 different priors
#Level 1 = Poisson process, Level 2 = Zero inflation
#Zero inflated poisson model
load(file="Routput/MCMC_Heated_Model_1.3.RData")

heated.total <- MCMC.heated
prior1.3 <- list(R = list(V = diag(2), nu = 1, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) # Chi Square Prior with SD 1

HW_model1.3 <- MCMCglmm(seed_pods ~ trait - 1 + plot,
                        random = ~idh(at.level(trait,1)):animal + idh(at.level(trait,2)):animal +
                                  idh(at.level(trait,1)):matID + idh(at.level(trait,2)):matID, 
                        ginverse = list(animal = Ainv), 
                        rcov = ~idh(trait):units,
                        family = "zipoisson", data = heated.total, prior=prior1.3, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable and have been removed. Use an informative prior

plot(HW_model1.3$VCV) # good trace plots
summary(HW_model1.3) # good effective sample size
autocorr.diag(HW_model1.3$VCV) # no autocorrelation
heidel.diag(HW_model1.3$VCV) # convergence success

#Model 1.3b
#Note: we fit an overdispersed Poisson using the simple family='poisson' model. Apparently, Kruuk et al. 2014b had no problems running this model!
load(file="Routput/MCMC_Heated_Model_1.3b.RData")

heated.total <- MCMC.heated
prior1.3b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model1.3b <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID, 
                         ginverse = list(animal = Ainv),
                         family = "poisson", data = heated.total, prior=prior1.3b, 
                         nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(HW_model1.3b$VCV) # good trace plots
summary(HW_model1.3b) # good effective sample size
autocorr.diag(HW_model1.3b$VCV) # no autocorrelation
heidel.diag(HW_model1.3b$VCV) # convergence success

HW_herit1.3b <- HW_model1.3b$VCV[, "animal"]/(HW_model1.3b$VCV[, "animal"] + HW_model1.3b$VCV[, "matID"] + HW_model1.3b$VCV[, "units"] )
mean(HW_herit1.3b)
HPDinterval(HW_herit1.3b)

yhat_H1.3b <- predict(HW_model1.3b, type = "terms")
mu_H1.3b <- mean(HW_model1.3b[["Sol"]][ , "(Intercept)"])
va_H1.3b <- mean(HW_model1.3b[["VCV"]][ , "animal"])
vp_H1.3b <- mean(rowSums(HW_model1.3b[["VCV"]]))


QG_H1.3b <- QGparams(predict=yhat_H1.3b, mu=mu_H1.3b, var.a=va_H1.3b, var.p=vp_H1.3b, model = "Poisson.log")
QG_H1.3b
va_H1.3b / vp_H1.3b

#














#### No. 2 Fecundity of Flowering Plants - Univariate | Poisson ####
#Note: For Fecundity, we only include plants that survived to reach flowering


#Calculating Variance for Fecundity for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1 & flower==1) %>%
  summarise(sd=sd(seed_pods), var=(sd)^2, mean=mean(seed_pods)) #var = 929.0913 *******

#Producing filtered data set for Fecundity
heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)

#From Committee Report #2: Fecundity follows a Poisson distribution

#Full model including Random: Va, Vm, and | Fixed: plot with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_2.1.RData")

heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior2.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model2.1 <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID +, 
                     ginverse = list(animal = Ainv),
                     family = "poisson", data = heated.Fecundity, prior = prior2.1, #poission distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(HW_model2.1$VCV) # 
summary(HW_model2.1) #
autocorr.diag(HW_model2.1$VCV) # 
heidel.diag(HW_model2.1$VCV) #

#Model 2.1b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Heated_Model_2.1b.RData")

heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior2.1b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model2.1b <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID +,
                      ginverse = list(animalDom = Dinv),
                      family = "poisson", data = heated.Fecundity, prior = prior2.1b, #poission distribution
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(HW_model2.1b, file="Routput/MCMC_Heated_Model_2.1b.RData")

plot(HW_model2.1b$VCV) # 
summary(HW_model2.1b) # 
autocorr.diag(HW_model2.1b$VCV) # 
heidel.diag(HW_model2.1b$VCV) # 

#Model 2.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_2.2.RData")

heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior2.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model2.2 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animalDom= Dinv),
                     family = "poisson", data = heated.Fecundity, prior = prior2.2,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(HW_model2.2, file="Routput/MCMC_Heated_Model_2.2.RData")

plot(HW_model2.2$VCV) # 
summary(HW_model2.2) # l
autocorr.diag(HW_model2.2$VCV) # 
heidel.diag(HW_model2.2$VCV) # 

#Model 2.2b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Heated_Model_2.2b.RData")

heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior2.2b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model2.2b <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "poisson", data = heated.Fecundity, prior = prior2.2b,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(HW_model2.2b$VCV) # good trace plots
summary(HW_model2.2b) # good effective sample size
autocorr.diag(HW_model2.2b$VCV) # no autocorrelation
heidel.diag(HW_model2.2b$VCV) # convergence success

HW_herit2.2b <- HW_model2.2b$VCV[, "animal"]/(HW_model2.2b$VCV[, "animal"] + HW_model2.2b$VCV[, "animalDom"] + HW_model2.2b$VCV[, "matID"] + HW_model2.2b$VCV[, "units"])
mean(HW_herit2.2b)
HPDinterval(HW_herit2.2b)

#Transformating data using "QGglmm" to get estimates of Va on data scale
yhat_H2.2b <- predict(HW_model2.2b, type = "terms")
mu_H2.2b <- mean(HW_model2.2b[["Sol"]][ , "(Intercept)"])
va_H2.2b <- mean(HW_model2.2b[["VCV"]][ , "animal"])
vp_H2.2b <- mean(rowSums(HW_model2.2b[["VCV"]]))

QG_H2.2b <- QGparams(predict=yhat_H2.2b, mu=mu_H2.2b, var.a=va_H2.2b, var.p=vp_H2.2b, model = "Poisson.log")
QG_H2.2b
HPDinterval(QG_H2.2b)
va_H2.2b / vp_H2.2b

#The phenotypic variance is so high for some reason... 

#Model 2.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_2.3.RData")

heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior2.3 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)))

HW_model2.3 <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "poisson", data = heated.Fecundity, prior = prior2.3,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use informative prior

plot(HW_model2.3$VCV) # good trace plots
summary(HW_model2.3) # poor effective sample size
autocorr.diag(HW_model2.3$VCV) # high autocorrelation
heidel.diag(HW_model2.3$VCV) # convergence success

#Model 2.3b attempted with Fisher Parameter Expanded Priors
load(file="Routput/MCMC_Heated_Model_2.3b.RData")

heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior2.3b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model2.3b <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "poisson", data = heated.Fecundity, prior = prior2.3b,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use informative prior

plot(HW_model2.3b$VCV) # good trace plots
summary(HW_model2.3b) # good effective sample size
autocorr.diag(HW_model2.3b$VCV) # no autocorrelation
heidel.diag(HW_model2.3b$VCV) # convergence success

HW_herit2.3b <- HW_model2.3b$VCV[, "animal"]/(HW_model2.3b$VCV[, "animal"] + 
                                              HW_model2.3b$VCV[, "matID"] + HW_model2.3b$VCV[, "units"])
mean(HW_herit2.3b)
HPDinterval(HW_herit2.3b)

#Transformating data using "QGglmm" to get estimates of Va on data scale
yhat_H2.3b <- predict(HW_model2.3b, type = "terms")
mu_H2.3b <- mean(HW_model2.3b[["Sol"]][ , "(Intercept)"])
va_H2.3b <- mean(HW_model2.3b[["VCV"]][ , "animal"])
vp_H2.3b <- mean(rowSums(HW_model2.3b[["VCV"]]))

QG_H2.3b <- QGparams(predict=yhat_H2.3b, mu=mu_H2.3b, var.a=va_H2.3b, var.p=vp_H2.3b, model = "Poisson.log")
QG_H2.3b
va_H2.3b / vp_H2.3b



#Model Comparison
HW_model2.2b$DIC
HW_model2.3b$DIC















#### No. 3 Survival to Flowering - Univariate | Bernoulli ####
#Note: For this trait, we include all plants from the experiment
#We use an uninformative prior suggested by in Animal tutorial papers for a Bernoulli distribution, which follows a X^2 distribution with df = 1.
#We apply a parameter extension with alpha.mu and alpha.V to allow for the use of the X^2 prior dist. 
#We also apply a probit link using family="ordinal" (this is switched to "threshold".. see below)


#Producing filtered data set for Survival
heated.Survival <- MCMC.heated

#From Committee Report #2, Survival follows a Bernoulli distribution
#Note: According to Pierre de Villemereuil, adding additional random factors is not recommended. 
#Further, it is required to fix residual variance (Vr) to 1 and estimate Va solely.
#We use a prior suggested by Villemereuil et al. 2012 for binary traits
#We change the family from "ordinal" to "threshold" as suggested by a forum post in 2017 by Jarrod Hadfield (https://stat.ethz.ch/pipermail/r-sig-mixed-models/2017q4/026115.html), and apply a truncation with trunc=TRUE

#Full model including Random: Va, Vm, and | Fixed: plot with Chi Sq df=1 prior
#Note, since autocorrelation is strong in MCMC w/ binary data, may need to multiply 'thin' and 'nitt' by 10.
load(file="Routput/MCMC_Heated_Model_3.1.RData")

heated.Survival <- MCMC.heated
prior3.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model3.1 <- MCMCglmm(flower ~ plot, random = ~animal + matID +, 
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = heated.Survival, prior = prior3.1, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(HW_model3.1$VCV) #
summary(HW_model3.1) #
autocorr.diag(HW_model3.1$VCV) #
heidel.diag(HW_model3.1$VCV) #


#Model 3.2 with additive, dominance and maternal effects + plot (fixed) using Chi Sq df=1 Prior
load(file="Routput/MCMC_Heated_Model_3.2.RData")

heated.Survival <- MCMC.heated
prior3.2 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model3.2 <- MCMCglmm(flower ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = heated.Survival, prior = prior3.2, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)


#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(HW_model3.2$VCV) # good trace plots
summary(HW_model3.2) # good effective sample size
autocorr.diag(HW_model3.2$VCV) # no autocorrelation
heidel.diag(HW_model3.2$VCV) # convergence success for all

HW_herit3.2 <- HW_model3.2$VCV[, "animal"]/(HW_model3.2$VCV[, "animal"] + HW_model3.2$VCV[, "animalDom"] + HW_model3.2$VCV[, "matID"] + 1)
mean(HW_herit3.2)
HPDinterval(HW_herit3.2)

yhat_H3.2 <- predict(HW_model3.2, type = "terms")
mu_H3.2 <- mean(HW_model3.2[["Sol"]][ , "(Intercept)"])
va_H3.2 <- mean(HW_model3.2[["VCV"]][ , "animal"])
vp_H3.2 <- mean(rowSums(HW_model3.2[["VCV"]]))

QG_H3.2 <- QGparams(predict=yhat_H3.2, mu=mu_H3.2, var.a=va_H3.2, var.p=vp_H3.2, model = "binom1.probit")
QG_H3.2
va_H3.2 / vp_H3.2


#Model 3.3 with animal + maternal effects
load(file="Routput/MCMC_Heated_Model_3.3.RData")

heated.Survival <- MCMC.heated
prior3.3 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model3.3 <- MCMCglmm(flower ~ plot, random = ~animal + matID, 
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = heated.Survival, prior = prior3.3, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model3.3$VCV) # good trace plots
summary(HW_model3.3) # good effective sample size
autocorr.diag(HW_model3.3$VCV) # no autocorrelation
heidel.diag(HW_model3.3$VCV) # convergence success
HW_herit3.3 <- HW_model3.3$VCV[, "animal"]/(HW_model3.3$VCV[, "animal"] + HW_model3.3$VCV[, "matID"] + 1)
mean(HW_herit3.3)
HPDinterval(HW_herit3.3)

yhat_H3.3 <- predict(HW_model3.3, type = "terms")
mu_H3.3 <- mean(HW_model3.3[["Sol"]][ , "(Intercept)"])
va_H3.3 <- mean(HW_model3.3[["VCV"]][ , "animal"])
vp_H3.3 <- mean(rowSums(HW_model3.3[["VCV"]]))

QG_H3.3 <- QGparams(predict=yhat_H3.3, mu=mu_H3.3, var.a=va_H3.3, var.p=vp_H3.3, model = "binom1.probit")
QG_H3.3
va_H3.3 / vp_H3.3

#Model Comparison
HW_model3.2$DIC
HW_model3.3$DIC















#### No. 4 Survival to Flowering + Fecundity - Bivariate ####
#Note: For this trait, we include all plants from the experiment
#Survival to Flowering = Bernoulli Dist | Fecundity = Poisson Dist

#Calculating Variance for Fecundity for the Non-informative Prior (equal variance)
MCMC.heated %>%
  summarise(mean=mean(seed_pods), sd=sd(seed_pods), var=(sd)^2) #mean = 9.261649
log(9.261649) # = 2.225882

MCMC.heated %>%
  summarise(sum.f=sum(flower), n=n(), mean=sum.f/n) #mean prob = 0.4423536
logit(0.4423536*(1-0.4423536))

#Bivariate Full Model
load(file="Routput/MCMC_Heated_Model_4.2.RData")

#Using prior similar to plasticity models
heated.Survival <- MCMC.heated
prior4.2 <- list(R = list(V = diag(2), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G3 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

HW_model4.2 <- MCMCglmm(cbind(flower, seed_pods) ~ trait + plot - 1,
                     random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                     ginverse = list(animal=Ainv, animalDom = Dinv), rcov = ~us(trait):units,
                     family = c("threshold", "poisson"), 
                     data = heated.Survival, prior = prior4.2, 
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(HW_model4.2, file="Routput/MCMC_Heated_Model_4.2.RData")

plot(HW_model4.2$VCV) #
summary(HW_model4.2) #
autocorr.diag(HW_model4.2$VCV) #


#Model 4.2b using alpha.V = diag(2) (basically at 1, suggested by David Aguirre)
#Did not run model 4.2 because expanded prior used in 4.2 ... this was discovered after Aguirre discussions
load(file="Routput/MCMC_Heated_Model_4.2b.RData")

heated.Survival <- MCMC.heated
prior4.2b <- list(R = list(V = diag(2), nu = 0, fix = 2),
                  G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*100),
                           G2 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*100), 
                           G3 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*100)))

HW_model4.2b <- MCMCglmm(cbind(flower, seed_pods) ~ trait + plot - 1,
                        random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                        ginverse = list(animal = Ainv, animalDom = Dinv), rcov = ~us(trait):units,
                        family = c("threshold", "poisson"), 
                        data = heated.Survival, prior = prior4.2b, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(HW_model4.2b, file="Routput/MCMC_Heated_Model_4.2b.RData")

plot(HW_model4.2b$VCV) #
summary(HW_model4.2b) #
autocorr.diag(HW_model4.2b$VCV) #



#Model 4.3 without dominance
load(file="Routput/MCMC_Heated_Model_4.3.RData")

heated.Survival <- MCMC.heated
prior4.3 <- list(R = list(V = diag(2), nu=2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

HW_model4.3 <- MCMCglmm(cbind(flower, seed_pods) ~ trait + plot - 1,
                        random = ~us(trait):animal + us(trait):matID,
                        ginverse = list(animal=Ainv), rcov = ~us(trait):units,
                        family = c("threshold", "poisson"), 
                        data = heated.Survival, prior = prior4.3, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(HW_model4.3$VCV) # trace plots don't look good for matID... but that's ok.
summary(HW_model4.3) # good effective sample size
autocorr.diag(HW_model4.3$VCV) # no autocorrelation
heidel.diag(HW_model4.3$VCV) # convergence fail for matID but whatever b/c not interested in it

#!!! 

#Calculating Genetic Correlations = Covariance b/w trait 1 & 2 / sqrt (var [trait 1] * var [trait 2])

# (1) Survival & Fecundity
gen.corrH4.3 <-HW_model4.3$VCV[,'traitflower:traitseed_pods.animal']/
  sqrt(HW_model4.3$VCV[,'traitflower:traitflower.animal']*HW_model4.3$VCV[,'traitseed_pods:traitseed_pods.animal']) 
mean(gen.corrH4.3) #Post Mean = 0.7220879
posterior.mode(gen.corrH4.3) #Post Mode = 0.7236313
HPDinterval(gen.corrH4.3) # Posterior 95% CI = (0.6951492, 0.7464828)

phen.corrH4.3 <- cor.test(heated.Survival$flower, heated.Survival$seed_pods, test="pearson")
phen.corrH4.3 # Corr = 0.4564391 | Posterior 95% CI =  (0.4291918, 0.4828564)









#### No. 5 Flowering Success - Univariate | Bernoulli ####
#Note: For this trait, we include only include plants that germinated
#We use an uninformative prior suggested by in Animal tutorial papers for a Bernoulli distribution, which follows a X^2 distribution with df = 1.
#We apply a parameter extension with alpha.mu and alpha.V to allow for the use of the X^2 prior dist. 
#We also apply a probit link using family="ordinal" (this is switched to "threshold".. see below)

#Producing filtered data set for Fecundity
heated.Flowering <- MCMC.heated %>% filter(germ==1)

#From Committee Report #2, Survival follows a Bernoulli distribution
#NEED TO CONFIRM THE PRIORS FOR BERNOULLI DISTRIBUTION

#Full model including Random: Va, Vm, and | Fixed: plot with Chi Sq df=1 prior
#Note, since autocorrelation is strong in MCMC w/ binary data, may need to multiply 'thin' and 'nitt' by 10.
load(file="Routput/MCMC_Heated_Model_5.1.RData")

heated.Flowering <- MCMC.heated %>% filter(germ==1)
prior5.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model5.1 <- MCMCglmm(flower ~ plot, random = ~animal + matID + ,
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = heated.Flowering, prior = prior5.1, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(HW_model5.1, file="Routput/MCMC_Heated_Model_5.1.RData")

plot(HW_model5.1$VCV) # 
summary(HW_model5.1) # 
autocorr.diag(HW_model5.1$VCV) # 
heidel.diag(HW_model5.1$VCV) # 

#Model with additive, dominance and maternal effects + plot (fixed) using Chi Sq df=1 Prior
load(file="Routput/MCMC_Heated_Model_5.2.RData")

heated.Flowering <- MCMC.heated %>% filter(germ==1)
prior5.2 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model5.2 <- MCMCglmm(flower ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = heated.Flowering, prior = prior5.2, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(HW_model5.2$VCV) # good trace plots
summary(HW_model5.2) # good effective sample size
autocorr.diag(HW_model5.2$VCV) # no autocorrelation
heidel.diag(HW_model5.2$VCV) # convergence success

HW_herit5.2 <- HW_model5.2$VCV[, "animal"]/(HW_model5.2$VCV[, "animal"] + HW_model5.2$VCV[, "animalDom"] + HW_model5.2$VCV[, "matID"] + 1)
mean(HW_herit5.2)
HPDinterval(HW_herit5.2)

yhat_H5.2 <- predict(HW_model5.2, type = "terms")
mu_H5.2 <- mean(HW_model5.2[["Sol"]][ , "(Intercept)"])
va_H5.2 <- mean(HW_model5.2[["VCV"]][ , "animal"])
vp_H5.2 <- mean(rowSums(HW_model5.2[["VCV"]]))

QG_H5.2 <- QGparams(predict=yhat_H5.2, mu=mu_H5.2, var.a=va_H5.2, var.p=vp_H5.2, model = "binom1.probit")
QG_H5.2
va_H5.2 / vp_H5.2

#Model with additive and maternal effects + fixed plot using a Chi Sq df=1 prior
load(file="Routput/MCMC_Heated_Model_5.3.RData")

heated.Flowering <- MCMC.heated %>% filter(germ==1)
prior5.3 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model5.3 <- MCMCglmm(flower ~ plot, random = ~animal + matID, 
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = heated.Flowering, prior = prior5.3, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model5.3$VCV) # good trace plots
summary(HW_model5.3) # good effective sample size
autocorr.diag(HW_model5.3$VCV) # no autocorrelation
heidel.diag(HW_model5.3$VCV) # convergence success

HW_herit5.3 <- HW_model5.3$VCV[, "animal"]/(HW_model5.3$VCV[, "animal"] + HW_model5.3$VCV[, "matID"] + 1)
mean(HW_herit5.3)
HPDinterval(HW_herit5.3)

yhat_H5.3 <- predict(HW_model5.3, type = "terms")
mu_H5.3 <- mean(HW_model5.3[["Sol"]][ , "(Intercept)"])
va_H5.3 <- mean(HW_model5.3[["VCV"]][ , "animal"])
vp_H5.3 <- mean(rowSums(HW_model5.3[["VCV"]]))

QG_H5.3 <- QGparams(predict=yhat_H5.3, mu=mu_H5.3, var.a=va_H5.3, var.p=vp_H5.3, model = "binom1.probit")
QG_H5.3
va_H5.3 / vp_H5.3

#Model Comparison
HW_model5.2$DIC
HW_model5.3$DIC















#### No. 6 Germination Success - Univariate | Bernoulli ####
#Note: For this trait, we include all plants from the study
#We use an uninformative prior suggested by in Animal tutorial papers for a Bernoulli distribution, which follows a X^2 distribution with df = 1.
#We apply a parameter extension with alpha.mu and alpha.V to allow for the use of the X^2 prior dist. 
#We also apply a probit link using family = "threshold"

#Producing filtered data set for Germination Success
heated.Germination <- MCMC.heated

#From Committee Report #2, Germination Success follows a Bernoulli distribution

#Full model including Random: Va, Vm, and | Fixed: plot with Chi Sq df=1 prior
#Note, since autocorrelation is strong in MCMC w/ binary data, may need to multiply 'thin' and 'nitt' by 10.
load(file="Routput/MCMC_Heated_Model_6.1.RData")

heated.Germination <- MCMC.heated
prior6.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model6.1 <- MCMCglmm(germ ~ plot, random = ~animal + matID +,
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = heated.Germination, prior = prior6.1, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(HW_model6.1, file="Routput/MCMC_Heated_Model_6.1.RData")

plot(HW_model6.1$VCV) #
summary(HW_model6.1) #
autocorr.diag(HW_model6.1$VCV) #
heidel.diag(HW_model6.1$VCV) #

#Model with additive, dominance and maternal effects + plot (fixed) using Chi Sq df=1 Prior
load(file="Routput/MCMC_Heated_Model_6.2.RData")

heated.Germination <- MCMC.heated
prior6.2 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model6.2 <- MCMCglmm(germ ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = heated.Germination, prior = prior6.2, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(HW_model6.2$VCV) # good trace plots
summary(HW_model6.2) # good effective sample size
autocorr.diag(HW_model6.2$VCV) # no autocorrelation
heidel.diag(HW_model6.2$VCV) # convergence success

HW_herit6.2 <- HW_model6.2$VCV[, "animal"]/(HW_model6.2$VCV[, "animal"] + HW_model6.2$VCV[, "animalDom"] + HW_model6.2$VCV[, "matID"] + 1)
mean(HW_herit6.2)
HPDinterval(HW_herit6.2)

yhat_H6.2 <- predict(HW_model6.2, type = "terms")
mu_H6.2 <- mean(HW_model6.2[["Sol"]][ , "(Intercept)"])
va_H6.2 <- mean(HW_model6.2[["VCV"]][ , "animal"])
vp_H6.2 <- mean(rowSums(HW_model6.2[["VCV"]]))

QG_H6.2 <- QGparams(predict=yhat_H6.2, mu=mu_H6.2, var.a=va_H6.2, var.p=vp_H6.2, model = "binom1.probit")
QG_H6.2
va_H6.2 / vp_H6.2


#Model with additive and maternal effects + fixed plot using a Chi Sq df=1 prior
load(file="Routput/MCMC_Heated_Model_6.3.RData")

heated.Germination <- MCMC.heated
prior6.3 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model6.3 <- MCMCglmm(germ ~ plot, random = ~animal + matID, 
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = heated.Germination, prior = prior6.3, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model6.3$VCV) # good trace plots
summary(HW_model6.3) # good effective sample size
autocorr.diag(HW_model6.3$VCV) # no autocorrelation
heidel.diag(HW_model6.3$VCV) # convergence success

HW_herit6.3 <- HW_model6.3$VCV[, "animal"]/(HW_model6.3$VCV[, "animal"] + HW_model6.3$VCV[, "matID"] + 1)
mean(HW_herit6.3)
HPDinterval(HW_herit6.3)

yhat_H6.3 <- predict(HW_model6.3, type = "terms")
mu_H6.3 <- mean(HW_model6.3[["Sol"]][ , "(Intercept)"])
va_H6.3 <- mean(HW_model6.3[["VCV"]][ , "animal"])
vp_H6.3 <- mean(rowSums(HW_model6.3[["VCV"]]))

QG_H6.3 <- QGparams(predict=yhat_H6.3, mu=mu_H6.3, var.a=va_H6.3, var.p=vp_H6.3, model = "binom1.probit")
QG_H6.3
va_H6.3 / vp_H6.3

#Model Comparison
HW_model6.2$DIC
HW_model6.3$DIC














#### No. 7 Seed Maturation Success - Univariate | Bernoulli ####
#Note: For this trait, we include only include plants that germinated and flowered
#We use an uninformative prior suggested by in Animal tutorial papers for a Bernoulli distribution, which follows a X^2 distribution with df = 1.
#We apply a parameter extension with alpha.mu and alpha.V to allow for the use of the X^2 prior dist. 
#We also apply a probit link using family = "threshold"


#Producing filtered data set for Seed Maturation Success
heated.SeedMaturation <- MCMC.heated %>% filter(germ==1, flower==1)

#From Committee Report #2, Seed Maturaiton Success follows a Bernoulli distribution

#Full model including Random: Va, Vm, and | Fixed: plot with Chi Sq df=1 prior
#Note, since autocorrelation is strong in MCMC w/ binary data, may need to multiply 'thin' and 'nitt' by 10.
load(file="Routput/MCMC_Heated_Model_7.1.RData")

heated.SeedMaturation <- MCMC.heated %>% filter(germ==1, flower==1)
prior7.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model7.1 <- MCMCglmm(seed ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = heated.SeedMaturation, prior = prior7.1, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(HW_model7.1, file="Routput/MCMC_Heated_Model_7.1.RData")

plot(HW_model7.1$VCV) # 
summary(HW_model7.1) # 
autocorr.diag(HW_model7.1$VCV) # 
heidel.diag(HW_model7.1$VCV) # 

#Model with additive, dominance and maternal effects + plot (fixed) using Chi Sq df=1 Prior
load(file="Routput/MCMC_Heated_Model_7.2.RData")

heated.SeedMaturation <- MCMC.heated %>% filter(germ==1, flower==1)
prior7.2 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model7.2 <- MCMCglmm(seed ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = heated.SeedMaturation, prior = prior7.2, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(HW_model7.2$VCV) # good trace plots
summary(HW_model7.2) #good effective sample size
autocorr.diag(HW_model7.2$VCV) # no autocorrelation
heidel.diag(HW_model7.2$VCV) # convergence success

HW_herit7.2 <- HW_model7.2$VCV[, "animal"]/(HW_model7.2$VCV[, "animal"] + HW_model7.2$VCV[, "animalDom"] + HW_model7.2$VCV[, "matID"] + 1)
mean(HW_herit7.2)
HPDinterval(HW_herit7.2)


yhat_H7.2 <- predict(HW_model7.2, type = "terms")
mu_H7.2 <- mean(HW_model7.2[["Sol"]][ , "(Intercept)"])
va_H7.2 <- mean(HW_model7.2[["VCV"]][ , "animal"])
vp_H7.2 <- mean(rowSums(HW_model7.2[["VCV"]]))

QG_H7.2 <- QGparams(predict=yhat_H7.2, mu=mu_H7.2, var.a=va_H7.2, var.p=vp_H7.2, model = "binom1.probit")
QG_H7.2
va_H7.2 / vp_H7.2

#Model with additive and maternal effects + fixed plot using a Chi Sq df=1 prior
load(file="Routput/MCMC_Heated_Model_7.3.RData")

heated.SeedMaturation <- MCMC.heated %>% filter(germ==1, flower==1)
prior7.3 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_model7.3 <- MCMCglmm(seed ~ plot, random = ~animal + matID, 
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = heated.SeedMaturation, prior = prior7.3, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model7.3$VCV) # good trace plots
summary(HW_model7.3) # good effective sample size
autocorr.diag(HW_model7.3$VCV) # no autocorrelation
heidel.diag(HW_model7.3$VCV) # convergence success

HW_herit7.3 <- HW_model7.3$VCV[, "animal"]/(HW_model7.3$VCV[, "animal"] + HW_model7.3$VCV[, "matID"] + 1)
mean(HW_herit7.3)
HPDinterval(HW_herit7.3)

yhat_H7.3 <- predict(HW_model7.3, type = "terms")
mu_H7.3 <- mean(HW_model7.3[["Sol"]][ , "(Intercept)"])
va_H7.3 <- mean(HW_model7.3[["VCV"]][ , "animal"])
vp_H7.3 <- mean(rowSums(HW_model7.3[["VCV"]]))

QG_H7.3 <- QGparams(predict=yhat_H7.3, mu=mu_H7.3, var.a=va_H7.3, var.p=vp_H7.3, model = "binom1.probit")
QG_H7.3
va_H7.3 / vp_H7.3

#Model Comparison
HW_model7.2$DIC
HW_model7.3$DIC


















#### No. 8 Leaf Number - Univariate | Gaussian ####
#Note: For leaf number, we only include plants that germinated
#We use an uninformative prior that divides the total phenotypic variance by the number of random effects included
#Variance in Total Fitness = squared standard deviation
#We treat plot as a random effect as from previous analysis (01_Data_exploration), we observe that there is a plot effect on this trait... There is likely some sampling error caused by plot due to the placement of the physical plots at KSR with respect to the tree line, etc so we want to give less statistical weight to the effect of plot. We also want to estimate the effect of plot 

#Preliminary GLMM analysis


#Calculating Variance for leaf number for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1) %>%
  summarise(sd=sd(leaf), var=(sd)^2) #var = 6.928249 *******

#Producing filtered data set for leaf number
heated.leaf <- MCMC.heated %>% filter(germ==1)

#From Committee Report #2: Leaf Number follows a Gaussian distribution

#Full model including 3 random effects: Va, Vd, and Vm + 2 fixed effects: plot and IGEs
load(file="Routput/MCMC_Heated_Model_8.1.RData")

heated.leaf <- MCMC.heated %>% filter(germ==1)
prior8.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model8.1 <- MCMCglmm(leaf ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = heated.leaf, prior = prior8.1, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(HW_Model8.1, file="Routput/MCMC_Heated_Model_8.1.RData")

plot(HW_model8.1$VCV) #
summary(HW_model8.1) #
autocorr.diag(HW_model8.1$VCV) #
heidel.diag(HW_model8.1$VCV) #


#Model 2.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_8.2.RData")

heated.leaf <- MCMC.heated %>% filter(germ==1)
prior8.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model8.2 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom= Dinv),
                     family = "gaussian", data = heated.leaf, prior = prior8.2, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(HW_model8.2$VCV) # poor trace plots
summary(HW_model8.2) # bad effective sample size
autocorr.diag(HW_model8.2$VCV) # autocorrelation
heidel.diag(HW_model8.2$VCV) # convergence fail

#Model 8.2b with a Fisher Expanded Parameter Prior
load(file="Routput/MCMC_Heated_Model_8.2b.RData")

heated.leaf <- MCMC.heated %>% filter(germ==1)
prior8.2b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model8.2b <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom= Dinv),
                     family = "gaussian", data = heated.leaf, prior = prior8.2b, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(HW_model8.2b$VCV) # good trace plots
summary(HW_model8.2b) # good effective sample size
autocorr.diag(HW_model8.2b$VCV) # no autocorrelation
heidel.diag(HW_model8.2b$VCV) # convergence success

HW_herit8.2b <- HW_model8.2b$VCV[, "animal"]/(HW_model8.2b$VCV[, "animal"] + HW_model8.2b$VCV[, "matID"] + HW_model8.2b$VCV[, "animalDom"] + HW_model8.2b$VCV[, "units"])
mean(HW_herit8.2b)
HPDinterval(HW_herit8.2b)



#Model 8.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_8.3.RData")

heated.leaf <- MCMC.heated %>% filter(germ==1)
prior8.3 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)))

HW_model8.3 <- MCMCglmm(leaf ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = heated.leaf, prior = prior8.3, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model8.3$VCV) # good trace plot
summary(HW_model8.3) # decent effective sample size
autocorr.diag(HW_model8.3$VCV) # no autocorrelation
heidel.diag(HW_model8.3$VCV) # convergence success

HW_herit8.3 <- HW_model8.3$VCV[, "animal"]/(HW_model8.3$VCV[, "animal"] + HW_model8.3$VCV[, "matID"] + HW_model8.3$VCV[, "units"])
mean(HW_herit8.3)
HPDinterval(HW_herit8.3)

#WE REPORT THE VALUES FROM 8.3

#Model 8.3b using Fisher Parameter Expanded Priors
load(file="Routput/MCMC_Heated_Model_8.3b.RData")

heated.leaf <- MCMC.heated %>% filter(germ==1)
prior8.3b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model8.3b <- MCMCglmm(leaf ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = heated.leaf, prior = prior8.3b, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model8.3b$VCV) # good trace plots
summary(HW_model8.3b) #  good effective sample size
autocorr.diag(HW_model8.3b$VCV) # no autocorrelation
heidel.diag(HW_model8.3b$VCV) # convergence success

HW_herit8.3b <- HW_model8.3b$VCV[, "animal"]/(HW_model8.3b$VCV[, "animal"] + HW_model8.3b$VCV[, "matID"] + HW_model8.3b$VCV[, "units"])
mean(HW_herit8.3b)
HPDinterval(HW_herit8.3b)


#Model Comparison
HW_model8.2b$DIC
HW_model8.3b$DIC














#### No. 9 Height - Univariate | Gaussian ####
#Note: For height, we only include plants that germinated
#We use an uninformative prior that divides the total phenotypic variance by the number of random effects included
#Variance in Total Fitness = squared standard deviation
#We treat plot as a random effect as from previous analysis (01_Data_exploration), we observe that there is a plot effect on this trait... There is likely some sampling error caused by plot due to the placement of the physical plots at KSR with respect to the tree line, etc so we want to give less statistical weight to the effect of plot. We also want to estimate the effect of plot 

#Calculating Variance for height for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1, !height==0) %>%
  summarise(sd=sd(height), var=(sd)^2) #var = 213.7212 *******


#Producing filtered data set for height
heated.height <- MCMC.heated %>% filter(germ==1, !height==0)

#From Committee Report #2: height follows a Gaussian distribution

#Full model including Random: Va, Vm, and | Fixed: plot with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_9.1.RData")

heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior9.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model9.1 <- MCMCglmm(height ~ plot, random = ~animal + matID +, 
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = heated.height, prior = prior9.1, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(HW_model9.1, file="Routput/MCMC_Heated_Model_9.1.RData")

plot(HW_model9.1$VCV) #
summary(HW_model9.1) #
autocorr.diag(HW_model9.1$VCV) #
heidel.diag(HW_model9.1$VCV) #


#Model 9.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_9.2.RData")

heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior9.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model9.2 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = heated.height, prior = prior9.2, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(HW_model9.2$VCV) # poor trace plots
summary(HW_model9.2) # poor effective sample size
autocorr.diag(HW_model9.2$VCV) # high autocorrelation
heidel.diag(HW_model9.2$VCV) # convergence fail

#Model 9.2b with a Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Heated_Model_9.2b.RData")

heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior9.2b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model9.2b <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = heated.height, prior = prior9.2b, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(HW_model9.2b$VCV) # ok trace plots...
summary(HW_model9.2b) # good effective sample size
autocorr.diag(HW_model9.2b$VCV) # no autocorrelation
heidel.diag(HW_model9.2b$VCV) # convergence success

HW_herit9.2b <- HW_model9.2b$VCV[, "animal"]/(HW_model9.2b$VCV[, "animal"] + HW_model9.2b$VCV[, "animalDom"] + 
                                              HW_model9.2b$VCV[, "matID"] + HW_model9.2b$VCV[, "units"])
mean(HW_herit9.2b)
HPDinterval(HW_herit9.2b)

#Testing for correlation between additive and dominance sampling
library(lmodel2)
animal.est <- HW_model9.2b$VCV[, "animal"][1:dim(HW_model9.2b$VCV)[1]]
animalDom.est <- HW_model9.2b$VCV[, "animalDom"][1:dim(HW_model9.2b$VCV)[1]]
mareg <- lmodel2(animal.est~animalDom.est)
x11(w = 8, h = 8)
plot(mareg, method = "MA",
     xlab = "animalDom estimates", ylab = "animal estimates",
     main = paste("Sampling correlation: ", round(mareg$r, 3), sep =""),
     sub = "Line represents the major axis regression") 

#Model 9.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_9.3.RData")

heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior9.3 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)))

HW_model9.3 <- MCMCglmm(height ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = heated.height, prior = prior9.3, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model9.3$VCV) # OK trace plots
summary(HW_model9.3) # good effective sample size
autocorr.diag(HW_model9.3$VCV) # no autocorrelation
heidel.diag(HW_model9.3$VCV) # convergence success 

#Model 9.3b with standard uninformative priors + increased thin and nitt
load(file="Routput/MCMC_Heated_Model_9.3b.RData")

heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior9.3b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model9.3b <- MCMCglmm(height ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = heated.height, prior = prior9.3b, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model9.3b$VCV) # good trace plots
summary(HW_model9.3b) #  good effective sample size
autocorr.diag(HW_model9.3b$VCV) # no autocorrelation
heidel.diag(HW_model9.3b$VCV) # convergence success

HW_herit9.3b <- HW_model9.3b$VCV[, "animal"]/(HW_model9.3b$VCV[, "animal"] + 
                                              HW_model9.3b$VCV[, "matID"] + HW_model9.3b$VCV[, "units"])
mean(HW_herit9.3b)
HPDinterval(HW_herit9.3b)

#Model Comparison
HW_model9.2b$DIC
HW_model9.3b$DIC










#### No. 10 Height*Leaf*Seed_Pod - Multivariate | Gaus,Gaus,Pois ####
#Note: For this G matrix, we only use plants that germinated

#Producing filtered data set for height
heated.trivar <- MCMC.heated %>% filter(germ==1, !height==0)

#Model 10.2 with 3 additive, dominance and maternal effects + fixed plot effects
#In the Punetes et al. (2016) paper, they use a prior with nu=4.001 for Inverse-Wishart

heated.trivar <- MCMC.heated %>% filter(germ==1, !height==0)
prior10.2 <- list(R=list(V=diag(3), nu=2.002), 
                  G = list(G1 = list(V=diag(3), nu = 2.002),
                           G2 = list(V=diag(3), nu = 2.002),
                           G3 = list(V=diag(3), nu = 2.002)))

HW_model10.2 <- MCMCglmm(cbind(height, leaf, seed_pods) ~ trait + plot - 1,
                      random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                      ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~us(trait):units,
                      family = c("gaussian","gaussian","poisson"), 
                      data = heated.trivar, prior = prior10.2,
                      nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(HW_model10.2, file="Routput/MCMC_Heated_Model_10.2.RData")

plot(HW_model10.2$VCV) #
summary(HW_model10.2) #
autocorr.diag(HW_model10.2$VCV) #
heidel.diag(HW_model10.2$VCV) #

#Model 10.2b with Fisher Parameter Expansion + Parallel Processing
load(file="Routput/MCMC_Heated_Model_10.2b.RData")

heated.trivar <- MCMC.heated %>% filter(germ==1, !height==0)
prior10.2b <- list(R=list(V=diag(3), nu=2.002), 
                  G = list(G1 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000),
                           G2 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000),
                           G3 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000)))


HW_model10.2b <- MCMCglmm(cbind(height, leaf, seed_pods) ~ trait + plot - 1,
           random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~us(trait):units,
           family = c("gaussian","gaussian","poisson"), 
           data = heated.trivar, prior = prior10.2b,
           nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(HW_model10.2b, file="Routput/MCMC_Heated_Model_10.2b.RData")

plot(HW_model10.2b$VCV) #
summary(HW_model10.2b) #
autocorr.diag(HW_model10.2b$VCV) #
heidel.diag(HW_model10.2b$VCV) #

#Model 10.3 with 3 Additive and Maternal Effects + plot
load(file="Routput/MCMC_Heated_Model_10.3.RData")

heated.trivar <- MCMC.heated %>% filter(germ==1, !height==0)
prior10.3 <- list(R=list(V=diag(3), nu=2.002), 
                  G = list(G1 = list(V=diag(3), nu = 2.002),
                           G2 = list(V=diag(3), nu = 2.002)))

HW_model10.3 <- MCMCglmm(cbind(height, leaf, seed_pods) ~ trait + plot - 1,
                      random = ~us(trait):animal + us(trait):matID,
                      ginverse = list(animal = Ainv), rcov=~us(trait):units,
                      family = c("gaussian","gaussian","poisson"), 
                      data = heated.trivar, prior = prior10.3,
                      nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model10.3$VCV) # fairly good trace plots
summary(HW_model10.3) # low effective sample size for matID
autocorr.diag(HW_model10.3$VCV) # autocorrelation in matID
heidel.diag(HW_model10.3$VCV) # convergence fail for some matID

#Model 10.3b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Heated_Model_10.3b.RData")

heated.trivar <- MCMC.heated %>% filter(germ==1, !height==0)
prior10.3b <- list(R=list(V=diag(3), nu=2.002), 
                  G = list(G1 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000),
                           G2 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000)))

HW_model10.3b <- MCMCglmm(cbind(height, leaf, seed_pods) ~ trait - 1,
                         random = ~us(trait):animal + us(trait):matID,
                         ginverse = list(animal = Ainv), rcov=~us(trait):units,
                         family = c("gaussian","gaussian","poisson"), 
                         data = heated.trivar, prior = prior10.3b,
                         nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model10.3b$VCV) # OK trace plots..
summary(HW_model10.3b) # good effective sample size
autocorr.diag(HW_model10.3b$VCV) # some autocorrelation in animal and matID (at cusp of 0.1)
autocorr.plot(HW_model10.3b$VCV)
heidel.diag(HW_model10.3b$VCV) #convergence fail for some

#We use this model for the Thesis.. Ambient uses inverse gamma priors but the difference is not marginal

#Calculating Genetic Correlations = Covariance b/w trait 1 & 2 / sqrt (var [trait 1] * var [trait 2])
#Trait Combinations: 
# (1) Height & Leaf Number
# (2) Height & Seed Pod Number
# (3) Leaf & Seed Pod NUmber

# (1) Height & Leaf Number
gen.corrH10.3b_1 <-HW_model10.3b$VCV[,'traitheight:traitleaf.animal']/
  sqrt(HW_model10.3b$VCV[,'traitheight:traitheight.animal']*HW_model10.3b$VCV[,'traitleaf:traitleaf.animal']) 
mean(gen.corrH10.3b_1) #Post Mean = 0.7997674
posterior.mode(gen.corrH10.3b_1) #Post Mode = 0.8384592
HPDinterval(gen.corrH10.3b_1) # Posterior 95% CI = (0.6036196, 0.982449)

phen.corrH10.3b_1 <- cor.test(heated.trivar$height, heated.trivar$leaf, test="pearson")
phen.corrH10.3b_1 # Corr = 0.6200271 | Posterior 95% CI =  (0.5870698, 0.6509326)


# (2) Height & Seed Pod Number
gen.corrH10.3b_2 <-HW_model10.3b$VCV[,'traitheight:traitseed_pods.animal']/
  sqrt(HW_model10.3b$VCV[,'traitheight:traitheight.animal']*HW_model10.3b$VCV[,'traitseed_pods:traitseed_pods.animal']) 
mean(gen.corrH10.3b_2) #Post Mean = 0.5455568
posterior.mode(gen.corrH10.3b_2) #Post Mode = 0.5712875
HPDinterval(gen.corrH10.3b_2) # Posterior 95% CI = (0.1832114, 0.9307405)

phen.corrH10.3b_2 <- cor.test(heated.trivar$height, heated.trivar$seed_pods, test="pearson")
phen.corrH10.3b_2 # Corr = 0.492376 | Poster 95% CI = (0.4520913, 0.5306561)


# (3) Leaf Number and Seed Pod Number
gen.corrH10.3b_3 <-HW_model10.3b$VCV[,'traitseed_pods:traitleaf.animal']/
  sqrt(HW_model10.3b$VCV[,'traitseed_pods:traitseed_pods.animal']*HW_model10.3b$VCV[,'traitleaf:traitleaf.animal']) 
mean(gen.corrH10.3b_3) #Post Mean = 0.8227794
posterior.mode(gen.corrH10.3b_3) #Post Mode = 0.901749
HPDinterval(gen.corrH10.3b_3) #Posterior 95% CI = (0.6473549, 0.9822817)

phen.corrH10.3b_3 <- cor.test(heated.trivar$seed_pods, heated.trivar$leaf, test="pearson")
phen.corrH10.3b_3 # Corr = 0.5733401| Posterior 95% CI = (0.5374894, 0.6071220)










#### No. 11 Stem Diameter - Univariate | Gaussian ####
#Note: For stem diameter, we only include plants that germinated
#We use an uninformative prior that divides the total phenotypic variance by the number of random effects included
#Variance in Total Fitness = squared standard deviation
#We treat plot as a random effect as from previous analysis (01_Data_exploration), we observe that there is a plot effect on this trait... There is likely some sampling error caused by plot due to the placement of the physical plots at KSR with respect to the tree line, etc so we want to give less statistical weight to the effect of plot. We also want to estimate the effect of plot 

#Calculating Variance for stem diameter for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1, !stem_diam==0) %>%
  summarise(mean=mean(stem_diam), sd=sd(stem_diam), var=(sd)^2) #var = 3.102643 *******

#Producing filtered data set for stem diameter
heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)

#From Committee Report #2: stem diameter follows a Polynomial distribution

#Full model including Random: Va, Vm, and | Fixed: plot with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_11.1.RData")

heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior11.1 <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 0.002),
                           G2 = list(V = 1, nu = 0.002), 
                           G3 = list(V = 1, nu = 0.002)))

HW_model11.1 <- MCMCglmm(stem_diam ~ plot, random = ~animal + matID +, 
                      ginverse = list(animal = Ainv),
                      family = "gaussian", data = heated.stem, prior = prior11.1, #gaussian distribution
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(HW_model11.1, file="Routput/MCMC_Heated_Model_11.1.RData")

plot(HW_model11.1$VCV) # 
summary(HW_model11.1) #
autocorr.diag(HW_model11.1$VCV) # 
heidel.diag(HW_model11.1$VCV) # 


#Model 11.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_11.2.RData")

heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior11.2 <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 0.002),
                            G2 = list(V = 1, nu = 0.002), 
                            G3 = list(V = 1, nu = 0.002)))

HW_model11.2 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                       ginverse = list(animal = Ainv, animalDom= Dinv),
                       family = "gaussian", data = heated.stem, prior = prior11.2, #gaussian distribution
                       nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable have been removed. Use an informative prior!

plot(HW_model11.2$VCV) # poor trace plots
summary(HW_model11.2) # poor effective sample size
autocorr.diag(HW_model11.2$VCV) # high autocorrelation
heidel.diag(HW_model11.2$VCV) # convergence fail

#Model 11.2b using a Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Heated_Model_11.2b.RData")

heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior11.2b <- list(R = list(V = 1, nu = 1),
                   G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                            G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                            G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model11.2b <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                       ginverse = list(animal = Ainv, animalDom= Dinv),
                       family = "gaussian", data = heated.stem, prior = prior11.2b, #gaussian distribution
                       nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable have been removed. Use an informative prior!

plot(HW_model11.2b$VCV) # good trace plots
summary(HW_model11.2b) # good effective sample size
autocorr.diag(HW_model11.2b$VCV) # no autocorrelation
heidel.diag(HW_model11.2b$VCV) # convergence success

HW_herit11.2b <- HW_model11.2b$VCV[, "animal"]/(HW_model11.2b$VCV[, "animal"] + HW_model11.2b$VCV[, "animalDom"] + HW_model11.2b$VCV[, "matID"] + HW_model11.2b$VCV[, "units"])
mean(HW_herit11.2b)
HPDinterval(HW_herit11.2b)


#Model 11.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_11.3.RData")

heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior11.3 <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 0.002),
                            G2 = list(V = 1, nu = 0.002)))

HW_model11.3 <- MCMCglmm(stem_diam ~ plot, random = ~animal + matID,
                       ginverse = list(animal = Ainv),
                       family = "gaussian", data = heated.stem, prior = prior11.3, #gaussian distribution
                       nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable have been removed. Use an informative prior!

plot(HW_model11.3$VCV) # trace plots not great for matID
summary(HW_model11.3) # low effective sample size
autocorr.diag(HW_model11.3$VCV) # high autocorrelation
heidel.diag(HW_model11.3$VCV) #  convergence fail

#Model 11.3b using a Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Heated_Model_11.3b.RData")

heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior11.3b <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                            G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model11.3b <- MCMCglmm(stem_diam ~ plot, random = ~animal + matID,
                       ginverse = list(animal = Ainv),
                       family = "gaussian", data = heated.stem, prior = prior11.3b, #gaussian distribution
                       nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable have been removed. Use an informative prior!

plot(HW_model11.3b$VCV) # good trace plots
summary(HW_model11.3b) # good effective sample size
autocorr.diag(HW_model11.3b$VCV) # no autocorrelation
heidel.diag(HW_model11.3b$VCV) # convergene success for all

HW_herit11.3b <- HW_model11.3b$VCV[, "animal"]/(HW_model11.3b$VCV[, "animal"] + 
                                                HW_model11.3b$VCV[, "matID"] + HW_model11.3b$VCV[, "units"])
mean(HW_herit11.3b)
HPDinterval(HW_herit11.3b)














#### No. 12 Five Trait - Multivariate | Ord, Gaus,Gaus, Ord, Pois ####
#Note: For this G matrix, we only use plants that germinated

#Producing filtered data set for height
heated.multi <- MCMC.heated %>% filter(germ==1)

#Model 12.2 with 3 additive, dominance and maternal effects + fixed plot effects
#We use a alpha.V value of 1 because we're using a prior that is suitable for an ordinal distribution
#Although family='ordinal' is available, we use family='threshold' because it is better at mixing and is suggested by Hadfield (2015)
prior12.2 <- list(R=list(V=diag(5), fix = 1), 
                  G = list(G1 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G2 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G3 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5))))

HW_model12.2 <- MCMCglmm(cbind(germ_census, leaf, height, flwr_census, seed_pods) ~ trait + plot - 1,
                         random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                         ginverse = list(animalDom = Dinv), rcov=~us(trait):units,
                         family = c("threshold", "gaussian","gaussian", "threshold", "poisson"), 
                         data = heated.multi, prior = prior12.2,
                         nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(HW_model12.2, file="Routput/MCMC_Heated_Model_12.2.RData")


plot(HW_model12.2$VCV) #
summary(HW_model12.2) #
autocorr.diag(HW_model12.2$VCV) #


#Model 12.2b with 3 additive, dominance and maternal effects + fixed plot effects
#We fix the residual covariance matrix at near 0 according to Hadfield (2015) Supplement 3
#Additionally, we also fix 
### STILL WORKING ON THIS

heated.multi <- MCMC.heated %>% filter(germ==1)
prior12.2b <- list(R=list(V=diag(5)*1e-15, fix = 1), 
                  G = list(G1 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G2 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G3 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5))))


HW_model12.2b <- MCMCglmm(cbind(germ_census, leaf, height, flwr_census, seed_pods) ~ trait + plot - 1,
                         random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                         ginverse = list(animalDom = Dinv), rcov=~us(trait):units,
                         family = c("threshold", "gaussian","gaussian", "threshold", "poisson"), 
                         data = heated.multi, prior = prior12.2b,
                         nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(HW_model12.2b, file="Routput/MCMC_Heated_Model_12.2b.RData")


plot(HW_model12.2b$VCV) #
summary(HW_model12.2b) #
autocorr.diag(HW_model12.2b$VCV) #

#Model 12.3 with 3 additive effects but no plot effect
heated.multi <- MCMC.heated %>% filter(germ==1)

prior12.3 <- list(R=list(V=diag(5), fix = 1), 
                  G = list(G1 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G2 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G3 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5))))

HW_model12.3 <- MCMCglmm(cbind(germ_census, leaf, height, flwr_census, seed_pods) ~ trait - 1,
                         random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                         ginverse = list(animalDom = Dinv), rcov=~us(trait):units,
                         family = c("threshold", "gaussian","gaussian", "threshold", "poisson"), 
                         data = heated.multi, prior = prior12.3,
                         nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(HW_model12.3, file="Routput/MCMC_Heated_Model_12.3.RData")

plot(HW_model12.3$VCV) #
summary(HW_model12.3) #
autocorr.diag(HW_model12.3$VCV) #















#### No. 13 Number of Flowering Clusters - Univariate | Poisson ####
#Note: For Flowering Clusters, we only include plants that survived to reach flowering


#Calculating Variance for Flowering Clusters for the Non-informative Prior (equal variance)
MCMC.heated %>%
  filter(germ==1 & flower==1) %>%
  summarise(mean=mean(flwr_clstr), sd=sd(flwr_clstr), var=(sd)^2) #var = 9.736969 

#Producing filtered data set for Flowering Clusters
heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)

#Full model including Random: Va, Vm, and | Fixed: plot with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_13.1.RData")

heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior13.1 <- list(R = list(V = 1, nu = 0.002), 
                  G = list(G1 = list(V = 1, nu = 0.002),
                           G2 = list(V = 1, nu = 0.002), 
                           G3 = list(V = 1, nu = 0.002)))

HW_model13.1 <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal=Ainv),
                     family = "poisson", data = heated.flwrclstr, prior = prior13.1, #poission distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(HW_model13.1, file="Routput/MCMC_Heated_Model_13.1.RData")

plot(HW_model13.1$VCV) # 
summary(HW_model13.1) # 
autocorr.diag(HW_model13.1$VCV) # 
heidel.diag(HW_model13.1$VCV) # 

#Model 13.1b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Heated_Model_13.1b.RData")

heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior13.1b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model13.1b <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal=Ainv),
                     family = "poisson", data = heated.flwrclstr, prior = prior13.1b, #poission distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(HW_model13.1b, file="Routput/MCMC_Heated_Model_13.1b.RData")

plot(HW_model13.1b$VCV) # 
summary(HW_model13.1b) # 
autocorr.diag(HW_model13.1b$VCV) # 
heidel.diag(HW_model13.1b$VCV) # 




#Model 13.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_13.2.RData")

heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior13.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

HW_model13.2 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "poisson", data = heated.flwrclstr, prior = prior13.2,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(HW_model13.2$VCV) # poor trace plot
summary(HW_model13.2) # poor effective sample size
autocorr.diag(HW_model13.2$VCV) # high autocorrelation
heidel.diag(HW_model13.2$VCV) # convergence fail

#Model 13.2b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Heated_Model_13.2b.RData")

heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior13.2b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model13.2b <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                      ginverse = list(animalDom= Dinv),
                      family = "poisson", data = heated.flwrclstr, prior = prior13.2b,
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)


#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effets, but use an informative prior!

plot(HW_model13.2b$VCV) # good trace plots
summary(HW_model13.2b) # medicore effective sample sizes
autocorr.diag(HW_model13.2b$VCV) # high autocorrelation
heidel.diag(HW_model13.2b$VCV) # convergence success

HW_herit13.2b <- HW_model13.2b$VCV[, "animal"]/(HW_model13.2b$VCV[, "animal"] + HW_model13.2b$VCV[, "animalDom"] + HW_model13.2b$VCV[, "matID"])
mean(HW_herit13.2b)
HPDinterval(HW_herit13.2b)

#Transformating data using "QGglmm" to get estimates of Va on data scale
yhat_H13.2b <- predict(HW_model13.2b, type = "terms")
mu_H13.2b <- mean(HW_model13.2b[["Sol"]][ , "(Intercept)"])
va_H13.2b <- mean(HW_model13.2b[["VCV"]][ , "animal"])
vp_H13.2b <- mean(rowSums(HW_model13.2b[["VCV"]]))

QG_H13.2b <- QGparams(predict=yhat_H13.2b, mu=mu_H13.2b, var.a=va_H13.2b, var.p=vp_H13.2b, model = "Poisson.log")
QG_H13.2b
va_H13.2b / vp_H13.2b






#Model 13.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Heated_Model_13.3.RData")

heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior13.3 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)))

HW_model13.3 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "poisson", data = heated.flwrclstr, prior = prior13.3,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model13.3$VCV) # good trace plots
summary(HW_model13.3) # OK effective sample size
autocorr.diag(HW_model13.3$VCV) # high autocorrelation
heidel.diag(HW_model13.3$VCV) # convergence success

#Model 13.3b attempted with Fisher Parameter Expanded Priors
load(file="Routput/MCMC_Heated_Model_13.3b.RData")

heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior13.3b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_model13.3b <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "poisson", data = heated.flwrclstr, prior = prior13.3b,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(HW_model13.3b$VCV) # good trace plots
summary(HW_model13.3b) # good effective sample size
autocorr.diag(HW_model13.3b$VCV) # no autocorrelation
heidel.diag(HW_model13.3b$VCV) # convergence success

HW_herit13.3b <- HW_model13.3b$VCV[, "animal"]/(HW_model13.3b$VCV[, "animal"] + HW_model13.3b$VCV[, "matID"] + HW_model13.3b$VCV[, "units"])
mean(HW_herit13.3b)
HPDinterval(HW_herit13.3b)

yhat_H13.3b <- predict(HW_model13.3b, type = "terms")
mu_H13.3b <- mean(HW_model13.3b[["Sol"]][ , "(Intercept)"])
va_H13.3b <- mean(HW_model13.3b[["VCV"]][ , "animal"])
vp_H13.3b <- mean(rowSums(HW_model13.3b[["VCV"]]))

QG_H13.3b <- QGparams(predict=yhat_H13.3b, mu=mu_H13.3b, var.a=va_H13.3b, var.p=vp_H13.3b, model = "Poisson.log")
QG_H13.3b
va_H13.3b / vp_H13.3b


#Model Comparison
HW_model13.2b$DIC
HW_model13.3b$DIC
