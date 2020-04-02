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

MCMC.ambient <- MCMC.2019 %>%
  filter(treatment=="A")



#'####################################################################'#
##############      MCMCglmm ANALYSIS FOR AMBIENT      ###############
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
MCMC.ambient$animalDom <- MCMC.ambient$animal














#### No. 1 Total Fitness - "Univariate" | Zero-Inflated Poisson ####
#We set two different priors for the zero inflation process, and the Poisson process. However, according to Hadfield (https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/018802.html), the genetic correlation between the two processes must be set to 0 and cannot be estimated. 

MCMC.ambient %>%
  summarise(mean=mean(seed_pods), sd=sd(seed_pods), var=(sd)^2) #var = 136.3176

#Model 1.2 with additive + dominance + maternal effects and plot fixed effects
#Zero-Inflated Poisson distribution + 2 different priors
#Level 1 = Poisson process, Level 2 = Zero inflation
load(file="Routput/MCMC_Ambient_Model_1.2.RData")
ambient.total <- MCMC.ambient

prior1.2 <- list(R = list(V = diag(2), nu = 1, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G5 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G6 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) # Chi Square Prior with SD 1

AM_model1.2 <- MCMCglmm(seed_pods ~ trait - 1 + plot,
                     random = ~idh(at.level(trait,1)):animal + idh(at.level(trait,2)):animal +
                       idh(at.level(trait,1)):animalDom + idh(at.level(trait,2)):animalDom + 
                       idh(at.level(trait,1)):matID + idh(at.level(trait,2)):matID, 
                     ginverse = list(animal=Ainv, animalDom= Dinv), 
                     rcov = ~idh(trait):units,
                     family = "zipoisson", data = ambient.total, prior=prior1.2, 
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr=TRUE)

#I don't even know how to diagnose this model.. also I forgot to set pr=TRUE so nothing was saved!

plot(AM_model1.2$VCV) # good trace plots
summary(AM_model1.2) # good effective sample size
autocorr.diag(AM_model1.2$VCV) # no autocorrelation
heidel.diag(AM_model1.2$VCV) # convergence success

#Model 1.2b
#Note: we fit an overdispersed Poisson using the simple family='poisson' model. Apparently, Kruuk et al. 2014b had no problems running this model!
load(file="Routput/MCMC_Ambient_Model_1.2b.RData")

ambient.total <- MCMC.ambient

prior1.2b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model1.2b <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID, 
                         ginverse = list(animal=Ainv, animalDom = Dinv),
                         family = "poisson", data = ambient.total, prior=prior1.2b, 
                         nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects could not be estimated and have been removed. Use singular=OK to sample these effects, but use an informative prior!

plot(AM_model1.2b$VCV) # good trace plots
summary(AM_model1.2b) # good effective sample size
autocorr.diag(AM_model1.2b$VCV) # no autocorrelation
heidel.diag(AM_model1.2b$VCV) # convergence success

#Transformating data using "QGglmm" to get estimates of Va on data scale
yhat_A1.2b <- predict(AM_model1.2b, type = "terms")
mu_A1.2b <- mean(AM_model1.2b[["Sol"]][ , "(Intercept)"])
va_A1.2b <- mean(AM_model1.2b[["VCV"]][ , "animal"])
vp_A1.2b <- mean(rowSums(AM_model1.2b[["VCV"]]))

QG_A1.2b <- QGparams(predict=yhat_A1.2b, mu=mu_A1.2b, var.a=va_A1.2b, var.p=vp_A1.2b, model = "Poisson.log")
QG_A1.2b
va_A1.2b / vp_A1.2b


#Model 1.3 with no dominance component
load(file="Routput/MCMC_Ambient_Model_1.3.RData")

ambient.total <- MCMC.ambient
prior1.3 <- list(R = list(V = diag(2), nu = 1, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), # Chi Square Prior with SD 1
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) # Chi Square Prior with SD 1

#Level 1 = Poisson process, Level 2 = Zero inflation
#Zero inflated poisson model
load(file="Routput/MCMC_Ambient_Model_1.3.RData")

AM_model1.3 <- MCMCglmm(seed_pods ~ trait - 1 + plot,
                        random = ~idh(at.level(trait,1)):animal + idh(at.level(trait,2)):animal + 
                                  idh(at.level(trait,1)):matID + idh(at.level(trait,2)):matID, 
                        ginverse = list(animal=Ainv), 
                        rcov = ~idh(trait):units,
                        family = "zipoisson", data = ambient.total, prior=prior1.3, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T)

plot(AM_model1.3$VCV) # Fairly good trace plots
summary(AM_model1.3) # good effective sample size
autocorr.diag(AM_model1.3$VCV) # no autocorrelation
heidel.diag(AM_model1.3$VCV) # convergence success

#Model 1.3b Fisher Parameter Expanded Prior
#Note: we fit an overdispersed Poisson using the simple family='poisson' model. Apparently, Kruuk et al. 2014b had no problems running this model!
load(file="Routput/MCMC_Ambient_Model_1.3b.RData")

ambient.total <- MCMC.ambient
prior1.3b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model1.3b <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID, 
                         ginverse = list(animal=Ainv),
                         family = "poisson", data = ambient.total, prior=prior1.3b, 
                         nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects could not be estimated and have been removed. Use singular=OK to sample these effects, but use an informative prior!

plot(AM_model1.3b$VCV) # good trace plots
summary(AM_model1.3b) # good effective sample size
autocorr.diag(AM_model1.3b$VCV) # no autocorrelation
heidel.diag(AM_model1.3b$VCV) # convergence success

AM_herit1.3b <- AM_model1.3b$VCV[, "animal"]/(AM_model1.3b$VCV[, "animal"] + AM_model1.3b$VCV[, "matID"] 
                                              + AM_model1.3b$VCV[, "units"])
mean(AM_herit1.3b)
HPDinterval(AM_herit1.3b)

yhat_A1.3b <- predict(AM_model1.3b, type = "terms")
mu_A1.3b <- mean(AM_model1.3b[["Sol"]][ , "(Intercept)"])
va_A1.3b <- mean(AM_model1.3b[["VCV"]][ , "animal"])
vp_A1.3b <- mean(rowSums(AM_model1.3b[["VCV"]]))


QG_A1.3b <- QGparams(predict=yhat_A1.3b, mu=mu_A1.3b, var.a=va_A1.3b, var.p=vp_A1.3b, model = "Poisson.log")
QG_A1.3b
va_A1.3b / vp_A1.3b







#### No. 2 Fecundity of Flowering Plants - Univariate | Poisson ####
#Note: For Fecundity, we only include plants that survived to reach flowering


#Calculating Variance for Fecundity for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1 & flower==1) %>%
  summarise(mean=mean(seed_pods), sd=sd(seed_pods), var=(sd)^2) #var = 244.9818 *******
log(12.12949)

#Producing filtered data set for Fecundity
ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)

#Full model including Random: Va, Vm, and | Fixed: plot with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_2.1.RData")

ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior2.1 <- list(R = list(V = 1, nu = 0.002), 
                  G = list(G1 = list(V = 1, nu = 0.002),
                           G2 = list(V = 1, nu = 0.002), 
                           G3 = list(V = 1, nu = 0.002)))

AM_model2.1 <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal=Ainv),
                     family = "poisson", data = ambient.Fecundity, prior = prior2.1, #poission distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(AM_model2.1, file="Routput/MCMC_Ambient_Model_2.1.RData")

plot(AM_model2.1$VCV) # 
summary(AM_model2.1) # 
autocorr.diag(AM_model2.1$VCV) # 
heidel.diag(AM_model2.1$VCV) # 

#Model 2.1b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Ambient_Model_2.1b.RData")

ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior2.1b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model2.1b <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal=Ainv),
                     family = "poisson", data = ambient.Fecundity, prior = prior2.1b, #poission distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(AM_model2.1b, file="Routput/MCMC_Ambient_Model_2.1b.RData")

plot(AM_model2.1b$VCV) # 
summary(AM_model2.1b) # 
autocorr.diag(AM_model2.1b$VCV) # 
heidel.diag(AM_model2.1b$VCV) # 




#Model 2.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_2.2.RData")

ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior2.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model2.2 <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "poisson", data = ambient.Fecundity, prior = prior2.2,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(AM_model2.2$VCV) # poor trace plots
summary(AM_model2.2) # low effective sample size
autocorr.diag(AM_model2.2$VCV) # autocorrelation
heidel.diag(AM_model2.2$VCV) # convergence failing

#Model 2.2b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Ambient_Model_2.2b.RData")

ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior2.2b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model2.2b <- MCMCglmm(seed_pods ~ plot, random = ~animal + animalDom + matID,
                      ginverse = list(animal=Ainv, animalDom= Dinv),
                      family = "poisson", data = ambient.Fecundity, prior = prior2.2b,
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effets, but use an informative prior!

plot(AM_model2.2b$VCV) # good trace plots
summary(AM_model2.2b) # good effetive sample size
autocorr.diag(AM_model2.2b$VCV) # no autocorrelation
heidel.diag(AM_model2.2b$VCV) # convergence success

AM_herit2.2b <- AM_model2.2b$VCV[, "animal"]/(AM_model2.2b$VCV[, "animal"] + AM_model2.2b$VCV[, "animalDom"] + AM_model2.2b$VCV[, "matID"])
mean(AM_herit2.2b)
HPDinterval(AM_herit2.2b)

#Transformating data using "QGglmm" to get estimates of Va on data scale
yhat_A2.2b <- predict(AM_model2.2b, type = "terms")
mu_A2.2b <- mean(AM_model2.2b[["Sol"]][ , "(Intercept)"])
va_A2.2b <- mean(AM_model2.2b[["VCV"]][ , "animal"])
vp_A2.2b <- mean(rowSums(AM_model2.2b[["VCV"]]))

QG_A2.2b <- QGparams(predict=yhat_A2.2b, mu=mu_A2.2b, var.a=va_A2.2b, var.p=vp_A2.2b, model = "Poisson.log")
QG_A2.2b
va_A2.2b / vp_A2.2b






#Model 2.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_2.3.RData")

ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior2.3 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)))

AM_model2.3 <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "poisson", data = ambient.Fecundity, prior = prior2.3,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(AM_model2.3$VCV) # OK trace plots
summary(AM_model2.3) # low effective sample size
autocorr.diag(AM_model2.3$VCV) # autocorrelation
heidel.diag(AM_model2.3$VCV) # convergence success

#Model 2.3b attempted with Fisher Parameter Expanded Priors
load(file="Routput/MCMC_Ambient_Model_2.3b.RData")

ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior2.3b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model2.3b <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "poisson", data = ambient.Fecundity, prior = prior2.3b,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(AM_model2.3b$VCV) # good trace plots
summary(AM_model2.3b) # good effective sample size
autocorr.diag(AM_model2.3b$VCV) # no autocorrelation
heidel.diag(AM_model2.3b$VCV) # convergence success

AM_herit2.3b <- AM_model2.3b$VCV[, "animal"]/(AM_model2.3b$VCV[, "animal"] + AM_model2.3b$VCV[, "matID"] + AM_model2.3b$VCV[, "units"])
mean(AM_herit2.3b)
HPDinterval(AM_herit2.3b)


yhat_A2.3b <- predict(AM_model2.3b, type = "terms")
mu_A2.3b <- mean(AM_model2.3b[["Sol"]][ , "(Intercept)"])
va_A2.3b <- mean(AM_model2.3b[["VCV"]][ , "animal"])
vp_A2.3b <- mean(rowSums(AM_model2.3b[["VCV"]]))

QG_A2.3b <- QGparams(predict=yhat_A2.3b, mu=mu_A2.3b, var.a=va_A2.3b, var.p=vp_A2.3b, model = "Poisson.log")
QG_A2.3b
va_A2.3b / vp_A2.3b


#Model Comparison
AM_model2.2b$DIC
AM_model2.3b$DIC
















#### No. 3 Survival to Flowering - Univariate | Bernoulli ####
#Note: For this trait, we include all plants from the experiment
#We use an uninformative prior suggested by in Animal tutorial papers for a Bernoulli distribution, which follows a X^2 distribution with df = 1.
#We apply a parameter extension with alpha.mu and alpha.V to allow for the use of the X^2 prior dist. 
#We also apply a probit link using family="ordinal" (this is switched to "threshold".. see below)

#Producing filtered data set for Survival
ambient.Survival <- MCMC.ambient

#Probability of Survival to Flowering
ambient.Survival %>% 
  summarise(tot.flwr=sum(flower), n=n(), prob=tot.flwr/n)

#From Committee Report #2, Survival follows a Bernoulli distribution
#Note: According to Pierre de Villemereuil, adding additional random factors is not recommended. 
#Further, it is required to fix residual variance (Vr) to 1 and estimate Va solely.
#We use a prior suggested by Villemereuil et al. 2012 for binary traits
#We change the family from "ordinal" to "threshold" as suggested by a forum post in 2017 by Jarrod Hadfield (https://stat.ethz.ch/pipermail/r-sig-mixed-models/2017q4/026115.html), and apply a truncation with trunc=TRUE

#Full model including Random: Va, Vm, and | Fixed: plot with Chi Sq df=1 prior
#Note, since autocorrelation is strong in MCMC w/ binary data, may need to multiply 'thin' and 'nitt' by 10.
load(file="Routput/MCMC_Ambient_Model_3.1.RData")

ambient.Survival <- MCMC.ambient
prior3.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model3.1 <- MCMCglmm(flower ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal=Ainv),
                     family = "threshold", data = ambient.Survival, prior = prior3.1, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(AM_model3.1, file="Routput/MCMC_Ambient_Model_3.1.RData")

plot(AM_model3.1$VCV) # 
summary(AM_model3.1) #
autocorr.diag(AM_model3.1$VCV) # 
heidel.diag(AM_model3.1$VCV) # 


#Model 3.2 with additive, dominance and maternal effects + plot (fixed) using Chi Sq df=1 Prior
load(file="Routput/MCMC_Ambient_Model_3.2.RData")

ambient.Survival <- MCMC.ambient
prior3.2 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model3.2 <- MCMCglmm(flower ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal=Ainv, animalDom = Dinv),
                     family = "threshold", data = ambient.Survival, prior = prior3.2, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)

plot(AM_model3.2$VCV) # good trace plots
summary(AM_model3.2) # good effective sample size
autocorr.diag(AM_model3.2$VCV) # no autocorrelation
heidel.diag(AM_model3.2$VCV) # convergence success

AM_herit3.2 <- AM_model3.2$VCV[, "animal"]/(AM_model3.2$VCV[, "animal"] + AM_model3.2$VCV[, "animalDom"] + AM_model3.2$VCV[, "matID"] + 1)
mean(AM_herit3.2)
HPDinterval(AM_herit3.2)

yhat_A3.2 <- predict(AM_model3.2, type = "terms")
mu_A3.2 <- mean(AM_model3.2[["Sol"]][ , "(Intercept)"])
va_A3.2 <- mean(AM_model3.2[["VCV"]][ , "animal"])
vp_A3.2 <- mean(rowSums(AM_model3.2[["VCV"]]))

QG_A3.2 <- QGparams(predict=yhat_A3.2, mu=mu_A3.2, var.a=va_A3.2, var.p=vp_A3.2, model = "binom1.probit")
QG_A3.2
va_A3.2 / vp_A3.2


#Model 3.3 with animal + maternal effects
load(file="Routput/MCMC_Ambient_Model_3.3.RData")

ambient.Survival <- MCMC.ambient
prior3.3 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model3.3 <- MCMCglmm(flower ~ plot, random = ~animal + matID, 
                     ginverse = list(animal=Ainv),
                     family = "threshold", data = ambient.Survival, prior = prior3.3, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)

#Some fixed effects were not estimable and removed. Use an informative prior.

plot(AM_model3.3$VCV) # good trace plot, although density seems bimodal
summary(AM_model3.3) # good effective sample size
autocorr.diag(AM_model3.3$VCV) # no autocorrelation
heidel.diag(AM_model3.3$VCV) # convergence success

AM_herit3.3 <- AM_model3.3$VCV[, "animal"]/(AM_model3.3$VCV[, "animal"] + AM_model3.3$VCV[, "matID"] + 1)
mean(AM_herit3.3)
HPDinterval(AM_herit3.3)

yhat_A3.3 <- predict(AM_model3.3, type = "terms")
mu_A3.3 <- mean(AM_model3.3[["Sol"]][ , "(Intercept)"])
va_A3.3 <- mean(AM_model3.3[["VCV"]][ , "animal"])
vp_A3.3 <- mean(rowSums(AM_model3.3[["VCV"]]))

QG_A3.3 <- QGparams(predict=yhat_A3.3, mu=mu_A3.3, var.a=va_A3.3, var.p=vp_A3.3, model = "binom1.probit")
QG_A3.3
va_A3.3 / vp_A3.3

#Model Comparison
AM_model3.2$DIC
AM_model3.3$DIC















#### No. 4 Survival to Flowering + Fecundity - Bivariate ####
#Note: For this trait, we include all plants from the experiment
#Survival to Flowering = Bernoulli Dist | Fecundity = Poisson Dist


#Calculating Variance for Fecundity for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  summarise(mean(seed_pods), sd=sd(seed_pods), var=(sd)^2) #mean = 4.985826
log(4.985826) # = 1.606599

MCMC.ambient %>%
  summarise(sum.f=sum(flower), n=n(), mean=sum.f/n) # mean prob = 0.41105
logit(0.41105*(1-0.41105)) # -1.141267

#Subset data
ambient.Survival <- MCMC.ambient

#Bivariate Full Model
load(file="Routput/MCMC_Ambient_Model_4.2.RData")

ambient.Survival <- MCMC.ambient
#prior4.2 <- list(R = list(V = diag(2), nu = 0.002),
                 #G = list(G1 = list(V = diag(2), nu = 1, alpha.mu = c(0,0), alpha.V = diag(2)*0.1),
                          #G2 = list(V = diag(2), nu = 1, alpha.mu = c(0,0), alpha.V = diag(2)*0.1), 
                          #G3 = list(V = diag(2), nu = 1, alpha.mu = c(0,0), alpha.V = diag(2)*0.1)))

#This is the prior I used for the plasticity models
#It's used for a binary + other distribution dataset
#Note that for binary distributions, the residual must be fixed to 1.
prior4.2 <- list(R = list(V = diag(2), nu=2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G3 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

AM_model4.2 <- MCMCglmm(cbind(flower, seed_pods) ~ trait + plot - 1,
                        random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                        ginverse = list(animal=Ainv, animalDom = Dinv), rcov = ~us(trait):units,
                        family = c("threshold", "poisson"), 
                        data = ambient.Survival, prior = prior4.2, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(AM_model4.2, file="Routput/MCMC_Ambient_Model_4.2.RData")

plot(AM_model4.2$VCV) #
summary(AM_model4.2) #
autocorr.diag(AM_model4.2$VCV) #
heidel.diag(AM_model4.2$VCV) #



#Model 4.2b --- Did not run because expanded priors used for model 4.2 anyway
#We use a Chi square distribution with df=1 prior for the binary trait (flower) and a parameter expanded prior for the Poisson process.
#Note that, it is unclear what to set nu for the Poisson process when the binary process requires nu = 1000
#We also use adivce from 5 years ago by Hadfield (see https://r-sig-mixed-models.r-project.narkive.com/Bi7GhoZ8/r-sig-me-mcmcglmm-bivariate-model-and-prior) to set the residual belief to 0.
load(file="Routput/MCMC_Ambient_Model_4.2.RData")

ambient.Survival <- MCMC.ambient
prior4.2b <- list(R = list(V = diag(2), nu = 0, fix = 2),
                 G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*100),
                          G2 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*100), 
                          G3 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*100)))

AM_model4.2b <- MCMCglmm(cbind(flower, seed_pods) ~ trait + plot - 1,
                        random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                        ginverse = list( animalDom = Dinv), rcov = ~us(trait):units,
                        family = c("threshold", "poisson"), 
                        data = ambient.Survival, prior = prior4.2b, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(AM_model4.2b, file="Routput/MCMC_Ambient_Model_4.2b.RData")

plot(AM_model4.2b$VCV) #
summary(AM_model4.2b) #
autocorr.diag(AM_model4.2b$VCV) #
heidel.diag(AM_model4.2b$VCV) #


#Model 4.3 without dominance variance
load(file="Routput/MCMC_Ambient_Model_4.3.RData")

ambient.Survival <- MCMC.ambient
prior4.3 <- list(R = list(V = diag(2), nu=2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

AM_model4.3 <- MCMCglmm(cbind(flower, seed_pods) ~ trait + plot - 1,
                        random = ~us(trait):animal + us(trait):matID,
                        ginverse = list(animal=Ainv), rcov = ~us(trait):units,
                        family = c("threshold", "poisson"), 
                        data = ambient.Survival, prior = prior4.3, 
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(AM_model4.3$VCV) # good trace plots (except matID)
summary(AM_model4.3) # good effective smaple size
autocorr.diag(AM_model4.3$VCV) # no autocorrelation
heidel.diag(AM_model4.3$VCV) # convergence issues for matID.. but not interested in it so its OK

#Calculating Genetic Correlations = Covariance b/w trait 1 & 2 / sqrt (var [trait 1] * var [trait 2])

# (1) Survival & Fecundity
gen.corrA4.3 <-AM_model4.3$VCV[,'traitflower:traitseed_pods.animal']/
  sqrt(AM_model4.3$VCV[,'traitflower:traitflower.animal']*AM_model4.3$VCV[,'traitseed_pods:traitseed_pods.animal']) 
mean(gen.corrA4.3) #Post Mean = 0.7317867
posterior.mode(gen.corrA4.3) #Post Mode = 0.729908
HPDinterval(gen.corrA4.3) # Posterior 95% CI = (0.7060835, 0.7554622)

phen.corrA4.3 <- cor.test(ambient.Survival$flower, ambient.Survival$seed_pods, test="pearson")
phen.corrA4.3 # Corr = 0.5112517 | Posterior 95% CI =  (0.4861976, 0.5354659)


















#### No. 5 Flowering Success - Univariate | Bernoulli ####
#Note: For this trait, we include only include plants that germinated
#We use an uninformative prior suggested by in Animal tutorial papers for a Bernoulli distribution, which follows a X^2 distribution with df = 1.
#We apply a parameter extension with alpha.mu and alpha.V to allow for the use of the X^2 prior dist. 
#We also apply a probit link using family="ordinal" (this is switched to "threshold".. see below)

#Producing filtered data set for Fecundity
ambient.Flowering <- MCMC.ambient %>% filter(germ==1)

#From Committee Report #2, Survival follows a Bernoulli distribution

#Full model including Random: Va, Vm, and | Fixed: plot with Chi Sq df=1 prior
#Note, since autocorrelation is strong in MCMC w/ binary data, may need to multiply 'thin' and 'nitt' by 10.
load(file="Routput/MCMC_Ambient_Model_5.1.RData")

ambient.Flowering <- MCMC.ambient %>% filter(germ==1)
prior5.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model5.1 <- MCMCglmm(flower ~ plot, random = ~animal + matID +, 
                     ginverse = list(animal=Ainv),
                     family = "threshold", data = ambient.Flowering, prior = prior5.1, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(AM_model5.1, file="Routput/MCMC_Ambient_Model_5.1.RData")

plot(AM_model5.1$VCV) # 
summary(AM_model5.1) # 
autocorr.diag(AM_model5.1$VCV) # 
heidel.diag(AM_model5.1$VCV) #


#Model with additive, dominance and maternal effects + plot (fixed) using Chi Sq df=1 Prior
load(file="Routput/MCMC_Ambient_Model_5.2.RData")

ambient.Flowering <- MCMC.ambient %>% filter(germ==1)
prior5.2 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model5.2 <- MCMCglmm(flower ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal=Ainv, animalDom = Dinv),
                     family = "threshold", data = ambient.Flowering, prior = prior5.2, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Warning: Some fixed effects not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(AM_model5.2$VCV) # good trace plots, but close to 0
summary(AM_model5.2) # good effective sample size
autocorr.diag(AM_model5.2$VCV) # no autocorrelation
heidel.diag(AM_model5.2$VCV) # convergence success

AM_herit5.2 <- AM_model5.2$VCV[, "animal"]/(AM_model5.2$VCV[, "animal"] + AM_model5.2$VCV[, "animalDom"] + AM_model5.2$VCV[, "matID"] + 1)
mean(AM_herit5.2)
HPDinterval(AM_herit5.2)


yhat_A5.2 <- predict(AM_model5.2, type = "terms")
mu_A5.2 <- mean(AM_model5.2[["Sol"]][ , "(Intercept)"])
va_A5.2 <- mean(AM_model5.2[["VCV"]][ , "animal"])
vp_A5.2 <- mean(rowSums(AM_model5.2[["VCV"]]))

QG_A5.2 <- QGparams(predict=yhat_A5.2, mu=mu_A5.2, var.a=va_A5.2, var.p=vp_A5.2, model = "binom1.probit")
QG_A5.2
va_A5.2 / vp_A5.2




#Model with additive and maternal effects + fixed plot using a Chi Sq df=1 prior
load(file="Routput/MCMC_Ambient_Model_5.3.RData")

ambient.Flowering <- MCMC.ambient %>% filter(germ==1)
prior5.3 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model5.3 <- MCMCglmm(flower ~ plot, random = ~animal + matID, 
                     ginverse = list(animal=Ainv),
                     family = "threshold", data = ambient.Flowering, prior = prior5.3, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects were not estimable and removed. Use an informative prior.

plot(AM_model5.3$VCV) # OK trace plot, not excellent
summary(AM_model5.3) # good effective sample size
autocorr.diag(AM_model5.3$VCV) # no autocorrelation
heidel.diag(AM_model5.3$VCV) # convergence success

AM_herit5.3 <- AM_model5.3$VCV[, "animal"]/(AM_model5.3$VCV[, "animal"] + AM_model5.3$VCV[, "matID"] + 1)
mean(AM_herit5.3)
HPDinterval(AM_herit5.3)

yhat_A5.3 <- predict(AM_model5.3, type = "terms")
mu_A5.3 <- mean(AM_model5.3[["Sol"]][ , "(Intercept)"])
va_A5.3 <- mean(AM_model5.3[["VCV"]][ , "animal"])
vp_A5.3 <- mean(rowSums(AM_model5.3[["VCV"]]))

QG_A5.3 <- QGparams(predict=yhat_A5.3, mu=mu_A5.3, var.a=va_A5.3, var.p=vp_A5.3, model = "binom1.probit")
QG_A5.3
va_A5.3 / vp_A5.3


#Model Comparison
AM_model5.2$DIC
AM_model5.3$DIC














#### No. 6 Germination Success - Univariate | Bernoulli ####
#Note: For this trait, we include all plants from the study
#We use an uninformative prior suggested by in Animal tutorial papers for a Bernoulli distribution, which follows a X^2 distribution with df = 1.
#We apply a parameter extension with alpha.mu and alpha.V to allow for the use of the X^2 prior dist. 
#We also apply a probit link using family = "threshold"

#Producing filtered data set for Germination Success
ambient.Germination <- MCMC.ambient

#From Committee Report #2, Germination Success follows a Bernoulli distribution

#Full model including Random: Va, Vm, and | Fixed: plot with Chi Sq df=1 prior
#Note, since autocorrelation is strong in MCMC w/ binary data, may need to multiply 'thin' and 'nitt' by 10.
load(file="Routput/MCMC_Ambient_Model_6.1.RData")

ambient.Germination <- MCMC.ambient
prior6.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model6.1 <- MCMCglmm(germ ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal=Ainv),
                     family = "threshold", data = ambient.Germination, prior = prior6.1, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(AM_model6.1, file="Routput/MCMC_Ambient_Model_6.1.RData")

plot(AM_model6.1$VCV) #
summary(AM_model6.1) #
autocorr.diag(AM_model6.1$VCV) #
heidel.diag(AM_model6.1$VCV) #

#Model with additive, dominance and maternal effects + plot (fixed) using Chi Sq df=1 Prior
load(file="Routput/MCMC_Ambient_Model_6.2.RData")

prior6.2 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model6.2 <- MCMCglmm(germ ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal=Ainv, animalDom = Dinv),
                     family = "threshold", data = ambient.Germination, prior = prior6.2, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects not estimable and have been removed, use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(AM_model6.2$VCV) # trace plots look OK
summary(AM_model6.2) # good effective sample size
autocorr.diag(AM_model6.2$VCV) # no autocorrelation
heidel.diag(AM_model6.2$VCV) # convergence success

AM_herit6.2 <- AM_model6.2$VCV[, "animal"]/(AM_model6.2$VCV[, "animal"] + AM_model6.2$VCV[, "animalDom"] + AM_model6.2$VCV[, "matID"] + 1)
mean(AM_herit6.2)
HPDinterval(AM_herit6.2)

yhat_A6.2 <- predict(AM_model6.2, type = "terms")
mu_A6.2 <- mean(AM_model6.2[["Sol"]][ , "(Intercept)"])
va_A6.2 <- mean(AM_model6.2[["VCV"]][ , "animal"])
vp_A6.2 <- mean(rowSums(AM_model6.2[["VCV"]]))

QG_A6.2 <- QGparams(predict=yhat_A6.2, mu=mu_A6.2, var.a=va_A6.2, var.p=vp_A6.2, model = "binom1.probit")
QG_A6.2
va_A6.2 / vp_A6.2



#Model with additive and maternal effects + fixed plot using a Chi Sq df=1 prior
load(file="Routput/MCMC_Ambient_Model_6.3.RData")

ambient.Germination <- MCMC.ambient
prior6.3 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model6.3 <- MCMCglmm(germ ~ plot, random = ~animal + matID, 
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = ambient.Germination, prior = prior6.3, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects were not estimable and removed. Use an informative prior.

plot(AM_model6.3$VCV) # OK trace plot
summary(AM_model6.3) # good effective sample size
autocorr.diag(AM_model6.3$VCV) # no autocorrelation
heidel.diag(AM_model6.3$VCV) # convergence success

AM_herit6.3 <- AM_model6.3$VCV[, "animal"]/(AM_model6.3$VCV[, "animal"] + AM_model6.3$VCV[, "matID"] + 1)
mean(AM_herit6.3)
HPDinterval(AM_herit6.3)

yhat_A6.3 <- predict(AM_model6.3, type = "terms")
mu_A6.3 <- mean(AM_model6.3[["Sol"]][ , "(Intercept)"])
va_A6.3 <- mean(AM_model6.3[["VCV"]][ , "animal"])
vp_A6.3 <- mean(rowSums(AM_model6.3[["VCV"]]))

QG_A6.3 <- QGparams(predict=yhat_A6.3, mu=mu_A6.3, var.a=va_A6.3, var.p=vp_A6.3, model = "binom1.probit")
QG_A6.3
va_A6.3 / vp_A6.3

#Model Comparison
AM_model6.2$DIC
AM_model6.3$DIC














#### No. 7 Seed Maturation Success - Univariate | Bernoulli ####
#Note: For this trait, we include only include plants that germinated and flowered
#We use an uninformative prior suggested by in Animal tutorial papers for a Bernoulli distribution, which follows a X^2 distribution with df = 1.
#We apply a parameter extension with alpha.mu and alpha.V to allow for the use of the X^2 prior dist. 
#We also apply a probit link using family = "threshold"


#Producing filtered data set for Seed Maturation Success
ambient.SeedMaturation <- MCMC.ambient %>% filter(germ==1, flower==1)

#From Committee Report #2, Seed Maturaiton Success follows a Bernoulli distribution

#Full model including Random: Va, Vm, and | Fixed: plot with Chi Sq df=1 prior
#Note, since autocorrelation is strong in MCMC w/ binary data, may need to multiply 'thin' and 'nitt' by 10.
load(file="Routput/MCMC_Ambient_Model_7.1.RData")

ambient.SeedMaturation <- MCMC.ambient %>% filter(germ==1, flower==1)
prior7.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model7.1 <- MCMCglmm(seed ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal=Ainv),
                     family = "threshold", data = ambient.SeedMaturation, prior = prior7.1, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(AM_model7.1, file="Routput/MCMC_Ambient_Model_7.1.RData")

plot(AM_model7.1$VCV) # 
summary(AM_model7.1) # 
autocorr.diag(AM_model7.1$VCV) # 
heidel.diag(AM_model7.1$VCV) #

#Model with additive, dominance and maternal effects + plot (fixed) using Chi Sq df=1 Prior
load(file="Routput/MCMC_Ambient_Model_7.2.RData")

ambient.SeedMaturation <- MCMC.ambient %>% filter(germ==1, flower==1)
prior7.2 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model7.2 <- MCMCglmm(seed ~ plot, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = ambient.SeedMaturation, prior = prior7.2, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects not estimable and have been removed, use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(AM_model7.2$VCV) # good trace plots
summary(AM_model7.2) # good effective sample size
autocorr.diag(AM_model7.2$VCV) # no autocorrelation
heidel.diag(AM_model7.2$VCV) # convergence success

AM_herit7.2 <- AM_model7.2$VCV[, "animal"]/(AM_model7.2$VCV[, "animal"] + AM_model7.2$VCV[, "animalDom"] + AM_model7.2$VCV[, "matID"] + 1)
mean(AM_herit7.2)
HPDinterval(AM_herit7.2)

yhat_A7.2 <- predict(AM_model7.2, type = "terms")
mu_A7.2 <- mean(AM_model7.2[["Sol"]][ , "(Intercept)"])
va_A7.2 <- mean(AM_model7.2[["VCV"]][ , "animal"])
vp_A7.2 <- mean(rowSums(AM_model7.2[["VCV"]]))

QG_A7.2 <- QGparams(predict=yhat_A7.2, mu=mu_A7.2, var.a=va_A7.2, var.p=vp_A7.2, model = "binom1.probit")
QG_A7.2
va_A7.2 / vp_A7.2



#Model with additive and maternal effects + fixed plot using a Chi Sq df=1 prior
load(file="Routput/MCMC_Ambient_Model_7.3.RData")

ambient.SeedMaturation <- MCMC.ambient %>% filter(germ==1, flower==1)
prior7.3 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_model7.3 <- MCMCglmm(seed ~ plot, random = ~animal + matID, 
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = ambient.SeedMaturation, prior = prior7.3, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#Some fixed effects were not estimable and removed. Use an informative prior.

plot(AM_model7.3$VCV) # OK trace plots, not great
summary(AM_model7.3) # good effective sample size
autocorr.diag(AM_model7.3$VCV) # no autocorrelation
heidel.diag(AM_model7.3$VCV) # convergence success for all

AM_herit7.3 <- AM_model7.3$VCV[, "animal"]/(AM_model7.3$VCV[, "animal"] + AM_model7.3$VCV[, "matID"] + 1)
mean(AM_herit7.3)
HPDinterval(AM_herit7.3)

yhat_A7.3 <- predict(AM_model7.3, type = "terms")
mu_A7.3 <- mean(AM_model7.3[["Sol"]][ , "(Intercept)"])
va_A7.3 <- mean(AM_model7.3[["VCV"]][ , "animal"])
vp_A7.3 <- mean(rowSums(AM_model7.3[["VCV"]]))

QG_A7.3 <- QGparams(predict=yhat_A7.3, mu=mu_A7.3, var.a=va_A7.3, var.p=vp_A7.3, model = "binom1.probit")
QG_A7.3
va_A7.3 / vp_A7.3

#Model Comparison
AM_model7.2$DIC
AM_model7.3$DIC













#### No. 8 Leaf Number - Univariate | Gaussian ####
#Note: For leaf number, we only include plants that germinated
#We use an uninformative prior that divides the total phenotypic variance by the number of random effects included
#Variance in Total Fitness = squared standard deviation
#We treat plot as a random effect as from previous analysis (01_Data_exploration), we observe that there is a plot effect on this trait... There is likely some sampling error caused by plot due to the placement of the physical plots at KSR with respect to the tree line, etc so we want to give less statistical weight to the effect of plot. We also want to estimate the effect of plot 


#Calculating Variance for leaf number for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1) %>%
  summarise(sd=sd(leaf), var=(sd)^2) #var = 5.709704 *******

#Producing filtered data set for leaf number
ambient.leaf <- MCMC.ambient %>% filter(germ==1)

#From Committee Report #2: Leaf Number follows a Gaussian distribution

#Full model including Random: Va, Vm, and | Fixed: plot with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_8.1.RData")

ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior8.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model8.1 <- MCMCglmm(leaf ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = ambient.leaf, prior = prior8.1, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(AM_model8.1, file="Routput/MCMC_Ambient_Model_8.1.RData")

plot(AM_model8.1$VCV) #  
summary(AM_model8.1) # 
autocorr.diag(AM_model8.1$VCV) # 
heidel.diag(AM_model8.1$VCV) # 


#Model 2.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_8.2.RData")

ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior8.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model8.2 <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.leaf, prior = prior8.2, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#AM_model8.2 used to have equal variance priors - redoing with Inverse Gamma Priors

plot(AM_model8.2$VCV) # poor trace plots
summary(AM_model8.2) # low effective sample size
autocorr.diag(AM_model8.2$VCV) # autocorrelation
heidel.diag(AM_model8.2$VCV) # convergence fail

#Model 8.2b with a Fisher Expanded Parameter Prior
load(file="Routput/MCMC_Ambient_Model_8.2b.RData")

ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior8.2b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model8.2b <- MCMCglmm(leaf ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.leaf, prior = prior8.2b, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and have been removed, use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(AM_model8.2b$VCV) # good trace plots
summary(AM_model8.2b) # good effective size
autocorr.diag(AM_model8.2b$VCV) # no autocorrelation
heidel.diag(AM_model8.2b$VCV) # convergence success

AM_herit8.2b <- AM_model8.2b$VCV[, "animal"]/(AM_model8.2b$VCV[, "animal"] + AM_model8.2b$VCV[, "animalDom"] + AM_model8.2b$VCV[, "matID"] + AM_model8.2b$VCV[, "units"])
mean(AM_herit8.2b)
HPDinterval(AM_herit8.2b)

#Model 8.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_8.3.RData")

ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior8.3 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)))

AM_model8.3 <- MCMCglmm(leaf ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = ambient.leaf, prior = prior8.3, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(AM_model8.3$VCV) # good trace plot, except maternal effects close to 0
summary(AM_model8.3) # decent effective sample size
autocorr.diag(AM_model8.3$VCV) # no autocorrelation
heidel.diag(AM_model8.3$VCV) # convergence success

#We report values from model 8.3, NOT 8.3b
AM_herit8.3 <- AM_model8.3$VCV[, "animal"]/(AM_model8.3$VCV[, "animal"] + AM_model8.3$VCV[, "matID"] + AM_model8.3$VCV[, "units"])
mean(AM_herit8.3)
HPDinterval(AM_herit8.3)


#Model 8.3b using Fisher Parameter Expanded Priors
load(file="Routput/MCMC_Ambient_Model_8.3b.RData")

ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior8.3b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model8.3b <- MCMCglmm(leaf ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = ambient.leaf, prior = prior8.3b, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(AM_model8.3b$VCV) # good trace plots (almost identical to 8.3)
summary(AM_model8.3b) # good effective sample size
autocorr.diag(AM_model8.3b$VCV) # no autocorrelation
heidel.diag(AM_model8.3b$VCV) # convergence success for all

#Model Comparison
AM_model8.2b$DIC
AM_model8.3b$DIC













#### No. 9 Height - Univariate | Gaussian ####
#Note: For height, we only include plants that germinated
#We use an uninformative prior that divides the total phenotypic variance by the number of random effects included
#Variance in Total Fitness = squared standard deviation
#We treat plot as a random effect as from previous analysis (01_Data_exploration), we observe that there is a plot effect on this trait... There is likely some sampling error caused by plot due to the placement of the physical plots at KSR with respect to the tree line, etc so we want to give less statistical weight to the effect of plot. We also want to estimate the effect of plot 

#Calculating Variance for height for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1, !height==0) %>%
  summarise(sd=sd(height), var=(sd)^2) #var = 147.7686 *******

#Producing filtered data set for height
ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)

#From Committee Report #2: height follows a Gaussian distribution

#Full model including Random: Va, Vm, and | Fixed: plot with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_9.1.RData")

ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior9.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model9.1 <- MCMCglmm(height ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = ambient.height, prior = prior9.1, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(AM_model9.1, file="Routput/MCMC_Ambient_Model_9.1.RData")

plot(AM_model9.1$VCV) #
summary(AM_model9.1) #
autocorr.diag(AM_model9.1$VCV) #
heidel.diag(AM_model9.1$VCV) # 


#Model 9.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_9.2.RData")

ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior9.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model9.2 <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.height, prior = prior9.2, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and have been removed, use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(AM_model9.2$VCV) # poor trace plots
summary(AM_model9.2) # poor effective sample size
autocorr.diag(AM_model9.2$VCV) # high autocorrelation
heidel.diag(AM_model9.2$VCV) # convergence fail

#Model 9.2b with a Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Ambient_Model_9.2b.RData")

ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior9.2b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model9.2b <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.height, prior = prior9.2b, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and have been removed, use singular.ok=TRUE to sample these effects, but use an informative prior!

plot(AM_model9.2b$VCV) # good trace plots
summary(AM_model9.2b) # high effective sample size
autocorr.diag(AM_model9.2b$VCV) # no autocorrelation
heidel.diag(AM_model9.2b$VCV) # convergence success

AM_herit9.2b <- AM_model9.2b$VCV[, "animal"]/(AM_model9.2b$VCV[, "animal"] + AM_model9.2b$VCV[, "animalDom"] + AM_model9.2b$VCV[, "matID"] + AM_model9.2b$VCV[, "units"])
mean(AM_herit9.2b)
HPDinterval(AM_herit9.2b)

library(lmodel2)
ID.est <- AM_model9.2b$VCV[, "animal"][1:dim(AM_model9.2b$VCV)[1]]
IDD.est <- AM_model9.2b$VCV[, "animalDom"][1:dim(AM_model9.2b$VCV)[1]]
mareg <- lmodel2(animal.est~animalDom.est)
x11(w = 8, h = 8)
plot(mareg, method = "MA", xlim = c(0,0.7), ylim = c(0,0.7),
       + xlab = "IDD estimates", ylab = "ID estimates",
       + main = paste("Sampling correlation: ", round(mareg$r, 3), sep =""),
       + sub = "Line represents the major axis regression") 

#Model 9.2c with a Fisher Parameter Expanded Prior and Maternal Group 2 ONLY
#Testing if we removing maternal group 2, we can get a stronger signal from 1st cousin relationships to estimate dominance
load(file="Routput/MCMC_Ambient_Model_9.2c.RData")

ambient.height <- MCMC.ambient %>% filter(germ==1 & !height==0 & mat_group=="2")
prior9.2c <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model9.2c <- MCMCglmm(height ~ plot, random = ~animal + animalDom + matID,
                         ginverse = list(animal = Ainv, animalDom= Dinv),
                         family = "gaussian", data = ambient.height, prior = prior9.2c, #gaussian distribution
                         nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(AM_model9.2b$VCV) # 
summary(AM_model9.2b) # 
autocorr.diag(AM_model9.2b$VCV) # 
heidel.diag(AM_model9.2b$VCV) # 


#Model 9.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_9.3.RData")

ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior9.3 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)))

AM_model9.3 <- MCMCglmm(height ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = ambient.height, prior = prior9.3, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(AM_model9.3$VCV) # good trace plot, with the exception of maternal effects
summary(AM_model9.3) # decent effective sample size
autocorr.diag(AM_model9.3$VCV) # some autocorrelation for additive effects
heidel.diag(AM_model9.3$VCV) # convergence fail for maternal effects

#Model 9.3b using Fisher Parameter Expanded Priors
load(file="Routput/MCMC_Ambient_Model_9.3b.RData")

ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior9.3b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model9.3b <- MCMCglmm(height ~ plot, random = ~animal + matID,
                      ginverse = list(animal = Ainv),
                      family = "gaussian", data = ambient.height, prior = prior9.3b, #gaussian distribution
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimable and removed. Use an informative prior!

plot(AM_model9.3b$VCV) # good trace plots
summary(AM_model9.3b) # good effective sample size
autocorr.diag(AM_model9.3b$VCV) # no autocorrelation
heidel.diag(AM_model9.3b$VCV) # convergence success
#WE REPORT VALUES FROM THIS MODEL

AM_herit9.3b <- AM_model9.3b$VCV[, "animal"]/(AM_model9.3b$VCV[, "animal"] + AM_model9.3b$VCV[, "matID"] + AM_model9.3b$VCV[, "units"])
mean(AM_herit9.3b)
HPDinterval(AM_herit9.3b)


#Model Comparison
AM_model9.2b$DIC
AM_model9.3b$DIC














#### No. 10 Height*Leaf*Seed_Pod - Multivariate | Gaus,Gaus,Pois ####
#Note: For this G matrix, we only use plants that germinated

#Producing filtered data set for height
ambient.trivar <- MCMC.ambient %>% filter(germ==1, !height==0)

#Model 10.2 with 3 additive, dominance and maternal effects + fixed plot effects
#In the Punetes et al. (2016) paper, they use a prior with nu=4.001
prior10.2 <- list(R=list(V=diag(3), nu=2.002), 
               G = list(G1 = list(V=diag(3), nu = 2.002),
                        G2 = list(V=diag(3), nu = 2.002),
                        G3 = list(V=diag(3), nu = 2.002)))

AM_model10.2 <- MCMCglmm(cbind(height, leaf, seed_pods) ~ trait + plot - 1,
                      random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                      ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~us(trait):units,
                      family = c("gaussian","gaussian","poisson"), 
                      data = ambient.trivar, prior = prior10.2,
                      nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(AM_model10.2, file="Routput/MCMC_Ambient_Model_10.2.RData")

plot(AM_model10.2$VCV) #
summary(AM_model10.2) #
autocorr.diag(AM_model10.2$VCV) #



#Model 10.2b with Fisher Parameter Expansion
load(file="Routput/MCMC_Ambient_Model_10.2b.RData")

ambient.trivar <- MCMC.ambient %>% filter(germ==1, !height==0)
prior10.2b <- list(R=list(V=diag(3), nu=2.002), 
                   G = list(G1 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000),
                            G2 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000),
                            G3 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000)))


AM_model10.2b <-  MCMCglmm(cbind(height, leaf, seed_pods) ~ trait + plot - 1,
           random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~us(trait):units,
           family = c("gaussian","gaussian","poisson"), 
           data = ambient.trivar, prior = prior10.2b,
           nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(AM_model10.2b, file="Routput/MCMC_Ambient_Model_10.2b.RData")

plot(AM_model10.2b$VCV) #
summary(AM_model10.2b) #
autocorr.diag(AM_model10.2b$VCV) #


#Model 10.3 with 3 Additive and Maternal Effects + plot
load(file="Routput/MCMC_Ambient_Model_10.3.RData")

ambient.trivar <- MCMC.ambient %>% filter(germ==1, !height==0)
prior10.3 <- list(R=list(V=diag(3), nu=2.002), 
                  G = list(G1 = list(V=diag(3), nu = 2.002),
                           G2 = list(V=diag(3), nu = 2.002)))

AM_model10.3 <- MCMCglmm(cbind(height, leaf, seed_pods) ~ trait + plot - 1,
                      random = ~us(trait):animal + us(trait):matID,
                      ginverse = list(animal=Ainv), rcov=~us(trait):units,
                      family = c("gaussian","gaussian","poisson"), 
                      data = ambient.trivar, prior = prior10.3,
                      nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

#Some fixed effects not estimable and removed. Use an informative prior!

plot(AM_model10.3$VCV) # good trace plots!
summary(AM_model10.3) # fairly good effective sample sizes
autocorr.diag(AM_model10.3$VCV) # no autocorrelation
autocorr.plot(AM_model10.3$VCV)
heidel.diag(AM_model10.3$VCV) #convergence success for all

#We use this model for the multivariate (10.3)

#Calculating Genetic Correlations = Covariance b/w trait 1 & 2 / sqrt (var [trait 1] * var [trait 2])
#Trait Combinations: 
  # (1) Height & Leaf Number
  # (2) Height & Seed Pod Number
  # (3) Leaf & Seed Pod NUmber

# (1) Height & Leaf Number
gen.corrA10.3_1 <-AM_model10.3$VCV[,'traitheight:traitleaf.animal']/
  sqrt(AM_model10.3$VCV[,'traitheight:traitheight.animal']*AM_model10.3$VCV[,'traitleaf:traitleaf.animal']) 
mean(gen.corrA10.3_1) #Post Mean = 0.783492
posterior.mode(gen.corrA10.3_1) #Post Mode = 0.8223059
HPDinterval(gen.corrA10.3_1) # Posterior 95% CI = (0.6145835, 0.9267772)

phen.corrA10.3_1 <- cor.test(ambient.trivar$height, ambient.trivar$leaf, test="pearson")
phen.corrA10.3_1 # Corr = 0.7735388 | Posterior 95% CI =  (0.7519433, 0.7934759)


# (2) Height & Seed Pod Number
gen.corrA10.3_2 <-AM_model10.3$VCV[,'traitheight:traitseed_pods.animal']/
  sqrt(AM_model10.3$VCV[,'traitheight:traitheight.animal']*AM_model10.3$VCV[,'traitseed_pods:traitseed_pods.animal']) 
mean(gen.corrA10.3_2) #Post Mean = 0.6083324
posterior.mode(gen.corrA10.3_2) #Post Mode = 0.6568689
HPDinterval(gen.corrA10.3_2) # Posterior 95% CI = (0.3465297, 0.8303014)

phen.corrA10.3_2 <- cor.test(ambient.trivar$height, ambient.trivar$seed_pods, test="pearson")
phen.corrA10.3_2 # Corr = 0.6101868 | Poster 95% CI = (0.5767315, 0.6415988)


# (3) Leaf Number and Seed Pod Number
gen.corrA10.3_3 <-AM_model10.3$VCV[,'traitseed_pods:traitleaf.animal']/
  sqrt(AM_modelA10.3$VCV[,'traitseed_pods:traitseed_pods.animal']*AM_model10.3$VCV[,'traitleaf:traitleaf.animal']) 
mean(gen.corrA10.3_3) #Post Mean = 0.5551225
posterior.mode(gen.corrA10.3_3) #Post Mode = 0.6065478
HPDinterval(gen.corrA10.3_3) #Posterior 95% CI = (0.2630008, 0.7921926)

phen.corrA10.3_3 <- cor.test(ambient.trivar$seed_pods, ambient.trivar$leaf, test="pearson")
phen.corrA10.3_3 # Corr = 0.5927418 | Posterior 95% CI = (0.5582001, 0.6252323)


#Model 10.3b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Ambient_Model_10.3b.RData")

ambient.trivar <- MCMC.ambient %>% filter(germ==1, !height==0)
prior10.3b <- list(R=list(V=diag(3), nu=2.002), 
                   G = list(G1 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000),
                            G2 = list(V=diag(3), nu = 1, alpha.mu = c(0,0,0), alpha.V = diag(3)*1000)))

AM_model10.3b <- MCMCglmm(cbind(height, leaf, seed_pods) ~ trait + plot - 1,
           random = ~us(trait):animal + us(trait):matID,
           ginverse = list(animal = Ainv), rcov=~us(trait):units,
           family = c("gaussian","gaussian","poisson"), 
           data = ambient.trivar, prior = prior10.3b,
           nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

#Some fixed effects not estimable and removed. Use an informative prior!

plot(AM_model10.3b$VCV) # trace plots not great for matID ... constricted by priors
summary(AM_model10.3b) # good effective sample sizes
autocorr.diag(AM_model10.3b$VCV) # no autocorrelation, but on the cusp
autocorr.plot(AM_model10.3b$VCV)
heidel.diag(AM_model10.3b$VCV) # convergence issues with matID

#Use model 10.3






#### No. 11 Stem Diameter - Univariate | Gaussian ####
#Note: For stem diameter, we only include plants that germinated
#We use an uninformative prior that divides the total phenotypic variance by the number of random effects included
#Variance in Total Fitness = squared standard deviation
#We treat plot as a random effect as from previous analysis (01_Data_exploration), we observe that there is a plot effect on this trait... There is likely some sampling error caused by plot due to the placement of the physical plots at KSR with respect to the tree line, etc so we want to give less statistical weight to the effect of plot. We also want to estimate the effect of plot 

#Calculating Variance for stem diameter for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1, !stem_diam==0) %>%
  summarise(mean=mean(stem_diam), sd=sd(stem_diam), var=(sd)^2) #var = 1.623994 *******

#Producing filtered data set for stem diameter
ambient.stem <- MCMC.ambient %>% filter(germ==1)
#From Committee Report #2: stem diameter follows a Polynomial distribution

#Full model including Random: Va, Vm, and | Fixed: plot with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_11.1.RData")

ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior11.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model11.1 <- MCMCglmm(stem_diam ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal = Ainv),
                     family = "gaussian", data = ambient.stem, prior = prior11.1, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(AM_model11.1, file="Routput/MCMC_Ambient_Model_11.1.RData")

plot(mode11.1$VCV) # 
summary(AM_model11.1) #
autocorr.diag(AM_model11.1$VCV) # 
heidel.diag(AM_model11.1$VCV) # 


#Model 11.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_11.2.RData")

ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior11.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model11.2 <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "gaussian", data = ambient.stem, prior = prior11.2, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable have been removed. Use an informative prior!

plot(AM_model11.2$VCV) # poor trace plots 
summary(AM_model11.2) # poor effective sample size
autocorr.diag(AM_model11.2$VCV) # high autocorrelation
heidel.diag(AM_model11.2$VCV) # convergence fail

#Model 11.2b using a Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Ambient_Model_11.2b.RData")

ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior11.2b <- list(R = list(V = 1, nu = 0.002),
                   G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                            G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                            G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model11.2b <- MCMCglmm(stem_diam ~ plot, random = ~animal + animalDom + matID,
                      ginverse = list(animal = Ainv, animalDom = Dinv),
                      family = "gaussian", data = ambient.stem, prior = prior11.2b, #gaussian distribution
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable have been removed. Use an informative prior!

plot(AM_model11.2b$VCV) # good trace plots
summary(AM_model11.2b) # good effective sample size
autocorr.diag(AM_model11.2b$VCV) # no autocorrelation
heidel.diag(AM_model11.2b$VCV) # convergence success

AM_herit11.2b <- AM_model11.2b$VCV[, "animal"]/(AM_model11.2b$VCV[, "animal"] + (AM_model11.2b$VCV[, "animalDom"]) + AM_model11.2b$VCV[, "matID"] + AM_model11.2b$VCV[, "units"])
mean(AM_herit11.2b)
HPDinterval(AM_herit11.2b)

#Model 11.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_11.3.RData")

ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior11.3 <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 0.002),
                           G2 = list(V = 1, nu = 0.002)))

AM_model11.3 <- MCMCglmm(stem_diam ~ plot, random = ~animal + matID,
                      ginverse = list(animal = Ainv),
                      family = "gaussian", data = ambient.stem, prior = prior11.3, #gaussian distribution
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable have been removed. Use an informative prior!

plot(AM_model11.3$VCV) # not great trace plot of animal
summary(AM_model11.3) # low effective sample size
autocorr.diag(AM_model11.3$VCV) # high autocorrelation
heidel.diag(AM_model11.3$VCV) # convergence success

#Model 11.3b using Fisher Parameter Expansion
#Note that parameter expansion is not ready for residual variance
load(file="Routput/MCMC_Ambient_Model_11.3b.RData")

ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior11.3b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model11.3b <- MCMCglmm(stem_diam ~ plot, random = ~animal + matID,
                      ginverse = list(animal = Ainv),
                      family = "gaussian", data = ambient.stem, prior = prior11.3b, #gaussian distribution
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable have been removed. Use an informative prior!

plot(AM_model11.3b$VCV) # better looking trace plot 
summary(AM_model11.3b) # good effective sample size
autocorr.diag(AM_model11.3b$VCV) # no autocorrelation
heidel.diag(AM_model11.3b$VCV) # convergence success

AM_herit11.3b <- AM_model11.3b$VCV[, "animal"]/(AM_model11.3b$VCV[, "animal"] + AM_model11.3b$VCV[, "matID"] + AM_model11.3b$VCV[, "units"])
mean(AM_herit11.3b)
HPDinterval(AM_herit11.3b)










#### No. 12 Five Trait - Multivariate | Ord, Gaus,Gaus, Ord, Pois ####
#Note: For this G matrix, we only use plants that germinated
#Traits: germination date (census), leaf #, height (cm), flowering date (census), seed pod #, relative fitness

#Producing filtered data set for height
ambient.multi <- MCMC.ambient %>% filter(germ==1)

#Model 12.2 with 3 additive, dominance and maternal effects + fixed plot effects
#We use a alpha.V value of 1 because we're using a prior that is suitable for an ordinal distribution
#Although family='ordinal' is available, we use family='threshold' because it is better at mixing and is suggested by Hadfield (2015)
ambient.multi <- MCMC.ambient %>% filter(germ==1)
prior12.2 <- list(R=list(V=diag(5), fix = 1), 
                  G = list(G1 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G2 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G3 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5))))

AM_model12.2 <- MCMCglmm(cbind(germ_census, leaf, height, flwr_census, seed_pods) ~ trait + plot - 1,
                         random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                         ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~us(trait):units,
                         family = c("threshold", "gaussian", "gaussian", "threshold", "poisson"), 
                         data = ambient.multi, prior = prior12.2,
                         nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(AM_model12.2, file="Routput/MCMC_Ambient_Model_12.2.RData")


plot(AM_model12.2$VCV) #
summary(AM_model12.2) #
autocorr.diag(AM_model12.2$VCV) #

#Model 12.3 with additive and maternal effects only
ambient.multi <- MCMC.ambient %>% filter(germ==1)

prior12.3 <- list(R=list(V=diag(5), fix = 1), 
                  G = list(G1 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G2 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5)),
                           G3 = list(V=diag(5), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(5))))

AM_model12.3 <- MCMCglmm(cbind(germ_census, leaf, height, flwr_census, seed_pods) ~ trait - 1,
                         random = ~us(trait):animal + us(trait):matID,
                         ginverse = list(animal = Ainv), rcov=~us(trait):units,
                         family = c("threshold", "gaussian","gaussian", "threshold", "poisson"), 
                         data = ambient.multi, prior = prior12.3,
                         nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(AM_model12.3, file="Routput/MCMC_Ambient_Model_12.3.RData")


#G-structure 3 is ill conditioned: use proper priors if you haven't or rescale data if you have

plot(AM_model12.3$VCV) #
summary(AM_model12.3) #
autocorr.diag(AM_model12.3$VCV) #



#Model 12.3b with 3 additive effects but no plot effect + scaled traits
ambient.multi <- MCMC.ambient %>% filter(germ==1)

#Only standarizing 3 traits - the census traits are not standardized
ambient.multi$leaf.st <- scale(ambient.multi$leaf)
ambient.multi$height.st <- scale(ambient.multi$height)
ambient.multi$seed_pods.st <- scale(ambient.multi$seed_pods)

ambient.multi %>% 
  summarise(germ.c.mean = mean(germ_census), flwr.c.mean = mean(germ_census),
            germ.c.sd = sd(germ_census), flwr.c.sd = sd(flwr_census))

hist(log(ambient.multi$seed_pods+1))
prior12.3b <- list(R=list(V=diag(6), fix = 1), 
                  G = list(G1 = list(V=diag(6), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(6)),
                           G2 = list(V=diag(6), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(6)),
                           G3 = list(V=diag(6), nu = 1, alpha.mu = c(0,0,0,0,0), alpha.V = diag(6))))

AM_model12.3b <- MCMCglmm(cbind(relative, germ_census, leaf.st, height.st, flwr_census, seed_pods.st) ~ trait - 1,
                         random = ~us(trait):animal + us(trait):animalDom + us(trait):matID,
                         ginverse = list(animalDom = Dinv), rcov=~us(trait):units,
                         family = c("threshold", "gaussian","gaussian", "threshold", "poisson"), 
                         data = ambient.multi, prior = prior12.3b,
                         nitt=2100000, thin=1000, burnin=100000, verbose = T, pr = TRUE) 

save(AM_model12.3b, file="Routput/MCMC_Ambient_Model_12.3b.RData")

plot(AM_model12.3b$VCV) #
summary(AM_model12.3b) #
autocorr.diag(AM_model12.3b$VCV) #












#### No. 13 Number of Flowering Clusters - Univariate | Poisson ####
#Note: For Flowering Clusters, we only include plants that survived to reach flowering


#Calculating Variance for Flowering Clusters for the Non-informative Prior (equal variance)
MCMC.ambient %>%
  filter(germ==1 & flower==1) %>%
  summarise(mean=mean(flwr_clstr), sd=sd(flwr_clstr), var=(sd)^2) #var = 2.535085 *******

#Producing filtered data set for Flowering Clusters
ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)

#Full model including Random: Va, Vm, and | Fixed: plot with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_13.1.RData")

ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior13.1 <- list(R = list(V = 1, nu = 0.002), 
                  G = list(G1 = list(V = 1, nu = 0.002),
                           G2 = list(V = 1, nu = 0.002), 
                           G3 = list(V = 1, nu = 0.002)))

AM_model13.1 <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal=Ainv),
                     family = "poisson", data = ambient.flwrclstr, prior = prior13.1, #poission distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(AM_model13.1, file="Routput/MCMC_Ambient_Model_13.1.RData")

plot(AM_model13.1$VCV) # 
summary(AM_model13.1) # 
autocorr.diag(AM_model13.1$VCV) # 
heidel.diag(AM_model13.1$VCV) # 

#Model 13.1b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Ambient_Model_13.1b.RData")

ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior13.1b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                          G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model13.1b <- MCMCglmm(seed_pods ~ plot, random = ~animal + matID + , 
                     ginverse = list(animal=Ainv),
                     family = "poisson", data = ambient.flwrclstr, prior = prior13.1b, #poission distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(AM_model13.1b, file="Routput/MCMC_Ambient_Model_13.1b.RData")

plot(AM_model13.1b$VCV) # 
summary(AM_model13.1b) # 
autocorr.diag(AM_model13.1b$VCV) # 
heidel.diag(AM_model13.1b$VCV) # 




#Model 13.2 Additive, Dominance, and Maternal Effects with inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_13.2.RData")

ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior13.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

AM_model13.2 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "poisson", data = ambient.flwrclstr, prior = prior13.2,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(AM_model13.2$VCV) # poor trace plot
summary(AM_model13.2) # poor effective sample size
autocorr.diag(AM_model13.2$VCV) # high autocorrelation
heidel.diag(AM_model13.2$VCV) # convergence fail

#Model 13.2b with Fisher Parameter Expanded Prior
load(file="Routput/MCMC_Ambient_Model_13.2b.RData")

ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior13.2b <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                           G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), 
                           G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model13.2b <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + animalDom + matID,
                      ginverse = list(animalDom= Dinv),
                      family = "poisson", data = ambient.flwrclstr, prior = prior13.2b,
                      nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effets, but use an informative prior!

hist(log(ambient.flwrclstr$flwr_clstr), breaks=30)
plot(AM_model13.2b$VCV) # OK trace plot
summary(AM_model13.2b) # moderate effective sample size
autocorr.diag(AM_model13.2b$VCV) # high autocorrelation
heidel.diag(AM_model13.2b$VCV) # convergence success

AM_herit13.2b <- AM_model13.2b$VCV[, "animal"]/(AM_model13.2b$VCV[, "animal"] + AM_model13.2b$VCV[, "animalDom"] + AM_model13.2b$VCV[, "matID"])
mean(AM_herit13.2b)
HPDinterval(AM_herit13.2b)

#Transformating data using "QGglmm" to get estimates of Va on data scale
yhat_A13.2b <- predict(AM_model13.2b, type = "terms")
mu_A13.2b <- mean(AM_model13.2b[["Sol"]][ , "(Intercept)"])
va_A13.2b <- mean(AM_model13.2b[["VCV"]][ , "animal"])
vp_A13.2b <- mean(rowSums(AM_model13.2b[["VCV"]]))

QG_A13.2b <- QGparams(predict=yhat_A13.2b, mu=mu_A13.2b, var.a=va_A13.2b, var.p=vp_A13.2b, model = "Poisson.log")
QG_A13.2b
va_A13.2b / vp_A13.2b






#Model 13.3 with additive + maternal effects + plot (fixed) + inverse gamma prior (0.001, 0.001)
load(file="Routput/MCMC_Ambient_Model_13.3.RData")

ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior13.3 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)))

AM_model13.3 <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "poisson", data = ambient.flwrclstr, prior = prior13.3,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#Some fixed effects not estimatable and removed. Use an informative prior!

plot(AM_model13.3$VCV) # good trace plots
summary(AM_model13.3) # poor effective sample size
autocorr.diag(AM_model13.3$VCV) # high autocorrelation for animal
heidel.diag(AM_model13.3$VCV) # convergence success

#Model 13.3b attempted with Fisher Parameter Expanded Priors
load(file="Routput/MCMC_Ambient_Model_13.3b.RData")

ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior13.3b <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                          G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_model13.3b <- MCMCglmm(flwr_clstr ~ plot, random = ~animal + matID,
                     ginverse = list(animal = Ainv),
                     family = "poisson", data = ambient.flwrclstr, prior = prior13.3b,
                     nitt = 3100000, thin = 1500, burnin = 100000, verbose = T, pr = TRUE)

#Increased length to 3 100 000 iterations + increased thinning interval to prevent autocorrelation
#If this doesn't solve the problem, then consider reducing alpha.V to <1000
#In the mean while going to report 
#Some fixed effects not estimatable and removed. Use an informative prior!

plot(AM_model13.3b$VCV) # good trace plot
summary(AM_model13.3b) # good effective sample size (but not exactly 2000 perfect)
autocorr.diag(AM_model13.3b$VCV) # high autocorrelation .. if increase nitt and thin, still autocorrelation but reduced
heidel.diag(AM_model13.3b$VCV) # convergence success

AM_herit13.3b <- AM_model13.3b$VCV[, "animal"]/(AM_model13.3b$VCV[, "animal"] + AM_model13.3b$VCV[, "matID"] + AM_model13.3b$VCV[, "units"])
mean(AM_herit13.3b)
HPDinterval(AM_herit13.3b)

yhat_A13.3b <- predict(AM_model13.3b, type = "terms")
mu_A13.3b <- mean(AM_model13.3b[["Sol"]][ , "(Intercept)"])
va_A13.3b <- mean(AM_model13.3b[["VCV"]][ , "animal"])
vp_A13.3b <- mean(rowSums(AM_model13.3b[["VCV"]]))

QG_A13.3b <- QGparams(predict=yhat_A13.3b, mu=mu_A13.3b, var.a=va_A13.3b, var.p=vp_A13.3b, model = "Poisson.log")
QG_A13.3b
va_A13.3b / vp_A13.3b


#Model Comparison
AM_model13.2b$DIC
AM_model13.3b$DIC