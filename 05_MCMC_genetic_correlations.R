#### PROJECT: Brassica rapa GxE Study (Data collected at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Cross-environment models to estimate Va, Vd, Vm and genetic correlation resampling

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Necessary Packages
library(tidyverse)
library(MCMCglmm)
library(QGglmm)
library(parallel)
library(MuMIn)
library(logitnorm)
library(nadiv)
library(lattice)
library(gtools)

#Importing Data using readr from tidyverse (Base R confuses the class for certain vectors). 
# d for double, f for factor, i for integer
col_types_list2 <- cols_only(posID = "d", individual = "d", plot = "f",
                             animal = "f", patID = "f", matID = "f", famID = "f",
                             treatment = col_factor(levels=c("A", "H")),
                             germ = "d", flower = "d", 
                             seed = "d", leaf = "d", flwr_clstr = "d", seed_pods = "d", seed_number = "d", 
                             height = "d", stem_diam = "d", germ_census = "d", flwr_census = "d")

MCMC.2019 <- read_csv("Rdata/MCMC_2019_cleaned.csv", col_names = TRUE, na = "NA",
                      col_types = col_types_list2)
spec(MCMC.2019)

ped <- read.csv("Rdata/heatarrays_animal_pedigree.csv") #note don't use readr to import the csv file.
for (x in 1:3) ped[, x] <- as.factor(ped[, x])

#Changing tibbles into dataframes because MCMCglmm cannot read tibbles
MCMC.2019 <- as.data.frame(MCMC.2019)
ped <- as.data.frame(ped)
lapply(ped, class)
lapply(MCMC.2019, class)

#Creating Dominance Matrix for Vd estimates
Ainv <- inverseA(ped[, 1:3])$Ainv
Dinv <- makeD(ped[, 1:3])$Dinv
MCMC.2019$animalDom <- MCMC.2019$animal

MCMC.ambient <- MCMC.2019 %>%
  filter(treatment=="A")
MCMC.heated <- MCMC.2019 %>%
  filter(treatment=="H")

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




#'###############################################################################'#
###############      CROSS ENVIRONMENT GENETIC CORRELATIONS      ################
#'###############################################################################'#

#Cross environment models (above) did not converge well and many models have autocorrelation.
#We instead estimate cross-environment additive genetic correlations using the separate environment models (scripts 03 and 04)
#We then check whether the cross-environment genetic correlation (of F2 parents since the only 1 individual in the F3 is expressed in one environment) is significantly different from 0. 

# ___ OTHER CHECKS TO DO ___

# (1) If coefficient of correlation b/w Ambient and Heated models = coefficient of correlation in full model
# (2) If coefficient of correlation b/w Ambient and Heated models = coefficient of correlation b/w models w/o dominance

## Trait 1: Survival Success (Binary) ####

  ## STEP 1 - Testing if r [individual models] = r [full models]]

#Genetic Correlation for models with Additive + Maternal Effects
load(file="Routput/Ambient/AM_model3.1.RData")
load(file="Routput/Heated/HW_model3.1.RData")
AM_3.1_BV <- AM_model3.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_3.1_BV <- HW_model3.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
r.3.1 <- 1:2000 
for (i in 1:2000) {
  r.3.1[i] <- (cov(AM_3.1_BV[i,], HW_3.1_BV[i,])) /
    sqrt(var(AM_3.1_BV[i,]) * var(HW_3.1_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.3.1), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.3.1)) # 95% interval
P.r.3.1 <- ifelse(r.3.1 <0, 0, 1) 
(1-sum(P.r.3.1)/2000)*2 # 


#Genetic correlation for full models
load(file="Routput/Full_Excl_Dom/Plasticity_1.5.RData")

#Extracting breeding values for each environment
f1.5_HW_BV <- (as.mcmc(PL_model1.5$Sol[,7799:7934])) # index [, 7799:7934] is HW BV 
f1.5_AM_BV <- (as.mcmc(PL_model1.5$Sol[,148:283])) #index [,148:283] is AM BV

#Resampling posterior to estimate correlation between BV from each environment
r.f1.5 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.f1.5[i] <- (cov(f1.5_AM_BV[i,], f1.5_HW_BV[i,])) /
    sqrt(var(f1.5_AM_BV[i,]) * var(f1.5_HW_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.f1.5), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.f1.5)) # 95% interval
P.r.f1.5 <- ifelse(r.f1.5 <0, 0, 1) 
(1-sum(P.r.f1.5)/2000)*2 # 

#Testing if different approaches to estimating the genetic correlation differs
P.diff1 <- ifelse(r.3.1 < r.f1.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff1)/2000)*2  # P = 0.905
    

    ## STEP 2 - Extracting r [individual models] including additive, dominance, maternal effects

#Genetic Correlation for models with Additive + Dominance + Maternal Effects
load(file="Routput/Ambient/AM_model3.0.RData")
load(file="Routput/Heated/HW_model3.0.RData")
AM_3.0_BV <- AM_model3.0$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_3.0_BV <- HW_model3.0$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.3.0 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.3.0[i] <- (cov(AM_3.0_BV[i,], HW_3.0_BV[i,])) /
    sqrt(var(AM_3.0_BV[i,]) * var(HW_3.0_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.3.0), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.3.0)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.3.0 <- ifelse(r.3.0 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.3.0)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 

# genetic variance for plasticity!
Va_plasticity3.0 <- numeric(2000)
for (i in 1:2000) {
  
  # this is the additive genetic variance for plasticity
  # it's the covariance of A and B divided by the square root of the product of their variances
  plasticity3.0 <- (AM_3.0_BV[i,] - HW_3.0_BV[i,]) 
  Va_plasticity3.0[i] <- var(plasticity3.0)
}
posterior.mode(as.mcmc(Va_plasticity3.0), na.rm=T) 
HPDinterval(as.mcmc(Va_plasticity3.0))
P_Va3.0 <- ifelse(Va_plasticity3.0 < 0, 0, 1)
(1-sum(P_Va3.0)/2000)*2


#Testing if including dominance changes anything to the genetic correlation (P value)
P.diff1b <- ifelse(r.3.0 > r.3.1, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff1)/2000)*2  # P = 0.792

#Plotting breeding values 
AM_3.0_BV <- colMeans(as.data.frame(AM_3.0_BV))
HW_3.0_BV <- colMeans(as.data.frame(HW_3.0_BV))
BV_surv <- bind_cols(AM_3.0_BV, HW_3.0_BV)
colnames(BV_surv)[1] <- "AM_surv_BV"
colnames(BV_surv)[2] <- "HW_surv_BV"
BV_surv %>% 
  ggplot(aes(x=AM_surv_BV, y=HW_surv_BV)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") + #1 to 1 line is dashed
  geom_smooth()

#










## Trait 2: Fecundity (Poisson) ####

    ## STEP 1 - Testing if r [individual models] = r [full models]]

#Genetic Correlation for models with Additive + Maternal Effects
load(file="Routput/Ambient/AM_model4.3.RData")
load(file="Routput/Heated/HW_model4.3.RData")
AM_4.3_BV <- AM_model4.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_4.3_BV <- HW_model4.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
r.4.3 <- 1:2000 
for (i in 1:2000) {
  r.4.3[i] <- (cov(AM_4.3_BV[i,], HW_4.3_BV[i,])) /
    sqrt(var(AM_4.3_BV[i,]) * var(HW_4.3_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.4.3), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.4.3)) # 95% interval
P.r.4.3 <- ifelse(r.4.3 <0, 0, 1) 
(1-sum(P.r.4.3)/2000)*2 # 


#Genetic correlation for full models
load(file="Routput/Full_Excl_Dom/Plasticity_2.5.RData")

#Extracting breeding values for each environment
f2.5_HW_BV <- (as.mcmc(PL_model2.5$Sol[,7799:7934])) # index [, 7799:7934] is HW BV 
f2.5_AM_BV <- (as.mcmc(PL_model2.5$Sol[,148:283])) #index [,148:283] is AM BV

#Resampling posterior to estimate correlation between BV from each environment
r.f2.5 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.f2.5[i] <- (cov(f2.5_AM_BV[i,], f2.5_HW_BV[i,])) /
    sqrt(var(f2.5_AM_BV[i,]) * var(f2.5_HW_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.f2.5), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.f2.5)) # 95% interval
P.r.f2.5 <- ifelse(r.f2.5 <0, 0, 1) 
(1-sum(P.r.f2.5)/2000)*2 # 

#Testing if different approaches to estimating the genetic correlation differs
P.diff2 <- ifelse(r.4.3 < r.f2.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff2)/2000)*2  # P = 0.915


    ## STEP 2 - Extracting r [individual models] including additive, dominance, maternal effects

#Genetic Correlation for models with Additive + Dominance + Maternal Effects
load(file="Routput/Ambient/AM_model4.1.RData")
load(file="Routput/Heated/HW_model4.1.RData")
AM_4.1_BV <- AM_model4.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_4.1_BV <- HW_model4.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file


#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.4.1 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.4.1[i] <- (cov(AM_4.1_BV[i,], HW_4.1_BV[i,])) /
    sqrt(var(AM_4.1_BV[i,]) * var(HW_4.1_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.4.1), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.4.1)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.4.1 <- ifelse(r.4.1 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.4.1)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 

# genetic variance for plasticity!
Va_plasticity4.1 <- numeric(2000)
for (i in 1:2000) {
  
  # this is the additive genetic variance for plasticity
  # it's the covariance of A and B divided by the square root of the product of their variances
  plasticity4.1 <- (AM_4.1_BV[i,] - HW_4.1_BV[i,]) 
  Va_plasticity4.1[i] <- var(plasticity4.1)
}
posterior.mode(as.mcmc(Va_plasticity4.1), na.rm=T) 
HPDinterval(as.mcmc(Va_plasticity4.1))
P_Va4.1 <- ifelse(Va_plasticity4.1 < 0, 0, 1)
(1-sum(P_Va4.1)/2000)*2

#Testing if including dominance changes anything to the genetic correlation (P value)
P.diff2b <- ifelse(r.4.1 > r.4.3, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff2b)/2000)*2  # P = **

#Plotting breeding values
AM_4.1_BV <- colMeans(as.data.frame(AM_4.1_BV))
HW_4.1_BV <- colMeans(as.data.frame(HW_4.1_BV))
BV_fecund <- bind_cols(AM_4.1_BV, HW_4.1_BV)
colnames(BV_fecund)[1] <- "AM_fecund_BV"
colnames(BV_fecund)[2] <- "HW_fecund_BV"

BV_fecund %>% 
  ggplot(aes(x=AM_fecund_BV, y=HW_fecund_BV)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") + #1 to 1 line is dashed
  geom_smooth()

#Converting BVs from latent to data scale using QGpsi (see Eqs 14 and 15 of de Villemereuil et al. 2016)
vp_AM4.1 <- mean(rowSums(AM_model4.1[["VCV"]])) #phenotypic variance
vp_HW4.1 <- mean(rowSums(HW_model4.1[["VCV"]])) #phenotypic variance
dinv <- function(x) {exp(x)} #derivative of inverse poisson log link
#Note derivative of inverse log link (exponential) is still exponential

#psi for linear function conversion from latent ot data scale
psi_AM4.1 <- QGpsi(mu = 0, var = vp_AM4.1, d.link.inv = dinv) 
psi_HW4.1 <- QGpsi(mu = 0, var = vp_HW4.1, d.link.inv = dinv)

#Applying psi
AM_4.1_BV <- AM_model4.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_4.1_BV <- HW_model4.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
AM_fecund_BV <- colMeans((psi_AM4.1 * AM_4.1_BV) - mean(psi_AM4.1 * AM_4.1_BV))
HW_fecund_BV <- colMeans((psi_HW4.1 * HW_4.1_BV) - mean(psi_HW4.1 * HW_4.1_BV))

BV_fecund <- bind_cols(AM_fecund_BV, HW_fecund_BV)
colnames(BV_fecund)[1] <- "AM_fecund_BV"
colnames(BV_fecund)[2] <- "HW_fecund_BV"

#Plotting mean breeding values again phenotypic mean
phen_fecund <- MCMC.2019 %>% group_by(treatment, famID, matID, patID) %>% summarise(mean_fecund=mean(seed_pods)) #calculating mean trait per treatment and parental ID
phen_fecund2 <- pivot_wider(phen_fecund, names_from = "treatment", values_from = "mean_fecund") #transposing so mean trait has col by treatment
BV_fecund$matID <- 136:271 #matID is the animal number of the parents
BV_fecund$matID <- as.factor(BV_fecund$matID)
BV_phen_fecund <- inner_join(phen_fecund2, BV_fecund, by = "matID") #combining BV (F2 parents) and phenotypic data (for F3 experimental), simultaenously removes rows of F2 parents without progeny in experiment
colnames(BV_phen_fecund)[4] <- "AM_fecund_mu"
colnames(BV_phen_fecund)[5] <- "HW_fecund_mu"


BV_phen_fecund %>% #ambient
  ggplot(aes(x=AM_fecund_BV, y=AM_fecund_mu)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept = 0, size = 1.2, linetype = "dashed") +
  geom_smooth(color="blue")

BV_phen_fecund %>% #heated
  ggplot(aes(x=HW_fecund_BV, y=HW_fecund_mu)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept = 0, size = 1.2, linetype = "dashed") +
  geom_smooth(color="red")

#








## Trait 3: Overwintering Survival / Germination Success (Binary) ####

    ## STEP 1 - Testing if r [individual models] = r [full models]]

#Genetic Correlation for models with Additive + Maternal Effects
load(file="Routput/Ambient/AM_model5.1.RData")
load(file="Routput/Heated/HW_model5.1.RData")
AM_5.1_BV <- AM_model5.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_5.1_BV <- HW_model5.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
r.5.1 <- 1:2000 
for (i in 1:2000) {
  r.5.1[i] <- (cov(AM_5.1_BV[i,], HW_5.1_BV[i,])) /
    sqrt(var(AM_5.1_BV[i,]) * var(HW_5.1_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.5.1), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.5.1)) # 95% interval
P.r.5.1 <- ifelse(r.5.1 <0, 0, 1) 
(1-sum(P.r.5.1)/2000)*2 # 


#Genetic correlation for full models
load(file="Routput/Full_Excl_Dom/Plasticity_3.5.RData")

#Extracting breeding values for each environment
f3.5_HW_BV <- (as.mcmc(PL_model3.5$Sol[,7799:7934])) # index [, 7799:7934] is HW BV 
f3.5_AM_BV <- (as.mcmc(PL_model3.5$Sol[,148:283])) #index [,148:283] is AM BV

#Resampling posterior to estimate correlation between BV from each environment
r.f3.5 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.f3.5[i] <- (cov(f3.5_AM_BV[i,], f3.5_HW_BV[i,])) /
    sqrt(var(f3.5_AM_BV[i,]) * var(f3.5_HW_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.f3.5), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.f3.5)) # 95% interval
P.r.f3.5 <- ifelse(r.f3.5 <0, 0, 1) 
(1-sum(P.r.f3.5)/2000)*2 # 

#Testing if different approaches to estimating the genetic correlation differs
P.diff3 <- ifelse(r.5.1 > r.f3.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff3)/2000)*2  # P = 0.995


    ## STEP 2 - Extracting r [individual models] including additive, dominance, maternal effects

#Genetic Correlation for models with Additive + Dominance + Maternal Effects
load(file="Routput/Ambient/AM_model5.0.RData")
load(file="Routput/Heated/HW_model5.0.RData")
AM_5.0_BV <- AM_model5.0$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_5.0_BV <- HW_model5.0$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.5.0 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.5.0[i] <- (cov(AM_5.0_BV[i,], HW_5.0_BV[i,])) /
    sqrt(var(AM_5.0_BV[i,]) * var(HW_5.0_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.5.0), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.5.0)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.5.0 <- ifelse(r.5.0 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.5.0)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 

# genetic variance for plasticity!
Va_plasticity5.0 <- numeric(2000)
for (i in 1:2000) {
  
  # this is the additive genetic variance for plasticity
  # it's the covariance of A and B divided by the square root of the product of their variances
  plasticity5.0 <- (AM_5.0_BV[i,] - HW_5.0_BV[i,]) 
  Va_plasticity5.0[i] <- var(plasticity5.0)
}
posterior.mode(as.mcmc(Va_plasticity5.0), na.rm=T) 
HPDinterval(as.mcmc(Va_plasticity5.0))
P_Va5.0 <- ifelse(Va_plasticity5.0 < 0, 0, 1)
(1-sum(P_Va5.0)/2000)*2

#Testing if including dominance changes anything to the genetic correlation (P value)
P.diff3b <- ifelse(r.5.0 > r.5.1, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff3b)/2000)*2  # P = 0.865

#Plotting breeding values 
AM_5.0_BV <- colMeans(as.data.frame(AM_5.0_BV))
HW_5.0_BV <- colMeans(as.data.frame(HW_5.0_BV))
BV_germ <- bind_cols(AM_5.0_BV, HW_5.0_BV)
colnames(BV_germ)[1] <- "AM_germ_BV"
colnames(BV_germ)[2] <- "HW_germ_BV"
BV_germ %>% 
  ggplot(aes(x=AM_germ_BV, y=HW_germ_BV)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") + #1 to 1 line is dashed
  geom_smooth()

#






## Trait 4: Spring-summer Survival / Flowering Success (Binary) ####

    ## STEP 1 - Testing if r [individual models] = r [full models]]

#Genetic Correlation for models with Additive + Maternal Effects
load(file="Routput/Ambient/AM_model6.1.RData")
load(file="Routput/Heated/HW_model6.1.RData")
AM_6.1_BV <- AM_model6.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_6.1_BV <- HW_model6.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
r.6.1 <- 1:2000 
for (i in 1:2000) {
  r.6.1[i] <- (cov(AM_6.1_BV[i,], HW_6.1_BV[i,])) /
    sqrt(var(AM_6.1_BV[i,]) * var(HW_6.1_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.6.1), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.6.1)) # 95% interval
P.r.6.1 <- ifelse(r.6.1 <0, 0, 1) 
(1-sum(P.r.6.1)/2000)*2 # 


#Genetic correlation for full models
load(file="Routput/Full_Excl_Dom/Plasticity_4.5.RData")

#Extracting breeding values for each environment
f4.5_HW_BV <- (as.mcmc(PL_model4.5$Sol[,7799:7934])) # index [, 7799:7934] is HW BV 
f4.5_AM_BV <- (as.mcmc(PL_model4.5$Sol[,148:283])) #index [,148:283] is AM BV

#Resampling posterior to estimate correlation between BV from each environment
r.f4.5 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.f4.5[i] <- (cov(f4.5_AM_BV[i,], f4.5_HW_BV[i,])) /
    sqrt(var(f4.5_AM_BV[i,]) * var(f4.5_HW_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.f4.5), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.f4.5)) # 95% interval
P.r.f4.5 <- ifelse(r.f4.5 <0, 0, 1) 
(1-sum(P.r.f4.5)/2000)*2 # 

#Testing if different approaches to estimating the genetic correlation differs
P.diff4 <- ifelse(r.6.1 < r.f4.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff4)/2000)*2  # P = 0.957


    ## STEP 2 - Extracting r [individual models] including additive, dominance, maternal effects

#Genetic Correlation for models with Additive + Dominance + Maternal Effects
load(file="Routput/Ambient/AM_model6.0.RData")
load(file="Routput/Heated/HW_model6.0.RData")
AM_6.0_BV <- AM_model6.0$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_6.0_BV <- HW_model6.0$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.6.0 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.6.0[i] <- (cov(AM_6.0_BV[i,], HW_6.0_BV[i,])) /
    sqrt(var(AM_6.0_BV[i,]) * var(HW_6.0_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.6.0), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.6.0)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.6.0 <- ifelse(r.6.0 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.6.0)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 

# genetic variance for plasticity!
Va_plasticity6.0 <- numeric(2000)
for (i in 1:2000) {
  
  # this is the additive genetic variance for plasticity
  # it's the covariance of A and B divided by the square root of the product of their variances
  plasticity6.0 <- (AM_6.0_BV[i,] - HW_6.0_BV[i,]) 
  Va_plasticity6.0[i] <- var(plasticity6.0)
}
posterior.mode(as.mcmc(Va_plasticity6.0), na.rm=T) 
HPDinterval(as.mcmc(Va_plasticity6.0))
P_Va6.0 <- ifelse(Va_plasticity6.0 < 0, 0, 1)
(1-sum(P_Va6.0)/2000)*2

#Testing if including dominance changes anything to the genetic correlation (P value)
P.diff4b <- ifelse(r.6.0 > r.6.1, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff4b)/2000)*2  # P = 0.792

#Plotting breeding values 
AM_6.0_BV <- colMeans(as.data.frame(AM_6.0_BV))
HW_6.0_BV <- colMeans(as.data.frame(HW_6.0_BV))
BV_flwr <- bind_cols(HW_6.0_BV, AM_6.0_BV)
colnames(BV_flwr)[1] <- "HW_flwr_BV"
colnames(BV_flwr)[2] <- "AM_flwr_BV"
BV_flwr %>% 
  ggplot(aes(x=AM_flwr_BV, y=HW_flwr_BV)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") + #1 to 1 line is dashed
  geom_smooth()

#






## Trait 5: Seed Maturation Success (Binary) ####

    ## STEP 1 - Testing if r [individual models] = r [full models]]

#Genetic Correlation for models with Additive + Maternal Effects
load(file="Routput/Ambient/AM_model7.1.RData")
load(file="Routput/Heated/HW_model7.1.RData")
AM_7.1_BV <- AM_model7.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_7.1_BV <- HW_model7.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
r.7.1 <- 1:2000 
for (i in 1:2000) {
  r.7.1[i] <- (cov(AM_7.1_BV[i,], HW_7.1_BV[i,])) /
    sqrt(var(AM_7.1_BV[i,]) * var(HW_7.1_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.7.1), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.7.1)) # 95% interval
P.r.7.1 <- ifelse(r.7.1 <0, 0, 1) 
(1-sum(P.r.7.1)/2000)*2 # 


#Genetic correlation for full models
load(file="Routput/Full_Excl_Dom/Plasticity_5.5.RData")

#Extracting breeding values for each environment
f5.5_HW_BV <- (as.mcmc(PL_model5.5$Sol[,7799:7934])) # index [, 7799:7934] is HW BV 
f5.5_AM_BV <- (as.mcmc(PL_model5.5$Sol[,148:283])) #index [,148:283] is AM BV

#Resampling posterior to estimate correlation between BV from each environment
r.f5.5 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.f5.5[i] <- (cov(f5.5_AM_BV[i,], f5.5_HW_BV[i,])) /
    sqrt(var(f5.5_AM_BV[i,]) * var(f5.5_HW_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.f5.5), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.f5.5)) # 95% interval
P.r.f5.5 <- ifelse(r.f5.5 >0, 0, 1) 
(1-sum(P.r.f5.5)/2000)*2 # 

#Testing if different approaches to estimating the genetic correlation differs
P.diff5 <- ifelse(r.5.1 < r.f5.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff5)/2000)*2  # P = 0.791


    ## STEP 2 - Extracting r [individual models] including additive, dominance, maternal effects

#Genetic Correlation for models with Additive + Dominance + Maternal Effects
load(file="Routput/Ambient/AM_model7.0.RData")
load(file="Routput/Heated/HW_model7.0.RData")
AM_7.0_BV <- AM_model7.0$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_7.0_BV <- HW_model7.0$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.7.0 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.7.0[i] <- (cov(AM_7.0_BV[i,], HW_7.0_BV[i,])) /
    sqrt(var(AM_7.0_BV[i,]) * var(HW_7.0_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.7.0), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.7.0)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.7.0 <- ifelse(r.7.0 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.7.0)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 

# genetic variance for plasticity!
Va_plasticity7.0 <- numeric(2000)
for (i in 1:2000) {
  
  # this is the additive genetic variance for plasticity
  # it's the covariance of A and B divided by the square root of the product of their variances
  plasticity7.0 <- (AM_7.0_BV[i,] - HW_7.0_BV[i,]) 
  Va_plasticity7.0[i] <- var(plasticity7.0)
}
posterior.mode(as.mcmc(Va_plasticity7.0), na.rm=T) 
HPDinterval(as.mcmc(Va_plasticity7.0))
P_Va7.0 <- ifelse(Va_plasticity7.0 < 0, 0, 1)
(1-sum(P_Va7.0)/2000)*2

#Testing if including dominance changes anything to the genetic correlation (P value)
P.diff5b <- ifelse(r.7.0 > r.f5.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff5b)/2000)*2  # P = 1.021

#Plotting breeding values 
AM_7.0_BV <- colMeans(as.data.frame(AM_7.0_BV))
HW_7.0_BV <- colMeans(as.data.frame(HW_7.0_BV))
BV_seed <- bind_cols(HW_7.0_BV, AM_7.0_BV)
colnames(BV_seed)[1] <- "HW_seed_BV"
colnames(BV_seed)[2] <- "AM_seed_BV"
BV_seed %>% 
  ggplot(aes(x=AM_seed_BV, y=HW_seed_BV)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") + #1 to 1 line is dashed
  geom_smooth()

#







## Trait 6: Leaf Number (Gaussian) ####

## STEP 1 - Testing if r [individual models] = r [full models]]

#Genetic Correlation for models with Additive + Maternal Effects
load(file="Routput/Ambient/AM_model8.3.RData")
load(file="Routput/Heated/HW_model8.3.RData")
AM_8.3_BV <- AM_model8.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_8.3_BV <- HW_model8.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
r.8.3 <- 1:2000 
for (i in 1:2000) {
  r.8.3[i] <- (cov(AM_8.3_BV[i,], HW_8.3_BV[i,])) /
    sqrt(var(AM_8.3_BV[i,]) * var(HW_8.3_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.8.3), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.8.3)) # 95% interval
P.r.8.3 <- ifelse(r.8.3 <0, 0, 1) 
(1-sum(P.r.8.3)/2000)*2 # 


#Genetic correlation for full models
load(file="Routput/Full_Excl_Dom/Plasticity_6.5.RData")

#Extracting breeding values for each environment
f6.5_HW_BV <- (as.mcmc(PL_model6.5$Sol[,7799:7934])) # index [, 7799:7934] is HW BV 
f6.5_AM_BV <- (as.mcmc(PL_model6.5$Sol[,148:283])) #index [,148:283] is AM BV

#Resampling posterior to estimate correlation between BV from each environment
r.f6.5 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.f6.5[i] <- (cov(f6.5_AM_BV[i,], f6.5_HW_BV[i,])) /
    sqrt(var(f6.5_AM_BV[i,]) * var(f6.5_HW_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.f6.5), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.f6.5)) # 95% interval
P.r.f6.5 <- ifelse(r.f6.5 <0, 0, 1) 
(1-sum(P.r.f6.5)/2000)*2 # 

#Testing if different approaches to estimating the genetic correlation differs
P.diff6 <- ifelse(r.8.3 < r.f6.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(1-sum(P.diff6)/2000)*2  # P =


    ## STEP 2 - Extracting r [individual models] including additive, dominance, maternal effects

#Genetic Correlation for models with Additive + Dominance + Maternal Effects
load(file="Routput/Ambient/AM_model8.1.RData")
load(file="Routput/Heated/HW_model8.1.RData")
AM_8.1_BV <- AM_model8.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_8.1_BV <- HW_model8.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.8.1 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.8.1[i] <- (cov(AM_8.1_BV[i,], HW_8.1_BV[i,])) /
    sqrt(var(AM_8.1_BV[i,]) * var(HW_8.1_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.8.1), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.8.1)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.8.1 <- ifelse(r.8.1 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.8.1)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 

# genetic variance for plasticity!
Va_plasticity8.1 <- numeric(2000)
for (i in 1:2000) {
  
  # this is the additive genetic variance for plasticity
  # it's the covariance of A and B divided by the square root of the product of their variances
  plasticity8.1 <- (AM_8.1_BV[i,] - HW_8.1_BV[i,]) 
  Va_plasticity8.1[i] <- var(plasticity8.1)
}
posterior.mode(as.mcmc(Va_plasticity8.1), na.rm=T) 
HPDinterval(as.mcmc(Va_plasticity8.1))
P_Va8.1 <- ifelse(Va_plasticity8.1 < 0, 0, 1)
(1-sum(P_Va8.1)/2000)*2

#Testing if including dominance changes anything to the genetic correlation (P value)
P.diff6b <- ifelse(r.8.1 > r.8.3, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff6b)/2000)*2  # P = **

#Plotting breeding values 
AM_8.1_BV <- colMeans(as.data.frame(AM_8.1_BV))
HW_8.1_BV <- colMeans(as.data.frame(HW_8.1_BV))
BV_leaf <- bind_cols(HW_8.1_BV, AM_8.1_BV)
colnames(BV_leaf)[1] <- "HW_leaf_BV"
colnames(BV_leaf)[2] <- "AM_leaf_BV"
BV_leaf %>% 
  ggplot(aes(x=AM_leaf_BV, y=HW_leaf_BV)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") + #1 to 1 line is dashed
  geom_smooth()















## Trait 7: Height (Gaussian) ####

## STEP 1 - Testing if r [individual models] = r [full models]]

#Genetic Correlation for models with Additive + Maternal Effects
load(file="Routput/Ambient/AM_model9.3.RData")
load(file="Routput/Heated/HW_model9.3.RData")
AM_9.3_BV <- AM_model9.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_9.3_BV <- HW_model9.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
r.9.3 <- 1:2000 
for (i in 1:2000) {
  r.9.3[i] <- (cov(AM_9.3_BV[i,], HW_9.3_BV[i,])) /
    sqrt(var(AM_9.3_BV[i,]) * var(HW_9.3_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.9.3), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.9.3)) # 95% interval
P.r.9.3 <- ifelse(r.9.3 <0, 0, 1) 
(1-sum(P.r.9.3)/2000)*2 # 


#Genetic correlation for full models
load(file="Routput/Full_Excl_Dom/Plasticity_7.5.RData")

#Extracting breeding values for each environment
f7.5_HW_BV <- (as.mcmc(PL_model7.5$Sol[,7799:7934])) # index [, 7799:7934] is HW BV 
f7.5_AM_BV <- (as.mcmc(PL_model7.5$Sol[,148:283])) #index [,148:283] is AM BV

#Resampling posterior to estimate correlation between BV from each environment
r.f7.5 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.f7.5[i] <- (cov(f7.5_AM_BV[i,], f7.5_HW_BV[i,])) /
    sqrt(var(f7.5_AM_BV[i,]) * var(f7.5_HW_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.f7.5), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.f7.5)) # 95% interval
P.r.f7.5 <- ifelse(r.f7.5 <0, 0, 1) 
(1-sum(P.r.f7.5)/2000)*2 # 

#Testing if different approaches to estimating the genetic correlation differs
P.diff7 <- ifelse(r.9.3 > r.f7.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff7)/2000)*2  # P = 0.969


    ## STEP 2 - Extracting r [individual models] including additive, dominance, maternal effects

#Genetic Correlation for models with Additive + Dominance + Maternal Effects
load(file="Routput/Ambient/AM_model9.1.RData")
load(file="Routput/Heated/HW_model9.1.RData")
AM_9.1_BV <- AM_model9.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_9.1_BV <- HW_model9.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.9.1 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.9.1[i] <- (cov(AM_9.1_BV[i,], HW_9.1_BV[i,])) /
    sqrt(var(AM_9.1_BV[i,]) * var(HW_9.1_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.9.1), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.9.1)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.9.1 <- ifelse(r.9.1 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.9.1)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 

# genetic variance for plasticity!
Va_plasticity9.1 <- numeric(2000)
for (i in 1:2000) {
  
  # this is the additive genetic variance for plasticity
  # it's the covariance of A and B divided by the square root of the product of their variances
  plasticity9.1 <- (AM_9.1_BV[i,] - HW_9.1_BV[i,]) 
  Va_plasticity9.1[i] <- var(plasticity9.1)
}
posterior.mode(as.mcmc(Va_plasticity9.1), na.rm=T) 
HPDinterval(as.mcmc(Va_plasticity9.1))
P_Va9.1 <- ifelse(Va_plasticity9.1 < 0, 0, 1)
(1-sum(P_Va9.1)/2000)*2

#Testing if including dominance changes anything to the genetic correlation (P value)
P.diff7b <- ifelse(r.9.1 > r.9.3, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff7b)/2000)*2  # P = 0.371

#Plotting breeding values 
AM_9.1_BV <- colMeans(as.data.frame(AM_9.1_BV))
HW_9.1_BV <- colMeans(as.data.frame(HW_9.1_BV))
BV_height <- bind_cols(HW_9.1_BV, AM_9.1_BV)
colnames(BV_height)[1] <- "HW_height_BV"
colnames(BV_height)[2] <- "AM_height_BV"
BV_height %>% 
  ggplot(aes(x=AM_height_BV, y=HW_height_BV)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") + #1 to 1 line is dashed
  geom_smooth()

#





## Trait 8: Stem Diameter (Gaussian) ####

    ## STEP 1 - Testing if r [individual models] = r [full models]]

#Genetic Correlation for models with Additive + Maternal Effects
load(file="Routput/Ambient/AM_model10.3.RData")
load(file="Routput/Heated/HW_model10.3.RData")
AM_10.3_BV <- AM_model10.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_10.3_BV <- HW_model10.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
r.10.3 <- 1:2000 
for (i in 1:2000) {
  r.10.3[i] <- (cov(AM_10.3_BV[i,], HW_10.3_BV[i,])) /
    sqrt(var(AM_10.3_BV[i,]) * var(HW_10.3_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.10.3), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.10.3)) # 95% interval
P.r.10.3 <- ifelse(r.10.3 <0, 0, 1) 
(1-sum(P.r.10.3)/2000)*2 # 


#Genetic correlation for full models
load(file="Routput/Full_Excl_Dom/Plasticity_9.5.RData")

#Extracting breeding values for each environment
f9.5_HW_BV <- (as.mcmc(PL_model9.5$Sol[,7799:7934])) # index [, 7799:7934] is HW BV 
f9.5_AM_BV <- (as.mcmc(PL_model9.5$Sol[,148:283])) #index [,148:283] is AM BV

#Resampling posterior to estimate correlation between BV from each environment
r.f9.5 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.f9.5[i] <- (cov(f9.5_AM_BV[i,], f9.5_HW_BV[i,])) /
    sqrt(var(f9.5_AM_BV[i,]) * var(f9.5_HW_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.f9.5), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.f9.5)) # 95% interval
P.r.f9.5 <- ifelse(r.f9.5 <0, 0, 1) 
(1-sum(P.r.f9.5)/2000)*2 # 


#Testing if different approaches to estimating the genetic correlation differs
P.diff8 <- ifelse(r.10.3 > r.f9.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff8)/2000)*2  # P =


    ## STEP 2 - Extracting r [individual models] including additive, dominance, maternal effects

#Genetic Correlation for models with Additive + Dominance + Maternal Effects
load(file="Routput/Ambient/AM_model10.1.RData")
load(file="Routput/Heated/HW_model10.1.RData")
AM_10.1_BV <- AM_model10.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_10.1_BV <- HW_model10.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.10.1 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.10.1[i] <- (cov(AM_10.1_BV[i,], HW_10.1_BV[i,])) /
    sqrt(var(AM_10.1_BV[i,]) * var(HW_10.1_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.10.1), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.10.1)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.10.1 <- ifelse(r.10.1 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.10.1)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 

# genetic variance for plasticity!
Va_plasticity10.1 <- numeric(2000)
for (i in 1:2000) {
  
  # this is the additive genetic variance for plasticity
  # it's the covariance of A and B divided by the square root of the product of their variances
  plasticity10.1 <- (AM_10.1_BV[i,] - HW_10.1_BV[i,]) 
  Va_plasticity10.1[i] <- var(plasticity10.1)
}
posterior.mode(as.mcmc(Va_plasticity10.1), na.rm=T) 
HPDinterval(as.mcmc(Va_plasticity10.1))
P_Va10.1 <- ifelse(Va_plasticity10.1 < 0, 0, 1)
(1-sum(P_Va10.1)/2000)*2

#Testing if including dominance changes anything to the genetic correlation (P value)
P.diff8b <- ifelse(r.10.1 > r.10.3, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff8b)/2000)*2  # P = **

#Plotting breeding values 
AM_10.1_BV <- colMeans(as.data.frame(AM_10.1_BV))
HW_10.1_BV <- colMeans(as.data.frame(HW_10.1_BV))
BV_stem <- bind_cols(HW_10.1_BV, AM_10.1_BV)
colnames(BV_stem)[1] <- "HW_stem_BV"
colnames(BV_stem)[2] <- "AM_stem_BV"
BV_stem %>% 
  ggplot(aes(x=AM_stem_BV, y=HW_stem_BV)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") + #1 to 1 line is dashed
  geom_smooth()

#





## Trait 9: Flowering Clusters (Poisson) ####

## STEP 1 - Testing if r [individual models] = r [full models]]

#Genetic Correlation for models with Additive + Maternal Effects
load(file="Routput/Ambient/AM_model11.3.RData")
load(file="Routput/Heated/HW_model11.3.RData")
AM_11.3_BV <- AM_model11.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_11.3_BV <- HW_model11.3$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
r.11.3 <- 1:2000 
for (i in 1:2000) {
  r.11.3[i] <- (cov(AM_11.3_BV[i,], HW_11.3_BV[i,])) /
    sqrt(var(AM_11.3_BV[i,]) * var(HW_11.3_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.11.3), na.rm=T) #point estimate of r
HPDinterval(as.mcmc(r.11.3)) # 95% interval
P.r.11.3 <- ifelse(r.11.3 <0, 0, 1) 
(1-sum(P.r.11.3)/2000)*2 # 


#Genetic correlation for full models
load(file="Routput/Full_Excl_Dom/Plasticity_9.5.RData")

#Extracting breeding values for each environment
f9.5_HW_BV <- (as.mcmc(PL_model9.5$Sol[,7799:7934])) # index [, 7799:7934] is HW BV 
f9.5_AM_BV <- (as.mcmc(PL_model9.5$Sol[,148:283])) #index [,148:283] is AM BV

#Resampling posterior to estimate correlation between BV from each environment
r.f9.5 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.f9.5[i] <- (cov(f9.5_AM_BV[i,], f9.5_HW_BV[i,])) /
    sqrt(var(f9.5_AM_BV[i,]) * var(f9.5_HW_BV[i,]))
}

#Testing if different approaches to estimating the genetic correlation differs
P.diff9 <- ifelse(r.11.3 > r.f9.5, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff9)/2000)*2  # P =


## STEP 2 - Extracting r [individual models] including additive, dominance, maternal effects

#Genetic Correlation for models with Additive + Dominance + Maternal Effects
load(file="Routput/Ambient/AM_model11.1-extend2.RData")
load(file="Routput/Heated/HW_model11.1-extend2.RData")
AM_11.1_BV <- AM_model11.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
HW_11.1_BV <- HW_model11.1$Sol[,142:277] #extracting breeding values of parents - parents are from 136 to 271 on pedigree file

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.11.1 <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.11.1[i] <- (cov(AM_11.1_BV[i,], HW_11.1_BV[i,])) /
    sqrt(var(AM_11.1_BV[i,]) * var(HW_11.1_BV[i,]))
}

#Testing if correlation is different from 0 
posterior.mode(as.mcmc(r.11.1), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.11.1)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.11.1 <- ifelse(r.11.1 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.11.1)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 

# genetic variance for plasticity!
Va_plasticity11.1 <- numeric(2000)
for (i in 1:2000) {
  
  # this is the additive genetic variance for plasticity
  # it's the covariance of A and B divided by the square root of the product of their variances
  plasticity11.1 <- (AM_11.1_BV[i,] - HW_11.1_BV[i,]) 
  Va_plasticity11.1[i] <- var(plasticity11.1)
}
posterior.mode(as.mcmc(Va_plasticity11.1), na.rm=T) 
HPDinterval(as.mcmc(Va_plasticity11.1))
P_Va11.1 <- ifelse(Va_plasticity11.1 < 0, 0, 1)
(1-sum(P_Va11.1)/2000)*2

#Testing if including dominance changes anything to the genetic correlation (P value)
P.diff2b <- ifelse(r.11.1 > r.11.3, 1, 0) # That is, knowing that the point estimate of r.Mod.2 > r.Mod.1, we call those whose r.Mod.2 < r.Mod.1, 1, else 0...
(sum(P.diff2)/2000)*2  # P = **

#Plotting breeding values 
AM_11.1_BV <- colMeans(as.data.frame(AM_11.1_BV))
HW_11.1_BV <- colMeans(as.data.frame(HW_11.1_BV))
BV_flwrclstr <- bind_cols(HW_11.1_BV, AM_11.1_BV)
colnames(BV_flwrclstr)[1] <- "HW_flwrclstr_BV"
colnames(BV_flwrclstr)[2] <- "AM_flwrclstr_BV"
BV_flwrclstr %>% 
  ggplot(aes(x=AM_flwrclstr_BV, y=HW_flwrclstr_BV)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") + #1 to 1 line is dashed
  geom_smooth()

#




## Trait 10: Genetic Correlation between Survival and Fecundity (Poisson) ####

## AMBIENT TREATMENT ##
load(file="Routput/Ambient/AM_model3.0.RData") # AM Survival
load(file="Routput/Ambient/AM_model4.1.RData") # AM Fecundity

AM_3.0_BV <- AM_model3.0$Sol[,13:7663] #extracting breeding values of entire pedigree
AM_4.1_BV <- AM_model4.1$Sol[,13:7663] #extracting breeding values of entire pedigree

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.AM_SxF <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.AM_SxF[i] <- (cov(AM_3.0_BV[i,], AM_4.1_BV[i,])) /
    sqrt(var(AM_3.0_BV[i,]) * var(AM_4.1_BV[i,]))
}

#Testing if correlation is different from 0 
mean(as.mcmc(r.AM_SxF), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.AM_SxF)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.AM_SxF <- ifelse(r.AM_SxF <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
p1 <- (1-sum(P.r.AM_SxF)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 
p1


## HEATED TREATMENT ##
load(file="Routput/Heated/HW_model3.0.RData") # HW Survival
load(file="Routput/Heated/HW_model4.1.RData") # HW Fecundity
HW_3.0_BV <- HW_model3.0$Sol[,13:7663] #extracting breeding values of entire pedigree
HW_4.1_BV <- HW_model4.1$Sol[,13:7663] #extracting breeding values of entire pedigree

#Pearson correlation between breeding values in AM and HW
#Resampling posterior to estimate correlation between BV from each environment
r.HW_SxF <- 1:2000 # set up empty vector to store each estimate of the correlation of the breeding values across environments
for (i in 1:2000) {
  
  # this is just a Pearson's correlation, spelled out in code
  # it's the covariance of A and B divided by the square root of the product of their variances
  r.HW_SxF[i] <- (cov(HW_3.0_BV[i,], HW_4.1_BV[i,])) /
    sqrt(var(HW_3.0_BV[i,]) * var(HW_4.1_BV[i,]))
}

#Testing if correlation is different from 0 
mean(as.mcmc(r.HW_SxF), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.HW_SxF)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.HW_SxF <- ifelse(r.HW_SxF <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
p2 <- (1-sum(P.r.HW_SxF)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 2000 estimates that were below 0, times 2. 
p2

## PHENOTYPIC CORRELATION ##

MCMC.ambient <- MCMC.2019 %>%
  filter(treatment=="A")
MCMC.heated <- MCMC.2019 %>%
  filter(treatment=="H")
heated.Survival <- MCMC.heated
heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
ambient.Survival <- MCMC.ambient
ambient.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)

phen.corrA2.1 <- cor.test(ambient.Survival$flower, ambient.Survival$seed_pods, method="pearson", exact=FALSE)
phen.corrA2.1 # Corr = 0.5112517 | Posterior 95% CI =  (0.4861976, 0.5354659)

phen.corrH2.1 <- cor.test(heated.Survival$flower, heated.Survival$seed_pods, test="pearson")
phen.corrH2.1 # Corr = 0.4564391 | Posterior 95% CI =  (0.4291918, 0.4828564)
