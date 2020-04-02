#### PROJECT: Brassica rapa GxE Study (Data collected by Cameron So 2019 at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Estimate fixed effect of treatment per trait
#### AUTHOR: Cameron So
#### DATE LAST MODIFIED: 2020/02/02

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Necessary Packages
library(MCMCglmm)
library(nadiv)
library(dplyr)
library(readr)

#Importing Data using readr from tidyverse (Base R confuses the class for certain vectors). 
# d for double, f for factor, i for integer
col_types_list2 <- cols_only(posID = "d", individual = "d", plot = col_factor(levels=c(1:12)),
                             animal = "f",
                             patID = "f", matID = "f", famID = col_factor(levels=c(1:62)),
                             treatment = col_factor(levels=c("A", "H")),
                             germ = "d", flower = "d",
                             seed = "d", leaf = "d", flwr_clstr = "d", seed_pods = "d",
                             height = "d", stem_diam = "d")
MCMC.2019 <- read_csv("Rdata/MCMC_2019_cleaned.csv", col_names = TRUE, na = "NA",
                      col_types=col_types_list2)
MCMC.2019$animal <- as.factor(MCMC.2019$animal)
ped <- read.csv("Rdata/heatarrays_animal_pedigree.csv") #note don't use readr to import the csv file.
for (x in 1:3) ped[, x] <- as.factor(ped[, x])

total.2019 <- as.data.frame(MCMC.2019)
ped <- as.data.frame(ped)

Ainv <- inverseA(ped[, 1:3])$Ainv
Dinv <- makeD(ped[, 1:3])$Dinv
total.2019$animalDom <- total.2019$animal

#'####################################################################'#
   ###########      MCMCglmm ANALYSIS FOR TRT EFFECT      ############
#'####################################################################'#



### No. 2 Fecundity of Flowering Plants - Univariate | Poisson ####
#Note: For Fecundity, we only include plants that survived to reach flowering
#We use an uninformative prior that divides the total phenotypic variance by the number of random effects included
#Variance in Total Fitness = squared standard deviation
#We treat plot as a random effect as from previous analysis (01_Data_exploration), we observe that there is a plot effect on this trait... There is likely some sampling error caused by plot due to the placement of the physical plots at KSR with respect to the tree line, etc so we want to give less statistical weight to the effect of plot. We also want to estimate the effect of plot 

#Calculating Variance for Fecundity for the Non-informative Prior (equal variance)
total.2019 %>%
  filter(germ==1 & flower==1) %>%
  summarise(mean=mean(seed_pods), sd=sd(seed_pods), var=(sd)^2) #var = 613.3011 *******


#Producing filtered data set for Fecundity
total.Fecundity <- total.2019 %>% filter(germ==1 & flower==1)

#From Committee Report #2: Fecundity follows a Poisson distribution

#Fecundity Total Model - Additive Effects Only
load(file="Routput/MCMC_Total_Model_3.5.RData")
plot(model3.5$VCV)
summary(model3.5)
autocorr.diag(model3.5)


#### No. 2 Fecundity of Flowering Plants - Univariate | Poisson ####

total.Fecundity <- total.2019 %>% filter(germ==1 & flower==1)

#var = 613.3011

prior2.1 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002),
                          G4 = list(V = 1, nu = 0.002),
                          G5 = list(V = 1, nu = 0.002)))

prior2.2 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002),
                          G4 = list(V = 1, nu = 0.002)))

prior2.3 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002), 
                          G3 = list(V = 1, nu = 0.002)))

prior2.4 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002),
                          G2 = list(V = 1, nu = 0.002)))

prior2.5 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 0.002)))


model2.1 <- MCMCglmm(seed_pods ~ treatment, random = ~animal + animalDom + matID + plot + posID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "poisson", data = total.Fecundity, prior = prior2.1, #poission distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE)

model2.2 <- MCMCglmm(seed_pods ~ treatment, random = ~animal + animalDom + matID + posID,
                     ginverse = list(animal= Ainv, animalDom= Dinv),
                     family = "poisson", data = total.Fecundity, prior = prior2.2,
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE)

model2.3 <- MCMCglmm(seed_pods ~ treatment, random = ~animal + animalDom + matID,
                     ginverse = list(animal= Ainv, animalDom= Dinv),
                     family = "poisson", data = total.Fecundity, prior = prior2.3,
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE)

model2.4 <- MCMCglmm(seed_pods ~ treatment, random = ~animal + animalDom,
                     ginverse = list(animal= Ainv, animalDom= Dinv),
                     family = "poisson", data = total.Fecundity, prior = prior2.4,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

model2.5 <- MCMCglmm(seed_pods ~ treatment, random = ~animal,
                     ginverse = list(animal= Ainv),
                     family = "poisson", data = total.Fecundity, prior = prior2.5,
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE)

model2.1_AC.Sol <- autocorr(model2.1$Sol)
model2.1_AC.VCV <- autocorr(model2.1$VCV)
model2.2_AC.Sol <- autocorr(model2.2$Sol)
model2.2_AC.VCV <- autocorr(model2.2$VCV)
model2.3_AC.Sol <- autocorr(model2.3$Sol)
model2.3_AC.VCV <- autocorr(model2.3$VCV)
model2.4_AC.Sol <- autocorr(model2.4$Sol)
model2.4_AC.VCV <- autocorr(model2.4$VCV)
model2.5_AC.Sol <- autocorr(model2.5$Sol)
model2.5_AC.VCV <- autocorr(model2.5$VCV)

#Note: The Total Model includes all data, and now the treatment is set as a fixed effect

#Fecundity Total Model 2.3 - Additive, Dominant, and Maternal Effects
#Model was re-run because of autocorrelation
load(file="Routput/MCMC_Total_Model_2.3.RData")
plot(model2.3$VCV) #terrible trace plot for additive and dominance
summary(model2.3) #low effective sample size for additive and dominance
autocorr.diag(model2.3$VCV) #high autocorrelation
heidel.diag(model2.3$VCV) # convergence success

#Fecundity Total Model 2.4 - Additive and Dominant Effects Only
#Model was re-run because of autocorrelation
load(file="Routput/MCMC_Total_Model_2.4.RData")
plot(model2.4$VCV) #terrible trace plot
summary(model2.4) #Half effective sample size obtained from initial 2000
autocorr.diag(model2.4$VCV) #high autocorrelation
heidel.diag(model2.4$VCV) # convergence success

#Fecundity Total Model 2.5 - Additive Effects Only
#Another model 2.5b was run to see the difference in autocorrelation, etc
load(file="Routput/MCMC_Total_Model_2.5.RData")
plot(model2.5$VCV) #Good trace plot
summary(model2.5) #Good effective sample size
autocorr.diag(model2.5$VCV) #low autocorrelation
heidel.diag(model2.5$VCV) # convergence success

#Another model 2.5b was run to see the difference in autocorrelation, etc if standard uninformative priors
load(file="Routput/MCMC_Total_Model_2.5b.RData")
plot(model2.5b$VCV) #Good trace plot
summary(model2.5b) #Good effective sample size
autocorr.diag(model2.5b$VCV) #autocorrelation
heidel.diag(model2.5b$VCV) # convergence success



#### No. 3 Survival to Flowering - Univariate | Bernoulli ####

total.Survival <- total.2019

prior3.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G5 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

prior3.2 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G4 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

prior3.3 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1), 
                          G3 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

prior3.4 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

prior3.5 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))



model3.1 <- MCMCglmm(flower ~ treatment, random = ~animal + animalDom + matID + plot + posID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = total.Survival, prior = prior3.1, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

model3.2 <- MCMCglmm(flower ~ treatment, random = ~animal + animalDom + matID + posID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = total.Survival, prior = prior3.2, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)

model3.3 <- MCMCglmm(flower ~ treatment, random = ~animal + animalDom + matID, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = total.Survival, prior = prior3.3, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)

model3.4 <- MCMCglmm(flower ~ treatment, random = ~animal + animalDom, 
                     ginverse = list(animal = Ainv, animalDom = Dinv),
                     family = "threshold", data = total.Survival, prior = prior3.4, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)

model3.5 <- MCMCglmm(flower ~ treatment, random = ~animal,
                     ginverse = list(animal = Ainv),
                     family = "threshold", data = total.Survival, prior = prior3.5, #Bernoulli distribution
                     nitt = 1100000, thin = 500, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)

model3.1_AC.Sol <- autocorr(model3.1$Sol)
model3.1_AC.VCV <- autocorr(model3.1$VCV)
model3.2_AC.Sol <- autocorr(model3.2$Sol)
model3.2_AC.VCV <- autocorr(model3.2$VCV)
model3.3_AC.Sol <- autocorr(model3.3$Sol)
model3.3_AC.VCV <- autocorr(model3.3$VCV)
model3.4_AC.Sol <- autocorr(model3.4$Sol)
model3.4_AC.VCV <- autocorr(model3.4$VCV)
model3.5_AC.Sol <- autocorr(model3.5$Sol)
model3.5_AC.VCV <- autocorr(model3.5$VCV)

#Survival Total Model 3.3 - Additive, Dominance and Maternal Effects Only
load(file="Routput/MCMC_Total_Model_3.3.RData")
plot(model3.3$VCV) #good trace plot
summary(model3.3) #good effective sample size, significant treatment effect
autocorr.diag(model3.3$VCV) # no autocorrelation
heidel.diag(model3.3$VCV) # convergence success

#Survival Total Model 3.4 - Additive and Dominance Effects Only
load(file="Routput/MCMC_Total_Model_3.4.RData")
plot(model3.4$VCV) #not great trace plot
summary(model3.4) #good effective sample size, significant treatment effect
autocorr.diag(model3.4$VCV) # no autocorrelation
heidel.diag(model3.4$VCV) # convergence success

#Survival Total Model 3.5 - Additive Effects Only
load(file="Routput/MCMC_Total_Model_3.5.RData")
plot(model3.5$VCV) #good trace plot
summary(model3.5) #good effective sample size, significant treatment effect
autocorr.diag(model3.5$VCV) #no autocorrelation
heidel.diag(model3.5$VCV) # convergence success
