#### PROJECT: Brassica rapa h2 comparison between environments
#### PURPOSE: Models to estimate Va, Vd, Vm + cross-environment Va, Vd, Vm + and genetic correlations

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
lapply(ped, class)
lapply(MCMC.2019, class)

MCMC.ambient <- MCMC.2019 %>%
  filter(treatment=="A")

MCMC.heated <- MCMC.2019 %>%
  filter(treatment=="H")

MCMC.ambient$animal <- as.factor(MCMC.ambient$animal) #double check
MCMC.heated$animal <- as.factor(MCMC.heated$animal) #double check
lapply(ped, class)
lapply(MCMC.ambient, class)
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
#############   'BROAD-SENSE' GENETIC VARIANCE MODEL    ##############
#'####################################################################'#



## Model 1 - Lifetime Fitness ####
#Zero-Inflated Poisson distribution + 2 different priors
#Level 1 = Poisson process, Level 2 = Zero inflation


#Ambient
ambient.total <- MCMC.ambient
prior1.0 <- list(R = list(V = diag(2), nu = 0.002, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))
                 ) # Chi Square Prior with SD 1

AM_Vg1.0 <- MCMCglmm(seed_pods ~ trait - 1 + plot, 
                     random = ~idh(at.level(trait,1)):famID + idh(at.level(trait,2)):famID,
                     family = "zipoisson", rcov = ~idh(trait):units, data = ambient.total, prior = prior1.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

# with covariance structure between the Zero inflation and poisson process
AM_Vg1.1 <- MCMCglmm(seed_pods ~ trait - 1 + plot, 
                     random = ~us(at.level(trait,1)):famID + us(at.level(trait,2)):famID,
                     family = "zipoisson", rcov = ~idh(trait):units, data = ambient.total, prior = prior1.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(AM_Vg1.0, file="Routput/Vg/AM_Vg1.0.RData")
save(AM_Vg1.0, file="Routput/Vg/AM_Vg1.1.RData")

load(file="Routput/Vg/AM_Vg1.0.RData")
plot(AM_Vg1.1$VCV) #  good trace plots
summary(AM_Vg1.0) # 
effectiveSize(AM_Vg1.0$VCV) #  good effective sample sizes
autocorr.diag(AM_Vg1.0$VCV) #  no autocorrelation
heidel.diag(AM_Vg1.0$VCV) # good convergence 

#Latent Scale Vg
mean(AM_Vg1.0$VCV[, "at.level(trait, 1).famID"]) #Vg poisson process
mean(AM_Vg1.0$VCV[, "at.level(trait, 2).famID"]) #Vg zero inflation process
HPDinterval(AM_Vg1.0$VCV[, "at.level(trait, 1).famID"])
HPDinterval(AM_Vg1.0$VCV[, "at.level(trait, 2).famID"])

mean(AM_Vg1.0[["VCV"]][, "units"]) #Vr 
HPDinterval(AM_Vg1.0[["VCV"]][, "units"])

mean(AM_Vg1.0[["VCV"]][ , "famID"]/(AM_Vg1.0[["VCV"]][, "famID"] + AM_Vg1.0[["VCV"]][, "units"])) #H2
HPDinterval(AM_Vg1.0[["VCV"]][, "famID"]/(AM_Vg1.0[["VCV"]][, "famID"] + AM_Vg1.0[["VCV"]][, "units"]))



#Heated
heated.total <- MCMC.heated
prior1.0 <- list(R = list(V = diag(2), nu = 0.002, fix=2),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000), # Expanded Fisher prior for Poisson process
                          G2 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))) # Chi Square Prior with SD 1

HW_Vg1.0 <- MCMCglmm(seed_pods ~ trait - 1 + plot, 
                     random = ~idh(at.level(trait,1)):famID + idh(at.level(trait,2)):famID,
                     family = "zipoisson", rcov = ~idh(trait):units, data = heated.total, prior = prior1.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

#with covariance s tructure using us
HW_Vg1.1 <- MCMCglmm(seed_pods ~ trait - 1 + plot, 
                     random = ~us(at.level(trait,1)):famID + us(at.level(trait,2)):famID,
                     family = "zipoisson", rcov = ~us(trait):units, data = heated.total, prior = prior1.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

save(HW_Vg1.0, file="Routput/Vg/HW_Vg1.0.RData")
save(HW_Vg1.1, file="Routput/Vg/HW_Vg1.1.RData")

load(file="Routput/Vg/HW_Vg1.0.RData")
plot(HW_Vg1.0$VCV) # good trace plots
summary(HW_Vg1.0) # 
effectiveSize(HW_Vg1.0$VCV) # good effective sample sizes
autocorr.diag(HW_Vg1.0$VCV) # no autocorrelation
heidel.diag(HW_Vg1.0$VCV) # good convergence 

#Latent Scale Vg
mean(HW_Vg1.0$VCV[, "at.level(trait, 1).famID"]) #Vg poisson process
mean(HW_Vg1.0$VCV[, "at.level(trait, 2).famID"]) #Vg zero inflation process
HPDinterval(HW_Vg1.0$VCV[, "at.level(trait, 1).famID"])
HPDinterval(HW_Vg1.0$VCV[, "at.level(trait, 2).famID"])

mean(HW_Vg1.0[["VCV"]][, "units"]) #Vr 
HPDinterval(HW_Vg1.0[["VCV"]][, "units"])

mean(HW_Vg1.0[["VCV"]][ , "famID"]/(HW_Vg1.0[["VCV"]][, "famID"] + HW_Vg1.0[["VCV"]][, "units"])) #H2
HPDinterval(HW_Vg1.0[["VCV"]][, "famID"]/(HW_Vg1.0[["VCV"]][, "famID"] + HW_Vg1.0[["VCV"]][, "units"]))






## Model 2 - Survival ####

#Ambient
ambient.Survival <- MCMC.ambient
prior2.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

AM_Vg2.0 <- MCMCglmm(flower ~ plot, random = ~famID,
                     family = "threshold", data = ambient.Survival, prior = prior2.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)
save(AM_Vg2.0, file="Routput/Vg/AM_Vg2.0.RData")

load(file="Routput/Vg/AM_Vg2.0.RData")
plot(AM_Vg2.0$VCV) # good trace plots
summary(AM_Vg2.0) # 
effectiveSize(AM_Vg2.0$VCV) # good effective sample sizes
autocorr.diag(AM_Vg2.0$VCV) # no autocorrelation
heidel.diag(AM_Vg2.0$VCV) # good convergence 

#Latent Scale Vg
mean(AM_Vg2.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(AM_Vg2.0[["VCV"]][ , "famID"]) 

mean(AM_Vg2.0[["VCV"]][, "units"]) #Vr 
HPDinterval(AM_Vg2.0[["VCV"]][, "units"])

mean(AM_Vg2.0[["VCV"]][ , "famID"]/(AM_Vg2.0[["VCV"]][, "famID"] + AM_Vg2.0[["VCV"]][, "units"])) #H2
HPDinterval(AM_Vg2.0[["VCV"]][, "famID"]/(AM_Vg2.0[["VCV"]][, "famID"] + AM_Vg2.0[["VCV"]][, "units"]))

#Heated
heated.Survival <- MCMC.heated
prior2.0 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

HW_Vg2.0 <- MCMCglmm(flower ~ plot, random = ~famID,
                     family = "threshold", data = heated.Survival, prior = prior2.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc=TRUE)
save(HW_Vg2.0, file="Routput/Vg/HW_Vg2.0.RData")

load(file="Routput/Vg/HW_Vg2.0.RData")
plot(HW_Vg2.0$VCV) # good trace plots
summary(HW_Vg2.0) # 
effectiveSize(HW_Vg2.0$VCV) #  good effective sample size
autocorr.diag(HW_Vg2.0$VCV) #  no autocorrelation
heidel.diag(HW_Vg2.0$VCV) # good convergence 

#Latent Scale Vg
mean(HW_Vg2.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(HW_Vg2.0[["VCV"]][ , "famID"])

mean(HW_Vg2.0[["VCV"]][, "units"]) #Vr 
HPDinterval(HW_Vg2.0[["VCV"]][, "units"])

mean(HW_Vg2.0[["VCV"]][ , "famID"]/(HW_Vg2.0[["VCV"]][, "famID"] + HW_Vg2.0[["VCV"]][, "units"])) #H2
HPDinterval(HW_Vg2.0[["VCV"]][, "famID"]/(HW_Vg2.0[["VCV"]][, "famID"] + HW_Vg2.0[["VCV"]][, "units"]))
#










## Model 3 - Fecundity ####

#Ambient
ambient.Fecundity <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior3.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_Vg3.0 <- MCMCglmm(seed_pods ~ plot, random = ~famID,
                     family = "poisson", data = ambient.Fecundity, prior = prior3.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
save(AM_Vg3.0, file="Routput/Vg/AM_Vg3.0.RData")

load(file="Routput/Vg/AM_Vg3.0.RData")
plot(AM_Vg3.0$VCV) # good trace plots (though skewed a bit)
summary(AM_Vg3.0) # 
effectiveSize(AM_Vg3.0$VCV) #  moderate effective sample sizes
autocorr.diag(AM_Vg3.0$VCV) #  no autocorrelation
heidel.diag(AM_Vg3.0$VCV) # good convergence 

#Latent Scale Vg
mean(AM_Vg3.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(AM_Vg3.0[["VCV"]][ , "famID"]) 

mean(AM_Vg3.0[["VCV"]][, "units"]) #Vr 
HPDinterval(AM_Vg3.0[["VCV"]][, "units"])

mean(AM_Vg3.0[["VCV"]][ , "famID"]/(AM_Vg3.0[["VCV"]][, "famID"] + AM_Vg3.0[["VCV"]][, "units"])) #H2
HPDinterval(AM_Vg3.0[["VCV"]][, "famID"]/(AM_Vg3.0[["VCV"]][, "famID"] + AM_Vg3.0[["VCV"]][, "units"]))

#Heated
heated.Fecundity <- MCMC.heated %>% filter(germ==1 & flower==1)
prior3.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_Vg3.0 <- MCMCglmm(seed_pods ~ plot, random = ~famID,
                     family = "poisson", data = heated.Fecundity, prior = prior3.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
save(HW_Vg3.0, file="Routput/Vg/HW_Vg3.0.RData")

load(file="Routput/Vg/HW_Vg3.0.RData")
plot(HW_Vg3.0$VCV) # good trace plots
summary(HW_Vg3.0) # 
effectiveSize(HW_Vg3.0$VCV) # good effective sHWple size
autocorr.diag(HW_Vg3.0$VCV) # no autocorrelation
heidel.diag(HW_Vg3.0$VCV) # convergence success

#Latent Scale Vg
mean(HW_Vg3.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(HW_Vg3.0[["VCV"]][ , "famID"])

mean(HW_Vg3.0[["VCV"]][, "units"]) #Vr 
HPDinterval(HW_Vg3.0[["VCV"]][, "units"])

mean(HW_Vg3.0[["VCV"]][ , "famID"]/(HW_Vg3.0[["VCV"]][, "famID"] + HW_Vg3.0[["VCV"]][, "units"])) #H2
HPDinterval(HW_Vg3.0[["VCV"]][, "famID"]/(HW_Vg3.0[["VCV"]][, "famID"] + HW_Vg3.0[["VCV"]][, "units"]))
#













## Model 4 - Leaf Number ####

#Ambient
ambient.leaf <- MCMC.ambient %>% filter(germ==1)
prior4.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_Vg4.0 <- MCMCglmm(leaf ~ plot, random = ~famID,
                     family = "gaussian", data = ambient.leaf, prior = prior4.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
save(AM_Vg4.0, file="Routput/Vg/AM_Vg4.0.RData")

load(file="Routput/Vg/AM_Vg4.0.RData")
plot(AM_Vg4.0$VCV) # good trace plots 
summary(AM_Vg4.0) # 
effectiveSize(AM_Vg4.0$VCV) # good effective sizes
autocorr.diag(AM_Vg4.0$VCV) # no autocorrelation
heidel.diag(AM_Vg4.0$VCV) # good convergence 

#Latent Scale Vg
mean(AM_Vg4.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(AM_Vg4.0[["VCV"]][ , "famID"]) 

mean(AM_Vg4.0[["VCV"]][, "units"]) #Vr 
HPDinterval(AM_Vg4.0[["VCV"]][, "units"])

mean(AM_Vg4.0[["VCV"]][ , "famID"]/(AM_Vg4.0[["VCV"]][, "famID"] + AM_Vg4.0[["VCV"]][, "units"])) #H2
HPDinterval(AM_Vg4.0[["VCV"]][, "famID"]/(AM_Vg4.0[["VCV"]][, "famID"] + AM_Vg4.0[["VCV"]][, "units"]))

#Heated
heated.leaf <- MCMC.heated %>% filter(germ==1)
prior4.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_Vg4.0 <- MCMCglmm(leaf ~ plot, random = ~famID,
                     family = "gaussian", data = heated.leaf, prior = prior4.0, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
save(HW_Vg4.0, file="Routput/Vg/HW_Vg4.0.RData")

load(file="Routput/Vg/HW_Vg4.0.RData")
plot(HW_Vg4.0$VCV) #good trace pots
summary(HW_Vg4.0) # 
effectiveSize(HW_Vg4.0$VCV) # good effective sample sizes
autocorr.diag(HW_Vg4.0$VCV) # no autocorrelation
heidel.diag(HW_Vg4.0$VCV) # good convergence 

#Latent Scale Vg
mean(HW_Vg4.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(HW_Vg4.0[["VCV"]][ , "famID"]) 

mean(HW_Vg4.0[["VCV"]][, "units"]) #Vr 
HPDinterval(HW_Vg4.0[["VCV"]][, "units"])

mean(HW_Vg4.0[["VCV"]][ , "famID"]/(HW_Vg4.0[["VCV"]][, "famID"] + HW_Vg4.0[["VCV"]][, "units"])) #H2
HPDinterval(HW_Vg4.0[["VCV"]][, "famID"]/(HW_Vg4.0[["VCV"]][, "famID"] + HW_Vg4.0[["VCV"]][, "units"]))
#








## Model 5 - Height ####

#Ambient
ambient.height <- MCMC.ambient %>% filter(germ==1, !height==0)
prior5.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_Vg5.0 <- MCMCglmm(height ~ plot, random = ~famID, 
                     family = "gaussian", data = ambient.height, prior = prior5.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
save(AM_Vg5.0, file="Routput/Vg/AM_Vg5.0.RData")

load(file="Routput/Vg/AM_Vg5.0.RData")
plot(AM_Vg5.0$VCV) #good trace plots
summary(AM_Vg5.0) # 
effectiveSize(AM_Vg5.0$VCV) # good effective sample sizes
autocorr.diag(AM_Vg5.0$VCV) # no autocorrelation
heidel.diag(AM_Vg5.0$VCV) # good convergence 

#Latent/Data Scale Vg
mean(AM_Vg5.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(AM_Vg5.0[["VCV"]][ , "famID"])

mean(AM_Vg5.0[["VCV"]][, "units"]) #Vr 
HPDinterval(AM_Vg5.0[["VCV"]][, "units"])

mean(AM_Vg5.0[["VCV"]][ , "famID"]/(AM_Vg5.0[["VCV"]][ , "famID"] + AM_Vg5.0[["VCV"]][, "units"])) #H2
HPDinterval(AM_Vg5.0[["VCV"]][, "famID"]/(AM_Vg5.0[["VCV"]][ , "famID"] + AM_Vg5.0[["VCV"]][, "units"]))

#Heated
heated.height <- MCMC.heated %>% filter(germ==1, !height==0)
prior5.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_Vg5.0 <- MCMCglmm(height ~ plot, random = ~famID,
                     family = "gaussian", data = heated.height, prior = prior5.0, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
save(HW_Vg5.0, file="Routput/Vg/HW_Vg5.0.RData")

load(file="Routput/Vg/HW_Vg5.0.RData")
plot(HW_Vg5.0$VCV) # good trace plots
summary(HW_Vg5.0) # 
effectiveSize(HW_Vg5.0$VCV) # good effective sample sizes
autocorr.diag(HW_Vg5.0$VCV) # no autocorrelation
heidel.diag(HW_Vg5.0$VCV) # good convergence issues

#Latent Scale Vg
mean(HW_Vg5.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(HW_Vg5.0[["VCV"]][ , "famID"]) 

mean(HW_Vg5.0[["VCV"]][, "units"]) #Vr 
HPDinterval(HW_Vg5.0[["VCV"]][, "units"])

mean(HW_Vg5.0[["VCV"]][ , "famID"]/(HW_Vg5.0[["VCV"]][ , "famID"] + HW_Vg5.0[["VCV"]][, "units"])) #H2
HPDinterval(HW_Vg5.0[["VCV"]][, "famID"]/(HW_Vg5.0[["VCV"]][ , "famID"] + HW_Vg5.0[["VCV"]][, "units"]))

#








## Model 6 - Stem Diameter ####

#Ambient
ambient.stem <- MCMC.ambient %>% filter(germ==1, !stem_diam==0)
prior6.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_Vg6.0 <- MCMCglmm(stem_diam ~ plot, random = ~famID, 
                     family = "gaussian", data = ambient.stem, prior = prior6.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
save(AM_Vg6.0, file="Routput/Vg/AM_Vg6.0.RData")

load(file="Routput/Vg/AM_Vg6.0.RData")
plot(AM_Vg6.0$VCV) # good trace plots
summary(AM_Vg6.0) # 
effectiveSize(AM_Vg6.0$VCV) # good effective sample sizes
autocorr.diag(AM_Vg6.0$VCV) # no autocorrelation
heidel.diag(AM_Vg6.0$VCV) # good convergence 

#Latent/Data Scale Vg
mean(AM_Vg6.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(AM_Vg6.0[["VCV"]][ , "famID"])

mean(AM_Vg6.0[["VCV"]][, "units"]) #Vr 
HPDinterval(AM_Vg6.0[["VCV"]][, "units"])

mean(AM_Vg6.0[["VCV"]][ , "famID"]/(AM_Vg6.0[["VCV"]][ , "famID"] + AM_Vg6.0[["VCV"]][, "units"])) #H2
HPDinterval(AM_Vg6.0[["VCV"]][, "famID"]/(AM_Vg6.0[["VCV"]][ , "famID"] + AM_Vg6.0[["VCV"]][, "units"]))

#Heated
heated.stem <- MCMC.heated %>% filter(germ==1, !stem_diam==0)
prior6.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_Vg6.0 <- MCMCglmm(stem_diam ~ plot, random = ~famID,
                     family = "gaussian", data = heated.stem, prior = prior6.0, #gaussian distribution
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
save(HW_Vg6.0, file="Routput/Vg/HW_Vg6.0.RData")

load(file="Routput/Vg/HW_Vg6.0.RData")
plot(HW_Vg6.0$VCV) # good trace plots
summary(HW_Vg6.0) # 
effectiveSize(HW_Vg6.0$VCV) # good effective sample sizes
autocorr.diag(HW_Vg6.0$VCV) # no autocorrelation
heidel.diag(HW_Vg6.0$VCV) # good convergence

#Latent Scale Vg
mean(HW_Vg6.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(HW_Vg6.0[["VCV"]][ , "famID"]) 

mean(HW_Vg6.0[["VCV"]][, "units"]) #Vr 
HPDinterval(HW_Vg6.0[["VCV"]][, "units"])

mean(HW_Vg6.0[["VCV"]][ , "famID"]/(HW_Vg6.0[["VCV"]][ , "famID"] + HW_Vg6.0[["VCV"]][, "units"])) #H2
HPDinterval(HW_Vg6.0[["VCV"]][, "famID"]/(HW_Vg6.0[["VCV"]][ , "famID"] + HW_Vg6.0[["VCV"]][, "units"]))

#







## Model 7 - Flowering clusters ####

#Ambient
ambient.flwrclstr <- MCMC.ambient %>% filter(germ==1 & flower==1)
prior7.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

AM_Vg7.0 <- MCMCglmm(flwr_clstr ~ plot, random = ~famID,
                     family = "poisson", data = ambient.flwrclstr, prior = prior7.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
#nitt = 8300000, thin = 1800, burnin = 200000,
save(AM_Vg7.0, file="Routput/Vg/AM_Vg7.0.RData")

load(file="Routput/Vg/AM_Vg7.0.RData")
plot(AM_Vg7.0$VCV) # good trace plots
summary(AM_Vg7.0) # 
effectiveSize(AM_Vg7.0$VCV) # good effective sample sizes
autocorr.diag(AM_Vg7.0$VCV) # no autocorrelation
heidel.diag(AM_Vg7.0$VCV) # good convergence 

#Latent Scale Vg
mean(AM_Vg7.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(AM_Vg7.0[["VCV"]][ , "famID"]) 

mean(AM_Vg7.0[["VCV"]][, "units"]) #Vr 
HPDinterval(AM_Vg7.0[["VCV"]][, "units"])

mean(AM_Vg7.0[["VCV"]][ , "famID"]/(AM_Vg7.0[["VCV"]][, "famID"] + AM_Vg7.0[["VCV"]][, "units"])) #H2
HPDinterval(AM_Vg7.0[["VCV"]][, "famID"]/(AM_Vg7.0[["VCV"]][, "famID"] + AM_Vg7.0[["VCV"]][, "units"]))


#Model 7.0 has nitt = 2.1 mil, thin = 1k, burnin = 100k
#Model 7.1 has nitt = 8.5 mil, thin = 4k, burnin = 400k
#Model 7.2 has nitt = 8.5 mil, thin = 3.2k, burnin = 500k
#Model 7.3 has nitt = 8.3 mil, thin = 1.8k, burnin = 200k

#Model 7.2 has the lowest autocorrelation (0.13.. which is OK)

#Heated
heated.flwrclstr <- MCMC.heated %>% filter(germ==1 & flower==1)
prior7.0 <- list(R = list(V = 1, nu = 0.002),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

HW_Vg7.0 <- MCMCglmm(flwr_clstr ~ plot, random = ~famID,
                     family = "poisson", data = heated.flwrclstr, prior = prior7.0,
                     nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
save(HW_Vg7.0, file="Routput/Vg/HW_Vg7.0.RData")

load(file="Routput/Vg/HW_Vg7.0.RData")
plot(HW_Vg7.0$VCV) # good trace plots
summary(HW_Vg7.0) # 
effectiveSize(HW_Vg7.0$VCV) # good effective sample sizes
autocorr.diag(HW_Vg7.0$VCV) # no autocorrelation
heidel.diag(HW_Vg7.0$VCV) # good convergence

#Latent Scale Vg
mean(HW_Vg7.0[["VCV"]][ , "famID"]) #Vg
HPDinterval(HW_Vg7.0[["VCV"]][ , "famID"])

mean(HW_Vg7.0[["VCV"]][, "units"]) #Vr 
HPDinterval(HW_Vg7.0[["VCV"]][, "units"])

mean(HW_Vg7.0[["VCV"]][ , "famID"]/(HW_Vg7.0[["VCV"]][, "famID"] + HW_Vg7.0[["VCV"]][, "units"])) #H2
HPDinterval(HW_Vg7.0[["VCV"]][, "famID"]/(HW_Vg7.0[["VCV"]][, "famID"] + HW_Vg7.0[["VCV"]][, "units"]))
#

