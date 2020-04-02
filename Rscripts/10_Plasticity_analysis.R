#### PROJECT: Brassica rapa GxE Study (Data collected by Cameron So 2019 at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Figures and models on plasticity for thesis and publication
#### AUTHOR: Cameron So
#### DATE LAST MODIFIED: 2020/02/02

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Necessary Packages
library(tidyverse)

#Importing data frame 'dat' from the cleaning process of 01_Data_exploration
load(file="Routput/heatarrays_data_explore.RData")

#'####################################################################'#
################      PLASTICITY-RELATED FIGURES      ################
#'####################################################################'#

#GxE Reaction Norms

#Lifetime Fitness
fam.absolute <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(seed_pods))
fam.absolute %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Absolute Fitness") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 40, 10)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Survival to Flowering
fam.surv <- dat %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(flower), n=n()) %>%
  mutate(fam.mean=(sum/n), se.p=(sqrt(fam.mean*(1-fam.mean))/n))
fam.surv %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Survival to Flowering") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 0.7, 0.1)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  ggtitle("A)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Fecundity
fam.fecund <- dat %>%
  filter(germ==1, flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(seed_pods))
fam.fecund %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Fecundity") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 60, 10)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


#Germination Success
fam.germ <- dat %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(germ), n=n()) %>%
  mutate(fam.mean=(sum/n), se.p=(sqrt(fam.mean*(1-fam.mean))/n))
fam.germ %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Germination Success") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 0.8, 0.1)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
       axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Flowering Success
fam.flower <- dat %>%
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(flower), n=n()) %>%
  mutate(fam.mean=(sum/n), se.p=(sqrt(fam.mean*(1-fam.mean))/n))
fam.flower %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Flowering Success") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0.5, 1, 0.1)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


#Leaf Number
fam.leaf <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(leaf))
fam.leaf %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Leaf Number") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 8, 2)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Height
fam.height <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(height))
fam.height %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Height (cm)") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 40, 5)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Flowering Clusters Number
fam.flwrclstr <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(flwr_clstr))
fam.flwrclstr %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Flowering Cluster Number") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 4, 1)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Germination Time by Census
fam.germtime <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(germ_census))
fam.germtime %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Germination Time by Census") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 6, 1)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Flowering Time by Census
fam.flwrtime <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(flwr_census))
fam.flwrtime %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Family Flowering Time by Census") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 6, 1)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Necessary Packages
library(tidyverse)
library(MCMCglmm)
library(MuMIn)
library(logitnorm)
library(nadiv)

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

lapply(ped, class)
lapply(MCMC.heated, class)

#'####################################################################'#
###############      PLASTICITY MODELS ON TRAITS      ################
#'####################################################################'#

#We adapt the framework in Arnold et al. 2019 for studying plasticity via reaction norms
#This is more exlusively mentioned in Nussey et al. 2007

#Covariance in residuals is not identifiable since only one observation has one treatment applied
#Especially in binary models, the residual variance can't be identified.. so fix = 1

## Trait 1: Survival Success (Binary) ####
load(file="Routput/Plasticity/Plasticity_1.1.RData")

plastic.survival <- total.2019
prior1.1 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25))),
                          G2 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25)))))

PL_model1.1 <- MCMCglmm(flower ~ treatment + plot - 1, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov = ~units,
                        family = "threshold", data = plastic.survival, prior = prior1.1, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

#I should probably attempt plot as a random effect, NESTED underneath treatment
#This probably looks like + idh(treatment):plot with a prior of V=diag(2), nu = 15, alpha.mu=c(0,0), alpha.V = diag(c(1.25,1.25))..
#Although, based on the course notes, I'm not sure if it is diag(6) based on the number of levels in the plot

plot(PL_model1.1$VCV) # OK trace plots.. extending past 0 though
summary(PL_model1.1) # good effective sample size
autocorr.diag(PL_model1.1$VCV) # no autocorrelation
heidel.diag(PL_model1.1$VCV) # convergence success


## Trait 1 without maternal effects
load(file="Routput/Plasticity/Plasticity_1.2.RData")

plastic.survival <- total.2019
prior1.2 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25)))))

PL_model1.2 <- MCMCglmm(flower ~ treatment + plot, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov = ~units,
                        family = "threshold", data = plastic.survival, prior = prior1.2, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(PL_model1.2, file="Routput/Plasticity/Plasticity_1.2.RData")

plot(PL_model1.2$VCV) # good trace plot
summary(PL_model1.2) # good effective sample size
autocorr.diag(PL_model1.2$VCV) # no autocorrelation
heidel.diag(PL_model1.2$VCV) # convergence success


gen.corrPL1.2 <-PL_model1.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model1.2$VCV[,'treatmentA:treatmentA.animal']*PL_model1.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL1.2) #Post Mean = 0.520711
posterior.mode(gen.corrPL1.2) #Post Mode = 0.6079715
HPDinterval(gen.corrPL1.2) # Posterior 95% CI = (0.2680, 0.7995)




















## Trait 2: Fecundity (Poisson) ####
load(file="Routput/Plasticity/Plasticity_2.1.RData")

plastic.fecundity <- total.2019 %>% filter(germ==1 & flower ==1)
prior2.1 <- list(R = list(V = diag(1),fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model2.1 <- MCMCglmm(seed_pods ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.fecundity, prior = prior2.1, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model2.1$VCV) # good trace plots
summary(PL_model2.1) # good effective sample size
autocorr.diag(PL_model2.1$VCV) # no autocorrelation
heidel.diag(PL_model2.1$VCV) # convergence fail for matID

## Trait 2 without maternal effects
load(file="Routput/Plasticity/Plasticity_2.2.RData")

plastic.fecundity <- total.2019 %>% filter(germ==1 & flower ==1)
prior2.2 <- list(R = list(V = diag(1),fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model2.2 <- MCMCglmm(seed_pods ~ treatment + plot, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.fecundity, prior = prior2.2, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model2.2$VCV) # good trace plots
summary(PL_model2.2) # good effective sample size
autocorr.diag(PL_model2.2$VCV) # no autocorrelation
heidel.diag(PL_model2.2$VCV) # convergence success

gen.corrPL2.2 <-PL_model2.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model2.2$VCV[,'treatmentA:treatmentA.animal']*PL_model2.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL2.2) #Post Mean = 0.01773589
posterior.mode(gen.corrPL2.2) #Post Mode = 0.01640536
HPDinterval(gen.corrPL2.2) # Posterior 95% CI = (-0.0413328, 0.07889832)










## Trait 3: Germination Success (Binary) ####
load(file="Routput/Plasticity/Plasticity_3.1.RData")

plastic.Germination <- total.2019
prior3.1 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25))),
                          G2 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25)))))

PL_model3.1 <- MCMCglmm(germ ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "threshold", data = plastic.Germination, prior = prior3.1, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model3.1$VCV) # good trace plots
summary(PL_model3.1) # good effective sample size
autocorr.diag(PL_model3.1$VCV) # no autocorrelation
heidel.diag(PL_model3.1$VCV) # convergence success

## Trait 3 without maternal effects
load(file="Routput/Plasticity/Plasticity_3.2.RData")

plastic.Germination <- total.2019
prior3.2 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25)))))

PL_model3.2 <- MCMCglmm(germ ~ treatment + plot, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "threshold", data = plastic.Germination, prior = prior3.2, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model3.2$VCV) # good trace plots
summary(PL_model3.2) # good effective sample size
autocorr.diag(PL_model3.2$VCV) # no autocorrelation
heidel.diag(PL_model3.2$VCV) # convergence success

gen.corrPL3.2 <-PL_model3.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model3.2$VCV[,'treatmentA:treatmentA.animal']*PL_model3.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL3.2) #Post Mean = 0.4559276
posterior.mode(gen.corrPL3.2) #Post Mode = 0.5192249
HPDinterval(gen.corrPL3.2) # Posterior 95% CI = (0.166213, 0.7288745)












## Trait 4: Flowering Success (Binary) ####
load(file="Routput/Plasticity/Plasticity_4.1.RData")

plastic.Flowering <- total.2019 %>% filter(germ==1)
prior4.1 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                   G = list(G1 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25))),
                            G2 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25)))))

PL_model4.1 <- MCMCglmm(flower ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "threshold", data = plastic.Flowering, prior = prior4.1, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model4.1$VCV) # OK trace plots
summary(PL_model4.1) # good effective sample size
autocorr.diag(PL_model4.1$VCV) # no autocorrelation
heidel.diag(PL_model4.1$VCV) # convergence fail

## Trait 4: Flowering Success (Binary) ####
load(file="Routput/Plasticity/Plasticity_4.2.RData")

plastic.Flowering <- total.2019 %>% filter(germ==1)
prior4.2 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25)))))

PL_model4.2 <- MCMCglmm(flower ~ treatment + plot, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "threshold", data = plastic.Flowering, prior = prior4.2, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model4.2$VCV) # OK trace plots not great
summary(PL_model4.2) # good effective sample size
autocorr.diag(PL_model4.2$VCV) # no autocorrelation
heidel.diag(PL_model4.2$VCV) # convergence fail

gen.corrPL4.2 <-PL_model4.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model4.2$VCV[,'treatmentA:treatmentA.animal']*PL_model4.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL4.2) #Post Mean = 0.0798662
posterior.mode(gen.corrPL4.2) #Post Mode = 0.07684446
HPDinterval(gen.corrPL4.2) # Posterior 95% CI = (-0.3531059, 0.5405406)









## Trait 5: Seed Maturation Success (Binary) ####
load(file="Routput/Plasticity/Plasticity_5.1.RData")

plastic.Seed <- total.2019 %>% filter(germ==1 & flower==1)
prior5.1 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                     G = list(G1 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25))),
                              G2 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25)))))

PL_model5.1 <- MCMCglmm(seed ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "threshold", data = plastic.Seed, prior = prior5.1, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model5.1$VCV) # poor trace plots
summary(PL_model5.1) # good effective sample size
autocorr.diag(PL_model5.1$VCV) #no autocorrelation
heidel.diag(PL_model5.1$VCV) # convergence fail

## Trait 5 without maternal effects
load(file="Routput/Plasticity/Plasticity_5.2.RData")

plastic.Seed <- total.2019 %>% filter(germ==1 & flower==1)
prior5.2 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25)))))

PL_model5.2 <- MCMCglmm(seed ~ treatment + plot, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "threshold", data = plastic.Seed, prior = prior5.2, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model5.2$VCV) # OK trace plots
summary(PL_model5.2) # good effective sample size
autocorr.diag(PL_model5.2$VCV) # no autocorrelation
heidel.diag(PL_model5.2$VCV) # convergence fail

gen.corrPL5.2 <-PL_model5.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model5.2$VCV[,'treatmentA:treatmentA.animal']*PL_model5.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL5.2) #Post Mean = 0.008757276
posterior.mode(gen.corrPL5.2) #Post Mode = 0.05864863
HPDinterval(gen.corrPL5.2) # Posterior 95% CI = (-0.4762743, 0.4903583)









## Trait 6: Leaf Number (Gaussian) ####
load(file="Routput/Plasticity/Plasticity_6.1.RData")

plastic.leaf <- total.2019 %>% filter(germ==1)
prior6.1 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model6.1 <- MCMCglmm(leaf ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.leaf, prior = prior6.1, #Gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model6.1$VCV) # good trace plots
summary(PL_model6.1) # good effective sample size
autocorr.diag(PL_model6.1$VCV) # no autocorrelation
heidel.diag(PL_model6.1$VCV) # convergence fail

#Trait 6 without maternal effects
load(file="Routput/Plasticity/Plasticity_6.2.RData")

plastic.leaf <- total.2019 %>% filter(germ==1)
prior6.2 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model6.2 <- MCMCglmm(leaf ~ treatment + plot, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.leaf, prior = prior6.2, #Gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model6.2$VCV) # good trace plots
summary(PL_model6.2) # good effective sample size
autocorr.diag(PL_model6.2$VCV) # no autocorrelation
heidel.diag(PL_model6.2$VCV) # convergence fail (but trace plots sufficient)

gen.corrPL6.2 <-PL_model6.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model6.2$VCV[,'treatmentA:treatmentA.animal']*PL_model6.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL6.2) #Post Mean = 0.008148733
posterior.mode(gen.corrPL6.2) #Post Mode = 0.01661873
HPDinterval(gen.corrPL6.2) # Posterior 95% CI = (-0.04875116, 0.07818992)












## Trait 7: Height (Gaussian) ####
load(file="Routput/Plasticity/Plasticity_7.1.RData")

plastic.height <- total.2019 %>% filter(germ==1, !height==0)
prior7.1 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model7.1 <- MCMCglmm(height ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.height, prior = prior7.1, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model7.1$VCV) # good trace plots.. except some bad
summary(PL_model7.1) # good effective sample size
autocorr.diag(PL_model7.1$VCV) # no autocorrelation
heidel.diag(PL_model7.1$VCV) # convergence fail


#Trait 7 without maternal effects
load(file="Routput/Plasticity/Plasticity_7.2.RData")

plastic.height <- total.2019 %>% filter(germ==1, !height==0)
prior7.2 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model7.2 <- MCMCglmm(height ~ treatment + plot, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.height, prior = prior7.2, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model7.2$VCV) # good trace plots
summary(PL_model7.2) # good effective sample size
autocorr.diag(PL_model7.2$VCV) # no autocorrelation
heidel.diag(PL_model7.2$VCV) # convergence fail (but trace plots sufficient)

gen.corrPL7.2 <-PL_model7.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model7.2$VCV[,'treatmentA:treatmentA.animal']*PL_model7.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL7.2) #Post Mean = 0.01237662
posterior.mode(gen.corrPL7.2) #Post Mode = 0.01734095
HPDinterval(gen.corrPL7.2) # Posterior 95% CI = (-0.06465876, 0.09503565)











## Trait 8: Flowering Clusters Number (Gaussian) ####
load(file="Routput/Plasticity/Plasticity_8.1.RData")

plastic.flwrclstr <- total.2019 %>% filter(germ==1)
prior8.1 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model8.1 <- MCMCglmm(flwr_clstr ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.flwrclstr, prior = prior8.1, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model8.1$VCV) # good trace plots
summary(PL_model8.1) # good effective sample size
autocorr.diag(PL_model8.1$VCV) # no autocorrelation
heidel.diag(PL_model8.1$VCV) # convergence fail

## Trait 8 without maternal effects
load(file="Routput/Plasticity/Plasticity_8.2.RData")

plastic.flwrclstr <- total.2019 %>% filter(germ==1)
prior8.2 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model8.2 <- MCMCglmm(flwr_clstr ~ treatment + plot, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.flwrclstr, prior = prior8.2, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model8.2$VCV) # trace plots don't look good
summary(PL_model8.2) # good effective sample size
autocorr.diag(PL_model8.2$VCV) # no autocorrelation
heidel.diag(PL_model8.2$VCV) # convergence fail

gen.corrPL8.2 <-PL_model8.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model8.2$VCV[,'treatmentA:treatmentA.animal']*PL_model8.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL8.2) #Post Mean = 0.003353434
posterior.mode(gen.corrPL8.2) #Post Mode = 0.008966134
HPDinterval(gen.corrPL8.2) # Posterior 95% CI = (-0.0576239, 0.066754)








## Trait 9: Stem Diameter (Gaussian) ####
load(file="Routput/Plasticity/Plasticity_9.1.RData")

plastic.stemdiam <- total.2019 %>% filter(germ==1, !stem_diam==0)
prior9.1 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model9.1 <- MCMCglmm(flwr_clstr ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.stemdiam, prior = prior9.1, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model9.1$VCV) # trace plots don't look great, although just concentrated
summary(PL_model9.1) # good effective sample size
autocorr.diag(PL_model9.1$VCV) # no autocorrelation
heidel.diag(PL_model9.1$VCV) # convergence fail

## Trait 9 without maternal effects
load(file="Routput/Plasticity/Plasticity_9.2.RData")

plastic.stemdiam <- total.2019 %>% filter(germ==1, !stem_diam==0)
prior9.2 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model9.2 <- MCMCglmm(flwr_clstr ~ treatment + plot, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.stemdiam, prior = prior9.2, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model9.2$VCV) # poor trace plots
summary(PL_model9.2) # good effective sample size
autocorr.diag(PL_model9.2$VCV) # no autocorrelation
heidel.diag(PL_model9.2$VCV) # convergence fail 

gen.corrPL9.2 <-PL_model9.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model9.2$VCV[,'treatmentA:treatmentA.animal']*PL_model9.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL9.2) #Post Mean = 0.001621548
posterior.mode(gen.corrPL9.2) #Post Mode = -0.01199429
HPDinterval(gen.corrPL9.2) # Posterior 95% CI = (-0.05717875, 0.0632061)















## Archive : brms code ####

#Assuming Ainv is created, we need the non-Inverse of the Additive Related Matrix for brms input
A <- solve(Ainv$Ainv)
rownames(A) <- rownames(Ainv$Ainv)

PL_model1.2 <- brm(flower ~ treatment + (1|animal) + (1|matID) + (1|treatment:plot), 
                   cov_ranef = list(animal = A), data = plastic.survival,
                   family = "binomial", cores = 4, chains = 4, iter = 2100000, warmup=100000, thin=1000,
                   #Can set seed = to some random number to make model reproducible 
)