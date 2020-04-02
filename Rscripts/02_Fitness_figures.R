#### PROJECT: Brassica rapa GxE Study (Data collected by Cameron So 2019 at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Produce fitness-related figures and models for thesis and publication
#### AUTHOR: Cameron So
#### DATE LAST MODIFIED: 2020/02/02

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Neces+sary Packages
library(gplots)
library(tidyverse)

#Importing data frame 'dat' from the cleaning process of 01_Data_exploration
load(file="Routput/heatarrays_data_explore.RData")

#Standardizing traits of interest for selection curves (sd=1) by treatment
#NOTE: THIS IS ALREADY INCORPORATED IN THE LOADED .RDATA FILE

dat <- dat %>%
  mutate(leaf.st = ifelse(treatment == "A", scale(leaf), scale(leaf)),
         flwr_clstr.st = ifelse(treatment == "A", scale(flwr_clstr), scale(flwr_clstr)),
         seed_pods.st = ifelse(treatment == "A", scale(seed_pods), scale(seed_pods)),
         height.st = ifelse(treatment == "A", scale(height), scale(height)),
         stem_diam.st = ifelse(treatment == "A", scale(stem_diam), scale(stem_diam)),
         germ_census.st = ifelse(treatment == "A", scale(germ_census), scale(germ_census)),
         flwr_census.st = ifelse(treatment == "A", scale(flwr_census), scale(flwr_census)))

#'####################################################################'#
#################      FITNESS-RELATED FIGURES      ##################
#'####################################################################'#

#Graphs 14 - Phenology on Germination Time vs Fitness (seed_pods) ####
dat$germ_census <- as.factor(dat$germ_census)
dat %>%
  filter(germ==1, !germ_census==0) %>%
  group_by(treatment, germ_census) %>%
  summarise(pod_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = germ_census, y = pod_avg, group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=pod_avg-se, ymax=pod_avg+se), width=0.5, size=0.5) + 
  geom_point(aes(color=treatment, size=n)) + 
  ggtitle("Germination Date vs Avg Seed Pod Number") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Census Week", y="Average Seed Pod Number")
#Germinating earlier gives a head start - it allows you to gather more resources to produce more pods


#Graphs 15 - Phenology on Flowering Time vs Fitness (seed_pods) ####
#Note: dates are converted in census # as data was collected weekly
dat$flwr_census <- as.factor(dat$flwr_census)
dat %>%
  filter(germ==1, !flwr_census==0) %>%
  group_by(treatment, flwr_census) %>%
  summarise(pod_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = flwr_census, y = pod_avg, group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=pod_avg-se, ymax=pod_avg+se), width=0.5, size=0.5) + 
  geom_point(aes(color=treatment, size=n)) + 
  ggtitle("First Flowering Date vs Avg Seed Pod Number") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Census Week", y="Average Seed Pod Number")

#Graphs 16 - Fitness curve using relative fitness on LEAF number ####
#To visualize whether there's directional (linear) or stabilizing (quadratic) selection acting, we use family means.

  #Method = linear model
dat %>%
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(leaf.st.avg=mean(leaf.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = leaf.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ggtitle("Selection Curve on Leaf Number") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Leaf Number (standardized)", y="Relative Fitness")

dat %>%
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(leaf.st.avg=mean(leaf.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = leaf.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() + 
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ggtitle("Selection Curve on Leaf Number") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Leaf Number (standardized)", y="Relative Fitness")

#Graphs 17 - Fitness curve using relative fitness on FLOWER CLUSTER number ####
#To visualize whether there's directional (linear) or stabilizing (quadratic) selection acting, we use family means.

#Method = linear model
dat %>%
  filter(germ==1 & flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(flwr_clstr.st.avg=mean(flwr_clstr.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = flwr_clstr.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ggtitle("Selection Curve on Flower Cluster Number") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Flower Cluster Number (standardized)", y="Relative Fitness")

dat %>%
  filter(germ==1 & flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(flwr_clstr.st.avg=mean(flwr_clstr.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = flwr_clstr.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() + 
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ggtitle("Selection Curve on Flower Cluster Number") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Flower Cluster Number (standardized)", y="Relative Fitness")

#Graphs 18 - Fitness curve using relative fitness on HEIGHT ####
#To visualize whether there's directional (linear) or stabilizing (quadratic) selection acting, we use family means.

#Method = linear model
dat %>%
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(height.st.avg=mean(height.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = height.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ggtitle("Selection Curve on Height") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Height (standardized)", y="Relative Fitness")

dat %>%
  filter(germ==1 & flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(height.st.avg=mean(height.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = height.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() + 
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ggtitle("Selection Curve on Height") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Height (standardized)", y="Relative Fitness")


#Graphs 19 - Fitness curve using relative fitness on GERMINATION TIME ####
#To visualize whether there's directional (linear) or stabilizing (quadratic) selection acting, we use family means.

#Method = linear model
dat %>%
  group_by(treatment, famID) %>%
  summarise(germ_census.st.avg=mean(germ_census.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = germ_census.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ggtitle("Selection Curve on Germination Time (by Census)") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Germination Time by Census (standardized)", y="Relative Fitness")

dat %>%
  group_by(treatment, famID) %>%
  summarise(germ_census.st.avg=mean(germ_census.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = germ_census.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ggtitle("Selection Curve on Germination Time (by Census)") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Germination Time by Census (standardized)", y="Relative Fitness")

#Graphs 20 - Fitness curve using relative fitness on FLOWERING TIME ####
#To visualize whether there's directional (linear) or stabilizing (quadratic) selection acting, we use family means.

#Method = linear model
dat %>%
  filter(germ==1 & flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(flwr_census.st.avg=mean(flwr_census.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = flwr_census.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ggtitle("Selection Curve on Flowering Time (by Census)") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Flowering Time by Census (standardized)", y="Relative Fitness")

dat %>%
  filter(germ==1 & flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(flwr_census.st.avg=mean(flwr_census.st), relative.avg=mean(relative), n=n()) %>%
  ggplot(aes(x = flwr_census.st.avg, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ggtitle("Selection Curve on Flowering Time (by Census)") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Flowering Time by Census (standardized)", y="Relative Fitness")


#'####################################################################'#
 ###############      SELECTION MODELS ON TRAITS      ################
#'####################################################################'#

#We employ regression models suggested in Stinchcombe et al. (2008) - 
#https://onlinelibrary-wiley-com.myaccess.library.utoronto.ca/doi/10.1111/j.1558-5646.2008.00449.x

#To test the effects on selection (relative fitness), we run univariate and multivariate analyses
#Univariate: (1) leaf number (2) height at flowering (3) flowering time (4) flowering cluster 
#Multivariate: All of the above as co-variates 
#We use family means to evaluate family-level selection
#We exclude non-germinated plants as we are interested in selection that is applied following germination, as there was treatment applied prior to germination

#All traits are standardized with mean 0 and sd = 0, separate by treatment
#Note, we employ MCMCglmm to standardize the methodology between analysis of genetic variance and selection

library(MCMCglmm)
library(lme4)

#### Family-Level Multivariate Model w/ Traits 1-5 ####

#Subsetting Datasets for both Treatments
ambient.Sel.F.Multi <- dat %>% 
  filter (treatment == "A", germ == 1) %>%
  group_by(famID, plot) %>%
  summarise(relative=mean(relative), leaf.st=mean(leaf.st), height.st=mean(height.st), 
            flwr_census.st=mean(flwr_census.st), flwr_clstr.st=mean(flwr_clstr.st),
            germ_census.st=mean(germ_census.st))
ambient.Sel.F.Multi <- as.data.frame(ambient.Sel.F.Multi)
lapply(ambient.Sel.F.Multi, class)

heated.Sel.F.Multi <- dat %>% 
  filter (treatment == "H", germ == 1) %>%
  group_by(famID, plot) %>%
  summarise(relative=mean(relative), leaf.st=mean(leaf.st), height.st=mean(height.st), 
            flwr_census.st=mean(flwr_census.st), flwr_clstr.st=mean(flwr_clstr.st),
            germ_census.st=mean(germ_census.st))
heated.Sel.F.Multi <- as.data.frame(heated.Sel.F.Multi)

#Verifying Distributions per Covariate Trait
par(mfrow=c(2,3))
hist(log(ambient.Sel.F.Multi$relative), breaks=30) #Gaussian if log transformed
hist(ambient.Sel.F.Multi$leaf.st, breaks=30) #Gaussian
hist(ambient.Sel.F.Multi$height.st, breaks=30) #Gaussian
hist(ambient.Sel.F.Multi$germ_census.st, breaks=30) #Gaussian .. or poisson
hist(ambient.Sel.F.Multi$flwr_census.st, breaks=30) #Gaussian
hist(ambient.Sel.F.Multi$flwr_clstr.st, breaks=30) #Poisson

hist(log(heated.Sel.F.Multi$relative), breaks=30) #Gaussian if log transformed
hist(heated.Sel.F.Multi$leaf.st, breaks=30) #Gaussian
hist(heated.Sel.F.Multi$height.st, breaks=30) #Gaussian
hist(heated.Sel.F.Multi$germ_census.st, breaks=30) #Gaussian .. or poisson
hist(heated.Sel.F.Multi$flwr_census.st, breaks=30) #Gaussian
hist(heated.Sel.F.Multi$flwr_clstr.st, breaks=30) #Poisson

#Prior for both treatments
prior1.1Sel <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 0.002)))

#Ambient Multivariate Model
AM_Sel.F.Model1 <- MCMCglmm(log(relative) ~ leaf.st + height.st + germ_census.st + flwr_census.st + flwr_clstr.st,
                          random = ~plot, family = "gaussian",
                          data = ambient.Sel.F.Multi, prior = prior1.1Sel,
                          nitt=100000, thin=100, burnin=10000, verbose = T, pr = TRUE) 

plot(AM_Sel.F.Model1) # good trace plots except 'Plot'
summary(AM_Sel.F.Model1) #good effective sample sizes
autocorr.diag(AM_Sel.F.Model1$Sol) #no autocorrelation
heidel.diag(AM_Sel.F.Model1$Sol) #Convergence for all random effects

#Heated Multivariate Model
HT_Sel.F.Model1 <- MCMCglmm(log(relative) ~ leaf.st + height.st + germ_census.st + flwr_census.st + flwr_clstr.st,
                          random = ~plot, family = "gaussian",
                          data = heated.Sel.F.Multi, prior = prior1.1Sel,
                          nitt=100000, thin=100, burnin=10000, verbose = T, pr = TRUE) 

plot(HT_Sel.F.Model1) # good trace plots except 'Plot'
summary(HT_Sel.F.Model1) #good effective sample sizes
autocorr.diag(HT_Sel.F.Model1$Sol) #no autocorrelation when thin = 1000
heidel.diag(HT_Sel.F.Model1$Sol) #Convergence for all random effects



#### Individual-Level Multivariate Model w/ Traits 1-5 ####

#Subsetting Datasets for both Treatments
ambient.Sel.I.Multi <- dat %>% 
  filter (treatment == "A", germ == 1)
ambient.Sel.I.Multi <- as.data.frame(ambient.Sel.I.Multi)

heated.Sel.I.Multi <- dat %>% 
  filter (treatment == "H", germ == 1)
heated.Sel.I.Multi <- as.data.frame(heated.Sel.I.Multi)

#Verifying Distributions per Covariate Trait
hist(log(ambient.Sel.I.Multi$relative), breaks=30) #Gaussian (because non-integers)
hist(ambient.Sel.I.Multi$leaf.st, breaks=30) #Gaussian
hist(ambient.Sel.I.Multi$height.st, breaks=30) #Gaussian
hist(ambient.Sel.I.Multi$germ_census.st, breaks=30) #Gaussian .. or poisson
hist(ambient.Sel.I.Multi$flwr_census.st, breaks=30) #Gaussian
hist(ambient.Sel.I.Multi$flwr_clstr.st, breaks=30) #Poisson

hist(log(heated.Sel.I.Multi$relative), breaks=30) #Gaussian (because non-integers)
hist(heated.Sel.I.Multi$leaf.st, breaks=30) #Gaussian
hist(heated.Sel.I.Multi$height.st, breaks=30) #Gaussian
hist(heated.Sel.I.Multi$germ_census.st, breaks=30) #Gaussian .. or poisson
hist(heated.Sel.I.Multi$flwr_census.st, breaks=30) #Gaussian
hist(heated.Sel.I.Multi$flwr_clstr.st, breaks=30) #Poisson

#Prior for both treatments ... reused from above
prior2.1Sel <- list(R = list(V = 1, nu = 0.002),
                  G = list(G1 = list(V = 1, nu = 0.002),
                           G2 = list(V = 1, nu = 0.002)))

#Ambient Multivariate Model
AM_Sel.I.Model2 <- MCMCglmm(log(relative) ~ leaf.st + height.st + germ_census.st + flwr_census.st + flwr_clstr.st,
                            random = ~plot + famID, family = "gaussian",
                            data = ambient.Sel.I.Multi, prior = prior2.1Sel,
                            nitt=100000, thin=100, burnin=10000, verbose = T, pr = TRUE) 

plot(AM_Sel.I.Model1) # good trace plots except 'Plot'
summary(AM_Sel.I.Model1) #good effective sample sizes
autocorr.diag(AM_Sel.I.Model1$Sol) #no autocorrelation
heidel.diag(AM_Sel.I.Model1$Sol) #Convergence for all random effects
gelman.diag()

#Heated Multivariate Model
HT_Sel.I.Model2 <- MCMCglmm(log(relative) ~ leaf.st + height.st + germ_census.st + flwr_census.st + flwr_clstr.st,
                            random = ~plot + famID, family = "gaussian",
                            data = heated.Sel.I.Multi, prior = prior2.1Sel,
                            nitt=100000, thin=100, burnin=10000, verbose = T, pr = TRUE) 

plot(HT_Sel.I.Model1) # good trace plots except 'Plot'
summary(HT_Sel.I.Model1) #good effective sample sizes
autocorr.diag(HT_Sel.I.Model1$Sol) #no autocorrelation when thin = 1000
heidel.diag(HT_Sel.I.Model1$Sol) #Convergence for all random effects


#'####################################################################'#
#############      CHANGES IN PHENOTYPIC MEAN MODELS      ##############
#'####################################################################'#

library(lme4)
library(tidyverse)
#library(gtools) #for logit, can't use car
library(car)
library(MASS)
library(vcd)
library(performance)
library(see)

## Trait 1: Lifetime Fitness ####
dat %>% group_by(treatment) %>% 
  summarise(mean=mean(seed_pods), st.dev=sd(seed_pods), var=(st.dev)^2) #Variance is much greater than the mean
phen1 <- as.data.frame(dat)

hist(phen1$seed_pods) 
mod1a <- lm(seed_pods ~ treatment, data = phen1)
qqnorm(rstandard(mod1a))
qqline(rstandard(mod1a)) #Residuals not normally distributed
gg_reshist(mod1a, bins = 30) #Residuals definitely not normally distributed
check_model(mod1a)

#Logistic regression with log link (family=poisson) is used, but not great
mod1 <- glmer(data = phen1, seed_pods ~ treatment + (1|treatment:plot) + (1|individual), family = "poisson")
summary(mod1)
exp(0.4591)
check_overdispersion(mod1) #No overdispersion apparently
check_zeroinflation(mod1) #Not zero-inflated??? I don't believe this but OK.
#The overdispersion and zero-inflation is removed when (1|individual) is added.
#Regardless, we can't use this model because the residuals are NOT normally distributed

# Need to run either zero-inflated Poisson model, ZI neg-binomial, or bootstrap
library(glmmTMB) #for Zero-inflated Poisson/Neg-Binomial model
# I use a ZI Neg-Bin model because of high variance relative to the mean
# I follow the tutorial: https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf
mod1b <- glmmTMB(seed_pods ~ treatment + (1|treatment:plot), data = phen1, 
                 ziformula = ~1, family = nbinom2)
summary(mod1b)
Anova(mod1b)
exp(0.6114) #1.84301 change due to treatment



## Trait 2: Survival ####
dat %>% group_by(treatment) %>% 
  summarise(sum=sum(flower), n=n(), prob=sum/n, var=(prob*(1-prob)), sd=sqrt(var))
phen2 <- dat

#Model: Logistic Regression using a logit link via family = "binomial"
mod2 <- glmer(data = phen2, flower ~ treatment + (1|treatment:plot), family = "binomial")
summary(mod2) #No significant effect of treatment
Anova(mod2)
#In Anova, type 2 = default, type 1 = assumed balanced data in levels, type 3 = deals best with imbalance
#type 1 is use in summary(mod2)
#Can change Anova to F -statistic than Wald chisquare





## Trait 3: Fecundity ####
dat %>% group_by(treatment) %>% filter(germ==1, flower==1) %>% 
  summarise(mean=mean(seed_pods), st.dev=sd(seed_pods), var=(st.dev)^2)
phen3 <- dat %>% filter(germ==1, flower==1)

hist(phen3$seed_pods, breaks=30) #Poisson distribution
mod3a <- lm(seed_pods ~ treatment, data = phen3)
qqnorm(rstandard(mod3a))
qqline(rstandard(mod3a)) #Residuals somewhat deviating from normal
gg_reshist(mod3a, bins = 30) #Residuals somewhat normal
check_model(mod3a)

#Model:
mod3 <- glmer(data = phen3, seed_pods ~ treatment + (1|treatment:plot) + (1|individual), family = "poisson")
check_overdispersion(mod3)
check_zeroinflation(mod3) #Probable zero-inflation
summary(mod3)
Anova(mod3)

#Need to run a zero-inflated model
mod3b <- glmmTMB(seed_pods ~ treatment + (1|treatment:plot) + (1|individual), data = phen3, 
                 ziformula = ~1, family = nbinom2)
summary(mod3b) #significant zero-inflation model
Anova(mod3b)
exp(0.6114) #0.4883 change due to treatment





## Trait 4: Germination Success ####
dat %>% group_by(treatment) %>% 
  summarise(sum=sum(germ), n=n(), prob=sum/n, var=(prob*(1-prob)), sd=sqrt(var))
phen4 <- dat

#Logistic Regression using a logit link via family = "binomial"
mod4 <- glmer(data = phen4, germ ~ treatment + (1|treatment:plot), family = "binomial")
summary(mod4) #No significant effect of treatment on germination
Anova(mod4)








## Trait 5: Flowering Success ####
dat %>% group_by(treatment) %>% filter(germ==1) %>% 
  summarise(sum=sum(flower), n=n(), prob=sum/n, var=(prob*(1-prob)), sd=sqrt(var))
phen5 <- dat %>% filter(germ==1)

#Logistic Regression using a logit link via family = "binomial"
mod5 <- glmer(data = phen5, flower ~ treatment + (1|treatment:plot), family = "binomial")
summary(mod5) #No significant effect of treatment on flowering success
Anova(mod5)








## Trait 6: Fruiting Success ####
dat %>% group_by(treatment) %>% filter(germ==1, flower==1) %>% 
  summarise(sum=sum(seed), n=n(), prob=sum/n, var=(prob*(1-prob)), sd=sqrt(var))
phen6 <- dat %>% filter(germ==1, flower==1)

#Logistic Regression using a logit link via family = "binomial"
mod6 <- glmer(data = phen6, seed ~ treatment + (1|treatment:plot), family = "binomial")
summary(mod6) #No significant effect of treatment on fruiting success
Anova(mod6)








## Trait 7: Leaf Number ####
dat %>% group_by(treatment) %>% filter(germ==1) %>% 
  summarise(mean=mean(leaf), st.dev=sd(leaf), var=(st.dev)^2)
phen7 <- dat %>% filter(germ==1)

hist(phen7$leaf, breaks=30) #Gaussian distribution
mod7a <- lm(leaf ~ treatment, data = phen7)
qqnorm(rstandard(mod7a))
qqline(rstandard(mod7a)) #Residuals fairly normal
gg_reshist(mod7a, bins = 30) #Residuals somewhat normal
check_model(mod7a)

#Linear Regression assuming a Gaussian distribution
mod7 <- lmer(data = phen7, leaf ~ treatment + (1|treatment:plot), REML=FALSE) #AIC = ***
check_heteroscedasticity(mod7) #data is heteroscedasticitic ... 
summary(mod7)
Anova(mod7) #for p values .. significant treatment effect



  
## Trait 8: Height ####
dat %>% group_by(treatment) %>% filter(germ==1, !height==0) %>% 
  summarise(mean=mean(height), st.dev=sd(height), var=(st.dev)^2)
phen8 <- dat %>% filter(germ==1, !height==0) #Removing 0s because previously NAs

hist(phen8$height, breaks=30) #Poisson distribution
mod8a <- lm(height ~ treatment, data = phen8)
qqnorm(rstandard(mod8a))
qqline(rstandard(mod8a)) #Residuals fairly normal
gg_reshist(mod8a, bins = 30) #Residuals somewhat normal
check_model(mod8a)

#Linear Regression assuming a Gaussian distribution
mod8 <- lmer(data = phen8, height ~ treatment + (1|treatment:plot), REML=FALSE) #AIC = ***
check_heteroscedasticity(mod8)
summary(mod8)
Anova(mod8) #for p values



## Trait 9: Flowering Clusters Number ####
dat %>% group_by(treatment) %>% filter(germ==1, flower==1) %>% 
  summarise(mean=mean(flwr_clstr), st.dev=sd(flwr_clstr), var=(st.dev)^2)
phen9 <- dat %>% filter(germ==1, flower==1)

hist(phen9$flwr_clstr, breaks=30) #Poisson distribution
mod9a <- lm(height ~ treatment, data = phen8)
qqnorm(rstandard(mod9a))
qqline(rstandard(mod9a)) #Residuals fairly normal
gg_reshist(mod9a, bins = 30) #Residuals somewhat normal
check_model(mod9a)

#GLMM Poisson model using a log link
mod9 <- glmer(data = phen9, flwr_clstr ~ treatment + (1|treatment:plot) + (1|individual), family = "poisson")
check_overdispersion(mod9) #Overdispersion - added + (1|individual)
check_zeroinflation(mod9) 
summary(mod9)
Anova(mod9) #for p values






## Trait 10: Stem Diameter ####
dat %>% group_by(treatment) %>% filter(germ==1, !stem_diam==0) %>% 
  summarise(mean=mean(stem_diam), mode=mode(stem_diam), st.dev=sd(stem_diam), var=(st.dev)^2)
phen10 <- dat %>% filter(germ==1, !stem_diam==0, !stem_diam==0.5) #Removing 0s because previously NAs

hist(phen10$stem_diam, breaks=30) #Poisson distribution
mod10a <- lm(height ~ treatment, data = phen10)
qqnorm(rstandard(mod10a))
qqline(rstandard(mod10a)) #Residuals fairly normal
gg_reshist(mod10a, bins = 30) #Residuals somewhat normal
check_model(mod10a)

mod10 <- lmer(data = phen10, stem_diam ~ treatment + (1|treatment:plot)) #AIC = ***
check_heteroscedasticity(mod10)
Anova(mod10) #for p values
summary(mod10)

