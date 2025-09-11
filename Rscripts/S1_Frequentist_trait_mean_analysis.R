#### PROJECT: Brassica rapa Va/W Study (Data collected at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Produce population+treatment level trait mean figures and models for manuscript

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Necessary Packages
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

#Graphs 21 - Fitness curve using relative fitness on SURVIVAL ####
#To visualize whether there's directional (linear) or stabilizing (quadratic) selection acting, we use family means.

bivariate <- dat %>%
  group_by(treatment, famID) %>%
  summarise(sum_flower=sum(flower), n=n(), mean_pods=mean(seed_pods)) %>%
  mutate(mean_flower=(sum_flower/n))

#Method = linear model
dat %>%
  group_by(treatment, famID) %>%
  summarise(sum_flower=sum(flower), n=n(), relative.avg=mean(relative)) %>%
  mutate(mean_flower=(sum_flower/n),
         mean_flower.st=ifelse(treatment == "A", scale(mean_flower), scale(mean_flower))) %>%
  ggplot(aes(x = mean_flower.st, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  ggtitle("Selection Curve on Flowering Time (by Census)") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Survival Success (standardized)", y="Relative Fitness")

dat %>%
  group_by(treatment, famID) %>%
  summarise(sum_flower=sum(flower), n=n(), relative.avg=mean(relative)) %>%
  mutate(mean_flower=(sum_flower/n),
         mean_flower.st=ifelse(treatment == "A", scale(mean_flower), scale(mean_flower))) %>%
  ggplot(aes(x = mean_flower.st, y = relative.avg, group = treatment, color = treatment)) + 
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  ggtitle("Selection Curve on Flowering Time (by Census)") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Survival Success by Census (standardized)", y="Relative Fitness")

#'##################################################################################'#
#############      CHANGES IN PHENOTYPIC MEAN MODELS + Vg ESTIMATES     ##############
#'##################################################################################'#

library(lme4)
library(tidyverse)
#library(gtools) #for logit, can't use car
library(car)
library(MASS)
library(vcd)
library(performance)
library(see)
library(merTools)
library(lmerTest)
library(DHARMa)
library(gplots)
library(lindia)
library(patchwork)

#Function to get 95% CI from lme4 models when ignoring random effects
easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link function
  linkinv <- family(model)$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low=pred0-crit*pred.se,
                conf.high=pred0+crit*pred.se))
}

#NOTE: We do not nest plot within treatment BECAUSE plot is numerically ordered from 1-12, where 6 are under ambient and 6 are under heated. 


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
mod1 <- glmer(data = phen1, seed_pods ~ treatment + (1|famID) + (1|plot), family = "poisson")
summary(mod1)
check_overdispersion(mod1) #Overdispersion apparently
check_zeroinflation(mod1) #Probable zero-inflation

# Need to run either zero-inflated Poisson model, ZI neg-binomial, or bootstrap
library(glmmTMB) #for Zero-inflated Poisson/Neg-Binomial model
# I use a ZI Neg-Bin model because of high variance relative to the mean
# I follow the tutorial: https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf
mod1b <- glmmTMB(seed_pods ~ treatment  + (1|famID) + (1|plot), data = phen1, 
                 ziformula = ~treatment + (1|famID) + (1|plot), family = nbinom2)
summary(mod1b)

  #Model Diagnostics
testDispersion(mod1b) #no dispersion
sim_mod1b <- simulateResiduals(fittedModel = mod1b) #calculating residuals using DHARMa package
plot(sim_mod1b, asFactor = T) #Change asFactor = T to scale across random effects
#No dispersion, outlier not concerning, within-group deviations is not too concerning either. 
plotResiduals(sim_mod1b, phen1$treatment, quantreg = T) #residuals plotted against treatment
testDispersion(sim_mod1b)

  #Test of treatment significance
Anova(mod1b, type = 2)

  #Model comparison with a GxE
mod1c <- glmmTMB(seed_pods ~ treatment  + (1|famID) + (1|plot) + (1|famID:treatment), data = phen1, 
                 ziformula = ~treatment + (1|famID) + (1|plot) + (1|famID:treatment), family = nbinom2)
anova(mod1b, mod1c) #mod1b slightly better

  #Getting estimated mean trait values per treatment
newdata1 <- data.frame(treatment = factor(c("A","A","A","A","A","A", "H","H","H","H","H","H")),
                       plot = factor(c("2", "3", "5", "8", "11", "12", "1", "4", "6", "7", "9", "10"))
                       ) #if want to include random effects
newdata2 <- data.frame(treatment = factor(levels(phen1$treatment), levels = levels(phen1$treatment)), plot = NA) #if want to exclude R.E
pred1a <- predict(mod1b, newdata2, type = "response", re.form = NA) #NULL to include all random effects, NA to exclude
pred1a
#Ambient estimate = 4.628692 including 0 inflation parameter
#Heated estimate = 8.384884 including 0 inflation parameter
exp(2.3662) #Alternatively, Ambient latent scale trait mean = 10.65682 (data collected from summary(mod1b))
exp(2.3662 + 0.5679) #Heated latent = 18.80457 when excluding 0 inflation parameter

  #Standard error estimates
pred1b <- predict(mod1b, newdata2, type = "response", se.fit = TRUE, re.form = NA)  #standard errors
pred1b
fit1 <- pred1b$fit
lwr1 <- fit1 - 1.96*pred1b$se.fit #95% confidence intervals (upper) .. first ambient, then heated
upr1 <- fit1 + 1.96*pred1b$se.fit #95% confidence intervals (lower) .. first ambient, then heated
lwr1
upr1






## Trait 2: Survival ####
dat %>% group_by(treatment) %>% 
  summarise(sum=sum(flower), n=n(), prob=sum/n, var=(prob*(1-prob)), sd=sqrt(var))
phen2 <- dat

#Model: Logistic Regression using a logit link via family = "binomial"
mod2a <- glmer(data = phen2, flower ~ treatment + (1|plot), family = "binomial")
summary(mod2a) #No significant effect of treatment
summary(residuals(mod2a))
plot(mod2a)
  
  #Test of treatment significance
Anova(mod2a)
#In Anova, type 2 = default, type 1 = assumed balanced data in levels, type 3 = deals best with imbalance
#type 1 is use in summary(mod2)
#Can change Anova to F -statistic than Wald chisquare

#Model comparison with a GxE
mod2b <- glmer(data = phen2, flower ~ treatment + (1|famID:treatment) + (1|plot), family = "binomial")
mod2c <- glm(data = phen2, flower ~ treatment, family = "binomial")
anova(mod2a, mod2b, mod2c) #mod2b slightly better, but we stick with mod2a as there's too many random effects to estimate with the data and downstream models run into problems (with estimating standard errors)

  #Model Diagnostics
testDispersion(mod2b)
sim_mod2 <- simulateResiduals(fittedModel = mod2) #calculating residuals using DHARMa package
plot(sim_mod2, asFactor = T) #Change asFactor = T to scale across random effects
#no issues with residuals, dispersion, etc
plotResiduals(sim_mod2, phen2$treatment, quantreg = T) #residuals plotted against treatment
testDispersion(sim_mod2)

  #Getting estimated values per treatment
newdata2 <- data.frame(treatment = factor(levels(phen1$treatment), levels = levels(phen1$treatment)), plot = NA)
predict(mod2b, newdata2, type = "response", re.form = NA) 
#Ambient estimate = 0.4077725
#Heated estimate = 0.4421335

  #Estimated (latent scale) values using plogis using manual calculation
summary(mod2b)
plogis(-0.3732) # Ambient probability estimate = 0.407768
plogis(-0.3732 + 0.1407) #Heated probability estimate = 0.4421354


  #Standard error estimates
newdata2 <- data.frame(treatment = factor(1:2, levels = 1:2, labels = levels(phen2$treatment)), plot=NA)
pred2 <- bootMer(mod2a, FUN=function(x) predict(x, re.form = NA, newdata = newdata2, type = "response"),
                 nsim = 1000) #standard errors
pred2 
easyPredCI(mod2a, newdata2) #95% confidence intervals

##  inverse link function (e.g. plogis() for binomial, beta;
##  exp() for Poisson, negative binomial



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

#Model 3 - glmm with Poisson
mod3a <- glmer(data = phen3, seed_pods ~ treatment + (1|plot), family = "poisson")
summary(mod3a)
Anova(mod3a)

  #Model Diagnostics
testDispersion(mod3a) # overdispersion
testZeroInflation(mod3a) #some zero inflation
sim_mod3a <- simulateResiduals(fittedModel = mod3a) #calculating residuals using DHARMa package
plot(sim_mod3a, asFactor = T) #Change asFactor = T to scale across random effects

#Model 3b - glmm with (1|individual) because of overdispersion
mod3b <- glmer(data = phen3, seed_pods ~ treatment + (1|plot) + (1|individual), family = "poisson")

  #Model Diagnostics
testDispersion(mod3b)
testZeroInflation(mod3b)
sim_mod3b <- simulateResiduals(fittedModel = mod3b) #calculating residuals using DHARMa package
plot(sim_mod3b, asFactor = T) #Change asFactor = T to scale across random effects

#Model 3c zero-inflated model
mod3c <- glmmTMB(seed_pods ~ treatment + (1|plot), data = phen3, 
                 ziformula = ~1, family = nbinom2)
summary(mod3c) #zero inflation contributes very little to the variance


#Model comparison with a GxE
mod3d <- glmmTMB(data = phen3, seed_pods ~ treatment + (1|famID:treatment) + (1|plot), ziformula = ~1, family = "poisson")
anova(mod3c, mod3d) # mod3c better

  #Model Diagnostics
testDispersion(mod3c) #no dispersion
testZeroInflation(mod3c) #some zero inflation
sim_mod3c <- simulateResiduals(fittedModel = mod3c) #calculating residuals using DHARMa package
plot(sim_mod3c, asFactor = T) #Change asFactor = T to scale across random effects
testOutliers(sim_mod3c, type="bootstrap", nBoot = nSim) #No significant outliers
#homogeneity not too concerning. Within margin of acceptance
plotResiduals(sim_mod3c, phen3$treatment, quantreg = T) #residuals plotted against treatment
testDispersion(sim_mod3c)

#CONCLUSION: we use model 3c because we deal with zero inflation appropriately in the data
#Test of significance
summary(mod3c)
Anova(mod3c)

  #Getting estimated values per treatment
newdata3 <- data.frame(treatment = factor(1:2, levels = 1:2, labels = levels(phen3$treatment)), plot=NA)
pred3 <- predict(mod3c, newdata3, type = "response", re.form = NA) #NULL to include all random effects
pred3
#Ambient estimate = 12.03740
#Heated estimate = 20.82984
exp(2.4880) #Alternatively, Ambient odds ratio = 11.06
exp(2.4880 + 0.5484) #Heated odds ratio = 20.30973

#Standard error estimates
pred3b <- predict(mod3c, newdata3, type = "response", se.fit = TRUE, re.form = NA)  #standard errors
pred3b
fit3 <- pred3c$fit
lwr3 <- fit3 - 1.96*pred3b$se.fit #95% confidence intervals (upper) .. first ambient, then heated
upr3 <- fit3 + 1.96*pred3b$se.fit #95% confidence intervals (lower) .. first ambient, then heated
lwr3
upr3









## Trait 4: Germination Success ####
dat %>% group_by(treatment) %>% 
  summarise(sum=sum(germ), n=n(), prob=sum/n, var=(prob*(1-prob)), sd=sqrt(var))
phen4 <- dat

#Logistic Regression using a logit link via family = "binomial"
mod4a <- glmer(data = phen4, germ ~ treatment + (1|plot), family = "binomial")
summary(mod4a) #No significant effect of treatment on germination

  #Test of significance
Anova(mod4a)

#Model comparison with a GxE
mod4b <- glmer(data = phen4, germ ~ treatment + (1|famID:treatment) + (1|plot), family = "binomial")
anova(mod4a, mod4b) # ______ slightly better

  #Model Diagnostics
testDispersion(mod4a)
sim_mod4 <- simulateResiduals(fittedModel = mod4a) #calculating residuals using DHARMa package
plot(sim_mod4a, asFactor = T) #Change asFactor = T to scale across random effects
plotResiduals(sim_mod4a, phen4$treatment, quantreg = T) #residuals plotted against treatment
testDispersion(sim_mod4a)

  #Getting estimated values per treatment
newdata4 <- data.frame(treatment = factor(levels(phen4$treatment), levels = levels(phen4$treatment)), plot = NA)
predict(mod4a, newdata4, type = "response", re.form = NA) 
#Ambient estimate = 0.4591922
#Heated estimate = 0.4723672

  #Estimated values using plogis
plogis(-0.16360) # Ambient probability estimate = 0.4104997
plogis(-0.16360 + 0.05295) #Heated probability estimate = 0.4435418

  #Standard error estimates
newdata4 <- data.frame(treatment = factor(levels(phen4$treatment), levels = levels(phen4$treatment)), plot = NA)
pred4 <- bootMer(mod4a, FUN=function(x) predict(x, re.form = NA, newdata = newdata4, type = "response"),
                 nsim = 1000) #standard errors with bootstrap 1000 simulations
pred4
easyPredCI(mod4a, newdata4) #95% confidence intervals

##  inverse link function (e.g. plogis() for binomial, beta;
##  exp() for Poisson, negative binomial






## Trait 5: Flowering Success ####
dat %>% group_by(treatment) %>% filter(germ==1) %>% 
  summarise(sum=sum(flower), n=n(), prob=sum/n, var=(prob*(1-prob)), sd=sqrt(var))
phen5 <- dat %>% filter(germ==1)

#Logistic Regression using a logit link via family = "binomial"
mod5a <- glmer(data = phen5, flower ~ treatment + (1|plot), family = "binomial")
summary(mod5a) #No significant effect of treatment on flowering success

  #Test of significance
Anova(mod5a)

#Model comparison with a GxE
mod5b <- glmer(data = phen5, flower ~ treatment + (1|famID:treatment) + (1|plot), family = "binomial")
anova(mod5a, mod5b) # ______ slightly better

  #Model Diagnostics
testDispersion(mod5) #no dispersion
sim_mod5 <- simulateResiduals(fittedModel = mod5) #calculating residuals using DHARMa package
plot(sim_mod5, asFactor = T) #Change asFactor = T to scale across random effects
plotResiduals(sim_mod5, phen5$treatment, quantreg = T) #residuals plotted against treatment
testDispersion(sim_mod5)

  #Getting estimated values per treatment
newdata5 <- data.frame(treatment = factor(levels(phen5$treatment), levels = levels(phen5$treatment)), plot = NA)
predict(mod5, newdata5, type = "response", re.form = NA) 
#Ambient estimate = 0.8970709
#Heated estimate = 0.9425588

  #Estimated values using plogis
plogis(2.1651) # Ambient probability estimate = 0.8970714
plogis(2.1651 + 0.6327) #Heated probability estimate = 0.9425568

exp(2.1651) #Alternatively, Ambient odds ratio = 8.715473
exp(2.1651 + 0.6327) #Heated odds ratio = 16.40851

  #Standard error estimates
newdata5 <- data.frame(treatment = factor(levels(phen5$treatment), levels = levels(phen5$treatment)), plot = NA)
pred5 <- bootMer(mod5, FUN=function(x) predict(x, re.form = NA, newdata = newdata5, type = "response"),
                 nsim = 1000) #standard errors with bootstrap 1000 simulations
pred5
easyPredCI(mod5, newdata5) #95% confidence intervals

##  inverse link function (e.g. plogis() for binomial, beta;
##  exp() for Poisson, negative binomial






## Trait 6: Fruiting Success ####
dat %>% group_by(treatment) %>% filter(germ==1, flower==1) %>% 
  summarise(sum=sum(seed), n=n(), prob=sum/n, var=(prob*(1-prob)), sd=sqrt(var))
phen6 <- dat %>% filter(germ==1, flower==1)

#Logistic Regression using a logit link via family = "binomial"
mod6a <- glmer(data = phen6, seed ~ treatment + (1|plot), family = "binomial")
summary(mod6a) #No significant effect of treatment on fruiting success
  
  #Test of significance
Anova(mod6a)

#Model comparison with a GxE
mod6b <- glmer(data = phen6, seed ~ treatment + (1|famID:treatment) + (1|plot), family = "binomial")
anova(mod6a, mod6b) # ______ slightly better

  #Model Diagnostics
testDispersion(mod6a)
sim_mod6 <- simulateResiduals(fittedModel = mod6a) #calculating residuals using DHARMa package
plot(sim_mod6, asFactor = T) #Change asFactor = T to scale across random effects
plotResiduals(sim_mod6, phen6$treatment, quantreg = T) #residuals plotted against treatment
testDispersion(sim_mod6)

  #Getting estimated values per treatment
newdata6 <- data.frame(treatment = factor(levels(phen6$treatment), levels = levels(phen6$treatment)), plot = NA)
predict(mod6a, newdata6, type = "response", re.form = NA) 
#Ambient estimate = 0.9929627
#Heated estimate = 0.9959487

  #Estimated values using plogis
plogis(4.9495) # Ambient probability estimate = 0.9929629
plogis(4.9495 + 0.5552) #Heated probability estimate = 0.9959489

exp(4.9495) #Alternatively, Ambient odds ratio = 141.1044
exp(4.9495 + 0.5552) #Heated odds ratio = 245.8447

  #Standard error estimates
newdata6 <- data.frame(treatment = factor(levels(phen6$treatment), levels = levels(phen6$treatment)), plot = NA)
pred6 <- bootMer(mod6a, FUN=function(x) predict(x, re.form = NA, newdata = newdata6, type = "response"),
                 nsim = 1000) #standard errors with bootstrap 1000 simulations
pred6
easyPredCI(mod6a, newdata6) #95% confidence intervals

##  inverse link function (e.g. plogis() for binomial, beta;
##  exp() for Poisson, negative binomial






## Trait 7: Leaf Number ####
dat %>% group_by(treatment) %>% filter(germ==1) %>% 
  summarise(mean=mean(leaf), st.dev=sd(leaf), var=(st.dev)^2)
phen7 <- dat %>% filter(germ==1)

hist(phen7$leaf, breaks=50) #Gaussian distribution
mod7a <- lm(leaf ~ treatment, data = phen7)
qqnorm(rstandard(mod7a))
qqline(rstandard(mod7a)) #Residuals fairly normal
gg_reshist(mod7a, bins = 30) #Residuals somewhat normal
check_model(mod7a)

#Linear Regression assuming a Gaussian distribution
mod7a <- lmer(data = phen7, leaf ~ treatment + (1|plot), REML=FALSE) #AIC = ***
check_heteroscedasticity(mod7a) #data is heteroscedasticitic ... 
summary(mod7a)
Anova(mod7a) #for p values .. significant treatment effect

#Model: Poisson model to see if better fit
mod7b <- glmer(data = phen7, leaf ~ treatment + (1|plot), family = "poisson")
check_overdispersion(mod7b) #no overdispersion
check_zeroinflation(mod7b) #Probable zero-inflation
summary(mod7b)
anova(mod7b, mod7a) #model with gaussian error structure fits better

#We use mod7a even though there's some heteroscedasticity.. residuals are not normally distributed but this is a minor issue as the histograms look OK. 
#There is some precision to be lost for derived p values. (can get lower than true values) and the standard errors. 

#Model comparison with a GxE
mod7c <- glmer(data = phen7, leaf ~ treatment + (1|famID:treatment) + (1|plot), family = "poisson")
anova(mod7b, mod7c) # mod7c only slightly better, so going to use 7b because don't want to lose power in detecting signal in famID

#Getting estimated values per treatment
newdata7 <- data.frame(treatment = factor(1:2, levels = 1:2, labels = levels(phen7$treatment)), plot=NA)
pred7 <- predict(mod7a, newdata7, type = "response", re.form = NA) #NULL to include all random effects
pred7
#Ambient estimate = 6.48819
#Heated estimate = 7.807218
6.4882 #Alternatively, extract from summary .. Ambient = 6.4882
6.4882+1.3190 #Heated = 7.8072

#Standard error estimates
summary(mod7a)











  
## Trait 8: Height ####
dat %>% group_by(treatment) %>% filter(germ==1, !height==0) %>% 
  summarise(mean=mean(height), st.dev=sd(height), var=(st.dev)^2)
phen8 <- dat %>% filter(germ==1, !height==0) #Removing 0s because previously NAs

hist(phen8$height, breaks=50) 
mod8a <- lm(height ~ treatment, data = phen8)
qqnorm(rstandard(mod8a))
qqline(rstandard(mod8a)) #Residuals fairly normal
gg_reshist(mod8a, bins = 30) #Residuals somewhat normal
check_model(mod8a)

#Linear Regression assuming a Gaussian distribution
mod8a <- lmer(data = phen8, height ~ treatment + (1|plot), REML=FALSE) #AIC = ***
check_heteroscedasticity(mod8a) #heteroscedastic but fairly normal Q-Q plot... so we choose to ignore this
summary(mod8a)
Anova(mod8a) #for p values

#Model comparison with a GxE
mod8b <- lmer(data = phen8, height ~ treatment + (1|famID:treatment) + (1|plot), REML=FALSE)
anova(mod8a, mod8b) # mod8b slightly better


#Getting estimated values per treatment
newdata8 <- data.frame(treatment = factor(1:2, levels = 1:2, labels = levels(phen8$treatment)), plot=NA)
pred8 <- predict(mod8b, newdata8, type = "response", re.form = NA) #NULL to include all random effects
pred8
#Ambient estimate = 30.61775
#Heated estimate = 37.27157
summary(mod8b)
30.61775 #Alternatively, extract from summary .. Ambient = 30.617
30.61775+6.654 #Heated = 37.27











## Trait 9: Flowering Clusters Number ####
dat %>% group_by(treatment) %>% filter(germ==1, flower==1) %>% 
  summarise(mean=mean(flwr_clstr), st.dev=sd(flwr_clstr), var=(st.dev)^2)
phen9 <- dat %>% filter(germ==1, flower==1)

hist(phen9$flwr_clstr, breaks=50) #Poisson distribution ... lots of 1s
mod9a <- lm(flwr_clstr ~ treatment, data = phen8)
qqnorm(rstandard(mod9a))
qqline(rstandard(mod9a)) #Residuals not normal
check_model(mod9a)

#GLMM Poisson model using a log link
mod9b <- glmer(data = phen9, flwr_clstr ~ treatment + (1|plot), family = "poisson")
check_overdispersion(mod9b) #Overdispersion - add + (1|individual)!

  #Model Diagnostics
testDispersion(mod9b) #still issues!
sim_mod9b <- simulateResiduals(fittedModel = mod9b)
plot(sim_mod9b, asFactor = T)

#GLMM Poisson model using log link + (1|individual) for decreasing dispersion
mod9c <- glmer(data = phen9, flwr_clstr ~ treatment + (1|plot) + (1|individual), family = "poisson")
check_overdispersion(mod9c) #no overdispersion

  #Model Diagnostics
testDispersion(mod9c)
sim_mod9c <- simulateResiduals(fittedModel = mod9c) #calculating residuals using DHARMa package
plot(sim_mod9c, asFactor = T) #Change asFactor = T to scale across random effects; #someone improved dispersion
testOutliers(sim_mod9c, type="bootstrap")
VarCorr(mod9c)
anova(mod9b, mod9c)

#GLMM Neg-Bin model using log link
mod9d <- glmmTMB(data = phen9, flwr_clstr ~ treatment + (1|plot) + (1|individual), family = "nbinom1") # 
sim_mod9d <- simulateResiduals(fittedModel = mod9d) #still overdispersion
plot(sim_mod9d, asFactor = T)
testDispersion(mod9d)
#Can't apply zero-inflation to the model because most values are 1 (and have no zeroes anyway)


#Zero truncated Poisson model
mod9e <- glmmTMB(flwr_clstr ~ treatment + (1|plot), ziformula = ~0, family = "poisson", data = phen9)
sim_mod9e <- simulateResiduals(fittedModel = mod9e) #still overdispersion
plot(sim_mod9e, asFactor = T)
testDispersion(mod9e)


  #Test of Significance
summary(mod9c)
Anova(mod9c)

  #Getting estimated values per treatment
newdata9 <- data.frame(treatment = factor(1:2, levels = 1:2, labels = levels(phen9$treatment)), plot=NA)
pred9 <- predict(mod9c, newdata9, type = "response", re.form = NA) #NULL to include all random effects
pred9
#Ambient estimate = 1.704752
#Heated estimate = 2.502231

  #Standard error estimates
newdata9 <- data.frame(treatment = factor(levels(phen9$treatment), levels = levels(phen9$treatment)), plot = NA)
pred9b <- bootMer(mod9b, FUN=function(x) predict(x, re.form = NA, newdata = newdata9, type = "response"),
                 nsim = 1000) #standard errors with bootstrap 1000 simulations
pred9b
easyPredCI(mod9b, newdata9) #95% confidence intervals















## Trait 10: Stem Diameter ####
dat %>% group_by(treatment) %>% filter(germ==1, !stem_diam==0) %>% 
  summarise(mean=mean(stem_diam), mode=mode(stem_diam), st.dev=sd(stem_diam), var=(st.dev)^2)
phen10 <- dat %>% filter(germ==1, !stem_diam==0) #Removing 0s because previously NAs

hist(phen10$stem_diam, breaks=30) 
#inverse Gaussian or gamma distribution? These aren't counts (non-integers), so can't use Poisson.
#Perhaps use log transform then use LMM
mod10a <- lm(height ~ treatment, data = phen10)
qqnorm(rstandard(mod10a))
qqline(rstandard(mod10a)) #Residuals fairly normal
gg_reshist(mod10a, bins = 30) #Residuals somewhat normal
check_model(mod10a)

mod10 <- lmer(data = phen10, log(stem_diam) ~ treatment + (1|plot)) 
summary(mod10)
check_model(mod10)
Anova(mod10) #for p values
exp(0.64513) #Ambient estimate = 1.906235
exp(0.64513 + 0.13868) #Heated estimate = 2.1898
exp(0.08181) #Ambient standard error = 1.08525
exp(0.08181 + 0.11573) #Heated standard error = 1.218402 

