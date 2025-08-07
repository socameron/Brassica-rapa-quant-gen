#### PROJECT: Brassica rapa Va/W Study
#### PURPOSE: Simple model to estimate broad-sense genetic variance

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Necessary Packages
library(tidyverse)
library(lme4)
#library(gtools) #for logit, can't use car
library(car)
library(MASS)
library(vcd)
library(performance)
library(see)
library(merTools)
library(lmerTest)
library(glmmTMB)
library(DHARMa)


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

#Importing data frame 'dat' from the cleaning process of 01_Data_exploration
load(file="Routput/heatarrays_data_explore.RData")
dat <- dat %>%
  mutate(leaf.st = ifelse(treatment == "A", scale(leaf), scale(leaf)),
         flwr_clstr.st = ifelse(treatment == "A", scale(flwr_clstr), scale(flwr_clstr)),
         seed_pods.st = ifelse(treatment == "A", scale(seed_pods), scale(seed_pods)),
         height.st = ifelse(treatment == "A", scale(height), scale(height)),
         stem_diam.st = ifelse(treatment == "A", scale(stem_diam), scale(stem_diam)),
         germ_census.st = ifelse(treatment == "A", scale(germ_census), scale(germ_census)),
         flwr_census.st = ifelse(treatment == "A", scale(flwr_census), scale(flwr_census)))


#'####################################################################'#
#############     BROAD-SENSE HERITABILITY MODELS        ##############
#'####################################################################'#

#NOTE: 
# (1) We exclude plot effect variance
# (2) Broad-sense heritability is defined as 2y/(2y + e) where
#y = genetic variance (random effect famID)
#e = residual variance (including environmental variance)
# (3) We extract genetic variance by each environment

## Trait 1: Leaf Number ####
dat %>% group_by(treatment) %>% filter(germ==1) %>% 
  summarise(mean=mean(leaf), st.dev=sd(leaf), var=(st.dev)^2, CV=st.dev/mean)

phen1 <- dat %>% filter(germ==1)
leaf_mod3 <- lmer(leaf ~ treatment + (1|famID) + (1|famID:treatment) + (1|plot), 
                  data = phen1)
leaf_mod2 <- update(leaf_mod3, .~ treatment + (1|famID) + (1|plot))
leaf_mod1 <- update(leaf_mod3, .~ treatment + (1|plot), data = phen1)

#Diagnostics
check_heteroscedasticity(leaf_mod3) #heteroscedastic
check_heteroscedasticity(leaf_mod2) #heteroscedastic
check_heteroscedasticity(leaf_mod1) #heteroscedastic

hist((phen1$leaf), breaks=50) #Gaussian distribution
mod1a <- lm(leaf ~ treatment, data = phen1)
plot(mod1a) #Residuals fairly normal, although some tailing
#We ignore tests of heteroscedasticity because of the Q-Q plots

#Model comparison
anova(leaf_mod3, leaf_mod2, leaf_mod1) #2nd model is significant - significant Vg
summary(leaf_mod2) #Vg
confint(leaf_mod2, oldNames = FALSE)
0.5846387^2 #2.5% Vg
0.9081174^2 #97.5% Vg

#Trait mean extraction by treatment
newdata1 <- data.frame(treatment = factor(1:2, levels = 1:2, labels = levels(phen1$treatment)), plot=NA)
pred1 <- predict(leaf_mod2, newdata1, type = "response", re.form = NA) #NULL to include all random effects
pred1
#Ambient estimate = 6.426610
#Heated estimate = 7.744185
6.4266 #Alternatively, extract from summary(leaf_mod2) .. Ambient = 6.4266
6.4266 + 1.3176 #Heated = 7.7442

#Standard error estimates
summary(leaf_mod2)

#P-value for treatment/environment
Anova(leaf_mod2)

#Vg by environment
#Ambient environment
AM_phen1 <- dat %>% filter(germ==1, treatment=="A")
AM_leaf <- lmer(leaf ~ plot + (1|famID), data = AM_phen1)
check_heteroscedasticity(AM_leaf) #homoscedastic
summary(AM_leaf)
AM_leaf_var <- as.data.frame(VarCorr(AM_leaf))
AM_leaf_Vg <- AM_leaf_var[1,'vcov']
AM_leaf_Vg
confint(AM_leaf, oldNames = FALSE)
0.48360798^2 #2.5% Vg CI
0.8397098^2 #97.5% Vg CI

AM_leaf_Vr <- AM_leaf_var[2,'vcov']
AM_leaf_Vr
AM_leaf_H2 <- AM_leaf_Vg / (AM_leaf_Vg + AM_leaf_Vr) #Broad-sense H2
AM_leaf_H2
(0.48360798^2) / ( (0.48360798^2) + (2.20007268^2) )#2.5% CI
(0.8397098^2) / ( (0.8397098^2) + (2.3619782^2) ) #97.5% CI

#Checking Vp
AM_phen1 %>% summarise(st.dev=sd(leaf), var=(st.dev)^2) #Raw Vp = 5.713137
(AM_leaf_Vg + AM_leaf_Vr) #Vp =  5.632654


#Heated environment
HW_phen1 <- dat %>% filter(germ==1, treatment=="H")
HW_leaf <- lmer(leaf ~ plot + (1|famID), data = HW_phen1)
plot(HW_leaf)
check_heteroscedasticity(HW_leaf) #heteroscedastic
hist((HW_phen1$leaf), breaks=50) # fairly normal
summary(HW_leaf)
HW_leaf_var <- as.data.frame(VarCorr(HW_leaf))
HW_leaf_Vg <- HW_leaf_var[1,'vcov']
HW_leaf_Vg
confint(HW_leaf, oldNames = FALSE)
0.6058223^2 #2.5% Vg
0.9990391^2 #97.5% Vg

HW_leaf_Vr <- HW_leaf_var[2,'vcov']
HW_leaf_H2 <- HW_leaf_Vg / (HW_leaf_Vg + HW_leaf_Vr) #Broad-sense H2
HW_leaf_H2
(0.6058223^2) / ( (0.6058223^2) + (2.3763832^2) )#2.5% CI
(0.9990391^2) / ( (0.9990391^2) + (2.5519727^2) ) #97.5% CI

#Checking Vp
HW_phen1 %>% summarise(st.dev=sd(leaf), var=(st.dev)^2) #Raw Vp = 6.931959
(HW_leaf_Vg + HW_leaf_Vr) #Vp = 6.701281
#









## Trait 2: Height ####
dat %>% group_by(treatment) %>% filter(germ==1, !height==0) %>% 
  summarise(mean=mean(height), st.dev=sd(height), var=(st.dev)^2, CV=st.dev/mean)

phen2 <- dat %>% filter(germ==1, !height==0)
height_mod3 <- lmer(height ~ treatment + (1|famID) + (1|famID:treatment) + (1|plot), 
                    data = phen2)
height_mod2 <- update(height_mod3, .~ treatment + (1|famID) + (1|plot))
height_mod1 <- update(height_mod3, .~ treatment + (1|plot), data = phen2)

#Diagnostics
check_heteroscedasticity(height_mod3) #heteroscedastic
check_heteroscedasticity(height_mod2) #heteroscedastic
check_heteroscedasticity(height_mod1) #heteroscedastic

hist(phen2$height, breaks=50) #Gaussian distribution
mod2a <- lm(height ~ treatment, data = phen2)
plot(mod2a) #Residuals fairly normal ... some tailing but better than transforming the data
#We ignore tests of heteroscedasticity because of the Q-Q plots

#Model comparison
anova(height_mod3, height_mod2, height_mod1) #3nd model is significant - significant Vg + GxE
summary(height_mod3) #Vg + Vgxe
confint(height_mod3, oldNames = FALSE)
1.335368^2 #2.5% Vg
4.550886^2 #97.5% Vg
2.487820^2 #2.5% Vgxe
5.022777^2 #97.5% Vgxe

#Trait mean extraction by treatment
newdata2 <- data.frame(treatment = factor(1:2, levels = 1:2, labels = levels(phen2$treatment)), 
                       plot=NA, famID=NA)
pred2 <- predict(height_mod3, newdata2, type = "response", re.form = NA) #NULL to include all random effects
pred2
#Field estimate = 49.05665
#Greenhouse estimate = 30.61226
49.06 #Alternatively, extract from summary .. Greenhouse = 49.06
49.06 + (-18.44) #Heated = 30.62

#Standard error of mean estimates
summary(height_mod3)

#P-value for treatment/environment
Anova(height_mod3)


#Vg by environment
#Greenhouse environment
AM_phen2 <- dat %>% filter(germ==1, !height==0, treatment=="G") #Removing 0s because previously NAs
AM_height <- lmer(height ~ plot + (1|famID), data = AM_phen2)
check_heteroscedasticity(AM_height)

#Model diagnostics
hist((AM_phen2$height), breaks=50) #left-skewed data
AM_phen2a <- lm(height ~ plot, data = AM_phen2) 
plot(AM_phen2a) #some tailing.. and transformations dont improve fit.
#We choose to ignore some of the violations... simply because we are only extracting the variance not doing any predictions/inferences

summary(AM_height)
AM_height_var <- as.data.frame(VarCorr(AM_height))
AM_height_Vg <- AM_height_var[1,'vcov']
AM_height_Vg
confint(AM_height, oldNames = FALSE)
4.3108965^2 #2.5% Vg
8.0905895^2 #97.5% Vg

AM_height_Vr <- AM_height_var[2,'vcov']
AM_height_H2 <- AM_height_Vg / (AM_height_Vg + AM_height_Vr) #Broad-sense H2
AM_height_H2
(4.3108965^2) / ( (4.3108965^2) + (14.8192568^2) )#2.5% CI
(8.09058950^2) / ( (8.09058950^2) + (16.61764062^2) ) #97.5% CI

#Checking Vp
AM_phen2 %>% summarise(st.dev=sd(height), var=(st.dev)^2) #Raw Vp = 287
(AM_height_Vg + AM_height_Vr) #Vp = 285

#Field environment
HW_phen2 <- dat %>% filter(germ==1, !height==0, treatment=="A") #Removing 0s because previously NAs
HW_height <- lmer(height ~ plot + (1|famID), data = HW_phen2)
check_heteroscedasticity(HW_height) #heteroscedastic

#Model diagnostics
hist(HW_phen2$height, breaks=50) #right-skewed data
AM_phen2a <- lm(height ~ plot, data = AM_phen2) 
plot(AM_phen2a) #some tailing.. log transformations removes heteroscedasticity but not severely different
#We choose to ignore some of the violations... simply because we are only extracting the variance not doing any predictions/inferences

summary(HW_height)
HW_height_var <- as.data.frame(VarCorr(HW_height))
HW_height_Vg <- HW_height_var[1,'vcov']
HW_height_Vg
confint(HW_height, oldNames = FALSE)
2.6423117^2 #2.5% Vg
4.498921^2 #97.5% Vg

HW_height_Vr <- HW_height_var[2,'vcov']
HW_height_H2 <- HW_height_Vg / (HW_height_Vg + HW_height_Vr) #Broad-sense H2
HW_height_H2
(2.6423117^2) / ( (2.6423117^2) + (11.0026985^2) )#2.5% CI
(4.498921^2) / ( (4.498921^2) + (11.854957^2) ) #97.5% CI

#Checking Vp
HW_phen2 %>% summarise(st.dev=sd(height), var=(st.dev)^2) #Raw Vp = 148
(HW_height_Vg + HW_height_Vr) #Vp = 143.1255
#










## Trait 3: Flowering Clusters Number ####
dat %>% group_by(treatment) %>% filter(germ==1, flower==1) %>% 
  summarise(mean=mean(flwr_clstr), st.dev=sd(flwr_clstr), var=(st.dev)^2, CV=st.dev/mean)

phen3 <- dat %>% filter(germ==1, flower==1)
flwrclstr_mod3 <- glmer(flwr_clstr ~ treatment + (1|famID) + (1|famID:treatment) + (1|plot), 
                        family = "poisson", data = phen3)
flwrclstr_mod2 <- update(flwrclstr_mod3, .~ treatment + (1|famID) + (1|plot))
flwrclstr_mod1 <- update(flwrclstr_mod3, .~ treatment + (1|plot))

#Diagnostics
check_overdispersion(flwrclstr_mod3) #no overdispersion
check_overdispersion(flwrclstr_mod2) #no overdispersion
check_overdispersion(flwrclstr_mod1) #marginal levels of dispersion

#Model Comparison
anova(flwrclstr_mod3, flwrclstr_mod2, flwrclstr_mod1) #3rd model is significant - significant Vg + GxE
summary(flwrclstr_mod3) #Vg + Vgxe
confint(flwrclstr_mod3, oldNames = FALSE)
0^2 #2.5% Vg
0.1451543^2 #97.5% Vg
0.07351783^2 #2.5% Vgxe
0.1853064^2 #97.5% Vgxe

#Treatment mean extraction by treatment
newdata3 <- data.frame(treatment = factor(1:2, levels = 1:2, labels = levels(phen3$treatment)), plot=NA)
pred3 <- predict(flwrclstr_mod3, newdata3, type = "response", re.form = NA) #NULL to include all random effects
pred3
#Field estimate = 1.883104
#Greenhouse estimate = 3.571158

#Standard error estimates
newdata3 <- data.frame(treatment = factor(levels(phen3$treatment), levels = levels(phen3$treatment)), plot = NA)
pred3b <- bootMer(flwrclstr_mod3, FUN=function(x) predict(x, re.form = NA, newdata = newdata3, type = "response"),
                  nsim = 1000) #standard errors with bootstrap 1000 simulations
pred3b #some sims didn't converge but was less than 100
easyPredCI(mod3b, newdata3) #95% confidence intervals

#P-value for treatment/environment
Anova(flwrclstr_mod3)

#Vg by environment
#Greenhouse environment
AM_phen3 <- dat %>% filter(germ==1, flower==1, treatment=="G")
AM_flwrclstr <- glmer(flwr_clstr ~ plot + (1|famID), family = "poisson", data = AM_phen3)
check_overdispersion(AM_flwrclstr) #no overdispersion
summary(AM_flwrclstr)
AM_flwrclstr_var <- as.data.frame(VarCorr(AM_flwrclstr))
AM_flwrclstr_Vg <- AM_flwrclstr_var[1,'vcov']
AM_flwrclstr_Vg
confint(AM_flwrclstr, oldNames = FALSE)
0.04714784^2 #2.5% Vg
0.18164403^2 #97.5% Vg

AM_phen3b <- AM_phen3 %>% 
  mutate(flwr_clstr1 = log(flwr_clstr+1)) %>%
  summarise(st.dev=sd(flwr_clstr1), var=(st.dev)^2)
AM_flwrclstr_Vp <- AM_phen3b[1,'var'] #Vr
AM_flwrclstr_Vp
AM_flwrclstr_H2 <- AM_flwrclstr_Vg / AM_flwrclstr_Vp #Broad-sense H2
AM_flwrclstr_H2


#Field environment
HW_phen3 <- dat %>% filter(germ==1, flower==1, treatment=="A")
HW_flwrclstr <- glmer(flwr_clstr ~ plot + (1|famID), family = "poisson", data = HW_phen3)
check_overdispersion(HW_flwrclstr) #overdispersion found, but very marginal.. so it is acceptable
summary(HW_flwrclstr)
HW_flwrclstr_var <- as.data.frame(VarCorr(HW_flwrclstr))
HW_flwrclstr_Vg <- HW_flwrclstr_var[1,'vcov']
HW_flwrclstr_Vg
confint(HW_flwrclstr, oldNames = FALSE)
0.12411269^2 #2.5% Vg
0.22867652^2 #97.5% Vg

HW_phen3b <- HW_phen3 %>% 
  mutate(flwr_clstr1 = log(flwr_clstr+1)) %>%
  summarise(st.dev=sd(flwr_clstr1), var=(st.dev)^2)
HW_flwrclstr_Vp <- HW_phen3b[1,'var'] #Vp
HW_flwrclstr_Vp
HW_flwrclstr_H2 <- HW_flwrclstr_Vg / HW_flwrclstr_Vp #Broad-sense H2
HW_flwrclstr_H2

#










## Trait 4: Fecundity ####
dat %>% group_by(treatment) %>% filter(germ==1, flower==1) %>% 
  summarise(mean=mean(seed_pods), st.dev=sd(seed_pods), var=(st.dev)^2, CV=st.dev/mean)

fHW_var <-dat %>% group_by(treatment, famID) %>%
  summarise(mean=mean(seed_pods), var=var(seed_pods))
fHW_var

phen4 <- dat %>% filter(germ==1, flower==1)
fecund_mod3 <- glmmTMB(seed_pods ~ treatment + (1|famID) + (1|famID:treatment) + (1|plot), 
                       family = nbinom2, data = phen4)
fecund_mod2 <- glmmTMB(seed_pods ~ treatment + (1|famID) + (1|plot), 
                       family = nbinom2, data = phen4)
fecund_mod1 <- glmmTMB(seed_pods ~ treatment + (1|plot), 
                       family = nbinom2, data = phen4)

#DHARMa diagnostics
sim_output3 <- simulateResiduals(fittedModel = fecund_mod3)
plot(sim_output3) # based on QQ plot residuals, no overdispersion detected
plotResiduals(sim_output3, phen4$treatment, quantreg = T)

sim_output2 <- simulateResiduals(fittedModel = fecund_mod2)
plot(sim_output2) # based on QQ plot residuals, no overdispersion detected

sim_output1 <- simulateResiduals(fittedModel = fecund_mod1)
plot(sim_output1)

testDispersion(fecund_mod3) # no overdispersion
testDispersion(fecund_mod2) # no overdispersion
testDispersion(fecund_mod1) # no overdispersion, but moderately underdispersed

testZeroInflation(fecund_mod3) # no zero inflation
testZeroInflation(fecund_mod2) # no zero inflation
testZeroInflation(fecund_mod3) # no zero inlation

#Model comparison - conservative because of underdispersion
anova(fecund_mod3, fecund_mod2) #3nd model is significant - significant Vg + GxE
summary(fecund_mod3) #Vg + Vgxe
confint(fecund_mod3)
0.04269934^2 #2.5% Vg
0.3520154^2 #97.5% Vg
0.11802783^2 #2.5% Vgxe
0.3129645^2 #97.5% Vgxe

#Extracting trait means by environment
newdata4 <- data.frame(treatment = factor(1:2, levels = 1:2, labels = levels(phen4$treatment)),
                       plot=NA, famID = NA)
pred4 <- predict(fecund_mod3, newdata4, type = "response", se.fit = TRUE, re.form = NA) #NULL to include all random effects
pred4
#Field estimate = 11.37696
#Greenhouse estimate = 59.72314
fit4 <- pred4$fit
lwr4 <- fit4 - 1.96*pred4$se.fit #95% confidence intervals (upper) .. first ambient, then heated
upr4 <- fit4 + 1.96*pred4$se.fit #95% confidence intervals (lower) .. first ambient, then heated
lwr4
upr4

#P-value for treatment/environment
Anova(fecund_mod3)

#Vg by environment
#Greenhouse environment
AM_phen4 <- dat %>% filter(germ==1, flower==1, treatment=="G")
AM_fecund <- glmmTMB(seed_pods ~ plot + (1|famID), 
                     ziformula = ~1, family = nbinom2, data = AM_phen4)

#Model diagnostics
hist(AM_phen4$seed_pods, breaks=50) #Poisson distribution with some inflation at the 0 end
testDispersion(AM_fecund) #no overdispersion
testZeroInflation(AM_fecund) #no zero-inflation

summary(AM_fecund) #famID variance = 0.006321
AM_fecund_var <- VarCorr(AM_fecund)
AM_fecund_Vg <- 0.079508^2
confint(AM_fecund) #outputs confidence intervals of famID st.dev
0.01950968^2 #95% confidence intervals, 2.5% = 0.0003806276
0.32401649^2 #97.5% = 0.1049867

#merTools::REsim(AM_fecund, n.sims=1000)
#cV <- ranef(AM_flwrclstr, condVar = TRUE) 

AM_phen4b <- AM_phen4 %>% 
  mutate(seed_pods1 = log(seed_pods+1)) %>%
  summarise(st.dev=sd(seed_pods), var=(st.dev)^2)
AM_fecund_Vp <- AM_phen4b[1,'var'] #Vp
AM_fecund_Vp
AM_fecund_H2 <- AM_fecund_Vg / AM_fecund_Vp #Broad-sense H2
AM_fecund_H2

#Field environment
HW_phen4 <- dat %>% filter(germ==1, flower==1, treatment=="A")
HW_fecund <- glmmTMB(seed_pods ~ plot + (1|famID), family = "nbinom2",
                     ziformula = ~1, data = HW_phen4)

#Model diagnostics
hist(HW_phen4$seed_pods, breaks=100) #Poisson distribution with some zero inflation
testDispersion(HW_fecund) #overdispersion at 1.4 ... this is not that acceptable at 1.1, but since we're not doing any significance testing, we'll take this.
testZeroInflation(HW_fecund) #somewhat underfitting the zeroes.

summary(HW_fecund)
HW_fecund_var <- VarCorr(HW_fecund)
HW_fecund_Vg <- 0.28844^2
confint(HW_fecund)
0.2188336^2 #2.5% Vg
0.3801931^2 #97.5% Vg

HW_phen4b <- HW_phen4 %>% 
  mutate(seed_pods1 = log(seed_pods+1)) %>%
  summarise(st.dev=sd(seed_pods), var=(st.dev)^2)
HW_fecund_Vp <- HW_phen4b[1,'var'] #Vp
HW_fecund_Vp
HW_fecund_H2 <- HW_fecund_Vg / HW_fecund_Vp #Broad-sense H2
HW_fecund_H2
#
