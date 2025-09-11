#### PROJECT: Brassica rapa GxE Study (Data collected at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Figures on plasticity and full models excluding dominance


## NOTICE ##
#We previously ran full models including dominance in MCMC_full_analysis, but these models failed to converge.
#In an earlier analysis, we also ran full models *excluding* dominance. To see if we can estimate genetic correlations from separate environment models (e.g from scripts 03 and 04), we test if genetic correlations agree from full models *excluding* dominance and from separate environment models. 
#This script just includes the full models excluding dominance. All genetic correlation estimations are found in 05_MCMC_genetic_correlations


#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Necessary Packages
library(tidyverse)
library(ggpubr)

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
  labs(x="Treatment", y="Mean Family Lifetime Fitness") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 40, 10)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.absolute2 <- fam.absolute %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.absolute2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5)  +
  labs(x="Ambient Absolute Fitness", y="Heated Absolute Fitness")



#Survival
fam.surv <- dat %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(flower), n=n())
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
fam.surv2 <- fam.surv %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.surv2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5)  +
  labs(x="Ambient Survival", y="Heated Survival")
  


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
fam.fecund2 <- fam.fecund %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.fecund2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient Fecundity", y="Heated Fecundity")



#Overwintering Survival to Germination (Germination Success)
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
fam.germ2 <- fam.germ %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.germ2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient Germination Success", y="Heated Germination Success")



#Spring-Summer Survival to Flowering (Flowering Success)
fam.flower <- dat %>%
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(flower), n=n()) %>%
  mutate(fam.mean=(sum/n), se.p=(sqrt(fam.mean*(1-fam.mean))/n))
fam.flower %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Spring-Summer Survival Success") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0.5, 1, 0.1)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.flower2 <- fam.flower %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.flower2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5)  +
  labs(x="Ambient Flowering Success", y="Heated Flowering Success")



#Fruiting Success
fam.seed <- dat %>%
  filter(germ==1, flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(seed), n=n()) %>%
  mutate(fam.mean=(sum/n), se.p=(sqrt(fam.mean*(1-fam.mean))/n))
fam.seed %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Fruiting Success") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0.9, 1, 0.02)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.seed2 <- fam.seed %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.seed2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient Fruiting Success", y="Heated Fruiting Success")



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
fam.flwrclstr2 <- fam.flwrclstr %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.flwrclstr2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient Flower Clusters", y="Heated Flower Clusters")



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
fam.leaf2 <- fam.leaf %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.leaf2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient Leaf", y="Heated Leaf")



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
fam.height2 <- fam.height %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.height2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient Height", y="Heated Height")



#Stem Diameter
fam.stem <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(stem_diam))
fam.stem %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Treatment", y="Mean Stem Diameter (mm)") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 3, 0.5)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.stem2 <- fam.stem %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.stem2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5)



# Additional Traits - Not included in Thesis #

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
############      "BROAD SENSE" GENETIC CORRELATIONS     #############
#'####################################################################'#

#Lifetime Fitness
fam.absolute2 <- pivot_wider(fam.absolute, id_cols = c(treatment, famID, fam.mean),
                             names_from = treatment, values_from = fam.mean) #Change data format
fam.absolute2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Lifetime Fitness")
qqnorm(fam.absolute2$A) 
qqline(fam.absolute2$A) #Ambient variable non-normal
hist(fam.absolute2$A, breaks=50)
qqnorm(fam.absolute2$H)
qqline(fam.absolute2$H) #Heated variable non-normal
hist(fam.absolute2$H, breaks=50)
shapiro.test(fam.absolute2$A) #Non-normal.. but heavily conservative test and sensitive to outliers
shapiro.test(fam.absolute2$H) #Non-normal.. but ""

cor.test(fam.absolute2$A, fam.absolute2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.absolute2$A, fam.absolute2$H, method="spearman") #Spearman ranked test (non-parametric)
#RESULT: Spearman's Rho = 0.6021254 .. relatively strong correlation


#Survival
fam.surv2 <- pivot_wider(fam.surv, id_cols = c(treatment, famID, fam.mean), 
                         names_from= treatment, values_from = fam.mean) #Change data format
fam.surv2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Survival")
qqnorm(fam.surv2$A) 
qqline(fam.surv2$A) #Ambient variable normal
hist(fam.surv2$A, breaks=50)
qqnorm(fam.surv2$H)
qqline(fam.surv2$H) #Heated variable non-normal slightly
hist(fam.surv2$H, breaks=50)
shapiro.test(fam.surv2$A) #Normal
shapiro.test(fam.surv2$H) #Normal

cor.test(fam.surv2$A, fam.surv2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.surv2$A, fam.surv2$H, method="spearman") #Spearman ranked test (non-parametric)
#RESULT: Pearson's r coeff = 0.580665 .. relatively strong correlation | R^2 = 0.3371736


#Fecundity
fam.fecund2 <- pivot_wider(fam.fecund, id_cols = c(treatment, famID, fam.mean), 
                           names_from= treatment, values_from = fam.mean) #Change data format
fam.fecund2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Fecundity")
qqnorm(fam.fecund2$A) 
qqline(fam.fecund2$A) #Ambient variable non-normal
hist(fam.fecund2$A, breaks=50)
qqnorm(fam.fecund2$H)
qqline(fam.fecund2$H) #Heated variable non-normal
hist(fam.fecund2$H, breaks=50)
shapiro.test(fam.fecund2$A) #Non-normal.. but heavily conservative test and sensitive to outliers
shapiro.test(fam.fecund2$H) #Non-normal.. but ""

cor.test(fam.fecund2$A, fam.fecund2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.fecund2$A, fam.fecund2$H, method="spearman") #Spearman ranked test (non-parametric)
#RESULT: Spearman's Rho = 0.518471 .. moderate correlation


#Overwinter Survival to Germination (Germination Success)
fam.germ2 <- pivot_wider(fam.germ, id_cols = c(treatment, famID, fam.mean), 
                         names_from= treatment, values_from = fam.mean) #Change data format
fam.germ2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Overwinter Survival")
qqnorm(fam.germ2$A) 
qqline(fam.germ2$A) #Ambient variable normal
hist(fam.germ2$A, breaks=50)
qqnorm(fam.germ2$H)
qqline(fam.germ2$H) #Heated variable normal
hist(fam.germ2$H, breaks=50)
shapiro.test(fam.germ2$A) #Normal.. but barely
shapiro.test(fam.germ2$H) #Normal

cor.test(fam.germ2$A, fam.germ2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.germ2$A, fam.germ2$H, method="spearman") #Spearman ranked test (non-parametric)
#RESULT: Pearson's r coeff = 0.4999616 .. weak correlation | R^2 = 0.2499616


#Spring-Summer Survival to Flowering (Flowering Success)
fam.flower2 <- pivot_wider(fam.flower, id_cols = c(treatment, famID, fam.mean), 
                         names_from= treatment, values_from = fam.mean) #Change data format
fam.flower2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Spring-Summer Survival")
qqnorm(fam.flower2$A) 
qqline(fam.flower2$A) #Ambient variable normal (barely)
hist(fam.flower2$A, breaks=50)
qqnorm(fam.flower2$H)
qqline(fam.flower2$H) #Heated variable non-normal
hist(fam.flower2$H, breaks=50)
shapiro.test(fam.flower2$A) #Non-normal.. but heavily conservative test and sensitive to outliers
shapiro.test(fam.flower2$H) #Non-normal.. but ""

cor.test(fam.flower2$A, fam.flower2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.flower2$A, fam.flower2$H, method="spearman", exact=FALSE) #Spearman ranked test (non-parametric), exact=FALSE because there are ranked ties in the heated (multiple 100% flowering success)
#RESULT: Spearman's Rho = 0.1565305 .. relatively weak correlation


#Seed Maturation Success
fam.seed2 <- pivot_wider(fam.seed, id_cols = c(treatment, famID, fam.mean), 
                         names_from= treatment, values_from = fam.mean) #Change data format
fam.seed2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Seed Maturation Success")
qqnorm(fam.seed2$A) 
qqline(fam.seed2$A) #Ambient variable non-normal
hist(fam.seed2$A, breaks=50)
qqnorm(fam.seed2$H)
qqline(fam.seed2$H) #Heated variable non-normal
hist(fam.seed2$H, breaks=50)
shapiro.test(fam.seed2$A) #Non-normal..
shapiro.test(fam.seed2$H) #Non-normal..

cor.test(fam.seed2$A, fam.seed2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.seed2$A, fam.seed2$H, method="spearman", exact=FALSE) #Spearman ranked test (non-parametric)
#RESULT: Spearman's Rho = 0.04257557.. no correlation


#Flowering Clusters
fam.flwrclstr2 <- pivot_wider(fam.flwrclstr, id_cols = c(treatment, famID, fam.mean), 
                              names_from= treatment, values_from = fam.mean) #Change data format
fam.flwrclstr2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Number of Flowering Clusters")
qqnorm(fam.flwrclstr2$A) 
qqline(fam.flwrclstr2$A) #Ambient variable non-normal
hist(fam.flwrclstr2$A, breaks=50)
qqnorm(fam.flwrclstr2$H)
qqline(fam.flwrclstr2$H) #Heated variable non-normal
hist(fam.flwrclstr2$H, breaks=50)
shapiro.test(fam.flwrclstr2$A) #Non-normal.. but heavily conservative test and sensitive to outliers
shapiro.test(fam.flwrclstr2$H) #Non-normal.. but ""

cor.test(fam.flwrclstr2$A, fam.flwrclstr2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.flwrclstr2$A, fam.flwrclstr2$H, method="spearman", exact=FALSE) #Spearman ranked test (non-parametric)
#RESULT: Spearman's Rho = 0.6244979 .. relatively strong correlation


#Leaf Number
fam.leaf2 <- pivot_wider(fam.leaf, id_cols = c(treatment, famID, fam.mean), 
                         names_from= treatment, values_from = fam.mean) #Change data format
fam.leaf2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Leaf Number")
qqnorm(fam.leaf2$A) 
qqline(fam.leaf2$A) #Ambient variable normal
hist(fam.leaf2$A, breaks=50)
qqnorm(fam.leaf2$H)
qqline(fam.leaf2$H) #Heated variable normal
hist(fam.leaf2$H, breaks=50)
shapiro.test(fam.leaf2$A) #Normal.. 
shapiro.test(fam.leaf2$H) #Normal.. 

cor.test(fam.leaf2$A, fam.leaf2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.leaf2$A, fam.leaf2$H, method="spearman") #Spearman ranked test (non-parametric)
#RESULT: Pearson's r coeff = 0.63964.. relatively strong correlation | R^2 = 0.4091393


#Height
fam.height2 <- pivot_wider(fam.height, id_cols = c(treatment, famID, fam.mean), 
                           names_from= treatment, values_from = fam.mean) #Change data format
fam.height2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Height")
qqnorm(fam.height2$A) 
qqline(fam.height2$A) #Ambient variable normal
hist(fam.height2$A, breaks=50)
qqnorm(fam.height2$H)
qqline(fam.height2$H) #Heated variable normal
hist(fam.height2$H, breaks=50)
shapiro.test(fam.height2$A) #Normal.. 
shapiro.test(fam.height2$H) #Normal.. 

cor.test(fam.height2$A, fam.height2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.height2$A, fam.height2$H, method="spearman") #Spearman ranked test (non-parametric)
#RESULT: Pearson's r coeff = 0.6418526 relatively strong correlation | R^2 = 0.4119748


#Stem Diameter
fam.stem2 <- pivot_wider(fam.stem, id_cols = c(treatment, famID, fam.mean), 
                         names_from= treatment, values_from = fam.mean) #Change data format
fam.stem2 %>% ggplot(aes(x=A, y=H)) + 
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  xlab("Phenotype - Ambient") + ylab("Phenotype - Heated") + ggtitle("Stem Diameter")
qqnorm(fam.stem2$A) 
qqline(fam.stem2$A) #Ambient variable normal
hist(fam.stem2$A, breaks=50)
qqnorm(fam.stem2$H)
qqline(fam.stem2$H) #Heated variable normal
hist(fam.stem2$H, breaks=50)
shapiro.test(fam.stem2$A) #Normal.. 
shapiro.test(fam.stem2$H) #Normal.. 

cor.test(fam.stem2$A, fam.stem2$H, method="pearson") #Pearson test (parametric)
cor.test(fam.stem2$A, fam.stem2$H, method="spearman") #Spearman ranked test (non-parametric)
#RESULT: Pearson's r coeff = 0.7094773 relatively strong correlation | R^2 = 0.503358

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
###############      FULL GENETIC MODELS ON TRAITS      ################
#'####################################################################'#

#Covariance in residuals is not identifiable since only one observation has one treatment applied
#Especially in binary models, the residual variance can't be identified.. so fix = 1 in R total_prior and rcov = ~units
#NOTE: **ALL OF THESE MODELS DID NOT CONVERGE WELL AND FAILED DIAGNOSTICS**. We keep the code here as a reference but they are not mentioned in the manuscript results. 

## Trait 1: Survival Success (Binary) ####
total_survival <- MCMC.2019
total_prior1.1 <- list(R = list(V = diag(1)*0.02, nu = 3, fix = 1),
                       G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)
                       ))

set.seed(1)
total_model1.1 <- mclapply(1:5, function(i) { 
  MCMCglmm(flower ~ plot:treatment + treatment, 
           random = ~us(treatment):animal + us(treatment):animalDom + us(treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov = ~units,
           family = "threshold", data = total_survival, prior = total_prior1.1, #Bernoulli distribution
           nitt = 500000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)
}, mc.cores = 5)

#Model 1.0 - Full model
load(file="Routput/Full/total_model1.0.RData")
t1.0_sol <- lapply(total_model1.0, function(m) m$Sol)
t1.0_sol <- do.call(mcmc.list, t1.0_sol)

t1.0_vcv <- lapply(total_model1.0, function(m) m$VCV)
t1.0_vcv <- do.call(mcmc.list, t1.0_vcv)
t1.0_vcv <- mcmc.stack(t1.0_vcv)
plot(t1.0_vcv) 
autocorr.diag(t1.0_vcv, relative=FALSE) 
heidel.diag(t1.0_vcv)
effectiveSize(t1.0_vcv) 

# Model 1.1 - Informative Priors
total_survival %>% group_by(treatment) %>%
  summarise(tot.flwr=sum(flower), n=n(), prob=tot.flwr/n, 
            var = n*prob*(1-prob), varlogit = logit(prob*(1-prob))) #varlogit = -1.42
#Removed all multipliers to V and alpha.V

load(file="Routput/Full/total_model1.1.RData")
t1.1_sol <- lapply(total_model1.1, function(m) m$Sol)
t1.1_sol <- do.call(mcmc.list, t1.1_sol)

t1.1_vcv <- lapply(total_model1.1, function(m) m$VCV)
t1.1_vcv <- do.call(mcmc.list, t1.1_vcv)
t1.1_vcv <- mcmc.stack(t1.1_vcv)
plot(t1.1_vcv) 
autocorr.diag(t1.1_vcv, relative=FALSE) 
heidel.diag(t1.1_vcv)
effectiveSize(t1.1_vcv) 

# Model 1.2 - Cross-environment Vd and Vm set to 0 (idh)
load(file="Routput/Full/total_model1.2.RData")
t1.2_sol <- lapply(total_model1.2, function(m) m$Sol)
t1.2_sol <- do.call(mcmc.list, t1.2_sol)

t1.2_vcv <- lapply(total_model1.2, function(m) m$VCV)
t1.2_vcv <- do.call(mcmc.list, t1.2_vcv)
t1.2_vcv <- mcmc.stack(t1.2_vcv)
plot(t1.2_vcv) 
autocorr.diag(t1.2_vcv, relative=FALSE) 
heidel.diag(t1.2_vcv)
effectiveSize(t1.2_vcv) 









## Trait 2: Fecundity (Poisson) ####
total_fecundity <- MCMC.2019 %>% filter(germ==1 & flower==1)
total_prior2.0 <- list(R = list(V = diag(2)*0.02, nu = 3),
                       G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)
                       ))

set.seed(1)
total_model2.0 <- mclapply(1:5, function(i) { 
  MCMCglmm(seed_pods ~ plot:treatment + treatment, 
           random = ~us(treatment):animal + us(treatment):animalDom + us(treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~idh(treatment):units,
           family = "poisson", data = total_fecundity, prior = total_prior2.0, #Poisson distribution
           nitt = 500000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
}, mc.cores = 5)

#Model 2.0 - Full model
load(file="Routput/Full/total_model2.0.2.RData")
t2.0_sol <- lapply(total_model2.0, function(m) m$Sol)
t2.0_sol <- do.call(mcmc.list, t2.0_sol)

t2.0_vcv <- lapply(total_model2.0, function(m) m$VCV)
t2.0_vcv <- do.call(mcmc.list, t2.0_vcv)
t2.0_vcv <- mcmc.stack(t2.0_vcv)
plot(t2.0_vcv) 
autocorr.diag(t2.0_vcv, relative=FALSE) 
heidel.diag(t2.0_vcv)
effectiveSize(t2.0_vcv)

#Model 2.1 - Informative Priors
total_fecundity %>% group_by(treatment) %>%
  summarise(sd=sd(seed_pods), var=(sd)^2, mean1=mean(seed_pods),
            logmu = log(mean1)) #log(mu) is ~2
#Had alpha.V=diag(2)*2 and all V = diag(2)

load(file="Routput/Full/total_model2.1.RData")
t2.1_sol <- lapply(total_model2.1, function(m) m$Sol)
t2.1_sol <- do.call(mcmc.list, t2.1_sol)

t2.1_vcv <- lapply(total_model2.1, function(m) m$VCV)
t2.1_vcv <- do.call(mcmc.list, t2.1_vcv)
t2.1_vcv <- mcmc.stack(t2.1_vcv)
plot(t2.1_vcv) 
autocorr.diag(t2.1_vcv, relative=FALSE) 
heidel.diag(t2.1_vcv)
effectiveSize(t2.1_vcv) 

# Model 2.2 - Cross-environment Vd and Vm set to 0 (idh)
load(file="Routput/Full/total_model2.2.RData")
t2.2_sol <- lapply(total_model2.2, function(m) m$Sol)
t2.2_sol <- do.call(mcmc.list, t2.2_sol)

t2.2_vcv <- lapply(total_model2.2, function(m) m$VCV)
t2.2_vcv <- do.call(mcmc.list, t2.2_vcv)
t2.2_vcv <- mcmc.stack(t2.2_vcv)
plot(t2.2_vcv) 
autocorr.diag(t2.2_vcv, relative=FALSE) 
heidel.diag(t2.2_vcv)
effectiveSize(t2.2_vcv)  










## Trait 3: Overwintering Survival / Germination Success (Binary) ####
total_germ <- MCMC.2019
total_prior3.1 <- list(R = list(V = diag(1), nu = 3, fix = 1),
                       G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G2 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G3 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)
                       ))

set.seed(1)
total_model3.1 <- mclapply(1:5, function(i) { 
  MCMCglmm(germ ~ plot:treatment + treatment, 
           random = ~us(treatment):animal + us(treatment):animalDom +  us(treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~units,
           family = "threshold", data = total_germ, prior = total_prior3.1, #Bernoulli distribution
           nitt = 500000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)
}, mc.cores = 5)

#Model 3.0 - Full model
load(file="Routput/Full/total_model3.0.RData")
t3.0_sol <- lapply(total_model3.0, function(m) m$Sol)
t3.0_sol <- do.call(mcmc.list, t3.0_sol)

t3.0_vcv <- lapply(total_model3.0, function(m) m$VCV)
t3.0_vcv <- do.call(mcmc.list, t3.0_vcv)
t3.0_vcv <- mcmc.stack(t3.0_vcv)
plot(t3.0_vcv) 
autocorr.diag(t3.0_vcv, relative=FALSE) 
heidel.diag(t3.0_vcv)
effectiveSize(t3.0_vcv) 

# Model 3.1 - Informative Priors
total_germ %>% group_by(treatment) %>%
  summarise(tot.germ=sum(germ), n=n(), prob=tot.germ/n, 
            var = n*prob*(1-prob), varlogit = logit(prob*(1-prob))) #varlogit = -1.10
#Removed all multipliers to V and alpha.V


load(file="Routput/Full/total_model3.1.RData")
t3.1_sol <- lapply(total_model3.1, function(m) m$Sol)
t3.1_sol <- do.call(mcmc.list, t3.1_sol)

t3.1_vcv <- lapply(total_model3.1, function(m) m$VCV)
t3.1_vcv <- do.call(mcmc.list, t3.1_vcv)
t3.1_vcv <- mcmc.stack(t3.1_vcv)
plot(t3.1_vcv) 
autocorr.diag(t3.1_vcv, relative=FALSE) 
heidel.diag(t3.1_vcv)
effectiveSize(t3.1_vcv) 

# Model 3.2 - Cross-environment Vd and Vm set to 0 (idh)
load(file="Routput/Full/total_model3.2.RData")
t3.2_sol <- lapply(total_model3.2, function(m) m$Sol)
t3.2_sol <- do.call(mcmc.list, t3.2_sol)

t3.2_vcv <- lapply(total_model3.2, function(m) m$VCV)
t3.2_vcv <- do.call(mcmc.list, t3.2_vcv)
t3.2_vcv <- mcmc.stack(t3.2_vcv)
plot(t3.2_vcv) 
autocorr.diag(t3.2_vcv, relative=FALSE) 
heidel.diag(t3.2_vcv)
effectiveSize(t3.2_vcv) # 











## Trait 4: Spring-summer Survival / Flowering Success (Binary) ####
total_flowering <- MCMC.2019 %>% filter(germ==1)
total_prior4.0 <- list(R = list(V = diag(1)*0.02, nu = 3, fix = 1),
                       G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)
                       ))

set.seed(1)
total_model4.0 <- mclapply(1:5, function(i) { 
  MCMCglmm(flower ~ plot:treatment + treatment, 
           random = ~us(treatment):animal + us(treatment):animalDom + us(treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~units,
           family = "threshold", data = total_flowering, prior = total_prior4.0, #Bernoulli distribution
           nitt = 500000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)
}, mc.cores = 5)

#Model 4.0 - Full model
load(file="Routput/Full/total_model4.0.RData")
t4.0_sol <- lapply(total_model4.0, function(m) m$Sol)
t4.0_sol <- do.call(mcmc.list, t4.0_sol)

t4.0_vcv <- lapply(total_model4.0, function(m) m$VCV)
t4.0_vcv <- do.call(mcmc.list, t4.0_vcv)
t4.0_vcv <- mcmc.stack(t4.0_vcv)
plot(t4.0_vcv) 
autocorr.diag(t4.0_vcv, relative=FALSE) 
heidel.diag(t4.0_vcv)
effectiveSize(t4.0_vcv) 

# Model 4.1 - Informative Priors
total_flowering %>% group_by(treatment) %>%
  summarise(tot.flower=sum(flower), n=n(), prob=tot.flower/n, 
            var = n*prob*(1-prob), varlogit = logit(prob*(1-prob))) #varlogit = -2.26

load(file="Routput/Full/total_model4.1.RData")
t4.1_sol <- lapply(total_model4.1, function(m) m$Sol)
t4.1_sol <- do.call(mcmc.list, t4.1_sol)

t4.1_vcv <- lapply(total_model4.1, function(m) m$VCV)
t4.1_vcv <- do.call(mcmc.list, t4.1_vcv)
t4.1_vcv <- mcmc.stack(t4.1_vcv)
plot(t4.1_vcv) 
autocorr.diag(t4.1_vcv, relative=FALSE) 
heidel.diag(t4.1_vcv)
effectiveSize(t4.1_vcv) 

# Model 4.2 - Cross-environment Vd and Vm set to 0 (idh)
load(file="Routput/Full/total_model4.2.RData")
t4.2_sol <- lapply(total_model4.2, function(m) m$Sol)
t4.2_sol <- do.call(mcmc.list, t4.2_sol)

t4.2_vcv <- lapply(total_model4.2, function(m) m$VCV)
t4.2_vcv <- do.call(mcmc.list, t4.2_vcv)
t4.2_vcv <- mcmc.stack(t4.2_vcv)
plot(t4.2_vcv) 
autocorr.diag(t4.2_vcv, relative=FALSE) 
heidel.diag(t4.2_vcv)
effectiveSize(t4.2_vcv) 











## Trait 5: Seed Maturation Success (Binary) ####
total_seed <- MCMC.2019 %>% filter(germ==1 & flower==1)
total_prior5.0 <- list(R = list(V = diag(1)*0.02, nu = 3, fix = 1),
                       G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)
                       ))

set.seed(1)
total_model5.0 <- mclapply(1:5, function(i) {
  MCMCglmm(seed ~ plot:treatment + treatment, 
           random = ~us(treatment):animal + us(treatment):animalDom + us(treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~units,
           family = "threshold", data = total_seed, prior = total_prior5.0, #Bernoulli distribution
           nitt = 500000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)
}, mc.cores = 5)

#Model 5.0 - Full model
load(file="Routput/Full/total_model5.0.RData")
t5.0_sol <- lapply(total_model5.0, function(m) m$Sol)
t5.0_sol <- do.call(mcmc.list, t5.0_sol)

t5.0_vcv <- lapply(total_model5.0, function(m) m$VCV)
t5.0_vcv <- do.call(mcmc.list, t5.0_vcv)
t5.0_vcv <- mcmc.stack(t5.0_vcv)
plot(t5.0_vcv) 
autocorr.diag(t5.0_vcv, relative=FALSE) 
heidel.diag(t5.0_vcv)
effectiveSize(t5.0_vcv) 

# Model 5.1 - Informative Priors
total_seed %>% group_by(treatment) %>%
  summarise(tot.seed=sum(seed), n=n(), prob=tot.seed/n, 
            var = n*prob*(1-prob), varlogit = logit(prob*(1-prob))) #varlogit = -1.42
#Removed all multipliers to V and alpha.V

load(file="Routput/Full/total_model5.1.RData")
t5.1_sol <- lapply(total_model5.1, function(m) m$Sol)
t5.1_sol <- do.call(mcmc.list, t5.1_sol)

t5.1_vcv <- lapply(total_model5.1, function(m) m$VCV)
t5.1_vcv <- do.call(mcmc.list, t5.1_vcv)
t5.1_vcv <- mcmc.stack(t5.1_vcv)
plot(t5.1_vcv) 
autocorr.diag(t5.1_vcv, relative=FALSE) 
heidel.diag(t5.1_vcv)
effectiveSize(t5.1_vcv) 

# Model 5.2 - Cross-environment Vd and Vm set to 0 (idh)
load(file="Routput/Full/total_model5.2.RData")
t5.2_sol <- lapply(total_model5.2, function(m) m$Sol)
t5.2_sol <- do.call(mcmc.list, t5.2_sol)

t5.2_vcv <- lapply(total_model5.2, function(m) m$VCV)
t5.2_vcv <- do.call(mcmc.list, t5.2_vcv)
t5.2_vcv <- mcmc.stack(t5.2_vcv)
plot(t5.2_vcv) 
autocorr.diag(t5.2_vcv, relative=FALSE) 
heidel.diag(t5.2_vcv)
effectiveSize(t5.2_vcv) 











## Trait 6: Leaf Number (Gaussian) ####
total_leaf <- MCMC.2019 %>% filter(germ==1)
total_prior6.0 <- list(R = list(V = diag(2)*0.02, nu = 3),
                       G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)
                       ))

set.seed(1)
total_model6.0 <- mclapply(1:5, function(i) {
  MCMCglmm(leaf ~ plot:treatment + treatment, 
           random = ~us(treatment):animal + us(treatment):animalDom + us(treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~idh(treatment):units,
           family = "gaussian", data = total_leaf, prior = total_prior6.0, #Gaussian distribution
           nitt = 500000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
}, mc.cores = 5)

#Model 6.0 - Full model
load(file="Routput/Full/total_model6.0.2.RData")
t6.0_sol <- lapply(total_model6.0, function(m) m$Sol)
t6.0_sol <- do.call(mcmc.list, t6.0_sol)

t6.0_vcv <- lapply(total_model6.0, function(m) m$VCV)
t6.0_vcv <- do.call(mcmc.list, t6.0_vcv)
t6.0_vcv <- mcmc.stack(t6.0_vcv)
plot(t6.0_vcv) 
autocorr.diag(t6.0_vcv, relative=FALSE) 
heidel.diag(t6.0_vcv)
effectiveSize(t6.0_vcv)  

# Model 6.1 - Informative Priors
total_leaf %>% group_by(treatment) %>%
  summarise(sd=sd(leaf), var=(sd)^2, mean1=mean(leaf))

load(file="Routput/Full/total_model6.1.RData")
t6.1_sol <- lapply(total_model6.1, function(m) m$Sol)
t6.1_sol <- do.call(mcmc.list, t6.1_sol)

t6.1_vcv <- lapply(total_model6.1, function(m) m$VCV)
t6.1_vcv <- do.call(mcmc.list, t6.1_vcv)
t6.1_vcv <- mcmc.stack(t6.1_vcv)
plot(t6.1_vcv) 
autocorr.diag(t6.1_vcv, relative=FALSE) 
heidel.diag(t6.1_vcv)
effectiveSize(t6.1_vcv) 

# Model 6.2 - Cross-environment Vd and Vm set to 0 (idh)
load(file="Routput/Full/total_model6.2.RData")
t6.2_sol <- lapply(total_model6.2, function(m) m$Sol)
t6.2_sol <- do.call(mcmc.list, t6.2_sol)

t6.2_vcv <- lapply(total_model6.2, function(m) m$VCV)
t6.2_vcv <- do.call(mcmc.list, t6.2_vcv)
t6.2_vcv <- mcmc.stack(t6.2_vcv)
plot(t6.2_vcv) 
autocorr.diag(t6.2_vcv, relative=FALSE) 
heidel.diag(t6.2_vcv)
effectiveSize(t6.2_vcv) 











## Trait 7: Height (Gaussian) ####
ggplot(total_height, aes(x=height)) + geom_histogram() + facet_grid(~treatment) + theme_bw()
total_height <- MCMC.2019 %>% filter(germ==1, !height==0)
total_prior7.0 <- list(R = list(V = diag(2)*0.02, nu = 3),
                       G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)
                       ))

set.seed(1)
total_model7.0 <- mclapply(1:5, function(i) {
  MCMCglmm(height ~ plot:treatment + treatment, 
           random = ~us(treatment):animal + us(treatment):animalDom + us(treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~idh(treatment):units,
           family = "gaussian", data = total_height, prior = total_prior7.0, #Gaussian distribution
           nitt = 500000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
}, mc.cores = 5)

#Model 7.0 - Full model
load(file="Routput/Full/total_model7.0.RData")
t7.0_sol <- lapply(total_model7.0, function(m) m$Sol)
t7.0_sol <- do.call(mcmc.list, t7.0_sol)

t7.0_vcv <- lapply(total_model7.0, function(m) m$VCV)
t7.0_vcv <- do.call(mcmc.list, t7.0_vcv)
t7.0_vcv <- mcmc.stack(t7.0_vcv)
plot(t7.0_vcv) 
autocorr.diag(t7.0_vcv, relative=FALSE) 
heidel.diag(t7.0_vcv)
effectiveSize(t7.0_vcv) 

# Model 7.1 - Informative priors
total_height%>% group_by(treatment) %>%
  summarise(sd=sd(height), var=(sd)^2) #var = 120 max *******

load(file="Routput/Full/total_model7.1.RData")
t7.1_sol <- lapply(total_model7.1, function(m) m$Sol)
t7.1_sol <- do.call(mcmc.list, t7.1_sol)

t7.1_vcv <- lapply(total_model7.1, function(m) m$VCV)
t7.1_vcv <- do.call(mcmc.list, t7.1_vcv)
t7.1_vcv <- mcmc.stack(t7.1_vcv)
plot(t7.1_vcv) 
autocorr.diag(t7.1_vcv, relative=FALSE) 
heidel.diag(t7.1_vcv)
effectiveSize(t7.1_vcv) 

# Model 7.2 - Cross-environment Vd and Vm set to 0 (idh)
load(file="Routput/Full/total_model7.2.RData")
t7.2_sol <- lapply(total_model7.2, function(m) m$Sol)
t7.2_sol <- do.call(mcmc.list, t7.2_sol)

t7.2_vcv <- lapply(total_model7.2, function(m) m$VCV)
t7.2_vcv <- do.call(mcmc.list, t7.2_vcv)
t7.2_vcv <- mcmc.stack(t7.2_vcv)
plot(t7.2_vcv) 
autocorr.diag(t7.2_vcv, relative=FALSE) 
heidel.diag(t7.2_vcv)
effectiveSize(t7.2_vcv)  











## Trait 8: Stem Diameter (Gaussian) ####
total_stemdiam <- MCMC.2019 %>% filter(germ==1, !stem_diam==0)
total_prior8.0 <- list(R = list(V = diag(2)*0.02, nu = 3),
                       G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)
                       ))

set.seed(1)
total_model8.0 <- mclapply(1:5, function(i) {
  MCMCglmm(flwr_clstr ~ plot:treatment + treatment, 
           random = ~us(treatment):animal + us(treatment):animalDom + us(treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~idh(treatment):units,
           family = "gaussian", data = total_stemdiam, prior = total_prior8.0, #Gaussian distribution
           nitt = 500000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
}, mc.cores = 5)

#Model 8.0 - Full model
load(file="Routput/Full/total_model8.0.RData") #total_model8.0.2 has longer nitt and nu = 2
t8.0_sol <- lapply(total_model8.0, function(m) m$Sol)
t8.0_sol <- do.call(mcmc.list, t8.0_sol)

t8.0_vcv <- lapply(total_model8.0, function(m) m$VCV)
t8.0_vcv <- do.call(mcmc.list, t8.0_vcv)
t8.0_vcv <- mcmc.stack(t8.0_vcv)
plot(t8.0_vcv) 
autocorr.plot(t8.0_vcv)
autocorr.diag(t8.0_vcv, relative=FALSE) 
heidel.diag(t8.0_vcv)
effectiveSize(t8.0_vcv) 

# Model 8.1 - Informative Priors
total_stemdiam%>% group_by(treatment) %>%
  summarise(sd=sd(stem_diam), var=(sd)^2) #var = 3.5 max *******

load(file="Routput/Full/total_model8.1.RData")
t8.1_sol <- lapply(total_model8.1, function(m) m$Sol)
t8.1_sol <- do.call(mcmc.list, t8.1_sol)

t8.1_vcv <- lapply(total_model8.1, function(m) m$VCV)
t8.1_vcv <- do.call(mcmc.list, t8.1_vcv)
t8.1_vcv <- mcmc.stack(t8.1_vcv)
plot(t8.1_vcv) 
autocorr.plot(t8.1_vcv)
autocorr.diag(t8.1_vcv, relative=FALSE) 
heidel.diag(t8.1_vcv)
effectiveSize(t8.1_vcv)  

# Model 8.2 - Cross-environment Vd and Vm set to 0 (idh)
load(file="Routput/Full/total_model8.2.RData")
t8.2_sol <- lapply(total_model8.2, function(m) m$Sol)
t8.2_sol <- do.call(mcmc.list, t8.2_sol)

t8.2_vcv <- lapply(total_model8.2, function(m) m$VCV)
t8.2_vcv <- do.call(mcmc.list, t8.2_vcv)
t8.2_vcv <- mcmc.stack(t8.2_vcv)
plot(t8.2_vcv) 
autocorr.plot(t8.2_vcv)
autocorr.diag(t8.2_vcv, relative=FALSE) 
heidel.diag(t8.2_vcv)
effectiveSize(t8.2_vcv) 











## Trait 9: Flowering Clusters Number (Poisson) ####
total_flwrclstr <- MCMC.2019 %>% filter(germ==1)
total_prior9.0 <- list(R = list(V = diag(2)*0.02, nu = 3),
                       G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                                G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)
                       ))

set.seed(1)
total_model9.0 <- mclapply(1:5, function(i) {
  MCMCglmm(flwr_clstr ~ plot:treatment + treatment,
           random = ~us(treatment):animal + us(treatment):animalDom + us(treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~idh(treatment):units,
           family = "poisson", data = total_flwrclstr, prior = total_prior9.0, #Poisson distribution
           nitt = 500000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)
}, mc.cores = 5)

#Model 9.0 - Full model
load(file="Routput/Full/total_model9.0.RData") #total_model9.0.2 has nu = 2 and it's very similar. 
t9.0_sol <- lapply(total_model9.0, function(m) m$Sol)
t9.0_sol <- do.call(mcmc.list, t9.0_sol)

t9.0_vcv <- lapply(total_model9.0, function(m) m$VCV)
t9.0_vcv <- do.call(mcmc.list, t9.0_vcv)
t9.0_vcv <- mcmc.stack(t9.0_vcv)
plot(t9.0_vcv) 
autocorr.diag(t9.0_vcv, relative=FALSE) 
heidel.diag(t9.0_vcv)
effectiveSize(t9.0_vcv) 

# Model 9.1 - Informative Prior
total_flwrclstr %>% group_by(treatment) %>%
  summarise(sd=sd(flwr_clstr), var=(sd)^2, mean1=mean(flwr_clstr),
            logmu = log(mean1)) #log(mu) is ~1
#Had alpha.V=diag(2)*2 and all V = diag(2)

load(file="Routput/Full/total_model9.1.RData")
t9.1_sol <- lapply(total_model9.1, function(m) m$Sol)
t9.1_sol <- do.call(mcmc.list, t9.1_sol)

t9.1_vcv <- lapply(total_model9.1, function(m) m$VCV)
t9.1_vcv <- do.call(mcmc.list, t9.1_vcv)
t9.1_vcv <- mcmc.stack(t9.1_vcv)
plot(t9.1_vcv) 
autocorr.diag(t9.1_vcv, relative=FALSE) 
heidel.diag(t9.1_vcv)
effectiveSize(t9.1_vcv) 


# Model 9.2 - Cross-environment Vd and Vm set to 0 (idh)
load(file="Routput/Full/total_model9.2.RData")
t9.2_sol <- lapply(total_model9.2, function(m) m$Sol)
t9.2_sol <- do.call(mcmc.list, t9.2_sol)

t9.2_vcv <- lapply(total_model9.2, function(m) m$VCV)
t9.2_vcv <- do.call(mcmc.list, t9.2_vcv)
t9.2_vcv <- mcmc.stack(t9.2_vcv)
plot(t9.2_vcv) 
autocorr.diag(t9.2_vcv, relative=FALSE) 
heidel.diag(t9.2_vcv)
effectiveSize(t9.2_vcv) 










## Trait 10: Lifetime Fitness (multivariate) (Gaussian) ####
total_lifetime <- MCMC.2019
total_prior10.0 <- list(R = list(V = diag(4)*0.002, nu = 5),
                        G = list(G1 = list(V = diag(4)*0.002, nu = 5, alpha.mu = c(0,0,0,0), alpha.V = diag(4)*1000),
                                 G2 = list(V = diag(4)*0.002, nu = 5, alpha.mu = c(0,0,0,0), alpha.V = diag(4)*1000),
                                 G3 = list(V = diag(4)*0.002, nu = 5, alpha.mu = c(0,0,0,0), alpha.V = diag(4)*1000)
                        ))

set.seed(1)
total_model10.0 <- mclapply(1:20, function(i) { 
  MCMCglmm(cbind(flower, seed_pods) ~ trait + plot:treatment + treatment - 1, 
           random = ~us(trait:treatment):animal + us(trait:treatment):animalDom + us(trait:treatment):matID, 
           ginverse = list(animal = Ainv, animalDom = Dinv), rcov=~idh(trait:treatment):units,
           family = c("threshold", "poisson"), data = total_lifetime, prior = total_prior10.0, #Bernoulli + Poisson
           nitt = 200000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)
}, mc.cores = 20)

#Model 10.0 - Full model
load(file="Routput/Full/total_model10.0.RData")
t10.0_sol <- lapply(total_model10.0, function(m) m$Sol)
t10.0_sol <- do.call(mcmc.list, t10.0_sol)

t10.0_vcv <- lapply(total_model10.0, function(m) m$VCV)
t10.0_vcv <- do.call(mcmc.list, t10.0_vcv)
t10.0_vcv <- mcmc.stack(t10.0_vcv)
plot(t10.0_vcv) 
autocorr.diag(t10.0_vcv, relative=FALSE) 
heidel.diag(t10.0_vcv)
effectiveSize(t10.0_vcv) 

# Model 10.1 - Informative Priors
load(file="Routput/Full/total_model9.1.RData")
t9.1_sol <- lapply(total_model9.1, function(m) m$Sol)
t9.1_sol <- do.call(mcmc.list, t9.1_sol)

t9.1_vcv <- lapply(total_model9.1, function(m) m$VCV)
t9.1_vcv <- do.call(mcmc.list, t9.1_vcv)
t9.1_vcv <- mcmc.stack(t9.1_vcv)
plot(t9.1_vcv) 
autocorr.plot(t9.1_vcv)
autocorr.diag(t9.1_vcv, relative=FALSE) 
heidel.diag(t9.1_vcv)
effectiveSize(t9.1_vcv)  

# Model 10.2 - Cross-environment Vd  and Vm set to 0 (idh)
load(file="Routput/Full/total_model9.2.RData")
t9.2_sol <- lapply(total_model9.2, function(m) m$Sol)
t9.2_sol <- do.call(mcmc.list, t9.2_sol)

t9.2_vcv <- lapply(total_model9.2, function(m) m$VCV)
t9.2_vcv <- do.call(mcmc.list, t9.2_vcv)
t9.2_vcv <- mcmc.stack(t9.2_vcv)
plot(t9.2_vcv) 
autocorr.plot(t9.2_vcv)
autocorr.diag(t9.2_vcv, relative=FALSE) 
heidel.diag(t9.2_vcv)
effectiveSize(t9.2_vcv) 


#'##################################################################################'#
###############      FULL MODELS EXCLUDING DOMINANCE ON TRAITS      ################
#'##################################################################################'#

#We previously ran full models including dominance but these models failed to converge.
#In an earlier analysis, we also ran full models *excluding* dominance. To see if we can estimate genetic correlations from separate environment models (e.g from scripts 03 and 04), we test if genetic correlations agree from full models *excluding* dominance and from separate environment models. 
#This script just includes the full models excluding dominance.


#Covariance in residuals is not identifiable since only one observation has one treatment applied
#Especially in binary models, the residual variance can't be identified.. so fix = 1 in R prior

## Trait 1: Survival Success (Binary) ####
load(file="Routput/Full_Excl_Dom/Plasticity_1.1.RData")

plastic.survival <- MCMC.2019
prior1.1 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25))),
                          G2 = list(V = diag(2), nu = 15, alpha.mu = c(0,0), alpha.V = diag(c(1.25, 1.25)))))

PL_model1.1 <- MCMCglmm(flower ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
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
load(file="Routput/Full_Excl_Dom/Plasticity_1.2.RData")

plastic.survival <- MCMC.2019
prior1.2 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model1.2 <- MCMCglmm(flower ~ plot:treatment + treatment, random = ~idh(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov = ~units,
                        family = "threshold", data = plastic.survival, prior = prior1.2, #Bernoulli distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

plot(PL_model1.2$VCV) # good trace plot
summary(PL_model1.2) # good effective sample size
autocorr.diag(PL_model1.2$VCV) # no autocorrelation
heidel.diag(PL_model1.2$VCV) # convergence success

#For nu = 3
posterior.mode(PL_model1.2$VCV[,'treatmentA:treatmentH.animal'])
HPDinterval(PL_model1.2$VCV[,'treatmentA:treatmentH.animal'])

gen.corrPL1.2 <-PL_model1.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model1.2$VCV[,'treatmentA:treatmentA.animal']*PL_model1.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL1.2) #Post Mean = 0.7852267
posterior.mode(gen.corrPL1.2) #Post Mode = 0.8276458
HPDinterval(gen.corrPL1.2) # Posterior 95% CI = (0.5547362, 0.979963)

#Random = idh() structure | Model 1.3 idh() | Model 1.4  but including maternal effects 
#Model 1.5 has R prior has V = diag(2) so rcov = ~idh(treatment):units .. nu = 3
#Model 1.6 has G prior with nu = 1000

plastic.survival <- MCMC.2019
prior1.6 <- list(R = list(V = diag(2), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1, 1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model1.6 <- MCMCglmm(flower ~ plot:treatment + treatment, random = ~idh(treatment):animal + idh(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov = ~idh(treatment):units,
                        family = "threshold", data = plastic.survival, prior = prior1.6, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Full_Excl_Dom/Plasticity_1.5.RData")
plot(PL_model1.5$VCV) # good trace plots
summary(PL_model1.5) # good effective sample size
autocorr.diag(PL_model1.5$VCV) # no autocorrelation
heidel.diag(PL_model1.5$VCV) #  convergence success

load(file="Routput/Full_Excl_Dom/Plasticity_1.6.RData")
plot(PL_model1.6$VCV) # good trace plots
summary(PL_model1.6) # good effective sample size
autocorr.diag(PL_model1.6$VCV) # no autocorrelation
heidel.diag(PL_model1.6$VCV) #  convergence success


#Pearson correlation between breeding values in AM and HW
PL_A1.3_BV <- (as.mcmc(PL_model1.3$Sol[,148:283])) #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
PL_H1.3_BV <- (as.mcmc(PL_model1.3$Sol[,7799:7934])) #extracting breeding values of parents - parents are from 136 to 271 on pedigree file
r.PL1.3 <- 1:2000 
for (i in 1:2000) {
  r.PL1.3[i] <- (cov(PL_A1.3_BV[i,], PL_H1.3_BV[i,])) /
    sqrt(var(PL_A1.3_BV[i,]) * var(PL_H1.3_BV[i,]))
}

posterior.mode(as.mcmc(r.PL1.3), na.rm=T) # point estimate of r (note any NAs, so that denominator can be adjusted below)
HPDinterval(as.mcmc(r.PL1.3)) # 95% HPD intervals of r. If these overlap zero then the correlation is not different from 0.
P.r.PL1.3 <- ifelse(r.PL1.3 <0, 0, 1) # Assuming that the point estimate is positive, we call those estimates that are <0, 0, else 1...
(1-sum(P.r.PL1.3)/2000)*2 # ...and this kind of a lazy way to get the P value, i.e. the proportion of those 1000 estimates that were below 0, times 2. 











## Trait 2: Fecundity (Poisson) ####
load(file="Routput/Full_Excl_Dom/Plasticity_2.1.RData")

plastic.fecundity <- MCMC.2019 %>% filter(germ==1 & flower ==1)
prior2.1 <- list(R = list(V = diag(1), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model2.1 <- MCMCglmm(seed_pods ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.fecundity, prior = prior2.1, #Poisson distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(PL_model2.1$VCV) # good trace plots
summary(PL_model2.1) # good effective sample size
autocorr.diag(PL_model2.1$VCV) # no autocorrelation
heidel.diag(PL_model2.1$VCV) # convergence fail for matID

## Trait 2 without maternal effects
load(file="Routput/Full_Excl_Dom/Plasticity_2.2.RData")

plastic.fecundity <- MCMC.2019 %>% filter(germ==1 & flower ==1)
prior2.2 <- list(R = list(V = diag(1), nu=0.002),
                 G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model2.2 <- MCMCglmm(seed_pods ~ plot:treatment + treatment, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.fecundity, prior = prior2.2, #Poisson distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE)

save(prior2.2, PL_model2.2, file="Routput/Plasticity-Rd6/Plasticity_2.2.RData")

plot(PL_model2.2$VCV) # good trace plots
summary(PL_model2.2) # good effective sample size
autocorr.diag(PL_model2.2$VCV) # no autocorrelation
autocorr.plot(PL_model2.2$VCV) #
heidel.diag(PL_model2.2$VCV) # convergence success


#With nu = 3
posterior.mode(PL_model2.2$VCV[,'treatmentA:treatmentH.animal']) #0.18965
HPDinterval(PL_model2.2$VCV[,'treatmentA:treatmentH.animal'])  # (0.09010, 0.3293)

gen.corrPL2.2 <-PL_model2.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model2.2$VCV[,'treatmentA:treatmentA.animal']*PL_model2.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL2.2) #Post Mean = 0.8728189
posterior.mode(gen.corrPL2.2) #Post Mode = 0.9252819
HPDinterval(gen.corrPL2.2) # Posterior 95% CI = (0.6895932, 0.9945704)

#Without trunc = TRUE, nu = 3
posterior.mode(PL_model2.2$VCV[,'treatmentA:treatmentH.animal']) # 0.176764
HPDinterval(PL_model2.2$VCV[,'treatmentA:treatmentH.animal']) # (0.08407403, 0.3242065)

gen.corrPL2.2 <-PL_model2.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model2.2$VCV[,'treatmentA:treatmentA.animal']*PL_model2.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL2.2) #Post Mean = 0.8668348
posterior.mode(gen.corrPL2.2) #Post Mode = 0.9225568
HPDinterval(gen.corrPL2.2) # Posterior 95% CI = (0.682938, 0.9950502)

#Random = idh() structure | Model 2.3 nu =3 | Model 2.4 nu=0.002 | Model 2.5 nu=0.002 but including maternal effects
#Model 2.6 has R prior has V = diag(2) so rcov = ~idh(treatment):units 
plastic.fecundity <- MCMC.2019 %>% filter(germ==1 & flower ==1)
prior2.6 <- list(R = list(V = diag(2), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))
                          ))

PL_model2.6 <- MCMCglmm(seed_pods ~ plot:treatment + plot, random = ~idh(treatment):animal + idh(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~idh(treatment):units,
                        family = "poisson", data = plastic.fecundity, prior = prior2.6, #Poisson distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Full_Excl_Dom/Plasticity_2.5.RData")
plot(PL_model2.5$VCV) # 
summary(PL_model2.5) # 
autocorr.diag(PL_model2.5$VCV) # 
heidel.diag(PL_model2.5$VCV) # 

load(file="Routput/Full_Excl_Dom/Plasticity_2.6.RData")
plot(PL_model2.6$VCV) # good trace plots
summary(PL_model2.6) # OK effective sample sizes
autocorr.diag(PL_model2.6$VCV) # moderate autocorrelation in Va
heidel.diag(PL_model2.6$VCV) # convergence success













## Trait 3: Germination Success (Binary) ####
load(file="Routput/Full_Excl_Dom/Plasticity_3.1.RData")

plastic.Germination <- MCMC.2019
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
load(file="Routput/Full_Excl_Dom/Plasticity_3.2.RData")

plastic.Germination <- MCMC.2019
prior3.2 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model3.2 <- MCMCglmm(germ ~ plot:treatment + treatment, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "threshold", data = plastic.Germination, prior = prior3.2, #Bernoulli distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(prior3.2, PL_model3.2, file="Routput/Plasticity-Rd6/Plasticity_3.2.RData")

plot(PL_model3.2$VCV) # good trace plots
summary(PL_model3.2) # good effective sample size
autocorr.diag(PL_model3.2$VCV) # no autocorrelation
autocorr.plot(PL_model3.2$VCV) #
heidel.diag(PL_model3.2$VCV) # convergence success

#With nu = 3
posterior.mode(PL_model3.2$VCV[,'treatmentA:treatmentH.animal'])
HPDinterval(PL_model3.2$VCV[,'treatmentA:treatmentH.animal'])

gen.corrPL3.2 <-PL_model3.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model3.2$VCV[,'treatmentA:treatmentA.animal']*PL_model3.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL3.2) #Post Mean = 0.708856
posterior.mode(gen.corrPL3.2) #Post Mode = 0.7674795
HPDinterval(gen.corrPL3.2) # Posterior 95% CI = (0.4017943, 0.9509564)

#Random = idh() structure | Model 3.3 idh() | Model 3.4  but including maternal effects
#Model 3.5 has R prior has V = diag(2) so rcov = ~idh(treatment):units 
#Model 3.6 has G prior with nu = 1000
plastic.Germination <- MCMC.2019
prior3.6 <- list(R = list(V = diag(2), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1, 1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model3.6 <- MCMCglmm(germ ~ plot:treatment + treatment, random = ~idh(treatment):animal + idh(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~idh(treatment):units,
                        family = "threshold", data = plastic.Germination, prior = prior3.6, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Full_Excl_Dom/Plasticity_3.5.RData")
plot(PL_model3.5$VCV) #
summary(PL_model3.5) # 
autocorr.diag(PL_model3.5$VCV) # 
heidel.diag(PL_model3.5$VCV) # 

load(file="Routput/Full_Excl_Dom/Plasticity_3.6.RData")
plot(PL_model3.6$VCV) # good trace plots
summary(PL_model3.6) # good effective sample size
autocorr.diag(PL_model3.6$VCV) # no autocorrelation
heidel.diag(PL_model3.6$VCV) # convergence success








## Trait 4: Flowering Success (Binary) ####
load(file="Routput/Full_Excl_Dom/Plasticity_4.1.RData")

plastic.Flowering <- MCMC.2019 %>% filter(germ==1)
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

## Trait 4 without maternal effects
load(file="Routput/Full_Excl_Dom/Plasticity_4.2.RData")

plastic.Flowering <- MCMC.2019 %>% filter(germ==1)
prior4.2 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model4.2 <- MCMCglmm(flower ~ plot:treatment + treatment, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "threshold", data = plastic.Flowering, prior = prior4.2, #Bernoulli distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(prior4.2, PL_model4.2, file="Routput/Plasticity-Rd6/Plasticity_4.2.RData")

plot(PL_model4.2$VCV) # OK trace plots
summary(PL_model4.2) # good effective sample size
autocorr.diag(PL_model4.2$VCV) # no autocorrelation
autocorr.plot(PL_model4.2$VCV) 
heidel.diag(PL_model4.2$VCV) # convergence success


#With nu = 3
posterior.mode(PL_model4.2$VCV[,'treatmentA:treatmentH.animal'])
HPDinterval(PL_model4.2$VCV[,'treatmentA:treatmentH.animal'])

gen.corrPL4.2 <-PL_model4.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model4.2$VCV[,'treatmentA:treatmentA.animal']*PL_model4.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL4.2) #Post Mean = 0.2077978
posterior.mode(gen.corrPL4.2) #Post Mode = 0.4965493
HPDinterval(gen.corrPL4.2) # Posterior 95% CI = (-0.5327902, 0.902755)

#Random = idh() structure | Model 4.3 idh() | Model 4.4  but including maternal effects
#Model 4.5 has R prior has V = diag(2) so rcov = ~idh(treatment):units 
#Model 4.6 has G prior with nu = 1000
plastic.Flowering <- MCMC.2019 %>% filter(germ==1)
prior4.6 <- list(R = list(V = diag(2), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1, 1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model4.6 <- MCMCglmm(flower ~ plot:treatment + treatment, random = ~idh(treatment):animal + idh(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~idh(treatment):units,
                        family = "threshold", data = plastic.Flowering, prior = prior4.6, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Full_Excl_Dom/Plasticity_4.5.RData")
plot(PL_model4.5$VCV) # 
summary(PL_model4.5) # 
autocorr.diag(PL_model4.5$VCV) # 
heidel.diag(PL_model4.5$VCV) # 

load(file="Routput/Full_Excl_Dom/Plasticity_4.6.RData")
plot(PL_model4.6$VCV) # good trace plots
summary(PL_model4.6) # good effective sample size
autocorr.diag(PL_model4.6$VCV) # no autocorrelation
heidel.diag(PL_model4.6$VCV) # convergence success











## Trait 5: Seed Maturation Success (Binary) ####
load(file="Routput/Full_Excl_Dom/Plasticity_5.1.RData")

plastic.Seed <- MCMC.2019 %>% filter(germ==1 & flower==1)
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

## Trait 5 without maternal effects | No maternal effects
load(file="Routput/Full_Excl_Dom/Plasticity_5.2.RData")

plastic.Seed <- MCMC.2019 %>% filter(germ==1 & flower==1)
prior5.2 <- list(R = list(V = diag(1), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model5.2 <- MCMCglmm(seed ~ plot:treatment + treatment, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "threshold", data = plastic.Seed, prior = prior5.2, #Bernoulli distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

save(prior5.2, PL_model5.2, file="Routput/Plasticity-Rd6/Plasticity_5.2.RData")

plot(PL_model5.2$VCV) # OK trace plots
summary(PL_model5.2) # good effective sample size
autocorr.diag(PL_model5.2$VCV) # no autocorrelation
heidel.diag(PL_model5.2$VCV) # convergence fail .. but this is expected  because there is almost NO variance in fruiting success between treatments... so genetic correlation would be essentially 0


#With nu = 3
posterior.mode(PL_model5.2$VCV[,'treatmentA:treatmentH.animal'])
HPDinterval(PL_model5.2$VCV[,'treatmentA:treatmentH.animal'])

gen.corrPL5.2 <-PL_model5.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model5.2$VCV[,'treatmentA:treatmentA.animal']*PL_model5.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL5.2) #Post Mean = 0.01037868
posterior.mode(gen.corrPL5.2) #Post Mode = 0.3027879
HPDinterval(gen.corrPL5.2) # Posterior 95% CI = (-0.761717, 0.8693172)

herit_PL5.2_A <-PL_model5.2$VCV[,'treatmentA:treatmentA.animal']/
  (PL_model5.2$VCV[,'treatmentA:treatmentA.animal'] + 1) 
mean(herit_PL5.2_A)

#Random = idh() structure | Model 5.3 idh() | Model 5.4  but including maternal effects
#Model 5.5 has R prior has V = diag(2) so rcov = ~idh(treatment):units 
#Model 5.6 has G prior with nu = 1000
plastic.Seed <- MCMC.2019 %>% filter(germ==1 & flower==1)
prior5.6 <- list(R = list(V = diag(2), nu = 2, fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1, 1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model5.6 <- MCMCglmm(seed ~ plot:treatment + treatment, random = ~idh(treatment):animal + idh(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~idh(treatment):units,
                        family = "threshold", data = plastic.Seed, prior = prior5.6, #Bernoulli distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE, trunc = TRUE)

load(file="Routput/Full_Excl_Dom/Plasticity_5.5.RData")
plot(PL_model5.5$VCV) # 
summary(PL_model5.5) # 
autocorr.diag(PL_model5.5$VCV) # 
heidel.diag(PL_model5.5$VCV) # 

load(file="Routput/Full_Excl_Dom/Plasticity_5.6.RData")
plot(PL_model5.6$VCV) # convergence success
summary(PL_model5.6) # high effective sample size
autocorr.diag(PL_model5.6$VCV) # no autocorrelation
heidel.diag(PL_model5.6$VCV) # convergence success















## Trait 6: Leaf Number (Gaussian) ####
load(file="Routput/Full_Excl_Dom/Plasticity_6.1.RData")


plastic.leaf <- MCMC.2019 %>% filter(germ==1)
prior6.1 <- list(R = list(V = diag(1), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model6.1 <- MCMCglmm(leaf ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.leaf, prior = prior6.1, #Gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(PL_model6.1$VCV) # good trace plots
summary(PL_model6.1) # good effective sample size
autocorr.diag(PL_model6.1$VCV) # no autocorrelation
heidel.diag(PL_model6.1$VCV) # convergence fail

#Trait 6 without maternal effects
load(file="Routput/Full_Excl_Dom/Plasticity_6.2.RData")


#Note for R elements, V=diag(1), and fix=1 because only 1 random effect.. although in G elements, crossing treatment with the random effect
plastic.leaf <- MCMC.2019 %>% filter(germ==1)
prior6.2 <- list(R = list(V = diag(1), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model6.2 <- MCMCglmm(leaf ~ plot:treatment + treatment, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.leaf, prior = prior6.2, #Gaussian distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE)

save(prior6.2, PL_model6.2, file="Routput/Plasticity-Rd6/Plasticity_6.2.RData")

plot(PL_model6.2$VCV) # good trace plots
summary(PL_model6.2) # good effective sample size
autocorr.diag(PL_model6.2$VCV) # no autocorrelation
heidel.diag(PL_model6.2$VCV) # convergence success


#With nu = 3
posterior.mode(PL_model6.2$VCV[,'treatmentA:treatmentH.animal'])
HPDinterval(PL_model6.2$VCV[,'treatmentA:treatmentH.animal'])

gen.corrPL6.2 <-PL_model6.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model6.2$VCV[,'treatmentA:treatmentA.animal']*PL_model6.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL6.2) #Post Mean = 0.896586
posterior.mode(gen.corrPL6.2) #Post Mode = 0.9411732
HPDinterval(gen.corrPL6.2) # Posterior 95% CI = (0.7732273, 0.9980959)

#With nu = 3 and without trunc=TRUE
posterior.mode(PL_model6.2$VCV[,'treatmentA:treatmentH.animal']) # 1.060199
HPDinterval(PL_model6.2$VCV[,'treatmentA:treatmentH.animal']) # (0.5878667, 1.61622)

gen.corrPL6.2 <-PL_model6.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model6.2$VCV[,'treatmentA:treatmentA.animal']*PL_model6.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL6.2) #Post Mean = 0.8948076
posterior.mode(gen.corrPL6.2) #Post Mode = 0.9441056
HPDinterval(gen.corrPL6.2) # Posterior 95% CI = (0.768444, 0.99639)

#Random = idh() structure | Model 6.3 nu =3 | Model 6.4 nu=0.002 | Model 6.5 nu=0.002 but including maternal effects
#Model 6.6 has R prior has V = diag(2) so rcov = ~idh(treatment):units 
plastic.leaf <- MCMC.2019 %>% filter(germ==1 & flower==1)
prior6.6 <- list(R = list(V = diag(2), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1, 1))),
                          G2 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model6.6 <- MCMCglmm(leaf ~ plot:treatment + treatment, random = ~idh(treatment):animal + idh(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~idh(treatment):units, 
                        family = "gaussian", data = plastic.leaf, prior = prior6.6, #Gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Full_Excl_Dom/Plasticity_6.5.RData")
plot(PL_model6.5$VCV) # 
summary(PL_model6.5) # 
autocorr.diag(PL_model6.5$VCV) # 
heidel.diag(PL_model6.5$VCV) # 

load(file="Routput/Full_Excl_Dom/Plasticity_6.6.RData")
plot(PL_model6.6$VCV) # good trace plots
summary(PL_model6.6) # good effective sample size
autocorr.diag(PL_model6.6$VCV) # no autocorrelation
heidel.diag(PL_model6.6$VCV) # convergence success










## Trait 7: Height (Gaussian) ####
load(file="Routput/Full_Excl_Dom/Plasticity_7.1.RData")

plastic.height <- MCMC.2019 %>% filter(germ==1, !height==0)
prior7.1 <- list(R = list(V = diag(1), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model7.1 <- MCMCglmm(height ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.height, prior = prior7.1, #Gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)


plot(PL_model7.1$VCV) # good trace plots.. except some bad
summary(PL_model7.1) # good effective sample size
autocorr.diag(PL_model7.1$VCV) # no autocorrelation
heidel.diag(PL_model7.1$VCV) # convergence fail


#Trait 7 without maternal effects
load(file="Routput/Full_Excl_Dom/Plasticity_7.2.RData")

plastic.height <- MCMC.2019 %>% filter(germ==1, !height==0)
prior7.2 <- list(R = list(V = diag(1), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model7.2 <- MCMCglmm(height ~ plot:treatment + treatment, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.height, prior = prior7.2, #Gaussian distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE)

save(prior7.2, PL_model7.2, file="Routput/Plasticity-Rd6/Plasticity_7.2.RData")

plot(PL_model7.2$VCV) # good trace plots
summary(PL_model7.2) # good effective sample size
autocorr.diag(PL_model7.2$VCV) # no autocorrelation
heidel.diag(PL_model7.2$VCV) # convergence success


#With nu = 3
posterior.mode(PL_model7.2$VCV[,'treatmentA:treatmentH.animal'])
HPDinterval(PL_model7.2$VCV[,'treatmentA:treatmentH.animal'])

gen.corrPL7.2 <-PL_model7.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model7.2$VCV[,'treatmentA:treatmentA.animal']*PL_model7.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL7.2) #Post Mean = 0.970785
posterior.mode(gen.corrPL7.2) #Post Mode = 0.9867147
HPDinterval(gen.corrPL7.2) # Posterior 95% CI = (0.9221566, 0.9995654)

#With nu = 3 and without trunc=TRUE
posterior.mode(PL_model7.2$VCV[,'treatmentA:treatmentH.animal']) # 27.33361
HPDinterval(PL_model7.2$VCV[,'treatmentA:treatmentH.animal']) # (17.86437, 45.06279)

gen.corrPL7.2 <-PL_model7.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model7.2$VCV[,'treatmentA:treatmentA.animal']*PL_model7.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL7.2) #Post Mean = 0.973086
posterior.mode(gen.corrPL7.2) #Post Mode = 0.9874767
HPDinterval(gen.corrPL7.2) # Posterior 95% CI = (0.9265667, 0.999596)


#Random idh() structure  | Model 7.3 nu =3 | Model 7.4 nu=0.002 | Model 7.5 nu=0.002 but including maternal effects
#Model 7.6 has R prior has V = diag(2) so rcov = ~idh(treatment):units 
plastic.height <- MCMC.2019 %>% filter(germ==1, !height==0)
prior7.6 <- list(R = list(V = diag(2), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1, 1))),
                          G2 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model7.6 <- MCMCglmm(height ~ plot:treatment + treatment, random = ~idh(treatment):animal + idh(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~idh(treatment):units,
                        family = "gaussian", data = plastic.height, prior = prior7.6, #Gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Full_Excl_Dom/Plasticity_7.5.RData")
plot(PL_model7.5$VCV) # 
summary(PL_model7.5) # 
autocorr.diag(PL_model7.5$VCV) # 
heidel.diag(PL_model7.5$VCV) # 


load(file="Routput/Full_Excl_Dom/Plasticity_7.6.RData")
plot(PL_model7.6$VCV) # good trace plots
summary(PL_model7.6) # good effective sample size
autocorr.diag(PL_model7.6$VCV) # no autocorrelation
heidel.diag(PL_model7.6$VCV) # convergence success (conservative warning on matID-A)







## Trait 8: Flowering Clusters Number (Gaussian) ####
load(file="Routput/Full_Excl_Dom/Plasticity_8.1.RData")

plastic.flwrclstr <- MCMC.2019 %>% filter(germ==1)
prior8.1 <- list(R = list(V = diag(1), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model8.1 <- MCMCglmm(flwr_clstr ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.flwrclstr, prior = prior8.1, #Poisson distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(PL_model8.1$VCV) # good trace plots
summary(PL_model8.1) # good effective sample size
autocorr.diag(PL_model8.1$VCV) # no autocorrelation
heidel.diag(PL_model8.1$VCV) # convergence fail

## Trait 8 without maternal effects
load(file="Routput/Full_Excl_Dom/Plasticity_8.2.RData")

plastic.flwrclstr <- MCMC.2019 %>% filter(germ==1)
prior8.2 <- list(R = list(V = diag(2), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model8.2 <- MCMCglmm(flwr_clstr ~ plot:treatment + treatment, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "poisson", data = plastic.flwrclstr, prior = prior8.2, #Poisson distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE)

save(prior8.2, PL_model8.2, file="Routput/Plasticity-Rd6/Plasticity_8.2.RData")

plot(PL_model8.2$VCV) # good trace plots
summary(PL_model8.2) # good effective sample size
autocorr.diag(PL_model8.2$VCV) # no autocorrelation
heidel.diag(PL_model8.2$VCV) # convergence success


#nu = 3
posterior.mode(PL_model8.2$VCV[,'treatmentA:treatmentH.animal'])
HPDinterval(PL_model8.2$VCV[,'treatmentA:treatmentH.animal'])

gen.corrPL8.2 <-PL_model8.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model8.2$VCV[,'treatmentA:treatmentA.animal']*PL_model8.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL8.2) #Post Mean = 0.846879
posterior.mode(gen.corrPL8.2) #Post Mode = 0.8945364
HPDinterval(gen.corrPL8.2) # Posterior 95% CI = (0.648065, 0.9963849)

#nu = 3 and without trunc=TRUE
posterior.mode(PL_model8.2$VCV[,'treatmentA:treatmentH.animal']) #0.1345883
HPDinterval(PL_model8.2$VCV[,'treatmentA:treatmentH.animal']) #(0.05413783, 0.2190489)

gen.corrPL8.2 <-PL_model8.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model8.2$VCV[,'treatmentA:treatmentA.animal']*PL_model8.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL8.2) #Post Mean = 0.8528467
posterior.mode(gen.corrPL8.2) #Post Mode = 0.9166738
HPDinterval(gen.corrPL8.2) # Posterior 95% CI = (0.653139, 0.9961332)

#Interestingly - what is the heritability within each treatment?

herit_PL8.2_A <-PL_model8.2$VCV[,'treatmentA:treatmentA.animal']/
  (PL_model8.2$VCV[,'treatmentA:treatmentA.animal'] + PL_model8.2$VCV[, 'units'])
mean(herit_PL8.2_A)

herit_PL8.2_H <-PL_model8.2$VCV[,'treatmentH:treatmentH.animal']/
  (PL_model8.2$VCV[,'treatmentH:treatmentH.animal'] + PL_model8.2$VCV[, 'units']) 
mean(herit_PL8.2_H)

#Random = idh() structure  | Model 8.3 nu =3 | Model 8.4 nu=0.002 | Model 8.5 nu=0.002 but including maternal effects
#Model 8.6 has R prior has V = diag(2) so rcov = ~idh(treatment):units 
#Model 8.6b (saved as 8.6b but named as 8.6) has nu = 0.002 in the G structure. Previously nu=3 for 8.6
plastic.flwrclstr <- MCMC.2019 %>% filter(germ==1)
prior8.6 <- list(R = list(V = diag(2), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1, 1))),
                          G2 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model8.6 <- MCMCglmm(flwr_clstr ~ plot:treatment + treatment, random = ~idh(treatment):animal + idh(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~idh(treatment):units,
                        family = "poisson", data = plastic.flwrclstr, prior = prior8.6, #Poisson distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Full_Excl_Dom/Plasticity_8.5.RData")
plot(PL_model8.5$VCV) # 
summary(PL_model8.5) # 
autocorr.diag(PL_model8.5$VCV) # 
heidel.diag(PL_model8.5$VCV) # 

load(file="Routput/Full_Excl_Dom/Plasticity_8.6.RData")
plot(PL_model8.6$VCV) # good trace plots
summary(PL_model8.6) # good effective sample size
autocorr.diag(PL_model8.6$VCV) # no autocorrelation
heidel.diag(PL_model8.6$VCV) # convergence success












## Trait 9: Stem Diameter (Gaussian) ####
load(file="Routput/Full_Excl_Dom/Plasticity_9.1.RData")

plastic.stemdiam <- MCMC.2019 %>% filter(germ==1, !stem_diam==0)
prior9.1 <- list(R = list(V = diag(1), fix = 1),
                 G = list(G1 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1))),
                          G2 = list(V = diag(2), nu = 1000, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model9.1 <- MCMCglmm(flwr_clstr ~ treatment + plot, random = ~us(treatment):animal + us(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.stemdiam, prior = prior9.1, #Gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

plot(PL_model9.1$VCV) # trace plots don't look great, although just concentrated
summary(PL_model9.1) # good effective sample size
autocorr.diag(PL_model9.1$VCV) # no autocorrelation
heidel.diag(PL_model9.1$VCV) # convergence fail

## Trait 9 without maternal effects
load(file="Routput/Full_Excl_Dom/Plasticity_9.2.RData")

plastic.stemdiam <- MCMC.2019 %>% filter(germ==1, !stem_diam==0)
prior9.2 <- list(R = list(V = diag(1), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 3, alpha.mu = c(0,0), alpha.V = diag(c(1,1)))))

PL_model9.2 <- MCMCglmm(flwr_clstr ~ plot:treatment + treatment, random = ~us(treatment):animal, 
                        ginverse = list(animal = Ainv), rcov=~units,
                        family = "gaussian", data = plastic.stemdiam, prior = prior9.2, #Gaussian distribution
                        nitt = 4100000, thin = 2000, burnin = 100000, verbose = T, pr = TRUE)

save(prior9.2, PL_model9.2, file="Routput/Plasticity-Rd6/Plasticity_9.2.RData")

plot(PL_model9.2$VCV) # good trace plots
summary(PL_model9.2) # good effective sample size
autocorr.diag(PL_model9.2$VCV) # low autocorrelation
heidel.diag(PL_model9.2$VCV) # convergence success


#For nu = 3
posterior.mode(PL_model9.2$VCV[,'treatmentA:treatmentH.animal'])
HPDinterval(PL_model9.2$VCV[,'treatmentA:treatmentH.animal'])

gen.corrPL9.2 <-PL_model9.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model9.2$VCV[,'treatmentA:treatmentA.animal']*PL_model9.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL9.2) #Post Mean = 0.9408211
posterior.mode(gen.corrPL9.2) #Post Mode = 0.9797814
HPDinterval(gen.corrPL9.2) # Posterior 95% CI = (0.8395064, 0.998984)

#For nu = 3
posterior.mode(PL_model9.2$VCV[,'treatmentA:treatmentH.animal']) #2.793273
HPDinterval(PL_model9.2$VCV[,'treatmentA:treatmentH.animal']) #(1.563941, 4.078929)

gen.corrPL9.2 <-PL_model9.2$VCV[,'treatmentA:treatmentH.animal']/
  sqrt(PL_model9.2$VCV[,'treatmentA:treatmentA.animal']*PL_model9.2$VCV[,'treatmentH:treatmentH.animal']) 
mean(gen.corrPL9.2) #Post Mean = 0.9434721
posterior.mode(gen.corrPL9.2) #Post Mode = 0.9740195
HPDinterval(gen.corrPL9.2) # Posterior 95% CI = (0.8513274, 0.9998809)

#Interestingly - what is the heritability within each treatment?

herit_PL9.2_A <-PL_model9.2$VCV[,'treatmentA:treatmentA.animal']/
  (PL_model9.2$VCV[,'treatmentA:treatmentA.animal'] + PL_model9.2$VCV[, 'units']) 
mean(herit_PL9.2_A)

herit_PL9.2_H <-PL_model9.2$VCV[,'treatmentH:treatmentH.animal']/
  (PL_model9.2$VCV[,'treatmentH:treatmentH.animal'] + PL_model9.2$VCV[, 'units']) 
mean(herit_PL9.2_H)


#Random = idh() structure  | Model 9.3 nu =3 | Model 9.4 nu=0.002 | Model 9.5 nu=0.002 but including maternal effects
#Model 9.6 has R prior has V = diag(2) so rcov = ~idh(treatment):units 
#Model 9.6b (saved as 9.6b but named as 9.6) has nu = 0.002 in the G structure. Previously nu=3 for 9.6
plastic.stemdiam <- MCMC.2019 %>% filter(germ==1, !stem_diam==0)
prior9.6 <- list(R = list(V = diag(2), nu = 0.002),
                 G = list(G1 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1, 1))),
                          G2 = list(V = diag(2), nu = 0.002, alpha.mu = c(0,0), alpha.V = diag(c(1, 1)))))

PL_model9.6 <- MCMCglmm(flwr_clstr ~ plot:treatment + treatment, random = ~idh(treatment):animal + idh(treatment):matID, 
                        ginverse = list(animal = Ainv), rcov=~idh(treatment):units,
                        family = "gaussian", data = plastic.stemdiam, prior = prior9.6, #Gaussian distribution
                        nitt = 2100000, thin = 1000, burnin = 100000, verbose = T, pr = TRUE)

load(file="Routput/Full_Excl_Dom/Plasticity_9.5.RData")
plot(PL_model9.5$VCV) # 
summary(PL_model9.5) # 
autocorr.diag(PL_model9.5$VCV) # 
heidel.diag(PL_model9.5$VCV) # 

load(file="Routput/Full_Excl_Dom/Plasticity_9.6.RData")
plot(PL_model9.6$VCV) # good trace plots
summary(PL_model9.6) # good effective sample size (decreased but ok)
autocorr.diag(PL_model9.6$VCV) # slight autocorrelation for matID-A, but still adequate
heidel.diag(PL_model9.6$VCV) # 


