#### PROJECT: Brassica rapa GxE Study (Data collected at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Flowering Time Analysis

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Installing and Loading Necessary Packages####
library(tidyverse)

#Loading cleaned dataset from 01_Data_exploration
load(file="Routput/heatarrays_data_explore.RData")
data1 <- dat

lapply(data1, class)

#'####################################################################'#
  ###############      DATASET PREPARATION & NOTES     ###############
#'####################################################################'#

#DATASET 1: data1 is the cleaned genernal dataset which includes:
      #germ_census, flwr_census, and flwrcompl_census contains the dates in census number format
      #germ_juln_date, flwr_juln_date, and flwr_compl_juln are the dates in julian format
#DATASET 2: data1a is for calculating cumulative germination, flowering and completed flowering over the season.
      #It includes all plants in the study and is a summarized version of data1
#DATASET 3: data1b is similar to 1a except it ONLY includes germinated plants
#DATASET 4: data_phen is for evaluating germination and flowering over time (NOT cumulative)



#Removing individuals that have seed date but no flwr date, etc
nrow(subset(data1, seed==1&flower==0)) # 0 errors
data1 <- data1[!(data1$seed==1&data1$flower==0),]

nrow(subset(data1, flwr_census==0&flwr_compl_juln>0)) # 2 errors
data1 <- data1[!(data1$flwr_juln_date==0&data1$flwr_compl_juln>0),]

#### Preparation for Flowering Time Analysis ####

#Selecting columns to clean dataset
data1 <- data1 %>%
  select(posID, individual, plot, matID, famID, treatment, germ, flower, seed,
         germ_juln_date, flwr_juln_date, flwr_compl_juln)#,
         #germ_census, flwr_census, flwrcompl_census)

#To create phenological values relative to germination time, add the following:  
  
  #mutate(flwr_census_relative = (as.numeric(flwr_census) - as.numeric(germ_census)),
        # flwrcompl_census_relative = (as.numeric(flwrcompl_census) - as.numeric(germ_census)))
#nrow(subset(data1, flwr_census_relative<0)) # 261 errors
#data1$flwr_census_relative[data1$flwr_census == 0] <- 0
#data1$flwrcompl_census_relative[data1$flwrcompl_census == 0] <- 0

#Creating Long-format Data for Graphing Purposes
data1a <- data1 %>%
  gather(key="variable", value="census", germ_juln_date:flwr_compl_juln)

summary.data <- data1a %>%
  group_by(treatment, census) %>%
  summarise(first.flwr=sum(variable=="flwr_juln_date"), #interchange with flwr_juln_relative or flwr_juln_date
            last.flwr=sum(variable=="flwr_compl_juln"), #interchange with 
            germ.sum=sum(variable=="germ_juln_date"))

data1 %>% #data1 = ALL plants
  group_by(treatment) %>%
  summarise(n=n())

data1b <- data1 %>% filter(germ==1)
data1b %>% #data1b = ONLY GERMINATING PLANTS
  group_by(treatment) %>%
  filter(germ==1) %>%
  summarise(n=n())

save(data1, file="Routput/Flowering_analysis_all_plants.RData")
save(data1a, file= "Routput/Flowering_analysis_all_plants_long.RData")
save(data1b, file="Routput/Flowering_analysis_germ_plants.RData")

load(file="Routput/Flowering_analysis_all_plants.RData")
load(file="Routput/Flowering_analysis_germ_plants.RData")

#Note, there are some germinated plants that DID NOT germinate.
#These plants are represented by census = 0

#'####################################################################'#
############       CUMULATIVE PHENOLOGICAL FIGURES       #############
#'####################################################################'#


#### Cumulative Percent Plants Germinating ####


#Values retrieved from the sum totals by Julian Calendar Day from summary.data
#Total values obtained from 01_Data_exploration tallies
germ.cum <- data.frame(treatment = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",
                                     "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"),
                       day = c(119,129,134,141,149,156,165,172,177,185,191,204,
                                  119,129,134,141,149,156,165,172,177,185,191,204),
                       value = c(602,607,209,110,40,12,6,1,1,0,0,0,
                                 1007,432,98,24,6,3,1,0,0,0,0,0),
                       total = c(3458, 3458, 3458, 3458, 3458, 3458,
                                 3458, 3458, 3458, 3458, 3458, 3458,
                                 3346, 3346, 3346, 3346, 3346, 3346,
                                 3346, 3346, 3346, 3346, 3346, 3346))

germ.cum <- germ.cum %>% group_by(treatment) %>% 
  mutate(cum.sum=cumsum(value), germ_cum_perc=100*(cum.sum/total))
lapply(germ.cum, class)

germ_cum_fig <- germ.cum %>%
  ggplot(aes(x=day, y=germ_cum_perc, group=treatment, color=treatment)) +
  geom_point(size=2, alpha = 0.75) +
  geom_line(size=1.2, alpha = 0.5) + 
  #scale_x_continuous(breaks=seq(1,12,1)) + #ONLY IF CENSUS
  scale_y_continuous(limits=c(0,100)) +
  labs(x="Julian calendar day", y="Cumulative % of plants germinated", tag = "A)")  +
  scale_color_manual("Environment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.text = element_text(size=10), legend.title = element_text(size=10))
germ_cum_fig
ggsave(plot = germ_cum_fig, filename = "Routput/Figures/fig_s5a.png", width = 8, height = 8, units ="cm")




#### Cumulative Percent Plants Reaching Flowering ####

#Values retrieved from the sum totals by Julian Calendar Day  from summary.data
#Total values obtained from 01_Data_exploration tallies
first.flwr <- data.frame(treatment = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",
                                   "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"),
                     day = c(119,129,134,141,149,156,165,172,177,185,191,204,
                                119,129,134,141,149,156,165,172,177,185,191,204),
                     value = c(0,0,0,0,0,3,284,539,372,180,35,8,
                               0,0,0,0,8,287,772,284,84,32,8,1),
                     total = c(1589, 1589, 1589, 1589, 1589, 1589, 
                               1589, 1589, 1589, 1589, 1589, 1589,
                               1576, 1576, 1576, 1576, 1576, 1576,
                               1576, 1576, 1576, 1576, 1576, 1576))

first.flwr <- first.flwr %>% group_by(treatment) %>% 
  mutate(cum.sum=cumsum(value), flwr_cum_perc=100*(cum.sum/total))
lapply(first.flwr, class)

first_flwr_cum_fig <- first.flwr %>%
  ggplot(aes(x=day, y=flwr_cum_perc, group=treatment, color=treatment)) +
  geom_point(size=2, alpha = 0.75) +
  geom_line(size=1.2, alpha = 0.5) + 
  #scale_x_continuous(breaks=seq(1,8,1)) + #ONLY IF CENSUS
  scale_y_continuous(limits=c(0,100)) +
  labs(x="Julian calendar day", y="Cumulative % of plants reached flowering", tag = "B)")  +
  scale_color_manual("Environment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.text = element_text(size=10), legend.title = element_text(size=10))
first_flwr_cum_fig
ggsave(plot = first_flwr_cum_fig, filename = "Routput/Figures/fig_s5b.png", width = 8, height = 8, units ="cm")


#### Cumulative Percent Plants Finishing Flowering ####

#Values retrieved from the sum totals by Julian Calendar Day  from summary.data
#Total values obtained from 01_Data_exploration tallies
last.flwr <- data.frame(treatment = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",
                                      "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"),
                         day = c(119,129,134,141,149,156,165,172,177,185,191,204,
                                    119,129,134,141,149,156,165,172,177,185,191,204),
                         value = c(0,0,0,0,0,0,0,17,223,775,309,88,
                                   0,0,0,0,0,0,1,356,578,466,55,19),
                         total = c(1589, 1589, 1589, 1589, 1589, 1589,
                                   1589, 1589, 1589, 1589, 1589, 1589,
                                   1576, 1576, 1576, 1576, 1576, 1576,
                                   1576, 1576, 1576, 1576, 1576, 1576))

last.flwr <- last.flwr %>% group_by(treatment) %>% 
  mutate(cum.sum=cumsum(value), flwr_cum_perc=100*(cum.sum/total))
lapply(last.flwr, class)

last_flwr_cum_fig <- last.flwr %>%
  ggplot(aes(x=day, y=flwr_cum_perc, group=treatment, color=treatment)) +
  geom_point(size=2, alpha = 0.75) +
  geom_line(size=1.2, alpha = 0.5) + 
  #scale_x_continuous(breaks=seq(1,8,1)) + #ONLY IF CENSUS
  scale_y_continuous(limits=c(0,100)) +
  labs(x="Julian calendar day", y="Cumulative % of plants completed flowering", tag = "C)")  +
  scale_color_manual("Environment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.text = element_text(size=10), legend.title = element_text(size=10))
last_flwr_cum_fig
ggsave(plot = last_flwr_cum_fig, filename = "Routput/Figures/fig_s5c.png", width = 8, height = 8, units ="cm")

#### Cumulative Percent Plants Reaching Flowering - Date Relative to Germination ####

#Values retrieved from the sum totals by Julian Calendar Day from summary.data
#Total values obtained from 01_Data_exploration tallies
first.flwr.relative <- data.frame(treatment = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",
                                       "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"),
                         day = c(119,129,134,141,149,156,165,172,177,185,191,204,
                                    119,129,134,141,149,156,165,172,177,185,191,204),
                         value = c(0,3,4,24,153,541,487,173,30,5,1,0,
                                   0,1,2,71,509,642,180,53,18,0,0,0),
                         total = c(1589, 1589, 1589, 1589, 1589, 1589, 
                                   1589, 1589, 1589, 1589, 1589, 1589,
                                   1576, 1576, 1576, 1576, 1576, 1576,
                                   1576, 1576, 1576, 1576, 1576, 1576))



first.flwr.relative <- first.flwr.relative %>% group_by(treatment) %>% 
  mutate(cum.sum=cumsum(value), flwr_cum_perc=100*(cum.sum/total))
lapply(first.flwr.relative, class)

first.flwr.relative %>%
  ggplot(aes(x=day, y=flwr_cum_perc, group=treatment, color=treatment)) +
  geom_point(size=4) +
  geom_line(size=1.2) + 
  #scale_x_continuous(breaks=seq(1,12,1)) + #ONLY IF CENSUS
  scale_y_continuous(limits=c(0,100)) +
  ggtitle("A)") +
  labs(x="Julian Calendar Day", y="Cumulative % of Plants Reached Flowering")  +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#### Cumulative Percent Plants Finishing Flowering - Date Relative to Germination ####

#Values retrieved from the sum totals by Julian Calendar Day from summary.data
#Total values obtained from 01_Data_exploration tallies
last.flwr.relative <- data.frame(treatment = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",
                                      "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"),
                        day = c(119,129,134,141,149,156,165,172,177,185,191,204,
                                   119,129,134,141,149,156,165,172,177,185,191,204),
                        value = c(0,0,2,4,17,61,208,512,453,131,23,0,
                                  0,0,1,0,12,121,463,514,310,40,10,0),
                        total = c(1589, 1589, 1589, 1589, 1589, 1589,
                                  1589, 1589, 1589, 1589, 1589, 1589,
                                  1576, 1576, 1576, 1576, 1576, 1576,
                                  1576, 1576, 1576, 1576, 1576, 1576))



last.flwr.relative <- last.flwr.relative %>% group_by(treatment) %>% 
  mutate(cum.sum=cumsum(value), flwr_cum_perc=100*(cum.sum/total))
lapply(last.flwr.relative, class)

last.flwr.relative %>%
  ggplot(aes(x=day, y=flwr_cum_perc, group=treatment, color=treatment)) +
  geom_point(size=4) +
  geom_line(size=1.2) + 
  #scale_x_continuous(breaks=seq(1,12,1)) + #ONLY IF CENSUS
  scale_y_continuous(limits=c(0,100)) +
  ggtitle("B)") +
  labs(x="Julian Calendar Day", y="Cumulative % of Plants Completed Flowering")  +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))









#'####################################################################'#
  ##################       PHENOLOGY FIGURES      ####################
#'####################################################################'#

#Using tidyverse to convert dataframe into long format
data_phen<- dat %>%
  select(posID, individual, plot, matID, patID, famID, mat_group, treatment,
         germ, flower,
         germ_census, flwr_census, flwrcompl_census, 
         leaf1, leaf2, leaf3, leaf4, leaf5, leaf6, leaf7, leaf8, leaf9,
         flwr_clstr1, bud_clstr1, flwr_clstr2, bud_clstr2, flwr_clstr3, bud_clstr3,
         flwr_clstr4, bud_clstr4, flwr_clstr5, bud_clstr5, flwr_clstr6, bud_clstr6,
         flwr_clstr7, bud_clstr7, flwr_clstr8, bud_clstr8, seed_pods, seed_dmg, seed_undmg) %>%
  gather("leaf1", "leaf2", "leaf3", "leaf4", "leaf5", "leaf6", "leaf7", "leaf8", "leaf9", "flwr_clstr1", "flwr_clstr2", "flwr_clstr3", "flwr_clstr4", "flwr_clstr5", "flwr_clstr6", "flwr_clstr7", "flwr_clstr8",
         "bud_clstr1", "bud_clstr2", "bud_clstr3", "bud_clstr4", "bud_clstr5", "bud_clstr6", "bud_clstr7", "bud_clstr8",
         key = "phenotype", value = "phen_num")  #This can be used to see the change across time for various phenotypes.

#Converting leaf number to Census Date
data_phen$census <- data_phen$phenotype #ignore the error, just some error b/c tibble doesn't know the class type
data_phen$census[data_phen$census=="leaf1"] <- 3
data_phen$census[data_phen$census=="leaf2"] <- 4
data_phen$census[data_phen$census=="leaf3"] <- 5
data_phen$census[data_phen$census=="leaf4"] <- 6
data_phen$census[data_phen$census=="leaf5"] <- 7
data_phen$census[data_phen$census=="leaf6"] <- 8
data_phen$census[data_phen$census=="leaf7"] <- 9
data_phen$census[data_phen$census=="leaf8"] <- 10
data_phen$census[data_phen$census=="leaf9"] <- 11

#Converting flower cluster number to census date
data_phen$census[data_phen$census=="flwr_clstr1"] <- 5
data_phen$census[data_phen$census=="flwr_clstr2"] <- 6
data_phen$census[data_phen$census=="flwr_clstr3"] <- 7
data_phen$census[data_phen$census=="flwr_clstr4"] <- 8
data_phen$census[data_phen$census=="flwr_clstr5"] <- 9
data_phen$census[data_phen$census=="flwr_clstr6"] <- 10
data_phen$census[data_phen$census=="flwr_clstr7"] <- 11
data_phen$census[data_phen$census=="flwr_clstr8"] <- 12

#Converting bud cluster number to census date
data_phen$census[data_phen$census=="bud_clstr1"] <- 5
data_phen$census[data_phen$census=="bud_clstr2"] <- 6
data_phen$census[data_phen$census=="bud_clstr3"] <- 7
data_phen$census[data_phen$census=="bud_clstr4"] <- 8
data_phen$census[data_phen$census=="bud_clstr5"] <- 9
data_phen$census[data_phen$census=="bud_clstr6"] <- 10
data_phen$census[data_phen$census=="bud_clstr7"] <- 11
data_phen$census[data_phen$census=="bud_clstr8"] <- 12

save(data_phen, file="Routput/heatarrays_data_phen.RData")

#### Graphs 21 - Flower Cluster No. vs Time ####
#Note: dates are converted in census # as data was collected weekly
#Subtracting flwr_clstr from bud_clstr = flower cluster number per census since data was collected as a cumulative number 
load(file="Routput/heatarrays_data_phen.RData")

data_phen$census <- as.factor(data_phen$census)
data_phen %>%
  filter(flower==1, str_detect(phenotype, "flwr_clstr")) %>%
  group_by(treatment, census) %>%
  summarise(flwr_avg=sum(phen_num, na.rm=TRUE), n=n(), sd=sd(flwr_avg), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = census, y = flwr_avg, group = treatment, color = treatment)) + 
  #geom_errorbar(aes(ymin=pod_avg-se, ymax=pod_avg+se), width=0.5, size=0.5) + 
  geom_smooth(se=F) +
  ggtitle("Flowering Clusters Observed in Summer 2019") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Census Week", y="Flower Cluster Number")

#### Graphs 22 - Growth Rate on Leaf No. vs Time ####
data_phen$census <- as.factor(data_phen$census)
data_phen %>%
  filter(germ==1, str_detect(phenotype, "leaf")) %>%
  group_by(treatment, census) %>%
  summarise(leaf_avg=sum(phen_num, na.rm=TRUE), n=n(), sd=sd(leaf_avg), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = census, y = leaf_avg, group = treatment, color = treatment)) + 
  #geom_errorbar(aes(ymin=pod_avg-se, ymax=pod_avg+se), width=0.5, size=0.5) + 
  geom_smooth(se=F) +
  ggtitle("Leaf Production Observed in Summer 2019") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Census Week", y="True Leaf Number")














#'####################################################################'#
 ##################       SURVIVAL ANALYSIS      #####################
#'####################################################################'#

## Packages required
library(survival)
library(survminer)
library(coxme)

## I employ untruncated, right-censored Cox models including random plot effects
## I use data1 which includes all flowering plants
load(file="Routput/Flowering_analysis_all_plants.RData")

#Setting Appropriate Classes to Numeric/Double as classes were altered during cleanup
lapply(data1, class)
data1$posID <- as.factor(data1$posID)
data1$individual <- as.integer(data1$individual)

#In this Survival Analysis, we give every individual a value which is the time to achieving some outcome 
  #(e.g first flower, finish flowering, or germinate).
#In the censor data, we give each individual that did NOT achieve some outcome (e.g y = 0) a censor value of 0, meaning it is "still alive' and has not progressed to the outcome

# SUMMARY MEDIAN PHENOLOGICAL DATES #

#Calculating median time to first flower, germination, and last flower treatment
data1 %>%
  filter(germ==1) %>%
  group_by(treatment) %>%
  summarise(first_flower = median(flwr_juln_date), germ = median(germ_juln_date), last_flower = median(flwr_compl_juln))



#### (NOT) Survival Analysis 1: Flowering Season Length ####

#Previous analysis is archived. See 01_Data_exploration for season length analysis as this is an incorrect approach. Requires a glmm. 
#Copy and pasted from 01_Data_exploration with some edits

library(lme4)
season.length <- data1 %>% 
  filter(flower==1, !flwr_compl_juln==0) %>%
  mutate(slength=(as.numeric(flwr_compl_juln) - as.numeric(flwr_juln_date)))

hist(season.length$slength[season.length$treatment=="H"], breaks=50)
qqnorm(season.length$slength[season.length$treatment=="H"])
qqline(season.length$slength[season.length$treatment=="H"]) #non-normal distribution
summary(model1a <- lmer(dat=dat1, first_flower ~ treatment + famID + (treatment|plot), REML=FALSE))
summary(model1b <- lmer(dat=dat1, first_flower ~ treatment * famID + (treatment|plot), REML=FALSE))


#






#### Survival Analysis 2: Flowering Time - Time to First Flower ####
surv.first <- data1 %>%
  filter(germ==1) %>%
  mutate(first_flower=(as.numeric(flwr_juln_date))) %>%
  mutate(censor= ifelse(first_flower <= 0, 0, 1)) %>%
  mutate(first_flower = ifelse(first_flower <= 0, 199, first_flower)) %>% #We give an arbitrary value of 199 to individuals with no first flower - this just indicates the last date to which we sampled the data
  mutate(treatment_B = ifelse(treatment == "H", 1, 0)) %>%
  select(individual, plot, treatment, matID, famID, treatment_B, first_flower, censor)

lapply(surv.first, class)
surv.first$famID <- as.factor(surv.first$famID) #
surv.first$plot <- as.factor(surv.first$plot) #
surv.first$treatment <- as.factor(surv.first$treatment)
hist(surv.first$first_flower[surv.first$treatment=="H"])
hist(surv.first$first_flower[surv.first$treatment=="A"])

#COX SURVIVAL MODELS
#Note, censor has 0 = did not flower (i.e still alive), 1 = flowered (i.e dead in Survival analysis). 
#Use coxph for models w/o a random effect, and coxme with random effects
cox.first_2A <- coxme(Surv(first_flower, censor, type="right") ~ treatment + famID + (1|treatment/plot), data=surv.first)
cox.first_2B <- coxme(Surv(first_flower, censor, type="right") ~ treatment * famID + (1|treatment/plot), data=surv.first)
cox.first_2C <- coxph(Surv(first_flower, censor, type="right") ~ treatment + famID + frailty(plot), data=surv.first)

#MODEL COMPARISON + RESULTS
anova(cox.first_2A, cox.first_2B) #Model 2B including interaction is better
summary(cox.first_2B) #Heated treatment DECREASES time to first flower by exp(2.0549) = 7.8 days
anova(cox.first_2B) #P value < 2.2e-16, Chisq = 431.87
survdiff(Surv(first_flower, censor, type="right") ~ treatment, data=surv.first) #Another method to check significance b/w treatments 

#NOTE: Because this is a SURVIVAL analysis
confint(cox.first_2B, level=0.95)
#Confidence Intervals: Lower 95% = exp(1.52756212) = 4.60, Upper 95% = exp(2.582433284) = 13.22

#Creating dataframe for plotting - not necessary though
cox.first.fit <- survfit(Surv(first_flower, censor, type="right") ~ treatment, data=surv.first)
ggsurvplot(cox.first.fit, data=surv.first)

#








#### Survival Analysis 3: Germination Time ####
surv.germ <- data1 %>%
  mutate(germ_time=as.numeric(germ_juln_date)) %>%
  mutate(censor=ifelse(germ_time <= 0, 0, 1)) %>%
  mutate(germ_time=ifelse(germ_time <= 0, 176, germ_time)) %>% #We give an arbitrary value of 176 to individuals with no germination time - this just indicates the last date to which we sampled the data
  mutate(treatment_B=ifelse(treatment == "H", 1, 0)) %>%
  select(individual, plot, treatment, matID, treatment_B, famID, germ_time, censor)

lapply(surv.germ, class)
surv.germ$famID <- as.factor(surv.germ$famID) #
surv.germ$plot <- as.factor(surv.germ$plot) #
surv.germ$treatment <- as.factor(surv.germ$treatment)
hist(surv.germ$germ_time[surv.germ$treatment=="H"])
hist(surv.germ$germ_time[surv.germ$treatment=="A"])

#COX SURVIVAL MODELS
#Note, censor has 0 = did not flower (i.e still alive), 1 = flowered (i.e dead in Survival analysis). 
cox.germ_3A <- coxme(Surv(germ_time, censor, type="right") ~ treatment + famID + (1|treatment/plot), data=surv.germ)
cox.germ_3B <- coxme(Surv(germ_time, censor, type="right") ~ treatment * famID + (1|treatment/plot), data=surv.germ)
cox.germ_3C <- coxph(Surv(germ_time, censor, type="right") ~ treatment + famID + frailty(plot), data=surv.germ)

#MODEL COMPARISON + RESULTS
#Same as 3, just using frailty to compare differences in coxph vs coxme and use for diagnostics
anova(cox.germ_3A, cox.germ_3B) # Model 3B is better
summary(cox.germ_3B) #Heated treatment DECREASES time to germination by exp(0.382295710) = 1.4656454 Julian days
#CAUTION: The p value from the summary is 0.2 but from the ANOVA, it is greater
anova(cox.germ_3B) # P value < 2.2e-16, Chisq = 190.898
survdiff(Surv(germ_time, censor, type="right") ~ treatment, data=surv.germ) #Another method to check significance b/w treatments 
confint(cox.germ_3B, level=0.95)
#Confidence Intervals: Lower 95% = exp(-0.200664902) = 0.818, Upper 95% = exp(0.96525632) = 2.62


#Creating dataframe for plotting - not necessary though
cox.germ.fit <- survfit(Surv(germ_time, censor, type="right") ~ treatment, data=surv.germ)
ggsurvplot(cox.germ.fit, data=surv.germ)
cox.germ.fit
#







#### Survival Analysis 4: Time to End of Flowering ####

#Unlike earlier graphs, I can make use of the flowering juln time without using census # as a dummy variable
#Treatment H = 1, A = 0 ... therefore comparison results are shown relative to the heated treatment
surv.last <- data1 %>% 
  filter(germ==1) %>%
  mutate(last_flower=as.numeric(flwr_compl_juln)) %>%
  mutate(censor= ifelse(last_flower <= 0, 0, 1)) %>%
  mutate(last_flower = ifelse(last_flower <= 0, 204, last_flower)) %>% #We give an arbitrary value of 36 to individuals with no season length - this just indicates the last date to which we sampled the data
  mutate(treatment_B = ifelse(treatment == "H", 1, 0)) %>%
  select(individual, plot, treatment, matID, famID, treatment_B, last_flower, censor)

lapply(surv.last, class)
surv.last$treatment <- as.factor(surv.last$treatment)
surv.last$famID <- as.factor(surv.last$famID) #Treating sibship is numeric so get generalized results
surv.last$plot <- as.factor(surv.last$plot) #same for plot effects
hist(surv.last$last_flower[surv.last$treatment=="H"])
hist(surv.last$last_flower[surv.last$treatment=="A"])

#COX SURVIVAL MODELS
#Note, censor has 0 = did not flower (i.e still alive), 1 = flowered (i.e dead in Survival analysis). 
#Use coxph for models w/o a random effect, and coxme with random effects
cox.time_4A <- coxme(Surv(last_flower, event=censor, type="right") ~ treatment + famID + (1|treatment/plot), data=surv.last)
cox.time_4B <- coxme(Surv(last_flower, event=censor, type="right") ~ treatment * famID + (1|treatment/plot), data=surv.last)
cox.time_4C <- coxph(Surv(last_flower, event=censor, type="right") ~ treatment + famID + frailty(plot), data=surv.last) #Same as 3, just using frailty to compare differences in coxph vs coxme and use for diagnostics

#RESULTS
anova(cox.time_4A, cox.time_4B) #Model 4B including interaction is better
summary(cox.time_4B) #Heated treatment DECREASES time to last flower by exp(1.076855830) = 2.935 Julian days
anova(cox.time_4B) #P value < 2.2e-16, Chisq = 380.08
survdiff(Surv(last_flower, censor, type="right") ~ treatment, data=surv.last) #Another method to check significance b/w treatments but does not take into account of family genetic effects
confint(cox.time_4B, level=0.95)
#Confidence Intervals: Lower 95% = exp(0.557251545) = 1.74, Upper 95% = exp(1.596460115) = 4.93

#Creating dataframe for plotting - not necessary though
cox.last.fit <- survfit(Surv(last_flower, censor, type="right") ~ treatment, data=surv.last)
ggsurvplot(cox.last.fit, data=surv.last)



#





















  





####OLD METHOD OF CALCULATING####
 
#Selecting columns needed for analysis (OLD)
data1 <- data1 %>%
  select(posID, individual, plot, matID, patID, famID, mat_group, treatment, 
       germ_census, flwr_census, flwrcompl_census, germ, flower, seed,
       flwr_clstr1, flwr_clstr2, flwr_clstr3, flwr_clstr4, flwr_clstr5, flwr_clstr6, flwr_clstr7, flwr_clstr8)
#### Percent Plants Reaching Flowering ####

#Totals: A=1589, H=1576
data1 %>% group_by(treatment) %>% 
  filter(!germ==0) %>%
  summarise(flower1 = sum(flwr_clstr1>0), n()) #A: 0 | H: 12

data1 %>% group_by(treatment) %>%
  filter(!germ==0, !flwr_clstr1>0) %>%
  summarise(flower2 = sum(flwr_clstr2>0), n()) #A: 4 | H: 291

data1 %>% group_by(treatment) %>%
  filter(!germ==0, !flwr_clstr1>0, !flwr_clstr2>0) %>%
  summarise(flower3 = sum(flwr_clstr3>0), n()) #A: 285 | H: 766

data1 %>% group_by(treatment) %>%
  filter(!germ==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0) %>%
  summarise(flower4 = sum(flwr_clstr4>0), n()) #A: 539 | H: 289

data1 %>% group_by(treatment) %>%
  filter(!germ==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0, !flwr_clstr4>0) %>%
  summarise(flower5 = sum(flwr_clstr5>0), n()) #A: 373 | H: 79

data1 %>% group_by(treatment) %>%
  filter(!germ==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0, !flwr_clstr4>0, !flwr_clstr5>0) %>%
  summarise(flower6 = sum(flwr_clstr6>0), n()) #A: 177 | H: 33

data1 %>% group_by(treatment) %>%
  filter(!germ==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0, !flwr_clstr4>0, !flwr_clstr5>0, !flwr_clstr6>0) %>%
  summarise(flower7 = sum(flwr_clstr7>0), n()) #A: 37 | H: 10

data1 %>% group_by(treatment) %>%
  filter(!germ==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0, !flwr_clstr4>0, !flwr_clstr5>0, !flwr_clstr6>0, !flwr_clstr7>0) %>%
  summarise(flower8 = sum(flwr_clstr8>0), n()) #A: 6 | H: 1

flower <- data.frame(treatment = c("A", "A", "A", "A", "A", "A", "A", "A",
                                   "H", "H", "H", "H", "H", "H", "H", "H"),
                     census = c(1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8),
                     value = c(0,4,285,539,373,177,37,6,
                               12,291,766,289,79,33,10,1),
                     total = c(1589, 1589, 1589, 1589, 1589, 1589, 1589, 1589,
                               1576, 1576, 1576, 1576, 1576, 1576, 1576, 1576))
flower <- flower %>% group_by(treatment) %>% 
  mutate(cum.sum=cumsum(value), flwr_cum_perc=100*(cum.sum/total))
lapply(flower, class)

flower %>%
  ggplot(aes(x=census, y=flwr_cum_perc, group=treatment, color=treatment)) +
  geom_point(size=4) +
  geom_line(size=1.2) + 
  scale_x_continuous(breaks=seq(1,8,1)) +
  scale_y_continuous(limits=c(0,100)) +
  labs(x="Census", y="Cumulative Percentage of Plants that Reached Flowered")  +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))





#### ARCHIVE: Survival Analysis 1: Flowering Season Length ####

## UPDATE + WARNING: This may be the incorrect model for flowering season length. It is probably better to use a GLMM

#Average season length by treatment
data1b %>%
  group_by(treatment) %>%
  summarise(mean_first_flwr=mean(flwr_juln_relative))

#Unlike earlier graphs, I can make use of the flowering juln time without using census # as a dummy variable
#Treatment H = 1, A = 0 ... therefore comparison results are shown relative to the heated treatment
surv.length <- data1 %>% 
  filter(germ==1) %>%
  mutate(season_length=(as.numeric(flwr_compl_juln) - as.numeric(flwr_juln_date))) %>%
  mutate(censor= ifelse(season_length <= 0, 0, 1)) %>%
  mutate(season_length = ifelse(season_length <= 0, 0, season_length)) %>% #We give an arbitrary value of 36 to individuals with no season length - this just indicates the last date to which we sampled the data
  mutate(treatment_B = ifelse(treatment == "H", 1, 0)) %>%
  select(individual, plot, treatment, matID, famID, treatment_B, season_length, censor)

lapply(surv.length, class)
surv.length$treatment <- as.factor(surv.length$treatment)
surv.length$famID <- as.factor(surv.length$famID) #Treating sibship is numeric so get generalized results
surv.length$plot <- as.factor(surv.length$plot) #same for plot effects
hist(surv.length$season_length[surv.length$treatment=="H"])
hist(surv.length$season_length[surv.length$treatment=="A"])

#COX SURVIVAL MODELS
#Note, censor has 0 = did not flower (i.e still alive), 1 = flowered (i.e dead in Survival analysis). 
#Use coxph for models w/o a random effect, and coxme with random effects
cox.time_1A <- coxme(Surv(season_length, event=censor, type="right") ~ treatment + famID + (1|treatment/plot), data=surv.length)
cox.time_1B <- coxme(Surv(season_length, event=censor, type="right") ~ treatment * famID + (1|treatment/plot), data=surv.length)
cox.time_1C <- coxph(Surv(season_length, event=censor, type="right") ~ treatment + famID + frailty(plot), data=surv.length) #Same as 3, just using frailty to compare differences in coxph vs coxme and use for diagnostics

#RESULTS
anova(cox.time_1A, cox.time_1B) #Model 1B including interaction is better
summary(cox.time_1B) #Heated treatment INCREASES flowering season length by exp(-0.094234881) = 0.910 (Julian days?)
anova(cox.time_1B) #P value = 0.0001717, Chisq = 14.118
survdiff(Surv(season_length, censor, type="right") ~ treatment, data=surv.length) #Another method to check significance b/w treatments but does not take into account of family genetic effects
confint(cox.time_1B, level=0.95)
#Confidence Intervals: Lower 95% = exp(-0.612346918) = 0.542, Upper 95% = exp(0.4238772) = 1.527

#Creating dataframe for plotting - not necessary though
cox.time.fit <- survfit(Surv(season_length, censor, type="right") ~ treatment, data=surv.length)
ggsurvplot(cox.time.fit, data=surv.length)
