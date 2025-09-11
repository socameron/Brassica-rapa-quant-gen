#### PROJECT: Brassica rapa Va/W Study (Data collected by Cameron So 2019 at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Clean data and observe distributions prior to MCMCglmm analysis
#### AUTHOR: Cameron So

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Loading Necessary Packages
library(gplots)
library(tidyverse)
library(zoo)
library(car)
library(multcompView)
#library(lsmeans) #lsmeans now integrated into emmeans
library(emmeans)
library(FSA)

#Importing Data using readr from tidyverse (Base R confuses the class for certain vectors). Excluding notes from the data import
# d for double, f for factor, i for integer
col_types_list <- cols_only(posID = "d", individual = "d", plot = "f", animal = "f", 
                            matID = "f", patID = "f", famID = "f",
                            mat_group = col_factor(levels=c("1", "2")), 
                            treatment = col_factor(levels=c("A", "H")),
                            double = "d", damage = "d", packID = "f",
                            germ_juln_date = "i", flwr_juln_date = "i", flwr_compl_juln = "i",
                            leaf1 = "d", leaf2 = "d", leaf3 = "d", leaf4 = "d", leaf5 = "d", 
                            leaf6 = "d", leaf7 = "d", leaf8 = "d", leaf9 = "d",
                            flwr_clstr1 = "d", flwr_clstr2 = "d",  flwr_clstr3 = "d", flwr_clstr4 = "d",
                            flwr_clstr5 = "d", flwr_clstr6 = "d",  flwr_clstr7 = "d", flwr_clstr8 = "d",
                            bud_clstr1 = "d", bud_clstr2 = "d", bud_clstr3 = "d", bud_clstr4 = "d",
                            bud_clstr5 = "d", bud_clstr6 = "d", bud_clstr7 = "d", bud_clstr8 = "d",
                            germ_bin = "d", flwr_bin = "d", seed_bin = "d",
                            seed_dmg = "d", seed_undmg = "d", seed_pods = "d", seed_number = "d",
                            height = "d", stem_diam = "d")
data <- read_csv("Rdata/heatarrays_brassica_2019_data.csv", col_names = TRUE, na = "NA",
                 col_types=col_types_list)
spec(data)
data <- as.data.frame(data)
lapply(data, class)
lapply(data, levels)

#Note: I am not using the calendar date data. However, if you wish to use that data, use package readr from tidyverse and function read_csv to parse the date columns correctly. 

#'####################################################################'#
######################      CLEANING DATASET      ######################
#'####################################################################'#

#Number of Positions wth Doubles
data %>% group_by(plot) %>% summarise(double=sum(double, na.rm=TRUE), n()) #229 doubles
  #Plot 7 has 79 doubles

#Number of accidentally damaged plants
data %>% summarise(damage=sum(damage, na.rm=TRUE), n()) #46 damaged

#From raw dataset, removing plant positions with double plants or accidental damage. 
data <- data %>%
  filter(is.na(double), is.na(damage)) 

#Double checking correct number of removals
data %>% summarise(n())#originally N=7270, now N=6995 (7270-229-46)
#7270 not 7260 but 10 planted additionally by accident (packet had low seed count)

#Number of individuals in families not to be included
data %>% summarise(exclude=sum(packID=="M1S20-14", packID=="M1S31-4", packID=="M1S33-13", packID=="M1S43-9", packID=="M2S28-18")) #186 individuals

#Removing families not to be included in the analysis
data <- data %>%
  filter(!packID=="M1S20-14", !packID=="M1S31-4", !packID=="M1S33-13", !packID=="M1S43-9", !packID=="M2S28-18")

#Double checking individuals have been removed
data %>% summarise(n()) #6809 individuals (6995-186)

#Converting NAs in data to 0
data[is.na(data)] <- 0

#Setting Appropriate Classes to Numeric/Double as classes were altered during cleanup
lapply(data, class)
data$posID <- as.factor(data$posID)
data$individual <- as.integer(data$individual)
data[,11:49] <- lapply(data[,11:49], as.numeric) #setting all neccessary parameters to class numeric

#Checking which flwr_clstr columns do not have corresponding data in the bud_clstr
nrow(subset(data, flwr_clstr1>1 & bud_clstr1==0))
nrow(subset(data, flwr_clstr2>1 & bud_clstr2==0))
nrow(subset(data, flwr_clstr3>1 & bud_clstr3==0))
nrow(subset(data, flwr_clstr4>1 & bud_clstr4==0))
nrow(subset(data, flwr_clstr5>1 & bud_clstr5==0))
nrow(subset(data, flwr_clstr6>1 & bud_clstr6==0))
nrow(subset(data, flwr_clstr7>1 & bud_clstr7==0))
nrow(subset(data, flwr_clstr8>1 & bud_clstr8==0))

#Alternatively, use a for loop
map2_int(data %>% select(starts_with("flwr_clstr")), 
         data %>% select(starts_with("bud_clstr")), 
         ~sum(.x  > 1 & .y == 0))  %>% unname() #No errors detected

#Grouping continuous traits into one parameter
leaf.col <- c("leaf1", "leaf2", "leaf3", "leaf4", "leaf5", "leaf6", "leaf7", "leaf8", "leaf9")
flwr.col <- c("flwr_clstr1", "flwr_clstr2", "flwr_clstr3", "flwr_clstr4", "flwr_clstr5", "flwr_clstr6", "flwr_clstr7", "flwr_clstr8")
bud.col <- c("bud_clstr1", "bud_clstr2", "bud_clstr3", "bud_clstr4", "bud_clstr5", "bud_clstr6", "bud_clstr7", "bud_clstr8")

#Collating continuous traits into one parameter (using the maximum value)
dat <- data %>%
  mutate(leaf=do.call(pmax, data[,leaf.col]),
         flwr_clstr=do.call(pmax, data[,flwr.col]),
         bud_clstr=do.call(pmax, data[,bud.col])) %>%
  mutate(germ=germ_bin, flower=flwr_bin, seed=seed_bin) %>% #renaming the fitness parameters for easier coding
  select(posID, animal, individual, plot, matID, patID, famID, mat_group, treatment,
         germ_juln_date, flwr_juln_date, flwr_compl_juln,
         leaf, flwr_clstr, bud_clstr, germ, flower, seed, seed_pods, seed_number, seed_dmg, seed_undmg, seed_number,
         height, stem_diam)

#Checking classes of the dataframe
lapply(dat, class)

#Looking for detection errors in data for binomial traits (germ, flower, seed)
nrow(subset(dat,germ==0&flower==0&seed==1)) # 0 errors
nrow(subset(dat,germ==0&flower==1&seed==0)) # 0 errors
nrow(subset(dat,germ==0&flower==1&seed==1)) # 0 errors


#Looking for inconsistencies in data for continuous traits (leaf, flwr_clstr, seed_pods)
nrow(subset(dat, leaf==0&flwr_clstr>0&seed_pods>0)) # 2 incidences
nrow(subset(dat, leaf==0&flwr_clstr==0&seed_pods>0)) # 3 incidences
nrow(subset(dat, flwr_clstr==0&seed_pods>0)) # 4 incidences
nrow(subset(dat, leaf==0&flwr_clstr>0)) # 2 incidences
nrow(subset(dat, leaf==0&seed_pods>0)) # 5 incidences
nrow(subset(dat, seed==1&flwr_clstr==0)) # 4 incidences 
nrow(subset(dat, seed_pods==0&seed_number>0)) #16 errors

#Evaluating individual errors
head(subset(dat, leaf==0&flwr_clstr>0&seed_pods>0)) # P7 372 & 392 are OK b/c no leaf data.
head(subset(dat, leaf==0&flwr_clstr==0&seed_pods>0)) #P2I2 needs to be checked in the collected bags
head(subset(dat, flwr_clstr==0&seed_pods>0)) #P3 441 and P7 157 to be deleted; #P7 170 & 329 no data 
head(subset(dat, leaf==0&flwr_clstr>0))
head(subset(dat, leaf==0&seed_pods>0))
head(subset(dat, seed_pods==0&seed_number>0)) # 16 errors - Check the lab data

#Fixing errors (removing 4 individuals)
dat <- dat[!(dat$flwr_clstr==0&dat$seed_pods>0),]
dat <- dat[!(dat$seed==1&dat$flwr_clstr==0),]
#When checking for errors again, only 2 errors because P7 372 & 392 remain in the dataset

#Looking for inconsistencies errors between binary and continuous traits (flwr_clstr & flower, and seed and seed_pods)
nrow(subset(dat, germ==0&flwr_clstr>0)) # 0 errors
nrow(subset(dat, germ==0&seed_pods>0)) # 0 errors
nrow(subset(dat, flower==0&flwr_clstr>0)) # 0 errors
nrow(subset(dat, flower==0&seed_pods>0)) # 0 errors
nrow(subset(dat, seed==0&seed_pods>0)) # 0 errors
nrow(subset(dat, seed==1&flwr_clstr==0)) # 0 errors
nrow(subset(dat, seed_pods==0&seed_number>0)) # 16 errors remain - but this will be removed in the second analysis for seed number only. 


#DATA PREP FOR GRAPHS 14-16 in 02_Fitness_figures ####
#Dates are converted into census number as data was collected weekly
#This provides some level of simplicity to the model

#Converting Germination Dates to Census Dates
dat$germ_census <- dat$germ_juln_date #ignore the error as a new column is created in this 'tibble'
dat$germ_census[dat$germ_census==119] <- 1
dat$germ_census[dat$germ_census==129] <- 2
dat$germ_census[dat$germ_census==134] <- 3
dat$germ_census[dat$germ_census==141] <- 4
dat$germ_census[dat$germ_census==147] <- 5
dat$germ_census[dat$germ_census==148] <- 5
dat$germ_census[dat$germ_census==149] <- 5
dat$germ_census[dat$germ_census==154] <- 6
dat$germ_census[dat$germ_census==155] <- 6
dat$germ_census[dat$germ_census==156] <- 6
dat$germ_census[dat$germ_census==162] <- 7
dat$germ_census[dat$germ_census==163] <- 7
dat$germ_census[dat$germ_census==165] <- 7
dat$germ_census[dat$germ_census==169] <- 8
dat$germ_census[dat$germ_census==170] <- 8
dat$germ_census[dat$germ_census==171] <- 8
dat$germ_census[dat$germ_census==172] <- 8
dat$germ_census[dat$germ_census==175] <- 9
dat$germ_census[dat$germ_census==176] <- 9
dat$germ_census[dat$germ_census==177] <- 9
dat$germ_census[dat$germ_census==183] <- 10
dat$germ_census[dat$germ_census==184] <- 10
dat$germ_census[dat$germ_census==185] <- 10
dat$germ_census[dat$germ_census==189] <- 11
dat$germ_census[dat$germ_census==190] <- 11
dat$germ_census[dat$germ_census==191] <- 11
dat$germ_census[dat$germ_census==199] <- 12
dat$germ_census[dat$germ_census==204] <- 12

#Converting First Flower Dates to Census #
dat$flwr_census <- dat$flwr_juln_date
dat$flwr_census[dat$flwr_census==119] <- 1
dat$flwr_census[dat$flwr_census==129] <- 2
dat$flwr_census[dat$flwr_census==134] <- 3
dat$flwr_census[dat$flwr_census==141] <- 4
dat$flwr_census[dat$flwr_census==147] <- 5
dat$flwr_census[dat$flwr_census==148] <- 5
dat$flwr_census[dat$flwr_census==149] <- 5
dat$flwr_census[dat$flwr_census==154] <- 6
dat$flwr_census[dat$flwr_census==155] <- 6
dat$flwr_census[dat$flwr_census==156] <- 6
dat$flwr_census[dat$flwr_census==162] <- 7
dat$flwr_census[dat$flwr_census==163] <- 7
dat$flwr_census[dat$flwr_census==165] <- 7
dat$flwr_census[dat$flwr_census==169] <- 8
dat$flwr_census[dat$flwr_census==170] <- 8
dat$flwr_census[dat$flwr_census==171] <- 8
dat$flwr_census[dat$flwr_census==172] <- 8
dat$flwr_census[dat$flwr_census==175] <- 9
dat$flwr_census[dat$flwr_census==176] <- 9
dat$flwr_census[dat$flwr_census==177] <- 9
dat$flwr_census[dat$flwr_census==183] <- 10
dat$flwr_census[dat$flwr_census==184] <- 10
dat$flwr_census[dat$flwr_census==185] <- 10
dat$flwr_census[dat$flwr_census==189] <- 11
dat$flwr_census[dat$flwr_census==190] <- 11
dat$flwr_census[dat$flwr_census==191] <- 11
dat$flwr_census[dat$flwr_census==199] <- 12
dat$flwr_census[dat$flwr_census==204] <- 12

dat$flwrcompl_census <- dat$flwr_juln_date
dat$flwrcompl_census[dat$flwrcompl_census == 119] <- 1
dat$flwrcompl_census[dat$flwrcompl_census == 129] <- 2
dat$flwrcompl_census[dat$flwrcompl_census == 134] <- 3
dat$flwrcompl_census[dat$flwrcompl_census == 141] <- 4
dat$flwrcompl_census[dat$flwrcompl_census == 147] <- 5
dat$flwrcompl_census[dat$flwrcompl_census == 148] <- 5
dat$flwrcompl_census[dat$flwrcompl_census == 149] <- 5
dat$flwrcompl_census[dat$flwrcompl_census == 154] <- 6
dat$flwrcompl_census[dat$flwrcompl_census == 155] <- 6
dat$flwrcompl_census[dat$flwrcompl_census == 156] <- 6
dat$flwrcompl_census[dat$flwrcompl_census == 162] <- 7
dat$flwrcompl_census[dat$flwrcompl_census == 163] <- 7
dat$flwrcompl_census[dat$flwrcompl_census == 165] <- 7
dat$flwrcompl_census[dat$flwrcompl_census == 169] <- 8
dat$flwrcompl_census[dat$flwrcompl_census == 170] <- 8
dat$flwrcompl_census[dat$flwrcompl_census == 171] <- 8
dat$flwrcompl_census[dat$flwrcompl_census == 172] <- 8
dat$flwrcompl_census[dat$flwrcompl_census == 175] <- 9
dat$flwrcompl_census[dat$flwrcompl_census == 176] <- 9
dat$flwrcompl_census[dat$flwrcompl_census == 177] <- 9
dat$flwrcompl_census[dat$flwrcompl_census == 183] <- 10
dat$flwrcompl_census[dat$flwrcompl_census == 184] <- 10
dat$flwrcompl_census[dat$flwrcompl_census == 185] <- 10
dat$flwrcompl_census[dat$flwrcompl_census == 189] <- 11
dat$flwrcompl_census[dat$flwrcompl_census == 190] <- 11
dat$flwrcompl_census[dat$flwrcompl_census == 191] <- 11
dat$flwrcompl_census[dat$flwrcompl_census == 199] <- 12
dat$flwrcompl_census[dat$flwrcompl_census == 204] <- 12

save(dat, file="Routput/heatarrays_data_explore.RData")

#'####################################################################'#
######################      DATA EXPLORATION      ######################
#'####################################################################'#

load(file="Routput/heatarrays_data_explore.RData")

dat %>%
  group_by(treatment) %>%
  summarise(n=n())

#Setting plot theme for ggplot
theme_set(theme_classic(base_size = 10))

#Tally and Proportions by Treatment for Binary Traits ####
#Consequent lifetime fitness traits are conditional on previous traits (e.g flowering success if germ=1)
#Also including z-test proportional tests

dat %>% #overwintering survival and germination
  group_by(treatment) %>%
  summarise(germ=sum(germ), n=n()) %>%
  mutate(germ_perc=(germ/n))
prop.test(x = c(1589, 1575), n=c(3457, 3348))

dat %>% #spring-summer survival
  group_by(treatment) %>%
  filter(germ==1) %>%
  summarise(flower=sum(flower), n=n()) %>%
  mutate(flwr_perc=(flower/n))
prop.test(x = c(1421, 1481), n=c(1589, 1575))

dat %>% #overall survival to flowering
  group_by(treatment) %>%
  summarise(surv=sum(flower), n=n()) %>%
  mutate(surv_perc=(surv/n), var=(surv_perc*n*(1-surv_perc)))
prop.test(x = c(1421, 1481), n=c(3456, 3346))

dat %>% #fruiting success
  group_by(treatment, famID) %>%
  filter(germ==1, flower==1) %>%
  summarise(seed=sum(seed), n=n()) %>%
  mutate(seed_perc=(seed/n))
prop.test(x = c(1411, 1475), n=c(1421, 1481))


#Tally and Averages by Treatment for Continuous Traits ####

a2 <- dat %>% #fecundity of flowering plants
  filter(!germ==0, !flower==0) %>%
  group_by(treatment) %>%
  summarise(pod_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n))), var=var(seed_pods))

a1 <- dat %>% #lifetime fitness
  group_by(treatment) %>%
  summarise(pod_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), var=var(seed_pods), se=(sd/(sqrt(n))), var=var(seed_pods))

dat %>% #flowering clusters
  filter(!germ==0) %>%
  group_by(treatment) %>%
  summarise(flower_avg=mean(flwr_clstr), n=n(), sd=sd(flwr_clstr), se=(sd/(sqrt(n))), var=var(flwr_clstr))

dat %>% #leaf number
  filter(!germ==0) %>%
  group_by(treatment) %>%
  summarise(leaf_avg=mean(leaf), n=n(), sd=sd(leaf), se=(sd/(sqrt(n))), var=var(leaf))

a3 <- dat %>% #height
  filter(!germ==0, !height==0) %>%
  group_by(treatment) %>%
  summarise(height_avg=mean(height), n=n(), sd=sd(height), se=(sd/(sqrt(n))), var=var(height))

dat %>% #stem diameter
  filter(!germ==0, !stem_diam==0) %>%
  group_by(treatment) %>%
  summarise(stem_avg=mean(stem_diam), n=n(), sd=sd(stem_diam), se=(sd/(sqrt(n))), var=var(stem_diam))


#Tally for Bird Damage ####
dat %>% 
  filter(!germ==0, !flower==0, !seed==0) %>%
  group_by(treatment, famID) %>%
  summarise(bird_plants=sum(seed_dmg>0), undmg_plants=sum(seed_undmg>=0 & seed_dmg==0), total_plants=sum(seed_pods>=0))

dat %>%
  filter(!germ==0, !flower==0, !seed==0) %>%
  group_by(treatment) %>%
  summarise(seeded_plants=sum(seed==1))

#Tally for Number of Maternal Sub-groups ####
dat %>% 
  group_by(mat_group, famID) %>%
  summarise(count = n_distinct(matID)) %>%
  summarise(mat_subgroups = sum(count)) #M1 has 39 mat-subgroups; M2 has 79 mat-subgroups

dat %>%
  group_by(mat_group) %>%
  summarise(sibships = n_distinct(famID)) #M1 has 22 sibships; M2 has 40 sibships


#Graphs 1 - Flowering Season Length ####
dat %>% #graphing by treatment only
  filter(!flwr_juln_date==0, !flwr_compl_juln==0, !flower==0, !germ==0) %>%
  mutate(flowering_time=(flwr_compl_juln - flwr_juln_date)) %>%
  filter(!flowering_time==0) %>% #Removed individuals with flowering time = 0 because we did not collect the flowering dates correctly in the field
  group_by(treatment) %>%
  ggplot(aes(x=treatment, y=flowering_time)) +
  geom_boxplot(aes(fill=treatment), lwd=1) +
  labs(x="Environment", y="Flowering Season Length (Days)") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0,40)) +
  scale_fill_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient, Heated")) + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm")) +
  coord_flip() +
  geom_jitter(aes(color=treatment), alpha=0.2, shape=16, position=position_jitter(0.2)) +
  scale_color_manual(values=c("dodgerblue4", "tomato4")) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=6,fill="darkgreen",alpha=0.5)
  
  #Flowering Time Summary
dat %>% #graphing by treatment only
  filter(!flwr_juln_date==0, !flwr_compl_juln==0, !flower==0, !germ==0) %>%
  mutate(flowering_time=(flwr_compl_juln - flwr_juln_date)) %>%
  filter(!flowering_time==0) %>%
  group_by(treatment) %>%
  summarise(flowering_time=median(flowering_time))

  #This put spacing between the xy labs and the graph
  #xlab("\nTreatment") + ylab("Flowering Season Length (Days)\n")
  
dat %>% #graphing by treatment and famID ... it's really messy! (62 different boxplots)
    filter(!flwr_juln_date==0, !flwr_compl_juln==0, !flower==0, !germ==0) %>%
    mutate(flowering_time=(flwr_compl_juln - flwr_juln_date)) %>%
    filter(!flowering_time==0) %>% #Removed individuals with flowering time = 0 because we did not collect the flowering dates correctly in the field
    filter(mat_group==2) %>% #maternal group 2 only
    group_by(treatment, famID) %>%
    ggplot(aes(x=famID, y=flowering_time), group=treatment) +
    geom_boxplot(aes(fill=treatment), lwd=1) +
    labs(x="Environment", y="Flowering season length (Days)") +
    scale_x_discrete(labels=c("Ambient", "Heated")) + 
    scale_fill_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient, Heated")) + 
    theme_bw() + theme(legend.position="none") +
    coord_flip()

### Flowering Date 
dat %>% #time to first flower after germination
  filter(!flwr_juln_date==0, !flwr_compl_juln==0, !flower==0, !germ==0) %>%
  mutate(first_flower=(flwr_juln_date - germ_juln_date)) %>%
  group_by(treatment) %>%
  ggplot(aes(x=treatment, y=first_flower)) +
  geom_boxplot(aes(fill=treatment), lwd=1) +
  labs(x="Environment", y="Days to first flower") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 80, 10)) +
  scale_fill_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient, Heated")) + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm")) +
  geom_jitter(aes(color=treatment), alpha=0.2, shape=16, position=position_jitter(0.2)) +
  scale_color_manual(values=c("dodgerblue4", "tomato4")) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=6,fill="darkgreen",alpha=0.5)  +
  coord_flip()


#Graphs 2 - Histograms for Continuous Traits (INCLUDING germ = 0) ####
#Including ENTIRE Dataset (where germ = 0)
dat %>% #flowering clusters
  ggplot(aes(x=flwr_clstr, fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of flowering clusters", y="Frequency")

dat %>% #lifetime fitness
  ggplot(aes(x=seed_pods, fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of seed pods", y="Frequency")

dat %>% #leaf number
  ggplot(aes(x=leaf, fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of leaves", y="Frequency") 

dat %>% #height
  ggplot(aes(x=height, fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Height (cm)", y="Frequency")

dat %>% #stem diameter
  ggplot(aes(x=log(stem_diam), fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Stem diameter (mm)", y="Frequency")

#Graphs 3 - Density plots for Continuous Traits (EXCLUDING germ = 0) ####
#Including SUBSET Dataset (EXCLUDING where germ = 0)

mu.flwr <- dat %>%
  filter(!germ==0) %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(flwr_clstr))
dat %>% #flowering clusters
  filter(!germ==0) %>%
  ggplot(aes(x=flwr_clstr, fill=treatment, color=treatment))+
  geom_density(adjust=2, alpha=0.3) + #NOTE, it would be multi-modal if adjust=1
  geom_vline(data=mu.flwr, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_x_continuous(breaks=seq(0, 35, 5)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Number of flowering clusters", y="Frequency") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

mu.pods <- dat %>%
  filter(!germ==0, !flower==0) %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(seed_pods))
dat %>% #fecundity of flowering plants
  filter(!germ==0, !flower==0) %>%
  ggplot(aes(x=seed_pods, fill=treatment, color=treatment)) +
  geom_density(alpha = 0.3) +
  geom_vline(data=mu.pods, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_x_continuous(breaks=seq(0, 300, 50)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Number of seed pods from flowering plants", y="Frequency")  + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

mu.leaf <- dat %>%
  filter(!germ==0) %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(leaf))
dat %>% #leaf number
  filter(!germ==0) %>%
  ggplot(aes(x=leaf, fill=treatment, color=treatment)) +
  geom_density(alpha=0.3) +
  geom_vline(data=mu.leaf, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_x_continuous(breaks=seq(0, 20, 2)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Number of leaves", y="Frequency")  + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

mu.height <- dat %>%
  filter(!germ==0, !height==0) %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(height))
dat %>% #height
  filter(!germ==0, !height==0) %>% #exlcuding height==0 because they were previously NA values
  ggplot(aes(x=height, fill=treatment, color=treatment)) +
  geom_density(adjust=0.5, alpha=0.3) +
  geom_vline(data=mu.height, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_x_continuous(breaks=seq(0, 100, 10)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Height (cm)", y="Frequency") + 
  theme_bw() + theme(legend.position="none") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

mu.stem <- dat %>%
  filter(!germ==0, !stem_diam==0) %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(stem_diam))
dat %>% #stem diameter
  filter(!germ==0, !stem_diam==0) %>% #excluding stem_diam==0 because they were previously NA values
  ggplot(aes(x=stem_diam, fill=treatment, color=treatment))+
  geom_density(adjust=1.1, alpha=0.3) +
  geom_vline(data=mu.stem, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_x_continuous(breaks=seq(0, 15, 5)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Stem diameter (mm)", y="Frequency") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

#Graphs 4 - Spatial (plot) averages for fitness components + lifetime fitness ####

dat %>% #survival
  group_by(treatment, plot) %>%
  summarise(survival_sum=sum(flower), n=n()) %>%
  mutate(survival_prob=survival_sum/n, se.p=(sqrt(survival_prob*(1-survival_prob))/n)) %>%
  ggplot(aes(x = plot, y = survival_prob, group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=survival_prob-se.p, ymax=survival_prob+se.p), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Plot", y="Fecundity of flowering plants")

dat %>% #fecundity
  filter(flower==1) %>%
  group_by(treatment, plot) %>%
  summarise(fecundity_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = plot, y = fecundity_avg, group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=fecundity_avg-se, ymax=fecundity_avg+se), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Plot", y="Fecundity of flowering plants")

dat %>% #lifetime fitness
  group_by(treatment, plot) %>%
  summarise(pods_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = plot, y = pods_avg, group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=pods_avg-se, ymax=pods_avg+se), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Plot", y="Lifetime fitness")


#Graphs 5 - Maternal Correlations for Survival & Reproductive Components (Binary) ####

germ_mat <- dat %>% #overwintering survival
  group_by(treatment, famID, matID) %>%
  summarise(sum=sum(germ), n=n()) %>%
  mutate(germ_by_fam=(sum/n), se.p=(sqrt(germ_by_fam*(1-germ_by_fam))/n)) %>%#se.p = standard error of proportions
  select(treatment, famID, matID, germ_by_fam, se.p) %>%
  ungroup() %>%
  group_by(treatment, famID) %>%
  mutate(matID = c("A", "B")[sequence(n())]) %>%
  pivot_wider(names_from = matID, values_from = c(germ_by_fam, se.p))
germ_mat %>% #Note: 12 removed values because no paired maternal group
  ggplot(aes(x=germ_by_fam_A, y=germ_by_fam_B, group=treatment, color=treatment)) +
  geom_point(size=1.5) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=germ_by_fam_A, ymin=(germ_by_fam_B-se.p_B), ymax=(germ_by_fam_B+se.p_B), width=0.01)) +
  geom_errorbarh(aes(y=germ_by_fam_B, xmin=(germ_by_fam_A-se.p_A), xmax=(germ_by_fam_A+se.p_A))) +
  scale_color_manual("Environment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  ggtitle(label="A)") +
  labs(x="Average overwintering survival - Maternal subgroup A", y="Average overwintering survival - Maternal Subgroup B")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

flwr_mat <- dat %>% #spring-summer survival
  filter(germ==1) %>%
  group_by(treatment, famID, matID) %>%
  summarise(sum=sum(flower), n=n()) %>%
  mutate(flwr_by_fam=(sum/n), se.p=(sqrt(flwr_by_fam*(1-flwr_by_fam))/n)) %>% #se.p = standard error of proportions
  select(treatment, famID, matID, flwr_by_fam, se.p) %>%
  ungroup() %>%
  group_by(treatment, famID) %>%
  mutate(matID = c("A", "B")[sequence(n())]) %>%
  pivot_wider(names_from = matID, values_from = c(flwr_by_fam, se.p))
flwr_mat %>% #Note: 12 removed values because no paired maternal group
  ggplot(aes(x=flwr_by_fam_A, y=flwr_by_fam_B, group=treatment, color=treatment)) +
  geom_point(size=1.5) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=flwr_by_fam_A, ymin=(flwr_by_fam_B-se.p_B), ymax=(flwr_by_fam_B+se.p_B), width=0.01)) +
  geom_errorbarh(aes(y=flwr_by_fam_B, xmin=(flwr_by_fam_A-se.p_A), xmax=(flwr_by_fam_A+se.p_A))) +
  scale_color_manual("Environment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0.5,1)) +
  scale_x_continuous(limits=c(0.5,1)) +
  ggtitle(label="B)") +
  labs(x="Average summer survival - Maternal subgroup A", y="Average summer survival - Maternal subgroup B")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

seed_mat <- dat %>% #fruiting success
  filter(germ==1, flower==1) %>%
  group_by(treatment, famID, matID) %>%
  summarise(sum=sum(seed), n=n()) %>%
  mutate(seed_by_fam=(sum/n), se.p=(sqrt(seed_by_fam*(1-seed_by_fam))/n)) %>% #se.p = standard error of proportions
  select(treatment, famID, matID, seed_by_fam, se.p) %>%
  ungroup() %>%
  group_by(treatment, famID) %>%
  mutate(matID = c("A", "B")[sequence(n())]) %>%
  pivot_wider(names_from = matID, values_from = c(seed_by_fam, se.p))
seed_mat %>% #Note: 12 removed values because no paired maternal group
  ggplot(aes(x=seed_by_fam_A, y=seed_by_fam_B, group=treatment, color=treatment)) +
  geom_point(size=5, alpha=0.5) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=seed_by_fam_A, ymin=(seed_by_fam_B-se.p_B), ymax=(seed_by_fam_B+se.p_B), width=0.01)) +
  geom_errorbarh(aes(y=seed_by_fam_B, xmin=(seed_by_fam_A-se.p_A), xmax=(seed_by_fam_A+se.p_A))) +
  scale_color_manual("Environment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0.8,1.1)) +
  scale_x_continuous(limits=c(0.8,1.1)) +
  ggtitle(label="C)") +
  labs(x="Average fruiting success - Maternal subgroup A", y="Average fruiting success - Maternal subgroup B")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Graphs 7 - Maternal Correlations for Fitness Components ####

surv_mat <- dat %>% #survival
  group_by(treatment, famID, matID) %>%
  summarise(sum=sum(flower), n=n()) %>%
  mutate(surv_by_fam=(sum/n), se.p=(sqrt(surv_by_fam*(1-surv_by_fam))/n)) %>%#se.p = standard error of proportions
  select(treatment, famID, matID, surv_by_fam, se.p) %>%
  ungroup() %>%
  group_by(treatment, famID) %>%
  mutate(matID = c("A", "B")[sequence(n())]) %>%
  pivot_wider(names_from = matID, values_from = c(surv_by_fam, se.p))
mat_surv_fig <- surv_mat %>% #Note: 12 removed values because no paired maternal group
  ggplot(aes(x=surv_by_fam_A, y=surv_by_fam_B, group=treatment, color=treatment)) +
  geom_point(size=2, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=surv_by_fam_A, ymin=(surv_by_fam_B-se.p_B), ymax=(surv_by_fam_B+se.p_B), width=0.01, alpha = 0.5)) +
  geom_errorbarh(aes(y=surv_by_fam_B, xmin=(surv_by_fam_A-se.p_A), xmax=(surv_by_fam_A+se.p_A), alpha=0.5)) +
  scale_color_manual("Environment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  labs(x="Mean survival - Maternal subgroup A", y="Mean survival - Maternal subgroup B", tag = "A)")   +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.text = element_text(size=10), legend.title = element_text(size=10))
mat_surv_fig
ggsave(plot = mat_surv_fig, filename = "Routput/Figures/fig_s5a.png", width = 8, height = 8, units ="cm")

fecundity_mat <- dat %>% #fecundity
  filter(flower==1) %>%
  group_by(treatment, famID, matID) %>%
  summarise(pods_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  select(treatment, famID, matID, pods_avg, se) %>%
  ungroup %>%
  group_by(treatment, famID) %>%
  mutate(matID = c("A", "B")[sequence(n())]) %>%
  pivot_wider(names_from = matID, values_from = c(pods_avg, se))
mat_fecund_fig <- fecundity_mat %>% #Note: 12 removed values because no paired maternal group
  ggplot(aes(x=pods_avg_A, y=pods_avg_B, group=treatment, color=treatment)) +
  geom_point(size=2, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=pods_avg_A, ymin=(pods_avg_B-se_B), ymax=(pods_avg_B+se_B), width=0.01, alpha=0.5)) +
  geom_errorbarh(aes(y=pods_avg_B, xmin=(pods_avg_A-se_A), xmax=(pods_avg_A+se_A), alpha = 0.5)) +
  scale_color_manual("Environment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Mean fecundity - Maternal subgroup A", y="Mean fecundity - Maternal subgroup B", tag = "B)") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.text = element_text(size=10), legend.title = element_text(size=10))
mat_fecund_fig
ggsave(plot = mat_fecund_fig, filename = "Routput/Figures/fig_s5b.png", width = 8, height = 8, units ="cm")

pods_mat <- dat %>% #lifetime fitness
  group_by(treatment, famID, matID) %>%
  summarise(pods_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  select(treatment, famID, matID, pods_avg, se) %>%
  ungroup %>%
  group_by(treatment, famID) %>%
  mutate(matID = c("A", "B")[sequence(n())]) %>%
  pivot_wider(names_from = matID, values_from = c(pods_avg, se))
mat_lifetime_fig <- pods_mat %>% #Note: 12 removed values because no paired maternal group
  ggplot(aes(x=pods_avg_A, y=pods_avg_B, group=treatment, color=treatment)) +
  geom_point(size=2, alpha=0.8) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=pods_avg_A, ymin=(pods_avg_B-se_B), ymax=(pods_avg_B+se_B), width=0.1, alpha=0.5)) +
  geom_errorbarh(aes(y=pods_avg_B, xmin=(pods_avg_A-se_A), xmax=(pods_avg_A+se_A), alpha=0.5)) +
  scale_color_manual("Environment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Mean lifetime fitness - Maternal subgroup A", y="Mean lifetime fitness - Maternal subgroup B", tag = "C)")  +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.text = element_text(size=10), legend.title = element_text(size=10))
mat_lifetime_fig
ggsave(plot = mat_lifetime_fig, filename = "Routput/Figures/fig_s5c.png", width = 8, height = 8, units ="cm")


#Graphs 8 - Genetic Correlations between Select Traits and Lifetime Fitness ####

dat %>% #Lifetime Fitness against Overwintering Survival - Family Level
  group_by(treatment, famID) %>%
  summarise(sum=sum(germ), n=n(), pods_avg=mean(seed_pods)) %>%
  mutate(germ_avg=(sum/n)) %>%
  ggplot(aes(x=germ_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Average germination success", y="Average lifetime fitness") +
  ggtitle("A)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

dat %>% #Lifetime Fitness against Spring-summer Survival - Family Level
  group_by(treatment, famID) %>%
  summarise(sum=sum(flower), n=n(), pods_avg=mean(seed_pods)) %>%
  mutate(flwr_avg=(sum/n)) %>%
  ggplot(aes(x=flwr_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Average flowering success", y="Average lifetime fitness") +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

dat %>% #Lifetime Fitness against Flowering Clusters - Family Level
  filter(!germ==0, !flower==0, !flwr_clstr==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), flwrclstr_avg=mean(flwr_clstr)) %>%
  ggplot(aes(x=flwrclstr_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Average flowering clusters", y="Average seed pods") +
  ggtitle("B)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

dat %>% #Lifetime Fitness against Leaf Number - Family Level
  filter(!germ==0, !leaf==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), leaf_avg=mean(leaf)) %>%
  ggplot(aes(x=leaf_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Average number of leaves", y="Average seed pods") +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

dat %>% #Lifetime Fitness against Height - Family Level
  filter(!germ==0, !height==0) %>% #filtering out height==0 because those were originally NAs
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), height_avg=mean(height)) %>%
  ggplot(aes(x=height_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Average height (cm)", y="Average seed pods") +
  ggtitle("B)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

dat %>% #Lifetime Fitness against Stem Diameter - Family Level
  filter(!germ==0, !stem_diam==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), stem_avg=mean(stem_diam)) %>%
  ggplot(aes(x=stem_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Average stem diameter (mm)", y="Average seed pods") +
  ggtitle("B)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


dat %>% #Lifetime Fitness against Germination Time - Family Level
  filter(!germ==0) %>%
  group_by(treatment, famID) %>%
  summarise(germ_census_avg=mean(as.numeric(germ_census)), pods_avg=mean(seed_pods)) %>%
  ggplot(aes(x=germ_census_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Mean germination date (census #)", y="Number of seed pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


dat %>% #Lifetime Fitness against Flowering Time - Family Level
  filter(!germ==0) %>%
  group_by(treatment, famID) %>%
  summarise(flwr_census_avg=mean(as.numeric(flwr_census)), pods_avg=mean(seed_pods)) %>%
  ggplot(aes(x=flwr_census_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Mean Flowering Date (census #)", y="Number of Seed Pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Graphs 9 - Lifetime Fitness Summary by Proportion of Total Plants ####
sum.fit <- dat %>%
  group_by(treatment) %>%
  summarise(Planted=n()/n(), Germinated=sum(germ)/n(), Flowered=sum(flower)/n(), Seeded=sum(seed)/n())
sum.fit <- gather(sum.fit, lifetime, measurement, Planted:Seeded, factor_key=TRUE)
sum.fit <- sum.fit %>%
  mutate(value= case_when(
    lifetime == "Planted" ~ 1,
    lifetime == "Germinated" ~ 2,
    lifetime == "Flowered" ~ 3,
    lifetime == "Seeded" ~ 4))
sum.fit %>%
  ggplot(aes(x=value, y=measurement, group=treatment)) +
  geom_point(aes(color=treatment), size=4) + 
  geom_line(aes(color=treatment), size=1) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 1, 0.1)) +
  scale_x_discrete(limits=c("Planted", "Germinated", "Flowered", "Seeded")) +
  labs(x="Life history event", y="Relative proportion") + 
  ggtitle("A)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#GxE Reaction Norms
mu.lifetime <- dat %>% #reused from above
  group_by(treatment) %>%
  summarise(grp.mean=mean(seed_pods))
fam.lifetime <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(seed_pods))
fam.lifetime %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean family lifetime fitness") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 40, 10)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  #ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#GxE using boxplot
fam.lifetime %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.5) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=6,fill="darkgreen",alpha=0.5) +
  labs(x="Environment", y="Mean family absolute fitness") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 40, 10)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Previously used the following for individual level boxplot
#geom_boxplot() +
  #geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.5) +

mu.lifetime <- dat %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(seed_pods))
dat %>% #seed pods
  ggplot(aes(x=seed_pods, fill=treatment, color=treatment)) +
  geom_density(alpha=0.3) +
  geom_vline(data=mu.lifetime, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  scale_x_continuous(limits=c(0,300)) + #The range is HIGH! From 0 to 300.. not going to show this graph
  labs(x="Total seed pods per plant", y="Frequency")  + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

#Graphs 10 - Density plots using Family Means for fitness components ####
#reused some code from above 

mu_lifetime <- dat %>% #Lifetime Fitness
  group_by(treatment) %>%
  summarise(grp.mean=mean(seed_pods))
fam_lifetime <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(seed_pods))
fam_lifetime_fig <- fam_lifetime %>% 
  ggplot(aes(x=fam.mean, fill=treatment, color=treatment)) +
  geom_density(alpha=0.3) +
  geom_vline(data=mu_lifetime, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_color_manual("", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual("", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_x_continuous(limits=c(0,40)) +
  labs(x="Lifetime fitness", y="", tag = ("C)"), color = "", fill = "") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.text = element_text(size=10), legend.title = element_text(size=10),
        legend.position  = c(0.7, 0.7)) +
  #guides(color = guide_legend(nrow = 1)) + #to put it on one horizontal row.. but it doesn't fit in saved picture
  theme(axis.title.x = ggtext::element_markdown())
fam_lifetime_fig
ggsave(plot = fam_lifetime_fig, filename = "Routput/Figures/fig_3c.tif", width = 8, height = 8, units ="cm")

mu_flwr <- dat %>% #Survival (to flowering)
  group_by(treatment) %>%
  summarise(sum=sum(flower), n=n(), flwr_avg=sum/n)
fam_flwr <- dat %>% 
  group_by(treatment, famID) %>%
  summarise(sum=sum(flower), n=n()) %>%
  mutate(flwr_by_fam=(sum/n), se.p=(sqrt(flwr_by_fam*(1-flwr_by_fam))/n)) #se.p = standard of proportions
fam_flwr_fig <- fam_flwr %>% 
  ggplot(aes(x=flwr_by_fam, fill=treatment, color=treatment)) +
  geom_density(alpha=0.3) +
  geom_vline(data=mu_flwr, aes(xintercept=flwr_avg, color=treatment), linetype="dashed", size=1) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  scale_x_continuous(limits=c(0,1)) +
  labs(x="Survivorship", y="Frequency", tag=("A)")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.position = "none", legend.text = element_text(size=10), legend.title = element_text(size=10))
ggsave(plot = fam_flwr_fig, filename = "Routput/Figures/fig_3a.tif", width = 8, height = 8, units ="cm")


mu_pods <- dat %>% #Fecundity (of flowering plants)
  filter(germ==1, flower==1) %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(seed_pods))
fam_pods <- dat %>%
  filter(germ==1, flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(seed_pods))
fam_pods_fig <- fam_pods %>% 
  ggplot(aes(x = fam.mean, fill = treatment, color = treatment)) + 
  geom_density(alpha=0.3) +
  geom_vline(data=mu_pods, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  scale_x_continuous(limits=c(0,65)) +
  labs(x="Fecundity", y="", tag = ("B)")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=10), axis.title = element_text(size=10),
        legend.position = "none", legend.text = element_text(size=10), legend.title = element_text(size=10))
ggsave(plot = fam_pods_fig, filename = "Routput/Figures/fig_3b.png", width = 8, height = 8, units ="cm")

fig_3 <- fam_flwr_fig + fam_pods_fig
ggsave(plot = fig_3, filename = "Routput/Figures/fig_3.tif", width = 16, height = 9.5, units = "cm")

#Graphs 11 - Reaction Norms of Family Means ####

fam.surv <- dat %>% #Survival (to flowering)
  group_by(treatment, famID) %>%
  summarise(sum=sum(flower), n=n()) %>%
  mutate(fam.mean=(sum/n), se.p=(sqrt(fam.mean*(1-fam.mean))/n))
reaction_surv_fig <- fam.surv %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean family survival to flowering", tag = "A)") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 0.7, 0.1)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.text = element_text(size=8), legend.title = element_text(size=8))
reaction_surv_fig
ggsave(plot = reaction_surv_fig, filename = "Routput/Figures/fig_s6a.png", width = 8, height = 8, units ="cm")

  #x-y plot
fam.surv2 <- fam.surv %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.surv2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5)  +
  labs(x="Ambient Survival", y="Heated Survival")



fam.fecund <- dat %>% #Fecundity (of flowering plants)
  filter(germ==1, flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(seed_pods))
reaction_fecund_fig <- fam.fecund %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean family fecundity", tag = "B)") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 60, 10)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.text = element_text(size=8), legend.title = element_text(size=8))
reaction_fecund_fig
ggsave(plot = reaction_fecund_fig, filename = "Routput/Figures/fig_s6b.png", width = 8, height = 8, units ="cm")

  #x-y plot
fam.fecund2 <- fam.fecund %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.fecund2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient Fecundity", y="Heated Fecundity")



fam.absolute <- dat %>% #Lifetime Fitness
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(seed_pods))
reaction_lifetime_fig <- fam.absolute %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean family lifetime fitness", tag = "C)") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 40, 10)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.text = element_text(size=8), legend.title = element_text(size=8))
reaction_lifetime_fig
ggsave(plot = reaction_lifetime_fig, filename = "Routput/Figures/fig_s6c.png", width = 8, height = 8, units ="cm")

  #x-y plot
fam.absolute2 <- fam.absolute %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.absolute2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5)  +
  labs(x="Ambient Absolute Fitness", y="Heated Absolute Fitness")

fam.germ <- dat %>% #Overwintering Survival to Germination (Germination Success)
  group_by(treatment, famID) %>%
  summarise(sum=sum(germ), n=n()) %>%
  mutate(fam.mean=(sum/n), se.p=(sqrt(fam.mean*(1-fam.mean))/n))
fam.germ %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean family germination success") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 0.8, 0.1)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.germ2 <- fam.germ %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.germ2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient germination success", y="Heated germination success")

fam.flower <- dat %>% #Spring-Summer Survival to Flowering (Flowering Success)
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(flower), n=n()) %>%
  mutate(fam.mean=(sum/n), se.p=(sqrt(fam.mean*(1-fam.mean))/n))
fam.flower %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean spring-summer survival success") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0.5, 1, 0.1)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.flower2 <- fam.flower %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.flower2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5)  +
  labs(x="Ambient flowering success", y="Heated flowering success")

fam.seed <- dat %>% #Fruiting Success
  filter(germ==1, flower==1) %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(seed), n=n()) %>%
  mutate(fam.mean=(sum/n), se.p=(sqrt(fam.mean*(1-fam.mean))/n))
fam.seed %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean fruiting success") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0.9, 1, 0.02)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.seed2 <- fam.seed %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.seed2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient fruiting Success", y="Heated fruiting success")

fam.flwrclstr <- dat %>% #Flowering Clusters Number
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(flwr_clstr))
fam.flwrclstr %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean family flowering cluster number") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 4, 1)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.flwrclstr2 <- fam.flwrclstr %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.flwrclstr2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient flower clusters", y="Heated flower clusters")

fam.leaf <- dat %>% #Leaf Number
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(leaf))
fam.leaf %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean family leaf number") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 8, 2)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.leaf2 <- fam.leaf %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.leaf2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient leaf", y="Heated leaf")


fam.height <- dat %>% #Height
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(height))
fam.height %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean family height (cm)") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 40, 5)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))
fam.height2 <- fam.height %>%
  pivot_wider(id_cols=famID, names_from=treatment, values_from=fam.mean) %>%
  select(famID, A, H)
fam.height2 %>%
  ggplot(aes(x=A, y=H)) + geom_point(size=2, alpha=0.5) +
  labs(x="Ambient height", y="Heated height")

#Stem Diameter
fam.stem <- dat %>%
  filter(germ==1) %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(stem_diam))
fam.stem %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment, group=famID)) + 
  geom_point(size=3, alpha=0.5) +
  geom_line(color="black") +
  labs(x="Environment", y="Mean stem diameter (mm)") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 3, 0.5)) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
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

#Graphs 12 - Evaluating Seed Pod as Proxy for Fitness ####

seed_proxy_fig <- dat %>% 
  filter(plot %in% c("1", "2", "10", "11"), !seed_pods==0) %>%
  ggplot(aes(x = seed_pods, y = seed_number, color = treatment)) +
  geom_point(size=3, alpha = 0.5) +
  geom_smooth(method="lm", formula=y~x, se=F) +
  scale_color_manual("Environment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  scale_x_continuous(breaks = seq(0, 300, 50)) +
  scale_y_continuous(breaks = seq(0, 3500, 500)) + 
  labs(x="Total seed pods per plant", y="Seed number per plant")  +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 1), axis.ticks = element_line(size = 1), 
        axis.ticks.length = unit(0.25, "cm"), plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=8), axis.title = element_text(size=8),
        legend.position = "none", legend.box = "vertical", 
        legend.text = element_text(size=8), legend.title = element_text(size=8))
seed_proxy_fig
ggsave(plot = seed_proxy_fig, filename = "Routput/Figures/fig_s3.png", width = 12, height = 9.5, units = "cm")

cor.test(dat$seed_pods[dat$plot=="1" | dat$plot=="2" | dat$plot=="10" | dat$plot=="11" | dat$treatment=="A"], 
         dat$seed_number[dat$plot=="1" | dat$plot=="2" | dat$plot=="10" | dat$plot=="11" | dat$treatment=="A"],
         method = "pearson")
Ar2 <- 0.7575373^2
Ar2

cor.test(dat$seed_pods[dat$plot=="1" | dat$plot=="2" | dat$plot=="10" | dat$plot=="11" | dat$treatment=="H"], 
         dat$seed_number[dat$plot=="1" | dat$plot=="2" | dat$plot=="10" | dat$plot=="11" | dat$treatment=="H"],
         method = "pearson")
Hr2 <- 0.5611888^2
Hr2

cor.test(dat$seed_pods[dat$plot=="1" | dat$plot=="2" | dat$plot=="10" | dat$plot=="11" ], 
         dat$seed_number[dat$plot=="1" | dat$plot=="2" | dat$plot=="10" | dat$plot=="11" ],
         method = "pearson")
Tr2 <- 0.911911^2
Tr2
