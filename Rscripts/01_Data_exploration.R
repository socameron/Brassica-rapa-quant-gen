#### PROJECT: Brassica rapa GxE Study (Data collected by Cameron So 2019 at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Clean data and observe distributions prior to aster analysis
#### AUTHOR: Cameron So
#### DATE LAST MODIFIED: 2019/09/24

.########################################################################.
 ##############      PACKAGE INSTALLATION AND IMPORT      ###############
.########################################################################.

#Installing and Loading Necessary Packages
install.packages("ggplot2")
install.packages("dplyr")
install.packages("lme4")
install.packages("gplots")
install.packages("tidyverse")
install.packages("zoo")
install.packages("car")
install.packages("multcompView")
install.packages("lsmeans")
install.packages("FSA")
install.packages("nlme")
install.packages("MASS")
install.packages("fitdistrplus")
library(ggplot2)
library(dplyr)
library(gplots)
library(tidyverse)
library(zoo)
library(lme4)
library(car)
library(multcompView)
#library(lsmeans) #lsmeans now integrated into emmeans
library(emmeans)
library(FSA)

#Importing Data using readr from tidyverse (Base R confuses the class for certain vectors). Excluding notes from the data import
# d for double, f for factor, i for integer
col_types_list <- cols_only(posID = "d", individual = "d", plot = col_factor(levels=c(1:12)),
                            matID = "f", patID = "f", famID = col_factor(levels=c(1:62)),
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
                            germ_bin = "?", flwr_bin = "?", seed_bin = "?",
                            seed_dmg = "d", seed_undmg = "d", seed_pods = "d",
                            height = "d", stem_diam = "d"
                            )
data <- read_csv("Rdata/heatarrays_brassica_2019_data.csv", col_names = TRUE, na = "NA", 
                 col_types=col_types_list)
lapply(data, class)

#Note: I am not using the calendar date data. However, if you wish to use that data, use package readr from tidyverse and function read_csv to parse the date columns correctly. 

.########################################################################.
 ######################      CLEANING DATASET      ######################
.########################################################################.

#Number of Positions wth Doubles
data %>% summarise(double=sum(double, na.rm=TRUE), n()) #229 doubles

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
data <- data %>% replace(., is.na(.), "0")

#Using tidyverse to convert dataframe into long format
data_phen<- data %>%
  select(posID, individual, plot, matID, patID, famID, mat_group, treatment,
         germ_juln_date, flwr_juln_date, flwr_compl_juln, 
         leaf1, leaf2, leaf3, leaf4, leaf5, leaf6, leaf7, leaf8, leaf9,
         flwr_clstr1, bud_clstr1, flwr_clstr2, bud_clstr2, flwr_clstr3, bud_clstr3,
         flwr_clstr4, bud_clstr4, flwr_clstr5, bud_clstr5, flwr_clstr6, bud_clstr6,
         flwr_clstr7, bud_clstr7, flwr_clstr8, bud_clstr8, seed_pods, seed_dmg, seed_undmg) %>%
  gather("leaf1", "leaf2", "leaf3", "leaf4", "leaf5", "leaf6", "leaf7", "leaf8", "leaf9", "flwr_clstr1", "flwr_clstr2", "flwr_clstr3", "flwr_clstr4", "flwr_clstr5", "flwr_clstr6", "flwr_clstr7", "flwr_clstr8",
         "bud_clstr1", "bud_clstr2", "bud_clstr3", "bud_clstr4", "bud_clstr5", "bud_clstr6", "bud_clstr7", "bud_clstr8",
         key = "phenotype", value = "phen_num")  #This can be used to see the change across time for various phenotypes.
         
data_leaf <- data %>%  
  select(posID, individual, plot, matID, patID, famID, mat_group, treatment, 
         germ_juln_date, flwr_juln_date, flwr_compl_juln, germ_bin, flwr_bin, seed_bin,
         leaf1, leaf2, leaf3, leaf4, leaf5, leaf6, leaf7, leaf8, leaf9) %>%
  gather(key="census", value="leaf", leaf1:leaf9)

data_flwrclstr <- data %>%
  select(posID, individual, plot, matID, patID, famID, mat_group, treatment, 
         germ_juln_date, flwr_juln_date, flwr_compl_juln, germ_bin, flwr_bin, seed_bin,
         flwr_clstr1, flwr_clstr2, flwr_clstr3, flwr_clstr4, flwr_clstr5, flwr_clstr6, flwr_clstr7, flwr_clstr8) %>%
  gather(key="census", value="flwr_clstr", flwr_clstr1:flwr_clstr8) #for flowering time analysis
         
#Setting Appropriate Classes to Numeric/Double as classes were altered during cleanup
lapply(data, class)
data$posID <- as.factor(data$posID)
data$individual <- as.integer(data$individual)
data[,11:48] <- lapply(data[,11:48], as.numeric) #setting all neccessary parameters to class numeric

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
  select(posID, individual, plot, matID, patID, famID, mat_group, treatment,
         germ_juln_date, flwr_juln_date, flwr_compl_juln,
         leaf, flwr_clstr, bud_clstr, germ, flower, seed, seed_pods, seed_dmg, seed_undmg, height, stem_diam)

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

#Evaluating individual errors
head(subset(dat, leaf==0&flwr_clstr>0&seed_pods>0)) # P7 372 & 392 are OK b/c no leaf data.
head(subset(dat, leaf==0&flwr_clstr==0&seed_pods>0)) #P2I2 needs to be checked in the collected bags
head(subset(dat, flwr_clstr==0&seed_pods>0)) #P3 441 and P7 157 to be deleted; #P7 170 & 329 no data 
head(subset(dat, leaf==0&flwr_clstr>0))
head(subset(dat, leaf==0&seed_pods>0))

#Fixing errors (removing 4 individuals)
dat <- dat[!(dat$flwr_clstr==0&dat$seed_pods>0),]
#When checking for errors again, only 2 errors because P7 372 & 392 remain in the dataset

#Looking for inconsistencies errors between binary and continuous traits (flwr_clstr & flower, and seed and seed_pods)
nrow(subset(dat, germ==0&flwr_clstr>0)) # 0 errors
nrow(subset(dat, germ==0&seed_pods>0)) # 0 errors
nrow(subset(dat, flower==0&flwr_clstr>0)) # 0 errors
nrow(subset(dat, flower==0&seed_pods>0)) # 0 errors
nrow(subset(dat, seed==0&seed_pods>0)) # 0 errors

.########################################################################.
 ######################      DATA EXPLORATION      ######################
.########################################################################.

#Setting plot theme for ggplot
theme_set(theme_classic(base_size = 10))

#Tally and Proportions by Treatment for Binary Traits (germ, flower, seed) ####
#Consequent lifetime fitness traits are conditional on previous traits (e.g flowering success if germ=1)
#Also including z-test proportional tests
dat %>% #germination
  group_by(treatment) %>%
  summarise(germ=sum(germ), n=n()) %>%
  mutate(germ_perc=(100*germ/n))
prop.test(x = c(1589, 1575), n=c(3457, 3348))

dat %>% #flowering success
  group_by(treatment) %>%
  filter(!germ==0) %>%
  summarise(flower=sum(flower), n=n()) %>%
  mutate(flwr_perc=(100*flower/n))
prop.test(x = c(1421, 1481), n=c(1589, 1575))

dat %>% #seed pod maturation
  group_by(treatment) %>%
  filter(!germ==0, !flower==0) %>%
  summarise(seed=sum(seed), n=n()) %>%
  mutate(seed_perc=(100*seed/n))
prop.test(x = c(1411, 1475), n=c(1421, 1481))


#Tally and Averages by Treatment for Continuous Traits (leaf, flwr_clstr, seed_pods) ####
dat %>% #flowering clusters
  filter(!germ==0) %>%
  group_by(treatment) %>%
  summarise(flower_avg=mean(flwr_clstr), n=n(), sd=sd(flwr_clstr), se=(sd/(sqrt(n))))

dat %>% #seed pod number
  filter(!germ==0, !flower==0) %>%
  group_by(treatment) %>%
  summarise(pod_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n))))

dat %>% #leaf number
  filter(!germ==0) %>%
  group_by(treatment) %>%
  summarise(leaf_avg=mean(leaf), n=n(), sd=sd(leaf), se=(sd/(sqrt(n))))

dat %>% #height
  filter(!germ==0, !height==0) %>%
  group_by(treatment) %>%
  summarise(height_avg=mean(height), n=n(), sd=sd(height), se=(sd/sqrt(n)))

dat %>% #stem diameter
  filter(!germ==0, !stem_diam==0) %>%
  group_by(treatment) %>%
  summarise(stem_avg=mean(stem_diam), n=n(), sd=sd(height), se=(sd/sqrt(n)))

dat %>% #absolute fitness
  group_by(treatment, famID) %>%
  summarise(pod_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n))))


#Tally for Bird Damage ####
dat %>% 
  filter(!germ==0, !flower==0, !seed==0) %>%
  group_by(treatment, famID) %>%
  summarise(bird_plants=sum(seed_dmg>0), undmg_plants=sum(seed_undmg>=0 & seed_dmg==0), total_plants=sum(seed_pods>=0))

dat %>%
  filter(!germ==0, !flower==0, !seed==0) %>%
  group_by(treatment) %>%
  summarise(seeded_plants=sum(seed==1))

# Graphs - Flowering Season Length ####
dat %>% #graphing by treatment only
  filter(!flwr_juln_date==0, !flwr_compl_juln==0, !flower==0, !germ==0) %>%
  mutate(flowering_time=(flwr_compl_juln - flwr_juln_date)) %>%
  filter(!flowering_time==0) %>% #Removed individuals with flowering time = 0 because we did not collect the flowering dates correctly in the field
  group_by(treatment) %>%
  ggplot(aes(x=treatment, y=flowering_time)) +
  geom_boxplot(aes(fill=treatment), lwd=1) +
  labs(x="Treatment", y="Flowering Season Length (Days)") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0,40)) +
  scale_fill_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient, Heated")) + 
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
    labs(x="Treatment", y="Flowering Season Length (Days)") +
    scale_x_discrete(labels=c("Ambient", "Heated")) + 
    scale_fill_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient, Heated")) + 
    theme_bw() + theme(legend.position="none") +
    coord_flip()

### Flowering Date 
dat %>% #time to first flower after germination
  filter(!flwr_juln_date==0, !flwr_compl_juln==0, !flower==0, !germ==0) %>%
  mutate(first_flower=(flwr_juln_date - germ_juln_date)) %>%
  group_by(treatment) %>%
  ggplot(aes(x=treatment, y=first_flower)) +
  geom_boxplot(aes(fill=treatment), lwd=1) +
  labs(x="Treatment", y="Days to First Flower") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 80, 10)) +
  scale_fill_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient, Heated")) + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm")) +
  geom_jitter(aes(color=treatment), alpha=0.2, shape=16, position=position_jitter(0.2)) +
  scale_color_manual(values=c("dodgerblue4", "tomato4")) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=6,fill="darkgreen",alpha=0.5)  +
  coord_flip()

### Maternal Sub-group Tally
dat %>% 
  group_by(mat_group, famID) %>%
  summarise(count = n_distinct(matID)) %>%
  summarise(mat_subgroups = sum(count)) #M1 has 39 mat-subgroups; M2 has 79 mat-subgroups

dat %>%
  group_by(mat_group) %>%
  summarise(sibships = n_distinct(famID)) #M1 has 22 sibships; M2 has 40 sibships


#Graphs - Histogram Distribution for Continuous Traits (INCLUDING germ = 0) ####
#Including ENTIRE Dataset (where germ = 0)
dat %>% #leaf number
  ggplot(aes(x=leaf, fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of Leaves", y="Frequency") 

dat %>% #flowering clusters
  ggplot(aes(x=flwr_clstr, fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of Flowering Clusters", y="Frequency")

dat %>% #seed pod
  ggplot(aes(x=seed_pods, fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of Seed Pods", y="Frequency")

dat %>% #height
  ggplot(aes(x=height, fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Height (cm)", y="Frequency")

dat %>% #stem diameter
  ggplot(aes(x=stem_diam, fill=treatment)) +
  geom_histogram(bins=30, position="dodge") +
  scale_fill_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Stem Diameter (mm)", y="Frequency")

#Graphs - Density Plots for Distribution for Continuous Traits (EXCLUDING germ = 0) ####
#Including SUBSET Dataset (EXCLUDING where germ = 0)
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
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Number of Leaves", y="Frequency")  + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

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
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Number of Flowering Clusters", y="Frequency") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

mu.pods <- dat %>%
  filter(!germ==0, !flower==0) %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(seed_pods))
dat %>% #seed pods
  filter(!germ==0, !flower==0) %>%
  ggplot(aes(x=seed_pods, fill=treatment, color=treatment)) +
  geom_density(alpha=0.3) +
  geom_vline(data=mu.pods, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_x_continuous(breaks=seq(0, 300, 50)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Number of Seed Pods from Flowering Plants", y="Frequency")  + 
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
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
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
  scale_x_continuous(breaks=seq(0, 40, 5)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Stem Diameter (mm)", y="Frequency") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

mu.area <- dat %>%
  filter(!germ==0, !stem_diam==0) %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(pi*((stem_diam/2)^2)))
dat %>% #cross sectional area
  filter(!germ==0, !stem_diam==0) %>% #excluding stem_diam==0 because they were previously NA values
  ggplot(aes(x=pi*((stem_diam/2)^2), fill=treatment, color=treatment))+
  geom_density(alpha=0.3) +
  geom_vline(data=mu.area, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  labs(x="Cross Sectional Area (mm)", y="Frequency") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

#Graphs - Trait distribution by PLOT for continuous Traits (leaf, flwr_clstr, seed_pods) ####
dat %>% #leaves
  group_by(treatment, plot) %>%
  summarise(leaf_avg=mean(leaf), n=n(), sd=sd(leaf), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = plot, y = leaf_avg, group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=leaf_avg-se, ymax=leaf_avg+se), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Plot", y="Average Number of Leaves")

dat %>% #flower clusters
  group_by(treatment, plot) %>%
  summarise(flwr_clstr_avg=mean(flwr_clstr), n=n(), sd=sd(flwr_clstr), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = plot, y = flwr_clstr_avg, group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=flwr_clstr_avg-se, ymax=flwr_clstr_avg+se), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Plot", y="Average Number of Flowering Clusters")

dat %>% #seed pods
  group_by(treatment, plot) %>%
  summarise(pods_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = plot, y = pods_avg, group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=pods_avg-se, ymax=pods_avg+se), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Plot", y="Average Number of Seed Pods")


#Graphs - Trait distribution by FAMILY/SIBSHIP for binary traits (germ, flower, seed) ####
germ_fam <- dat %>% #germination
  group_by(treatment, famID) %>%
  summarise(sum=sum(germ), n=n()) %>%
  mutate(germ_by_fam=(sum/n), se.p=(sqrt(germ_by_fam*(1-germ_by_fam))/n)) #se.p = standard error of proportions

  #Calculating average for each treatment
  germ_summary <- dat %>%
    group_by(treatment) %>%
    summarise(sum=sum(germ), n=n(), germ_avg=sum/n)
    germ_A <- 0.4596471
    germ_H <- 0.4704301
  
  #Plotting in descending rank order of germination by family with AMBIENT as base
  germ_fam %>%
    group_by(famID) %>%
    mutate(germ_by_fam2 = if_else(treatment=="A", 1, 0)*germ_by_fam) %>% 
    ungroup %>%
    ggplot(aes(x=reorder(famID, desc(germ_by_fam2)), y=germ_by_fam, group=treatment, color=treatment)) +
    geom_point(size=1.5) +
    geom_hline(yintercept=germ_A, color="steelblue2") +
    geom_hline(yintercept=germ_H, color="tomato2") +
    geom_errorbar(aes(x=famID, ymin=(germ_by_fam-se.p), ymax=(germ_by_fam+se.p), width=0.5)) + 
    geom_smooth(se=F) +
    scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
    labs(x="Sibship", y="Average Germination Success")  +
    theme_bw() + theme(legend.position="none") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.5), 
          axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
    theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))
  
flwr_fam <- dat %>% #flowering success
    group_by(treatment, famID) %>%
    summarise(sum=sum(flower), n=n()) %>%
    mutate(flwr_by_fam=(sum/n), se.p=(sqrt(flwr_by_fam*(1-flwr_by_fam))/n)) #se.p = standard error of proportions  
  
  #Calculating average for each treatment
  flwr_summary <- dat %>%
    group_by(treatment) %>%
    summarise(sum=sum(flower), n=n(), flwr_avg=sum/n)
    flwr_A <- 0.4110500
    flwr_H <- 0.4423536
  
  #Plotting in descending rank order of flowering success by family with AMBIENT as base
  flwr_fam %>%
    group_by(famID) %>%
    mutate(flwr_by_fam2 = if_else(treatment=="A", 1, 0)*flwr_by_fam) %>% 
    ungroup %>%
    ggplot(aes(x=reorder(famID, desc(flwr_by_fam2)), y=flwr_by_fam, group=treatment, color=treatment)) +
    geom_point(size=1.5) +
    geom_hline(yintercept=flwr_A, color="steelblue2") +
    geom_hline(yintercept=flwr_H, color="tomato2") +
    geom_errorbar(aes(x=famID, ymin=(flwr_by_fam-se.p), ymax=(flwr_by_fam+se.p), width=0.5)) + 
    geom_smooth(se=F) +
    scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
    labs(x="Sibship", y="Average Flowering Success")  +
    theme_bw() + theme(legend.position="none") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
    theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

  
seed_fam <- dat %>% #seed maturation success
    group_by(treatment, famID) %>%
    summarise(sum=sum(seed), n=n()) %>%
    mutate(seed_by_fam=(sum/n), se.p=(sqrt(seed_by_fam*(1-seed_by_fam))/n)) #se.p = standard error of proportions  
  
  #Calculating average for each treatment
  seed_summary <- dat %>%
    group_by(treatment) %>%
    summarise(sum=sum(seed), n=n(), seed_avg=sum/n)
  seed_A <- 0.4084466
  seed_H <- 0.4411589
  
  #Plotting in descending rank order of seed maturation success by family with AMBIENT as base
  seed_fam %>%
    group_by(famID) %>%
    mutate(seed_by_fam2 = if_else(treatment=="A", 1, 0)*seed_by_fam) %>% 
    ungroup %>%
    ggplot(aes(x=reorder(famID, desc(seed_by_fam2)), y=seed_by_fam, group=treatment, color=treatment)) +
    geom_point(size=1.5) +
    geom_hline(yintercept=seed_A, color="steelblue2") +
    geom_hline(yintercept=seed_H, color="tomato2") +
    geom_errorbar(aes(x=famID, ymin=(seed_by_fam-se.p), ymax=(seed_by_fam+se.p), width=0.5)) + 
    geom_smooth(se=F) +
    scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
    labs(x="Sibship", y="Average Seed Maturation Success")  +
    theme_bw() + theme(legend.position="none") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.5), 
          axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
    theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))
  
  
#Graphs - Trait distributions by FAMILY/SIBSHIP for continuous traits (leaf, flwr_clstr, seed_pods) ####
  
leaf1 <- dat %>% #calculating leaf average per treatment
    group_by(treatment) %>%
    summarise(leaf_avg=mean(leaf))
    leaf_A <- 2.986404
    leaf_H <- 3.638590
  
   dat %>% #leaf number
    group_by(treatment, famID) %>%
    summarise(leaf_avg=mean(leaf), n=n(), sd=sd(leaf), se=(sd/(sqrt(n)))) %>%
    mutate(leaf_avg2 = if_else(treatment=="A", 1, 0)*leaf_avg) %>% 
    ungroup %>%
    ggplot(aes(x=reorder(famID, desc(leaf_avg2)), y=leaf_avg, 
             group=treatment, color=treatment)) +
    geom_errorbar(aes(ymin=leaf_avg-se, ymax=leaf_avg+se), width=0.2, size=0.5) + 
    geom_point(aes(color=treatment)) + 
    geom_smooth(se=F) + 
    geom_hline(yintercept=leaf_A, color="steelblue2") +
    geom_hline(yintercept=leaf_H, color="tomato2") +
    scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
    labs(x="Sibship", y="Average Leaf Number")  +
    theme_bw() + theme(legend.position="none") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.5), 
          axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
    theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))


flwr_clstr1 <- dat %>% #calculating flower cluster average per treatment
    group_by(treatment) %>%
    summarise(flwrclstr_avg=mean(flwr_clstr))
    flwrclstr_A <- 0.753559
    flwrclstr_H <- 1.145618

  
    dat %>% #flower clusters
      group_by(treatment, famID) %>%
      summarise(flwr_clstr_avg=mean(flwr_clstr), n=n(), sd=sd(flwr_clstr), se=(sd/(sqrt(n)))) %>%
      mutate(flwr_clstr_avg2 = if_else(treatment=="A", 1, 0)*flwr_clstr_avg) %>% 
      ungroup %>%
      ggplot(aes(x=reorder(famID, desc(flwr_clstr_avg2)), y=flwr_clstr_avg, 
                 group=treatment, color=treatment)) + 
      geom_errorbar(aes(ymin=flwr_clstr_avg-se, ymax=flwr_clstr_avg+se), width=0.2, size=0.5) + 
      geom_point(aes(color=treatment)) + 
      geom_smooth(se=F) +
      geom_hline(yintercept=flwrclstr_A, color="steelblue2") +
      geom_hline(yintercept=flwrclstr_H, color="tomato2") +
      scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
      labs(x="Sibship", y="Average Flowering Clusters Number")  +
      theme_bw() + theme(legend.position="none") + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.line = element_line(colour = "black", size = 1.5), 
            axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
      theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

seed_pods1 <- dat %>% #calculating seed pod average per treatment
    group_by(treatment) %>%
    summarise(pods_avg=mean(seed_pods))
    pods_A <- 4.985826
    pods_H <- 9.261649

    dat %>% #seed pods
      group_by(treatment, famID) %>%
      summarise(pods_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
      mutate(pods_avg2 = if_else(treatment=="A", 1, 0)*pods_avg) %>% 
      ungroup %>%
      ggplot(aes(x = reorder(famID, desc(pods_avg2)), y = pods_avg,
                 group = treatment, color = treatment)) + 
      geom_errorbar(aes(ymin=pods_avg-se, ymax=pods_avg+se), width=0.2, size=0.5) + 
      geom_point(aes(color=treatment)) + 
      geom_smooth(se=F) +
      geom_hline(yintercept=pods_A, color="steelblue2") +
      geom_hline(yintercept=pods_H, color="tomato2") +
      scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
      labs(x="Sibship", y="Average Seed Pod Number")  +
      theme_bw() + theme(legend.position="none") + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.line = element_line(colour = "black", size = 1.5), 
            axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
      theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

#Excluding individuals that did NOT produce any seed pods at all
seed_pods2 <- dat %>% #calculating flower cluster average per treatment
    group_by(treatment) %>%
    filter(seed==1) %>%
    summarise(pods_avg=mean(seed_pods))
    pods_A2 <- 12.20680
    pods_H2 <- 20.99391

    dat %>% #seed pods, excluding individuals that did NOT produce seed pods
      group_by(treatment, famID) %>%
      filter(seed==1) %>%
      summarise(pods_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
      mutate(pods_avg2 = if_else(treatment=="A", 1, 0)*pods_avg) %>% 
      ungroup %>%
      ggplot(aes(x = reorder(famID, desc(pods_avg2)), y = pods_avg,
                 group = treatment, color = treatment)) + 
      geom_errorbar(aes(ymin=pods_avg-se, ymax=pods_avg+se), width=0.2, size=0.5) + 
      geom_point(aes(color=treatment)) + 
      geom_smooth(se=F) +
      geom_hline(yintercept=pods_A2, color="steelblue2") +
      geom_hline(yintercept=pods_H2, color="tomato2") +
      scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
      labs(x="Sibship", y="Average Seed Pod Number")  +
      theme_bw() + theme(legend.position="none") + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.line = element_line(colour = "black", size = 1.5), 
            axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
      theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

    
height1 <- dat %>% #calculating height average per treatment
    group_by(treatment) %>%
    summarise(height_avg=mean(height))
    height_A <- 12.84842
    height_H <- 15.76912
    
    dat %>% #height
      group_by(treatment, famID) %>%
      summarise(height_avg=mean(height), n=n(), sd=sd(height), se=(sd/(sqrt(n)))) %>%
      mutate(height_avg2 = if_else(treatment=="A", 1, 0)*height_avg) %>% 
      ungroup %>%
      ggplot(aes(x = reorder(famID, desc(height_avg2)), y = height_avg,
                 group = treatment, color = treatment)) + 
      geom_errorbar(aes(ymin=height_avg-se, ymax=height_avg+se), width=0.2, size=0.5) + 
      geom_point(aes(color=treatment)) + 
      geom_smooth(se=F) +
      geom_hline(yintercept=height_A, color="steelblue2") +
      geom_hline(yintercept=height_H, color="tomato2") +
      scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
      labs(x="Sibship", y="Average Height (cm)")  +
      theme_bw() + theme(legend.position="none") + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.line = element_line(colour = "black", size = 1.5), 
            axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
      theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))
    
stem_diam1 <- dat %>% #calculating stem diameter average per treatment
    group_by(treatment) %>%
    summarise(stem_avg=mean(stem_diam))
    stem_A <- 0.9224761
    stem_H <- 1.1099164
    
    dat %>% #stem diameter
      group_by(treatment, famID) %>%
      summarise(stem_avg=mean(stem_diam), n=n(), sd=sd(stem_diam), se=(sd/(sqrt(n)))) %>%
      mutate(stem_avg2 = if_else(treatment=="A", 1, 0)*stem_avg) %>% 
      ungroup %>%
      ggplot(aes(x = reorder(famID, desc(stem_avg2)), y = stem_avg,
                 group = treatment, color = treatment)) + 
      geom_errorbar(aes(ymin=stem_avg-se, ymax=stem_avg+se), width=0.2, size=0.5) + 
      geom_point(aes(color=treatment)) + 
      geom_smooth(se=F) +
      geom_hline(yintercept=stem_A, color="steelblue2") +
      geom_hline(yintercept=stem_H, color="tomato2") +
      scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
      labs(x="Sibship", y="Average Stem Diameter (mm)")  +
      theme_bw() + theme(legend.position="none") + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.line = element_line(colour = "black", size = 1.5), 
            axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
      theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))    
    
    
    
    
###!!!!NOTE: Graphs INCLUDE maternal subgroups !!!!####

#Graphs - Trait distribution by FAMILY/SIBSHIP for BINARY traits (germ, flower, seed) ####
germ_fam2 <- dat %>% #germination
  group_by(treatment, famID, matID) %>%
  summarise(sum=sum(germ), n=n()) %>%
  mutate(germ_by_fam=(sum/n), se.p=(sqrt(germ_by_fam*(1-germ_by_fam))/n)) #se.p = standard error of proportions

  #Plotting in descending rank order of germination by family with AMBIENT as base
  germ_fam2 %>%
  group_by(famID) %>%
    mutate(germ_by_fam2 = if_else(treatment=="A", 1, 0)*germ_by_fam) %>% 
    ungroup %>%
    ggplot(aes(x=reorder(famID, desc(germ_by_fam2)), y=germ_by_fam, group=treatment, color=treatment)) +
    geom_point(size=1.5) +
    geom_hline(yintercept=germ_A, color="steelblue2") +
    geom_hline(yintercept=germ_H, color="tomato2") +
    geom_errorbar(aes(x=famID, ymin=(germ_by_fam-se.p), ymax=(germ_by_fam+se.p), width=0.5)) + 
    geom_smooth(se=F) +
    scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
    labs(x="Sibship", y="Average Germination Success")  +
    theme_bw() + theme(legend.position="none") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.5), 
          axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
    theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

flwr_fam2 <- dat %>% #flowering success
  group_by(treatment, famID, matID) %>%
  summarise(sum=sum(flower), n=n()) %>%
  mutate(flwr_by_fam=(sum/n), se.p=(sqrt(flwr_by_fam*(1-flwr_by_fam))/n)) #se.p = standard error of proportions  

  #Plotting in descending rank order of flowering success by family with AMBIENT as base
  flwr_fam2 %>%
    group_by(famID) %>%
    mutate(flwr_by_fam2 = if_else(treatment=="A", 1, 0)*flwr_by_fam) %>% 
    ungroup %>%
    ggplot(aes(x=reorder(famID, desc(flwr_by_fam2)), y=flwr_by_fam, group=treatment, color=treatment)) +
    geom_point(size=1.5) +
    geom_hline(yintercept=flwr_A, color="steelblue2") +
    geom_hline(yintercept=flwr_H, color="tomato2") +
    geom_errorbar(aes(x=famID, ymin=(flwr_by_fam-se.p), ymax=(flwr_by_fam+se.p), width=0.5)) + 
    geom_smooth(se=F) +
    scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
    labs(x="Sibship", y="Average Flowering Success")  +
    theme_bw() + theme(legend.position="none") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.5), 
          axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
    theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))


seed_fam2 <- dat %>% #seed maturation success
  group_by(treatment, famID, matID) %>%
  summarise(sum=sum(seed), n=n()) %>%
  mutate(seed_by_fam=(sum/n), se.p=(sqrt(seed_by_fam*(1-seed_by_fam))/n)) #se.p = standard error of proportions  

  #Plotting in descending rank order of seed maturation success by family with AMBIENT as base
  seed_fam2 %>%
    group_by(famID) %>%
    mutate(seed_by_fam2 = if_else(treatment=="A", 1, 0)*seed_by_fam) %>% 
    ungroup %>%
    ggplot(aes(x=reorder(famID, desc(seed_by_fam2)), y=seed_by_fam, group=treatment, color=treatment)) +
    geom_point(size=1.5) +
    geom_hline(yintercept=seed_A, color="steelblue2") +
    geom_hline(yintercept=seed_H, color="tomato2") +
    geom_errorbar(aes(x=famID, ymin=(seed_by_fam-se.p), ymax=(seed_by_fam+se.p), width=0.5)) + 
    geom_smooth(se=F) +
    scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
    labs(x="Sibship", y="Average Seed Maturation Success")  +
    theme_bw() + theme(legend.position="none") + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
    theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

#Graphs - Trait distributions by FAMILY/SIBSHIP for CONTINUOUS traits (leaf, flwr_clstr, seed_pods) ####

dat %>% #leaf number
  group_by(treatment, famID, matID) %>%
  summarise(leaf_avg=mean(leaf), n=n(), sd=sd(leaf), se=(sd/(sqrt(n)))) %>%
  mutate(leaf_avg2 = if_else(treatment=="A", 1, 0)*leaf_avg) %>% 
  ungroup %>%
  ggplot(aes(x=reorder(famID, desc(leaf_avg2)), y=leaf_avg, 
             group=treatment, color=treatment)) +
  geom_errorbar(aes(ymin=leaf_avg-se, ymax=leaf_avg+se), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  geom_smooth(se=F) + 
  geom_hline(yintercept=leaf_A, color="steelblue2") +
  geom_hline(yintercept=leaf_H, color="tomato2") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Sibship", y="Average Leaf Number")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))


dat %>% #flower clusters
  group_by(treatment, famID, matID) %>%
  summarise(flwr_clstr_avg=mean(flwr_clstr), n=n(), sd=sd(flwr_clstr), se=(sd/(sqrt(n)))) %>%
  mutate(flwr_clstr_avg2 = if_else(treatment=="A", 1, 0)*flwr_clstr_avg) %>% 
  ungroup %>%
  ggplot(aes(x=reorder(famID, desc(flwr_clstr_avg2)), y=flwr_clstr_avg, 
             group=treatment, color=treatment)) + 
  geom_errorbar(aes(ymin=flwr_clstr_avg-se, ymax=flwr_clstr_avg+se), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  geom_smooth(se=F) +
  geom_hline(yintercept=flwrclstr_A, color="steelblue2") +
  geom_hline(yintercept=flwrclstr_H, color="tomato2") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Sibship", y="Average Flowering Clusters Number")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

dat %>% #seed pods
  group_by(treatment, famID, matID) %>%
  summarise(pods_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  mutate(pods_avg2 = if_else(treatment=="A", 1, 0)*pods_avg) %>% 
  ungroup %>%
  ggplot(aes(x = reorder(famID, desc(pods_avg2)), y = pods_avg,
             group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=pods_avg-se, ymax=pods_avg+se), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  geom_smooth(se=F) +
  geom_hline(yintercept=pods_A, color="steelblue2") +
  geom_hline(yintercept=pods_H, color="tomato2") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Sibship", y="Average Seed Pod Number")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

dat %>% #seed pods, excluding individuals that did NOT produce seed pods
  group_by(treatment, famID, matID) %>%
  filter(seed==1) %>%
  summarise(pods_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  mutate(pods_avg2 = if_else(treatment=="A", 1, 0)*pods_avg) %>% 
  ungroup %>%
  ggplot(aes(x = reorder(famID, desc(pods_avg2)), y = pods_avg,
             group = treatment, color = treatment)) + 
  geom_errorbar(aes(ymin=pods_avg-se, ymax=pods_avg+se), width=0.2, size=0.5) + 
  geom_point(aes(color=treatment)) + 
  geom_smooth(se=F) +
  geom_hline(yintercept=pods_A2, color="steelblue2") +
  geom_hline(yintercept=pods_H2, color="tomato2") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Sibship", y="Average Seed Pod Number")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

#Graphs - Maternal Correlations for BINARY traits (germ, flower, seed) ####
germ_mat <- dat %>% #germination
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
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  ggtitle(label="A)") +
  labs(x="Average Germination Success - Maternal Subgroup A", y="Average Germination Success - Maternal Subgroup B")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

flwr_mat <- dat %>% #flowering success
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
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  ggtitle(label="B)") +
  labs(x="Average Flowering Success - Maternal Subgroup A", y="Average Flowering Success - Maternal Subgroup B")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


seed_mat <- dat %>% #seed maturation success
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
  geom_point(size=1.5) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=seed_by_fam_A, ymin=(seed_by_fam_B-se.p_B), ymax=(seed_by_fam_B+se.p_B), width=0.01)) +
  geom_errorbarh(aes(y=seed_by_fam_B, xmin=(seed_by_fam_A-se.p_A), xmax=(seed_by_fam_A+se.p_A))) +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,1)) +
  ggtitle(label="C)") +
  labs(x="Average Seed Maturation Success - Maternal Subgroup A", y="Average Seed Maturation Success - Maternal Subgroup B")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Graphs - Maternal Correlations for CONTINUOUS traits (leaf, flwr_clstr, seed_pods) ####

leaf_mat <- dat %>% #leaf number
  group_by(treatment, famID, matID) %>%
  summarise(leaf_avg=mean(leaf), n=n(), sd=sd(leaf), se=(sd/(sqrt(n)))) %>%
  select(treatment, famID, matID, leaf_avg, se) %>%
  ungroup %>%
  group_by(treatment, famID) %>%
  mutate(matID = c("A", "B")[sequence(n())]) %>%
  pivot_wider(names_from = matID, values_from = c(leaf_avg, se))

leaf_mat %>% #Note: 12 removed values because no paired maternal group
  ggplot(aes(x=leaf_avg_A, y=leaf_avg_B, group=treatment, color=treatment)) +
  geom_point(size=1.5) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=leaf_avg_A, ymin=(leaf_avg_B-se_B), ymax=(leaf_avg_B+se_B), width=0.01)) +
  geom_errorbarh(aes(y=leaf_avg_B, xmin=(leaf_avg_A-se_A), xmax=(leaf_avg_A+se_A))) +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_y_continuous(limits=c(0,10)) +
  scale_x_continuous(limits=c(0,10)) +
  labs(x="Average Leaf Number - Maternal Subgroup A", y="Average Seed Maturation Success - Maternal Subgroup B")  +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


flwrclstr_mat <- dat %>% #flower clusters
  group_by(treatment, famID, matID) %>%
  summarise(flwr_clstr_avg=mean(flwr_clstr), n=n(), sd=sd(flwr_clstr), se=(sd/(sqrt(n)))) %>%
  select(treatment, famID, matID, flwr_clstr_avg, se) %>%
  ungroup %>%
  group_by(treatment, famID) %>%
  mutate(matID = c("A", "B")[sequence(n())]) %>%
  pivot_wider(names_from = matID, values_from = c(flwr_clstr_avg, se))

flwrclstr_mat %>% #Note: 12 removed values because no paired maternal group
  ggplot(aes(x=flwr_clstr_avg_A, y=flwr_clstr_avg_B, group=treatment, color=treatment)) +
  geom_point(size=1.5) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=flwr_clstr_avg_A, ymin=(flwr_clstr_avg_B-se_B), ymax=(flwr_clstr_avg_B+se_B), width=0.01)) +
  geom_errorbarh(aes(y=flwr_clstr_avg_B, xmin=(flwr_clstr_avg_A-se_A), xmax=(flwr_clstr_avg_A+se_A))) +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Average Flower Cluster Number - Maternal Subgroup A", y="Average Flower Cluster Number - Maternal Subgroup B")  +
  ggtitle("A)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

pods_mat <- dat %>% #seed pods
  group_by(treatment, famID, matID) %>%
  summarise(pods_avg=mean(seed_pods), n=n(), sd=sd(seed_pods), se=(sd/(sqrt(n)))) %>%
  select(treatment, famID, matID, pods_avg, se) %>%
  ungroup %>%
  group_by(treatment, famID) %>%
  mutate(matID = c("A", "B")[sequence(n())]) %>%
  pivot_wider(names_from = matID, values_from = c(pods_avg, se))

pods_mat %>% #Note: 12 removed values because no paired maternal group
  ggplot(aes(x=pods_avg_A, y=pods_avg_B, group=treatment, color=treatment)) +
  geom_point(size=1.5) +
  geom_abline(slope=1, intercept=0, size=1.2, linetype="dashed") +
  geom_errorbar(aes(x=pods_avg_A, ymin=(pods_avg_B-se_B), ymax=(pods_avg_B+se_B), width=0.01)) +
  geom_errorbarh(aes(y=pods_avg_B, xmin=(pods_avg_A-se_A), xmax=(pods_avg_A+se_A))) +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Average Seed Pod Number - Maternal Subgroup A", y="Average Seed Pod Number - Maternal Subgroup B")  +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


#Graphs - Genetic Correlations between Select Traits and Seed Pod Number (Fitness) ####

#Seed Pod Number against Germination Success - Family Level ONLY
dat %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(germ), n=n(), pods_avg=mean(seed_pods)) %>%
  mutate(germ_avg=(sum/n)) %>%
  ggplot(aes(x=germ_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Average Germination Success", y="Average Seed Pods") +
  ggtitle("A)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Flowering Success - FAMILY LEVEL ONLY
dat %>%
  group_by(treatment, famID) %>%
  summarise(sum=sum(flower), n=n(), pods_avg=mean(seed_pods)) %>%
  mutate(flwr_avg=(sum/n)) %>%
  ggplot(aes(x=flwr_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Average Flowering Success", y="Average Seed Pods") +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


#Seed Pod Number against Leaf Number - Individual Level
dat %>%
  filter(!germ==0) %>%
  ggplot(aes(x=leaf, y=seed_pods, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of Leaves", y="Number of Seed Pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Leaf Number - Family Level
dat %>%
  filter(!germ==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), leaf_avg=mean(leaf)) %>%
  ggplot(aes(x=leaf_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Average Number of Leaves", y="Average Seed Pods") +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Budding Branches - Individual Level
dat %>%
  filter(!germ==0, !bud_clstr==0) %>%
  ggplot(aes(x=bud_clstr, y=seed_pods, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of Budding Branches", y="Number of Seed Pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Budding Branches - Family Level
dat %>%
  filter(!germ==0, !bud_clstr==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), buds_avg=mean(bud_clstr)) %>%
  ggplot(aes(x=buds_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Average Budding Branches", y="Average Seed Pods") +
  ggtitle("B)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


#Seed Pod Number against Flowering Clusters - Individual Level
dat %>%
  filter(!germ==0, !flower==0) %>%
  ggplot(aes(x=flwr_clstr, y=seed_pods, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of Flowering Clusters", y="Number of Seed Pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Flowering Clusters - Family Level
dat %>%
  filter(!germ==0, !flower==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), flwrclstr_avg=mean(flwr_clstr)) %>%
  ggplot(aes(x=flwrclstr_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Average Flowering Clusters", y="Average Seed Pods") +
  ggtitle("B)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Height - Individual Level
dat %>%
  filter(!germ==0, !height==0) %>%
  ggplot(aes(x=height, y=seed_pods, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Height (cm)", y="Number of Seed Pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Height - Family Level
dat %>%
  filter(!germ==0, !height==0) %>% #filtering out height==0 because those were originally NAs
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), height_avg=mean(height)) %>%
  ggplot(aes(x=height_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Average Height (cm)", y="Average Seed Pods") +
  ggtitle("B)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Stem Diameter - Individual Level
dat %>%
  filter(!germ==0, !stem_diam==0) %>%
  ggplot(aes(x=stem_diam, y=seed_pods, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Stem Diameter (mm)", y="Number of Seed Pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Stem Diameter - Family Level
dat %>%
  filter(!germ==0, !stem_diam==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), stem_avg=mean(stem_diam)) %>%
  ggplot(aes(x=stem_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Average Stem Diameter (mm)", y="Average Seed Pods") +
  ggtitle("B)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


#Graphs - Life-time Fitness Summary by Proportion of Total Plants ####
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
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 1, 0.1)) +
  scale_x_discrete(limits=c("Planted", "Germinated", "Flowered", "Seeded")) +
  labs(x="Life History Event", y="Relative Proportion") + 
  ggtitle("A)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Plot seed pods for each treatment
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
  labs(x="Treatment", y="Mean Family Absolute Fitness") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 40, 10)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#GxE
fam.lifetime %>%
  ggplot(aes(x=treatment, y=fam.mean, color=treatment)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.5) +
  stat_summary(fun.y=mean, geom="point", shape=23, size=6,fill="darkgreen",alpha=0.5) +
  labs(x="Treatment", y="Mean Family Absolute Fitness") +
  scale_x_discrete(labels=c("Ambient", "Heated")) +
  scale_y_continuous(breaks=seq(0, 40, 10)) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
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
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  scale_x_continuous(limits=c(0,300)) + #The range is HIGH! From 0 to 300.. not going to show this graph
  labs(x="Total Seed Pods per Plant", y="Frequency")  + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

####Graphs - Density plots using FAMILY MEANS for fitness traits####
mu.lifetime <- dat %>% #reused from above
  group_by(treatment) %>%
  summarise(grp.mean=mean(seed_pods))
fam.lifetime <- dat %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(seed_pods))
fam.lifetime %>% #lifetime fitness (seed pods)
  ggplot(aes(x=fam.mean, fill=treatment, color=treatment)) +
  geom_density(alpha=0.3) +
  geom_vline(data=mu.lifetime, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  scale_x_continuous(limits=c(0,40)) +
  labs(x="Mean Family Absolute Fitness (W)", y="Frequency")  + 
  ggtitle("C)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

germ_fam <- dat %>% #reused from above
  group_by(treatment, famID) %>%
  summarise(sum=sum(germ), n=n()) %>%
  mutate(germ_by_fam=(sum/n), se.p=(sqrt(germ_by_fam*(1-germ_by_fam))/n)) #se.p = standard error of proportions
germ_summary <- dat %>%
  group_by(treatment) %>%
  summarise(sum=sum(germ), n=n(), germ_avg=sum/n)
germ_fam %>% #germination success of all plants
  ggplot(aes(x=germ_by_fam, fill=treatment, color=treatment)) +
  geom_density(alpha=0.3) +
  geom_vline(data=germ_summary, aes(xintercept=germ_avg, color=treatment), linetype="dashed", size=1) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  scale_x_continuous(limits=c(0,1)) +
  labs(x="Mean Family Survivorship", y="Frequency")  + 
  ggtitle("A)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

mu.pods <- dat %>% #reused from above 
  filter(!germ==0, !flower==0) %>%
  group_by(treatment) %>%
  summarise(grp.mean=mean(seed_pods))
fam.pods <- dat %>%
  filter(!germ==0, !flower==0) %>%
  group_by(treatment, famID) %>%
  summarise(fam.mean=mean(seed_pods))
fam.pods %>% #seed pods
  ggplot(aes(x = fam.mean, fill = treatment, color = treatment)) + 
  geom_density(alpha=0.3) +
  geom_vline(data=mu.pods, aes(xintercept=grp.mean, color=treatment), linetype="dashed", size=1) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  scale_fill_manual(values=c("dodgerblue2", "tomato2")) +
  scale_x_continuous(limits=c(0,65)) +
  labs(x="Mean Family Fecundity", y="Frequency")  +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm")) +
  theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))

#GLM + GLMM testing mean differences caused by TREATMENT ONLY for Continuous Traits ####
#Using package "Car" #Random Effect: plot | Fixed Effects: treatment + famID
#Only height and stem diameter can be evaluted. Flwr_clstr, seed_pods, and leaf are poisson distributed and to include random effects, must be analyzed using MCMCglmm. 

#NOTE: flowering time and first flower analysis are outdated and have been replaced by survival Cox analyses
#Time to First Flower using package lme4
library(lme4)
dat1 <- dat %>% 
  filter(!flwr_juln_date==0, !flwr_compl_juln==0, !flower==0) %>%
  mutate(first_flower=(flwr_juln_date - germ_juln_date))
hist(dat1$first_flower) #fairly normal (note Shapiro-Wilk not good for dealing with big sample sizes)
qqnorm(dat1$first_flower)
qqline(dat1$first_flower)
summary(model1a <- lmer(dat=dat1, first_flower ~ treatment + famID + (1|plot), REML=FALSE))
summary(model1b <- lmer(dat=dat1, first_flower ~ treatment*famID + (1|plot), REML=FALSE))
anova(model1a, model1b) #model1a  lower AIC/BIC value and selected
summary(model1a)
Anova(model1a)

#Flowering Season Length using package nlme
library(lme4)
library(MASS)
dat2 <- dat %>% 
  filter(!flwr_juln_date==0, !flwr_compl_juln==0, !flower==0) %>%
  mutate(flowering_time=(flwr_compl_juln - flwr_juln_date)) %>%
  filter(!flowering_time==0)
dat2 %>% 
  ggplot(aes(x=flowering_time, fill=treatment)) + 
  geom_histogram(bins=30, position="dodge") #tri-modal distribution
  #It's trimodal because census were taken every week.. so flowering_time could follow a different distribution if you account for days in between censuses.
qqnorm(dat2$flowering_time)
qqline(dat2$flowering_time) #non-normal distribution
model2 <- glmmPQL(dat=dat2, flowering_time ~ treatment + famID, ~1 | plot, family = gaussian(link = ""))

#Height
dat3 <- dat %>% filter(!germ==0, !height==0) #Removing 0s because previously NAs
hist(dat3$height) #fairly normal (note Shapiro-Wilk not good for dealing with big sample sizes)
qqnorm(dat3$height)
qqline(dat3$height)
model3a <- glm(dat=dat3, height ~ treatment) #AIC = 23126
model3b <- glm(dat=dat3, height ~ treatment + famID) #AIC = 22945
model3c <- glm(dat=dat3, height ~ treatment * famID) #AIC = 23019
model3d <- lmer(dat=dat3, height ~ treatment + famID + (1|plot), REML=FALSE) #AIC 22526 ***
model3e <- lmer(dat=dat3, height ~ treatment * famID + (1|plot), REML=FALSE) #AIC = 22603
anova(model3a, model3b) #model1a lower AIC/BIC value and selected
summary(model3d)
Anova(model3d) #for p values


#Stem Diameter
library(fitdistrplus) #using this package to discover the distribution type
dat4 <- dat %>% filter(!germ==0, !stem_diam==0) #Removing 0s because previously NAs
plotdist(dat4$stem_diam, histo = TRUE, demp = TRUE)
descdist(dat4$stem_diam, boot=1000) #potentially gamma or exponential distribution
dist1a <- fitdist(dat4$stem_diam, "gamma")
dist1b <- fitdist(dat4$stem_diam, "lnorm")
par(mfrow = c(2, 2))
plot.legend <- c("gamma", "lognormal")
denscomp(list(dist1a, dist1b), legendtext = plot.legend)
qqcomp(list(dist1a, dist1b), legendtext = plot.legend)
cdfcomp(list(dist1a, dist1b), legendtext = plot.legend)
ppcomp(list(dist1a, dist1b), legendtext = plot.legend) #It's not great, but I can fit the log-normal distribution.
par(mfrow= c(1,1))

  #Using PQLglmm
model4a <- glm(data=dat4, log(stem_diam) ~ treatment) #AIC = 5122.7
model4b <- glm(data=dat4, log(stem_diam) ~ treatment + famID) ## Preferred | AIC = 5036.2
model4c <- glm(data=dat4, log(stem_diam) ~ treatment * famID) #AIC = 5105.6
anova(model4a, model4b, model4c)
library(nlme)
model4d <- glmmPQL(data=dat4, stem_diam ~ treatment + famID , ~1 | plot, family = gaussian(link="log"), verbose=FALSE)
summary(model4d)

#Leaf Number
dat5 <- dat %>% filter(!germ==0)
hist(dat5$leaf) #fairly normal (note Shapiro-Wilk not good for dealing with big sample sizes)
qqnorm(dat3$leaf)
qqline(dat3$leaf)
model5a <- glm(dat=dat5, leaf ~ treatment) #AIC = 14815
model5b <- glm(dat=dat5, leaf ~ treatment + famID) #AIC = 14594 **
model5c <- glm(dat=dat5, leaf ~ treatment*famID) #AIC = 14650
model5d <- lmer(dat=dat5, leaf ~ treatment + famID + (1|plot), REML=FALSE) #AIC = 14541 *****
model5e <- lmer(dat=dat5, leaf ~ treatment*famID + (1|plot), REML=FALSE)
anova(model5a, model5b, model5c) #model with treatment + famID preferred
anova(model5d, model5e) #model5d preferred with random effects
summary(model5d)
Anova(model5d)


#LMM (Frequentist) and MCMCglmm (Bayesian) models for Poisson distributed Traits ####


  #Seed Pod Number
dat6 <- dat %>% filter(!germ==0, !flower==0) #excluding plants didn't germinate and flower BECAUSE want to understand differences in seed pod production between treatments
model6a <- glm(seed_pods ~ treatment, family="poisson", data=dat6)
model6b <- glm(seed_pods ~ treatment + famID, family="poisson", data=dat6) #AIC value much lower in model6b.. AIC = 65472
model6c <- glm(seed_pods ~ treatment * famID, family="poisson", data=dat6) #AIC = 64319
model6d <- glmer(seed_pods ~ treatment + famID + (1|plot), family="poisson", data=dat6) #AIC much lower in model 6e.. AIC = 61075.3
model6e <- glmer(seed_pods ~ treatment * famID + (1|plot), family="poisson", data=dat6) #AIC = 60015
model6f <- hurdle(seed_pods ~ treatment * famID | flwr_clstr, dist="poisson", link="logit", data=dat6) #using pscl
summary(model6e) #Laplace approximation *** used this because Poisson distributed and less zero inflated
summary(model6f) #hurdle model 
vuong(model6c, model6f)

  #Life time fitness measured as Seed Pod Number
dat7 <- dat
hist(dat7$seed_pods) #Highly zero inflated poisson distribution
qqnorm(dat7$seed_pods)
qqline(dat7$seed_pods)
plotdist(dat7$seed_pods, histo = TRUE, demp = TRUE)
descdist(dat7$seed_pods, boot=1000) #potentially gamma or exponential distribution
dist7a <- fitdist(dat7$seed_pods, "pois")
summary(dist7a)
plotdist(dist7a)
qqcomp(dist7a)
cdfcomp(dist7a)

library(pscl) #LIFETIME FITNESS
model7a <- glmer(seed_pods ~ treatment * famID + (1|plot), family="poisson", data=dat7) #using lme4
  #AIC = 137385.9, treatmentH = 0.7193, SE = 0.131, z = 5.465, P< 2e-16
model7b <- hurdle(seed_pods ~ treatment * famID | germ + flwr_clstr, dist="poisson", link="logit", data=dat7) #using zero inflated
model7c <- hurdle(seed_pods ~ treatment + famID | germ + flwr_clstr, dist="poisson", link="logit", data=dat7) #using zero inflated
model7d <- glm(seed_pods ~ treatment * famID, family="poisson", data=dat7) #normal Poisson model
  #treatmentH = 0.494, se = 0.0609, z = 8.121, P= 4.64e-16
  #zero hurdle model has est 22.38, se=484, z=0.046, P=0.963
summary(model7a) #Laplace approximation
summary(model7b) #hurdle model **** used this
vuong(model7d, model7b) #model comparison: model with hurdle is more significant
vuong(model7b, model7c) #model with interaction is significant


  #Flowering Clusters number
dat8 <- dat %>% filter(!germ==0) #excluding plants didn't germinate BECAUSE want to understand differences in seed pod production between treatments
model8a <- glm(flwr_clstr ~ treatment, family="poisson", data=dat8) #AIC = 13059
model8b <- glm(flwr_clstr ~ treatment + famID, family="poisson", data=dat8) #AIC = 12574
model8c <- glm(flwr_clstr ~ treatment * famID, family="poisson", data=dat8) #AIC = 12567
model8d <- glmer(flwr_clstr ~ treatment + famID + (1|plot), family="poisson", data=dat8) #AIC 12789
model8e <- glmer(flwr_clstr ~ treatment * famID + (1|plot), family="poisson", data=dat8) #FAILED TO CONVERGE
model8f <- hurdle(flwr_clstr ~ treatment * famID | flwr_clstr, dist="poisson", link="logit", data=dat8) #using zero inflated
summary(model8f) #hurdle model
summary(model8c)
vuong(model8c, model8f) #although hurdle model better, used Poisson because its a subsample and not that zero-inflated

  #germination success
dat9 <- dat
model9a <- glm(germ ~ treatment, family = "binomial", data=dat9) #AIC = 9403.5
model9b <- glm(germ ~ treatment + famID, family = "binomial", data=dat9) #AIC = 9333.3 ***
model9c <- glm(germ ~ treatment * famID, family = "binomial", data=dat9) #AIC = 9385
model9d <- glmer(germ ~ treatment + famID + (1|plot), family = "binomial", data=dat9) #failed to converge
model9e <- glmer(germ ~ treatment * famID + (1|plot), family = "binomial", data=dat9) #failed to converge
anova(model9b, model9a)
summary(model9b)

  #flowering success - reusing dat8
model10a <- glm(flower ~ treatment, family = "binomial", data=dat8) #AIC = 1788.7
model10b <- glm(flower ~ treatment + famID, family = "binomial", data=dat8) #AIC = 1818.5
model10c <- glm(flower ~ treatment * famID, family = "binomial", data=dat8) #AIC = 1859.5
model10d <- glmer(flower ~ treatment + famID + (1|plot), family = "binomial", data=dat8) #failed to converge
model10e <- glmer(flower ~ treatment * famID + (1|plot), family = "binomial", data=dat8) #failed to converge
anova(model10a, model10b, model10c)
summary(model10a)

  #seed maturation success - reusing dat6
model11a <- glm(seed ~ treatment, family = "binomial", data=dat6) #AIC = 201.14 ***
model11b <- glm(seed ~ treatment + famID, family = "binomial", data=dat6) #AIC = 269.67
model11c <- glm(seed ~ treatment * famID, family = "binomial", data=dat6) #373.7
model11d <- glmer(seed ~ treatment + famID + (1|plot), family = "binomial", data=dat6) #failed to converge
model11e <- glmer(seed ~ treatment * famID + (1|plot), family = "binomial", data=dat6) #failed to converge
summary(model11a)



## MCMCglmm ... not used
library(MCMCglmm)
prior6 = list(R = list(V = 1, n = 0, fix = 1), G = list(G1 = list(V = 1, n = 1),
              G2 = list(V = 1, n = 1), G3 = list(V = 1, n = 1), G4 = list(V = 1, n = 1),
              G5 = list(V = 1, n = 1)))
set.seed(45)
MCMC <- MCMCglmm(profit ~ 1, random = ~year + farmer + place + gen + district,
                 data = farmers, family = "categorical", prior = prior, verbose = FALSE)
summary(MCMC)

#### Calculating Raw Estimate of Adaptive Potential ####

#Ambient Treatment
dat %>% 
  group_by(treatment, famID) %>%
  summarise(fitness.fam=mean(seed_pods)) %>%
  ungroup() %>%
  group_by(treatment) %>%
  summarise(fitness.trt=mean(fitness.fam), variance=var(fitness.fam), adapt.pot=variance/fitness.trt)
  



.########################################################################.
 ################      DATA PREP FOR ASTER ANALYSIS      ################
.########################################################################.

#Selecting variables needed for aster: posID, individual, plot, matID, patID, famID, treatment, germ, flower, seed, leaf, flwr_clstr, seed_pods
aster.2019<- dat %>%
  select(posID, individual, plot, matID, patID, famID, treatment, germ, flower, seed, leaf, flwr_clstr, seed_pods)

#Checking for number of plants with aborted flowers
nrow(subset(aster.2019, flower==1 & seed_pods==0)) #182 plants flowered but did not produce seed pods

#Writing combined dataframe into a .csv file
write.csv(aster.2019,"Rdata/aster_2019_cleaned.csv",row.names = FALSE)

  
.########################################################################.
 ###################      ADDITIONAL SPATIAL MAPS     ###################
.########################################################################.

#Germination heatmap using base R
#WARNING: the heat maps do not work because after removing the doubles and accidentally damaged plants, the total number of individuals in the study is 549 and does not fit into the grid of a 20x30 space (unfilled spaces). Thus, there is an error created. 
p1_Germ <- dat %>% #p1 stands for plot 1
  filter(plot==1, posID<601) %>% #Note that the extra 5+ plants were omitted
  select(germ)
p1_Germ <- matrix(p1_Germ$germ, nrow=20, ncol=30, byrow=TRUE)
heatmap(p1_Germ, Rowv=NA, Colv=NA) #note that heatmaps are NOT useful for binary data, such as germination data. It's better for count data
heatmap.2(p1_Germ, Rowv=FALSE, Colv=FALSE, dendrogram='none', cellnote=p1_Germ, notecol="black", trace='none', key=FALSE, lwid = c(.01,.99),lhei = c(.01,.99), margins = c(4,4)) #this uses package gplots

#Running loops for plots 1 to 12
list.germ <- data %>% #this method splits the dataset into a list
  filter(posID < 601) %>%
  group_split(plot) %>%
  map(. %>% select(germ_bin))

unique_plot <- unique(dataset$plot)  #this method uses for loops to create a list
plot_list <- list(length = length(unique_plot))
for(i in seq_along(unique_plot)) {
  plot_list[[i]] <- dataset %>%
    filter(plot == unique_plot[i], pos_ID<601) %>%
    select(germ_bin)
}

lapply(split(dataset, data$plot), function(x) #this method does not use a list
  subset(x, pos_ID < 601, select = germ_bin, drop = FALSE))
  
#Germination heatmap using ggplot2 - Need to provide x and y values to each
ggplot(plot1, aes(X, Y, z=Z)) + geom_tile(aes(fill=Z)) + theme_bw() #this isn't quite complete because you need to assign x and y values for each plant individual.. which I will do late



.########################################################################.
 #######################      ADDED FUNCTIONS     #######################
.########################################################################.


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}