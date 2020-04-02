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
library(ggplot2)
library(dplyr)
library(gplots)
library(tidyverse)

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
dat %>% #germination
  group_by(treatment) %>%
  summarise(germ=sum(germ), n=n()) %>%
  mutate(germ_perc=(100*germ/n))

dat %>% #flowering success
  group_by(treatment) %>%
  filter(!germ==0) %>%
  summarise(flower=sum(flower), n=n()) %>%
  mutate(flwr_perc=(100*flower/n))

dat %>% #seed pod maturation
  group_by(treatment) %>%
  filter(!germ==0, !flower==0) %>%
  summarise(seed=sum(seed), n=n()) %>%
  mutate(seed_perc=(100*seed/n))


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

#Tally for Bird Damage ####
dat %>% 
  filter(!germ==0, !flower==0, !seed==0) %>%
  group_by(treatment, famID) %>%
  summarise(bird_plants=sum(seed_dmg>0), undmg_plants=sum(seed_undmg>=0 & seed_dmg==0), total_plants=sum(seed_pods>=0))

dat %>%
  filter(!germ==0, !flower==0, !seed==0) %>%
  group_by(treatment) %>%
  summarise(seeded_plants=sum(seed==1))

.########################################################################.
################      DATA PREP FOR ASTER ANALYSIS      ################
.########################################################################.

#Selecting variables needed for aster: posID, individual, plot, matID, patID, famID, treatment, germ, flower, seed, leaf, flwr_clstr, seed_pods
aster.2019<- dat %>%
  select(posID, individual, plot, animal, matID, patID, famID, treatment, germ, flower, seed, leaf, flwr_clstr, seed_pods)

#Checking for number of plants with aborted flowers
nrow(subset(aster.2019, flower==1 & seed_pods==0)) #182 plants flowered but did not produce seed pods

#Writing combined dataframe into a .csv file
write.csv(aster.2019,"Rdata/aster_2019_cleaned.csv",row.names = FALSE)
