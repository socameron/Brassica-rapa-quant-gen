#### PROJECT: Brassica rapa Va/W Study (Data collected at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Clean data prior to MCMCglmm analysis - Note: this code is already provided in 01_Data_exploration

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Installing and Loading Necessary Packages
library(gplots)
library(tidyverse)

#Importing Data using readr from tidyverse (Base R confuses the class for certain vectors). Excluding notes from the data import
# d for double, f for factor, i for integer
col_types_list <- cols_only(posID = "d", individual = "d", plot = "f",
			                      animal = "f", matID = "f", patID = "f", famID = "f",
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
                            height = "d", stem_diam = "d",
			                      posX = "d", posY = "d",
			                      germ_raster = "d", germ_neigh = "d")
data <- read_csv("Rdata/heatarrays_brassica_2019_data.csv", 
                 col_types=col_types_list, col_names = TRUE, na = "NA")
spec(data)
data <- as.data.frame(data)
lapply(data, class)
lapply(data, levels)

#Note: I am not using the calendar date data. However, if you wish to use that data, use package readr from tidyverse and function read_csv to parse the date columns correctly. 

#'####################################################################'#
######################      CLEANING DATASET      ######################
#'####################################################################'#

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
data[is.na(data)] <- 0

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
  select(posID, individual, plot, animal, matID, patID, famID, mat_group, treatment,
         germ_juln_date, flwr_juln_date, flwr_compl_juln,
         leaf, flwr_clstr, bud_clstr, germ, flower, seed, seed_pods, seed_dmg, seed_undmg, seed_number,
         height, stem_diam,
         germ_neigh)

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
head(subset(dat, seed==1&flwr_clstr==0)) # P1 133, P4 180, and P5 369 were still flowering. P7 170 did not record and flowering data, but seed pod data. Remove from analysis.
head(subset(dat, seed_pods==0&seed_number>0)) 

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
nrow(subset(dat, seed_pods==0&seed_number>0)) # 16 errors remain - see below for removal in another dataset file. 

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

#'####################################################################'#
##############      DATA PREP FOR MCMCglmm ANALYSIS      ###############
#'####################################################################'#

#Selecting variables needed for MCMCglmm:
MCMC.2019 <- dat %>%
  select(posID, individual, animal, plot, matID, patID, famID, mat_group, treatment, germ, flower, seed, leaf, flwr_clstr, seed_pods, seed_number, height, stem_diam, germ_census, flwr_census)

#Checking for number of plants with aborted flowers
nrow(subset(MCMC.2019, flower==1 & seed_pods==0)) #182 plants flowered but did not produce seed pods

#Writing combined dataframe into a .csv file
write.csv(MCMC.2019,"Rdata/MCMC_2019_cleaned.csv", row.names = FALSE)

