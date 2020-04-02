#### PROJECT: Brassica rapa GxE Study (Data collected by Cameron So 2019 at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Flowering Time Analysis
#### AUTHOR: Cameron So
#### DATE LAST MODIFIED: 2020/02/02

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

#Installing and Loading Necessary Packages####
library(ggplot2)
library(dplyr)
library(gplots)
library(tidyverse)
library(zoo)
library(car)
library(coxme)

col_types_list1 <- cols_only(posID = "d", individual = "d", plot = col_factor(levels=c(1:12)),
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
data1 <- read_csv("Rdata/heatarrays_brassica_2019_data.csv", col_names = TRUE, na = "NA", 
                 col_types=col_types_list1)
lapply(data1, class)

#From raw dataset, removing plant positions with double plants or accidental damage. 
data1 <- data1 %>%
  filter(is.na(double), is.na(damage)) 

#Removing families not to be included in the analysis
data1 <- data1 %>%
  filter(!packID=="M1S20-14", !packID=="M1S31-4", !packID=="M1S33-13", !packID=="M1S43-9", !packID=="M2S28-18")

#Converting NAs in data to 0
data1 <- data1 %>% replace(., is.na(.), "0")

#'####################################################################'#
  ############      CUMULATIVE PHENOLOGICAL FIGURES      #############
#'####################################################################'#

#Converting flowering juln dates to census number
#Note: data1a is for calculating cumulative germination, flowering and completed flowering over the season.
#Note: data_phen is for evaluating germination and flowering over time (NOT cumulative)
  #census 1-4
data1$flwr_juln_date[data1$flwr_juln_date == "119"] <- 1
data1$flwr_juln_date[data1$flwr_juln_date == "129"] <- 2
data1$flwr_juln_date[data1$flwr_juln_date == "134"] <- 3
data1$flwr_juln_date[data1$flwr_juln_date == "141"] <- 4
  #census 5
data1$flwr_juln_date[data1$flwr_juln_date == "147"] <- 5
data1$flwr_juln_date[data1$flwr_juln_date == "148"] <- 5
data1$flwr_juln_date[data1$flwr_juln_date == "149"] <- 5
  #census 6
data1$flwr_juln_date[data1$flwr_juln_date == "154"] <- 6
data1$flwr_juln_date[data1$flwr_juln_date == "155"] <- 6
data1$flwr_juln_date[data1$flwr_juln_date == "156"] <- 6
  #census 7
data1$flwr_juln_date[data1$flwr_juln_date == "162"] <- 7
data1$flwr_juln_date[data1$flwr_juln_date == "163"] <- 7
data1$flwr_juln_date[data1$flwr_juln_date == "165"] <- 7
  #census 8
data1$flwr_juln_date[data1$flwr_juln_date == "169"] <- 8
data1$flwr_juln_date[data1$flwr_juln_date == "170"] <- 8
data1$flwr_juln_date[data1$flwr_juln_date == "171"] <- 8
data1$flwr_juln_date[data1$flwr_juln_date == "172"] <- 8
  #census 9
data1$flwr_juln_date[data1$flwr_juln_date == "175"] <- 9
data1$flwr_juln_date[data1$flwr_juln_date == "176"] <- 9
data1$flwr_juln_date[data1$flwr_juln_date == "177"] <- 9
  #census 10
data1$flwr_juln_date[data1$flwr_juln_date == "183"] <- 10
data1$flwr_juln_date[data1$flwr_juln_date == "184"] <- 10
data1$flwr_juln_date[data1$flwr_juln_date == "185"] <- 10
  #census 11
data1$flwr_juln_date[data1$flwr_juln_date == "189"] <- 11
data1$flwr_juln_date[data1$flwr_juln_date == "190"] <- 11
data1$flwr_juln_date[data1$flwr_juln_date == "191"] <- 11
  #census 12
data1$flwr_juln_date[data1$flwr_juln_date == "199"] <- 12
data1$flwr_juln_date[data1$flwr_juln_date == "204"] <- 12

#Converting finished flowering juln dates to census number

#census 1-4
data1$flwr_compl_juln[data1$flwr_compl_juln == "119"] <- 1
data1$flwr_compl_juln[data1$flwr_compl_juln == "129"] <- 2
data1$flwr_compl_juln[data1$flwr_compl_juln == "134"] <- 3
data1$flwr_compl_juln[data1$flwr_compl_juln == "141"] <- 4
#census 5
data1$flwr_compl_juln[data1$flwr_compl_juln == "147"] <- 5
data1$flwr_compl_juln[data1$flwr_compl_juln == "148"] <- 5
data1$flwr_compl_juln[data1$flwr_compl_juln == "149"] <- 5
#census 6
data1$flwr_compl_juln[data1$flwr_compl_juln == "154"] <- 6
data1$flwr_compl_juln[data1$flwr_compl_juln == "155"] <- 6
data1$flwr_compl_juln[data1$flwr_compl_juln == "156"] <- 6
#census 7
data1$flwr_compl_juln[data1$flwr_compl_juln == "162"] <- 7
data1$flwr_compl_juln[data1$flwr_compl_juln == "163"] <- 7
data1$flwr_compl_juln[data1$flwr_compl_juln == "165"] <- 7
#census 8
data1$flwr_compl_juln[data1$flwr_compl_juln == "169"] <- 8
data1$flwr_compl_juln[data1$flwr_compl_juln == "170"] <- 8
data1$flwr_compl_juln[data1$flwr_compl_juln == "171"] <- 8
data1$flwr_compl_juln[data1$flwr_compl_juln == "172"] <- 8
#census 9
data1$flwr_compl_juln[data1$flwr_compl_juln == "175"] <- 9
data1$flwr_compl_juln[data1$flwr_compl_juln == "176"] <- 9
data1$flwr_compl_juln[data1$flwr_compl_juln == "177"] <- 9
#census 10
data1$flwr_compl_juln[data1$flwr_compl_juln == "183"] <- 10
data1$flwr_compl_juln[data1$flwr_compl_juln == "184"] <- 10
data1$flwr_compl_juln[data1$flwr_compl_juln == "185"] <- 10
#census 11
data1$flwr_compl_juln[data1$flwr_compl_juln == "189"] <- 11
data1$flwr_compl_juln[data1$flwr_compl_juln == "190"] <- 11
data1$flwr_compl_juln[data1$flwr_compl_juln == "191"] <- 11
#census 12
data1$flwr_compl_juln[data1$flwr_compl_juln == "199"] <- 12
data1$flwr_compl_juln[data1$flwr_compl_juln == "204"] <- 12

#Converting germination juln dates to census number

#census 1-4
data1$germ_juln_date[data1$germ_juln_date == "119"] <- 1
data1$germ_juln_date[data1$germ_juln_date == "129"] <- 2
data1$germ_juln_date[data1$germ_juln_date == "134"] <- 3
data1$germ_juln_date[data1$germ_juln_date == "141"] <- 4
#census 5
data1$germ_juln_date[data1$germ_juln_date == "147"] <- 5
data1$germ_juln_date[data1$germ_juln_date == "148"] <- 5
data1$germ_juln_date[data1$germ_juln_date == "149"] <- 5
#census 6
data1$germ_juln_date[data1$germ_juln_date == "154"] <- 6
data1$germ_juln_date[data1$germ_juln_date == "155"] <- 6
data1$germ_juln_date[data1$germ_juln_date == "156"] <- 6
#census 7
data1$germ_juln_date[data1$germ_juln_date == "162"] <- 7
data1$germ_juln_date[data1$germ_juln_date == "163"] <- 7
data1$germ_juln_date[data1$germ_juln_date == "165"] <- 7
#census 8
data1$germ_juln_date[data1$germ_juln_date == "169"] <- 8
data1$germ_juln_date[data1$germ_juln_date == "170"] <- 8
data1$germ_juln_date[data1$germ_juln_date == "171"] <- 8
data1$germ_juln_date[data1$germ_juln_date == "172"] <- 8
#census 9
data1$germ_juln_date[data1$germ_juln_date == "175"] <- 9
data1$germ_juln_date[data1$germ_juln_date == "176"] <- 9
data1$germ_juln_date[data1$germ_juln_date == "177"] <- 9
#census 10
data1$germ_juln_date[data1$germ_juln_date == "183"] <- 10
data1$germ_juln_date[data1$germ_juln_date == "184"] <- 10
data1$germ_juln_date[data1$germ_juln_date == "185"] <- 10
#census 11
data1$germ_juln_date[data1$germ_juln_date == "189"] <- 11
data1$germ_juln_date[data1$germ_juln_date == "190"] <- 11
data1$germ_juln_date[data1$germ_juln_date == "191"] <- 11
#census 12
data1$germ_juln_date[data1$germ_juln_date == "199"] <- 12
data1$germ_juln_date[data1$germ_juln_date == "204"] <- 12

#Removing individuals that have seed date but no flwr date, etc
nrow(subset(data1, seed_bin==1&flwr_bin==0)) # 3 errors
data1 <- data1[!(data1$seed_bin==1&data1$flwr_bin==0),]

nrow(subset(data1, flwr_juln_date==0&flwr_compl_juln>0)) # 2 errors
data1 <- data1[!(data1$flwr_juln_date==0&data1$flwr_compl_juln>0),]
data1b <- data1 

#### Preparation for Flowering Time Analysis ####

#Changing dataset into long format for flowering + germination by time analysis
data1 <- data1 %>%
  filter(!germ_bin==0) %>%
  select(posID, individual, plot, matID, famID, treatment, germ_bin, flwr_bin, seed_bin,
         germ_juln_date, flwr_juln_date, flwr_compl_juln) %>%
  mutate(flwr_juln_relative = (as.numeric(flwr_juln_date) - as.numeric(germ_juln_date)),
         flwr_compl_relative = (as.numeric(flwr_compl_juln) - as.numeric(germ_juln_date)))

nrow(subset(data1, flwr_juln_relative<0)) # 262 errors
data1$flwr_juln_relative[data1$flwr_juln_date == 0] <- 0
data1$flwr_compl_relative[data1$flwr_compl_juln == 0] <- 0

#Creating Long-format Data for Graphing Purposes
data1a <- data1 %>%
  gather(key="variable", value="census", germ_juln_date:flwr_compl_relative)

summary.data <- data1a %>%
  group_by(treatment, census) %>%
  summarise(first.flwr=sum(variable=="flwr_juln_relative"), #interchange with flwr_juln_relative or flwr_juln_date 
            last.flwr=sum(variable=="flwr_compl_relative"), #interchange
            germ.sum=sum(variable=="germ_juln_date"))

data1 %>% #data1 = ALL plants
  group_by(treatment) %>%
  summarise(n=n())

data1b %>% #data1b = ONLY GERMINATING PLANTS
  group_by(treatment) %>%
  summarise(n=n())

save(data1, file="Routput/Flowering_analysis_all_plants.RData")
save(data1b, file="Routput/Flowering_analysis_germ_plants.RData")

load(file="Routput/Flowering_analysis_all_plants.RData")
load(file="Routput/Flowering_analysis_germ_plants.RData")

#Note, there are some germinated plants that DID NOT germinate.
#These plants are represented by census = 0

#### Cumulative Percent Plants Germinating ####

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
#Values retrieved from the sum totals by Julian Calendar Day, and Total values obtained from 01_Data_exploration tallies

germ.cum <- germ.cum %>% group_by(treatment) %>% 
  mutate(cum.sum=cumsum(value), germ_cum_perc=100*(cum.sum/total))
lapply(germ.cum, class)

germ.cum %>%
  ggplot(aes(x=day, y=germ_cum_perc, group=treatment, color=treatment)) +
  geom_point(size=4) +
  geom_line(size=1.2) + 
  #scale_x_continuous(breaks=seq(1,12,1)) + #ONLY IF CENSUS
  scale_y_continuous(limits=c(0,100)) +
  ggtitle("A)") +
  labs(x="Julian Calendar Day", y="Cumulative % of Plants Germinated")  +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#### Cumulative Percent Plants Reaching Flowering ####

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

#Values retrieved from the sum totals by Julian Calendar Day, and Total values obtained from 01_Data_exploration tallies

first.flwr <- first.flwr %>% group_by(treatment) %>% 
  mutate(cum.sum=cumsum(value), flwr_cum_perc=100*(cum.sum/total))
lapply(first.flwr, class)

first.flwr %>%
  ggplot(aes(x=day, y=flwr_cum_perc, group=treatment, color=treatment)) +
  geom_point(size=4) +
  geom_line(size=1.2) + 
  #scale_x_continuous(breaks=seq(1,8,1)) + #ONLY IF CENSUS
  scale_y_continuous(limits=c(0,100)) +
  ggtitle("A)") +
  labs(x="Julian Calendar Day", y="Cumulative % of Plants Reached Flowering")  +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#### Cumulative Percent Plants Finishing Flowering ####

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

#Values retrieved from the sum totals by Julian Calendar Day, and Total values obtained from 01_Data_exploration tallies

last.flwr <- last.flwr %>% group_by(treatment) %>% 
  mutate(cum.sum=cumsum(value), flwr_cum_perc=100*(cum.sum/total))
lapply(last.flwr, class)

last.flwr %>%
  ggplot(aes(x=day, y=flwr_cum_perc, group=treatment, color=treatment)) +
  geom_point(size=4) +
  geom_line(size=1.2) + 
  #scale_x_continuous(breaks=seq(1,8,1)) + #ONLY IF CENSUS
  scale_y_continuous(limits=c(0,100)) +
  ggtitle("B)") +
  labs(x="Julian Calendar Day", y="Cumulative % of Plants Completed Flowering")  +
  scale_color_manual("Treatment", values=c("steelblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))










#### Cumulative Percent Plants Reaching Flowering - Date Relative to Germination ####

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

#Values retrieved from the sum totals by Julian Calendar Day, and Total values obtained from 01_Data_exploration tallies

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

#Values retrieved from the sum totals by Julian Calendar Day, and Total values obtained from 01_Data_exploration tallies

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
  ###############      PHENOLOGY OVER THE SEASON      ################
#'####################################################################'#

#Using tidyverse to convert dataframe into long format
data_phen<- data %>%
  select(posID, individual, plot, matID, patID, famID, mat_group, treatment,
         germ_bin, flwr_bin,
         germ_juln_date, flwr_juln_date, flwr_compl_juln, 
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

#Graphs 21 - Flower Cluster No. vs Time ####
#Note: dates are converted in census # as data was collected weekly
#Subtracting flwr_clstr from bud_clstr = flower cluster number per census since data was collected as a cumulative number 
load(file="Routput/heatarrays_data_phen.RData")

data_phen$census <- as.factor(data_phen$census)
data_phen %>%
  filter(flwr_bin==1, str_detect(phenotype, "flwr_clstr")) %>%
  group_by(treatment, census) %>%
  summarise(flwr_avg=sum(phen_num, na.rm=TRUE), n=n(), sd=sd(flwr_avg), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = census, y = flwr_avg, group = treatment, color = treatment)) + 
  #geom_errorbar(aes(ymin=pod_avg-se, ymax=pod_avg+se), width=0.5, size=0.5) + 
  geom_smooth(se=F) +
  ggtitle("Flowering Clusters Observed in Summer 2019") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Census Week", y="Flower Cluster Number")

#Graphs 22 - Growth Rate on Leaf No. vs Time ####
data_phen$census <- as.factor(data_phen$census)
data_phen %>%
  filter(germ_bin==1, str_detect(phenotype, "leaf")) %>%
  group_by(treatment, census) %>%
  summarise(leaf_avg=sum(phen_num, na.rm=TRUE), n=n(), sd=sd(leaf_avg), se=(sd/(sqrt(n)))) %>%
  ggplot(aes(x = census, y = leaf_avg, group = treatment, color = treatment)) + 
  #geom_errorbar(aes(ymin=pod_avg-se, ymax=pod_avg+se), width=0.5, size=0.5) + 
  geom_smooth(se=F) +
  ggtitle("Leaf Production Observed in Summer 2019") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels = c("Ambient", "Heated")) +
  labs(x="Census Week", y="True Leaf Number")


#### Failure Time Analysis AKA Survival Analysis (Cox) ####
library(survival)
library(survminer)
library(coxme)
view(eortc)

## We employ untruncated, right-censored Cox models including random plot effects

#Setting Appropriate Classes to Numeric/Double as classes were altered during cleanup
lapply(data1, class)
data1$posID <- as.factor(data1$posID)
data1$individual <- as.integer(data1$individual)
data1[,10:47] <- lapply(data1[,10:47], as.numeric) #setting all neccessary parameters to class numeric




## Flowering Schedule ##
#Unlike earlier graphs, I can make use of the flowering juln time without using census # as a dummy variable
#Treatment H = 1, A = 0 ... therefore comparison results are shown relative to the heated treatment
data1 %>%
  group_by(treatment) %>%
  summarise(mean_first_flwr=mean(flwr_juln_relative))

surv.time <- data1 %>% 
  filter(!germ_bin==0) %>%
  mutate(flowering_time=(as.numeric(flwr_compl_relative) - as.numeric(flwr_juln_relative))) %>%
  mutate(uncensor= ifelse(flowering_time <= 0, 0, 1)) %>%
  mutate(flowering_time = ifelse(flowering_time <= 0, 0, flowering_time)) %>%
  mutate(treatment = ifelse(treatment == "H", 1, 0)) %>%
  select(individual, plot, treatment, matID, famID, treatment, flowering_time, uncensor)

lapply(surv.time, class)
surv.time$treatment <- as.factor(surv.time$treatment)
#surv.time$famID <- as.numeric(surv.time$famID) #Treating sibship is numeric so get generalized results
#surv.time$plot <- as.numeric(surv.time$plot) #same for plot effects

#Note, uncensor has 0 = did not flower, 1 = flowered. 
cox.time1 <- coxph(Surv(flowering_time, uncensor, type="right") ~ treatment, data=surv.time)
cox.time2 <- coxph(Surv(flowering_time, uncensor, type="right") ~ treatment + famID, data=surv.time)
cox.time2a <- coxph(Surv(flowering_time, uncensor, type="right") ~ treatment * famID, data=surv.time)
cox.time3 <- coxme(Surv(flowering_time, uncensor, type="right") ~ treatment + famID + (1|plot), data=surv.time)
cox.time3a <- coxme(Surv(flowering_time, uncensor, type="right") ~ treatment * famID + (1|plot), data=surv.time)
cox.time3b <- coxph(Surv(flowering_time, uncensor, type="right") ~ treatment + famID + frailty(plot), data=surv.time) #Same as 3, just using frailty to compare differences in coxph vs coxme and use for diagnostics
anova(cox.time1, cox.time2)  #just famID not significant
anova(cox.time1, cox.time2a) #interaction with famID not significant
anova(cox.time1, cox.time3) #including random effects significant
anova(cox.time3, cox.time3a) #Interaction term not significant for random effect model
print(cox.time3)
summary(cox.time3a) #Heated treatment INCREASES flowering time by 47% (or a factor of exp(coef)= 0.53)
anova(cox.time3a)
#NOTE: Because this is a SURVIVAL analysis, it would have been REDUCED survival by 47%, but we treated flowering AS survival, so it flips.
confint(cox.time3a, level=0.95)
#Confidence Intervals: Lower 95% = 0.75, Upper 95% = 0.816, Hazard Ratio = 0.81

#Testing Cox assumptions
cox.zph(cox.time3a) #no significant for covariates -- GOOD .. using cox.time2 for covariate b/w treatment and famID.. but they shouldn't covary regardless because quite different values!
ggcoxdiagnostics(cox.time3a, type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_bw()) #GOOD..? except TRT???
lme(flowering_time ~ treatment + famID, random = ~1|plot, data=surv.time, method="ML") #NOT FINISHED - need to run diagnostics on residuals 

#Creating dataframe for plotting - not necessary though
cox.time2.fit <- survfit(Surv(flowering_time, uncensor, type="right") ~ treatment, data=surv.time)
ggsurvplot(cox.time2.fit, data=surv.time)




 ## First Flower ##
surv.first <- data1 %>% #time to first flower after germination
  filter(!germ_bin==0) %>%
  mutate(first_flower=(as.numeric(flwr_juln_date) - as.numeric(germ_juln_date))) %>%
  mutate(uncensor= ifelse(first_flower <= 0, 0, 1)) %>%
  mutate(first_flower = ifelse(first_flower <= 0, 0, first_flower)) %>%
  mutate(treatment = ifelse(treatment == "H", 1, 0)) %>%
  select(individual, plot, treatment, matID, famID, treatment, first_flower, uncensor)

lapply(surv.first, class)
surv.first$famID <- as.factor(surv.first$famID) #
surv.first$plot <- as.factor(surv.first$plot) #
surv.first$treatment <- as.factor(surv.first$treatment)

#Note, uncensor has 0 = did not flower, 1 = flowered. 
cox.first1 <- coxph(Surv(first_flower, uncensor, type="right") ~ treatment, data=surv.first)
cox.first2 <- coxph(Surv(first_flower, uncensor, type="right") ~ treatment + famID, data=surv.first)
cox.first2a <- coxph(Surv(first_flower, uncensor, type="right") ~ treatment * famID, data=surv.first)
cox.first3 <- coxme(Surv(first_flower, uncensor, type="right") ~ treatment + famID + (1|plot), data=surv.first)
cox.first3a <-coxme(Surv(first_flower, uncensor, type="right") ~ treatment*famID + (1|plot), data=surv.first)
cox.first3b <- coxph(Surv(first_flower, uncensor, type="right") ~ treatment + famID + frailty(plot), data=surv.first) #Same as 3, just using frailty to compare differences in coxph vs coxme and use for diagnostics
anova(cox.first1, cox.first2) #including famID important
anova(cox.first2, cox.first2a, cox.first3a) #interaction term not important
anova(cox.first2, cox.first3) #model 2 significant, random effect not neccessary.
summary(cox.first3a) #Heated treatment INCREASES flowering time by 21% (or a factor of exp(coef)= 0.79)
#NOTE: Because this is a SURVIVAL analysis, it would have been REDUCED survival by 21%, but we treated flowering AS survival, so it flips.
#Confidence Intervals: Lower 95% = 0.75, Upper 95% = 0.816, Hazard Ratio = 0.81
anova(cox.first3a)
#
















  





####OLD METHOD OF CALCULATING####
 
#Selecting columns needed for analysis (OLD)
data1 <- data1 %>%
  select(posID, individual, plot, matID, patID, famID, mat_group, treatment, 
       germ_juln_date, flwr_juln_date, flwr_compl_juln, germ_bin, flwr_bin, seed_bin,
       flwr_clstr1, flwr_clstr2, flwr_clstr3, flwr_clstr4, flwr_clstr5, flwr_clstr6, flwr_clstr7, flwr_clstr8)
#### Percent Plants Reaching Flowering ####

#Totals: A=1589, H=1576
data1 %>% group_by(treatment) %>% 
  filter(!germ_bin==0) %>%
  summarise(flower1 = sum(flwr_clstr1>0), n()) #A: 0 | H: 12

data1 %>% group_by(treatment) %>%
  filter(!germ_bin==0, !flwr_clstr1>0) %>%
  summarise(flower2 = sum(flwr_clstr2>0), n()) #A: 4 | H: 291

data1 %>% group_by(treatment) %>%
  filter(!germ_bin==0, !flwr_clstr1>0, !flwr_clstr2>0) %>%
  summarise(flower3 = sum(flwr_clstr3>0), n()) #A: 285 | H: 766

data1 %>% group_by(treatment) %>%
  filter(!germ_bin==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0) %>%
  summarise(flower4 = sum(flwr_clstr4>0), n()) #A: 539 | H: 289

data1 %>% group_by(treatment) %>%
  filter(!germ_bin==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0, !flwr_clstr4>0) %>%
  summarise(flower5 = sum(flwr_clstr5>0), n()) #A: 373 | H: 79

data1 %>% group_by(treatment) %>%
  filter(!germ_bin==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0, !flwr_clstr4>0, !flwr_clstr5>0) %>%
  summarise(flower6 = sum(flwr_clstr6>0), n()) #A: 177 | H: 33

data1 %>% group_by(treatment) %>%
  filter(!germ_bin==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0, !flwr_clstr4>0, !flwr_clstr5>0, !flwr_clstr6>0) %>%
  summarise(flower7 = sum(flwr_clstr7>0), n()) #A: 37 | H: 10

data1 %>% group_by(treatment) %>%
  filter(!germ_bin==0, !flwr_clstr1>0, !flwr_clstr2>0, !flwr_clstr3>0, !flwr_clstr4>0, !flwr_clstr5>0, !flwr_clstr6>0, !flwr_clstr7>0) %>%
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
