#### PROJECT: Climate analysis study (Data collected by KSR Weather Station at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Clean data and observe distributions over summer 2019
#### AUTHOR: Cameron So
#### DATE LAST MODIFIED: 2019/09/16

.########################################################################.
##############      PACKAGE INSTALLATION AND IMPORT      ###############
.########################################################################.

#Installing and Loading Necessary Packages
library(lme4)
library(gplots)
library(tidyverse)
library(lubridate)


#Importing Data - NOTE: The KSR climate data scripts are archived below and are no longer used as the data is unreliable.
climate <- read.csv("Rdata/climate_combined_data.csv")
lapply(climate, class)

climate <- climate %>%
  gather(key="variable", value="value", Temp:Total_Precip, factor_key=TRUE)
climate$Month <- factor(climate$Month, levels=c("April", "May", "June", "July", "August"))



lapply(climate, class)

.########################################################################.
######################      DATA EXPLORATION      #######################
.########################################################################.
Months <- c("April", "May", "June", "July", "August")

#Graphs: North York Temperature 2019 vs St. Georges 30-YR Normal
climate %>%
  filter(!variable=="Total_Precip" & !variable=="Temp_StDev") %>%
  ggplot(aes(x=D_Month, y=value, linetype=Station)) +
  geom_line(aes(color=variable), size=1.5) +
  scale_linetype_manual("Station", values=c("solid", "dotted"), labels=c("North York (ON)", "Saint-Georges (QC)")) +
  scale_color_manual("Temperature", values=c("gray0", "dodgerblue3", "firebrick3"), labels=c("Daily Average", "Minimum", "Maximum")) +
  guides(color=guide_legend("Temperature"), linetype=FALSE) + #Remove colour guide
  labs(y="Temperature (°C)") +
  scale_x_continuous("Month", labels=c("April", "May", "June", "July", "August")) +
  theme_bw() + theme(legend.title=element_blank(), legend.position="top", legend.box="vertical") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

#Graphs: North York Precipitation 2019 vs St. Georges 30-YR Normal
climate %>%
  filter(variable=="Total_Precip") %>%
  ggplot(aes(x=D_Month, y=value, fill=Station)) +
  geom_bar(color="black", position="dodge", stat="identity", size=1.2) +
  labs(y="Monthly Precipitation (mm)") +
  scale_x_continuous("Month", breaks=seq(4,8,1), labels=c("April", "May", "June", "July", "August")) + 
  scale_fill_manual("Station", values=c("steelblue3", "slategray1"), 
                    labels=c("North York (ON)", "Saint-Georges (QC)")) +
  guides(color=FALSE) + #Remove colour guide 
  labs(x="Month", y="Monthly Precipitation (mm)") +
  theme_bw() + theme(legend.title=element_blank(), legend.position="top") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

#


























#### Archive ####


ksr.climate <- read_csv("Rdata/heatarrays_climate_2019_data.csv", col_names = TRUE, na = "NAN")
ksr.climate$Station <- as.factor(ksr.climate$Station)
spec(ksr.climate) #or lapply(climate, class) works too

#This dataset includes temp + max/min and precipitation only
ksr.clim2 <- read_csv("Rdata/heatarrays_climate_fig_data.csv", col_names = TRUE, na = "NA")
ksr.clim2 <- as.data.frame(ksr.clim2)
ksr.clim2 <- na.omit(ksr.clim2)
ksr.clim2$Station <- as.factor(ksr.clim2$Station)
ksr.clim2$Timestamp <- as.Date(ksr.clim2$Timestamp)
lapply(ksr.clim2, class)
#spec(ksr.clim2) #or lapply(climate, class) works too if not a tibble

king.climate <- read_csv("Rdata/climate_normal_king_data.csv", col_names = TRUE)
king.climate$Station <- as.factor(king.climate$Station)
spec(king.climate)

stgeorge.climate <- read_csv("Rdata/climate_normal_stgeorges_data.csv", col_names = TRUE)
stgeorge.climate$Station <- as.factor(stgeorge.climate$Station)
spec(stgeorge.climate)


#Note: Plot 12 = Heated, Plot 10 = Ambient

#Calculating average climate data for KSR
ksr2 <- ksr.clim2 %>%
  group_by(Month=floor_date(Timestamp, "month"), Station) %>%
  summarise(Temp_KSR=mean(Temp, na.rm=TRUE), Temp_Min=mean(Temp_Min),
            Temp_Max=mean(Temp_Max), Precip_KSR=sum(Precipitation, na.rm=TRUE))
ksr2$Month <- months(as.Date(ksr2$Month))

#Graphs: Station Temperature vs Array Temps
ksr2$Month <- factor(ksr2$Month, levels=c("May", "June", "July", "August"))
ylim.primKSR <- c(0,100)
ylim.secKSR <- c(0,35)
bK <- diff(ylim.primKSR)/diff(ylim.secKSR)
aK <- bK*(ylim.primKSR[1] - ylim.secKSR[1])

ggplot(ksr2) +
  geom_bar(aes(x=Month, y=Precip_KSR), alpha=0.7, stat="identity", fill="steelblue1", color="steelblue4") +
  geom_line(aes(x=Month, y=aK + Temp_KSR*bK, group=1, colour="black"), size=2, linetype="solid") +
  geom_line(aes(x=Month, y=aK + Temp_Min*bK, group=1, color="dodgerblue4"), size=2, linetype="dashed") +
  geom_line(aes(x=Month, y=aK + Temp_Max*bK, group=1,  color="red4"), size=2, linetype="dashed") +
  ggtitle("A)", ) +
  scale_color_manual(values = c("black", "dodgerblue4", "red4"),
                     labels =c("Daily Avg", "Daily Avg Min", "Daily Avg Max"),
                     name = "") +
  scale_y_continuous("Precipitation (mm)", sec.axis = sec_axis( ~ (. -aK)/bK, name = "Temperature (°C)")) +
  labs(x="Date", y="Daily Average Temperature (°C)") +
  theme_bw() + theme(legend.position="top") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))

stgeorge.climate$Month <- factor(stgeorge.climate$Month, levels=c("May", "June", "July", "August"))
ylim.primStG <- c(0,140)
ylim.secStG <- c(0,35)

bStG <- diff(ylim.primStG)/diff(ylim.secStG)
aStG <- bStG*(ylim.primStG[1] - ylim.secStG[1])
ggplot(stgeorge.climate) +
  geom_bar(aes(x=Month, y=Precip_StG), alpha=0.7, stat="identity", fill="steelblue1", color="steelblue4") +
  geom_line(aes(x=Month, y=aStG + Temp_StG*bStG, color="black", group=1), size=2, linetype="solid") +
  geom_line(aes(x=Month, y=aStG + TempMin_StG*bStG, color="dodgerblue4",group=1), size=2, linetype="dashed") +
  geom_line(aes(x=Month, y=aStG + TempMax_StG*bStG, color="red4", group=1), size=2, linetype="dashed") +
  scale_color_manual(values = c("black", "dodgerblue4", "red4"),
                     labels =c("Daily Avg", "Daily Avg Min", "Daily Avg Max"),
                     name = "") +
  scale_y_continuous("Precipitation (mm)", sec.axis = sec_axis( ~ (. -aStG)/bStG, name = "Temperature (°C)")) +
  labs(x="Date", y="Daily Average Temperature (°C)") +
  theme_bw() + theme(legend.position="top") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))


  
#Organizing data into long format for available use
ksr <- ksr.climate %>%
  group_by(Date, Station) %>%
  summarise(Temp_KSR=mean(Temperature, na.rm=TRUE), RH_KSR=mean(RH), 
            SoilMoist_H=mean(Soil_Water_Content_12),
            Precip_KSR=sum(Precipitation, na.rm=TRUE), Heat_Index=mean(Heat_Index), Wind_Chill=mean(Wind_Chill),
            Temp_A=mean(Average_Reference_Temp, na.rm=TRUE), Temp_H=mean(Average_Heated_Temp, na.rm=TRUE)) %>%
  gather(key="variable", value="value", Temp_KSR:Temp_H, factor_key=TRUE)

#Transforming datasets into long format
#NOTE: NOT NECESSARY TO PRODUCE THE CLIMATE GRAPHS FOR KSR & ST. GEORGE
ksr2 <- ksr2 %>%
  gather(key="variable", value="value", Temp_KSR:Precip_KSR, factor_key=TRUE)
ksr2 <- as.data.frame(ksr2)

king <- king.climate %>%
  gather(key="variable", value="value", Temp_King:Precip_King, factor_key=TRUE)

stgeorge <- stgeorge.climate %>%
  gather(key="variable", value="value", Temp_StG:Precip_StG, factor_key=TRUE)
stgeorge <- as.data.frame(stgeorge)

climate <- rbind(ksr2, stgeorge) #Note, this does not work with tibbles. Otherwise, to use tibbles, use dplyr::bind_rows()
climate



  #Graphs: Station Temperature vs Array Temps
  climate %>%
    filter(variable=="Temp_KSR" | variable=="Temp_King") %>%
    ggplot(aes(x=Date, y=value, group=variable)) +
    geom_line(aes(linetype=variable, color=variable), size=1.2) +
    scale_linetype_manual("Station", values=c("twodash", "solid"), labels=c("King City 30-Yr Normal", "KSR")) +
    scale_color_manual("Station", values=c("grey5", "firebrick1"), labels=c("King City 30-Yr Normal", "KSR")) +
    labs(x="Date", y="Daily Average Temperature (°C)") +
    theme_bw() + theme(legend.position="top") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.5), 
          axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))
  
  #Graphs: Precipitation
  climate %>%
    filter(variable=="Precip_KSR" | variable=="Precip_King") %>%
    ggplot(aes(x=Date, y=value, group=variable)) +
    geom_line(aes(linetype=variable, color=variable), size=1.2) +
    scale_linetype_manual("Station", values=c("twodash", "solid"), labels=c("King City 30-Yr Normal", "KSR")) +
    scale_color_manual("Station", values=c("grey5", "steelblue2"), labels=c("King City 30-Yr Normal", "KSR")) +
    labs(x="Date", y="Daily Precipitation (mm)") +
    theme_bw() + theme(legend.position="top") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line = element_line(colour = "black", size = 1.5), 
          axis.ticks=element_line(size=1.5), axis.ticks.length=unit(0.2, "cm"))
  
  
  
