#### PROJECT: Brassica rapa Indirect Genetic Effects (Data collected by Cameron So 2019 at Koffler Scientific Reserve, King City, ON)
#### PURPOSE: Assess Indirect Genetic Effects (IGE) using heatmaps and update data with competitor values
#### AUTHOR: Cameron So
#### DATE LAST MODIFIED: 2020/01/28

#'####################################################################'#
##############      PACKAGE INSTALLATION AND IMPORT      ###############
#'####################################################################'#

library(tidyverse)

col_types_list <- cols_only(posID = "d", individual = "d", plot = col_factor(levels=c(1:12)),
                            animal = "d", 
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
                            height = "d", stem_diam = "d",
                            posX = "d", posY = "d", germ_raster = "d", cone_vol = "d"
)

data <- read_csv("Rdata/heatarrays_brassica_2019_data.csv", col_names = TRUE, na = "NA", 
                 col_types=col_types_list)

lapply(data, class)


#'####################################################################'#
 #######################      SPATIAL MAPS     #######################
#'####################################################################'#

#Spatial Maps using ggplot2 - NOTE: the package 'raster' also creates heatmaps and be visualized below.

#Select Columns of Data Desired
heatlist.germ <- subset(data, select = c(posX, posY, germ_bin, flwr_bin, seed_bin, plot))

#Split dataframe into multiple based on plot grouping, and create a list
heatlist.germ <- split(heatlist.germ, f=as.factor(heatlist.germ$plot))

#List of all of heatmaps: - NOTE: There is an EXCELLENT explanation of how to create functions in the thread I created on Stack Overflow

heatlist.germ <- heatlist.germ 
plot_data_fcn <- function (heatlist.germ) {
  ggplot(heatlist.germ, aes(x=posX, y=posY, fill=germ_bin)) + 
    geom_tile(aes(fill=germ_bin)) + 
    geom_text(aes(label=germ_bin)) +
    scale_fill_gradient(low = "gray90", high="darkolivegreen4") +
    ggtitle(aes(label=plot)) +
    scale_x_continuous("Position X", breaks=seq(1,30)) +
    scale_y_continuous("Position Y (REVERSED)", breaks=seq(1,20))
  }

heatlist.test <- lapply(heatlist.germ, plot_data_fcn)

#To access a heat map, just insert the plot number into the bracketed area
plot(heatlist.test[[1]])

#'####################################################################'#
 ##################      RASTER - IGE EXTRACTION     #################
#'####################################################################'#

#Note, this section is the script to extract the neighbourhood data (germ_neigh, cone_vol) that was attached to the dataset for analyses in scripts 01 and 02. 

## Extraction of Sum of Germinated Neighbours ##

library(raster)

#Select Columns of Data Desired
out <- subset(data, select = c(posX, posY, germ_bin, flwr_bin, seed_bin, plot))

#Split dataframe into multiple based on plot grouping, and create a list
out <- split(out, f=as.factor(out$plot))

#Convert each df into a Raster
out <- lapply(out, rasterFromXYZ)

#Renaming dataframes in the list for clarity
names(out) <- c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12")

#To access the dataframe individually, use the following
out[["p1"]] #or out[[1]]

#Applying multi-Focal function on each rasterBrick in the list
#Note: Layer 1=germ, Layer 2=flower, Layer 3=Seed Mat, Layer 4=Plot (ignore)
out <- lapply(names(out), function(x) multiFocal(out[[x]], w=matrix(1,3,3), fun=sum, na.rm=TRUE, pad=TRUE))

#To plot a raster layer of sum of germinated neighbours, try
plot(out[[1]])
text(out[[1]])

#Extracting Values
plot1_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot1_germ, file="Rdata/Spatial_data/Germ_spatial_p1.csv")

plot2_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot2_germ, file="Rdata/Spatial_data/Germ_spatial_p2.csv")

plot3_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot3_germ, file="Rdata/Spatial_data/Germ_spatial_p3.csv")

plot4_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot4_germ, file="Rdata/Spatial_data/Germ_spatial_p4.csv")

plot5_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot5_germ, file="Rdata/Spatial_data/Germ_spatial_p5.csv")

plot6_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot6_germ, file="Rdata/Spatial_data/Germ_spatial_p6.csv")

plot7_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot7_germ, file="Rdata/Spatial_data/Germ_spatial_p7.csv")

plot8_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot8_germ, file="Rdata/Spatial_data/Germ_spatial_p8.csv")

plot9_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot9_germ, file="Rdata/Spatial_data/Germ_spatial_p9.csv")

plot10_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot10_germ, file="Rdata/Spatial_data/Germ_spatial_p10.csv")

plot11_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot11_germ, file="Rdata/Spatial_data/Germ_spatial_p11.csv")

plot12_germ <- as.data.frame(out[[1]]$layer.1, xy=TRUE)
write.csv(plot12_germ, file="Rdata/Spatial_data/Germ_spatial_p12.csv")








## Extraction of Sum of Cone Volumes per Germinated Neighbours ##

library(raster)

#Select Columns of Data Desired
out2 <- subset(data, select = c(posX, posY, cone_vol, plot))

#Split dataframe into multiple based on plot grouping, and create a list
out2 <- split(out2, f=as.factor(out2$plot))

#Convert each df into a Raster
out2 <- lapply(out2, rasterFromXYZ)

#Renaming dataframes in the list for clarity
names(out2) <- c("p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10", "p11", "p12")

#To access the dataframe individually, use the following
out2[["p1"]] #or out[[1]]

#Applying multi-Focal function on each rasterBrick in the list
#Note: Layer 1=germ, Layer 2=flower, Layer 3=Seed Mat, Layer 4=Plot (ignore)
out2 <- lapply(names(out2), function(x) multiFocal(out2[[x]], w=matrix(1,3,3), fun=sum, na.rm=TRUE, pad=TRUE))

#To plot a raster layer of sum of germinated neighbours, try
plot(out2[[1]])
text(out2[[1]])

#Extracting Values
plot1_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot1_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p1.csv")

plot2_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot2_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p2.csv")

plot3_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot3_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p3.csv")

plot4_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot4_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p4.csv")

plot5_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot5_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p5.csv")

plot6_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot6_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p6.csv")

plot7_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot7_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p7.csv")

plot8_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot8_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p8.csv")

plot9_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot9_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p9.csv")

plot10_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot10_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p10.csv")

plot11_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot11_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p11.csv")

plot12_coneVol <- as.data.frame(out2[[1]]$layer.1, xy=TRUE)
write.csv(plot12_coneVol, file="Rdata/Spatial_data/ConeVol_spatial_p12.csv")

#'####################################################################'#
  ############      CORRELATION WITH EXTRACTED DATA     ##############
#'####################################################################'#

load(file="Routput/heatarrays_data_explore.RData")

#Converting all values below zero to zero for neighbourhood volume
dat$cone_vol_neigh[dat$cone_vol_neigh < 0] <- 0

#Number of Neighbours against Seed Pod Number - INDIVIDUAL LEVEL
dat %>%
  filter(!germ==0) %>%
  ggplot(aes(x=as.factor(germ_neigh), y=seed_pods, color=as.factor(treatment))) +
  geom_boxplot(aes(color=treatment), alpha=0.5) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Number of Germinated Neighbours", y="Number of Seed Pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Seed Pod Number against Seed Pod Number - FAMILY level
dat %>%
  filter(!germ==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), germ_neigh_avg=mean(germ_neigh)) %>%
  ggplot(aes(x=germ_neigh_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Average Number of Germinated Neighbours", y="Average Seed Pods") +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Number of Neighbours against Height - INDIVIDUAL LEVEL
dat %>%
  filter(!germ==0) %>%
  ggplot(aes(x=as.factor(germ_neigh), y=height, color=treatment)) +
  geom_boxplot(aes(color=treatment), alpha=0.5) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="NNumber of Germinated Neighbours", y="Height") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Number of Neighbours against Height - FAMILY level
dat %>%
  filter(!germ==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), germ_neigh_avg=mean(germ_neigh)) %>%
  ggplot(aes(x=germ_neigh_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Average Number of Germinated Neighbours", y="Average Height") +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))














#Neighbourhood Volume (cm^2) against Seed Pod Number - INDIVIDUAL LEVEL
dat %>%
  filter(!germ==0) %>%
  ggplot(aes(x=cone_vol_neigh, y=seed_pods, color=as.factor(treatment))) +
  geom_point(aes(color=treatment), alpha=0.5) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Neighbourhood Volume (cm^2)", y="Number of Seed Pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Neighbourhood Volume (cm^2) against Seed Pod Number - FAMILY level
dat %>%
  filter(!germ==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), cone_vol_neigh_avg=mean(cone_vol_neigh)) %>%
  ggplot(aes(x=cone_vol_neigh_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Cumulative Neighbourhood Volume by Family (cm^2)", y="Average Seed Pods") +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))


#Neighbourhood Volume (cm^2) against Height - INDIVIDUAL LEVEL
dat %>%
  filter(!germ==0) %>%
  ggplot(aes(x=cone_vol_neigh, y=height, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x="Neighbourhood Volume (cm^2)", y="Number of Seed Pods") +
  ggtitle("A)") + 
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))

#Neighbourhood Volume (cm^2) against Height - FAMILY level
dat %>%
  filter(!germ==0) %>%
  group_by(treatment, famID) %>%
  summarise(pods_avg=mean(seed_pods), cone_vol_neigh_avg=mean(cone_vol_neigh)) %>%
  ggplot(aes(x=cone_vol_neigh_avg, y=pods_avg, color=treatment)) +
  geom_point(aes(color=treatment), alpha=0.5) +
  geom_smooth(se=F, method="glm") +
  scale_color_manual("Treatment", values=c("dodgerblue2", "tomato2"), labels=c("Ambient", "Heated")) +
  labs(x= "Cumulative Neighbourhood Volume (cm^2)", y="Average Height (cm)") +
  ggtitle("B)") +
  theme_bw() + theme(legend.position="none") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black", size = 1.5), 
        axis.ticks=element_line(size=1), axis.ticks.length=unit(0.2, "cm"))






















### Function: multi-Focal on a rasterBrick ###
#We use a 3x3 grid for focal analysis
multiFocal <- function(x, w=matrix(1, nr=3, nc=3), ...) {
  
  
  if(is.character(x)) {
    x <- brick(x)
  }
  # The function to be applied to each individual layer
  fun <- function(ind, x, w, ...){
    focal(x[[ind]], w=w, ...)
  }
  
  
  n <- seq(nlayers(x))
  list <- lapply(X=n, FUN=fun, x=x, w=w, ...)
  
  
  out <- stack(list)
  return(out)
}
