##
# D. montana project clinal populations
# PCA/QBGLM data analysis
# Last Modified: Feb 2017
##

# Clean environment
rm(list = ls(all = TRUE))

#
setwd("~/Desktop/Data/RData")
# Load libraries
require(raster)
library(FactoMineR)
library(ggplot2)
library(ggmap)
#citation("ggmap")
source("RScripts/ggplot_theme.R")

#----------------------------------------------#
# Get climate data for BayeScEnv and Quasi GLM
#----------------------------------------------#
# first load WC bio variables at the resolution of 2.5 deg
biod <- getData("worldclim", var="bio", res=2.5)
#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)) ###
#BIO3 = Isothermality (BIO2/BIO7) (* 100) ###
#BIO4 = Temperature Seasonality (standard deviation *100) ###
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter ###
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter)

tmind <- getData("worldclim", var="tmin", res=2.5)

tmaxd <- getData("worldclim", var="tmax", res=2.5)

precd <- getData("worldclim", var="prec", res=2.5)


# read csv file with geographic coordinates
geod<-read.table("~/Desktop/Data/montana_project/clinal_population_genomics/montana_pop_info.csv", header=T, stringsAsFactors=F,sep=",")
names(geod)
geod$Notes<-gsub("\n"," ",geod$Notes)
geod$Notes<-gsub("\r"," ",geod$Notes)
head(geod)
# Draw map
NA_W_map<-get_map(location=c(lon=mean(geod$Long[geod$Region=="NoAm"]),
                             lat=mean(geod$Lat[geod$Region=="NoAm"])),
                  zoom=3,maptype = "toner-lite",source = "stamen")
ggmap(NA_W_map)+
  geom_point(data=geod[geod$Region=="NoAm" & geod$PoolGen == 1,],
             aes(x=Long,y=Lat),fill="Blue",size=5,shape=21,alpha=1/2)+
  geom_point(data=geod[geod$Region=="NoAm" & geod$PoolGen == 0,],
             aes(x=Long,y=Lat),fill="Red",size=5,shape=24,alpha=1/2)+
  geom_label(data=geod[geod$Region=="NoAm"& geod$PoolGen == 1,],
             aes(x = Long+c(-10,-10,-10,-14),
                 y=Lat+c(0,0,0,-2),label=LocationName),size = 7)+
  xlab("Longitude")+
  ylab("Latitude")+my.theme+theme(axis.text = element_text(size=14),
                                  text = element_text(size = 14))

Eur_N_map<-get_map(location=c(lon=mean(geod$Long[geod$Region=="Eur"]),
                             lat=mean(geod$Lat[geod$Region=="Eur"])),
                  zoom=5,maptype = "toner-lite",source = "stamen")
ggmap(Eur_N_map)+
  geom_point(data=geod[geod$Region=="Eur" & geod$PoolGen == 1,],
             aes(x=Long,y=Lat),fill="Blue",size = 5,shape=21,alpha=1/2)+
  geom_point(data=geod[geod$Region=="Eur"& geod$PoolGen == 0,],
             aes(x=Long,y=Lat),fill="Red",size = 5,shape=22,alpha=1/2)+
  geom_label(data=geod[geod$Region=="Eur" & geod$PoolGen == 1,],
             aes(x = Long+c(3,3),y=Lat,label=LocationName),size = 7)+
  xlab("Longitude")+
  ylab("Latitude")+my.theme+theme(axis.text = element_text(size=14),
                                  text = element_text(size = 14))

# extact for each coordinate bio clim variables
# extract(biovariable, c(long,lat))
bio<-extract(biod, geod[,c(5,4)])
tmin<-extract(tmind, geod[,c(5,4)])
tmax<-extract(tmaxd, geod[,c(5,4)])
precd<-extract(precd, geod[,c(5,4)])

# create a full dataset
geod<-cbind(geod,bio,tmin,tmax,precd)

# save into external file
#write.table(geod,
#            file="~/Desktop/Data/montana_project/clinal_population_genomics#/TableS1_montana_pop_info_bioclimdat.csv",
#            sep=",", row.names=FALSE ,quote=FALSE)

#######
# Run PCA for all populations
#######
#to get variance and mean for PCA
geo<-geod[,11:65]
names(geo)
rownames(geo)<-geod[,3]
# Run PCA
pca<-PCA(geo)
# Plots
pca_dat<-as.data.frame(pca$ind$coord)
pca_dat$loc<-rownames(pca_dat)
#Modify position of "Seward" label
pca_dat_labs<-pca_dat[,c(1,2,6)]
pca_dat_labs[1,2]<-pca_dat_labs[1,2]-0.3
pca_dat_labs$reg <- geod$Region
pca_dat_labs$PoolGen <- geod$PoolGen
head(pca_dat_labs)

ggplot()+
  geom_point(data=pca_dat_labs[pca_dat_labs$PoolGen == 0 &
                                 pca_dat_labs$reg=="NoAm",],
             aes(Dim.1,Dim.2),shape = 24, fill = "red",size = 3)+
  geom_point(data=pca_dat_labs[pca_dat_labs$PoolGen == 0 &
                                 pca_dat_labs$reg=="Eur",],
             aes(Dim.1,Dim.2),shape = 22, fill = "red",size = 3)+
  
  geom_point(data=pca_dat_labs[pca_dat_labs$PoolGen == 1,],
             aes(Dim.1,Dim.2),shape = 21, fill = "blue",size = 3)+

  geom_text(data=pca_dat_labs[pca_dat_labs$PoolGen==1,],
             aes(Dim.1-0.5,Dim.2+0.45,label=loc),size = 4,colour="blue")+
  geom_text(data=pca_dat_labs[pca_dat_labs$PoolGen==0 &
                                pca_dat_labs$reg=="Eur",],
            aes(Dim.1-c(0,-0.5,-0.8,2.5),
                Dim.2+c(0.45,-0.45,0.45,0),
                label=loc),size = 4,colour="red")+
  geom_text(data=pca_dat_labs[pca_dat_labs$PoolGen==0 &
                                pca_dat_labs$reg=="NoAm",],
            aes(Dim.1-c(0,0,0,0,0,0,0,
                        0,0,-2.5,0,1,0,0),
                Dim.2+c(0.45,0.45,0.45,-0.45,0.45,0.45,0.45,
                        0.45,0.45,0.45,0.45,0.45,0.45,-0.45),
                label=loc),size = 4,colour="red")+
  
  scale_x_continuous(limits = c(-15,15))+
  xlab(paste("PC1 (",round(pca$eig$`percentage of variance`[1],2),"%)"))+
  ylab(paste("PC2 (",round(pca$eig$`percentage of variance`[2],2),"%)"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  my.theme+theme(axis.text = element_text(size = 14),
                 axis.title = element_text(size = 14))

pca_dat_labs2<-pca_dat[,c(3,4,6)]
pca_dat_labs2[1,2]<-pca_dat_labs[1,2]-0.3
pca_dat_labs2$reg <- geod$Region
pca_dat_labs2$PoolGen <- geod$PoolGen

ggplot()+
  geom_point(data=pca_dat_labs2[pca_dat_labs2$PoolGen == 0 &
                                 pca_dat_labs$reg=="NoAm",],
             aes(Dim.3,Dim.4),shape = 24, fill = "red",size = 3)+
  geom_point(data=pca_dat_labs2[pca_dat_labs2$PoolGen == 0 &
                                 pca_dat_labs$reg=="Eur",],
             aes(Dim.3,Dim.4),shape = 22, fill = "red",size = 3)+
  
  geom_point(data=pca_dat_labs2[pca_dat_labs2$PoolGen == 1,],
             aes(Dim.3,Dim.4),shape = 21, fill = "blue",size = 3)+
  geom_label(data=pca_dat_labs2[pca_dat_labs2$PoolGen == 1,],
             aes(Dim.3+c(1,1.5,1,1,-1,1),Dim.4+c(0.4,0.1,0.3,0.3,0.3,0.3),label=loc),
             size=4)+
  scale_x_continuous(limits = c(-10,10))+
  xlab(paste("PC3 (",round(pca$eig$`percentage of variance`[3],2),"%)"))+
  ylab(paste("PC4 (",round(pca$eig$`percentage of variance`[4],2),"%)"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  my.theme+theme(axis.text = element_text(size = 14),
                 axis.title = element_text(size = 14))

# Make a "Scree" plot
pca_eig_dat<-as.data.frame(pca$eig)
pca_eig_dat$comp<-rownames(pca_eig_dat)
pca_eig_dat$comp <- gsub("comp ","",pca_eig_dat$comp)
pca_eig_dat$comp <- as.numeric(pca_eig_dat$comp)
pca_eig_dat$comp<-factor(pca_eig_dat$comp,
                         labels=paste("PC",seq(1,23),sep=""))

plot1<-ggplot(data=pca_eig_dat)+
  geom_bar(aes(comp,eigenvalue),stat = "identity",width = 0.5)+
  geom_hline(yintercept = 1,linetype="dashed",colour="blue")+
  xlab("")+
  ylab("Eigenvalues")+
  my.theme+theme(axis.text = element_text(size = 14),
                 axis.text.x = element_text(angle=45,
                                            hjust=1,vjust=1),
                 axis.title = element_text(size = 14))
plot(plot1)
# First 4 vectors have eigenvalues > 1
plot1_ins<-ggplot(data=pca_eig_dat)+
  geom_bar(aes(comp,`cumulative percentage of variance`),
           stat = "identity",width = 0.5)+
  geom_hline(yintercept = 98,linetype="dashed",colour="blue")+
  xlab("")+
  ylab("Cumulative Variance\nExplained (%)")+
  my.theme+theme(axis.text = element_text(size = 14),
                 axis.text.x = element_text(size=8,angle=45,
                                            hjust=1,vjust=1),
                 axis.title = element_text(size = 14))
plot(plot1_ins)

# make plot1_ins an "inset" in plot1
# Specify position of plot2 (in percentages of plot1)
# This is in the top left and 25% width and 25% height
xleft   = 0.25
xright  = 0.95
ybottom = 0.45
ytop    = 0.95 

# Calculate position in plot1 coordinates
# Extract x and y values from plot1
l1 = ggplot_build(plot1)
x1 = l1$layout$panel_ranges[[1]]$x.range[1]
x2 = l1$layout$panel_ranges[[1]]$x.range[2]
y1 = l1$layout$panel_ranges[[1]]$y.range[1]
y2 = l1$layout$panel_ranges[[1]]$y.range[2]
xdif = x2-x1
ydif = y2-y1
xmin  = x1 + (xleft*xdif)
xmax  = x1 + (xright*xdif)
ymin  = y1 + (ybottom*ydif)
ymax  = y1 + (ytop*ydif) 

# Get plot2 and make grob
g2 = ggplotGrob(plot1_ins)
plot3 = plot1 + annotation_custom(grob = g2, 
                                  xmin=xmin, xmax=xmax, 
                                  ymin=ymin, ymax=ymax)
plot(plot3)


# First 4 Principal Components have eigenvalues > 1
# Run PCA again keeping only the first 4 dimensions.
pca<-PCA(geo,ncp = 4)
pca$var$cor
cat(capture.output(pca$var$coord),
    file="~/Desktop/Data/montana_project/clinal_population_genomics/montana_clinalpops-GIS_PCA.rot",
    sep="")
pca$ind
pca$var$cor

# Describe the first 4 PCA components
dimdesc(pca,axes = c(1,2,3,4))
# pca$ind$coord
# Ztransform the PCA results
pca$ind$coord[,1]<-scale(x = pca$ind$coord[,1],center = TRUE,scale=TRUE)
pca$ind$coord[,2]<-scale(x = pca$ind$coord[,2],center = TRUE,scale=TRUE)
pca$ind$coord[,3]<-scale(x = pca$ind$coord[,3],center = TRUE,scale=TRUE)
pca$ind$coord[,4]<-scale(x = pca$ind$coord[,4],center = TRUE,scale=TRUE)
pca$ind$coord
plot(pca$ind$coord[,1],pca$ind$coord[,2])
# Bind the PCA results to the main data
geod[,c(1,3)]
colnames(pca$ind$coord)<-gsub("Dim.","PC",colnames(pca$ind$coord))
bio.data<-cbind(geod,pca$ind$coord)

write.table(bio.data,
            file="~/Desktop/Data/montana_project/clinal_population_genomics/TableS1_montana_pop_info_bioclimdat.csv",
            sep="\t", row.names=FALSE ,quote=FALSE)

# Read in table S1
bio.data<-read.table("~/Desktop/Data/montana_project/clinal_population_genomics/TableS1_montana_pop_info_bioclimdat.csv",
                     header=TRUE,sep=",")
# Write PC1 and PC2 in the order in which they are run in mpileup
names(bio.data)
order<-c("Oulanka","Korpilahti","Ashford",
        "Seward","Crested Butte","Terrace")
match(order,bio.data$LocationName)
Pool_seq_alt_PCdata<-bio.data[match(order,bio.data$LocationName),
                          c(3,6,66,67)]
Pool_seq_alt_lat<-bio.data[match(order,bio.data$LocationName),c(3,4,6)]

# Write to file
write.table(Pool_seq_alt_PCdata,
            file="~/Desktop/Data/montana_project/clinal_population_genomics/results/climPCA/montana_poolpops_PC1PC2.csv",
            sep=",", row.names=FALSE ,quote=FALSE)
write.table(Pool_seq_alt_lat,
            file="~/Desktop/Data/montana_project/clinal_population_genomics/results/climPCA/montana_poolpops_lat.csv",
            sep=",", row.names=FALSE ,quote=FALSE)



####
# Correlation between altitude and PC1 and PC2
####
bio.data[bio.data$LocationName=="Crested Butte",c(6,66,67)]

cor.test(bio.data$Alt_m,bio.data$PC1)
ggplot()+
  geom_point(data=bio.data,aes(Alt_m,PC1))+
  xlab("Altitude (m)")+
  ylab("PC1")+
  geom_label(data=bio.data[bio.data$LocationName=="Crested Butte",],
             aes(Alt_m-200,PC1+0.1,label=LocationName))+
  geom_label(aes(x=2500,y=1,label="cor = -0.398\np = 0.092"))+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text = element_text(size=10))

cor.test(bio.data$Alt_m,bio.data$PC2)
ggplot()+
  geom_point(data=bio.data,aes(Alt_m,PC2))+
  xlab("Altitude (m)")+
  ylab("PC2")+
  geom_label(data=bio.data[bio.data$LocationName=="Crested Butte",],
             aes(Alt_m-200,PC2+0.1,label=LocationName))+
  geom_label(aes(x=2500,y=1,label="cor = -0.578\np = 0.0096"))+
  my.theme+theme(axis.title = element_text(size=10),
                 axis.text = element_text(size=10))

#######
# Run PCA analysis on only the 6 populations for pool seq
#######
geod_pool <- geod[geod$PoolGen == 1,] 

#to get variance and mean for PCA
geo<-geod_pool[,11:65]
rownames(geo)<-geod_pool[,3]
# Run PCA
pca<-PCA(geo)
# Plots
pca_dat<-as.data.frame(pca$ind$coord)
pca_dat$loc<-rownames(pca_dat)
#Modify position of "Seward" label
pca_dat_labs<-pca_dat[,c(1,2,6)]
pca_dat_labs[1,2]<-pca_dat_labs[1,2]-0.3
pca_dat_labs$reg <- geod_pool$Region
pca_dat_labs$PoolGen <- geod_pool$PoolGen

ggplot()+
  geom_point(data=pca_dat_labs,
             aes(Dim.1,Dim.2),size = 3)+
  geom_label(data=pca_dat_labs,
             aes(Dim.1-0.5,Dim.2+0.45,label=loc),size = 4)+
  scale_x_continuous(limits = c(-10,10))+
  xlab(paste("PC1 (",round(pca$eig$`percentage of variance`[1],2),"%)"))+
  ylab(paste("PC2 (",round(pca$eig$`percentage of variance`[2],2),"%)"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  my.theme+theme(axis.text = element_text(size = 14),
                 axis.title = element_text(size = 14))

pca_dat_labs2<-pca_dat[,c(3,4,6)]
pca_dat_labs2[1,2]<-pca_dat_labs[1,2]-0.3
pca_dat_labs2$reg <- geod$Region
pca_dat_labs2$PoolGen <- geod$PoolGen

ggplot()+
  geom_point(data=pca_dat_labs2,
             aes(Dim.3,Dim.4),size = 3)+
  geom_label(data=pca_dat_labs2,
             aes(Dim.3+0.1,Dim.4+0.15,label=loc),
             size=4)+
  scale_x_continuous(limits = c(-4,4))+
  xlab(paste("PC3 (",round(pca$eig$`percentage of variance`[3],2),"%)"))+
  ylab(paste("PC4 (",round(pca$eig$`percentage of variance`[4],2),"%)"))+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  my.theme+theme(axis.text = element_text(size = 14),
                 axis.title = element_text(size = 14))

# Make a "Scree" plot
pca_eig_dat<-as.data.frame(pca$eig)
pca_eig_dat$comp<-rownames(pca_eig_dat)
pca_eig_dat$comp <- gsub("comp ","",pca_eig_dat$comp)
pca_eig_dat$comp <- as.numeric(pca_eig_dat$comp)

plot1<-ggplot(data=pca_eig_dat)+
  geom_bar(aes(comp,eigenvalue),stat = "identity",width = 0.5)+
  geom_hline(yintercept = 1,linetype="dashed",colour="blue")+
  xlab("")+
  ylab("Eigenvalues")+
  my.theme+theme(axis.text = element_text(size = 10),
                 axis.title = element_text(size = 10))
plot(plot1)
pca_eig_dat$comp<-factor(pca_eig_dat$comp,
                         labels=c("PC1","PC2","PC3","PC4","PC5"))
# First 4 vectors have eigenvalues > 1
plot1_ins<-ggplot(data=pca_eig_dat)+
  geom_bar(aes(comp,`cumulative percentage of variance`),
           stat = "identity",width = 0.5)+
#  geom_hline(yintercept = 80,linetype="dashed",colour="blue")+
  xlab("")+
  ylab("Cumulative Variance\nExplained (%)")+
  my.theme+theme(axis.text = element_text(size = 10),
                 axis.title = element_text(size = 10))
plot(plot1_ins)

# make plot1_ins an "inset" in plot1
# Specify position of plot2 (in percentages of plot1)
# This is in the top left and 25% width and 25% height
xleft   = 0.4
xright  = 0.95
ybottom = 0.55
ytop    = 0.95 

# Calculate position in plot1 coordinates
# Extract x and y values from plot1
l1 = ggplot_build(plot1)
x1 = l1$layout$panel_ranges[[1]]$x.range[1]
x2 = l1$layout$panel_ranges[[1]]$x.range[2]
y1 = l1$layout$panel_ranges[[1]]$y.range[1]
y2 = l1$layout$panel_ranges[[1]]$y.range[2]
xdif = x2-x1
ydif = y2-y1
xmin  = x1 + (xleft*xdif)
xmax  = x1 + (xright*xdif)
ymin  = y1 + (ybottom*ydif)
ymax  = y1 + (ytop*ydif) 

# Get plot2 and make grob
g2 = ggplotGrob(plot1_ins)
plot3 = plot1 + annotation_custom(grob = g2, 
                                  xmin=xmin, xmax=xmax, 
                                  ymin=ymin, ymax=ymax)
plot(plot3)


# First 4 Principal Components have eigenvalues > 1
# Run PCA again keeping only the first 4 dimensions.
pca<-PCA(geo,ncp = 4)
pca$var$cor
cat(capture.output(pca$var$coord),
    file="~/Desktop/Data/montana_project/clinal_population_genomics/montana_clinalpops-GIS_PCA_popsub.rot",
    sep="")
pca$ind
pca$var$cor
# Describe the first 4 PCA components
dimdesc(pca,axes = c(1,2,3,4))


