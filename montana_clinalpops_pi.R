##
# D. montana project clinal populations
# QBGLM data analysis
# Last Modified: Feb 2017
##
# Clean environment
rm(list = ls(all = TRUE))

#
setwd("~/Desktop/Data/RData/montana_clines/")
# Load libraries
library(ggplot2)
library(qvalue)
library(reshape)
library(reshape2)
se <- function(x) {sqrt(var(x, na.rm = TRUE))/sqrt(length(x))}
#citation("ggmap")
source("~/Desktop/Data/RData/RScripts/ggplot_theme.R")
options(scipen=10)

#
# Read in the pi data.
# ####
Oulanka_pi<-read.table("Oulanka.pi",header=FALSE)
# Latitude: 66.67
colnames(Oulanka_pi)<-c("chr","winpos","NrSNPs","PropGdReads","pi")
Oulanka_pi$pi<-as.numeric(as.character(Oulanka_pi$pi))
Oulanka_pi$pop<-rep("Oulanka",nrow(Oulanka_pi))
head(Oulanka_pi)
mean(Oulanka_pi$pi,na.rm = TRUE)
# Pi distribution
ggplot()+geom_histogram(data=Oulanka_pi,aes(pi))


Korpilahti_pi<-read.table("Korpilahti.pi",header=FALSE)
# Latitude: 62.33
colnames(Korpilahti_pi)<-c("chr","winpos","NrSNPs","PropGdReads","pi")
Korpilahti_pi$pi<-as.numeric(as.character(Korpilahti_pi$pi))
Korpilahti_pi$pop<-rep("Korpilahti",nrow(Korpilahti_pi))
head(Korpilahti_pi)
mean(Korpilahti_pi$pi,na.rm = TRUE)


Ashford_pi<-read.table("Ashford.pi",header=FALSE)
# Latitude: 46.75
colnames(Ashford_pi)<-c("chr","winpos","NrSNPs","PropGdReads","pi")
Ashford_pi$pi<-as.numeric(as.character(Ashford_pi$pi))
Ashford_pi$pop<-rep("Ashford",nrow(Ashford_pi))
head(Ashford_pi)
mean(Ashford_pi$pi,na.rm = TRUE)


Seward_pi<-read.table("Seward.pi",header=FALSE)
# Latitude: 60.15
colnames(Seward_pi)<-c("chr","winpos","NrSNPs","PropGdReads","pi")
Seward_pi$pi<-as.numeric(as.character(Seward_pi$pi))
Seward_pi$pop<-rep("Seward",nrow(Seward_pi))
head(Seward_pi)
mean(Seward_pi$pi,na.rm = TRUE)

Crested_Butte_pi<-read.table("Crested_Butte.pi",header=FALSE)
# Latitude: 38.9
colnames(Crested_Butte_pi)<-c("chr","winpos","NrSNPs","PropGdReads","pi")
Crested_Butte_pi$pi<-as.numeric(as.character(Crested_Butte_pi$pi))
Crested_Butte_pi$pop<-rep("Crested Butte",nrow(Crested_Butte_pi))
head(Crested_Butte_pi)
mean(Crested_Butte_pi$pi,na.rm = TRUE)
# Pi distribution
ggplot()+geom_histogram(data=Crested_Butte_pi,aes(pi))

Terrace_pi<-read.table("Terrace.pi",header=FALSE)
# Latitude: 54.45
colnames(Terrace_pi)<-c("chr","winpos","NrSNPs","PropGdReads","pi")
Terrace_pi$pi<-as.numeric(as.character(Terrace_pi$pi))
Terrace_pi$pop<-rep("Terrace",nrow(Terrace_pi))
head(Terrace_pi)
mean(Terrace_pi$pi,na.rm = TRUE)

# ####

# Plot mean pi as a function of latitude
pidat<-data.frame(pi=c(mean(Oulanka_pi$pi,na.rm = TRUE),
                       mean(Korpilahti_pi$pi,na.rm = TRUE),
                       mean(Ashford_pi$pi,na.rm = TRUE),
                       mean(Seward_pi$pi,na.rm = TRUE),
                       mean(Terrace_pi$pi,na.rm = TRUE),
                       mean(Crested_Butte_pi$pi,na.rm = TRUE)),
                  piSE=c(se(Oulanka_pi$pi),
                         se(Korpilahti_pi$pi),
                         se(Ashford_pi$pi),
                         se(Seward_pi$pi),
                         se(Terrace_pi$pi),
                         se(Crested_Butte_pi$pi)),
                  pop=c("Oulanka","Korpilahti","Ashford",
                        "Seward","Terrace","Crested Butte"),
                  lat=c(66.67,62.33,46.75,60.15,54.45,38.9))
pidat
ggplot()+geom_point(data=pidat,aes(lat,pi),size=2)+
  geom_errorbar(data=pidat,aes(x = lat,ymin=pi-piSE,ymax=pi+piSE),
                colour="red",width=1)+
  geom_label(data=pidat,aes(x = lat,y=pi+0.0003,label=pop))+
  xlab("Latitude")+
  ylab(expression(pi))+
  my.theme+theme(axis.text = element_text(size=10),
                 axis.title.y = element_text(size=14))

# Plot pi box plots
all_pi_dat<-rbind(Oulanka_pi,Korpilahti_pi,Ashford_pi,
                  Seward_pi,Terrace_pi,Crested_Butte_pi)
all_pi_dat$pop<-factor(all_pi_dat$pop,
                       levels = c("Oulanka","Korpilahti",
                                  "Seward","Terrace",
                                  "Ashford","Crested Butte"))
str(all_pi_dat)

ggplot()+
  geom_boxplot(data=all_pi_dat,aes(pop,pi))+
  xlab("")+
  ylab(expression(pi))+
  ggtitle("A)")+
  my.theme+theme(axis.text = element_text(size=10),
                 axis.title.y = element_text(size=14),
                 plot.title = element_text(size=14))
ggsave("~/Desktop/Data/montana_project/clinal_population_genomics/results/plots/montana_clines_pi.png",device = "png",dpi=600)


