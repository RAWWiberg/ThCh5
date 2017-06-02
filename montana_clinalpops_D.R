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
Oulanka_D<-read.table("Oulanka.D",header=FALSE)
# Latitude: 66.67
colnames(Oulanka_D)<-c("chr","pos","NrSNPs","PropGdReads","D")
Oulanka_D$D<-as.numeric(as.character(Oulanka_D$D))
Oulanka_D$pop<-rep("Oulanka",nrow(Oulanka_D))
# Remove the "size" from scaffold names
Oulanka_D$chr<-gsub("-size.*","",Oulanka_D$chr)
head(Oulanka_D)
summary(Oulanka_D$D,na.rm = TRUE)
quantile(Oulanka_D$D,probs = 0.05,na.rm = TRUE)
# D distribution
ggplot()+geom_histogram(data=Oulanka_D,aes(D))+ggtitle("Oulanka")

Korpilahti_D<-read.table("Korpilahti.D",header=FALSE)
# Latitude: 62.33
colnames(Korpilahti_D)<-c("chr","pos","NrSNPs","PropGdReads","D")
Korpilahti_D$D<-as.numeric(as.character(Korpilahti_D$D))
Korpilahti_D$pop<-rep("Korpilahti",nrow(Korpilahti_D))
# Remove the "size" from scaffold names
Korpilahti_D$chr<-gsub("-size.*","",Korpilahti_D$chr)
head(Korpilahti_D)
summary(Korpilahti_D$D,na.rm = TRUE)
quantile(Korpilahti_D$D,probs = 0.05,na.rm = TRUE)
# D distribution
ggplot()+geom_histogram(data=Korpilahti_D,aes(D))+ggtitle("Korpilahti")

Seward_D<-read.table("Seward.D",header=FALSE)
# Latitude: 60.15
colnames(Seward_D)<-c("chr","pos","NrSNPs","PropGdReads","D")
Seward_D$D<-as.numeric(as.character(Seward_D$D))
Seward_D$pop<-rep("Seward",nrow(Seward_D))
# Remove the "size" from scaffold names
Seward_D$chr<-gsub("-size.*","",Seward_D$chr)
head(Seward_D)
summary(Seward_D$D,na.rm = TRUE)
quantile(Seward_D$D,probs = 0.05,na.rm = TRUE)
# D distribution
ggplot()+geom_histogram(data=Seward_D,aes(D))+ggtitle("Seward")

Ashford_D<-read.table("Ashford.D",header=FALSE)
# Latitude: 46.75
colnames(Ashford_D)<-c("chr","pos","NrSNPs","PropGdReads","D")
Ashford_D$D<-as.numeric(as.character(Ashford_D$D))
Ashford_D$pop<-rep("Ashford",nrow(Ashford_D))
# Remove the "size" from scaffold names
Ashford_D$chr<-gsub("-size.*","",Ashford_D$chr)
head(Ashford_D)
summary(Ashford_D$D,na.rm = TRUE)
quantile(Ashford_D$D,probs = 0.05,na.rm = TRUE)
# D distribution
ggplot()+geom_histogram(data=Ashford_D,aes(D))+ggtitle("Ashford")

Terrace_D<-read.table("Terrace.D",header=FALSE)
# Latitude: 54.45
colnames(Terrace_D)<-c("chr","pos","NrSNPs","PropGdReads","D")
Terrace_D$D<-as.numeric(as.character(Terrace_D$D))
Terrace_D$pop<-rep("Terrace",nrow(Terrace_D))
# Remove the "size" from scaffold names
Terrace_D$chr<-gsub("-size.*","",Terrace_D$chr)
head(Terrace_D)
summary(Terrace_D$D,na.rm = TRUE)
quantile(Terrace_D$D,probs = 0.05,na.rm = TRUE)
# D distribution
ggplot()+geom_histogram(data=Terrace_D,aes(D))+ggtitle("Terrace")


Crested_Butte_D<-read.table("Crested_Butte.D",header=FALSE)
# Latitude: 38.9
colnames(Crested_Butte_D)<-c("chr","pos","NrSNPs","PropGdReads","D")
Crested_Butte_D$D<-as.numeric(as.character(Crested_Butte_D$D))
Crested_Butte_D$pop<-rep("Crested Butte",nrow(Crested_Butte_D))
# Remove the "size" from scaffold names
Crested_Butte_D$chr<-gsub("-size.*","",Crested_Butte_D$chr)
head(Crested_Butte_D)
summary(Crested_Butte_D$D,na.rm = TRUE)
quantile(Crested_Butte_D$D,probs = 0.05,na.rm = TRUE)
# D distribution
ggplot()+geom_histogram(data=Crested_Butte_D,aes(D))+ggtitle("Crested Butte")










# Convert the D data to linkage group data
# Scaffold sizes from genome freeze v.1.4
scaffdat<-read.table("montana_scaffsizes.tab",header=FALSE,sep=",")
colnames(scaffdat)<-c("chr","size")
nrow(scaffdat)
head(scaffdat[scaffdat$chr == "scaffold3615",])
# Load the v1 genome linkage groups
LGdat1<-read.table("markers_and_positions_v2.1.csv",header=TRUE,sep="\t")
colnames(LGdat1)<-c("scaff","linkage_gr","lkg_pos","lg","lgmrk")
head(LGdat1[LGdat1$scaff=="scaffold3615",])
head(LGdat1)
# Load the v1.4 genome linkage groups
LGdat1.4<-read.table("vennys_map_info_v1.4.txt",header=TRUE,sep="\t")
LGdat1.4$scaf_name<-gsub("-size.*","",LGdat1.4$scaf_name)
colnames(LGdat1.4)<-c("scaff","linkage_gr","size")
LGdat1.4$linkage_gr<-paste("Chr ",LGdat1.4$linkage_gr,sep="")
LGdat1.4$lkg_pos<-vector(length=nrow(LGdat1.4))
# Add scaffold positions on LGs to LGdat1.4
for(row in 1:nrow(LGdat1.4)){
  scf<-LGdat1.4$scaff[row]
  LGdat1.4$lkg_pos[row]<-LGdat1$lkg_pos[LGdat1$scaff==scf]
}
head(LGdat1.4[LGdat1.4$scaff=="scaffold3615",])
LGdat1.4<-LGdat1.4[order(LGdat1.4$linkage_gr,LGdat1.4$lkg_pos),]

## Add linkage groups and positions to the data set
head(LGdat1.4)
head(scaffdat)
head(Oulanka_D)
coord_converter<-function(scaff_data,LGdata,sizedata){
  scaff_data$lg<-vector(length=nrow(scaff_data))
  scaff_data$lgpos<-vector(length=nrow(scaff_data))
  scaff_data$scfsize<-vector(length=nrow(scaff_data))
  scaff_data_new<-data.frame()
  for(scaff in unique(scaff_data$chr)){
    scaff_data_sub<-scaff_data[scaff_data$chr==scaff,]
    lg<-as.character(LGdata$linkage_gr[LGdata$scaff==scaff])
    # If lg is null, the scaffold is not anchored to a linkage group.
    # Call it "Unplaced"
    if(identical(lg,character(0))){
      lg<-"Unplaced"
      #print(lg)
      scaff_data_sub$lg<-rep(lg,nrow(scaff_data_sub))
      scaff_data_sub$scfsize<-rep(scaffdat$size[scaffdat$chr==scaff],
                                 nrow(scaff_data_sub))
      scaff_data_new<-rbind(scaff_data_new,scaff_data_sub)
    }else{
      scaff_data_sub$lg<-rep(lg,nrow(scaff_data_sub))
      scaff_data_sub$lgpos<-rep(LGdata$lkg_pos[LGdata$scaff==scaff],
                              nrow(scaff_data_sub))
      scaff_data_sub$scfsize<-rep(LGdata$size[LGdata$scaff==scaff],
                              nrow(scaff_data_sub))
      scaff_data_new<-rbind(scaff_data_new,scaff_data_sub)
    }
  }
  # Subset to "placed" scaffolds
  scaff_data_new_pl<-scaff_data_new[grep("Unplaced",scaff_data_new$lg,invert=TRUE),]
  scaff_data_new_pl<-scaff_data_new_pl[
    order(scaff_data_new_pl$lg,scaff_data_new_pl$lgpos,scaff_data_new_pl$pos),]
  head(scaff_data_new_pl)
  ## Add an absoulte position for each linkage group
  scaff_data_new_pl$lg_win_pos<-vector(length=nrow(scaff_data_new_pl))
  scaff_data<-data.frame()
  for(LG in unique(scaff_data_new_pl$lg)){
    # Restart totpos at 0
    totpos<-0
    scaff_data_new_sub2<-scaff_data_new_pl[scaff_data_new_pl$lg == LG,]
    # Last scaffold length is 0 for the first scaffold on the LG
    last_scaff_l<-0
    # Looping through all the scaffolds on the LG ensures we can see the 
    # regions where we don't have any SNPs
    for(scaff in LGdata$scaff[LGdata$linkage_gr==LG]){
      scaff_data_new_sub3<-scaff_data_new_sub2[scaff_data_new_sub2$chr == scaff,]
      # Get scaffold size from the scaffdat table
      curr_scaff_l<-scaffdat$size[scaffdat$chr == scaff]
      #print(scaff)
      #print(curr_scaff_l)
      # Get a new position for each SNP by:
      # lg_snp_pos <- startpos + last_scaff_l + scaff_snp_pos
      scaff_data_new_sub3$lg_win_pos<-scaff_data_new_sub3$pos+totpos
      scaff_data<-rbind(scaff_data,scaff_data_new_sub3)
      last_scaff_l<-curr_scaff_l
      totpos <- totpos+last_scaff_l
    }
  }
  scaff_data$lg<-factor(scaff_data$lg)
  return(scaff_data)
}
Oulanka_D_LG<-coord_converter(Oulanka_D,LGdat1.4,scaffdat)
Korpilahti_D_LG<-coord_converter(Korpilahti_D,LGdat1.4,scaffdat)
Seward_D_LG<-coord_converter(Seward_D,LGdat1.4,scaffdat)
Crested_Butte_D_LG<-coord_converter(Crested_Butte_D,LGdat1.4,scaffdat)
Ashford_D_LG<-coord_converter(Ashford_D,LGdat1.4,scaffdat)
Terrace_D_LG<-coord_converter(Terrace_D,LGdat1.4,scaffdat)

allTajD<-rbind(Oulanka_D_LG,
           Korpilahti_D_LG,
           Seward_D_LG,
           Terrace_D_LG,
           Ashford_D_LG,
           Crested_Butte_D_LG)
allTajD$pop<-factor(allTajD$pop,levels=c("Oulanka","Korpilahti","Seward",
                                 "Terrace","Ashford","Crested Butte"))
allTajD$chr<-factor(allTajD$chr,levels=as.character(unique(allTajD$chr)))
  
D_plot<-ggplot()+
  geom_boxplot(data=allTajD,
            aes(pop,D))+
  ylab(expression(paste("Tajima's ",italic(D))))+
  xlab("")+
  ggtitle("B)")
D_plot+my.theme+theme(axis.text = element_text(size=10),
                       axis.title.y = element_text(size=14),
                      plot.title = element_text(size=14))
ggsave("~/Desktop/Data/montana_project/clinal_population_genomics/results/plots/montana_clines_TajD.png",device = "png",dpi=600)

# Plot TajD overall
D_plot<-ggplot()+
  geom_line(data=allTajD,
            aes(lg_win_pos,D,colour=chr,linetype=pop))+
  geom_hline(yintercept = 0,colour="grey")+
  scale_color_manual("",values=rep(c("black","red"),
                                   length.out=length(unique(all$chr))),
                     guide=FALSE)+
  scale_linetype_discrete("")+
  ylab(expression(paste("Tajima's ",italic(D))))+
  xlab("Position (kb)")+
  #ylim(-2,2)+
  facet_grid(lg~.)+
  my.theme+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=10),
        strip.text.y = element_text(size=10,angle=0),
        plot.title = element_text(hjust=0,size=10))
D_plot


# Get the locations of interesting genes
# See the script: montana_clinalpops_BSCAN.R
genedat

# Plot Tajima's D for the gene +/- 200kb
# Subset
allTajD[allTajD$lg=="Chr 5" &
          allTajD$lg_win_pos < 3600000 &
          allTajD$lg_win_pos > 3100000,]

# Collect the gene +/- 100 kb regions
tajD_data<-data.frame()
for(gene in genedat$geneid){
  LG<-genedat$lg[genedat$geneid==gene]
  SPOS<-genedat$g_lgspos[genedat$geneid==gene]
  EPOS<-genedat$g_lgepos[genedat$geneid==gene]
  sub<-allTajD[allTajD$lg == LG & 
                 allTajD$lg_win_pos > (SPOS - 200000) &
                 allTajD$lg_win_pos < (EPOS + 200000),]
  sub$geneid<-rep(paste(LG,": ",gene,sep=""),nrow(sub))
  tajD_data<-rbind(tajD_data,sub)
}
head(tajD_data)

D_plot<-ggplot()+
  geom_line(data=tajD_data,
                     aes(lg_win_pos/1000000,
                         D,colour=chr,linetype=pop))+
  geom_hline(yintercept = 0,colour="grey")+
  scale_color_manual("",values=rep(c("black","red"),
                                   length.out=length(unique(tajD_data$chr))),
                     guide=FALSE)+
  scale_linetype_discrete("")+
  ylab(expression(paste("Tajima's ",italic(D))))+
  xlab("Position (Mb)")+
  #ylim(-2,2)+
  facet_wrap(~geneid,scales = "free_x")+
  my.theme+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=10),
        strip.text.x = element_text(size=12,angle=0),
        plot.title = element_text(hjust=0,size=10),
        legend.position = "top",
        legend.text = element_text(size=10))
D_plot

head(genedat)
genedat2<-genedat
genedat2$geneid<-paste(genedat2$lg,": ",genedat2$geneid,sep="")
D_plot+geom_rect(data=genedat2,
                 aes(ymin=0.8,ymax=1,
                     xmin=(g_lgspos/1000000),
                     xmax=(g_lgepos/1000000)),fill="darkgreen")





