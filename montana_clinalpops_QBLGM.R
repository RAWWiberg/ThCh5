##
# D. montana project clinal populations
# QBGLM data analysis
# Last Modified: Feb 2017
##
# Clean environment
rm(list = ls(all = TRUE))

#
setwd("~/Desktop/Data/RData/montana_clines/")
dir<-"~/Desktop/Data/montana_project/clinal_population_genomics/results/syncfiles/"
# Load libraries
library(ggplot2)
library(qvalue)
library(reshape)
library(reshape2)
#citation("ggmap")
source("~/Desktop/Data/RData/RScripts/ggplot_theme.R")
options(scipen=10)
#------------------------------#
# Check coverage distributions
#------------------------------#
# ##########
# Read in data
# OULANKA
oulanka_cov<-read.table(
  "Oulanka_srt_covhist_genome.cov",
                        header = FALSE,sep="\t")
colnames(oulanka_cov)<-c("level","coverage","count","total","prop")
oulanka_cov$data<-rep("Oulanka",nrow(oulanka_cov))
# Remove anything with < perc10 and > 10xN in coverage (10 * 50 = 500)
oulanka_cov <- oulanka_cov[oulanka_cov$coverage < 500 &
                             oulanka_cov$coverage > 0,]
# Calculate the cumulative probability
oulanka_cov$cumprop <- cumsum(oulanka_cov$prop)
# What is the 90th percentile:
perc90<-max(oulanka_cov$coverage[oulanka_cov$cumprop < 0.9])
# What is the 10th percentile:
perc10<-max(oulanka_cov$coverage[oulanka_cov$cumprop < 0.1])
# What is the 5th percentile:
perc5<-max(oulanka_cov$coverage[oulanka_cov$cumprop < 0.05])
perc90
perc10
perc5

# KORPILAHTI
korpilahti_cov<-read.table(
  "Korpilahti_srt_covhist_genome.cov",
                        header = FALSE,sep="\t")
colnames(korpilahti_cov)<-c("level","coverage","count","total","prop")
korpilahti_cov$data<-rep("Korpilahti",nrow(korpilahti_cov))
# Remove anything with < 0 and > 10xN in coverage (10 * 50 = 500)
korpilahti_cov <- korpilahti_cov[korpilahti_cov$coverage < 500 &
                                korpilahti_cov$coverage > 0,]
# Calculate the cumulative probability
korpilahti_cov$cumprop <- cumsum(korpilahti_cov$prop)
# What is the 90th percentile:
perc90<-max(korpilahti_cov$coverage[korpilahti_cov$cumprop < 0.9])
# What is the 10th percentile:
perc10<-max(korpilahti_cov$coverage[korpilahti_cov$cumprop < 0.1])
# What is the 5th percentile:
perc5<-max(korpilahti_cov$coverage[korpilahti_cov$cumprop < 0.05])
perc90
perc10
perc5

# ASHFORD
ashford_cov<-read.table(
    "Ashford_srt_covhist_genome.cov",
                        header = FALSE,sep="\t")
colnames(ashford_cov)<-c("level","coverage","count","total","prop")
ashford_cov$data<-rep("Ashford",nrow(ashford_cov))
# Remove anything with < 0 and > 10xN in coverage (10 * 50 = 500)
ashford_cov <- ashford_cov[ashford_cov$coverage < 500 &
                             ashford_cov$coverage > 0,]
# Calculate the cumulative probability
ashford_cov$cumprop <- cumsum(ashford_cov$prop)
# What is the 90th percentile:
perc90<-max(ashford_cov$coverage[ashford_cov$cumprop < 0.9])
# What is the 10th percentile:
perc10<-max(ashford_cov$coverage[ashford_cov$cumprop < 0.1])
# What is the 5th percentile:
perc5<-max(ashford_cov$coverage[ashford_cov$cumprop < 0.05])
perc90
perc10
perc5

# CRESTED BUTTE
crested_butte_cov<-read.table(
  "Crested_Butte_srt_covhist_genome.cov",
                        header = FALSE,sep="\t")
colnames(crested_butte_cov)<-c("level","coverage","count","total","prop")
crested_butte_cov$data<-rep("Crested Butte",nrow(crested_butte_cov))
# Remove anything with < 0 and > 10xN in coverage (10 * 50 = 500)
crested_butte_cov <- crested_butte_cov[crested_butte_cov$coverage < 500 &
                                         crested_butte_cov$coverage > 0,]
# Calculate the cumulative probability
crested_butte_cov$cumprop <- cumsum(crested_butte_cov$prop)
# What is the 90th percentile:
perc90<-max(crested_butte_cov$coverage[crested_butte_cov$cumprop < 0.9])
# What is the 10th percentile:
perc10<-max(crested_butte_cov$coverage[crested_butte_cov$cumprop < 0.1])
# What is the 5th percentile:
perc5<-max(crested_butte_cov$coverage[crested_butte_cov$cumprop < 0.05])
perc90
perc10
perc5

# TERRACE
terrace_cov<-read.table(
  "Terrace_srt_covhist_genome.cov",
                        header = FALSE,sep="\t")
colnames(terrace_cov)<-c("level","coverage","count","total","prop")
terrace_cov$data<-rep("Terrace",nrow(terrace_cov))
# Remove anything with < 0 and > 10xN in coverage (10 * 50 = 500)
terrace_cov <- terrace_cov[terrace_cov$coverage < 500 &
                             terrace_cov$coverage > 0,]
# Calculate the cumulative probability
terrace_cov$cumprop <- cumsum(terrace_cov$prop)
# What is the 90th percentile:
perc90<-max(terrace_cov$coverage[terrace_cov$cumprop < 0.9])
# What is the 10th percentile:
perc10<-max(terrace_cov$coverage[terrace_cov$cumprop < 0.1])
# What is the 5th percentile:
perc5<-max(terrace_cov$coverage[terrace_cov$cumprop < 0.05])
perc90
perc10
perc5

# SEWARD
seward_cov<-read.table(
  "Seward_subs_srt_covhist_genome.cov",
                        header = FALSE,sep="\t")
colnames(seward_cov)<-c("level","coverage","count","total","prop")
seward_cov$data<-rep("Seward",nrow(seward_cov))
# Remove anything with < 0 and > 10xN in coverage (10 * 50 = 500)
seward_cov <- seward_cov[seward_cov$coverage < 500 &
                           seward_cov$coverage > 0,]
# Calculate the cumulative probability
seward_cov$cumprop <- cumsum(seward_cov$prop)
# What is the 90th percentile:
perc90<-max(seward_cov$coverage[seward_cov$cumprop < 0.9])
# What is the 10th percentile:
perc10<-max(seward_cov$coverage[seward_cov$cumprop < 0.1])
# What is the 5th percentile:
perc5<-max(seward_cov$coverage[seward_cov$cumprop < 0.05])
perc90
perc10
perc5


# PLOTTING
all_cov<-rbind(oulanka_cov,
               korpilahti_cov,
               crested_butte_cov,
               terrace_cov,
               ashford_cov,
               seward_cov)

# Plot the cumulative proportion coverage
cov_hist<-ggplot()+
  geom_line(data=all_cov[all_cov$coverage > 0,],
           aes(x=coverage,y=cumprop),
           stat="identity")+
  #  scale_x_continuous(breaks=seq(0,400,10))+
  ylab("Proportion of Sites")+
  xlab("Coverage")+
  facet_grid(data~.)

cov_hist+my.theme+theme(axis.text.x = element_text(size = 12,
                                                   angle=45,
                                                   vjust = 1,
                                                   hjust = 1),
                        axis.text.y = element_text(size = 12),
                        strip.text.y=element_text(size=12,angle=0))

# Plot the coverage histograms
cov_hist<-ggplot()+
  geom_line(data=all_cov[all_cov$coverage > 0,],
            aes(x=coverage,y=prop),
            stat="identity")+
  #  scale_x_continuous(breaks=seq(0,400,10))+
  ylab("Proportion of Sites")+
  xlab("Coverage")+
  facet_grid(data~.)

cov_hist+my.theme+theme(axis.text.x = element_text(size = 12,
                                                   angle=45,
                                                   vjust = 1,
                                                   hjust = 1),
                        axis.text.y = element_text(size = 12),
                        strip.text.y=element_text(size=12,angle=0))


# Get percentiles for the aggregated distribution
aggr_cov<-aggregate(all_cov$count,list(coverage=all_cov$coverage),sum)
aggr_cov$prop <- aggr_cov$x/sum(as.numeric(aggr_cov$x))
aggr_cov$data<-rep("Aggregate",nrow(aggr_cov))

colnames(aggr_cov)<-c("coverage","count","prop","data")

aggr_cov$prop2<-aggr_cov$count/sum(as.numeric(aggr_cov$count))
aggr_cov$cumprop <- cumsum(aggr_cov$prop2)

# 10th percentile
max(aggr_cov$coverage[aggr_cov$cumprop < 0.1])

# 90th percentile
max(aggr_cov$coverage[aggr_cov$cumprop < 0.9])
#
vlin<-data.frame(xint1=max(aggr_cov$coverage[aggr_cov$cumprop < 0.1]),
                 xint2=max(aggr_cov$coverage[aggr_cov$cumprop < 0.9]),
                 data="Aggregate")

head(aggr_cov[c(1,2,3,4)])
head(all_cov[c(2,3,5,6)])
all_cov <- rbind(all_cov[c(2,3,5,6)],
                 aggr_cov[c(1,2,3,4)])
all_cov$data<-factor(all_cov$data,
                     levels=c("Oulanka","Korpilahti","Seward",
                              "Terrace","Ashford","Crested Butte",
                              "Aggregate"))

cov_hist<-ggplot()+
  geom_bar(data=all_cov[all_cov$coverage > 0,],
           aes(x=coverage,y=prop),
           stat="identity")+
  scale_y_continuous(breaks=c(0,0.010,0.020))+
  scale_x_continuous(breaks=c(0,50,100,150),limits=c(0,155))+
  geom_vline(data=vlin,aes(xintercept = xint1),colour="red",linetype="dashed")+
  geom_vline(data=vlin,aes(xintercept = xint2),colour="red",linetype="dashed")+
  ylab("Proportion of Sites")+
  xlab("Coverage")+
  facet_grid(data~.)

cov_hist+my.theme+theme(axis.text.x = element_text(size = 12,
                                                   angle=45,
                                                   vjust = 1,
                                                   hjust = 1),
                        axis.title = element_text(size=10),
                        axis.text.y = element_text(size = 10),
                        strip.text.y=element_text(size=12,angle=0))


# MEAN COVERAGE PER SCAFFOLD
pops<-c("Ashford","Crested_Butte","Seward","Terrace","Korpilahti","Oulanka")
all_covperscaff<-data.frame()
# Read in coverage data
for(pop in pops){
 covperscaff<-read.csv(
   paste(pop,"_srt_mcov_perscaff.tab",sep=""),sep="\t",header=FALSE)
 colnames(covperscaff)<-c("scaff","cov")
 covperscaff$pop<-rep(pop,nrow(covperscaff))
 print(covperscaff[covperscaff$scaff == "genome",])
 covperscaff<-covperscaff[covperscaff$scaff != "genome",]
 covperscaff$scaffsize<-vector(length = nrow(covperscaff))
 covperscaff$scaffsize <- covperscaff$scaffsize<-as.data.frame(
   do.call("rbind",
           strsplit(as.character(covperscaff$scaff),"-")))[,2]
 covperscaff$scaffsize<-gsub("size","",covperscaff$scaffsize)
 covperscaff$scaff<-gsub("-size.*","",covperscaff$scaff)
 covperscaff$pop<-rep(pop,nrow(covperscaff))
 all_covperscaff<-rbind(all_covperscaff,covperscaff)
}
head(all_covperscaff)
str(all_covperscaff)
all_covperscaff$size<-as.numeric(all_covperscaff$scaffsize)
all_covperscaff[is.na(all_covperscaff$size),]

# Read in X scaff data
xscaffs<-read.table("Xonly.csv",header=TRUE,sep=",")
colnames(xscaffs)<-c("scaff","pos")
head(xscaffs)
xscaff_nms<-unique(xscaffs$scaff)


# Add chromosome type to data
all_covperscaff_aut<-all_covperscaff[
  which(!(all_covperscaff$scaff %in% xscaff_nms)),]
all_covperscaff_aut$type<-rep("Aut",nrow(all_covperscaff_aut))
all_covperscaff_x<-all_covperscaff[
  which(all_covperscaff$scaff %in% xscaff_nms),]
all_covperscaff_x$type<-rep("X",nrow(all_covperscaff_x))

all_covperscaff<-rbind(all_covperscaff_aut,all_covperscaff_x)
str(all_covperscaff)

# Remove anything with cov > maxcov and < mincov
maxcov <- 500
mincov <- 0
minscaffsize <- 10000
all_covperscaff_sub<-all_covperscaff[
  all_covperscaff$cov <= maxcov & 
    all_covperscaff$cov > mincov &
    all_covperscaff$scaffsize > minsize,]

str(all_covperscaff_sub)
tail(all_covperscaff_sub)

ggplot()+
  geom_point(data=all_covperscaff_sub,
             aes(size/10000,cov,colour=type),
             alpha=1/2)+
  facet_grid(pop~.)

tapply(all_covperscaff_sub$cov,
       INDEX = list(all_covperscaff_sub$type,
                    all_covperscaff_sub$pop),mean)

ggplot()+
  geom_histogram(data=all_covperscaff_sub,
    aes(cov,fill=type),
    position="identity",alpha=1/2)+
  facet_wrap(pop~type,ncol=2,scales="free_y")



# ##########



#------------------------------#
# Load QBGLM data
#------------------------------#
v<-"altpc1pc2"
test<-"qbglm"
scale<-"neff"
minscaffsize<-10000
qbglm_dat<-read.table(paste("montana_clines",
                            "_",v,
                            "_",scale,
                            "_",test,
                            "_minc10_cov37-106_regrouped.rout",
                            sep=""),
                  header = TRUE)
qbglm_dat$chr<-gsub("-size.*","",qbglm_dat$chr)
head(qbglm_dat)

#------------------------------#
# Load the scaffold size data  #
# and the linkage group data   #
#------------------------------#
# Scaffold sizes from genome freeze v.1.4
scaffdat<-read.table("montana_scaffsizes.tab",header=FALSE,sep=",")
colnames(scaffdat)<-c("chr","size")
nrow(scaffdat)
# Load the v1 genome linkage groups
LGdat1<-read.table("markers_and_positions_v2.1.csv",header=TRUE,sep="\t")
colnames(LGdat1)<-c("scaff","linkage_gr","lkg_pos","lg","lgmrk")
# Load the v1.4 genome linkage groups
LGdat1.4<-read.table("vennys_map_info_v1.4.txt",header=TRUE,sep="\t")
LGdat1.4$scaf_name<-gsub("-size.*","",LGdat1.4$scaf_name)
colnames(LGdat1.4)<-c("scaff","linkage_gr","size")
LGdat1.4$linkage_gr<-paste("Chr ",LGdat1.4$linkage_gr,sep="")
LGdat1.4$lkg_pos<-vector(length=nrow(LGdat1.4))
head(LGdat1.4)

# Add scaffold positions on LGs to LGdat1.4
for(row in 1:nrow(LGdat1.4)){
  scf<-LGdat1.4$scaff[row]
  LGdat1.4$lkg_pos[row]<-LGdat1$lkg_pos[LGdat1$scaff==scf]
}
LGdat1.4<-LGdat1.4[order(LGdat1.4$linkage_gr,LGdat1.4$lkg_pos),]
unique(LGdat1.4$linkage_gr)
# Get names of scaffolds mapped to the X-chromosome (Venera/Darren)
xscaffs<-LGdat1.4[LGdat1.4$linkage_gr=="Chr X",c(1,4)]
colnames(xscaffs)<-c("scaff","pos")
head(xscaffs)
xscaff_nms<-unique(xscaffs$scaff)
# Get names of scaffolds mapped to the Autosomes (Venera/Darren)
mapped_aut<-LGdat1.4[grep("Chr X",LGdat1.4$linkage_gr,invert=TRUE),c(1,4)]
colnames(mapped_aut)<-c("scaff","pos")
head(mapped_aut)
mapped_aut_nms<-unique(mapped_aut$scaff)

# Histogram of scaffold sizes
scaff_l_hist<-ggplot()+
  geom_histogram(data=qbglm_dat[
    qbglm_dat$scaff %in% xscaff_nms,],aes(scaffsize/1000))+
  xlab("Scaffold Size (kb)")+
  ylab("Count")
scaff_l_hist+my.theme+theme(axis.text.y = element_text(size=10),
                           axis.text.x = element_text(size=10),
                           axis.title = element_text(size=10))
summary(qbglm_dat$scaffsize/1000)

# Remove the scaffolds that are smaller than minscaffsize
qbglm_dat<-qbglm_dat[qbglm_dat$scaffsize>=minscaffsize,]
# Subset to the scaffolds that are mapped to autosomal linkage groups
qbglm_dat_aut<-qbglm_dat[which(qbglm_dat$chr %in% mapped_aut_nms),]
nrow(qbglm_dat_aut)
# Subset to the scaffolds that are mapped to X chromosome linkage groups
qbglm_dat_x<-qbglm_dat[which((qbglm_dat$chr %in% xscaff_nms)),]
nrow(qbglm_dat_x)

nrow(qbglm_dat)
nrow(qbglm_dat_aut)
nrow(qbglm_dat_x)
# Write the .sync file for the X and Autosome sscaffolds which
# are longer than minscaffsize
head(qbglm_dat_x)
write.table(qbglm_dat_aut[,c(16,2,3,4,5,6,7,8,9)],
            paste(dir,"montana_clines_AutLG",
                  "Mins",minscaffsize,".sync",sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

write.table(qbglm_dat_x[,c(16,2,3,4,5,6,7,8,9)],
            paste(dir,"montana_clines_XLG",
                  "Mins",minscaffsize,".sync",sep=""),
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

length(unique(qbglm_dat$chr))
nrow(qbglm_dat)
head(qbglm_dat_aut)

# Draw q-q plots
# Alt
#####
# Distribution of pvalues
pval_hist<-ggplot()+
  geom_histogram(data=qbglm_dat_aut,aes(alt_p))+
  xlab("P value")+
  ylab("Count")
pval_hist<-pval_hist+my.theme+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=10))
#pval_hist+theme(axis.text.y=element_text(size=10))
## Q-Q plot of p-values
o<-sort(qbglm_dat$alt_p,decreasing=FALSE)
pdat<-data.frame(o=o,
                 e=(1:length(o))/length(o))
q_q_plot<-ggplot()+
  geom_point(data=pdat,aes(-log10(e),-log10(o)),size=0.2)+
  geom_line(data=pdat,aes(-log10(e),-log10(e)),
            linetype="dashed",colour="red")+
  xlab("Expected -log10(p-values)")+
  ylab("Observed -log10(p-values)")+
  my.theme+
  ggtitle("A)")+
  theme(axis.text = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12))
#q_q_plot

# Make Q-Q plot with p-value distribution inset.
xleft   = 0.05
xright  = 0.35
ybottom = 0.55
ytop    = 0.95 

# Calculate position in plot1 coordinates
# Extract x and y values from q_q_plot
l1 = ggplot_build(q_q_plot)
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
g2 = ggplotGrob(pval_hist)
plot3 = q_q_plot + annotation_custom(grob = g2, 
                                     xmin=xmin, xmax=xmax, 
                                     ymin=ymin, ymax=ymax)
plot(plot3)

ggsave("~/Desktop/Data/montana_project/clinal_population_genomics/results/plots/qbglm_qqplot_alt",device = "png",dpi=600)
#####
# PC1
#####
# Distribution of pvalues
pval_hist<-ggplot()+
  geom_histogram(data=qbglm_dat_aut,aes(pc1_p))+
  xlab("P value")+
  ylab("Count")
pval_hist<-pval_hist+my.theme+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=10))
#pval_hist+theme(axis.text.y=element_text(size=10))
## Q-Q plot of p-values
o<-sort(qbglm_dat$pc1_p,decreasing=FALSE)
pdat<-data.frame(o=o,
                 e=(1:length(o))/length(o))
q_q_plot<-ggplot()+
  geom_point(data=pdat,aes(-log10(e),-log10(o)),size=0.2)+
  geom_line(data=pdat,aes(-log10(e),-log10(e)),
            linetype="dashed",colour="red")+
  xlab("Expected -log10(p-values)")+
  ylab("Observed -log10(p-values)")+
  my.theme+
  ggtitle("B)")+
  theme(axis.text = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12))
#q_q_plot

# Make Q-Q plot with p-value distribution inset.
xleft   = 0.05
xright  = 0.35
ybottom = 0.55
ytop    = 0.95 

# Calculate position in plot1 coordinates
# Extract x and y values from q_q_plot
l1 = ggplot_build(q_q_plot)
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
g2 = ggplotGrob(pval_hist)
plot3 = q_q_plot + annotation_custom(grob = g2, 
                                     xmin=xmin, xmax=xmax, 
                                     ymin=ymin, ymax=ymax)
plot(plot3)

ggsave("~/Desktop/Data/montana_project/clinal_population_genomics/results/plots/qbglm_qqplot_pc1",device = "png",dpi=600)

#####
# PC2
#####
# Distribution of pvalues
pval_hist<-ggplot()+
  geom_histogram(data=qbglm_dat_aut,aes(pc2_p))+
  xlab("P value")+
  ylab("Count")
pval_hist<-pval_hist+my.theme+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=10))
#pval_hist+theme(axis.text.y=element_text(size=10))
## Q-Q plot of p-values
o<-sort(qbglm_dat$pc2_p,decreasing=FALSE)
pdat<-data.frame(o=o,
                 e=(1:length(o))/length(o))
q_q_plot<-ggplot()+
  geom_point(data=pdat,aes(-log10(e),-log10(o)),size=0.2)+
  geom_line(data=pdat,aes(-log10(e),-log10(e)),
            linetype="dashed",colour="red")+
  xlab("Expected -log10(p-values)")+
  ylab("Observed -log10(p-values)")+
  my.theme+
  ggtitle("C)")+
  theme(axis.text = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12))
#q_q_plot

# Make Q-Q plot with p-value distribution inset.
xleft   = 0.05
xright  = 0.35
ybottom = 0.55
ytop    = 0.95 

# Calculate position in plot1 coordinates
# Extract x and y values from q_q_plot
l1 = ggplot_build(q_q_plot)
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
g2 = ggplotGrob(pval_hist)
plot3 = q_q_plot + annotation_custom(grob = g2, 
                                  xmin=xmin, xmax=xmax, 
                                  ymin=ymin, ymax=ymax)
plot(plot3)
ggsave("~/Desktop/Data/montana_project/clinal_population_genomics/results/plots/qbglm_qqplot_pc2",device = "png",dpi=600)
#####



## Add linkage groups and positions to the data set
# Also subsets to the "placed" scaffolds
# ####
head(qbglm_dat_aut)
qbglm_dat_aut$lg<-vector(length=nrow(qbglm_dat_aut))
qbglm_dat_aut$lgpos<-vector(length=nrow(qbglm_dat_aut))
qbglm_dat_aut$scfsize<-vector(length=nrow(qbglm_dat_aut))
qbglm_dat_aut_new<-data.frame()
for(scaff in unique(qbglm_dat_aut$chr)){
  qbglm_dat_aut_sub<-qbglm_dat_aut[qbglm_dat_aut$chr==scaff,]
  lg<-as.character(LGdat1.4$linkage_gr[LGdat1.4$scaff==scaff])
  # If lg is null, the scaffold is not anchored to a linkage group.
  # Call it "Unplaced"
  if(identical(lg,character(0))){
    lg<-"Unplaced"
    #print(lg)
    qbglm_dat_aut_sub$lg<-rep(lg,nrow(qbglm_dat_aut_sub))
    qbglm_dat_aut_sub$scfsize<-rep(scaffdat$size[scaffdat$chr==scaff],
                              nrow(qbglm_dat_aut_sub))
    qbglm_dat_aut_new<-rbind(qbglm_dat_aut_new,qbglm_dat_aut_sub)
  }else{
    qbglm_dat_aut_sub$lg<-rep(lg,nrow(qbglm_dat_aut_sub))
    qbglm_dat_aut_sub$lgpos<-rep(LGdat1.4$lkg_pos[LGdat1.4$scaff==scaff],
                            nrow(qbglm_dat_aut_sub))
    qbglm_dat_aut_sub$scfsize<-rep(LGdat1.4$size[LGdat1.4$scaff==scaff],
                              nrow(qbglm_dat_aut_sub))
    qbglm_dat_aut_new<-rbind(qbglm_dat_aut_new,qbglm_dat_aut_sub)
  }
}
rm(qbglm_dat_aut_sub)
head(qbglm_dat_aut_new)
summary(qbglm_dat_aut_new)

## Subset to only placed scaffolds.
# Order the placed dataset by LG, position on LG, position in scaffold
qbglm_dat_aut_new<-qbglm_dat_aut_new[
  grep("Unplaced",qbglm_dat_aut_new$lg,invert=TRUE),]
qbglm_dat_aut_new<-qbglm_dat_aut_new[
  order(qbglm_dat_aut_new$lg,qbglm_dat_aut_new$lgpos,qbglm_dat_aut_new$pos),]
head(qbglm_dat_aut_new)
qbglm_dat_aut_new$lg<-factor(qbglm_dat_aut_new$lg)

nrow(qbglm_dat_aut_new)
summary(qbglm_dat_aut_new)

## Add an absolute position for each linkage group
qbglm_dat_aut_new$lg_snp_pos<-vector(length=nrow(qbglm_dat_aut_new))
qbglm_dat_aut<-data.frame()
for(LG in unique(qbglm_dat_aut_new$lg)){
  # Restart totpos at 0
  totpos<-0
  qbglm_dat_aut_new_sub2<-qbglm_dat_aut_new[qbglm_dat_aut_new$lg == LG,]
  # Last scaffold length is 0 for the first scaffold on the LG
  last_scaff_l<-0
  # Looping through all the scaffolds on the LG ensures we can see the 
  # regions where we don't have any SNPs
  for(scaff in LGdat1.4$scaff[LGdat1.4$linkage_gr==LG]){
    qbglm_dat_aut_new_sub3<-qbglm_dat_aut_new_sub2[qbglm_dat_aut_new_sub2$chr == scaff,]
    # Get scaffold size from the scaffdat table
    curr_scaff_l<-scaffdat$size[scaffdat$chr == scaff]
    #print(scaff)
    #print(curr_scaff_l)
    # Get a new position for each SNP by:
    # lg_snp_pos <- startpos + last_scaff_l + scaff_snp_pos
    qbglm_dat_aut_new_sub3$lg_snp_pos<-qbglm_dat_aut_new_sub3$pos+totpos
    qbglm_dat_aut<-rbind(qbglm_dat_aut,qbglm_dat_aut_new_sub3)
    last_scaff_l<-curr_scaff_l
    totpos <- totpos+last_scaff_l
  }
}
rm(qbglm_dat_aut_new_sub2)
rm(qbglm_dat_aut_new_sub3)
qbglm_dat_aut$lg<-factor(qbglm_dat_aut$lg)
# ####


# Print a subset of SNPs to a file for plotting
# ####
subs<-qbglm_dat[sample(seq(1,nrow(qbglm_dat_aut)),size = 100,replace = TRUE),]
head(subs)
write.table(subs[,1:9],
            file = paste(dir,"qbglm_",v,
                         "_",scale,
                         "_",test,"_subs.sync",sep=""),
            row.names=FALSE,col.names=FALSE,quote=FALSE,sep = "\t")
# ####

# Convert p-values to to q-values
write.table(qbglm_dat_aut$pc1_p,
            file = paste("qbglm_",v,
                         "_",scale,
                         "_",test,"_pvals.list",sep=""),
            row.names=FALSE,col.names=FALSE,quote=FALSE)
pvals<-scan(paste("qbglm_",v,"_",scale,"_",test,"_pvals.list",sep=""))

qobj<-qvalue(pvals,pi0.method = "smoother")
plot(qobj)
qsummary(qobj)
qbglm_dat_aut$pc1_qval<-qobj$qvalues
min(qbglm_dat_aut$pc1_qval)
min(qbglm_dat_aut$pc1_p)

# How many SNPS
nrow(qbglm_dat_aut)
# How many with q-vals < 0.05
nrow(qbglm_dat_aut[qbglm_dat_aut$qval < 0.05,])
nrow(qbglm_dat_aut[qbglm_dat_aut$qval < 0.01,])
min(qbglm_dat_aut$qval)
# How many below bonferroni 
bonf<-0.05/nrow(qbglm_dat_aut)
nrow(qbglm_dat_aut[qbglm_dat_aut$pval < bonf,])

# pvalues vs. qvalues
qval_v_pval<-ggplot()+
  geom_point(data=qbglm_dat_aut,aes(pc1_p,pc1_qval),
             alpha=1/4,size=0.1)+
  geom_line(aes(seq(0,1,0.001),seq(0,1,0.001)),linetype="dashed",colour="red")+
  xlim(0,1)+
  ylim(0,1)+
  xlab("p-values")+
  ylab("q-values")+
  ggtitle("D)")
qval_v_pval+my.theme+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12))
ggsave("~/Desktop/Data/montana_project/clinal_population_genomics/results/plots/qbglm_q_vs_p_pc1.png",device = "png",dpi=600)

# Plot p-values on linkage groups
# p-values
head(qbglm_dat_aut)
manhplot<-ggplot()+
  geom_point(data=qbglm_dat_aut,
             aes(lg_snp_pos/1000000,-log10(pc1_p),colour=chr),
             size=0.5) +
  xlab("Position (Mb)")+
  ylab("-log10(p-value)")+
  ggtitle("Results for pc1")+
  geom_hline(yintercept = -log10(0.05),colour="red",linetype="dashed")+
  geom_hline(yintercept = -log10(0.01),colour="blue",linetype="dashed")+
#  geom_hline(yintercept = -log10(bonf),colour="green",linetype="dashed")+
  scale_colour_manual(guide=FALSE,
                      values=rep(c("grey","black"),
                                 length.out=length(unique(qbglm_dat_aut$chr))))+
  facet_grid(lg~.)
manhplot+my.theme+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(size=12),
        strip.text.y = element_text(angle=0,size=10))
ggsave("~/Desktop/Data/montana_project/clinal_population_genomics/results/plots/qbglm_pc1_manhplot_pvals.png",device = "png",dpi=600)

# Q-values
manhplot<-ggplot()+
  geom_point(data=qbglm_dat_aut,
             aes(lg_snp_pos/1000000,-log10(pc1_qval),colour=chr),
             size=0.5) +
  xlab("Position (Mb)")+
  ylab("-log10(q-value)")+
  ggtitle("Results for pc1")+
  geom_hline(yintercept = -log10(0.05),colour="red",linetype="dashed")+
  geom_hline(yintercept = -log10(0.01),colour="blue",linetype="dashed")+
  #  geom_hline(yintercept = -log10(bonf),colour="green",linetype="dashed")+
  scale_colour_manual(guide=FALSE,
                      values=rep(c("grey","black"),
                                 length.out=length(unique(qbglm_dat_aut$chr))))+
  facet_grid(lg~.)

manhplot+my.theme+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        plot.title = element_text(size=12),
        strip.text.y = element_text(angle=0,size=10))
ggsave("~/Desktop/Data/montana_project/clinal_population_genomics/results/plots/qbglm_pc1_manhplot_qvals.png",device = "png",dpi=600)


