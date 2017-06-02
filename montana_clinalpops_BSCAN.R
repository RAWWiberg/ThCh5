##
# D. montana project clinal populations
# BSCAN data analysis
# Last Modified: Feb 2017
##
# Clean environment
rm(list = ls(all = TRUE))

pc<-"pc1"
#
setwd(paste("~/Desktop/Data/RData/montana_clines/bscan_",pc,sep=""))
# Load libraries
library(ggplot2)
library(qvalue)
library(reshape)
library(reshape2)
#citation("ggmap")
source("~/Desktop/Data/RData/RScripts/ggplot_theme.R")
options(scipen=10)
#Define Functions
l_uniq<-function(x){length(unique(x))}
#------------------------------#
# Load BAYESCENV data
#------------------------------#
bscanout_dat<-read.table("montana_clines_Aut_all_bscanin_fst_sub.txt",header = TRUE)
colnames(bscanout_dat)<-c("snp","g_qval","g")
bscansnps_dat<-read.table("montana_clines_Aut_all_snps.txt",header=FALSE)
colnames(bscansnps_dat)<-c("chr","pos","snp")

# Combine the data
bscandat<-cbind(bscansnps_dat[,c(1,2)],bscanout_dat)
head(bscandat)
# How many SNPs
nrow(bscandat)


#------------------------------#
# Load the scaffold size data  #
# and the linkage group data   #
#------------------------------#
# Scaffold sizes from genome freeze v.1.4
scaffdat<-read.table("../montana_scaffsizes.tab",header=FALSE,sep=",")
colnames(scaffdat)<-c("chr","size")
nrow(scaffdat)
head(scaffdat[scaffdat$chr == "scaffold3615",])
# Load the v1 genome linkage groups
LGdat1<-read.table("../markers_and_positions_v2.1.csv",header=TRUE,sep="\t")
colnames(LGdat1)<-c("scaff","linkage_gr","lkg_pos","lg","lgmrk")
head(LGdat1[LGdat1$scaff=="scaffold3615",])
head(LGdat1)

# Load the v1.4 genome linkage groups
LGdat1.4<-read.table("../vennys_map_info_v1.4.txt",header=TRUE,sep="\t")
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
head(LGdat1.4)

#---------------------------------#
# BAYSCENV SUMMARIES AND ANALYSIS #
#---------------------------------#

# What number of scaffolds are < 10kb
length(scaffdat$size[scaffdat$size>10000])/length(scaffdat$size)
length(scaffdat$size[scaffdat$size<10000])/length(scaffdat$size)

sum(scaffdat$size[scaffdat$size>10000])/sum(scaffdat$size)

# Subset to scaffolds > 10kb
scaffdat_sub<-scaffdat[scaffdat$size>10000,]
# What number of scaffolds are not on the linkage groups
length(which(scaffdat_sub$chr %in% LGdat1.4$scaff))
length(which(!(scaffdat_sub$chr %in% LGdat1.4$scaff)))

# what is the total length of scaffolds on linkage groups
totl<-sum(scaffdat_sub$size)
sum(scaffdat_sub$size[which(scaffdat_sub$chr %in% LGdat1.4$scaff)])/totl
sum(scaffdat_sub$size[which(scaffdat_sub$chr %in% LGdat1.4$scaff)])/1000000
sum(scaffdat_sub$size[which(!(scaffdat_sub$chr %in% LGdat1.4$scaff))])/totl
sum(scaffdat_sub$size[which(!(scaffdat_sub$chr %in% LGdat1.4$scaff))])/1000000
# How long is each lg
tapply(LGdat1.4$size,INDEX = list(LGdat1.4$linkage_gr),sum)/1000000

nrow(bscandat)
length(unique(bscandat$chr))
length(unique(LGdat1.4$scaff))
length(which(bscandat$chr %in% LGdat1.4$scaff))

## Add linkage groups and positions to the data set
head(bscandat)
bscandat$lg<-vector(length=nrow(bscandat))
bscandat$lgpos<-vector(length=nrow(bscandat))
bscandat$scfsize<-vector(length=nrow(bscandat))
bscandat_new<-data.frame()
for(scaff in unique(bscandat$chr)){
  bscandat_sub<-bscandat[bscandat$chr==scaff,]
  lg<-as.character(LGdat1.4$linkage_gr[LGdat1.4$scaff==scaff])
  # If lg is null, the scaffold is not anchored to a linkage group.
  # Call it "Unplaced"
  if(identical(lg,character(0))){
    lg<-"Unplaced"
    #print(lg)
    bscandat_sub$lg<-rep(lg,nrow(bscandat_sub))
    bscandat_sub$scfsize<-rep(scaffdat$size[scaffdat$chr==scaff],
                              nrow(bscandat_sub))
    bscandat_new<-rbind(bscandat_new,bscandat_sub)
  }else{
    bscandat_sub$lg<-rep(lg,nrow(bscandat_sub))
    bscandat_sub$lgpos<-rep(LGdat1.4$lkg_pos[LGdat1.4$scaff==scaff],
                          nrow(bscandat_sub))
    bscandat_sub$scfsize<-rep(LGdat1.4$size[LGdat1.4$scaff==scaff],
                              nrow(bscandat_sub))
    bscandat_new<-rbind(bscandat_new,bscandat_sub)
  }
}
rm(bscandat_sub)
head(bscandat_new)
summary(bscandat_new)

## Order scaffolds wich are not on linkage groups by size
bscandat_new_unpl<-bscandat_new[grep("Unplaced",bscandat_new$lg),]
bscandat_new_unpl<-bscandat_new_unpl[
  order(-bscandat_new_unpl$scfsize,bscandat_new_unpl$chr,bscandat_new_unpl$pos),]
head(bscandat_new_unpl)
bscandat_new_unpl$lg<-factor(bscandat_new_unpl$lg)

nrow(bscandat_new_unpl)
summary(bscandat_new_unpl)

## Order the placed dataset by LG, position on LG, position in scaffold
bscandat_new_pl<-bscandat_new[grep("Unplaced",bscandat_new$lg,invert=TRUE),]
bscandat_new_pl<-bscandat_new_pl[
  order(bscandat_new_pl$lg,bscandat_new_pl$lgpos,bscandat_new_pl$pos),]
head(bscandat_new_pl)
bscandat_new_pl$lg<-factor(bscandat_new_pl$lg)

nrow(bscandat_new_pl)
summary(bscandat_new_pl)

# What number of scaffolds are "placed"
tapply(bscandat_new_pl$scfsize,
       INDEX = list(bscandat_new_pl$lg),FUN=l_uniq)

# How many scaffolds are "unplaced"
tapply(bscandat_new_unpl$scfsize,
       INDEX = list(bscandat_new_unpl$lg),FUN=l_uniq)

## Add an absoulte position for each linkage group
bscandat_new_pl$lg_snp_pos<-vector(length=nrow(bscandat_new_pl))
bscandat<-data.frame()
for(LG in unique(bscandat_new$lg)){
  # Restart totpos at 0
  totpos<-0
  bscandat_new_sub2<-bscandat_new_pl[bscandat_new_pl$lg == LG,]
  # Last scaffold length is 0 for the first scaffold on the LG
  last_scaff_l<-0
  # Looping through all the scaffolds on the LG ensures we can see the 
  # regions where we don't have any SNPs
  for(scaff in LGdat1.4$scaff[LGdat1.4$linkage_gr==LG]){
    bscandat_new_sub3<-bscandat_new_sub2[bscandat_new_sub2$chr == scaff,]
    # Get scaffold size from the scaffdat table
    curr_scaff_l<-scaffdat$size[scaffdat$chr == scaff]
    #print(scaff)
    #print(curr_scaff_l)
    # Get a new position for each SNP by:
    # lg_snp_pos <- startpos + last_scaff_l + scaff_snp_pos
    bscandat_new_sub3$lg_snp_pos<-bscandat_new_sub3$pos+totpos
    bscandat<-rbind(bscandat,bscandat_new_sub3)
    last_scaff_l<-curr_scaff_l
    totpos <- totpos+last_scaff_l
  }
}
bscandat$lg<-factor(bscandat$lg)
rm(bscandat_new_sub2)
rm(bscandat_new_sub3)

head(bscandat)
str(bscandat)
nrow(bscandat)
summary(bscandat)
length(unique(bscandat$chr))


## Add an AbsPos to the "unplaced" data
head(bscandat_new_unpl)
bscandat_new_unpl$lg_snp_pos<-vector(length=nrow(bscandat_new_unpl))
totpos<-0
last_scaff_l<-0
bscandat_unpl<-data.frame()
for(scaff in unique(bscandat_new_unpl$chr)){
  bscandat_unpl_sub<-bscandat_new_unpl[bscandat_new_unpl$chr==scaff,]
  curr_scaff_l<-scaffdat$size[scaffdat$chr==scaff]
  # Get a new position for each SNP by:
  # lg_snp_pos <- startpos + last_scaff_l + scaff_snp_pos
  bscandat_unpl_sub$lg_snp_pos<-bscandat_unpl_sub$pos+totpos
  bscandat_unpl<-rbind(bscandat_unpl,bscandat_unpl_sub)
  last_scaff_l<-curr_scaff_l
  totpos <- totpos+last_scaff_l
}
head(bscandat_unpl)
length(unique(bscandat_unpl$chr))

# Write as table: UNPLACED SNPS
setwd(paste("~/Desktop/Data/RData/montana_clines/bscan_",pc,sep=""))
write.table(bscandat_unpl,
            paste("montana_clines_Aut_",pc,"_bscanin_fst_sub_unpl.txt",
                           sep=""),
            quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")

# Write as table: PLACED SNPS
write.table(bscandat,paste("montana_clines_Aut_",pc,"_bscanin_fst_sub_lg.txt",
                           sep=""),
            quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")


#----------------#
# PLOTTING       #
#----------------#
###
# PC1
#####
pc<-"pc1"
setwd(paste("~/Desktop/Data/RData/montana_clines/bscan_",pc,sep=""))
# Read table
bscandatpc1<-read.table(paste("montana_clines_Aut_",
                           pc,
                           "_bscanin_fst_sub_lg.txt",
                           sep=""),
                     header=TRUE,sep="\t")
bscandat_unpl_pc1<-read.table(paste("montana_clines_Aut_",
                                pc,
                                "_bscanin_fst_sub_unpl.txt",
                                sep=""),
                          header=TRUE,sep="\t")
head(bscandatpc1)
# Plotting: PLACED SNPS
manhplot<-ggplot()+
  geom_point(data=bscandatpc1,
             aes(lg_snp_pos/1000000,-log10(g_qval),colour=chr),
             size=0.2)+
  geom_hline(yintercept = -log10(0.05),colour="red")+
  scale_colour_manual("",
                      values=rep(c("black","grey"),
                                 length.out=length(unique(bscandatpc1$chr))),
                      guide=FALSE)+
  xlab("Position (Mb)")+
  ylab(expression(paste("-log10(q-value) of ",italic(g),sep="")))+
  facet_grid(lg~.)+
  ggtitle(paste("Results for ",pc,sep=""))+
  my.theme+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=10),
        strip.text.y = element_text(size=10,angle=0),
        plot.title = element_text(hjust=0,size=10))
manhplot

# How many "top" SNPs are there 
qval_thr<-0.05
nrow(bscandatpc1[bscandatpc1$g_qval < qval_thr,]) # Mapped

# How many SNPs on each lg
lg_snps<-tapply(bscandatpc1$scfsize,
                INDEX = list(bscandatpc1$lg),FUN=length)

# How many top SNPs on each lg
lg_top_snps<-tapply(bscandatpc1$scfsize[bscandatpc1$g_qval < qval_thr],
                    INDEX = list(bscandatpc1$lg[bscandatpc1$g_qval < qval_thr]),
                    FUN=length)
chisq.test(lg_snps)
chisq.test(lg_top_snps)
chisq.test(lg_top_snps,p = lg_snps/nrow(bscandatpc1))
str(chisq.test(lg_top_snps,p = lg_snps/nrow(bscandatpc1)))

# How long is each lg
lg_ls<-tapply(LGdat1.4$size,
              INDEX = list(LGdat1.4$linkage_gr),sum)[c(1,2,3,4)]/1000000
# Correlations?
cor.test(lg_snps,lg_ls,method="spearman")
cor.test(lg_top_snps,lg_ls,method="spearman")
plot(lg_ls,lg_snps)
plot(lg_ls,lg_top_snps)

nrow(bscandatpc1_unpl[bscandatpc1_unpl$g_qval < qval_thr,]) # UnMapped

head(bscandatpc1[bscandatpc1$g_qval < qval_thr,])

outdir<-paste("~/Desktop/Data/montana_project/clinal_population_genomics/results/bscanout/bscan_",pc,"/",sep="")
# Write the top SNPs as a table
#PLACED
write.table(bscandatpc1[bscandatpc1$g_qval < qval_thr,c(1,2,2)],
            paste(outdir,"montana_clines_Aut_",pc,"_bscanin_topSNPs.tab",
                  sep=""),sep="\t",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
#UNPLACED
write.table(bscandatpc1_unpl[bscandatpc1_unpl$g_qval < qval_thr,c(1,2,2)],
            paste(outdir,"montana_clines_Aut_",pc,"_bscanin_topSNPs_unpl.tab",
                  sep=""),sep="\t",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
#####


###
# PC2
######
pc<-"pc2"
setwd(paste("~/Desktop/Data/RData/montana_clines/bscan_",pc,sep=""))
# Read table
bscandatpc2<-read.table(paste("montana_clines_Aut_",
                              pc,
                              "_bscanin_fst_sub_lg.txt",
                              sep=""),
                        header=TRUE,sep="\t")
bscandat_unpl_pc2<-read.table(paste("montana_clines_Aut_",
                                    pc,
                                    "_bscanin_fst_sub_unpl.txt",
                                    sep=""),
                              header=TRUE,sep="\t")
head(bscandatpc2)
# Plotting: PLACED SNPS
manhplot<-ggplot()+
  geom_point(data=bscandatpc2,
             aes(lg_snp_pos/1000000,-log10(g_qval),colour=chr),
             size=0.2)+
  geom_hline(yintercept = -log10(0.05),colour="red")+
  scale_colour_manual("",
                      values=rep(c("black","grey"),
                                 length.out=length(unique(bscandatpc2$chr))),
                      guide=FALSE)+
  xlab("Position (Mb)")+
  ylab(expression(paste("-log10(q-value) of ",italic(g),sep="")))+
  facet_grid(lg~.)+
  ggtitle(paste("Results for ",pc,sep=""))+
  my.theme+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=10),
        strip.text.y = element_text(size=10,angle=0),
        plot.title = element_text(hjust=0,size=10))
manhplot

# How many "top" SNPs are there 
qval_thr<-0.05
nrow(bscandatpc2[bscandatpc2$g_qval < qval_thr,]) # Mapped

# How many SNPs on each lg
lg_snps<-tapply(bscandatpc2$scfsize,
                INDEX = list(bscandatpc2$lg),FUN=length)

# How many top SNPs on each lg
lg_top_snps<-tapply(bscandatpc2$scfsize[bscandatpc2$g_qval < qval_thr],
                    INDEX = list(bscandatpc2$lg[bscandatpc2$g_qval < qval_thr]),
                    FUN=length)
chisq.test(lg_snps)
chisq.test(lg_top_snps)
chisq.test(lg_top_snps,p = lg_snps/nrow(bscandatpc2))
str(chisq.test(lg_top_snps,p = lg_snps/nrow(bscandatpc1)))

# How long is each lg
lg_ls<-tapply(LGdat1.4$size,
              INDEX = list(LGdat1.4$linkage_gr),sum)[c(1,2,3,4)]/1000000
# Correlations?
cor.test(lg_snps,lg_ls,method="spearman")
cor.test(lg_top_snps,lg_ls,method="spearman")
plot(lg_ls,lg_snps)
plot(lg_ls,lg_top_snps)

nrow(bscandatpc2_unpl[bscandatpc2_unpl$g_qval < qval_thr,]) # UnMapped

head(bscandatpc2[bscandatpc2$g_qval < qval_thr,])

outdir<-paste("~/Desktop/Data/montana_project/clinal_population_genomics/results/bscanout/bscan_",pc,"/",sep="")
# Write the top SNPs as a table
#PLACED
write.table(bscandatpc2[bscandatpc2$g_qval < qval_thr,c(1,2,2)],
            paste(outdir,"montana_clines_Aut_",pc,"_bscanin_topSNPs.tab",
                  sep=""),sep="\t",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
#UNPLACED
write.table(bscandatpc2_unpl[bscandatpc2_unpl$g_qval < qval_thr,c(1,2,2)],
            paste(outdir,"montana_clines_Aut_",pc,"_bscanin_topSNPs_unpl.tab",
                  sep=""),sep="\t",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
#####


#----------------#
# CLOSEST GENES  #
#----------------#
# ####
outdir<-paste("~/Desktop/Data/montana_project/clinal_population_genomics/results/bscanout/",sep="")
# Load the closest genes table
clgenepc1<-read.table(paste(outdir,
                         "bscan_pc1/montana_clines_Aut_pc1_bscanin_topSNPs_closestgenes.tab",
                         sep=""),
                   header=FALSE,sep="\t")
colnames(clgenepc1)<-c("scaff","s_spos","s_epos",
                    "gscaff","g_spos","g_epos","name","sym","dist")
clgenepc2<-read.table(paste(outdir,
                            "bscan_pc2/montana_clines_Aut_pc2_bscanin_topSNPs_closestgenes.tab",
                            sep=""),
                      header=FALSE,sep="\t")
colnames(clgenepc2)<-c("scaff","s_spos","s_epos",
                       "gscaff","g_spos","g_epos","name","sym","dist")

head(clgenepc1)
head(clgenepc2)
# Remove those with no closest SNP
clgenepc1<-clgenepc1[clgenepc1$name != ".",]
clgenepc2<-clgenepc2[clgenepc2$name != ".",]

# How many
nrow(clgenepc1)
nrow(clgenepc2)

# How many within 1Mb
nrow(clgenepc1[abs(clgenepc1$dist)<1000000,])
clgenepc1_1Mb<-clgenepc1[abs(clgenepc1$dist)<1000000,]
nrow(clgenepc2[abs(clgenepc2$dist)<1000000,])
clgenepc2_1Mb<-clgenepc2[abs(clgenepc2$dist)<1000000,]

# How many within 1kb
nrow(clgenepc1[abs(clgenepc1$dist)<1000,])
clgenepc1_1kb<-clgenepc1[abs(clgenepc1$dist)<1000,]
nrow(clgenepc2[abs(clgenepc2$dist)<1000,])
clgenepc2_1kb<-clgenepc2[abs(clgenepc2$dist)<1000,]


# How many have FB ids
nrow(clgenepc1_1Mb[grep("FBgn",clgenepc1_1Mb$name),])
nrow(clgenepc2_1Mb[grep("FBgn",clgenepc2_1Mb$name),])

nrow(clgenepc1_1kb[grep("FBgn",clgenepc1_1kb$name),])
nrow(clgenepc2_1kb[grep("FBgn",clgenepc2_1kb$name),])

# Print FB IDs as lists for FlyBase
write.table(clgene_1Mb$name[grep("FBgn",clgene_1Mb$name)],
            paste(outdir,"montana_clines_Aut_",
                  pc,
                  "_bscanin_topSNPs_closestgenes_1Mb_FBid.list",
                  sep=""),sep="\t",
            quote=FALSE,col.names=FALSE,row.names=FALSE)

write.table(clgene_1kb$name[grep("FBgn",clgene_1kb$name)],
            paste(outdir,"montana_clines_Aut_",
                  pc,
                  "_bscanin_topSNPs_closestgenes_1kb_FBid.list",
                  sep=""),sep="\t",
            quote=FALSE,col.names=FALSE,row.names=FALSE)

# Load the converted table from FlyBase
pc1<-read.table("~/Desktop/Data/montana_project/clinal_population_genomics/results/genes/closest_genes_PC1_unique.tab",header=TRUE,sep="\t")
pc2<-read.table("~/Desktop/Data/montana_project/clinal_population_genomics/results/genes/closest_genes_PC2_unique.tab",header=TRUE,sep="\t")
both<-read.table("~/Desktop/Data/montana_project/clinal_population_genomics/results/genes/closest_genes_PC1PC2_common.tab",header=TRUE,sep="\t")
all<-rbind(pc1,pc2,both)


# Load Parker et al., 2015 Table S1
parker2015<-read.table("Parker_2015_Table_S1.csv",header=TRUE,sep=",")
# Remove any genes that are not Dvir orthologs
parker2015<-parker2015[grep("Dvir",parker2015$gene_name),]

head(parker2015)
# How many are DE in both
unique(parker2015$sp)
nrow(parker2015[parker2015$sp=="Both species",])
# How many genes occur in the Parker et al., 2015 set of genes that are 
# DE in both Dvir and Dmon or just Dmon
parker2015$Dmel_orth[parker2015$sp %in% c("Both species","D. montana") & 
                       parker2015$gene_name %in% pc1$SYMBOL]
pc1_o<-length(
  parker2015$Dmel_orth[parker2015$sp %in% c("Both species","D. montana") & 
                       parker2015$gene_name %in% pc1$SYMBOL])

parker2015$Dmel_orth[parker2015$sp %in% c("Both species","D. montana") & 
                       parker2015$gene_name %in% pc2$SYMBOL]
pc2_o<-length(
  parker2015$Dmel_orth[parker2015$sp %in% c("Both species","D. montana") & 
                       parker2015$gene_name %in% pc2$SYMBOL])

parker2015$Dmel_orth[parker2015$sp %in% c("Both species","D. montana") & 
                       parker2015$gene_name %in% both$SYMBOL]
pcboth_o<-length(
  parker2015$Dmel_orth[parker2015$sp %in% c("Both species","D. montana") & 
                       parker2015$gene_name %in% both$SYMBOL])

pc1_o
pc2_o
pcboth_o
n_overlaps<-pc1_o + pc2_o + pcboth_o
# How many expected by chance:
dmon_genes<-read.table("~/Desktop/Data/montana_project/genome/montana_genome/D_montana_genes.list",header=FALSE)
colnames(dmon_genes)<-c("gene")
head(dmon_genes)
# Bootstrap

n_genes<-nrow(pc1)+nrow(pc2)+nrow(both)
(n_overlaps/n_genes)*100

nboots<-1000
bootsamps<-vector(length=nboots)
for(boot in 1:nboots){
  genes<-sample(dmon_genes$gene,size = n_genes,replace = FALSE)
  bootsamps[boot]<-length(
    parker2015$Dmel_orth[parker2015$sp %in% c("Both species","D. montana") & 
                           parker2015$gene_name %in% genes])
}
# 'p-value'
(length(bootsamps[bootsamps>n_overlaps])/length(bootsamps))
hist(bootsamps)

# Load Kankare et al., 2016 diapause Table S1
kankare2016<-read.table("Kankare_2016_Table_S1.csv",header=TRUE,sep=",")
kankare2016$gene_name<-gsub("Dvir_","Dvir\\\\",kankare2016$gene_name)
nrow(kankare2016)
# Remove any genes that are not Dvir orthologs
kankare2016<-kankare2016[grep("Dvir",kankare2016$gene_name),]
nrow(kankare2016)

unique(kankare2016$gene_name[kankare2016$gene_name %in% pc1$SYMBOL][
  order(kankare2016$gene_name[kankare2016$gene_name %in% pc1$SYMBOL])])
pc1_o<-length(unique(kankare2016$gene_name[kankare2016$gene_name %in% pc1$SYMBOL]))

unique(kankare2016$gene_name[kankare2016$gene_name %in% pc2$SYMBOL][
  order(kankare2016$gene_name[kankare2016$gene_name %in% pc2$SYMBOL])])
pc2_o<-length(unique(kankare2016$gene_name[kankare2016$gene_name %in% pc2$SYMBOL]))

unique(kankare2016$gene_name[kankare2016$gene_name %in% both$SYMBOL][
  order(kankare2016$gene_name[kankare2016$gene_name %in% both$SYMBOL])])
pcboth_o<-length(unique(kankare2016$gene_name[kankare2016$gene_name %in% both$SYMBOL]))

pc1_o
pc2_o
pcboth_o
n_overlaps<-pc1_o + pc2_o + pcboth_o
n_genes<-nrow(pc1)+nrow(pc2)+nrow(both)
(n_overlaps/n_genes)*100

nboots<-1000
bootsamps<-vector(length=nboots)
for(boot in 1:nboots){
  genes<-sample(dmon_genes$gene,size = n_genes,replace = FALSE)
  bootsamps[boot]<-length(
    unique(kankare2016$gene_name[kankare2016$gene_name %in% genes]))
}
# 'p-value'
(length(bootsamps[bootsamps>n_overlaps])/length(bootsamps))
hist(bootsamps)


# Load Parker et al., 2016 changing light Table S1
parker2016<-read.table("Parker_2016_Table_S1.csv",header=TRUE,sep=",")
parker2016$gene_name<-gsub("Dvir_","Dvir\\\\",parker2016$gene_name)
nrow(parker2016)
# Remove any genes that are not Dvir orthologs
parker2016<-parker2016[grep("Dvir",parker2016$gene_name),]
head(parker2016)
nrow(parker2016)

unique(parker2016$gene_name[parker2016$gene_name %in% pc1$SYMBOL])
pc1_o<-length(unique(parker2016$gene_name[parker2016$gene_name %in% pc1$SYMBOL]))

unique(parker2016$gene_name[parker2016$gene_name %in% pc2$SYMBOL])
pc2_o<-length(unique(parker2016$gene_name[parker2016$gene_name %in% pc2$SYMBOL]))

unique(parker2016$gene_name[parker2016$gene_name %in% both$SYMBOL])
pcboth_o<-length(unique(parker2016$gene_name[parker2016$gene_name %in% both$SYMBOL]))

pc1_o
pc2_o
pcboth_o
n_overlaps<-pc1_o + pc2_o + pcboth_o
n_genes<-nrow(pc1)+nrow(pc2)+nrow(both)
(n_overlaps/n_genes)*100

nboots<-1000
bootsamps<-vector(length=nboots)
for(boot in 1:nboots){
  genes<-sample(dmon_genes$gene,size = n_genes,replace = FALSE)
  bootsamps[boot]<-length(
    unique(parker2016$gene_name[parker2016$gene_name %in% genes]))
}
# 'p-value'
(length(bootsamps[bootsamps>n_overlaps])/length(bootsamps))
hist(bootsamps)

# ####




# Plot the locations of important genes on the plot
# Genes which are DE in at least two of the previous RNAseq studies 
# and also occur in this study.
#####
# Common genes
commongenes<-c("FBgn0203225",
               "FBgn0204964","FBgn0205258",
               "FBgn0205992","FBgn0207876")
geneids<-c("CG42313","Pex12",
           "CG9008","Yp3","Slik")
commongenedat<-data.frame(lg=vector(length=length(commongenes)),
                    geneid=geneids,
                    gene=commongenes,
                    g_spos=vector(length=length(commongenes)),
                    g_epos=vector(length=length(commongenes)),
                    gscaff=vector(length=length(commongenes)),
                    g_lgspos=vector(length=length(commongenes)),
                    g_lgepos=vector(length=length(commongenes)))


# PC1
pc1_uniquegenes<-c("FBgn0205243","FBgn0207689")
pc1_geneids<-c("GJ1807","Inos")
pc1_genedat<-data.frame(lg=vector(length=length(pc1_uniquegenes)),
                        geneid=pc1_geneids,
                        gene=pc1_uniquegenes,
                        g_spos=vector(length=length(pc1_uniquegenes)),
                        g_epos=vector(length=length(pc1_uniquegenes)),
                        gscaff=vector(length=length(pc1_uniquegenes)),
                        g_lgspos=vector(length=length(pc1_uniquegenes)),
                        g_lgepos=vector(length=length(pc1_uniquegenes)))

# PC2
pc2_uniquegenes<-c("FBgn0199132","FBgn0203500",
                   "FBgn0204709","FBgn0210867")
pc2_geneids<-c("Ltn1","GJ16316","vri","Sap47")
pc2_genedat<-data.frame(lg=vector(length=length(pc2_uniquegenes)),
                        geneid=pc2_geneids,
                        gene=pc2_uniquegenes,
                        g_spos=vector(length=length(pc2_uniquegenes)),
                        g_epos=vector(length=length(pc2_uniquegenes)),
                        gscaff=vector(length=length(pc2_uniquegenes)),
                        g_lgspos=vector(length=length(pc2_uniquegenes)),
                        g_lgepos=vector(length=length(pc2_uniquegenes)))

# Combine with common genes
genedat<-rbind(commongenedat,pc1_genedat,pc2_genedat)
genedat$dat<-c(rep("Common",nrow(commongenedat)),
               rep("PC1",nrow(pc1_genedat)),
               rep("PC2",nrow(pc2_genedat)))


# Gather plotting data for arrows: PC1
i<-1
for(gene in genedat$gene){
  # if PC1 or Common:
  if(genedat$dat[genedat$gene==gene] =="PC1" | 
     genedat$dat[genedat$gene==gene] =="Common"){
    # Get the positions of SNPs near gene
    # If there are more than one SNP, take the closest.
    s_spos<-min(clgenepc1_1Mb$s_spos[clgenepc1_1Mb$name == gene])
    s_scaff<-clgenepc1_1Mb$scaff[clgenepc1_1Mb$name == gene]
    
    # Get the distance of the gene from the SNP
    # If there are more than one SNP, take the closest.
    dist<-clgenepc1_1Mb$dist[clgenepc1_1Mb$name == gene]
    
    # Get the length of the gene
    g_l<-clgenepc1_1Mb$g_epos[
      clgenepc1_1Mb$name == gene]-clgenepc1_1Mb$g_spos[
        clgenepc1_1Mb$name == gene]
    
    # Get the scaffold and scaffold pos of the gene
    genedat$g_spos[i]<-clgenepc1_1Mb$g_spos[clgenepc1_1Mb$name == gene]
    genedat$g_epos[i]<-clgenepc1_1Mb$g_epos[clgenepc1_1Mb$name == gene]
    genedat$gscaff[i]<-as.character(
      clgenepc1_1Mb$gscaff[clgenepc1_1Mb$name == gene])
    
    # Get the linkage group and linkage group Pos of the SNP
    lg<-bscandatpc1$lg[bscandatpc1$chr == as.character(s_scaff) & 
                      bscandat$pos == s_spos]
    lgabspos<-bscandatpc1$lg_snp_pos[bscandatpc1$chr == as.character(s_scaff) & 
                                       bscandatpc1$pos == s_spos]
    genedat$lg[i]<-as.character(lg)
    genedat$g_lgspos[i]<-lgabspos+dist
    genedat$g_lgepos[i]<-lgabspos+dist+g_l
    i<-i+1
  }else if(genedat$dat[genedat$gene==gene] == "PC2"){
    # Get the positions of SNPs near gene
    # If there are more than one SNP, take the closest.
    s_spos<-min(clgenepc2_1Mb$s_spos[clgenepc2_1Mb$name == gene])
    s_scaff<-clgenepc2_1Mb$scaff[clgenepc2_1Mb$name == gene]
    
    # Get the distance of the gene from the SNP
    # If there are more than one SNP, take the closest.
    dist<-clgenepc2_1Mb$dist[clgenepc2_1Mb$name == gene]
    
    # Get the length of the gene
    g_l<-clgenepc2_1Mb$g_epos[
      clgenepc2_1Mb$name == gene]-clgenepc2_1Mb$g_spos[
        clgenepc2_1Mb$name == gene]
    
    # Get the scaffold and scaffold pos of the gene
    genedat$g_spos[i]<-clgenepc2_1Mb$g_spos[clgenepc2_1Mb$name == gene]
    genedat$g_epos[i]<-clgenepc2_1Mb$g_epos[clgenepc2_1Mb$name == gene]
    genedat$gscaff[i]<-as.character(
      clgenepc2_1Mb$gscaff[clgenepc2_1Mb$name == gene])
    
    # Get the linkage group and linkage group Pos of the SNP
    lg<-bscandatpc2$lg[bscandatpc2$chr == as.character(s_scaff) & 
                      bscandatpc2$pos == s_spos]
    lgabspos<-bscandatpc2$lg_snp_pos[bscandatpc2$chr == as.character(s_scaff) & 
                                    bscandatpc2$pos == s_spos]
    genedat$lg[i]<-as.character(lg)
    genedat$g_lgspos[i]<-lgabspos+dist
    genedat$g_lgepos[i]<-lgabspos+dist+g_l
    i<-i+1
  }
}
genedat


# Plot manhattan plot with arrows and labels for interesting genes: PC1
manhplot + ggtitle("") +
  geom_segment(data=genedat[genedat$dat == "PC1" & genedat$dat == "Common",],
               aes(yend=3,y=3.5,
                   xend=(g_lgspos/1000000),
                   x=(g_lgspos/1000000)),
               arrow = arrow(length = unit(0.1,"cm")))+
  geom_text(data=genedat,aes(y=3.6,x=(g_lgepos/1000000),label=geneid),size=2)

#####





###
# Other interesting genes from this study
#####
commongenes<-c("FBgn0202755","FBgn0204659")
geneids<-c("Hsp60C","mio")
commongenedat<-data.frame(lg=vector(length=length(commongenes)),
                          geneid=geneids,
                          gene=commongenes,
                          g_spos=vector(length=length(commongenes)),
                          g_epos=vector(length=length(commongenes)),
                          gscaff=vector(length=length(commongenes)),
                          g_lgspos=vector(length=length(commongenes)),
                          g_lgepos=vector(length=length(commongenes)))

# PC1
pc1_uniquegenes<-c("FBgn0022834","FBgn0015211","FBgn0199202",
                   "FBgn0198702","FBgn0208901")
pc1_geneids<-c("tim","gl","DnaJ-1","bora","so")
pc1_genedat<-data.frame(lg=vector(length=length(pc1_uniquegenes)),
                        geneid=pc1_geneids,
                        gene=pc1_uniquegenes,
                        g_spos=vector(length=length(pc1_uniquegenes)),
                        g_epos=vector(length=length(pc1_uniquegenes)),
                        gscaff=vector(length=length(pc1_uniquegenes)),
                        g_lgspos=vector(length=length(pc1_uniquegenes)),
                        g_lgepos=vector(length=length(pc1_uniquegenes)))

# PC2
pc2_uniquegenes<-c("FBgn0198685","FBgn0201989","FBgn0208150","FBgn0210387")
pc2_geneids<-c("Clk","nocte","inaC","slmb")
pc2_genedat<-data.frame(lg=vector(length=length(pc2_uniquegenes)),
                        geneid=pc2_geneids,
                        gene=pc2_uniquegenes,
                        g_spos=vector(length=length(pc2_uniquegenes)),
                        g_epos=vector(length=length(pc2_uniquegenes)),
                        gscaff=vector(length=length(pc2_uniquegenes)),
                        g_lgspos=vector(length=length(pc2_uniquegenes)),
                        g_lgepos=vector(length=length(pc2_uniquegenes)))

# Combine with common genes
genedat<-rbind(commongenedat,pc1_genedat,pc2_genedat)
genedat$dat<-c(rep("Common",nrow(commongenedat)),
               rep("PC1",nrow(pc1_genedat)),
               rep("PC2",nrow(pc2_genedat)))


# Gather plotting data
i<-1
for(gene in genedat$gene){
  # if PC1 or Common:
  if(genedat$dat[genedat$gene==gene] =="PC1" | 
     genedat$dat[genedat$gene==gene] =="Common"){
    # Get the positions of SNPs near gene
    # If there are more than one SNP, take the closest.
    s_spos<-min(clgenepc1_1Mb$s_spos[clgenepc1_1Mb$name == gene])
    s_scaff<-clgenepc1_1Mb$scaff[clgenepc1_1Mb$name == gene]
    
    # Get the distance of the gene from the SNP
    # If there are more than one SNP, take the closest.
    dist<-clgenepc1_1Mb$dist[clgenepc1_1Mb$name == gene]
    
    # Get the length of the gene
    g_l<-clgenepc1_1Mb$g_epos[
      clgenepc1_1Mb$name == gene]-clgenepc1_1Mb$g_spos[
        clgenepc1_1Mb$name == gene]
    
    # Get the scaffold and scaffold pos of the gene
    genedat$g_spos[i]<-clgenepc1_1Mb$g_spos[clgenepc1_1Mb$name == gene]
    genedat$g_epos[i]<-clgenepc1_1Mb$g_epos[clgenepc1_1Mb$name == gene]
    genedat$gscaff[i]<-as.character(
      clgenepc1_1Mb$gscaff[clgenepc1_1Mb$name == gene])
    
    # Get the linkage group and linkage group Pos of the SNP
    lg<-bscandatpc1$lg[bscandatpc1$chr == as.character(s_scaff) & 
                         bscandat$pos == s_spos]
    lgabspos<-bscandatpc1$lg_snp_pos[bscandatpc1$chr == as.character(s_scaff) & 
                                       bscandatpc1$pos == s_spos]
    genedat$lg[i]<-as.character(lg)
    genedat$g_lgspos[i]<-lgabspos+dist
    genedat$g_lgepos[i]<-lgabspos+dist+g_l
    i<-i+1
  }else if(genedat$dat[genedat$gene==gene] == "PC2"){
    # Get the positions of SNPs near gene
    # If there are more than one SNP, take the closest.
    s_spos<-min(clgenepc2_1Mb$s_spos[clgenepc2_1Mb$name == gene])
    s_scaff<-clgenepc2_1Mb$scaff[clgenepc2_1Mb$name == gene]
    
    # Get the distance of the gene from the SNP
    # If there are more than one SNP, take the closest.
    dist<-clgenepc2_1Mb$dist[clgenepc2_1Mb$name == gene]
    
    # Get the length of the gene
    g_l<-clgenepc2_1Mb$g_epos[
      clgenepc2_1Mb$name == gene]-clgenepc2_1Mb$g_spos[
        clgenepc2_1Mb$name == gene]
    
    # Get the scaffold and scaffold pos of the gene
    genedat$g_spos[i]<-clgenepc2_1Mb$g_spos[clgenepc2_1Mb$name == gene]
    genedat$g_epos[i]<-clgenepc2_1Mb$g_epos[clgenepc2_1Mb$name == gene]
    genedat$gscaff[i]<-as.character(
      clgenepc2_1Mb$gscaff[clgenepc2_1Mb$name == gene])
    
    # Get the linkage group and linkage group Pos of the SNP
    lg<-bscandatpc2$lg[bscandatpc2$chr == as.character(s_scaff) & 
                         bscandatpc2$pos == s_spos]
    lgabspos<-bscandatpc2$lg_snp_pos[bscandatpc2$chr == as.character(s_scaff) & 
                                       bscandatpc2$pos == s_spos]
    genedat$lg[i]<-as.character(lg)
    genedat$g_lgspos[i]<-lgabspos+dist
    genedat$g_lgepos[i]<-lgabspos+dist+g_l
    i<-i+1
  }
}
genedat
#####



