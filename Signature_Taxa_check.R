#-- This script is to check the Distribution PTR and Relative abundances of Signature taxa (Reported in Nature 2019 paper) in all samples 

Colors=c("darkred", "steelblue2", "orange")
library(ggplot2)
library(egg)

###########################################
#Load PTRs Calculated on Reference Genomes
###########################################
Ref_PTR_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/PTR_analysis/PTR_meta.txt", sep="\t", header = T, row.names =1)
Ref_PTR <- Ref_PTR_meta[,1:188]

#--meta
meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)
#--- If Needed to fix the meta
source('~/Box/David_Casero_Lab/HMP_IBD_Project/Common_scripts/fix_metaata.R')
meta_new <- fix_metadata(meta) #-- Check If required otherwise- Its fine to use the Original Metadata Only

#Merge new metadata
Ref_PTR_meta <- merge(Ref_PTR, meta, by=0, all=F)
rownames(Ref_PTR_meta) <- Ref_PTR_meta$Row.names; Ref_PTR_meta$Row.names <- NULL

########################################################
#Load Relative Abundances Calculated on Reference Genomes
#########################################################
Ref_Relab_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/Relab_analysis/Relabun_met.txt", sep="\t", header = T, row.names =1)
#-- Get same References Bacterial genomes (188) for which PTRs were Used
bac_keeps <- c("X1297617.4.patric","X349741.6.patric","X484020.3.patric","X568816.4.patric","X657309.4.patric","X657319.15.patric","X709991.3.patric","X890402.3.patric","X1073375.3.patric","X1121101.3.patric","X1128111.23.patric","X1203465.3.patric","X1236512.3.patric","X1262754.3.patric","X1262776.3.patric","X1262777.3.patric","X1262889.3.patric","X1262947.3.patric","X1262949.3.patric","X1262979.3.patric","X1262989.3.patric","X1262992.3.patric","X1263070.3.patric","X1263072.3.patric","X1263073.3.patric","X1263074.4.patric","X1263105.3.patric","X1263107.3.patric","X1432052.9.patric","X1522.9.patric","X1650661.3.patric","X1673721.3.patric","X169435.7.patric","X1697787.3.patric","X1897011.3.patric","X1952006.3.patric","X208479.8.patric","X246787.5.patric","X2562617057.img","X2562617183.img","X261299.4.patric","X2654588179.img","X28116.19.patric","X33043.5.patric","X39488.4.patric","X445970.5.patric","X445974.6.patric","X457421.5.patric","X460384.4.patric","X46503.4.patric","X469591.4.patric","X469610.4.patric","X483215.6.patric","X545696.5.patric","X573.1549.patric","X59620.13.patric","X59620.24.patric","X59620.49.patric","X649756.4.patric","X665951.3.patric","X702446.3.patric","X742725.3.patric","X742821.3.patric","X765821.5.patric","X765821.8.patric","X817.56.patric","X927665.4.patric","X999410.3.patric","ERS235525_24.hgm","ERS235534_25.hgm","ERS235560_81.hgm","ERS235564_44.hgm","ERS396406_17.hgm","ERS396428_22.hgm","ERS396436_36.hgm","ERS396465_42.hgm","ERS396493_5.hgm","ERS396509_31.hgm","ERS396528_23.hgm","ERS473051_14.hgm","ERS473088_19.hgm","ERS473122_41.hgm","ERS473180_17.hgm","ERS473219_31.hgm","ERS473221_34.hgm","ERS473240_3.hgm","ERS473296_9.hgm","ERS473331_31.hgm","ERS473343_37.hgm","ERS473381_6.hgm","ERS473414_25.hgm","ERS537180_31.hgm","ERS537183_49.hgm","ERS537188_55.hgm","ERS537194_38.hgm","ERS537197_59.hgm","ERS537219_61.hgm","ERS537221_53.hgm","ERS537229_8.hgm","ERS537233_63.hgm","ERS537236_4.hgm","ERS537240_18.hgm","ERS537246_12.hgm","ERS537247_14.hgm","ERS537252_4.hgm","ERS537292_27.hgm","ERS537293_53.hgm","ERS537304_27.hgm","ERS537306_47.hgm","ERS537338_17.hgm","ERS537352_6.hgm","ERS537358_50.hgm","ERS537369_59.hgm","ERS537392_14.hgm","ERS608491_27.hgm","ERS608493_46.hgm","ERS608495_69.hgm","ERS608502_34.hgm","ERS608507_35.hgm","ERS608510_48.hgm","ERS608519_39.hgm","ERS608530_15.hgm","ERS608548_42.hgm","ERS608550_64.hgm","ERS608550_86.hgm","ERS608606_84.hgm","ERS631833_15.hgm","SRS015663_20.hgm","SRS017307_13.hgm","SRS019582_79.hgm","SRS019910_14.hgm","SRS020233_5.hgm","SRS020328_21.hgm","SRS021219_12.hgm","SRS024075_3.hgm","SRS049402_2.hgm","SRS050026_15.hgm","SRS050752_17.hgm","SRS056519_39.hgm","SRS075878_45.hgm","SRS076929_37.hgm","SRS077502_20.hgm","SRS077552_27.hgm","SRS077849_38.hgm","SRS098061_38.hgm","SRS098717_4.hgm","SRS104400_3.hgm","SRS142712_3.hgm","SRS143991_10.hgm","SRS144714_48.hgm","SRS147652_2.hgm","SRS1596777_76.hgm","SRS1596798_70.hgm","SRS1596798_81.hgm","SRS1596800_42.hgm","SRS1596807_67.hgm","SRS1596842_68.hgm","SRS1596849_88.hgm","SRS1596857_138.hgm","SRS1596872_171.hgm","SRS1719091_5.hgm","SRS1719444_21.hgm","SRS1719567_3.hgm","SRS1735456_15.hgm","SRS1735557_8.hgm","SRS1735679_4.hgm","SRS1735703_14.hgm","SRS259497_2.hgm","SRS294840_55.hgm","SRS294930_28.hgm","SRS294956_23.hgm","SRS294989_46.hgm","SRS295008_15.hgm","SRS295010_88.hgm","SRS295012_18.hgm","SRS311548_12.hgm","SRS332210_5.hgm","SRS372948_2.hgm","SRS373206_4.hgm","SRS475574_7.hgm","SRS475693_1.hgm","SRS476023_14.hgm","SRS476024_18.hgm","SRS476037_7.hgm","SRS476121_4.hgm","SRS476326_100.hgm","SRS820594_22.hgm","SRS820605_29.hgm")
Ref_Relab <- Ref_Relab_meta[bac_keeps]

#--meta
meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)
#--- If Needed to fix the meta
source('~/Box/David_Casero_Lab/HMP_IBD_Project/Common_scripts/fix_metaata.R')
meta_new <- fix_metadata(meta) #-- Check If required otherwise- Its fine to use the Original Metadata Only

#Merge new metadata
Ref_Relab_meta <- merge(Ref_Relab, meta, by=0, all=F)
rownames(Ref_Relab_meta) <- Ref_Relab_meta$Row.names; Ref_Relab_meta$Row.names <- NULL


#---- Now make density plot of Signature Taxa
#ERS537229_8.hgm - Faecalibacterium prausnitzii_A
# ERS396428_22.hgm - Roseburia hominis
# X665951.3 - Ruminococcus torques
#X1073375.3.patric - Ruminococcus gnavus
# X2562617057 - Escherichia coli

#-- Faecali
#Fecali <- c("ERS537229_8.hgm","diagnosis")
#Fecali_PTR <- Ref_PTR_meta[Fecali]
#Fecali_PTR <- subset(Fecali_PTR, ERS537229_8.hgm != 0) #-Remove samples in which Growth rates were not estimated
jpeg("SignatureTaxa/FaecalibacteriumPTR.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_PTR_meta, aes(x=sqrt(ERS537229_8.hgm), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Faecalibacterium prausnitzii PTR") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

#Fecali_Relab <- Ref_Relab_meta[Fecali]
#Fecali_Relab <- subset(Fecali_Relab, ERS537229_8.hgm != 0) #-Remove samples in which Abudance were not estimated
jpeg("SignatureTaxa/Faecalibacterium_Relab.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_Relab_meta, aes(x=sqrt(ERS537229_8.hgm), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Faecalibacterium prausnitzii Relab") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

#-- Roseburia
#rosubria <- c("ERS396428_22.hgm","diagnosis")
#rosubria_PTR <- Ref_PTR_meta[rosubria]
#rosubria_PTR <- subset(rosubria_PTR, ERS396428_22.hgm != 0) #-Remove samples in which Growth rates were not estimated
jpeg("SignatureTaxa/RoseburiaPTR.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_PTR_meta, aes(x=sqrt(ERS396428_22.hgm), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Roseburia hominis PTR") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()
jpeg("SignatureTaxa/Roseburia_Relab.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_Relab_meta, aes(x=sqrt(ERS396428_22.hgm), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Roseburia hominis relab") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

#-- Ruminococcus torques
jpeg("SignatureTaxa/R.torquesPTR.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_PTR_meta, aes(x=sqrt(X665951.3.patric), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Ruminococcus torques PTR") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()
jpeg("SignatureTaxa/R.torques_Relab.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_Relab_meta, aes(x=sqrt(X665951.3.patric), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Ruminococcus torques relab") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

#-- Ruminococcus gnavus
jpeg("SignatureTaxa/R.gnavusPTR.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_PTR_meta, aes(x=sqrt(X1073375.3.patric), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Ruminococcus gnavus PTR") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()
jpeg("SignatureTaxa/R.gnavus_Relab.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_Relab_meta, aes(x=sqrt(X1073375.3.patric), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Ruminococcus gnavus relab") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

#-- Escherichia coli
jpeg("SignatureTaxa/E.coliPTR.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_PTR_meta, aes(x=sqrt(X2562617057.img), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Escherichia coli PTR") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()
jpeg("SignatureTaxa/E.coli_Relab.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=Ref_Relab_meta, aes(x=sqrt(X2562617057.img), color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  #geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Escherichia coli relab") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()
