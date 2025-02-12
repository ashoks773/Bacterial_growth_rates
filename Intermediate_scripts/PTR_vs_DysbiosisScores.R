#--- To Find out highly Correlated bacterial growth rates with Dysbiosis Scores for Each Patients
#source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/dysbiosis_Score_allset.R", echo=TRUE)
#- Computing dysbiosis scores is a one time Task Load directly for All samples

#--- Load AllSets Disease Activity
AllSets_Activity_Index <- read.csv(file="../Disease_Activity/AllSets_Activity_IndexUpdated.txt", sep = "\t", row.names = 1, header = T)

#Load meta
meta <- read.csv(file="~/Work/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)

AllSets_Activity_Indexmeta <- merge(AllSets_Activity_Index, meta, by=0, all=F)
rownames(AllSets_Activity_Indexmeta) <- AllSets_Activity_Indexmeta$Row.names; AllSets_Activity_Indexmeta$Row.names <- NULL

#---- Load RefGenomes PTRs
Ref_PTR_meta <- read.csv("~/Work/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/PTR_analysis/PTR_meta.txt", sep="\t", header = T, row.names =1)
Ref_PTR <- Ref_PTR_meta[,1:188]

#Merge new metadata
#Ref_PTR_meta <- merge(Ref_PTR, meta, by=0, all=F)
#rownames(Ref_PTR_meta) <- Ref_PTR_meta$Row.names; Ref_PTR_meta$Row.names <- NULL

#---- Load RefGenomes Relative abundances
Ref_Relab_meta <- read.csv("~/Work/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/Relab_analysis/Relabun_met.txt", sep="\t", header = T, row.names =1)
#-- Get same References Bacterial genomes (188) for which PTRs were Used
#bac_keeps <- c("X1297617.4.patric","X349741.6.patric","X484020.3.patric","X568816.4.patric","X657309.4.patric","X657319.15.patric","X709991.3.patric","X890402.3.patric","X1073375.3.patric","X1121101.3.patric","X1128111.23.patric","X1203465.3.patric","X1236512.3.patric","X1262754.3.patric","X1262776.3.patric","X1262777.3.patric","X1262889.3.patric","X1262947.3.patric","X1262949.3.patric","X1262979.3.patric","X1262989.3.patric","X1262992.3.patric","X1263070.3.patric","X1263072.3.patric","X1263073.3.patric","X1263074.4.patric","X1263105.3.patric","X1263107.3.patric","X1432052.9.patric","X1522.9.patric","X1650661.3.patric","X1673721.3.patric","X169435.7.patric","X1697787.3.patric","X1897011.3.patric","X1952006.3.patric","X208479.8.patric","X246787.5.patric","X2562617057.img","X2562617183.img","X261299.4.patric","X2654588179.img","X28116.19.patric","X33043.5.patric","X39488.4.patric","X445970.5.patric","X445974.6.patric","X457421.5.patric","X460384.4.patric","X46503.4.patric","X469591.4.patric","X469610.4.patric","X483215.6.patric","X545696.5.patric","X573.1549.patric","X59620.13.patric","X59620.24.patric","X59620.49.patric","X649756.4.patric","X665951.3.patric","X702446.3.patric","X742725.3.patric","X742821.3.patric","X765821.5.patric","X765821.8.patric","X817.56.patric","X927665.4.patric","X999410.3.patric","ERS235525_24.hgm","ERS235534_25.hgm","ERS235560_81.hgm","ERS235564_44.hgm","ERS396406_17.hgm","ERS396428_22.hgm","ERS396436_36.hgm","ERS396465_42.hgm","ERS396493_5.hgm","ERS396509_31.hgm","ERS396528_23.hgm","ERS473051_14.hgm","ERS473088_19.hgm","ERS473122_41.hgm","ERS473180_17.hgm","ERS473219_31.hgm","ERS473221_34.hgm","ERS473240_3.hgm","ERS473296_9.hgm","ERS473331_31.hgm","ERS473343_37.hgm","ERS473381_6.hgm","ERS473414_25.hgm","ERS537180_31.hgm","ERS537183_49.hgm","ERS537188_55.hgm","ERS537194_38.hgm","ERS537197_59.hgm","ERS537219_61.hgm","ERS537221_53.hgm","ERS537229_8.hgm","ERS537233_63.hgm","ERS537236_4.hgm","ERS537240_18.hgm","ERS537246_12.hgm","ERS537247_14.hgm","ERS537252_4.hgm","ERS537292_27.hgm","ERS537293_53.hgm","ERS537304_27.hgm","ERS537306_47.hgm","ERS537338_17.hgm","ERS537352_6.hgm","ERS537358_50.hgm","ERS537369_59.hgm","ERS537392_14.hgm","ERS608491_27.hgm","ERS608493_46.hgm","ERS608495_69.hgm","ERS608502_34.hgm","ERS608507_35.hgm","ERS608510_48.hgm","ERS608519_39.hgm","ERS608530_15.hgm","ERS608548_42.hgm","ERS608550_64.hgm","ERS608550_86.hgm","ERS608606_84.hgm","ERS631833_15.hgm","SRS015663_20.hgm","SRS017307_13.hgm","SRS019582_79.hgm","SRS019910_14.hgm","SRS020233_5.hgm","SRS020328_21.hgm","SRS021219_12.hgm","SRS024075_3.hgm","SRS049402_2.hgm","SRS050026_15.hgm","SRS050752_17.hgm","SRS056519_39.hgm","SRS075878_45.hgm","SRS076929_37.hgm","SRS077502_20.hgm","SRS077552_27.hgm","SRS077849_38.hgm","SRS098061_38.hgm","SRS098717_4.hgm","SRS104400_3.hgm","SRS142712_3.hgm","SRS143991_10.hgm","SRS144714_48.hgm","SRS147652_2.hgm","SRS1596777_76.hgm","SRS1596798_70.hgm","SRS1596798_81.hgm","SRS1596800_42.hgm","SRS1596807_67.hgm","SRS1596842_68.hgm","SRS1596849_88.hgm","SRS1596857_138.hgm","SRS1596872_171.hgm","SRS1719091_5.hgm","SRS1719444_21.hgm","SRS1719567_3.hgm","SRS1735456_15.hgm","SRS1735557_8.hgm","SRS1735679_4.hgm","SRS1735703_14.hgm","SRS259497_2.hgm","SRS294840_55.hgm","SRS294930_28.hgm","SRS294956_23.hgm","SRS294989_46.hgm","SRS295008_15.hgm","SRS295010_88.hgm","SRS295012_18.hgm","SRS311548_12.hgm","SRS332210_5.hgm","SRS372948_2.hgm","SRS373206_4.hgm","SRS475574_7.hgm","SRS475693_1.hgm","SRS476023_14.hgm","SRS476024_18.hgm","SRS476037_7.hgm","SRS476121_4.hgm","SRS476326_100.hgm","SRS820594_22.hgm","SRS820605_29.hgm")
bac_keeps <- colnames(Ref_PTR_meta[,1:188]) #-- Same genome which were used for PTR
Ref_Relab <- Ref_Relab_meta[bac_keeps]

#---Load Ref. Genome ID and Taxa Names
Refgid_gname <- read.csv("../PTR_analysis_on_IGG_Ref/RefGid_Gnames.txt", sep = "\t", header = T, row.names = 1)

###################################################################################
# Correlation between Dysbiosis Scores and RefGenomes PTRs and Relative Abundances
##################################################################################
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3001_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3003_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3006_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3002_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3008_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3009_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3004_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3010_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3005_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3012_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3013_ds.R", echo=TRUE) #C3013
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3015_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3016_ds.R", echo=TRUE)

source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3017_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3011_ds.R", echo=TRUE) #C3011
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3022_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3023_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3021_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3027_ds.R", echo=TRUE) #C3027
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3028_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3029_ds.R", echo=TRUE)

source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3030_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3031_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3032_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3034_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3035_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/C3037_ds.R", echo=TRUE)

source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/E5001_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/E5004_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/E5009_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/E5013_ds.R", echo=TRUE)

source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4001_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4004_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4006_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4007_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4008_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4009_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4010_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4013_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4014_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4015_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4016_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4017_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4018_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4019_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4020_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4022_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4023_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4024_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4027_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4028_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4030_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4031_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4032_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4035_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4038_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4039_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4040_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4042_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4043_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4044_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/H4045_ds.R", echo=TRUE)

source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2008_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2014_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2021_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2025_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2026_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2027_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2028_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2034_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2039_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2041_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2042_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2047_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2048_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2060_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2061_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2064_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2068_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2069_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2071_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2072_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2075_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2077_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2079_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2083_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2084_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2085_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2097_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/M2103_ds.R", echo=TRUE)

source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6005_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6009_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6010_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6012_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6013_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6014_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6016_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6017_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6018_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6024_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6025_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6028_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6033_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6035_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6037_ds.R", echo=TRUE)
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts/pid_scripts/P6038_ds.R", echo=TRUE)

#--- Try Loop that is not working
#for (i in names(unique_pid)){
for (i in unique_pid){
  C3001_diseaseAI <- subset(AllSets_Activity_Indexmeta, Participant.ID==i)
  #range_d <- range (list(C3001_diseaseAI[,1:5])) #-- Only for Bins PTR and Bins Relabundances
  
  C3001_ref_PTR <- subset(Ref_PTR_meta, pid==i)
  C3001_ref_PTR_filtered <- dropspc(C3001_ref_PTR[,1:188], 4) #-- Remove RefGenomes which is not present in Even 50% of the samples
  
  C3001_ref_relab <- subset(Ref_Relab_meta, pid==i)
  C3001_ref_relab_filtered <- dropspc(C3001_ref_relab[,1:188], 4) 
  
  C3001_ref_PTR_cor_dysbiosis <- cor(C3001_diseaseAI$Ref_PTR, as.matrix(C3001_ref_PTR_filtered[,1:ncol(C3001_ref_PTR_filtered)]),  method="pearson")
  C3001_ref_PTR_cor_dysbiosis <- data.frame(t(C3001_ref_PTR_cor_dysbiosis))
  C3001_ref_relab_cor_dysbiosis <- cor(C3001_diseaseAI$Ref_Relab, as.matrix(C3001_ref_relab_filtered[,1:ncol(C3001_ref_relab_filtered)]),  method="pearson")
  C3001_ref_relab_cor_dysbiosis <- data.frame(t(C3001_ref_relab_cor_dysbiosis))
  
  Ref_Genome_PTR_relab_Corrs <- merge(C3001_ref_PTR_cor_dysbiosis, C3001_ref_relab_cor_dysbiosis, by=0, all=T)
  rownames(Ref_Genome_PTR_relab_Corrs) <- Ref_Genome_PTR_relab_Corrs$Row.names; Ref_Genome_PTR_relab_Corrs$Row.names <- NULL
  colnames(Ref_Genome_PTR_relab_Corrs) <- c("Dysbiosis score (RefGenome PTRs)", "Dysbiosis score (RefGenome Relabun)")
  
  #--- Fix Id names in Corr File
  rownames(Ref_Genome_PTR_relab_Corrs) <- gsub(".hgm", "", rownames(Ref_Genome_PTR_relab_Corrs), fixed=TRUE)
  rownames(Ref_Genome_PTR_relab_Corrs) <- gsub(".patric", "", rownames(Ref_Genome_PTR_relab_Corrs), fixed=TRUE)
  rownames(Ref_Genome_PTR_relab_Corrs) <- gsub(".img", "", rownames(Ref_Genome_PTR_relab_Corrs), fixed=TRUE)
  rownames(Ref_Genome_PTR_relab_Corrs) <- gsub("X", "", rownames(Ref_Genome_PTR_relab_Corrs), fixed=TRUE)
  #--- Give IDs a Name
  Ref_Genome_PTR_relab_Corrs_TaxaNames <- merge(Ref_Genome_PTR_relab_Corrs, Refgid_gname, by=0, all=F)
  rownames(Ref_Genome_PTR_relab_Corrs_TaxaNames) <- Ref_Genome_PTR_relab_Corrs_TaxaNames$Row.names; Ref_Genome_PTR_relab_Corrs_TaxaNames$Row.names <- NULL
  
  #write.table(Ref_Genome_PTR_relab_Corrs_TaxaNames, "PatientSpceific_Ref_taxa_dysbiosis_scores/C3001_Ref_Genome_PTR_relab_Corrs.txt", sep = "\t")
  #filename <- paste(i, ".txt", sep="")
  #write.table(filename, col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
  write.table(unique_pid[i], file = paste(names(unique_pid[i]), ".txt", sep = ""), col.names= TRUE, sep = "\t", quote=FALSE)
  #write.csv2(get(unique_pid[i]),paste0("PatientSpceific_Ref_taxa_dysbiosis_scores/", unique_pid[i],".csv"),row.names = FALSE)
}
