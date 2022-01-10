PTR_meta_no_filtered <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta_no_filtered.txt", sep="\t", header = T, row.names =1)
C3001 <- subset(PTR_meta_no_filtered, pid=="C3001")
#-- Remove Bin which is not present in Even 50% of the samples
C3001_filtered <- dropspc(C3001[,1:1503], 4)
C3001_bins <- colnames(C3001_filtered) #These are the Bins which are representative of Sample C3001
#These are the 70 Bins which is consistently present in All samples
bins_keeps <- c("Bin104","Bin1049","Bin12","Bin1258","Bin1265","Bin1291","Bin1393","Bin1479","Bin151","Bin1549","Bin1555","Bin1560","Bin1633","Bin1677","Bin168","Bin1712","Bin175","Bin1811","Bin1818","Bin1898","Bin1908","Bin1962","Bin1966","Bin197","Bin1970","Bin2018","Bin2019","Bin2023","Bin2057","Bin2075","Bin2106","Bin2132","Bin2212","Bin2332","Bin2333","Bin2338","Bin2365","Bin2434","Bin2480","Bin2483","Bin2582","Bin282","Bin2837","Bin3052","Bin3091","Bin3140","Bin3317","Bin3329","Bin3340","Bin3347","Bin3362","Bin3587","Bin3685","Bin3733","Bin3746","Bin3838","Bin3842","Bin387","Bin390","Bin3955","Bin4039","Bin4103","Bin418","Bin421","Bin584","Bin620","Bin627","Bin645","Bin687","Bin780")

my_list <- list(C3001_bins, bins_keeps)
C3001_uniq_list <- unique(unlist(my_list)) #-- Total Unique Bins including the Bins specific to Sample C3001
#--- Final set
binsPTR <- PTR_meta_no_filtered[C3001_uniq_list]
#--meta
meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)

#Merge new metadata
binsPTR_meta <- merge(binsPTR, meta, by=0, all=F)
rownames(binsPTR_meta) <- binsPTR_meta$Row.names; binsPTR_meta$Row.names <- NULL

#--- No need to Calculate Relative abundances on binsPTR data
library(vegan)
binsPTR_D <- as.matrix(vegdist(binsPTR_meta[,1:ncol(binsPTR)], method="bray"))

# The "reference set" of inactive samples
binsPTR_ref_set <- (binsPTR_meta$diagnosis == "nonIBD") & (binsPTR_meta$week_num >= 20)

# Calculate the activity index
binsPTR_meta$activity_index <- sapply(seq_along(binsPTR_ref_set), function(i)
  median(binsPTR_D[i, binsPTR_ref_set & (binsPTR_meta$site_sub_coll != binsPTR_meta$site_sub_coll[i])]))

# Threshold activity
binsPTR_disease_activity_threshold <- quantile(binsPTR_meta$activity_index[binsPTR_meta$diagnosis=="nonIBD"], 0.9)

binsPTR_activity_index <- data.frame (binsPTR_meta$activity_index)
rownames(binsPTR_activity_index) <- row.names(binsPTR_meta)
library(tidyr)
binsPTR_activity_index <- binsPTR_activity_index %>% drop_na()

############################################
#----- Calculate on Bins Relative Abundances
BinsRelab_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/RelativeAbun_analysis/Relabun_meta_no_filtered.txt", sep="\t", header = T, row.names =1)
#--C3001_uniq_list contains 70 Bins + Represetative bins of C3001
BinsRelab <- BinsRelab_meta[C3001_uniq_list]

#Merge new metadata
BinsRelab_meta <- merge(BinsRelab, meta, by=0, all=F)
rownames(BinsRelab_meta) <- BinsRelab_meta$Row.names; BinsRelab_meta$Row.names <- NULL

#--- No need to Calculate Relative abundances on PTR data
library(vegan)
binsRelab_D <- as.matrix(vegdist(BinsRelab_meta[,1:ncol(BinsRelab)], method="bray"))

# The "reference set" of inactive samples
binsRelab_ref_set <- (BinsRelab_meta$diagnosis == "nonIBD") & (BinsRelab_meta$week_num >= 20)

# Calculate the activity index
BinsRelab_meta$activity_index <- sapply(seq_along(binsRelab_ref_set), function(i)
  median(binsRelab_D[i, binsRelab_ref_set & (BinsRelab_meta$site_sub_coll != BinsRelab_meta$site_sub_coll[i])]))

# Threshold activity
BinsRelab_disease_activity_threshold <- quantile(BinsRelab_meta$activity_index[BinsRelab_meta$diagnosis=="nonIBD"], 0.9)

binsRelab_activity_index <- data.frame (BinsRelab_meta$activity_index)
rownames(binsRelab_activity_index) <- row.names(BinsRelab_meta)
library(tidyr)
binsRelab_activity_index <- binsRelab_activity_index %>% drop_na()

############################################
#----- Calculate on Reference Genome PTRs
Ref_PTR_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/PTR_analysis/PTR_meta.txt", sep="\t", header = T, row.names =1)
Ref_PTR <- Ref_PTR_meta[,1:188]

#Merge new metadata
Ref_PTR_meta <- merge(Ref_PTR, meta, by=0, all=F)
rownames(Ref_PTR_meta) <- Ref_PTR_meta$Row.names; Ref_PTR_meta$Row.names <- NULL

#--- No need to Calculate Relative abundances on PTR data
library(vegan)
Ref_PTR_D <- as.matrix(vegdist(Ref_PTR_meta[,1:ncol(Ref_PTR)], method="bray"))

# The "reference set" of inactive samples
Ref_PTR_ref_set <- (Ref_PTR_meta$diagnosis == "nonIBD") & (Ref_PTR_meta$week_num >= 20)

# Calculate the activity index
Ref_PTR_meta$activity_index <- sapply(seq_along(Ref_PTR_ref_set), function(i)
  median(Ref_PTR_D[i, Ref_PTR_ref_set & (Ref_PTR_meta$site_sub_coll != Ref_PTR_meta$site_sub_coll[i])]))

Ref_PTR_activity_index <- data.frame (Ref_PTR_meta$activity_index)
rownames(Ref_PTR_activity_index) <- row.names(Ref_PTR_meta)
library(tidyr)
Ref_PTR_activity_index <- Ref_PTR_activity_index %>% drop_na()

#########################################################
#----- Calculate on Reference Genome Relative Abundances
Ref_Relab_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/Relab_analysis/Relabun_met.txt", sep="\t", header = T, row.names =1)
#-- Get same References Bacterial genomes (188) for which PTRs were Used
#bac_keeps <- c("X1297617.4.patric","X349741.6.patric","X484020.3.patric","X568816.4.patric","X657309.4.patric","X657319.15.patric","X709991.3.patric","X890402.3.patric","X1073375.3.patric","X1121101.3.patric","X1128111.23.patric","X1203465.3.patric","X1236512.3.patric","X1262754.3.patric","X1262776.3.patric","X1262777.3.patric","X1262889.3.patric","X1262947.3.patric","X1262949.3.patric","X1262979.3.patric","X1262989.3.patric","X1262992.3.patric","X1263070.3.patric","X1263072.3.patric","X1263073.3.patric","X1263074.4.patric","X1263105.3.patric","X1263107.3.patric","X1432052.9.patric","X1522.9.patric","X1650661.3.patric","X1673721.3.patric","X169435.7.patric","X1697787.3.patric","X1897011.3.patric","X1952006.3.patric","X208479.8.patric","X246787.5.patric","X2562617057.img","X2562617183.img","X261299.4.patric","X2654588179.img","X28116.19.patric","X33043.5.patric","X39488.4.patric","X445970.5.patric","X445974.6.patric","X457421.5.patric","X460384.4.patric","X46503.4.patric","X469591.4.patric","X469610.4.patric","X483215.6.patric","X545696.5.patric","X573.1549.patric","X59620.13.patric","X59620.24.patric","X59620.49.patric","X649756.4.patric","X665951.3.patric","X702446.3.patric","X742725.3.patric","X742821.3.patric","X765821.5.patric","X765821.8.patric","X817.56.patric","X927665.4.patric","X999410.3.patric","ERS235525_24.hgm","ERS235534_25.hgm","ERS235560_81.hgm","ERS235564_44.hgm","ERS396406_17.hgm","ERS396428_22.hgm","ERS396436_36.hgm","ERS396465_42.hgm","ERS396493_5.hgm","ERS396509_31.hgm","ERS396528_23.hgm","ERS473051_14.hgm","ERS473088_19.hgm","ERS473122_41.hgm","ERS473180_17.hgm","ERS473219_31.hgm","ERS473221_34.hgm","ERS473240_3.hgm","ERS473296_9.hgm","ERS473331_31.hgm","ERS473343_37.hgm","ERS473381_6.hgm","ERS473414_25.hgm","ERS537180_31.hgm","ERS537183_49.hgm","ERS537188_55.hgm","ERS537194_38.hgm","ERS537197_59.hgm","ERS537219_61.hgm","ERS537221_53.hgm","ERS537229_8.hgm","ERS537233_63.hgm","ERS537236_4.hgm","ERS537240_18.hgm","ERS537246_12.hgm","ERS537247_14.hgm","ERS537252_4.hgm","ERS537292_27.hgm","ERS537293_53.hgm","ERS537304_27.hgm","ERS537306_47.hgm","ERS537338_17.hgm","ERS537352_6.hgm","ERS537358_50.hgm","ERS537369_59.hgm","ERS537392_14.hgm","ERS608491_27.hgm","ERS608493_46.hgm","ERS608495_69.hgm","ERS608502_34.hgm","ERS608507_35.hgm","ERS608510_48.hgm","ERS608519_39.hgm","ERS608530_15.hgm","ERS608548_42.hgm","ERS608550_64.hgm","ERS608550_86.hgm","ERS608606_84.hgm","ERS631833_15.hgm","SRS015663_20.hgm","SRS017307_13.hgm","SRS019582_79.hgm","SRS019910_14.hgm","SRS020233_5.hgm","SRS020328_21.hgm","SRS021219_12.hgm","SRS024075_3.hgm","SRS049402_2.hgm","SRS050026_15.hgm","SRS050752_17.hgm","SRS056519_39.hgm","SRS075878_45.hgm","SRS076929_37.hgm","SRS077502_20.hgm","SRS077552_27.hgm","SRS077849_38.hgm","SRS098061_38.hgm","SRS098717_4.hgm","SRS104400_3.hgm","SRS142712_3.hgm","SRS143991_10.hgm","SRS144714_48.hgm","SRS147652_2.hgm","SRS1596777_76.hgm","SRS1596798_70.hgm","SRS1596798_81.hgm","SRS1596800_42.hgm","SRS1596807_67.hgm","SRS1596842_68.hgm","SRS1596849_88.hgm","SRS1596857_138.hgm","SRS1596872_171.hgm","SRS1719091_5.hgm","SRS1719444_21.hgm","SRS1719567_3.hgm","SRS1735456_15.hgm","SRS1735557_8.hgm","SRS1735679_4.hgm","SRS1735703_14.hgm","SRS259497_2.hgm","SRS294840_55.hgm","SRS294930_28.hgm","SRS294956_23.hgm","SRS294989_46.hgm","SRS295008_15.hgm","SRS295010_88.hgm","SRS295012_18.hgm","SRS311548_12.hgm","SRS332210_5.hgm","SRS372948_2.hgm","SRS373206_4.hgm","SRS475574_7.hgm","SRS475693_1.hgm","SRS476023_14.hgm","SRS476024_18.hgm","SRS476037_7.hgm","SRS476121_4.hgm","SRS476326_100.hgm","SRS820594_22.hgm","SRS820605_29.hgm")
bac_keeps <- colnames(Ref_PTR) #-- Same genome which were used for PTR
Ref_Relab <- Ref_Relab_meta[bac_keeps]

#Merge new metadata
Ref_Relab_meta <- merge(Ref_Relab, meta, by=0, all=F)
rownames(Ref_Relab_meta) <- Ref_Relab_meta$Row.names; Ref_Relab_meta$Row.names <- NULL

#--- No need to Calculate Relative abundances on PTR data
library(vegan)
Ref_Relab_D <- as.matrix(vegdist(Ref_Relab_meta[,1:ncol(Ref_Relab)], method="bray"))

# The "reference set" of inactive samples
Ref_Relab_ref_set <- (Ref_Relab_meta$diagnosis == "nonIBD") & (Ref_Relab_meta$week_num >= 20)

# Calculate the activity index
Ref_Relab_meta$activity_index <- sapply(seq_along(Ref_Relab_ref_set), function(i)
  median(Ref_Relab_D[i, Ref_Relab_ref_set & (Ref_Relab_meta$site_sub_coll != Ref_Relab_meta$site_sub_coll[i])]))

Ref_Relab_activity_index <- data.frame (Ref_Relab_meta$activity_index)
rownames(Ref_Relab_activity_index) <- row.names(Ref_Relab_meta)
library(tidyr)
Ref_Relab_activity_index <- Ref_Relab_activity_index %>% drop_na()

#########################################################
#----- Calculate on MetaPhln Relative Abundances
#bugs <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/hmp2_analysis/taxonomic_profiles.tsv", sep="\t", row.names = 1, header = T)
#bugs <- data.frame(t(bugs))
#colnames(bugs)
#--- Get Bacteria and Till species level annotation
#bugs <- bugs[, grepl("k__Bacteria", colnames(bugs))]
#bugs <- bugs[, grepl("s__", colnames(bugs))]
#-- All Bugs should be Present in Atleast 5% of the samples (82 out of 1638 samples)
#bugs_filtered <- dropspc(bugs, 82)
#write.table(bugs_filtered, file= "~/Box/David_Casero_Lab/HMP_IBD_Project/hmp2_analysis/bugs_filtered.txt", sep="\t")
bugs_filtered <- read.csv(file= "~/Box/David_Casero_Lab/HMP_IBD_Project/hmp2_analysis/bugs_filtered.txt", sep="\t", row.names = 1, header = T)

bugs_meta <- merge(bugs_filtered, meta, by=0, all=F)
rownames(bugs_meta) <- bugs_meta$Row.names; bugs_meta$Row.names <- NULL

library(vegan)
bugs_relab <- decostand(bugs_meta[,1:ncol(bugs_filtered)], method = "total")*100
metaphln_D <- as.matrix(vegdist(bugs_relab, method="bray"))

# The "reference set" of inactive samples
metaphln_ref_set <- (bugs_meta$diagnosis == "nonIBD") & (bugs_meta$week_num >= 20)

# Calculate the activity index
bugs_meta$activity_index <- sapply(seq_along(metaphln_ref_set), function(i)
  median(metaphln_D[i, metaphln_ref_set & (bugs_meta$site_sub_coll != bugs_meta$site_sub_coll[i])]))

metaphln_activity_index <- data.frame (bugs_meta$activity_index)
rownames(metaphln_activity_index) <- row.names(bugs_meta)
library(tidyr)
metaphln_activity_index <- metaphln_activity_index %>% drop_na() #--For 5 samples there were no AI calculated and Hene Removed


#------- Activity Index using All Five Genome set and measures has been Calculated
#--- Merge All Five
ab <- merge(metaphln_activity_index, binsPTR_activity_index, by=0, all=F)
rownames(ab) <- ab$Row.names; ab$Row.names <- NULL
abc <- merge(ab, binsRelab_activity_index, by=0, all=F) # Add third One.
rownames(abc) <- abc$Row.names; abc$Row.names <- NULL
abcd <- merge(abc, Ref_PTR_activity_index, by=0, all=F) # Add Fourth One.
rownames(abcd) <- abcd$Row.names; abcd$Row.names <- NULL
AllSets_Activity_Index <- merge(abcd, Ref_Relab_activity_index, by=0, all=F) # Add Last one and make a Final File
rownames(AllSets_Activity_Index) <- AllSets_Activity_Index$Row.names; AllSets_Activity_Index$Row.names <- NULL

#-- Change column names
colnames(AllSets_Activity_Index) <- c("metaphln", "binsPTR", "binsRelab", "Ref_PTR", "Ref_Relab")

#----> Make Plot for C3001

#Merge metadata with Disease activity Index
AllSets_Activity_Indexmeta <- merge(AllSets_Activity_Index, meta, by=0, all=F)
rownames(AllSets_Activity_Indexmeta) <- AllSets_Activity_Indexmeta$Row.names; AllSets_Activity_Indexmeta$Row.names <- NULL

C3001_diseaseAI <- subset(AllSets_Activity_Indexmeta, Participant.ID=="C3001")
range_d <- range (list(C3001_diseaseAI[,1:5])) #-- Only for Bins PTR and Bins Relabundances

#--- Plot Dysbiosis scores calculated on All five measures
#jpeg("PatientSpceific_output/C3001_Dysbiosis_score.jpg", height = 5, width = 7, units = 'in', res = 600)
#plot(C3001_diseaseAI$visit_num, C3001_diseaseAI$metaphln, type = "b", frame = FALSE, pch = 19, 
    # col = "darkgray", xlab = "visits", ylab = "Dysbiosis Score", ylim = c(0.45, range_d[2]+0.1)) 
#lines(C3001_diseaseAI$visit_num, C3001_diseaseAI$binsPTR, pch = 19, col = "darkgoldenrod4", type = "b", lty = 1)
#lines(C3001_diseaseAI$visit_num, C3001_diseaseAI$binsRelab, pch = 19, col = "Maroon", type = "b", lty = 1)
#lines(C3001_diseaseAI$visit_num, C3001_diseaseAI$Ref_PTR, pch = 19, col = "Magenta", type = "b", lty = 1)
#lines(C3001_diseaseAI$visit_num, C3001_diseaseAI$Ref_Relab, pch = 19, col = "black", type = "b", lty = 1)
# Add a legend to the plot
#legend("topright", legend=c("MetaPhlAn_Relabun","MAGs_PTRs","MAGs_Relabun","Ref.Genome_PTRs","Ref.Genome_Relabun"),
    #   col=c("darkgray", "darkgoldenrod4", "Maroon", "Magenta", "black"),
     #  lty = 0.5, cex=0.5, ncol = 5:5, pch=19)
#dev.off ()

#--- Plot Dysbiosis scores calculated only on binsPTRs and bins Relative abundances
range_d <- range (list(C3001_diseaseAI[,2:3]))
jpeg("PatientSpceific_output/C3001_Dysbiosis_score.jpg", height = 5, width = 7, units = 'in', res = 600)
plot(C3001_diseaseAI$visit_num, C3001_diseaseAI$binsPTR, type = "b", frame = FALSE, pch = 19, 
     col = "darkgoldenrod4", xlab = "visits", ylab = "Dysbiosis Score", ylim = c(range_d[1]-0.2, range_d[2]+0.1)) 
lines(C3001_diseaseAI$visit_num, C3001_diseaseAI$binsRelab, pch = 19, col = "Maroon", type = "b", lty = 1)
# Add a legend to the plot
legend("topright", legend=c("MAGs_PTRs","MAGs_Relabun"),
       col=c("darkgoldenrod4", "Maroon"),
       lty = 0.5, cex=0.5, ncol = 2:2, pch=19)
abline(h=binsPTR_disease_activity_threshold, col="darkgoldenrod4", lwd=1, lty=2)
abline(h=BinsRelab_disease_activity_threshold, col="Maroon", lwd=1, lty=2)
dev.off ()

######################################################
#------- Now Plot log2 PTR of Bins belonging to C3001
######################################################
C3001_filtered_meta <- merge(C3001_filtered, meta, by=0, all=F)
rownames(C3001_filtered_meta) <- C3001_filtered_meta$Row.names; C3001_filtered_meta$Row.names <- NULL
C3001binslist <- colnames(C3001_filtered)
range_ptr = range (list (C3001_filtered))

jpeg("PatientSpceific_output/C3001_logPTR.jpg", height = 5, width = 7, units = 'in', res = 600)
# Create a first line
plot(C3001_filtered_meta$visit, C3001_filtered_meta[,C3001binslist[[1]]], type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "visits", ylab = "log2(PTR)", ylim = c(range_ptr [1], range_ptr [2] + 0.7))
# Add other lines
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,2], pch = 19, col = "blue", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,3], pch = 19, col = "orange", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,4], pch = 19, col = "darkred", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,5], pch = 19, col = "cyan", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,6], pch = 19, col = "darkgray", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,7], pch = 19, col = "Magenta", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,8], pch = 19, col = "Maroon", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,9], pch = 19, col = "limegreen", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,10], pch = 19, col = "darkgreen", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,11], pch = 19, col = "purple", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,12], pch = 19, col = "navy", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,13], pch = 19, col = "darkgoldenrod4", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta[,14], pch = 19, col = "black", type = "b", lty = 1)
# Add a legend to the plot
legend("topright", legend=C3001binslist,
       col=c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"), 
       lty = 1, cex=0.8, ncol = 5:5, pch=19)
#legend("topright", legend=c("Bin1049", "Bin12", "Bin18", "Bin2018", "Bin21", "Bin23", "Bin2362", "Bin24", "Bin2402", "Bin25", "Bin3897", "Bin659", "Bin687", "Bin7"),
# col=c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"), 
# lty = 1, cex=0.8, ncol = 5:5, pch=19)
dev.off ()

################################################################
#------- Now Plot Relative Abundance of Bins belonging to C3001
###############################################################
BinsRelab_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/RelativeAbun_analysis/Relabun_meta_no_filtered.txt", sep="\t", header = T, row.names =1)
BinsRelab <- BinsRelab_meta[C3001binslist] #-- 14 Bins for which PTR was calculated in 50% of the samples

BinsRelab_meta <- merge(BinsRelab, meta, by=0, all=F)
rownames(BinsRelab_meta) <- BinsRelab_meta$Row.names; BinsRelab_meta$Row.names <- NULL

C3001_relab_filtered_meta <- subset(BinsRelab_meta, Participant.ID =="C3001")
range <- range(list(C3001_relab_filtered_meta[,1:ncol(BinsRelab)]))


#---- Sel relative abundance of only Bins for Which Bacterial Growth Rate were calculated (at least in 50% of the samples)

jpeg("PatientSpceific_output/C3001_Relab.jpg", height = 5, width = 7, units = 'in', res = 600)
# Create a first line
plot(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,1], type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "visits", ylab = "Relative Abundance", ylim = c(range[1], range[2]+0.15))
# Add other lines
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,2], pch = 19, col = "blue", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,3], pch = 19, col = "orange", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,4], pch = 19, col = "darkred", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,5], pch = 19, col = "cyan", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,6], pch = 19, col = "darkgray", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,7], pch = 19, col = "Magenta", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,8], pch = 19, col = "Maroon", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,9], pch = 19, col = "limegreen", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,10], pch = 19, col = "darkgreen", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,11], pch = 19, col = "purple", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,12], pch = 19, col = "navy", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,13], pch = 19, col = "darkgoldenrod4", type = "b", lty = 1)
lines(C3001_relab_filtered_meta$visit, C3001_relab_filtered_meta[,14], pch = 19, col = "black", type = "b", lty = 1)
# Add a legend to the plot
legend("topright", legend=C3001binslist,
       col=c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"), 
       lty = 1, cex=0.8, ncol = 5:5, pch=19)
dev.off ()


################################################################
#------- Now Plot Sequencing Depth belonging to C3001
###############################################################
#C3001_filtered_meta has been generated by comibining Filtered samples with Metadata
range_depth = range (C3001_filtered_meta$reads_filtered)

jpeg("PatientSpceific_output/C3001_filteredReads_depth.jpg", height = 5, width = 7, units = 'in', res = 600)
# Create a first line
plot(C3001_filtered_meta$visit, C3001_filtered_meta$reads_filtered, type = "b", frame = FALSE, pch = 19, 
     col = "darkgray", xlab = "visits", ylab = "Sequencing Depth (Filtered Reads)", ylim = c(range_depth [1], range_depth [2]))
dev.off ()


#################################################################
#-- Correlation between Sequencing Depth and Bins PTR
#################################################################
C3001_corr <- cor(C3001_filtered_meta[,1:ncol(C3001_filtered)],C3001_filtered_meta$reads_filtered, method = "spearman")
cor_range <- range(C3001_corr)
jpeg("PatientSpceific_output/C3001_ptr_corr_depth.jpg", height = 5, width = 7, units = 'in', res = 600)
plot(C3001_corr, pch=19, cex=1.5, 
     col=c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"),
     xlab = "", ylab = "Spearmean's Corr", ylim = c(cor_range[1], cor_range[2]+0.1),
     main = "Correlation between sequencing depth vs Bins PTRs")
text(C3001_corr, labels=rownames(C3001_corr), cex=0.8, font=2, pos=3)
dev.off ()

#####################################################################
#-- Correlation between Sequencing Depth and Bins Relative Abundance
#####################################################################
C3001_rel_corr <- cor(C3001_relab_filtered_meta[,1:ncol(BinsRelab)],C3001_relab_filtered_meta$reads_filtered, method = "spearman")
rel_cor_range <- range(C3001_rel_corr)
jpeg("PatientSpceific_output/C3001_relab_corr_depth.jpg", height = 5, width = 7, units = 'in', res = 600)
plot(C3001_rel_corr, pch=19, cex=1.5, 
     col=c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"),
     xlab = "", ylab = "Spearmean's Corr", ylim = c(rel_cor_range[1], rel_cor_range[2]+0.1),
     main = "Correlation between sequencing depth vs Bins Relabun")
text(C3001_rel_corr, labels=rownames(C3001_rel_corr), cex=0.8, font=2, pos=3)
dev.off ()

############################################################################
# Correlation between Dysbiosis Scores and bins PTRs and Relative Abundances
############################################################################
#----- Correlation of Dysbiosis Score calculated on PTRs with Bins PTRs and on relative abundannces with Bins Relabtive abundnaces
cor_ptr <- cor(C3001_diseaseAI$binsPTR, as.matrix(C3001_filtered[,1:ncol(C3001_filtered)]),  method="pearson")
cor_relab <- cor(C3001_diseaseAI$binsRelab, as.matrix(C3001_relab_filtered_meta[,1:ncol(BinsRelab)]), method="pearson")
corr_new <- rbind (cor_ptr, cor_relab)
row.names(corr_new) <- c("Dysbiosis score (Bins PTRs)", "Dysbiosis score (Bins Relabun)")
#corrplot::corrplot(corr_new)
jpeg("PatientSpceific_output/C3001_binPTR_Corr_dysbiosis.jpg", height = 5, width = 7, units = 'in', res = 600)
barplot(corr_new[1,], angle = 45, col=c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", 
                                        "darkgreen", "purple", "navy", "darkgoldenrod4", "black"), ylim = c(-1, 1), 
        ylab = "Correlation: dysbiosis scores and MAGs PTRs")
dev.off ()
#-- in PTR - Bin8, Bin24, and Bin25 Showd significant correlations
jpeg("PatientSpceific_output/C3001_binRelab_Corr_dysbiosis.jpg", height = 5, width = 7, units = 'in', res = 600)
barplot(corr_new[2,], angle = 45, col=c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", 
                                        "darkgreen", "purple", "navy", "darkgoldenrod4", "black"), ylim = c(-1, 1), 
        ylab = "Correlation: dysbiosis scores and MAGs Relabun")
dev.off ()

############################################################################
# Correlation between Dysbiosis Scores and RefGenomes PTRs and Relative Abundances
############################################################################
C3001_ref_PTR <- subset(Ref_PTR_meta, Participant.ID=="C3001")
C3001_ref_PTR_filtered <- dropspc(C3001_ref_PTR[,1:188], 4) #-- Remove RefGenomes which is not present in Even 50% of the samples

C3001_ref_relab <- subset(Ref_Relab_meta, Participant.ID=="C3001")
C3001_ref_relab_filtered <- dropspc(C3001_ref_relab[,1:188], 4) 

C3001_ref_PTR_cor_dysbiosis <- cor(C3001_diseaseAI$binsPTR, as.matrix(C3001_ref_PTR_filtered[,1:ncol(C3001_ref_PTR_filtered)]),  method="pearson")
C3001_ref_PTR_cor_dysbiosis <- data.frame(t(C3001_ref_PTR_cor_dysbiosis))
C3001_ref_relab_cor_dysbiosis <- cor(C3001_diseaseAI$binsPTR, as.matrix(C3001_ref_relab_filtered[,1:ncol(C3001_ref_relab_filtered)]),  method="pearson")
C3001_ref_relab_cor_dysbiosis <- data.frame(t(C3001_ref_relab_cor_dysbiosis))

Ref_Genome_PTR_relab_Corrs <- merge(C3001_ref_PTR_cor_dysbiosis, C3001_ref_relab_cor_dysbiosis, by=0, all=T)
rownames(Ref_Genome_PTR_relab_Corrs) <- Ref_Genome_PTR_relab_Corrs$Row.names; Ref_Genome_PTR_relab_Corrs$Row.names <- NULL
colnames(Ref_Genome_PTR_relab_Corrs) <- c("Dysbiosis score (RefGenome PTRs)", "Dysbiosis score (RefGenome Relabun)")

#---Load Ref. Genome ID and Taxa Names
Refgid_gname <- read.csv("../PTR_analysis_on_IGG_Ref/RefGid_Gnames.txt", sep = "\t", header = T, row.names = 1)

rownames(Ref_Genome_PTR_relab_Corrs) <- gsub(".hgm", "", rownames(Ref_Genome_PTR_relab_Corrs), fixed=TRUE)
rownames(Ref_Genome_PTR_relab_Corrs) <- gsub(".patric", "", rownames(Ref_Genome_PTR_relab_Corrs), fixed=TRUE)
rownames(Ref_Genome_PTR_relab_Corrs) <- gsub(".img", "", rownames(Ref_Genome_PTR_relab_Corrs), fixed=TRUE)
rownames(Ref_Genome_PTR_relab_Corrs) <- gsub("X", "", rownames(Ref_Genome_PTR_relab_Corrs), fixed=TRUE)

Ref_Genome_PTR_relab_Corrs_TaxaNames <- merge(Ref_Genome_PTR_relab_Corrs, Refgid_gname, by=0, all=F)
rownames(Ref_Genome_PTR_relab_Corrs_TaxaNames) <- Ref_Genome_PTR_relab_Corrs_TaxaNames$Row.names; Ref_Genome_PTR_relab_Corrs_TaxaNames$Row.names <- NULL

write.table(Ref_Genome_PTR_relab_Corrs_TaxaNames, "PatientSpceific_output/C3001_Ref_Genome_PTR_relab_Corrs.txt", sep = "\t")

#-- Corr r>0.6 for PTR
colnames(Ref_Genome_PTR_relab_Corrs_TaxaNames) <- c("PTR", "Relab", "taxa")
taxa_C3001_PTR <- subset(Ref_Genome_PTR_relab_Corrs_TaxaNames, abs(PTR) > 0.6)
taxa_C3001_relab <- subset(Ref_Genome_PTR_relab_Corrs_TaxaNames, abs(Relab) > 0.6)
