#############################################
#3--Activity Index on Bins Relative Abundances
#############################################
#-- Load Bins Relative abundances to Calculate Activity Index
#Calculated These Indexes on Same number of Bins for which PTRs were calcualted in at least 5% of samples that means 70 Bins
BinsRelab_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/RelativeAbun_analysis/Relab_meta.txt", sep="\t", header = T, row.names =1)
bins_keeps <- c("Bin104","Bin1049","Bin12","Bin1258","Bin1265","Bin1291","Bin1393","Bin1479","Bin151","Bin1549","Bin1555","Bin1560","Bin1633","Bin1677","Bin168","Bin1712","Bin175","Bin1811","Bin1818","Bin1898","Bin1908","Bin1962","Bin1966","Bin197","Bin1970","Bin2018","Bin2019","Bin2023","Bin2057","Bin2075","Bin2106","Bin2132","Bin2212","Bin2332","Bin2333","Bin2338","Bin2365","Bin2434","Bin2480","Bin2483","Bin2582","Bin282","Bin2837","Bin3052","Bin3091","Bin3140","Bin3317","Bin3329","Bin3340","Bin3347","Bin3362","Bin3587","Bin3685","Bin3733","Bin3746","Bin3838","Bin3842","Bin387","Bin390","Bin3955","Bin4039","Bin4103","Bin418","Bin421","Bin584","Bin620","Bin627","Bin645","Bin687","Bin780")
BinsRelab <- BinsRelab_meta[bins_keeps]

#--meta
meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)

#Merge new metadata
BinsRelab_meta <- merge(BinsRelab, meta, by=0, all=F)
rownames(BinsRelab_meta) <- BinsRelab_meta$Row.names; BinsRelab_meta$Row.names <- NULL

#--- No need to Calculate Relative abundances on PTR data
library(vegan)
binsRelab_D <- as.matrix(vegdist(BinsRelab_meta[,1:70], method="bray"))

# The "reference set" of inactive samples
binsRelab_ref_set <- (BinsRelab_meta$diagnosis == "nonIBD") & (BinsRelab_meta$week_num >= 20)

# Calculate the activity index
BinsRelab_meta$activity_index <- sapply(seq_along(binsRelab_ref_set), function(i)
  median(binsRelab_D[i, binsRelab_ref_set & (BinsRelab_meta$site_sub_coll != BinsRelab_meta$site_sub_coll[i])]))

# Threshold activity
BinsRelab_disease_activity_threshold <- quantile(BinsRelab_meta$activity_index[BinsRelab_meta$diagnosis=="nonIBD"], 0.9)
BinsRelab_eubiosis_lower_threshold <- quantile(BinsRelab_meta$activity_index[BinsRelab_meta$diagnosis=="nonIBD"], 0.1)
BinsRelab_meta$active <- BinsRelab_meta$activity_index >= BinsRelab_disease_activity_threshold

binsRelab_activity_index <- data.frame (BinsRelab_meta$activity_index)
rownames(binsRelab_activity_index) <- row.names(BinsRelab_meta)
library(tidyr)
binsRelab_activity_index <- binsRelab_activity_index %>% drop_na()