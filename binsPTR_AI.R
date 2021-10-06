##############################
#2.--Activity Index on Bins binsPTRs
#############################
#Calculated These Indexes on Final BINs which are present consistently across samples.
#binsPTR_meta_no_filtered <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta_no_filtered.txt", sep="\t", header = T, row.names =1)
binsPTR_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.txt", sep="\t", header = T, row.names =1) #70 Bins present in atleast 5% of the samples
binsPTR <- binsPTR_meta[,1:70]
#binsPTR <- binsPTR_meta_no_filtered[,1:1503]

#--meta
meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)

#Merge new metadata
binsPTR_meta <- merge(binsPTR, meta, by=0, all=F)
rownames(binsPTR_meta) <- binsPTR_meta$Row.names; binsPTR_meta$Row.names <- NULL

#--- No need to Calculate Relative abundances on binsPTR data
library(vegan)
binsPTR_D <- as.matrix(vegdist(binsPTR_meta[,1:70], method="bray"))
#bins_D <- as.matrix(vegdist(binsPTR_meta[,1:1503], method="bray"))

# The "reference set" of inactive samples
binsPTR_ref_set <- (binsPTR_meta$diagnosis == "nonIBD") & (binsPTR_meta$week_num >= 20)

# Calculate the activity index
binsPTR_meta$activity_index <- sapply(seq_along(binsPTR_ref_set), function(i)
  median(binsPTR_D[i, binsPTR_ref_set & (binsPTR_meta$site_sub_coll != binsPTR_meta$site_sub_coll[i])]))

# Threshold activity
binsPTR_disease_activity_threshold <- quantile(binsPTR_meta$activity_index[binsPTR_meta$diagnosis=="nonIBD"], 0.9)
binsPTR_eubiosis_lower_threshold <- quantile(binsPTR_meta$activity_index[binsPTR_meta$diagnosis=="nonIBD"], 0.1)
binsPTR_meta$active <- binsPTR_meta$activity_index >= binsPTR_disease_activity_threshold

binsPTR_activity_index <- data.frame (binsPTR_meta$activity_index)
rownames(binsPTR_activity_index) <- row.names(binsPTR_meta)
library(tidyr)
binsPTR_activity_index <- binsPTR_activity_index %>% drop_na()
