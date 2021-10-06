########################################################################
#4--Activity Index on PTRs Calculated on Referene Genomes
########################################################################
#Calculated These Indexes on Final BINs which are present consistently across samples.
Ref_PTR_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/PTR_analysis/PTR_meta.txt", sep="\t", header = T, row.names =1)
Ref_PTR <- Ref_PTR_meta[,1:188]

#--meta
meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)

#Merge new metadata
Ref_PTR_meta <- merge(Ref_PTR, meta, by=0, all=F)
rownames(Ref_PTR_meta) <- Ref_PTR_meta$Row.names; Ref_PTR_meta$Row.names <- NULL

#--- No need to Calculate Relative abundances on PTR data
library(vegan)
Ref_PTR_D <- as.matrix(vegdist(Ref_PTR_meta[,1:188], method="bray"))

# The "reference set" of inactive samples
Ref_PTR_ref_set <- (Ref_PTR_meta$diagnosis == "nonIBD") & (Ref_PTR_meta$week_num >= 20)

# Calculate the activity index
Ref_PTR_meta$activity_index <- sapply(seq_along(Ref_PTR_ref_set), function(i)
  median(Ref_PTR_D[i, Ref_PTR_ref_set & (Ref_PTR_meta$site_sub_coll != Ref_PTR_meta$site_sub_coll[i])]))

# Threshold activity
Ref_PTR_disease_activity_threshold <- quantile(Ref_PTR_meta$activity_index[Ref_PTR_meta$diagnosis=="nonIBD"], 0.9)
Ref_PTR_eubiosis_lower_threshold <- quantile(Ref_PTR_meta$activity_index[Ref_PTR_meta$diagnosis=="nonIBD"], 0.1)
Ref_PTR_meta$active <- Ref_PTR_meta$activity_index >= Ref_PTR_disease_activity_threshold

Ref_PTR_activity_index <- data.frame (Ref_PTR_meta$activity_index)
rownames(Ref_PTR_activity_index) <- row.names(Ref_PTR_meta)
library(tidyr)
Ref_PTR_activity_index <- Ref_PTR_activity_index %>% drop_na()