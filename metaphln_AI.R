#-- This script will be Used to Calculate Activity Index on Relative Abundances of All Bugs- Calculated by Authors in Nature.
#-- Check in our own Script
bugs <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/hmp2_analysis/taxonomic_profiles.tsv", sep="\t", row.names = 1, header = T)
bugs <- data.frame(t(bugs))
colnames(bugs)
#--- Get Bacteria and Till species level annotation
bugs <- bugs[, grepl("k__Bacteria", colnames(bugs))]
bugs <- bugs[, grepl("s__", colnames(bugs))]
#-- All Bugs should be Present in Atleast 5% of the samples (82 out of 1638 samples)
bugs_filtered <- dropspc(bugs, 82)

#--meta
meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)

bugs_meta <- merge(bugs_filtered, meta, by=0, all=F)
rownames(bugs_meta) <- bugs_meta$Row.names; bugs_meta$Row.names <- NULL

library(vegan)
#bugs_relab <- decostand(bugs_meta[,1:975], method = "total")*100 #I haven't used all Bugs
bugs_relab <- decostand(bugs_meta[,1:314], method = "total")*100
metaphln_D <- as.matrix(vegdist(bugs_relab, method="bray"))

# The "reference set" of inactive samples
metaphln_ref_set <- (bugs_meta$diagnosis == "nonIBD") & (bugs_meta$week_num >= 20)

# Calculate the activity index
bugs_meta$activity_index <- sapply(seq_along(metaphln_ref_set), function(i)
  median(metaphln_D[i, metaphln_ref_set & (bugs_meta$site_sub_coll != bugs_meta$site_sub_coll[i])]))

# Threshold activity
disease_activity_threshold <- quantile(bugs_meta$activity_index[bugs_meta$diagnosis=="nonIBD"], 0.9)
eubiosis_lower_threshold <- quantile(bugs_meta$activity_index[bugs_meta$diagnosis=="nonIBD"], 0.1)
bugs_meta$active <- bugs_meta$activity_index >= disease_activity_threshold

metaphln_activity_index <- data.frame (bugs_meta$activity_index)
rownames(metaphln_activity_index) <- row.names(bugs_meta)
library(tidyr)
metaphln_activity_index <- metaphln_activity_index %>% drop_na() #--For 5 samples there were no AI calculated and Hene Removed

