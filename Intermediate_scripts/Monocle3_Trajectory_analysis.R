# For Monocle3 installation
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
 #                      'HDF5Array', 'terra', 'ggrastr'))

#remotes::install_github("rspatial/terra")
#renv::install("github::rspatial/terra")
#devtools::install_github('cole-trapnell-lab/monocle3')
#devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")

library(monocle3)
#-- Example Set
#expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds"))
#cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_colData.rds"))
#gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_rowData.rds"))

#- Get metadata
PTR_meta <- read.csv (file="~/Work/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.txt", sep="\t", row.names =1, header = T)
sample_metadata <- PTR_meta[,71:90]
#sample_metadata <- new("AnnotatedDataFrame", data = sample_metadata) # Don't use this step in Monocle3

#- Get Ref Gene Counts
Ref_counts <- read.csv(file="../RefGenome_Counts_output.csv", row.names = 1, header = T, sep = ",")
# Fix names
rownames(Ref_counts) <- gsub(".hgm", "", rownames(Ref_counts), fixed=TRUE)
rownames(Ref_counts) <- gsub(".patric", "", rownames(Ref_counts), fixed=TRUE)
rownames(Ref_counts) <- gsub(".img", "", rownames(Ref_counts), fixed=TRUE)
rownames(Ref_counts) <- gsub("X", "", rownames(Ref_counts), fixed=TRUE)
samples <- row.names(PTR_meta)
Ref_counts_Sel <- Ref_counts[,samples] #-- to keep same samples as we have in the metadata file

#-- Check Sparsity(i.e. the percentage of 0 values) of a Dataframe
sum(Ref_counts_Sel == 0)/(dim(Ref_counts_Sel)[1]*dim(Ref_counts_Sel)[2])
#[1] 0.9082249: 90%


#- Get gene name File (here it is taxa names)
Refgid_gname_Org <- read.csv("../../PTR_analysis_on_IGG_Ref/RefGid_Gnames.txt", sep = "\t", header = T, row.names = 1)
Refgid_gname_Ref_counts <- merge(Refgid_gname_Org, Ref_counts_Sel, by=0, all=F)
rownames(Refgid_gname_Ref_counts) <- Refgid_gname_Ref_counts$Row.names; Refgid_gname_Ref_counts$Row.names <- NULL
Refgid_gname <- Refgid_gname_Ref_counts[,1] #-- Just to get same names as we have in Expression matrix
Refgid_gname <- data.frame(Refgid_gname)
row.names(Refgid_gname) <- row.names(Refgid_gname_Ref_counts)
colnames(Refgid_gname) <- "gene_short_name"
#-- If taxa names to be made Unique
#write.table(Refgid_gname, file="Refgid_gname.txt", sep="\t")
#_ These steps to reanme duplicate Taxa Names
Refgid_gname_uniq <- make.unique(Refgid_gname$gene_short_name)
Refgid_gname_uniq <- data.frame(Refgid_gname_uniq)
Refgid_gname_updated <- cbind(Refgid_gname, Refgid_gname_uniq)
Refgid_gname_updated <- data.frame(Refgid_gname_updated[,2])
row.names(Refgid_gname_updated) <- row.names(Refgid_gname)
colnames(Refgid_gname_updated) <- "gene_short_name"
#gene_annotation <- new("AnnotatedDataFrame", data = Refgid_gname_updated) # Don't use this Step in Monocle3


Ref_counts_Dataset <- new_cell_data_set(as.matrix(Ref_counts_Sel),
                         cell_metadata = sample_metadata,
                         gene_metadata = Refgid_gname_updated)


# Pre-process the data
Ref_counts_Dataset <- preprocess_cds(Ref_counts_Dataset, num_dim = 100)
#It's a good idea to check that you're using enough PCs to capture most of the variation in taxa abundance across 
#all samples in the data set. You can look at the fraction of variation explained by each PC using plot_pc_variance_explained():
plot_pc_variance_explained(Ref_counts_Dataset)

#Reduce dimensionality and visualize the cells
Ref_counts_Dataset <- reduce_dimension(Ref_counts_Dataset)
plot_cells(Ref_counts_Dataset)
plot_cells(Ref_counts_Dataset, color_cells_by="diagnosis")
plot_cells(Ref_counts_Dataset, genes=c("Roseburia intestinalis", "Subdoligranulum sp.", "Flavonifractor plautii"))

#Reduce dimensions using tSNE
Ref_counts_Dataset <- reduce_dimension(Ref_counts_Dataset, reduction_method="tSNE")
plot_cells(Ref_counts_Dataset, reduction_method="tSNE", color_cells_by="diagnosis")


#Group cells into clusters
#Grouping cells into clusters is an important step in identifying the cell types represented in your data. Monocle uses a technique called community detection to group cells. This approach was introduced by Levine et al as part of the phenoGraph algorithm. You can cluster your cells using the cluster_cells() function, like this:
  
Ref_counts_Dataset <- cluster_cells(Ref_counts_Dataset, resolution=1e-5)
plot_cells(Ref_counts_Dataset)

plot_cells(Ref_counts_Dataset, color_cells_by="diagnosis", group_cells_by="partition")


#--- -Trajectories
Ref_counts_Dataset <- new_cell_data_set(as.matrix(Ref_counts_Sel),
                                        cell_metadata = sample_metadata,
                                        gene_metadata = Refgid_gname)
#Ref_counts_Dataset <- detectGenes(Ref_counts_Dataset, min_expr = 0.1)

Ref_counts_Dataset <- preprocess_cds(Ref_counts_Dataset, num_dim = 50)

Ref_counts_Dataset <- reduce_dimension(Ref_counts_Dataset)
plot_cells(Ref_counts_Dataset, label_groups_by_cluster=FALSE,  color_cells_by = "diagnosis")
