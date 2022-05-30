setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts")
library (vegan)
library (ape)
library (variancePartition)

############*********************
#-- Step1: Ordination Analysis
############*********************
#####################################--------------------On MAGs PTRs
PTR_meta <- read.csv (file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.txt", sep="\t", row.names =1, header = T)

PTR_dist <-vegdist(PTR_meta[,1:70], "bray")
PTR_dist_pcoa <- pcoa(PTR_dist)
pcoa_scores <- data.frame(PTR_dist_pcoa$vectors[,1:2])
pcoa_scores <- cbind(pcoa_scores, PTR_meta[,71:90])

Color <- c("darkred", "steelblue2", "orange")
#Color <- c(#8b0000, #5cacee, #FFA500)
jpeg("../Combined_Final_Figures/MAG_PTR_PCoA.jpg", height = 5, width = 6, units = 'in', res = 600)
ggplot(pcoa_scores, aes(x=Axis.1, y=Axis.2, colour=diagnosis)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + 
  xlab("PCo1 - 17.2") + ylab("PCo2 - 11.1") + ggtitle("") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
#One measure of multivariate dispersion (variance) for a group of samples is to calculate the average distance of group members to the group centroid or spatial median 
#(both referred to as 'centroid' from now on unless stated otherwise) in multivariate space. 
beta_Disper <- betadisper(PTR_dist, PTR_meta$diagnosis, type = c("median","centroid"), bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)

jpeg("../Combined_Final_Figures/MAG_PTR_CenteroidDist.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(beta_Disper$distances ~ PTR_meta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="Distance from centroid", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8), ylim = c(0, 1))
#stripchart(beta_Disper$distances ~ PTR_meta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
kruskalmc(beta_Disper$distances ~ PTR_meta$diagnosis)
#--- All PERMANOVA Results can be seen in "Step1_ptr_analysis.R" at Original location ~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis

#####################################--------------------On MAGs Relative Abundances
Relabun_meta <- read.csv (file="~/Box/David_Casero_Lab/HMP_IBD_Project/RelativeAbun_analysis/Relab_meta.txt", sep="\t", row.names =1, header = T)
Relab_dist <-vegdist(Relabun_meta[,1:803], "bray")
Relab_dist_pcoa <- pcoa(Relab_dist)
pcoa_scores <- data.frame(Relab_dist_pcoa$vectors[,1:2])
pcoa_scores <- cbind(pcoa_scores, Relabun_meta[,804:823])

Color <- c("darkred", "steelblue2", "orange")
jpeg("../Combined_Final_Figures/MAG_Relab_PCoA.jpg", height = 5, width = 6, units = 'in', res = 600)
ggplot(pcoa_scores, aes(x=Axis.1, y=Axis.2, colour=diagnosis)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + 
  xlab("PCo1 - 6.1") + ylab("PCo2 - 5.3") + ggtitle("") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
#One measure of multivariate dispersion (variance) for a group of samples is to calculate the average distance of group members to the group centroid or spatial median 
#(both referred to as 'centroid' from now on unless stated otherwise) in multivariate space. 
beta_Disper <- betadisper(Relab_dist, Relabun_meta$diagnosis, type = c("median","centroid"), bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)

jpeg("../Combined_Final_Figures/MAG_Relab_CenteroidDist.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(beta_Disper$distances ~ Relabun_meta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="Distance from centroid", cex.axis=0.8, 
        par(cex.lab=0.8), par(cex.axis=0.8), ylim = c(0.4, 0.8))
#stripchart(beta_Disper$distances ~ Relabun_meta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
kruskalmc(beta_Disper$distances ~ Relabun_meta$diagnosis)
#--- All PERMANOVA Results can be seen in "Complete_Analysis.R" at Original location ~/Box/David_Casero_Lab/HMP_IBD_Project/RelativeAbun_analysis


#####################################--------------------On RefGenome PTRs
RefGenome_PTR_meta <- read.csv (file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/PTR_analysis/PTR_meta.txt", sep="\t", row.names =1, header = T)

RefGenome_PTR_dist <-vegdist(RefGenome_PTR_meta[,1:188], "bray")
RefGenome_PTR_dist_pcoa <- pcoa(RefGenome_PTR_dist)
pcoa_scores <- data.frame(RefGenome_PTR_dist_pcoa$vectors[,1:2])
pcoa_scores <- cbind(pcoa_scores, RefGenome_PTR_meta[,189:208])

Color <- c("darkred", "steelblue2", "orange")
jpeg("../Combined_Final_Figures/RefGenome_PTR_PCoA.jpg", height = 5, width = 6, units = 'in', res = 600)
ggplot(pcoa_scores, aes(x=Axis.1, y=Axis.2, colour=diagnosis)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + 
  xlab("PCo1 - 15.9") + ylab("PCo2 - 8.6") + ggtitle("") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
#One measure of multivariate dispersion (variance) for a group of samples is to calculate the average distance of group members to the group centroid or spatial median 
#(both referred to as 'centroid' from now on unless stated otherwise) in multivariate space. 
beta_Disper <- betadisper(RefGenome_PTR_dist, RefGenome_PTR_meta$diagnosis, type = c("median","centroid"), bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)
jpeg("../Combined_Final_Figures/RefGenome_PTR_CenteroidDist.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(beta_Disper$distances ~ RefGenome_PTR_meta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="Distance from centroid", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8), ylim = c(0, 1))
#stripchart(beta_Disper$distances ~ RefGenome_PTR_meta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
kruskalmc(beta_Disper$distances ~ RefGenome_PTR_meta$diagnosis)
#--- All PERMANOVA Results can be seen in "PTR_analysis.R" at Original location ~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref

#####################################--------------------On RefGenome Relative Abundances
RefGenome_Relabun_meta <- read.csv (file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/Relab_analysis/Relabun_met.txt", sep="\t", row.names =1, header = T)
RefGenome_Relabun_dist <-vegdist(RefGenome_Relabun_meta[,1:419], "bray")
RefGenome_Relabun_dist_pcoa <- pcoa(RefGenome_Relabun_dist)
pcoa_scores <- data.frame(RefGenome_Relabun_dist_pcoa$vectors[,1:2])
pcoa_scores <- cbind(pcoa_scores, RefGenome_Relabun_meta[,420:439])

Color <- c("darkred", "steelblue2", "orange")
jpeg("../Combined_Final_Figures/RefGenome_Relab_PCoA.jpg", height = 5, width = 6, units = 'in', res = 600)
ggplot(pcoa_scores, aes(x=Axis.1, y=Axis.2, colour=diagnosis)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + 
  xlab("PCo1 - 14.2") + ylab("PCo2 - 9.6") + ggtitle("") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
#One measure of multivariate dispersion (variance) for a group of samples is to calculate the average distance of group members to the group centroid or spatial median 
#(both referred to as 'centroid' from now on unless stated otherwise) in multivariate space. 
beta_Disper <- betadisper(RefGenome_Relabun_dist, RefGenome_Relabun_meta$diagnosis, type = c("median","centroid"), bias.adjust = FALSE, sqrt.dist = FALSE, add = FALSE)

jpeg("../Combined_Final_Figures/RefGenome_Relab_CentroidDist.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(beta_Disper$distances ~ RefGenome_Relabun_meta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="Distance from centroid", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8), ylim = c(0, 1))
#stripchart(beta_Disper$distances ~ RefGenome_Relabun_meta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
kruskalmc(beta_Disper$distances ~ RefGenome_Relabun_meta$diagnosis)
#--- All PERMANOVA Results can be seen in "Relative_abund_ana.R" at Original location ~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref

#-- PERANOVA Figures has been Created in "~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/Permanova_Variance.xlsx
#--- If any additional PERMANOVA Needs to be Added

##################################################***********************
#-- Step2: Calculate Dysbiosys Score and Use them For Variance Partition
##################################################***********************
#####################################-------------------- MAGs PTRs
setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")
#---- Calculate Disease Activity Index on All sets - Bins, RefGenomes - PTRs and Relative abundnaces, and Metaphln relab also
source("../Disease_Activity/Scripts/binsPTR_AI.R")  #Calculated on PTR estimates of 70 Bins (present in atleast 5% of the samples)
#-- Load BinPTR_metaAgain (Metadata used in Disease Activity file is not formatted and include weekNumber and site_sub_coll - Important for calculations)
binsPTR_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.txt", sep="\t", header = T, row.names =1) #70 Bins present in atleast 5% of the samples
#-- Now merge Activity (Dysbiosys Score with Original File)
binsPTR_meta_DysbiosisScores <- merge(binsPTR_meta, binsPTR_activity_index, by=0, all=F)
rownames(binsPTR_meta_DysbiosisScores) < binsPTR_meta_DysbiosisScores$Row.names; binsPTR_meta_DysbiosisScores$Row.names <- NULL
#write.table(binsPTR_meta_DysbiosisScores, file= "~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta_dysbiosisScores.txt", sep="\t")

#--VariancePartiton (Run on HPC Cluster)
# scp -r ~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta_dysbiosisScores.txt sharmaa4@csclprd3-s001v:/home/sharmaa4/IBD_datasets/HMP2/Analysis_R/
# qsub -q 1tb.q -cwd -o $PWD -e $PWD -l mem_free=400G -pe mpich 12 run_R.sh
# 
library (variancePartition)
PTR_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta_dysbiosisScores.txt", sep="\t", header = T, row.names =1)
PTR_filtered <- data.frame(t(PTR_meta[,1:70]))
meta_filtered <- PTR_meta[,71:91] #91 Column is Dysbiosis Score
model_new <- ~ binsPTR_meta.activity_index + (1|diagnosis)+(1|sex)+(1|site)+(1|alchohol)+(1|antibiotics)+(1|immunosuppressants)+(1|chemotherapy)+(1|bowel_surgery)+(1|vegetables)+(1|probiotics)+(1|yogurt)+(1|fruits)+(1|pid)

var_part_model <- fitExtractVarPartModel(PTR_filtered, model_new,meta_filtered)
#save(var_part_model, file="../Combined_Final_Models/MAG_PTRs_VarPart.Rdata")

dim(var_part_model)
write.table(var_part_model, file="../Combined_Final_Models/MAG_PTR_varPartModel.txt", sep = "\t")
#sort the terms of the model by average variance explained across all genes, so when we plot they will be sorted by overall importance:
vp2 <- sortCols(var_part_model)
# Violin plot
jpeg("../Combined_Final_Figures/MAG_PTR_varPartModel.jpg", height = 7, width = 7, units = 'in', res = 600)
plotVarPart( vp2  ,label.angle = 90)
dev.off ()
#sort genes based on variance explained by Condition (#Here it is a treatment)
head(var_part_model[order(var_part_model$diagnosis, decreasing=TRUE),])

#####################################-------------------- MAGs Relative Abundnaces
source("../Disease_Activity/Scripts/binsRelab_AI.R") #Calculated on Relative abundance estimates of Above 70 Bins 
#-- Load Bins Relab Metadata fiel Again (Metadata used in Disease Activity file is not formatted and include weekNumber and site_sub_coll - Important for calculations)
BinsRelab_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/RelativeAbun_analysis/Relab_meta.txt", sep="\t", header = T, row.names =1) #70 Bins present in atleast 5% of the samples
#-- Now merge Activity (Dysbiosys Score with Original File)
binsRelab_meta_DysbiosisScores <- merge(BinsRelab_meta, binsRelab_activity_index, by=0, all=F)
rownames(binsRelab_meta_DysbiosisScores) < binsRelab_meta_DysbiosisScores$Row.names; binsRelab_meta_DysbiosisScores$Row.names <- NULL
#write.table(binsRelab_meta_DysbiosisScores, file= "~/Box/David_Casero_Lab/HMP_IBD_Project/RelativeAbun_analysis/binsRelab_meta_DysbiosisScores.txt", sep="\t")

library (variancePartition)
Relabun_meta <- read.csv("~/Box/David_Casero_Lab/HMP_IBD_Project/RelativeAbun_analysis/binsRelab_meta_DysbiosisScores.txt", sep="\t", header = T, row.names =1)
Relabun_filtered <- data.frame(t(Relabun_meta[,1:803]))
meta_filtered <- Relabun_meta[,804:824] #824 Column is Dysbiosis Score

model_new <- ~ BinsRelab_meta.activity_index + (1|diagnosis)+(1|sex)+(1|site)+(1|alchohol)+(1|antibiotics)+(1|immunosuppressants)+(1|chemotherapy)+(1|bowel_surgery)+(1|vegetables)+(1|probiotics)+(1|yogurt)+(1|fruits)+(1|pid)
var_part_model_new <- fitExtractVarPartModel(Relabun_filtered, model_new,meta_filtered)
save(var_part_model_new, file="../Combined_Final_Models/MAG_Relabun_VarPart.Rdata")
#-----
# Load the model generated on Cluster
load ("var_part_model_new.Rdata")
dim(var_part_model_new)
write.table(var_part_model_new, file="Total_analysis/Relabun_varPart_model_new.txt", sep = "\t")
#sort the terms of the model by average variance explained across all genes, so when we plot they will be sorted by overall importance:
vp2 <- sortCols(var_part_model_new)
# Violin plot
jpeg("Total_analysis/Relabun_variance_partision_new.jpg", height = 7, width = 7, units = 'in', res = 600)
plotVarPart( vp2  ,label.angle = 90)
dev.off ()
#sort genes based on variance explained by Condition (#Here it is a treatment)
head(var_part_model_new[order(var_part_model_new$diagnosis, decreasing=TRUE),])



##################################################***********************
#-- Step3: For Dysbiosis Plots
##################################################***********************
#setwd("~/Work/David_Casero_Lab/HMP_IBD_Project/Disease_Activity/")
#-- All Dysbiosis Scores Along with their associations with Calprotectin  Level are present in 
# load ("~/Work/David_Casero_Lab/HMP_IBD_Project/Disease_Activity/disease_activity.RData")
# "Script Disease_activity.R"
# "Script Dysbiotic_Stats.R" for getting Survival Rates using different Measures


##################################################***********************
#-- Step4: Try Unsupervised Clustering on Ref Genomes PTRs
##################################################***********************
source("~/Work/David_Casero_Lab/HMP_IBD_Project/Disease_Activity/Scripts/RefPTR_AI.R")

Ref_PTR_dist_bray <- vegdist(Ref_PTR_meta[,1:188], "bray")
Ref_PTR_dist_bray <- as.matrix (Ref_PTR_dist_bray)
Ref_PTR_dist_pcoa <- pcoa(Ref_PTR_dist_bray)
Ref_PTR_pcoa_scores <- data.frame(Ref_PTR_dist_pcoa$vectors[,1:2])
Ref_PTR_pcoa_scores <- cbind(Ref_PTR_pcoa_scores, Ref_PTR_meta[,189:678])

jpeg("../Combined_Final_Figures/NumberOfOptimumClusters.jpg", height = 4, width = 5, units = 'in', res = 600)
fviz_nbclust(Ref_PTR_dist_bray, kmeans, method = "wss") + geom_vline(xintercept = 9, linetype = 2)
dev.off ()

library(cluster)
pam <- pam(Ref_PTR_dist_bray, k=10)
clusters <- pam$clustering
clusters <- data.frame(clusters)

Ref_PTR_pcoa_scores_metaClusters <- merge(clusters, Ref_PTR_pcoa_scores, by=0, all=F)
rownames(Ref_PTR_pcoa_scores_metaClusters) <- Ref_PTR_pcoa_scores_metaClusters$Row.names; Ref_PTR_pcoa_scores_metaClusters$Row.names <- NULL

Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "1"] <- "Cluster1"
Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "2"] <- "Cluster2"
Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "3"] <- "Cluster3"
Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "4"] <- "Cluster4"
Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "5"] <- "Cluster5"
Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "6"] <- "Cluster6"
Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "7"] <- "Cluster7"
Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "8"] <- "Cluster8"
Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "9"] <- "Cluster9"
Ref_PTR_pcoa_scores_metaClusters$clusters[Ref_PTR_pcoa_scores_metaClusters$clusters == "10"] <- "Cluster10"

#write.table(Ref_PTR_pcoa_scores_metaClusters, file="../Combined_Final_Figures/Ref_PTRs_Clusters.txt", sep="\t")

jpeg("../Combined_Final_Figures/Bray_kmean_MetabotypesSplot.jpg", height = 4, width = 5, units = 'in', res = 600)
s.class(Ref_PTR_dist_pcoa$vectors, fac= as.factor(Ref_PTR_pcoa_scores_metaClusters$clusters), 
        col = c("forestgreen", "wheat3", "red", "goldenrod4", "orange", "gray32", "darkgray", "blue", "khaki4", "black"), 
        label = c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8","Cluster9","Cluster10"))
dev.off ()

#-- Find Indval Scores Attach Metadata with Clusters
Clusters <- Ref_PTR_pcoa_scores_metaClusters$clusters
Clusters <- data.frame(Clusters)
rownames(Clusters) <- row.names(Ref_PTR_pcoa_scores_metaClusters)

Ref_PTR_Clusters <- merge(Clusters, Ref_PTR_meta[,1:188], by=0, all=F)
rownames(Ref_PTR_Clusters) <- Ref_PTR_Clusters$Row.names; Ref_PTR_Clusters$Row.names <- NULL

Axis_Scores <- Ref_PTR_pcoa_scores_metaClusters[,2:3]
Ref_PTR_Clusters_AxisScore <- merge(Axis_Scores, Ref_PTR_Clusters, by=0, all=F)
rownames(Ref_PTR_Clusters_AxisScore) <- Ref_PTR_Clusters_AxisScore$Row.names; Ref_PTR_Clusters_AxisScore$Row.names <- NULL

library (labdsv)
#Genus2 <- Otu_table_Clusters[ , which(!apply(Otu_table_Clusters==0,2,all))]
#Genus2[is.na(Genus2)] <- 0
iva <- indval(Ref_PTR_Clusters_AxisScore[,4:191], Ref_PTR_Clusters_AxisScore$Clusters)
gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(Ref_PTR_Clusters_AxisScore[,4:191]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
write.table (indvalsummary, file="../Combined_Final_Figures/Ref_PTR_Clusters_IndVal.txt", sep = "\t")

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#--- Cluster2 (Faecalibacterium HGM13285) Majorly nonIBD
jpeg("../Combined_Final_Figures/Clusters/Cluster2_SRS017307_13.hgm.jpg", height = 4, width = 5.5, units = 'in', res = 600)
ggplot(Ref_PTR_Clusters_AxisScore, aes(x=Axis.1, y=Axis.2, colour=SRS017307_13.hgm)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

#--- Cluster3 (Bacteroides rodentium) Mainly CD
jpeg("../Combined_Final_Figures/Clusters/Cluster3_X1236512.3.patric.jpg", height = 4, width = 5.5, units = 'in', res = 600)
ggplot(Ref_PTR_Clusters_AxisScore, aes(x=Axis.1, y=Axis.2, colour=X1236512.3.patric)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

#--- Cluster5 (Roseburia HGM12473) Mainly CD
jpeg("../Combined_Final_Figures/Clusters/Cluster5_SRS295008_15.hgm.jpg", height = 4, width = 5.5, units = 'in', res = 600)
ggplot(Ref_PTR_Clusters_AxisScore, aes(x=Axis.1, y=Axis.2, colour=SRS295008_15.hgm)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

#--- Cluster7 (Ruminococcus torques) # CD and UC
jpeg("../Combined_Final_Figures/Clusters/Cluster7_X665951.3.patric.jpg", height = 4, width = 5.5, units = 'in', res = 600)
ggplot(Ref_PTR_Clusters_AxisScore, aes(x=Axis.1, y=Axis.2, colour=X665951.3.patric)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

#--- Cluster8 (Oscillibacter sp.) # majorly nonIBD
jpeg("../Combined_Final_Figures/Clusters/Cluster8_SRS1596798_70.hgm.jpg", height = 4, width = 6, units = 'in', res = 600)
ggplot(Ref_PTR_Clusters_AxisScore, aes(x=Axis.1, y=Axis.2, colour=SRS1596798_70.hgm)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

#--- Cluster9 (Oscillospiraceae sp.) # majorly nonIBD
jpeg("../Combined_Final_Figures/Clusters/Cluster9_ERS631833_15.hgm.jpg", height = 4, width = 6, units = 'in', res = 600)
ggplot(Ref_PTR_Clusters_AxisScore, aes(x=Axis.1, y=Axis.2, colour=ERS631833_15.hgm)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

#--- Cluster10 (Phascolarctobacterium sp.) # CD and UC only
jpeg("../Combined_Final_Figures/Clusters/Cluster10_ERS537306_47.hgm.jpg", height = 4, width = 6, units = 'in', res = 600)
ggplot(Ref_PTR_Clusters_AxisScore, aes(x=Axis.1, y=Axis.2, colour=ERS537306_47.hgm)) + geom_point(size=3) + theme_bw() + scale_colour_gradientn(colours = myPalette(100)) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()


##################################################***********************
#-- Step4: Try Unsupervised Clustering on Ref Genomes Relative Abundances
##################################################***********************

