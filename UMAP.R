#---- Load data
#setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")


#--------------------------------------- ****************** ---------------------------------------
#---- RefGenome PTR
load ("~/Work/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/PTR_analysis/PTR_meta.RData")

#---- Get the Distance Matrix First
RefGenome_PTR_dist <- vegdist(PTR_meta[,1:188], method = "jaccard")
RefGenome_PTR_dist <- as.matrix(RefGenome_PTR_dist)
write.table(RefGenome_PTR_dist, file="RefGenome_PTR_dist.tsv", sep = "\t")
#- Vi remove " and put one tab in the first line
#conda activate qiime2-2021.4
#qiime tools import --input-path RefGenome_PTR_dist.tsv --output-path RefGenome_PTR_dist.qza --type DistanceMatrix

#write.table (PTR_meta[,189:208], file= "RefGenome_PTR_meta.tsv", sep = "\t")
# Format as per the Qiime and convert it in QZA

#qiime umap embed --i-distance-matrix RefGenome_PTR_dist.qza --p-n-neighbors 500 --o-umap RefGenome_PTR_umap.qza 
#qiime emperor plot --i-pcoa RefGenome_PTR_umap.qza --m-metadata-file RefGenome_PTR_meta.tsv --o-visualization RefGenome_PTR_umap-emperor.qzv

RefGenome_PTR_ord <- read.csv(file="RefGenome_PTR_umap/data/ordination_edited.txt", sep = "\t", row.names = 1, header = T)
RefGenome_PTR_ord_meta <- merge(RefGenome_PTR_ord, PTR_meta[,189:208], by=0, all=F)
rownames(RefGenome_PTR_ord_meta) <- RefGenome_PTR_ord_meta$Row.names; RefGenome_PTR_ord_meta$Row.names <- NULL

#--- Plot UMAP
Color <- c("darkred", "steelblue2", "orange")
jpeg("RefFGenomePTR_umap.jpg", height = 5, width = 6, units = 'in', res = 600)
ggplot(RefGenome_PTR_ord_meta, aes(x=UMAP1, y=UMAP2, colour=diagnosis)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + 
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

jpeg("RefGenome_PTR_UMAP1.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(RefGenome_PTR_ord_meta$UMAP1 ~ RefGenome_PTR_ord_meta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="UMAP1", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
#stripchart(beta_Disper$distances ~ PTR_meta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
kruskalmc(RefGenome_PTR_ord_meta$UMAP1 ~ RefGenome_PTR_ord_meta$diagnosis)

aggregate(RefGenome_PTR_ord_meta$UMAP1, list(RefGenome_PTR_ord_meta$diagnosis), FUN=median)
#Group.1          x
#1      CD  0.4514606
#2  nonIBD -2.1518843
#3      UC  0.3148153

#--------------------------------------- ****************** ---------------------------------------
#---- RefGenome Relative Abundances
load ("~/Work/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/Relab_analysis/Relabun_meta.RData")

#---- Get the Distance Matrix First
RefGenome_relab_dist <- vegdist(Relabun_meta[,1:419], method = "jaccard")
RefGenome_relab_dist <- as.matrix(RefGenome_relab_dist)
write.table(RefGenome_relab_dist, file="RefGenome_Relab_dist.tsv", sep = "\t")
#- Vi remove " and put one tab in the first line

#qiime tools import --input-path RefGenome_Relab_dist.tsv --output-path RefGenome_Relab_dist.qza --type DistanceMatrix 
#qiime umap embed --i-distance-matrix RefGenome_Relab_dist.qza --p-n-neighbors 500 --o-umap RefGenome_Relab_umap.qza

RefGenome_Relab_ord <- read.csv(file="RefGenome_Relab_Umap/data/ordination_edited.txt", sep = "\t", row.names = 1, header = T)
RefGenome_Relab_ord_meta <- merge(RefGenome_Relab_ord, Relabun_meta[,420:439], by=0, all=F)
rownames(RefGenome_Relab_ord_meta) <- RefGenome_Relab_ord_meta$Row.names; RefGenome_Relab_ord_meta$Row.names <- NULL

#--- Plot UMAP
Color <- c("darkred", "steelblue2", "orange")
jpeg("RefFGenomeRelab_umap.jpg", height = 5, width = 6, units = 'in', res = 600)
ggplot(RefGenome_Relab_ord_meta, aes(x=UMAP1, y=UMAP2, colour=diagnosis)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + 
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

jpeg("RefGenome_Relab_UMAP1.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(RefGenome_Relab_ord_meta$UMAP1 ~ RefGenome_Relab_ord_meta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="UMAP1", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
#stripchart(beta_Disper$distances ~ PTR_meta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
kruskalmc(RefGenome_Relab_ord_meta$UMAP1 ~ RefGenome_Relab_ord_meta$diagnosis)

aggregate(RefGenome_Relab_ord_meta$UMAP1, list(RefGenome_Relab_ord_meta$diagnosis), FUN=median)
#Group.1         x
#1      CD -0.649589
#2  nonIBD  1.147530
#3      UC -0.564796

#--------------------------------------- ****************** ---------------------------------------
#---- Bins (MAGs PTRs)
load ("~/Work/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.RData")

#---- Get the Distance Matrix First
Bins_PTR_dist <- vegdist(PTR_meta[,1:70], method = "jaccard")
Bins_PTR_dist <- as.matrix(Bins_PTR_dist)
write.table(Bins_PTR_dist, file="Bins_PTR_dist.tsv", sep = "\t")
#- Vi remove " and put one tab in the first line
#qiime tools import --input-path Bins_PTR_dist.tsv --output-path Bins_PTR_dist.qza --type DistanceMatrix
#qiime umap embed --i-distance-matrix Bins_PTR_dist.qza --p-n-neighbors 500 --o-umap Bins_PTR_umap.qza

Bins_PTR_ord <- read.csv(file="Bins_PTR_Umap/data/ordination_edited.txt", sep = "\t", row.names = 1, header = T)
Bins_PTR_ord_meta <- merge(Bins_PTR_ord, PTR_meta[,71:90], by=0, all=F)
rownames(Bins_PTR_ord_meta) <- Bins_PTR_ord_meta$Row.names; Bins_PTR_ord_meta$Row.names <- NULL

#--- Plot UMAP
Color <- c("darkred", "steelblue2", "orange")
jpeg("Bins_PTR_umap.jpg", height = 5, width = 6, units = 'in', res = 600)
ggplot(Bins_PTR_ord_meta, aes(x=UMAP1, y=UMAP2, colour=diagnosis)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + 
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

jpeg("Bins_PTR_UMAP1.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(Bins_PTR_ord_meta$UMAP1 ~ Bins_PTR_ord_meta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="UMAP1", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
#stripchart(beta_Disper$distances ~ PTR_meta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
kruskalmc(Bins_PTR_ord_meta$UMAP1 ~ Bins_PTR_ord_meta$diagnosis)

aggregate(Bins_PTR_ord_meta$UMAP1, list(Bins_PTR_ord_meta$diagnosis), FUN=median)
#Group.1          x
#1      CD  0.4714638
#2  nonIBD -1.2624006
#3      UC -0.5986924

#--------------------------------------- ****************** ---------------------------------------
#---- Bins (MAGs Relab)
load ("~/Work/David_Casero_Lab/HMP_IBD_Project/RelativeAbun_analysis/Relabun_meta.RData")

#---- Get the Distance Matrix First
Bins_Relab_dist <- vegdist(Relabun_meta[,1:803], method = "jaccard")
Bins_Relab_dist <- as.matrix(Bins_Relab_dist)
write.table(Bins_Relab_dist, file="Bins_Relab_dist.tsv", sep = "\t")

#- Vi remove " and put one tab in the first line
#qiime tools import --input-path Bins_Relab_dist.tsv --output-path Bins_Relab_dist.qza --type DistanceMatrix
#qiime umap embed --i-distance-matrix Bins_Relab_dist.qza --p-n-neighbors 500 --o-umap Bins_Relab_umap.qza

Bins_Relab_ord <- read.csv(file="Bins_relab_Umap/data/ordination_edited.txt", sep = "\t", row.names = 1, header = T)
Bins_Relab_ord_meta <- merge(Bins_Relab_ord, Relabun_meta[,804:823], by=0, all=F)
rownames(Bins_Relab_ord_meta) <- Bins_Relab_ord_meta$Row.names; Bins_Relab_ord_meta$Row.names <- NULL

#--- Plot UMAP
Color <- c("darkred", "steelblue2", "orange")
jpeg("Bins_Relab_umap.jpg", height = 5, width = 6, units = 'in', res = 600)
ggplot(Bins_Relab_ord_meta, aes(x=UMAP1, y=UMAP2, colour=diagnosis)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + 
  xlab("UMAP1") + ylab("UMAP2") + ggtitle("") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +
  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

jpeg("Bins_Relab_UMAP1.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(Bins_Relab_ord_meta$UMAP1 ~ Bins_Relab_ord_meta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="UMAP1", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
#stripchart(beta_Disper$distances ~ PTR_meta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
kruskalmc(Bins_Relab_ord_meta$UMAP1 ~ Bins_Relab_ord_meta$diagnosis)

aggregate(Bins_Relab_ord_meta$UMAP1, list(Bins_Relab_ord_meta$diagnosis), FUN=median)



#------------------ Tried but Not worked
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("M3C")
library(M3C)
umap(PTR_sample_filtered_filtered)
umap(PTR_meta[,1:188],labels=as.factor(PTR_meta$diagnosis),controlscale=TRUE,scale=3, colvec=c("skyblue", "gold", "blue"))


################ Not Used
library(umap)
RefPTR.data = PTR_meta[,1:188]
RefPTR.labels = PTR_meta[, "diagnosis"]

RefPTR.umap = umap(RefPTR.data)

RefPTR.umap
head(RefPTR.umap$layout, 3)

plot.function(RefPTR.umap, RefPTR.labels)
