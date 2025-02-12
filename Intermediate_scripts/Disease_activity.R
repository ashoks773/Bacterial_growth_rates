#---- Calculate Disease Activity Index on All sets - Bins, RefGenomes - PTRs and Relative abundnaces, and Metaphln relab also
setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/Disease_Activity")

source("Scripts/metaphln_AI.R") # Calculated on the relative abundances Calculated by Authors using Metaphln ** Downloaded from iHMP website **
source("Scripts/binsPTR_AI.R")  #Calculated on PTR estimates of 70 Bins (present in atleast 5% of the samples)
source("Scripts/binsRelab_AI.R") #Calculated on Relative abundance estimates of Above 70 Bins 
source("Scripts/RefPTR_AI.R") #Caclcualted on PTR estimates of 188 Reference Genomes (present in atleast 5% of the samples)
source("Scripts/RefRelab_AI.R") #Caclcualted on PTR estimates of 188 Reference Genomes (present in atleast 5% of the samples)

#------ Activity Indexes
#1. metaphln_activity_index: Calculated on the relative abundances Calculated by Authors using Metaphln
#2. binsPTR_activity_index: Calculated on PTR estimates of 70 Bins (present in atleast 5% of the samples)
#3. binsRelab_activity_index: Calculated on Relative abundance estimates of Above 70 Bins
#4. Ref_PTR_activity_index: Caclcualted on PTR estimates of 188 Reference Genomes (present in atleast 5% of the samples)
#5. Ref_Relab_activity_index: Caclcualted on Relative abundance estimates of above 188 Reference Genomes

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
write.table(AllSets_Activity_Index, file="~/Box/David_Casero_Lab/HMP_IBD_Project/Disease_Activity/AllSets_Activity_Index.txt", sep = "\t")
save.image("~/Box/David_Casero_Lab/HMP_IBD_Project/Disease_Activity/disease_activity.RData")

#---- Disease Activity and Disease Diagnosis
load ("~/Box/David_Casero_Lab/HMP_IBD_Project/Disease_Activity/disease_activity.RData")

#--meta
meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)
#--- If Needed to fix the meta
source('~/Box/David_Casero_Lab/HMP_IBD_Project/Common_scripts/fix_metaata.R')
meta_new <- fix_metadata(meta) #-- Check If required otherwise- Its fine to use the Original Metadata Only

#Merge metadata with Disease activity Index
AllSets_Activity_Indexmeta <- merge(AllSets_Activity_Index, meta, by=0, all=F)
rownames(AllSets_Activity_Indexmeta) <- AllSets_Activity_Indexmeta$Row.names; AllSets_Activity_Indexmeta$Row.names <- NULL


jpeg("Figures/metaphln_AI.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(AllSets_Activity_Indexmeta$metaphln ~ AllSets_Activity_Indexmeta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="Disease Activity (MetaPhlAn)", 
        cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(AllSets_Activity_Indexmeta$metaphln ~ AllSets_Activity_Indexmeta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'gray78')
dev.off ()
kruskalmc(AllSets_Activity_Indexmeta$metaphln ~ AllSets_Activity_Indexmeta$diagnosis)


jpeg("Figures/binPTR_AI.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(AllSets_Activity_Indexmeta$binsPTR ~ AllSets_Activity_Indexmeta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="Disease Activity (bins PTRs)", 
        cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(AllSets_Activity_Indexmeta$binsPTR ~ AllSets_Activity_Indexmeta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'gray78')
dev.off ()
kruskalmc(AllSets_Activity_Indexmeta$binsPTR ~ AllSets_Activity_Indexmeta$diagnosis)
#Comparisons
#obs.dif critical.dif difference
#CD-nonIBD 197.89920     60.45665       TRUE
#CD-UC      14.98562     60.50859      FALSE
#nonIBD-UC 182.91358     67.30298       TRUE

jpeg("Figures/binRelab_AI.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(AllSets_Activity_Indexmeta$binsRelab ~ AllSets_Activity_Indexmeta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="Disease Activity (bins Relative Abundance)", 
        cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(AllSets_Activity_Indexmeta$binsRelab ~ AllSets_Activity_Indexmeta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'gray78')
dev.off ()
kruskalmc(AllSets_Activity_Indexmeta$binsRelab ~ AllSets_Activity_Indexmeta$diagnosis)
#obs.dif critical.dif difference
#CD-nonIBD 187.51873     60.45665       TRUE
#CD-UC      23.17384     60.50859      FALSE
#nonIBD-UC 164.34488     67.30298       TRUE

jpeg("Figures/Ref_PTR_AI.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(AllSets_Activity_Indexmeta$Ref_PTR ~ AllSets_Activity_Indexmeta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="Disease Activity (Ref Genome PTRs)", 
        cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(AllSets_Activity_Indexmeta$Ref_PTR ~ AllSets_Activity_Indexmeta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'gray78')
dev.off ()
kruskalmc(AllSets_Activity_Indexmeta$Ref_PTR ~ AllSets_Activity_Indexmeta$diagnosis)
#obs.dif critical.dif difference
#CD-nonIBD 227.59738     60.45665       TRUE
#CD-UC      30.64264     60.50859      FALSE
#nonIBD-UC 196.95474     67.30298       TRUE

jpeg("Figures/Ref_Relab_AI.jpg", height = 5, width = 3, units = 'in', res = 600)
boxplot(AllSets_Activity_Indexmeta$Ref_Relab ~ AllSets_Activity_Indexmeta$diagnosis, col=c("darkred", "steelblue2", "orange"), ylab="Disease Activity (Ref Genome Relative Abundance)", 
        cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(AllSets_Activity_Indexmeta$Ref_Relab ~ AllSets_Activity_Indexmeta$diagnosis, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'gray78')
dev.off ()
kruskalmc(AllSets_Activity_Indexmeta$Ref_Relab ~ AllSets_Activity_Indexmeta$diagnosis)
#obs.dif critical.dif difference
#CD-nonIBD 137.48724     60.45665       TRUE
#CD-UC      25.23325     60.50859      FALSE
#nonIBD-UC 112.25399     67.30298       TRUE


##########################
#---- Density plots -----
#########################

Colors=c("darkred", "steelblue2", "orange")
library(ggplot2)
library(egg)

jpeg("Figures/metaphln_Density.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=AllSets_Activity_Indexmeta, aes(x=metaphln, color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  geom_vline(xintercept=disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Dysbiosis score (MetaPhlAn Relabun)") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

jpeg("Figures/binPTR_Density.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=AllSets_Activity_Indexmeta, aes(x=binsPTR, color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  geom_vline(xintercept=binsPTR_disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Dysbiosis score (MAG PTR)") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

jpeg("Figures/binRelab_Density.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=AllSets_Activity_Indexmeta, aes(x=binsRelab, color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  geom_vline(xintercept=BinsRelab_disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Dysbiosis score (MAG Relabun)") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

jpeg("Figures/Ref_PTR_Density.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=AllSets_Activity_Indexmeta, aes(x=Ref_PTR, color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  geom_vline(xintercept=Ref_PTR_disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Dysbiosis score (Ref.Genome PTR)") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

jpeg("Figures/Ref_Relab_Density.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data=AllSets_Activity_Indexmeta, aes(x=Ref_Relab, color=diagnosis, fill=diagnosis)) +
  geom_density(alpha=0.3, size=0.25) +
  geom_vline(xintercept=Ref_Relab_disease_activity_threshold, size=0.5) +
  scale_color_manual(values=Colors) +
  scale_fill_manual(values=Colors) +
  xlab("Dysbiosis score (Ref.Genome Relabun)") + ylab("Density") +
  guides(alpha="none", fill="none", color="none")
dev.off ()

##################
#---Find variable Importance
#################
set.seed(7)
diagnosis <- data.frame (AllSets_Activity_Indexmeta$diagnosis)
colnames(diagnosis) <- "diagnosis"
Activity_Index <- AllSets_Activity_Indexmeta[,1:5] #Remove metaphlan
Activity_Index_diagnosis <- cbind (Activity_Index, diagnosis)

#-- Random Forest
RF1 <- randomForest(factor(diagnosis)~., data=Activity_Index_diagnosis, ntree = 500, importance = TRUE, do.trace = 100, cv.fold = 10, na.action = na.omit)
jpeg("Figures/RF_varImp.jpg", height = 5, width = 4, units = 'in', res = 600)
varImpPlot(RF1, type=1)
dev.off ()
imp <- importance(RF1, type=1)
write.table(imp, file="Results/rf_importance10.txt", sep="\t")

setwd("Results/")
ldf <- list() # creates a list
listcsv <- dir(pattern = "*.txt") # creates the list of all the csv files in the directory
for (k in 1:length(listcsv)){
  ldf[[k]] <- read.csv(listcsv[k], sep="\t", row.names = 1, header = T)
}
str(ldf[[1]]) 
importance_10times <- cbind (ldf[[1]], ldf[[2]], ldf[[3]], ldf[[4]], ldf[[5]], ldf[[6]], ldf[[7]], ldf[[8]], ldf[[9]], ldf[[10]])
importance_10times <- data.frame(t(importance_10times))

library(earth)
regressor <- earth(factor(diagnosis)~., data=Activity_Index_diagnosis) # build model
ev <- evimp (regressor) # estimate variable importance
plot (ev)

library(caret)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
regressor<- train(factor(diagnosis)~., data=Activity_Index_diagnosis, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(regressor, scale=FALSE)

library(randomForest)
library(DALEX)
library("breakDown")

regressor <- randomForest(factor(diagnosis)~., data=Activity_Index_diagnosis, importance=TRUE) # fit the random forest with default parameter
# Variable importance with DALEX
explained_rf <- explain(regressor, data=Activity_Index_diagnosis, y=Activity_Index_diagnosis$diagnosis)
# Get the variable importances
varimps = variable_dropout(explained_rf, type='raw')
print(varimps)
plot(varimps)


fit <- lm(metaphln~diagnosis, data=Activity_Index_diagnosis)


