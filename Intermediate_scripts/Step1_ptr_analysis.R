##################################################################################################################
# ----- This is the Prelimany distribution Check of Bins - in total, diagnosis specific and Patient specific ----
#################################################################################################################

#########################
#@ Data Filtering #######
#########################

setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")
#-- Load PTR values for all samples
PTR <- read.csv(file="PTR_output.csv", row.names=1, header=T)
PTR[is.na(PTR)] <- 0 # Replace NAs with Zero

PTR.sums <- colSums(PTR)
PTR_sample_filtered <- PTR[ , which(PTR.sums > 0)]
#-In 23 samples no bacterial growth rates were observed and hence Removed
#CSM79HPU, CSM7KOOR, HSM67VDZ,HSM6XRR9, HSM6XRS8, HSM7J4KA, HSM7J4LH, HSM7J4LP, HSM7J4QJ, HSMA33IO, HSMA33MX, HSMA33N4, MSM5LLDU, MSM6J2HB, MSM9VZLZ, PSM6XBSU.TR, PSM6XBTD, PSM7J129, PSM7J16H, PSM7J177, PSM7J1AW, PSM7J1BH, PSMA266C
PTR_sample_filtered <- data.frame(t(PTR_sample_filtered))
#-- Remove Bins - which is not present in at least 5% of the samples (66 samples out of 1315)
PTR_sample_filtered_filtered <- dropspc(PTR_sample_filtered, 66)
#write.table (PTR_sample_filtered_filtered, file="PTR_output_filetered.txt", sep="\t")

#--- Load Metadata - In this metafile 300 samples are extra which were not analyzed 
meta <- read.csv(file="hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)
pid <- meta$Participant.ID
visit <- meta$visit_num
age <- meta$consent_age
sex <- meta$sex
site <- meta$site_name
#--- Dietary data
alchohol <- meta$Alcohol..beer..brandy..spirits..hard.liquor..wine..aperitif..etc..
vegetables <- meta$Vegetables..salad..tomatoes..onions..greens..carrots..peppers..green.beans..etc.
probiotics <- meta$Probiotic
diarymilk <- meta$Dairy..milk..cream..ice.cream..cheese..cream.cheese.
yogurt <- meta$Yogurt.or.other.foods.containing.active.bacterial.cultures..kefir..sauerkraut.
fruits <- meta$Fruits..no.juice...Apples..raisins..bananas..oranges..strawberries..blueberries
#--This one was only measured for UC patients and hence won't be included in complete comparisons
diagnosis <- meta$diagnosis
sccai <- meta$sccai
fecalcal <- meta$fecalcal
meta_new <- cbind(pid, visit, age, sex, site, alchohol, vegetables, probiotics, diarymilk, yogurt, fruits, diagnosis, sccai, fecalcal)
meta_new <- data.frame(meta_new)
rownames(meta_new) <- row.names(meta)

PTR_meta <- merge(PTR_sample_filtered_filtered, meta_new, by=0, all=F)
rownames(PTR_meta) <- PTR_meta$Row.names; PTR_meta$Row.names <- NULL

#---- No Taxa Filtering. Specifically will be Used for Patient Specific Analysis
PTR_meta_no_filtered <- merge(PTR_sample_filtered, meta_new, by=0, all=F)
rownames(PTR_meta_no_filtered) <- PTR_meta_no_filtered$Row.names; PTR_meta_no_filtered$Row.names <- NULL

#--Save this data to be used in Other Scripts- If need to add any metadata add save it again
#save.image("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.RData")
#write.table (PTR_meta, file="PTR_meta.txt", sep="\t")

###########################################################
##### Part1: Total Analysis################################
###########################################################
#-- This analysis is on total 70 Bins across all samples
library (vegan)
library (ape)
library (variancePartition)

setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")
load ("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.RData")

#----- Microbial growth rates are varying in response to which Variable- SampleSpecific, age, sex, visit, diagnosis 
PTR_dist <-vegdist(PTR_meta[,1:70], "bray")
PTR_dist_pcoa <- pcoa(PTR_dist)

adonis2(PTR_dist ~ PTR_meta$diagnosis)
#PTR_meta$diagnosis    2    10.28 0.01853 12.384  0.001 ***
adonis2(PTR_dist ~ PTR_meta$sex)
#PTR_meta$sex    1     4.90 0.00883 11.699  0.001 ***
adonis2(PTR_dist ~ PTR_meta$visit)
#PTR_meta$visit   23     5.95 0.01072 0.608      1
adonis2(PTR_dist ~ PTR_meta$pid)
#PTR_meta$pid  105   418.44 0.75385 35.263  0.001 ***

#------ Variance Particiting Analysis
library (variancePartition)
PTR_filtered <- data.frame(t(PTR_meta[,1:70]))
meta_filtered <- PTR_meta[,71:82] #SCCAI was not included

model <- ~ (1|diagnosis)+(1|age)+(1|sex)+(1|site)+(1|alchohol)+(1|vegetables)+(1|probiotics)+(1|diarymilk)+(1|yogurt)+(1|fruits)+(1|visit)+(1|pid)
var_part_model <- fitExtractVarPartModel(PTR_filtered, model,meta_filtered)
save(var_part_model, file="Total_analysis/var_part_model.Rdata")

dim(var_part_model)
write.table(var_part_model, file="Total_analysis/varPart_model.txt", sep = "\t")
#sort the terms of the model by average variance explained across all genes, so when we plot they will be sorted by overall importance:
vp2 <- sortCols(var_part_model)
# Violin plot
jpeg("Total_analysis/variance_partision.jpg", height = 7, width = 7, units = 'in', res = 600)
plotVarPart( vp2  ,label.angle = 90)
dev.off ()
#sort genes based on variance explained by Condition (#Here it is a treatment)
head(var_part_model[order(var_part_model$diagnosis, decreasing=TRUE),])

#---- Check the Growth Rate changes across patients, during diffent visits in CD, UC and nonIBD groups
visit_order <- c('4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '16', '18', '19', '20', '21', '22', '23', '25', '26', '27', '28', '29', '30')
jpeg("Total_analysis/Total_diagnosis_visit_Bin1049.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = PTR_meta, aes(factor(visit, level = visit_order), y = Bin1049, color = diagnosis, group=diagnosis)) +       
  #ggplot(data = CD_filtered_meta, aes(x = visit, y = Bin1049, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1049")
dev.off ()
jpeg("Total_analysis/Total_diagnosis_PID_Bin1049.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = PTR_meta, aes(x=pid, y = Bin1049, color = diagnosis, group=diagnosis)) +       
  #ggplot(data = CD_filtered_meta, aes(x = visit, y = Bin1049, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1049")
dev.off ()

#########################################################################
##### Part2: Diagnosis Specific Analysis ################################
#########################################################################
setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")
load ("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.RData")

#---- Diagnosis Specific Filtering
################################################################################
#------------------------------- CD ----------------------------------------
################################################################################
CD <- subset(PTR_meta, diagnosis=="CD")
#-- Remove Bin which is not present in Even 25% of the samples
CD_filtered <- dropspc(CD[,1:70], 148)
#-- Now to remove samples in Which non of these 11 Bins is present
CD_filtered <- data.frame(t(CD_filtered))
CD.sums <- colSums(CD_filtered)
CD_sample_filtered <- CD_filtered[ , which(CD.sums > 0)]
#-In 69 Samples None of these 11 Bins present - So removed
CD_sample_filtered <- data.frame(t(CD_sample_filtered))
#--- Finally for "522" samples Now add Metadata
CD_filtered_meta <- merge(CD_sample_filtered, meta_new, by=0, all=F)
rownames(CD_filtered_meta) <- CD_filtered_meta$Row.names; CD_filtered_meta$Row.names <- NULL
#write.table(CD_filtered_meta, file="DiagnosisSpecific_output/CD_filtered_PTR.txt", sep="\t")

#-------- Plot Each Bin to check their Growth rate disribution across Patients
visit_order <- c('4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '16', '18', '19', '20', '21', '22', '23', '25', '26', '27', '28', '29', '30')
jpeg("DiagnosisSpecific_output/CD_Bin1049.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1049, color = pid, group=pid)) +       
  #ggplot(data = CD_filtered_meta, aes(x = visit, y = Bin1049, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1049")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin12.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin12, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin12")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin168.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin168, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin168")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin1908.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1908, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1908")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin1970.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1970, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1970")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin2338.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin2338, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin2338")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin2480.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin2480, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin2480")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin3362.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin3362, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin3362")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin3733.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin3733, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin3733")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin3746.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin3746, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin3746")
dev.off ()

jpeg("DiagnosisSpecific_output/CD_Bin387.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = CD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin387, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin387")
dev.off ()


#-------- Plot variance of each Bin Across patients with Chron's Disease
#-- Calculate Varyance - Max Variance is 0.25
CD_bins_details <- sapply(CD_filtered_meta[,1:11], function(x) c(sum=sum(x), mean=mean(x), var=var(x), sd=sd(x)))
CD_bins_details <- data.frame(t(CD_bins_details))
jpeg("DiagnosisSpecific_output/CD_var.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (CD_bins_details$var, xlab = "Bins", ylab="Variance (log2 (PTR)- CD patients", pch=19, cex=1.5, ylim = c(0,0.25),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(CD_bins_details$var, labels=rownames(CD_bins_details), cex=0.7, font=1, pos=3)
dev.off ()

#Calcualte Index of dispersion: Like variance-to-mean ratio (VMR)
CD_bins_details_dispersionIndex <- CD_bins_details %>% mutate(Index_of_dis = var/mean)
jpeg("DiagnosisSpecific_output/CD_Index_of_Dispersion.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (CD_bins_details_dispersionIndex$Index_of_dis, xlab = "Bins", ylab="Index of Dispersion (log2 (PTR)- CD patients", pch=19, cex=1.5, ylim = c(0.2,0.6),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(CD_bins_details_dispersionIndex$Index_of_dis, labels=rownames(CD_bins_details_dispersionIndex), cex=0.7, font=1, pos=3)
dev.off ()

# Coefficient of variation (CV), also known as relative standard deviation (RSD), is a standardized measure of dispersion of a probability distribution or frequency distribution
CD_bins_details_dispersionIndexco_of_variation <- CD_bins_details_dispersionIndex %>% mutate(Covariation = sd/mean)
jpeg("DiagnosisSpecific_output/CD_coefficient_of_Variation.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (CD_bins_details_dispersionIndexco_of_variation$Covariation, xlab = "Bins", ylab="Coefficient of variation (log2 (PTR)- CD patients", pch=19, cex=1.5, ylim = c(0.5,2),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(CD_bins_details_dispersionIndexco_of_variation$Covariation, labels=rownames(CD_bins_details_dispersionIndexco_of_variation), cex=0.7, font=1, pos=3)
dev.off ()

#-Alternative Way is to Create a HeatMap (which will Show -scores of each Bin)
jpeg("DiagnosisSpecific_output/CD_heatmap_Zscores.jpg", height = 4, width = 5, units = 'in', res = 600)
aheatmap(CD_filtered_meta[,1:11], scale = "row", color = "-RdYlBu2:100", Colv = FALSE)
dev.off ()

#-------- Plot Z-scores of each Bin Across patients with Chron's Disease
#Bin1049_zscores <- CD_filtered_meta %>% mutate(Bin1049_zscore = (Bin1049 - mean(Bin1049))/sd(Bin1049))
#Bin12_zscores <- Bin1049_zscores %>% mutate(Bin12_zscore = (Bin12 - mean(Bin12))/sd(Bin12))
#Bin168_zscores <- Bin12_zscores %>% mutate(Bin168_zscore = (Bin168 - mean(Bin168))/sd(Bin168))
#Bin1908_zscores <- Bin168_zscores %>% mutate(Bin1908_zscore = (Bin1908 - mean(Bin1908))/sd(Bin1908))
#Bin1970_zscores <- Bin1908_zscores %>% mutate(Bin1970_zscore = (Bin1970 - mean(Bin1970))/sd(Bin1970))
#Bin2338_zscores <- Bin1970_zscores %>% mutate(Bin2338_zscore = (Bin2338 - mean(Bin2338))/sd(Bin2338))
#Bin2480_zscores <- Bin2338_zscores %>% mutate(Bin2480_zscore = (Bin2480 - mean(Bin2480))/sd(Bin2480))
#Bin3362_zscore <- Bin2480_zscores %>% mutate(Bin3362_zscore = (Bin3362 - mean(Bin3362))/sd(Bin3362))
#Bin3733_zscore <- Bin3362_zscore %>% mutate(Bin3733_zscore = (Bin3733 - mean(Bin3733))/sd(Bin3733))
#Bin3746_zscore <- Bin3733_zscore %>% mutate(Bin3746_zscore = (Bin3746 - mean(Bin3746))/sd(Bin3746))
#CD_All_bins_zscores <- Bin3746_zscore %>% mutate(Bin387_zscore = (Bin387 - mean(Bin387))/sd(Bin387))

#CD_zscores_summary <- sapply(CD_All_bins_zscores[,26:36], function(x) c(sum=sum(x), mean=mean(x), var=var(x), sd=sd(x)))
#CD_zscores_summary <- data.frame(t(CD_zscores_summary))
#plot (CD_zscores_summary$mean, xlab = "Bins", ylab="mean Z-scores (log2 (PTR)- CD patients", pch=19, cex=1.5,
   #   col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
#text(CD_zscores_summary$mean, labels=rownames(CD_bins_details), cex=0.7, font=1, pos=3)

################################################################################
#------------------------------- UC ----------------------------------------
################################################################################
UC <- subset(PTR_meta, diagnosis=="UC")
#-- Remove Bin which is not present in Even 25% of the samples
UC_filtered <- dropspc(UC[,1:70], 91)
#-- Now to remove samples in Which non of these 11 Bins is present
UC_filtered <- data.frame(t(UC_filtered))
UC.sums <- colSums(UC_filtered)
UC_sample_filtered <- UC_filtered[ , which(UC.sums > 0)]
#-In 5 Samples None of these 11 Bins present - So removed
UC_sample_filtered <- data.frame(t(UC_sample_filtered))
#--- Finally for "357" samples Now add Metadata
UC_filtered_meta <- merge(UC_sample_filtered, meta_new, by=0, all=F)
rownames(UC_filtered_meta) <- UC_filtered_meta$Row.names; UC_filtered_meta$Row.names <- NULL
#write.table(UC_filtered_meta, file="DiagnosisSpecific_output/UC_filtered_PTR.txt", sep="\t")

#-------- Plot Each Bin to check their Growth rate disribution across Patients
visit_order <- c('4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '16', '18', '19', '20', '21', '22', '23', '25', '26', '27', '28', '29', '30')
jpeg("DiagnosisSpecific_output/UC_Bin1049.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1049, color = pid, group=pid)) +       
  #ggplot(data = UC_filtered_meta, aes(x = visit, y = Bin1049, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1049")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin12.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin12, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin12")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin168.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin168, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin168")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin1908.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1908, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1908")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin1970.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1970, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1970")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin2338.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin2338, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin2338")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin2480.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin2480, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin2480")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin3362.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin3362, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin3362")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin3733.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin3733, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin3733")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin3746.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin3746, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin3746")
dev.off ()

jpeg("DiagnosisSpecific_output/UC_Bin387.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = UC_filtered_meta, aes(factor(visit, level = visit_order), y = Bin387, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin387")
dev.off ()

#-------- Plot variance of each Bin Across Ulcerative Colitis patients 
#-- Calculate Varyance - Max Variance is 0.25
UC_bins_details <- sapply(UC_filtered_meta[,1:11], function(x) c(sum=sum(x), mean=mean(x), var=var(x), sd=sd(x)))
UC_bins_details <- data.frame(t(UC_bins_details))
jpeg("DiagnosisSpecific_output/UC_var.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (UC_bins_details$var, xlab = "Bins", ylab="Variance (log2 (PTR)- UC patients", pch=19, cex=1.5, ylim = c(0,0.25),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(UC_bins_details$var, labels=rownames(UC_bins_details), cex=0.7, font=1, pos=3)
dev.off ()

#Calcualte Index of dispersion: Like variance-to-mean ratio (VMR)
UC_bins_details_dispersionIndex <- UC_bins_details %>% mutate(Index_of_dis = var/mean)
jpeg("DiagnosisSpecific_output/UC_Index_of_Dispersion.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (UC_bins_details_dispersionIndex$Index_of_dis, xlab = "Bins", ylab="Index of Dispersion (log2 (PTR)- UC patients", pch=19, cex=1.5, ylim = c(0.2,0.6),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(UC_bins_details_dispersionIndex$Index_of_dis, labels=rownames(UC_bins_details_dispersionIndex), cex=0.7, font=1, pos=3)
dev.off ()

# Coefficient of variation (CV), also known as relative standard deviation (RSD), is a standardized measure of dispersion of a probability distribution or frequency distribution
UC_bins_details_dispersionIndexco_of_variation <- UC_bins_details_dispersionIndex %>% mutate(Covariation = sd/mean)
jpeg("DiagnosisSpecific_output/UC_coefficient_of_Variation.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (UC_bins_details_dispersionIndexco_of_variation$Covariation, xlab = "Bins", ylab="Coefficient of variation (log2 (PTR)- CD patients", pch=19, cex=1.5, ylim = c(0.5,2),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(UC_bins_details_dispersionIndexco_of_variation$Covariation, labels=rownames(UC_bins_details_dispersionIndexco_of_variation), cex=0.7, font=1, pos=3)
dev.off ()

#-Alternative Way is to Create a HeatMap (which will Show -scores of each Bin)
jpeg("DiagnosisSpecific_output/UC_heatmap_Zscores.jpg", height = 4, width = 5, units = 'in', res = 600)
aheatmap(UC_filtered_meta[,1:11], scale = "row", color = "-RdYlBu2:100", Colv = FALSE)
dev.off ()


################################################################################
#------------------------------- nonIBD ----------------------------------------
################################################################################
nonIBD <- subset(PTR_meta, diagnosis=="nonIBD")
#-- Remove Bin which is not present in Even 25% of the samples
nonIBD_filtered <- dropspc(nonIBD[,1:70], 91)
#-- Now to remove samples in Which non of these 12 Bins is present
nonIBD_filtered <- data.frame(t(nonIBD_filtered))
nonIBD.sums <- colSums(nonIBD_filtered)
nonIBD_sample_filtered <- nonIBD_filtered[ , which(nonIBD.sums > 0)]
#-In 13 Samples None of these 12 Bins present - So removed
nonIBD_sample_filtered <- data.frame(t(nonIBD_sample_filtered))
#--- Finally for "349" samples Now add Metadata
nonIBD_filtered_meta <- merge(nonIBD_sample_filtered, meta_new, by=0, all=F)
rownames(nonIBD_filtered_meta) <- nonIBD_filtered_meta$Row.names; nonIBD_filtered_meta$Row.names <- NULL
#write.table(nonIBD_filtered_meta, file="DiagnosisSpecific_output/nonIBD_filtered_PTR.txt", sep="\t")

#-------- Plot Each Bin to check their Growth rate disribution across Patients
visit_order <- c('4', '5', '6', '7', '8', '9', '11', '12', '13', '14', '15', '16', '18', '19', '20', '21', '22', '23', '25', '26', '27', '28', '29', '30')
jpeg("DiagnosisSpecific_output/nonIBD_Bin1049.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1049, color = pid, group=pid)) +       
  #ggplot(data = nonIBD_filtered_meta, aes(x = visit, y = Bin1049, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1049")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin12.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin12, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin12")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin168.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin168, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin168")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin1811.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1811, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1811")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin1908.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1908, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1908")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin1970.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin1970, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin1970")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin2338.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin2338, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin2338")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin2480.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin2480, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin2480")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin3362.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin3362, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin3362")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin3733.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin3733, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin3733")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin3746.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin3746, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin3746")
dev.off ()

jpeg("DiagnosisSpecific_output/nonIBD_Bin387.jpg", height = 5, width = 7, units = 'in', res = 600)
ggplot(data = nonIBD_filtered_meta, aes(factor(visit, level = visit_order), y = Bin387, color = pid, group=pid)) +       
  geom_line() + geom_point(size = 2, shape = 21) + theme_bw() + 
  ylab(paste("log2 (PTR)")) + ggtitle("Bin387")
dev.off ()

#-------- Plot variance of each Bin Across nonIBD patients
#-- Calculate Varyance - Max Variance is 0.25
nonIBD_bins_details <- sapply(nonIBD_filtered_meta[,1:12], function(x) c(sum=sum(x), mean=mean(x), var=var(x), sd=sd(x)))
nonIBD_bins_details <- data.frame(t(nonIBD_bins_details))
jpeg("DiagnosisSpecific_output/nonIBD_var.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (nonIBD_bins_details$var, xlab = "Bins", ylab="Variance (log2 (PTR)- nonIBD subjects", pch=19, cex=1.5, ylim = c(0,0.25),
      col = c("red", "blue", "orange", "black", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(nonIBD_bins_details$var, labels=rownames(nonIBD_bins_details), cex=0.7, font=1, pos=3)
dev.off ()

#Calcualte Index of dispersion: Like variance-to-mean ratio (VMR)
nonIBD_bins_details_dispersionIndex <- nonIBD_bins_details %>% mutate(Index_of_dis = var/mean)
jpeg("DiagnosisSpecific_output/nonIBD_Index_of_Dispersion.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (nonIBD_bins_details_dispersionIndex$Index_of_dis, xlab = "Bins", ylab="Index of Dispersion (log2 (PTR)- nonIBD", pch=19, cex=1.5, ylim = c(0.15,0.6),
      col = c("red", "blue", "orange", "black", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(nonIBD_bins_details_dispersionIndex$Index_of_dis, labels=rownames(nonIBD_bins_details_dispersionIndex), cex=0.7, font=1, pos=3)
dev.off ()

# Coefficient of variation (CV), also known as relative standard deviation (RSD), is a standardized measure of dispersion of a probability distribution or frequency distribution
nonIBD_bins_details_dispersionIndexIndexco_of_variation <- nonIBD_bins_details_dispersionIndex %>% mutate(Covariation = sd/mean)
jpeg("DiagnosisSpecific_output/nonIBD_coefficient_of_Variation.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (nonIBD_bins_details_dispersionIndexIndexco_of_variation$Covariation, xlab = "Bins", ylab="Coefficient of variation (log2 (PTR)- nonIBD", pch=19, cex=1.5, ylim = c(0.5,2),
      col = c("red", "blue", "orange", "black", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(nonIBD_bins_details_dispersionIndexIndexco_of_variation$Covariation, labels=rownames(nonIBD_bins_details_dispersionIndexIndexco_of_variation), cex=0.7, font=1, pos=3)
dev.off ()

#-Alternative Way is to Create a HeatMap (which will Show -scores of each Bin)
jpeg("DiagnosisSpecific_output/nonIBD_heatmap_Zscores.jpg", height = 4, width = 5, units = 'in', res = 600)
aheatmap(nonIBD_filtered_meta[,1:12], scale = "row", color = "-RdYlBu2:100", Colv = FALSE)
dev.off ()

#########################################################################
##### Part3: Patient Specific Analysis ################################
#########################################################################
setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")
load ("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.RData")
#---- Patient Specific Filtering - For this No Filtered data will be Used.
C3001 <- subset(PTR_meta_no_filtered, pid=="C3001")
#-- Remove Bin which is not present in Even 50% of the samples
C3001_filtered <- dropspc(C3001[,1:1503], 4)
#write.table(C3001_filtered, "PatientSpceific_output/C3001_PTR_filtered.txt", sep = "\t")
C3001_filtered_meta <- cbind (C3001_filtered, C3001[,1504:1517])


jpeg("PatientSpceific_output/C3001.jpg", height = 5, width = 7, units = 'in', res = 600)
# Create a first line
plot(C3001_filtered_meta$visit, C3001_filtered_meta$Bin1049, type = "b", frame = FALSE, pch = 19, 
     col = "red", xlab = "visits", ylab = "log2(PTR)", ylim = c(0, 3))
# Add other lines
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin12, pch = 19, col = "blue", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin18, pch = 19, col = "orange", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin2018, pch = 19, col = "darkred", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin21, pch = 19, col = "cyan", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin23, pch = 19, col = "darkgray", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin2362, pch = 19, col = "Magenta", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin24, pch = 19, col = "Maroon", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin2402, pch = 19, col = "limegreen", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin25, pch = 19, col = "darkgreen", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin3897, pch = 19, col = "purple", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin659, pch = 19, col = "navy", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin687, pch = 19, col = "darkgoldenrod4", type = "b", lty = 1)
lines(C3001_filtered_meta$visit, C3001_filtered_meta$Bin7, pch = 19, col = "black", type = "b", lty = 1)
# Add a legend to the plot
legend("topright", legend=c("Bin1049", "Bin12", "Bin18", "Bin2018", "Bin21", "Bin23", "Bin2362", "Bin24", "Bin2402", "Bin25", "Bin3897", "Bin659", "Bin687", "Bin7"),
       col=c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"), 
       lty = 1, cex=0.8, ncol = 5:5, pch=19)
dev.off ()

#-- Calculate Varyance of Each Variable
C3001_bins_details <- sapply(C3001_filtered_meta[,1:14], function(x) c(sum=sum(x), mean=mean(x), var=var(x), sd=sd(x)))
C3001_bins_details <- data.frame(t(C3001_bins_details))
jpeg("PatientSpceific_output/C3001_var.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (C3001_bins_details$var, xlab = "Bins", ylab="Variance (log2 (PTR) among different visits", pch=19, cex=1.5, ylim = c(0,1.1),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"))
text(C3001_bins_details$var, labels=rownames(C3001_bins_details), cex=0.7, font=1, pos=3)
dev.off ()

#Calcualte Index of dispersion: Like variance-to-mean ratio (VMR)
C3001_bins_details_dispersionIndex <- C3001_bins_details %>% mutate(Index_of_dis = var/mean)
jpeg("PatientSpceific_output/C3001_Index_of_Dispersion.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (C3001_bins_details_dispersionIndex$Index_of_dis, xlab = "Bins", ylab="Index of Dispersion (log2 (PTR)- C3001", pch=19, cex=1.5, ylim = c(0.0005,0.75),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"))
text(C3001_bins_details_dispersionIndex$Index_of_dis, labels=rownames(C3001_bins_details_dispersionIndex), cex=0.7, font=1, pos=3)
dev.off ()

# Coefficient of variation (CV), also known as relative standard deviation (RSD), is a standardized measure of dispersion of a probability distribution or frequency distribution
C3001_bins_details_dispersionIndexco_of_variation <- C3001_bins_details_dispersionIndex %>% mutate(Covariation = sd/mean)
jpeg("PatientSpceific_output/C3001_coefficient_of_Variation.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (C3001_bins_details_dispersionIndexco_of_variation$Covariation, xlab = "Bins", ylab="Coefficient of variation (log2 (PTR)- C3001", pch=19, cex=1.5, ylim = c(0.03,0.75),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"))
text(C3001_bins_details_dispersionIndexco_of_variation$Covariation, labels=rownames(C3001_bins_details_dispersionIndexco_of_variation), cex=0.7, font=1, pos=3)
dev.off ()

#-Alternative Way is to Create a HeatMap (which will Show -scores of each Bin)
jpeg("PatientSpceific_output/C3001_heatmap_Zscores.jpg", height = 2, width = 5, units = 'in', res = 600)
aheatmap(C3001_filtered_meta[,1:14], scale = "row", color = "-RdYlBu2:100", Colv = FALSE)
dev.off ()


###################################################################
#------- Select Three Patients one from Each UC, CD and nonIBD

C3005 <- subset(PTR_meta_no_filtered, pid=="C3005") #For nonIBD
C3005_filtered <- dropspc(C3005[,1:1503], 4)
C3010 <- subset(PTR_meta_no_filtered, pid=="C3010") #For UC
C3010_filtered <- dropspc(C3010[,1:1503], 6)
M2041 <- subset(PTR_meta_no_filtered, pid=="M2041") #For CD
M2041_filtered <- dropspc(M2041[,1:1503], 7)

##---- Merge Three Patients
library(plyr)
df_list <- list(C3005_filtered, C3010_filtered, M2041_filtered)
for(i in 1:length(df_list)){
  colnames(df_list[[i]]) <- paste0( names(df_list)[i], "_", colnames(df_list[[i]]) )
  df_list[[i]]$ROWNAMES  <- rownames(df_list[[i]])
}
Merged_Df <- join_all( df_list, by="ROWNAMES", type="full" )
rownames(Merged_Df) <- Merged_Df$ROWNAMES; Merged_Df$ROWNAMES <- NULL

Merged_Df[is.na(Merged_Df)] <- 0 # Replace NAs with Zero


##-- Remove Bins Which are not consistently present in All three Patients.
ThreePID <- rbind(H4039, M2064, M2025)

#-- Remove Bin which is not present in Even 50% of the samples
ThreePID_filtered <- dropspc(ThreePID[,1:1503], 16)
ThreePID_filtered_meta <- cbind (ThreePID_filtered, ThreePID[,1514:1517])
