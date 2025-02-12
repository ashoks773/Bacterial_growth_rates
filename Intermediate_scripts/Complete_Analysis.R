library (vegan)
library (ape)
library (variancePartition)

setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")
load ("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.RData")

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



#------ Density Plots of Bins- UC, CD, and nonIBD
########
#Bin645
########
#--- Remove samples which has Zero Growth Rates
PTR_meta_bin645 <- subset(PTR_meta, Bin645 > 0)
a <- ggplot(PTR_meta_bin645, aes(x = Bin645))

library("dplyr")
mu <- PTR_meta_bin645 %>% 
  group_by(diagnosis) %>%
  summarise(grp.mean = mean(Bin645))
mu

# Change line color by sex
a + geom_density(aes(color = diagnosis)) +
  scale_color_manual(values = c("darkred", "cyan2", "orange"))
# Change fill color by sex and add mean line
# Use semi-transparent fill: alpha = 0.4
a + geom_density(aes(fill = diagnosis), alpha = 0.4) +
  geom_vline(aes(xintercept = grp.mean, color = diagnosis),
             data = mu, linetype = "dashed") +
  scale_color_manual(values = c("darkred", "cyan2", "orange"))+
  scale_fill_manual(values = c("darkred", "cyan2", "orange"))

