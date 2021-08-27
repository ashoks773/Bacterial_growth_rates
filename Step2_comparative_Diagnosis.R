######################################################################################################
# ------- This scipt is to compare Bins across CD, UC, nonIBD and check their Distributions ----------
######################################################################################################

####################
#Comprative Analysis
####################
setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")
load ("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.RData")

#################################################################################################
#---- To get Bins which are Significantly Discriminating among three Groups: CD, UC, and nonIBD
################################################################################################
library (labdsv)

iva <- indval(PTR_meta[,1:70], PTR_meta$diagnosis)

gr <- iva$maxcls[iva$pval<=0.05]
iv <- iva$indcls[iva$pval<=0.05]
pv <- iva$pval[iva$pval<=0.05]
fr <- apply(PTR_meta[,1:70]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
write.table (indvalsummary, file="Bins_Diagnosis_indvalsummary.txt", sep = "\t")

##########
#-- For CD
##########
#-----Bin645
PTR_meta_bin645 <- subset(PTR_meta, Bin645 > 0)
boxplot(PTR_meta_bin645$Bin645 ~ PTR_meta_bin645$diagnosis)
kruskalmc(PTR_meta_bin645$Bin645, PTR_meta_bin645$diagnosis)

#-----Bin2019
PTR_meta_bin2019 <- subset(PTR_meta, Bin2019 > 0)
boxplot(PTR_meta_bin2019$Bin2019 ~ PTR_meta_bin2019$diagnosis)
kruskalmc(PTR_meta_bin2019$Bin2019, PTR_meta_bin2019$diagnosis)

##########
#-- For UC
##########
#-----Bin3140
PTR_meta_bin3140 <- subset(PTR_meta, Bin3140 > 0)
boxplot(PTR_meta_bin3140$Bin3140 ~ PTR_meta_bin3140$diagnosis)
kruskalmc(PTR_meta_bin3140$Bin3140, PTR_meta_bin3140$diagnosis)

#----Bin1510
PTR_meta_bin1510 <- subset(PTR_meta, Bin1510 > 0)
boxplot(PTR_meta_bin1510$Bin1510 ~ PTR_meta_bin1510$diagnosis)
kruskalmc(PTR_meta_bin1510$Bin1510, PTR_meta_bin1510$diagnosis)


#####################################
#------ If need to make density plot
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

