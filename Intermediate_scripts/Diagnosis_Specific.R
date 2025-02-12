setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")
load ("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.RData")

#---- Diagnosis Specific Filtering
################################################################################
#------------------------------- CD ----------------------------------------
################################################################################
CD <- subset(PTR_meta, diagnosis=="CD")
#-- Remove Bin which is not present in Even 25% of the samples
CD_filtered <- dropspc(CD[,1:70], 148)
CD_filtered_meta <- cbind(CD_filtered, CD[,71:84])
write.table(CD_filtered, file="DiagnosisSpecific_output/CD_filtered_PTR.txt", sep="\t")

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
CD_bins_details <- sapply(CD_filtered_meta[,1:11], function(x) c(sum=sum(x), var=var(x), sd=sd(x)))
CD_bins_details <- data.frame(t(CD_bins_details))
jpeg("DiagnosisSpecific_output/CD_var.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (CD_bins_details$var, xlab = "Bins", ylab="Variance (log2 (PTR)- CD patients", pch=19, cex=1.5, ylim = c(0,0.25),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(CD_bins_details$var, labels=rownames(CD_bins_details), cex=0.7, font=1, pos=3)
dev.off ()


################################################################################
#------------------------------- UC ----------------------------------------
################################################################################
UC <- subset(PTR_meta, diagnosis=="UC")
#-- Remove Bin which is not present in Even 25% of the samples
UC_filtered <- dropspc(UC[,1:70], 91)
UC_filtered_meta <- cbind(UC_filtered, UC[,71:84])
write.table(UC_filtered, file="DiagnosisSpecific_output/UC_filtered_PTR.txt", sep="\t")

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


#-------- Plot variance of each Bin Across patients with Chron's Disease
#-- Calculate Varyance - Max Variance is 0.25
UC_bins_details <- sapply(UC_filtered_meta[,1:11], function(x) c(sum=sum(x), var=var(x), sd=sd(x)))
UC_bins_details <- data.frame(t(UC_bins_details))
jpeg("DiagnosisSpecific_output/UC_var.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (UC_bins_details$var, xlab = "Bins", ylab="Variance (log2 (PTR)- UC patients", pch=19, cex=1.5, ylim = c(0,0.25),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(UC_bins_details$var, labels=rownames(UC_bins_details), cex=0.7, font=1, pos=3)
dev.off ()


################################################################################
#------------------------------- nonIBD ----------------------------------------
################################################################################
nonIBD <- subset(PTR_meta, diagnosis=="nonIBD")
#-- Remove Bin which is not present in Even 25% of the samples
nonIBD_filtered <- dropspc(nonIBD[,1:70], 91)
nonIBD_filtered_meta <- cbind(nonIBD_filtered, nonIBD[,71:84])
write.table(nonIBD_filtered, file="DiagnosisSpecific_output/nonIBD_filtered_PTR.txt", sep="\t")

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


#-------- Plot variance of each Bin Across patients with Chron's Disease
#-- Calculate Varyance - Max Variance is 0.25
nonIBD_bins_details <- sapply(nonIBD_filtered_meta[,1:12], function(x) c(sum=sum(x), var=var(x), sd=sd(x)))
nonIBD_bins_details <- data.frame(t(nonIBD_bins_details))
jpeg("DiagnosisSpecific_output/nonIBD_var.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (nonIBD_bins_details$var, xlab = "Bins", ylab="Variance (log2 (PTR)- nonIBD subjects", pch=19, cex=1.5, ylim = c(0,0.25),
      col = c("red", "blue", "orange", "black", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple"))
text(nonIBD_bins_details$var, labels=rownames(nonIBD_bins_details), cex=0.7, font=1, pos=3)
dev.off ()

