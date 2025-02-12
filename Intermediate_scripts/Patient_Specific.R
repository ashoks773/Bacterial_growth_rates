setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis")
load ("~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.RData")
#---- Patient Specific Filtering
C3001 <- subset(PTR_meta_no_filtered, pid=="C3001")
#-- Remove Bin which is not present in Even 50% of the samples
C3001_filtered <- dropspc(C3001[,1:1503], 4)
write.table(C3001_filtered, "PatientSpceific_output/C3001_PTR_filtered.txt", sep = "\t")
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
C3001_bins_details <- sapply(C3001_filtered_meta[,1:14], function(x) c(sum=sum(x), var=var(x), sd=sd(x)))
C3001_bins_details <- data.frame(t(C3001_bins_details))
jpeg("PatientSpceific_output/C3001_var.jpg", height = 5, width = 7, units = 'in', res = 600)
plot (C3001_bins_details$var, xlab = "Bins", ylab="Variance (log2 (PTR) among different visits", pch=19, cex=1.5, ylim = c(0,1.1),
      col = c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black"))
text(C3001_bins_details$var, labels=rownames(C3001_bins_details), cex=0.7, font=1, pos=3)
dev.off ()


#-------
C3001 <- subset(PTR_meta_no_filtered, pid=="C3003")
C3001_filtered <- dropspc(C3001[,1:1503], 4)
C3001_filtered_meta <- cbind(C3001_filtered, C3001[,1504:1517])









##################
#---- Using GGPLOT
##################
colcode <- c("red", "blue", "orange", "darkred", "cyan", "darkgray", "Magenta", "Maroon", "limegreen", "darkgreen", "purple", "navy", "darkgoldenrod4", "black")
plot <- ggplot(C3001_bins_details, aes(x = row.names(C3001_bins_details), y = var, colour=row.names(C3001_bins_details))) + geom_point() + theme_bw()
plot1 <- plot + scale_color_manual(values=colcode)
countries_sp +
  geom_label_repel(aes(label = row.names(C3001_bins_details)), size = 2)

####################
#----- Using GGPLOT
###################
C3001_filtered_visit <- C3001_filtered[,c(1:14,16)]
C3001_filtered_melted <- melt(C3001_filtered_visit, id.vars = "visit")
library(ggplot2)
library(directlabels)

ggplot(C3001_filtered_melted, aes(x = visit, y = value, group = variable, colour = variable)) + 
  geom_line() +
  scale_colour_discrete(guide = 'none') + theme_bw()+
  scale_x_discrete(expand=c(0, 1)) +
  geom_dl(aes(label = variable), method = list(dl.combine("first.points", "last.points")), cex = 0.8)

p + geom_label_repel()
