#setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/Combined_Final_Scripts")

########################################################################################
#----- Part1- Check correlations between MAGs and their represetative Reference Genomes
########################################################################################

#--- Load MAG PTRs
MAG_PTR_meta <- read.csv (file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/PTR_meta.txt", row.names = 1, header = T, sep="\t")
#--- Load Reference Genomes PTRs
Ref_PTR_meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis_on_IGG_Ref/PTR_analysis/PTR_meta.txt", row.names = 1, header = T, sep="\t")

#--- Merge Both
MAG_Ref_Merged <- merge(MAG_PTR_meta, Ref_PTR_meta, by=0, all=F)
rownames(MAG_Ref_Merged) <- MAG_Ref_Merged$Row.names; MAG_Ref_Merged$Row.names <- NULL

#--- 49 MAGs have 95% Similar ANI with 49 Ref Genomes (~/Box/David_Casero_Lab/HMP_IBD_Project/Representative_Genomes or On Cluster)
#-- To do one to one Correlations we to Make the Same One to one Matrix
#- Get Ref PTRs ()
keep_Ref <- c("X1262754.3.patric","X1262776.3.patric","X1262777.3.patric","X1262992.3.patric","X1263070.3.patric","X1263107.3.patric","X1297617.4.patric","X1432052.9.patric","X1650661.3.patric","X208479.8.patric","X2562617057.img","X2562617183.img","X2654588179.img","X39488.4.patric","X59620.13.patric","X702446.3.patric","X742725.3.patric","X765821.8.patric","SRS259497_2.hgm","SRS1735557_8.hgm","SRS1596800_42.hgm","ERS608550_64.hgm","ERS396465_42.hgm","ERS631833_15.hgm","ERS235534_25.hgm","ERS608519_39.hgm","SRS1735679_4.hgm","SRS077552_27.hgm","ERS537197_59.hgm","SRS077502_20.hgm","ERS537292_27.hgm","SRS144714_48.hgm","ERS396493_5.hgm","SRS1596849_88.hgm","SRS021219_12.hgm","SRS294930_28.hgm","SRS142712_3.hgm","ERS396528_23.hgm","SRS050026_15.hgm","SRS475574_7.hgm","ERS473331_31.hgm","ERS608530_15.hgm","SRS019910_14.hgm","ERS537229_8.hgm","SRS1596798_81.hgm","SRS076929_37.hgm","ERS537221_53.hgm","ERS608493_46.hgm","SRS1719091_5.hgm")
Ref_PTR_sel <- MAG_Ref_Merged[keep_Ref]
#- Get MAG PTRs
keep_MAG <- c("Bin1291","Bin627","Bin584","Bin1962","Bin1555","Bin1633","Bin390","Bin175","Bin687","Bin197","Bin2019","Bin1549","Bin12","Bin2483","Bin3347","Bin168","Bin620","Bin2365","Bin104","Bin1265","Bin1393","Bin1479","Bin1560","Bin1677","Bin1712","Bin1811","Bin1898","Bin1908","Bin1966","Bin1970","Bin2023","Bin2057","Bin2075","Bin2132","Bin2212","Bin2333","Bin2480","Bin2582","Bin282","Bin3317","Bin3340","Bin3587","Bin3685","Bin3733","Bin3842","Bin3955","Bin418","Bin421","Bin780")
MAG_PTR_sel <- MAG_Ref_Merged[keep_MAG]

######################################################################
#-------- One to one Correaltions ** Correlation between Two matrices
#install.packages("lineup")
library (lineup)
r <- corbetw2mat(MAG_PTR_sel, Ref_PTR_sel)
r_allpairs <- corbetw2mat(MAG_PTR_sel, Ref_PTR_sel, what="bestpairs", corthresh=0.6)
r_bestright <- corbetw2mat(MAG_PTR_sel, Ref_PTR_sel, what="bestright")

#jpeg("Correlations.jpg", height = 6, width = 6, units = 'in', res = 300)
#plot (r_bestright[,1:2], col="darkred", pch=19)
#dev.off ()
#--- Or try this
#mapply(cor,as.data.frame(Ref_PTR_sel),as.data.frame(MAG_PTR_sel))

hist(r_bestright$cor, col=rgb(1,0,0,0.5), xlab="Correlation", ylab="Number of Bins", main="distribution of correlation between MAGs and Reference Genomes")

######################################################################
#----Select Random Reference Genomes and do correlations - Repeat it 10 Times
Ref_PTR <- data.frame(t(MAG_Ref_Merged[,91:278]))
random1 <- Ref_PTR[sample(nrow(Ref_PTR), 49), ]
random1 <- data.frame(t(random1))
random_cor1 <- corbetw2mat(MAG_PTR_sel, random1)
random_cor1 <- data.frame(random_cor1)
#-- Repeat this step 10 times

#random_corrs <- cbind(random_cor1, random_cor2, random_cor3, random_cor4, random_cor5, random_cor6, random_cor7, random_cor8, random_cor9, random_cor10)
#write.table(random_corrs, file="random_corrs_iter10.txt", sep="\t")
random_corrs <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/Representative_Genomes/random_corrs_iter10.txt", sep = "\t", header = T, row.names = 1)

mean_random_corrs <- rowMeans(random_corrs)
mean_random_corrs <- data.frame(mean_random_corrs)

hist(mean_random_corrs$mean_random_corrs, col=rgb(0,0,1,0.5), xlab="Correlation", ylab="Number of Bins", main="distribution of correlation between MAGs and randomly sampled Reference Genomes")

#-- Plot Them Together
jpeg("../Combined_Final_Figures/Correlations_Histogram.jpg", height = 6, width = 8, units = 'in', res = 300)
set.seed(1)

par(
  mfrow=c(1,2),
  mar=c(4,4,1,0)
)
# First distribution
hist(r_bestright$cor, col=rgb(1,0.5,0,0.5), xlab="Correlation", ylab="Number of Bins", 
     main="distribution of correlation")
lines(x = density(x = r_bestright$cor), col = "red", lty= 5, lwd=3)
# Second Distribution
hist(mean_random_corrs$mean_random_corrs, col=rgb(0,0,1,0.5), xlab="Correlation", ylab="Number of Bins", 
     main="")
lines(x = density(x = r_bestright$cor), col = "red", lty=5, lwd=3)
# Add legend
legend("topright", legend=c("Pairwise Corr","Randomly sampled"), col=c(rgb(1,0,0,0.5), 
                                                                       rgb(0,0,1,0.5)), pt.cex=2, pch=15 )
dev.off ()
