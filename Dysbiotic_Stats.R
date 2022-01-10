#-- Dysbiotic Stats
#---- Calculate Disease Activity Index on All sets - Bins, RefGenomes - PTRs and Relative abundnaces, and Metaphln relab also
setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/Disease_Activity")

source("Scripts/metaphln_AI.R") # Calculated on the relative abundances Calculated by Authors using Metaphln ** Downloaded from iHMP website **
#source("Scripts/binsPTR_AI.R")  #Calculated on PTR estimates of 70 Bins (present in atleast 5% of the samples)
#source("Scripts/binsRelab_AI.R") #Calculated on Relative abundance estimates of Above 70 Bins 
source("Scripts/RefPTR_AI.R") #Caclcualted on PTR estimates of 188 Reference Genomes (present in atleast 5% of the samples)
#source("Scripts/RefRelab_AI.R") #Caclcualted on Relative abundance of 188 Reference Genomes (present in atleast 5% of the samples)

# Dysbiosis stats
source('~/Box/David_Casero_Lab/HMP_IBD_Project/hmp2_analysis/pcl_utils.R')
# https://bitbucket.org/biobakery/hmp2_analysis/src/master/disease_activity/src/define_disease_activity.r

all_dysbiotic <- c()
dysbiosis_dt <- data.frame()
for (sub in Ref_PTR_meta$site_sub_coll) {
  Ref_PTR_meta.sub <- Ref_PTR_meta %>% pcl.filter.s(site_sub_coll==sub) %>% pcl.sort.s(collection)
  #Ref_PTR_meta.sub <- subset(Ref_PTR_meta, site_sub_coll==sub)
  
  act_change <- diff(Ref_PTR_meta.sub$active)
  change_wk <- Ref_PTR_meta.sub$week_num[act_change!=0]
  change_delta <- diff(change_wk)
  
  dysbiosis_dt <- rbind(dysbiosis_dt, data.frame(
    active = (seq_along(change_delta) + Ref_PTR_meta.sub$active[1]) %% 2 == 1,
    delta = change_delta,
    censored = rep(F, length(change_delta)),
    diagnosis = rep(Ref_PTR_meta.sub$diagnosis[1], length(change_delta))
  ))
  if (length(change_wk) >= 1) {
    dysbiosis_dt <- rbind(dysbiosis_dt, data.frame(
      active = Ref_PTR_meta.sub$active[Ref_PTR_meta.sub$ns],
      delta = Ref_PTR_meta.sub$week_num[Ref_PTR_meta.sub$ns] - change_wk[length(change_wk)],
      censored = T,
      diagnosis = Ref_PTR_meta.sub$diagnosis[1]
    ))
  }
  
  if (all(Ref_PTR_meta.sub$active)) {
    all_dysbiotic <- c(all_dysbiotic, sub)
  }
}

AllSets_Activity_Index <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/Disease_Activity/AllSets_Activity_Index.txt", sep = "\t", row.names = 1, header = T)
# Load Meta
meta <- read.csv(file="~/Box/David_Casero_Lab/HMP_IBD_Project/PTR_analysis/hmp2_metadata_mgx_edited.csv", header = T, row.names = 1)

AllSets_Activity_Index_meta <- merge(AllSets_Activity_Index, meta, by=0, all=F)
rownames(AllSets_Activity_Index_meta) <- AllSets_Activity_Index_meta$Row.names; AllSets_Activity_Index_meta$Row.names <- NULL

#--- Plot Dysbiosis scores calculated only on binsPTRs and bins Relative abundances
range_d <- range (list(C3001_diseaseAI[,2:3]))
jpeg("PatientSpceific_output/C3001_Dysbiosis_score.jpg", height = 5, width = 7, units = 'in', res = 600)
plot(C3001_diseaseAI$visit_num, C3001_diseaseAI$binsPTR, type = "b", frame = FALSE, pch = 19, 
     col = "darkgoldenrod4", xlab = "visits", ylab = "Dysbiosis Score", ylim = c(range_d[1]-0.2, range_d[2]+0.1)) 
lines(C3001_diseaseAI$visit_num, C3001_diseaseAI$binsRelab, pch = 19, col = "Maroon", type = "b", lty = 1)
# Add a legend to the plot
legend("topright", legend=c("MAGs_PTRs","MAGs_Relabun"),
       col=c("darkgoldenrod4", "Maroon"),
       lty = 0.5, cex=0.5, ncol = 2:2, pch=19)
abline(h=binsPTR_disease_activity_threshold, col="darkgoldenrod4", lwd=1, lty=2)
abline(h=BinsRelab_disease_activity_threshold, col="Maroon", lwd=1, lty=2)
dev.off ()

#-- Take only CD and UC patients
CD <- subset(AllSets_Activity_Index_meta, diagnosis == "CD")
UC <- subset(AllSets_Activity_Index_meta, diagnosis == "UC")
alL_CD_UC <- rbind(CD, UC)

# https://stackoverflow.com/questions/3777174/plotting-two-variables-as-lines-using-ggplot2-on-the-same-graph
#g <- alL_CD_UC %>%
 # gather(key,value, metaphln, Ref_PTR) %>%
 ## ggplot(aes(x=visit_num, y=value, colour=key)) + 
 # geom_point() + 
 # labs(title = "", subtitle = "", y = "Dysbiosis score", x = "") + theme_bw() +
  #facet_grid(. ~ Participant.ID) +
 # geom_line()
# Change line type and color  geom_hline(data = lt,
#p <- g + geom_hline(yintercept=c(0.84, 0.75), linetype="dashed", color = "black")


pdf("CD_UC_Dysbiosis9.pdf")
alL_CD_UC %>%
  gather(key,value, metaphln, Ref_PTR) %>%
  ggplot(aes(x=visit_num, y=value, colour=key)) + 
  geom_point() + 
  labs(title = "", subtitle = "", y = "Dysbiosis score", x = "") + theme_bw() +
  facet_wrap_paginate(~ Participant.ID, ncol = 3, nrow = 3, page = 9) + geom_line() +
  geom_hline(yintercept=c(0.84, 0.75), linetype="dashed", color = "black")
dev.off ()

#---