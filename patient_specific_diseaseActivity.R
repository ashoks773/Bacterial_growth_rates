#########################################################################
##### Patient Specific Analysis ################################
#########################################################################
#---- For DysBiosys Scores:
#-- We can not use Dysbiosys scores calculated on only 70 Bins- Because these common Bins might not be a reprentative of Each Patient
# Take 70 Bins and Bins which are specific to C3001; Repeat same for All samples

#--- Only thing needs to be changed is the number of lines in the last two plot based on Number of Bins present in 50% of the samples
#-- And box size based on the 

#-- For these many samples we need to Repeat this Process
sort (unique (bugs_meta$Participant.ID))
#[1] "C3001" "C3002" "C3003" "C3004" "C3005" "C3006" "C3008" "C3009" "C3010" "C3011" "C3012" "C3013" "C3015" "C3016" "C3017" "C3021" "C3022"
#[18] "C3023" "C3027" "C3028" "C3029" "C3030" "C3031" "C3032" "C3034" "C3035" "C3037" "E5001" "E5004" "E5009" "E5013" "H4001" "H4004" "H4006"
#[35] "H4007" "H4008" "H4009" "H4010" "H4013" "H4014" "H4015" "H4016" "H4017" "H4018" "H4019" "H4020" "H4022" "H4023" "H4024" "H4027" "H4028"
#[52] "H4030" "H4031" "H4032" "H4035" "H4038" "H4039" "H4040" "H4042" "H4043" "H4044" "H4045" "M2008" "M2014" "M2021" "M2025" "M2026" "M2027"
#[69] "M2028" "M2034" "M2039" "M2041" "M2042" "M2047" "M2048" "M2060" "M2061" "M2064" "M2068" "M2069" "M2071" "M2072" "M2075" "M2077" "M2079"
#[86] "M2083" "M2084" "M2085" "M2097" "M2103" "P6005" "P6009" "P6010" "P6012" "P6013" "P6014" "P6016" "P6017" "P6018" "P6024" "P6025" "P6028"
#[103] "P6033" "P6035" "P6037" "P6038"


#-- Make a Unique List
setwd("~/Box/David_Casero_Lab/HMP_IBD_Project/Disease_Activity")
source("Scripts/C3001_Dysbiosis_and_bins.R") # 14 Bins are present in 50% of the samples so check you are using 14
source("Scripts/C3002_Dysbiosis_and_bins.R") # 6 Bins are present in 50% of the samples so check you are using 6
source("Scripts/C3003_Dysbiosis_and_bins.R") # 19 Bins are present in 50% of the samples so check you are using 19
