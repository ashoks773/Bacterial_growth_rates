#!/bin/bash
# By: Ashok K. Sharma;  
# Date 12th July 2021 

#----------------##########------------------------------------------
# This step is to calcuate bacterial Growth rates using CoPTR method
# This program will be run for All final selected bins (>90% completeness, <5% contamination, and less then 175 Scaffolds per 1Mbps) from all 106 individuals 
# Rename bins to keep them in a single folder and keep the map files in Individuals folders. 

cd ~/home/sharmaa4/IBD_datasets/HMP2/Bacterial_replication/Co_assembly
source /hpc/apps/miniconda/3/etc/profile.d/conda.sh
conda activate coptr

mkdir CoPTR_map_outputs
#--Make sure the Bins (MAGs) clustered by 95% ANI and select one representative genome per species
coptr index --bt2-threads 12 --bt2-packed Total_Selected_Bins CoPTR_map_outputs/bins_indexed
mkdir CoPTR_map_outputs/bam
#--Make Sure your paired end reads ends with *_1.fastq.gz and *_2.fastq.gz
coptr map --threads 12 --paired CoPTR_map_outputs/bins_indexed /home/sharmaa4/IBD_datasets/HMP2/test_fastq CoPTR_map_outputs/bam
mkdir CoPTR_map_outputs/coverage-maps
coptr extract CoPTR_map_outputs/bam CoPTR_map_outputs/coverage-maps
coptr estimate CoPTR_map_outputs/coverage-maps PTR_output --min-reads 5000
conda deactivate


