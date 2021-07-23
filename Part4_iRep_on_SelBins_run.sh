#!/bin/bash
# By: Ashok K. Sharma; 
# Date 12th July 2021

# This program will be run for All final selected bins (>90% completeness, <5% contamination, and less then 175 Scaffolds per 1Mbps) from all 106 individuals
# Rename bins to keep them in a single folder and keep the map files in Individuals folders.

module load bowtie/2.3.2

#-- Change directory
cd ~/IBD_datasets/HMP2/Bacterial_replication/Co_assembly

#- In this folder put all renamed bins from individual specific co-assemblies
ls Total_Selected_Bins/*fa > Bins_list.txt

#---- Build bowtie2 Index for each Bin
mkdir Allsamples_mapped_bins_SAM
mkdir Bins_growth_rates

for i in $(cat Bins_list.txt)
do
bid=$(echo "$i" | sed "s/Total_Selected_Bins\///")
bid_name=$(echo "$bid" | sed "s/.fa//")
echo "$bid_name"
	
	bowtie2-build "$i" "$i"
	#-- change with file name which contain all samples
	for j in $(cat check.txt)
	do
	pid=$(echo "$j" | sed "s/_sample_ID.txt//")
        pid_name=$(echo "$pid" | sed "s/Files_listing\///")
	echo "$pid_name"
	
		for k in $(cat "$j")
		do
			bowtie2 -x "$i" -1 /home/sharmaa4/IBD_datasets/HMP2/WGS_rawReads/"$k"_R1.fastq.gz -2 /home/sharmaa4/IBD_datasets/HMP2/WGS_rawReads/"$k"_R2.fastq.gz -S Allsamples_mapped_bins_SAM/"$bid_name"_"$pid_name"_"$k".sam --reorder
			~/Softwares/shrinksam-master/shrinksam -i Allsamples_mapped_bins_SAM/"$bid_name"_"$pid_name"_"$k".sam -k Allsamples_mapped_bins_SAM/"$bid_name"_"$pid_name"_"$k"-shrunk.sam
			rm Allsamples_mapped_bins_SAM/"$bid_name"_"$pid_name"_"$k".sam
		done
	done

module load python3/3.8.0
iRep -f "$i" -s Allsamples_mapped_bins_SAM/"$bid_name"*sam -o Bins_growth_rates/"$bid_name".iRep

done


