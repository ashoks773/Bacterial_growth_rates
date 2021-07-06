#!/bin/bash

cd ~/IBD_datasets/HMP2/Bacterial_replication/Co_assembly

ls ../../WGS_rawReads/*fastq.gz | cut -d"_" -f2 | cut -d"/" -f2 | sort -u > SID.txt
grep -f SID.txt ../../hmp2_metadata.csv | grep "metagenomics" | cut -d"," -f3 > PID.txt
grep -f SID.txt ../../hmp2_metadata.csv | grep "metagenomics" | cut -d"," -f2,3 > SID_PID.txt

for i in $(cat PID.txt)
 do
  grep "$i" SID_PID.txt | cut -d"," -f1 > "$i"_sample_ID.txt
 done
mkdir Files_listing
mv *sample_ID.txt Files_listing
ls Files_listing/*sample_ID.txt > sample_list.txt



