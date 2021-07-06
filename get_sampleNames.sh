#!/bin/bash

cd ~/IBD_datasets/HMP2/Bacterial_replication/Co_assembly

ls ../../WGS_rawReads/*fastq.gz | cut -d"_" -f2 | cut -d"/" -f2 | sort -u > SID.txt
ls ../../WGS_rawReads/*_TR*.fastq.gz | cut -d"_" -f2,3 | cut -d"/" -f2 | sort -u >> SID.txt

grep -f SID.txt ../../hmp2_metadata.csv | grep "metagenomics" | cut -d"," -f3 > PID.txt
grep -f SID.txt ../../hmp2_metadata.csv | grep "metagenomics" | cut -d"," -f2,3 > SID_PID.txt

#--- Samples ending with _P were present in 1338 files (Downloaded from HMP website)
grep -v "_P" SID_PID.txt > SID_PID_filtered.txt

for i in $(cat PID.txt)
 do
  grep "$i" SID_PID_filtered.txt | cut -d"," -f1 > "$i"_sample_ID.txt
 done
mkdir Files_listing
mv *sample_ID.txt Files_listing
ls Files_listing/*sample_ID.txt > sample_list.txt



