#!/bin/bash

#--- List.txt file will be generated which will contain all the File names with Patient IDs. In each patient ID files names - all the names belonging to that patient will bewritten. I have used the following Script file for that purpose.
# sh get_sampleNames.sh

module load bowtie/2.3.2
module load samtools/1.9

cd ~/IBD_datasets/HMP2/Bacterial_replication/Co_assembly

#*************** Start for loop
for i in $(cat sample_list.txt)
do
pid=$(echo "$i" | sed "s/_sample_ID.txt//")
pid_name=$(echo "$pid" | sed "s/Files_listing\///")
echo "$pid_name"
 
#***** Start for loop 
   for j in $(cat "$i")
    do
	cp /home/sharmaa4/IBD_datasets/HMP2/WGS_rawReads/"$j"_*.fastq.gz .
    done
#*****End for loop

   cat *_R1.fastq.gz > "$pid_name"_R1_Merged.fastq.gz
   cat *_R2.fastq.gz > "$pid_name"_R2_Merged.fastq.gz

###############################
#----- Step1: Generate Contigs
###############################
~/Softwares/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 "$pid_name"_R1_Merged.fastq.gz -2 "$pid_name"_R2_Merged.fastq.gz  -m 0.99 --mem-flag 0 -o "$pid_name"_MEGAHIT_assembly -t 12

##############################################
#----- Step2: Generate Contig db using Anvi'o
##############################################
source /hpc/apps/miniconda/3/etc/profile.d/conda.sh
conda activate anvio-7
anvi-script-reformat-fasta "$pid_name"_MEGAHIT_assembly/final.contigs.fa -o "$pid_name"_MEGAHIT_assembly/contigs.fa -l 1000 --simplify-name

anvi-gen-contigs-database -f "$pid_name"_MEGAHIT_assembly/contigs.fa -o "$pid_name"_MEGAHIT_assembly/contigs.db -n "my metagenome" -T 12
anvi-run-hmms -c "$pid_name"_MEGAHIT_assembly/contigs.db -I Campbell_et_al -T 12
#anvi-run-hmms -c "$pid_name"_MEGAHIT_assembly/contigs.db -I Bacteria_71 -T 12
#anvi-run-scg-taxonomy -c "$pid_name"_MEGAHIT_assembly/contigs.db -T 12
anvi-setup-ncbi-cogs -T 12 --just-do-it
anvi-run-ncbi-cogs -c "$pid_name"_MEGAHIT_assembly/contigs.db  -T 12
#--- Assign taxonomy with tool called Centrifuge
anvi-get-sequences-for-gene-calls -c "$pid_name"_MEGAHIT_assembly/contigs.db -o "$pid_name"_MEGAHIT_assembly/gene_calls.fa
#-- Use Kaiju to import taxonomy
~/Softwares/kaiju/bin/kaiju -t ~/Softwares/kaiju/kaijudb/nodes.dmp -f ~/Softwares/kaiju/kaijudb/kaiju_db.fmi -i "$pid_name"_MEGAHIT_assembly/gene_calls.fa -o "$pid_name"_MEGAHIT_assembly/gene_calls_nr.out -z 16 -v
~/Softwares/kaiju/bin/kaiju-addTaxonNames -t ~/Softwares/kaiju/kaijudb/nodes.dmp -n ~/Softwares/kaiju/kaijudb/names.dmp -i "$pid_name"_MEGAHIT_assembly/gene_calls_nr.out -o "$pid_name"_MEGAHIT_assembly/gene_calls_nr.names -r superkingdom,phylum,order,class,family,genus,species

anvi-import-taxonomy-for-genes -i "$pid_name"_MEGAHIT_assembly/gene_calls_nr.names -c "$pid_name"_MEGAHIT_assembly/contigs.db -p ~/Softwares/kaiju/bin/kaiju --just-do-it
centrifuge -f -x /home/sharmaa4/Softwares/centrifuge/indices/p_compressed+h+v "$pid_name"_MEGAHIT_assembly/gene_calls.fa -S "$pid_name"_MEGAHIT_assembly/centrifuge_hits.tsv
#Check if this command not work then remove pid name before the centrifuge report.tsv file
anvi-import-taxonomy-for-genes -c "$pid_name"_MEGAHIT_assembly/contigs.db -i "$pid_name"_MEGAHIT_assembly/"$pid_name"_centrifuge_report.tsv "$pid_name"_MEGAHIT_assembly/centrifuge_hits.tsv -p centrifuge
#rm -r "$pid_name"_MEGAHIT_assembly/SAM/*ANVIO_PROFILE
conda deactivate

#########################################################
#----Step3: Map reads on contigs to Generat Sam/Bam Files
########################################################

bowtie2-build "$pid_name"_MEGAHIT_assembly/contigs.fa "$pid_name"_MEGAHIT_assembly/contigs.fa

mkdir "$pid_name"_MEGAHIT_assembly/SAM
#ls *_R1.fastq.gz | cut -d"_" -f1 > "$pid_name"_SRR.txt
for k in $(cat "$i")
 do
  bowtie2 -q -x "$pid_name"_MEGAHIT_assembly/contigs.fa -1 ~/IBD_datasets/HMP2/WGS_rawReads/"$k"_R1.fastq.gz -2 ~/IBD_datasets/HMP2/WGS_rawReads/"$k"_R2.fastq.gz --no-unal -p 12 -S "$pid_name"_MEGAHIT_assembly/SAM/"$k".sam
  samtools view -b -o "$pid_name"_MEGAHIT_assembly/SAM/"$k"_raw.bam "$pid_name"_MEGAHIT_assembly/SAM/"$k".sam
  samtools sort -o "$pid_name"_MEGAHIT_assembly/SAM/"$k".bam "$pid_name"_MEGAHIT_assembly/SAM/"$k"_raw.bam
  samtools index "$pid_name"_MEGAHIT_assembly/SAM/"$k".bam
  samtools view -h SAM/"$k".bam > SAM/"$k"_sort.sam
 done

#---- Make coverage files for MaxBin and MetaBat if needed in the future
module load java/1.8.0_221
mkdir "$pid_name"_MEGAHIT_assembly/Abundance
for m in $(cat "$i")
 do
~/Softwares/bbmap/pileup.sh in="$pid_name"_MEGAHIT_assembly/SAM/"$m".sam out="$pid_name"_MEGAHIT_assembly/Abundance/"$m".cov.txt
awk '{print $1"\t"$5}' "$pid_name"_MEGAHIT_assembly/Abundance/"$m".cov.txt | grep -v '^#' > "$pid_name"_MEGAHIT_assembly/Abundance/"$m".abundance.txt
done
ls "$pid_name"_MEGAHIT_assembly/Abundance/*abundance.txt > "$pid_name"_MEGAHIT_assembly/abund_list_file

#---Now remove sam and raw bam files
#rm "$pid_name"_MEGAHIT_assembly/SAM/*sam
rm "$pid_name"_MEGAHIT_assembly/SAM/*_raw.bam

module load metabat/2.12.1
jgi_summarize_bam_contig_depths "$pid_name"_MEGAHIT_assembly/SAM/*bam --outputDepth "$pid_name"_MEGAHIT_assembly/depth.txt

###########################################################
#----Step4: Run Anvi'o to generate proifles for each sample
###########################################################
source /hpc/apps/miniconda/3/etc/profile.d/conda.sh
conda activate anvio-7
for n in $(cat "$i")
 do
   anvi-profile -i "$pid_name"_MEGAHIT_assembly/SAM/"$n".bam -c "$pid_name"_MEGAHIT_assembly/contigs.db -T 12
 done
anvi-merge "$pid_name"_MEGAHIT_assembly/SAM/*ANVIO_PROFILE/PROFILE.db -o "$pid_name"_MEGAHIT_assembly/SAMPLES-MERGED -c "$pid_name"_MEGAHIT_assembly/contigs.db
#rm -r "$pid_name"_MEGAHIT_assembly/SAM/*ANVIO_PROFILE

###########################################################################################
#----Step5: Use profile_db and contig_db to do Binning and remove redundancy using DAS_Tool
###########################################################################################
#---- Run Binning and Refinement using Anvio
profile_db="$pid_name"_MEGAHIT_assembly/SAMPLES-MERGED/PROFILE.db
contig_db="$pid_name"_MEGAHIT_assembly/contigs.db
 
mkdir "$pid_name"_MEGAHIT_assembly/Binning_methods
bin_out="$pid_name"_MEGAHIT_assembly/Binning_methods
 
threads=12
 
#Clustering with metabat2
anvi-cluster-contigs -p $profile_db -c $contig_db -C metabat --driver metabat2 -T $threads --just-do-it
anvi-summarize -p $profile_db -c $contig_db -C metabat -o $bin_out/Metabat_bins
 
#Clustering with concoct
anvi-cluster-contigs -p $profile_db -c $contig_db -C concoct --driver concoct -T $threads --just-do-it
anvi-summarize -p $profile_db -c $contig_db -C concoct -o $bin_out/Concoct_bins
 
#Clustering with maxbin2
anvi-cluster-contigs -p $profile_db -c $contig_db -C maxbin --driver maxbin2 -T $threads --just-do-it
anvi-summarize -p $profile_db -c $contig_db -C maxbin -o $bin_out/Maxbin_bins
 
#Consensus binning with DASTOOL
anvi-cluster-contigs -p $profile_db -c $contig_db --driver dastool -S metabat,concoct,maxbin -C dastool --just-do-it -T $threads --search_engine diamond
 
#Obtain bins as fasta files
anvi-summarize -p $profile_db -c $contig_db -C dastool -o $bin_out/dastool_out
 
#Copy all bin fasta files for downstream processing such as quality estimation, taxonomic classification, and functional annotation
#Make directory for fasta files for easy access
mkdir $bin_out/fasta_files
cp $bin_out/dastool_out/*/*fasta $bin_out/fasta_files
conda deactivate

#rm *R*_Merged.fastq.gz
#---Important to remove all R1 and R2 files belong to First Patient ID. Otherwise these R1 and R2 will be merged with samples files of the Second patient.
rm *_R1.fastq.gz
rm *_R2.fastq.gz
done
#************** End for Loop
