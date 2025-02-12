#!/bin/bash

#By: Ashok K. Sharma
#Date: 7th July 2021

#https://educe-ubc.github.io/MICB405/introduction-to-anvio.html
#https://astrobiomike.github.io/metagenomics/metagen_anvio
#https://merenlab.org/2016/06/18/importing-taxonomy/
#https://github.com/infphilo/centrifuge
#https://merenlab.org/tutorials/infant-gut/#inferring-taxonomy-for-metagenomes

#--- Anvio installtion - One time Task
#source /hpc/apps/miniconda/3/etc/profile.d/conda.sh
#conda activate anvio-7
# Metabat2
#conda install -c bioconda metabat2
# Maxbin2
#conda install -c bioconda maxbin2
# Concoct
#conda install -c bioconda concoct
# DAS Tool- if DAS tools need to be run with Diamond
#conda install -c bioconda das_tool
#conda install -c bioconda diamond
#- If needed then install these two also
# CheckM
#conda install -c bioconda checkm-genome
# GTDBTK
#conda install -c bioconda gtdbtk 
#conda deactivate

#--- Load Modules
#module load bowtie/2.3.2
#module load samtools/1.9

cd ~/IBD_datasets/HMP2/Bacterial_replication/Co_assembly
source /hpc/apps/miniconda/3/etc/profile.d/conda.sh
conda activate anvio-7

#######################################################################
#--- Step5- DAS tool was not working so run it Here
#- Run DAStools seperately
#--- Prepare Scaffolds2bin files for all three Binning Methods
#######################################################################
for i in $(cat sample_list.txt)
do
pid=$(echo "$i" | sed "s/_sample_ID.txt//")
pid_name=$(echo "$pid" | sed "s/Files_listing\///")
echo "$pid_name"

#-- Make contig info file for MaxBin
mkdir "$pid_name"_MEGAHIT_assembly/Binning_methods/Maxbin_fasta
cp "$pid_name"_MEGAHIT_assembly/Binning_methods/Maxbin_bins/bin_by_bin/*/*fa "$pid_name"_MEGAHIT_assembly/Binning_methods/Maxbin_fasta
grep ">" "$pid_name"_MEGAHIT_assembly/Binning_methods/Maxbin_fasta/*fa > "$pid_name"_MEGAHIT_assembly/Binning_methods/names
cut -d":" -f1 "$pid_name"_MEGAHIT_assembly/Binning_methods/names > "$pid_name"_MEGAHIT_assembly/Binning_methods/file1
cut -d":" -f1 "$pid_name"_MEGAHIT_assembly/Binning_methods/names | cut -d"/" -f4 > "$pid_name"_MEGAHIT_assembly/Binning_methods/file1
sed -i 's/.fa//g' "$pid_name"_MEGAHIT_assembly/Binning_methods/file1
cut -d":" -f2 "$pid_name"_MEGAHIT_assembly/Binning_methods/names > "$pid_name"_MEGAHIT_assembly/Binning_methods/file2
sed -i 's/>//g' "$pid_name"_MEGAHIT_assembly/Binning_methods/file2
paste "$pid_name"_MEGAHIT_assembly/Binning_methods/file2 "$pid_name"_MEGAHIT_assembly/Binning_methods/file1 > "$pid_name"_MEGAHIT_assembly/Binning_methods/maxbin2_Scaffolds2bin_format.tsv

#-- Make contig info file for Metabat2
mkdir "$pid_name"_MEGAHIT_assembly/Binning_methods/Metabat_fasta
cp "$pid_name"_MEGAHIT_assembly/Binning_methods/Metabat_bins/bin_by_bin/*/*fa "$pid_name"_MEGAHIT_assembly/Binning_methods/Metabat_fasta
grep ">" "$pid_name"_MEGAHIT_assembly/Binning_methods/Metabat_fasta/*fa > "$pid_name"_MEGAHIT_assembly/Binning_methods/names
cut -d":" -f1 "$pid_name"_MEGAHIT_assembly/Binning_methods/names > "$pid_name"_MEGAHIT_assembly/Binning_methods/file1
cut -d":" -f1 "$pid_name"_MEGAHIT_assembly/Binning_methods/names | cut -d"/" -f4 > "$pid_name"_MEGAHIT_assembly/Binning_methods/file1
sed -i 's/.fa//g' "$pid_name"_MEGAHIT_assembly/Binning_methods/file1
cut -d":" -f2 "$pid_name"_MEGAHIT_assembly/Binning_methods/names > "$pid_name"_MEGAHIT_assembly/Binning_methods/file2
sed -i 's/>//g' "$pid_name"_MEGAHIT_assembly/Binning_methods/file2
paste "$pid_name"_MEGAHIT_assembly/Binning_methods/file2 "$pid_name"_MEGAHIT_assembly/Binning_methods/file1 > "$pid_name"_MEGAHIT_assembly/Binning_methods/metabat_Scaffolds2bin_format.tsv

#-- Make contig info file for Concoct
mkdir "$pid_name"_MEGAHIT_assembly/Binning_methods/Concoct_fasta
cp "$pid_name"_MEGAHIT_assembly/Binning_methods/Concoct_bins/bin_by_bin/*/*fa "$pid_name"_MEGAHIT_assembly/Binning_methods/Concoct_fasta
grep ">" "$pid_name"_MEGAHIT_assembly/Binning_methods/Concoct_fasta/*fa > "$pid_name"_MEGAHIT_assembly/Binning_methods/names
cut -d":" -f1 "$pid_name"_MEGAHIT_assembly/Binning_methods/names > "$pid_name"_MEGAHIT_assembly/Binning_methods/file1
cut -d":" -f1 "$pid_name"_MEGAHIT_assembly/Binning_methods/names | cut -d"/" -f4 > "$pid_name"_MEGAHIT_assembly/Binning_methods/file1
sed -i 's/.fa//g' "$pid_name"_MEGAHIT_assembly/Binning_methods/file1
cut -d":" -f2 "$pid_name"_MEGAHIT_assembly/Binning_methods/names > "$pid_name"_MEGAHIT_assembly/Binning_methods/file2
sed -i 's/>//g' "$pid_name"_MEGAHIT_assembly/Binning_methods/file2
paste "$pid_name"_MEGAHIT_assembly/Binning_methods/file2 "$pid_name"_MEGAHIT_assembly/Binning_methods/file1 > "$pid_name"_MEGAHIT_assembly/Binning_methods/concoct_Scaffolds2bin_format.tsv

#--- Now run DAS Tools
DAS_Tool -i "$pid_name"_MEGAHIT_assembly/Binning_methods/metabat_Scaffolds2bin_format.tsv,"$pid_name"_MEGAHIT_assembly/Binning_methods/concoct_Scaffolds2bin_format.tsv,"$pid_name"_MEGAHIT_assembly/Binning_methods/maxbin2_Scaffolds2bin_format.tsv -l metabat2,concoct,maxbin2 -c "$pid_name"_MEGAHIT_assembly/contigs.fa -o "$pid_name"_MEGAHIT_assembly/Binning_methods/DASTool_Run1 --search_engine diamond -t 12 --write_bins

done
conda deactivate

###############################################################################################
#----Step6: Assess the quality of final non-redundant bins obtained from DAS tool using CheckM 
##############################################################################################
#-- If checkm is not installed in anvi'o then use it normally
module load python3/3.8.0
module load prodigal/2.6.3
module load pplacer/1.1.19
module load hmmer/3.3.2

for k in $(cat sample_list.txt)
do
pid=$(echo "$k" | sed "s/_sample_ID.txt//")
pid_name=$(echo "$pid" | sed "s/Files_listing\///")
echo "$pid_name"

~/.local/bin/checkm tree "$pid_name"_MEGAHIT_assembly/Binning_methods/DASTool_Run1_DASTool_bins "$pid_name"_MEGAHIT_assembly/Refined_bin_checkM -x .fa -t 24
~/.local/bin/checkm tree_qa "$pid_name"_MEGAHIT_assembly/Refined_bin_checkM -f "$pid_name"_MEGAHIT_assembly/marker_file_refined
~/.local/bin/checkm lineage_set "$pid_name"_MEGAHIT_assembly/Refined_bin_checkM "$pid_name"_MEGAHIT_assembly/marker_file_refined
~/.local/bin/checkm analyze "$pid_name"_MEGAHIT_assembly/marker_file_refined "$pid_name"_MEGAHIT_assembly/Binning_methods/DASTool_Run1_DASTool_bins "$pid_name"_MEGAHIT_assembly/Refined_bin_checkM_analyse -x .fa -t 24
~/.local/bin/checkm qa "$pid_name"_MEGAHIT_assembly/marker_file_refined "$pid_name"_MEGAHIT_assembly/Refined_bin_checkM_analyse -o 1 -t 24 > "$pid_name"_MEGAHIT_assembly/Refined_bin_Stats.txt
done

###############################################
# Step7: Taxonomic Assignment using CAT/BAT 
###############################################
#-- One time Run
#https://github.com/dutilh/CAT
#cd ~/Databases
#wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20210107.tar.gz
#tar -zxvf CAT_prepare_20210107.tar.gz

for j in $(cat sample_list.txt)
do
pid=$(echo "$j" | sed "s/_sample_ID.txt//")
pid_name=$(echo "$pid" | sed "s/Files_listing\///")
echo "$pid_name"

mkdir "$pid_name"_MEGAHIT_assembly/Refined_bin_Taxonomy
cd ~/IBD_datasets/HMP2/Bacterial_replication/Co_assembly/"$pid_name"_MEGAHIT_assembly/Refined_bin_Taxonomy

~/Softwares/CAT-master/CAT_pack/CAT bins -b ../Binning_methods/DASTool_Run1_DASTool_bins --bin_suffix .fa -d ~/Databases/CAT_prepare_20210107/2021-01-07_CAT_database -t ~/Databases/CAT_prepare_20210107/2021-01-07_taxonomy --path_to_diamond ~/Databases/CAT_prepare_20210107/Diamond_2.0.6/diamond
~/Softwares/CAT-master/CAT_pack/CAT add_names -i out.BAT.bin2classification.txt -o out.BAT.bin2classification_Names.txt -t ~/Databases/CAT_prepare_20210107/2021-01-07_taxonomy --only_official --exclude_scores
~/Softwares/CAT-master/CAT_pack/CAT summarise -i out.BAT.bin2classification_Names.txt -o BAT_summary
done

###################################################################################
# *** Not Used *** : Step7: Taxonomic Assignment using CAT/BAT or GTDB-Tk in Anvi'o 
###################################################################################
#-- For CAT/BAT check MAG_pipeline.sh or Replicaiton_rate_new.sh Script files
#-- For Anvi'o https://educe-ubc.github.io/MICB405/processing-mags.html#taxonomic-classification-using-gtdb-tk

#Create GTDB-Tk output directory
#mkdir ~/creepingFat/Data/Genome_reconstruction/GTDBTK_Output
#Create variable for output directory
#gtdbtk_out=~/CreepingFat/Data/Genome_reconstruction/GTDBTK_Output

#Establish bins’ fasta files directory 
#in_fasta=$bin_out/fasta_files

#Export GTDB-Tk data path (Check this data- where it is downloaded and installed)
#export GTDBTK_DATA_PATH=~/path_to_GTDB_data

#Operational parameters – change as required on your system
#threads=48
#pplacer_threads=24

#Run GTDB-Tk
#gtdbtk classify_wf -x fa --genome_dir $in_fasta \
 #   --out_dir $gtdbtk_out \
 #   --cpus $threads \
  #  --pplacer_cpus $pplacer_threads
 
