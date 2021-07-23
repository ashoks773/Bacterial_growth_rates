#!/bin/bash
# By: Ashok K. Sharma; 
# Date 12th July 2021

cd ~/IBD_datasets/HMP2/Bacterial_replication/Co_assembly/Total_Bins

###############################################################################################
#----Step6: Assess the quality of final non-redundant bins obtained from DAS tool using CheckM 
##############################################################################################
#-- If checkm is not installed in anvi'o then use it normally
module load python3/3.8.0
module load prodigal/2.6.3
module load pplacer/1.1.19
module load hmmer/3.3.2

#-- Bins from all patients were collected and renamed using the following scripts
# sh metabat_bins.sh
# sh maxbin_bins.sh
# sh concoct_bins.sh
# sh refined_bins.sh

#-- For Metabat
~/.local/bin/checkm tree Metabat_all_Bins Metabat_all_Bins_Checkm -x .fa -t 24
~/.local/bin/checkm tree_qa Metabat_all_Bins_Checkm -f marker_metabat
~/.local/bin/checkm lineage_set Metabat_all_Bins_Checkm marker_metabat
~/.local/bin/checkm analyze marker_metabat Metabat_all_Bins Metabat_all_Bins_Checkm_analyse -x .fa -t 24
~/.local/bin/checkm qa marker_metabat Metabat_all_Bins_Checkm_analyse -o 1 -t 24 Metabat_bins_Stats.txt

#-- For MaxBin
~/.local/bin/checkm tree Maxbin_all_Bins Maxbin_all_Bins_Checkm -x .fa -t 24
~/.local/bin/checkm tree_qa Maxbin_all_Bins_Checkm -f marker_maxbin
~/.local/bin/checkm lineage_set Maxbin_all_Bins_Checkm marker_maxbin
~/.local/bin/checkm analyze marker_maxbin Maxbin_all_Bins Maxbin_all_Bins_Checkm_analyse -x .fa -t 24
~/.local/bin/checkm qa marker_maxbin Maxbin_all_Bins_Checkm_analyse -o 1 -t 24 Maxbin_bins_Stats.txt

#-- For Concoct
~/.local/bin/checkm tree Concoct_all_Bins Concoct_all_Bins_Checkm -x .fa -t 24
~/.local/bin/checkm tree_qa Concoct_all_Bins_Checkm -f marker_concoct
~/.local/bin/checkm lineage_set Concoct_all_Bins_Checkm marker_concoct
~/.local/bin/checkm analyze marker_concoct Concoct_all_Bins Concoct_all_Bins_Checkm_analyse -x .fa -t 24
~/.local/bin/checkm qa marker_concoct Concoct_all_Bins_Checkm_analyse -o 1 -t 24 Concoct_bins_Stats.txt

#-- For Merged Refined Bins DASTools
~/.local/bin/checkm tree Refined_Bins Refined_Bins_Checkm -x .fa -t 24
~/.local/bin/checkm tree_qa Refined_Bins_Checkm -f marker_refined
~/.local/bin/checkm lineage_set Refined_Bins_Checkm marker_refined
~/.local/bin/checkm analyze marker_refined Refined_Bins Refined_Bins_Checkm_analyse -x .fa -t 24
~/.local/bin/checkm qa marker_refined Refined_Bins_Checkm_analyse -o 1 -t 24 Refined_bins_Stats.txt

#----- Assign taxonomy using CAT/BAT
mkdir Refined_Bins_Taxonomy
cd ~/IBD_datasets/HMP2/Bacterial_replication/Co_assembly/Total_Bins/Refined_Bins_Taxonomy

~/Softwares/CAT-master/CAT_pack/CAT bins -b Refined_Bins --bin_suffix .fa -d ~/Databases/CAT_prepare_20210107/2021-01-07_CAT_database -t ~/Databases/CAT_prepare_20210107/2021-01-07_taxonomy --path_to_diamond ~/Databases/CAT_prepare_20210107/Diamond_2.0.6/diamond
~/Softwares/CAT-master/CAT_pack/CAT add_names -i out.BAT.bin2classification.txt -o out.BAT.bin2classification_Names.txt -t ~/Databases/CAT_prepare_20210107/2021-01-07_taxonomy --only_official --exclude_scores
~/Softwares/CAT-master/CAT_pack/CAT summarise -i out.BAT.bin2classification_Names.txt -o BAT_summary

rm out.BAT.concatenated.alignment.diamond

