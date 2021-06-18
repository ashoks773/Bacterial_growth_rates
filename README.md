# Bacterial_growth_rates
Computational pipeline to calculate bacterial growth rates from metagenomic samples
By: Ashok K. Sharma

Any questions about this pipeline please mail to *ashoks773@gmail.com*

## Pacakges to be installed:
* Check the following modules at your HPC and install them manually if they are not available.
  *  **Kneaddata:** https://huttenhower.sph.harvard.edu/kneaddata/ (If reads are not quality filtered)
  *  **Trimmomatic:** http://www.usadellab.org/cms/?page=trimmomatic
  *  **Bowtie2:** https://www.metagenomics.wiki/tools/bowtie2/install
  *  **Samtools:** https://gist.github.com/adefelicibus/f6fd06df1b4bb104ceeaccdd7325b856
  *  **Megahit:** https://github.com/voutcn/megahit
  *  **MetaBAT2:** https://bitbucket.org/berkeleylab/metabat/src/master/
  *  **Prodigal:** https://github.com/hyattpd/Prodigal/wiki/installation
  *  **pplacer:** https://github.com/matsen/pplacer
  *  **hmmer:** http://hmmer.org/documentation.html
  *  **CheckM:** https://github.com/Ecogenomics/CheckM
  *  **DAS_Tools:** https://github.com/cmks/DAS_Tool
  *  **CAT and BAT:** https://github.com/dutilh/CAT
  *  **Phython3:** https://www.python.org/downloads/
  *  **Java:** https://docs.oracle.com/en/java/javase/11/install/installation-jdk-linux-platforms.html#GUID-737A84E4-2EFF-4D38-8E60-3E29D1B884B8
  *  **Perl:** https://www.perl.org/get.html
  
  ## Step1: Co-assembly: to generate contigs/scaffolds using MEGAHIT
  Before running this step, raw metagenomic reads need to be processed for **1.Qaulity filtering: Using Trimmomatic** and **2. Host contamination removal: using  Bowtie2 and Samtools**. Or kneaddata https://huttenhower.sph.harvard.edu/kneaddata/ can be used for this purpose. Kneaddata can be run using the following command. 
``` r
kneaddata --input ${SAMPLE}_1.fastq.gz --input ${SAMPLE}_2.fastq.gz -db /home/sharmaa4/Databases/Human_Genome/human_GRCh38.p13 --output /home/sharmaa4/IBD_datasets/HMP2/step1_WGS_rawReads_QC/${SAMPLE}_kneaddata_output --bypass-trf
```
Output of kneaddata: will generate multiple files per sample. We need to use the following four files per sample: **paired_1.fastq, paired_2.fastq, unmatched_1.fastq, unmatched_2.fastq**. All four files can be merged to generate one file for each sample if someone want to use Humann3 (https://huttenhower.sph.harvard.edu/humann/) or MetaPhlAn 3.0 (https://huttenhower.sph.harvard.edu/metaphlan/) for functional and taxonomic classification respectively.
**Or** paired end files can be generated by concatenating paired_1.fastq and unmatched_1.fastq and paired_2.fastq and unmatched_2.fastq files. 
``` r
cat ${SAMPLE}_paired_1.fastq ${SAMPLE}_unmatched_1.fastq > ${SAMPLE}_1.fastq
cat ${SAMPLE}_paired_2.fastq ${SAMPLE}_unmatched_2.fastq > ${SAMPLE}_2.fastq
gzip *fastq
```

**We are starting from here:** For this tutorial, quality filtered paired end reads (_1.fastq.gz and _2.fastq.gz_) will be used. We have total 20 fastq.gz files for 10 Samples. This directly will be used for downstream processing.
``` r
module load bowtie/2.3.2
module load samtools/1.9

cd ~/IBD_datasets/HMP2/Bacterial_replication/Sample_Dataset/

~/Softwares/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 Samples/SRR5936130_1.fastq.gz,Samples/SRR5936131_1.fastq.gz,Samples/SRR5936132_1.fastq.gz,Samples/SRR5936133_1.fastq.gz,Samples/SRR5936134_1.fastq.gz,Samples/SRR5936135_1.fastq.gz,Samples/SRR5936136_1.fastq.gz,Samples/SRR5936137_1.fastq.gz,Samples/SRR5936138_1.fastq.gz,Samples/SRR5936139_1.fastq.gz -2 Samples/SRR5936130_2.fastq.gz,Samples/SRR5936131_2.fastq.gz,Samples/SRR5936132_2.fastq.gz,Samples/SRR5936133_2.fastq.gz,Samples/SRR5936134_2.fastq.gz,Samples/SRR5936135_2.fastq.gz,Samples/SRR5936136_2.fastq.gz,Samples/SRR5936137_2.fastq.gz,Samples/SRR5936138_2.fastq.gz,Samples/SRR5936139_2.fastq.gz -m 100000000000 -o MEGAHIT_assembly -t 120
``` 
Or concatenate files before running MEGAHIT
``` r
cat *_1.fastq.gz > allsample_1.fastq.gz
cat *_2.fastq.gz > allsample_2.fastq.gz
~/Softwares/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 allsample_1.fastq.gz -2 allsample_2.fastq.gz -m 100000000000 -o MEGAHIT_assembly -t 120
``` 

# Step2: Alignmnet of reads on assembled contigs
This step includes indexing of contigs files, followed by mapping of reads to generate sam and bam files which needs to be sorted and indexed to be used for other packages such as MetaBAT2 and MaxBin (for construction of MAGs) and DEMIC (to calculate Bacterial growth rates)
``` r
bowtie2-build MEGAHIT_assembly/final.contigs.fa MEGAHIT_assembly/final.contigs.fa
mkdir SAM # Make a folder name SAM to store alignment outputs
cd Samples
for file in *_1.fastq.gz;
do
  SAMPLE=$(echo ${file} | sed "s/_1.fastq.gz//")
  echo ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz
  bowtie2 -q -x MEGAHIT_assembly/final.contigs.fa -1 ${SAMPLE}_1.fastq.gz -2 ${SAMPLE}_2.fastq.gz --no-unal -p 120 -S SAM/${SAMPLE}.sam
  samtools view -b -o SAM/${SAMPLE}.bam SAM/${SAMPLE}.sam
  samtools sort -o SAM/${SAMPLE}_sort.bam SAM/${SAMPLE}.bam
  samtools index SAM/${SAMPLE}_sort.bam
  samtools view -h SAM/${SAMPLE}_sort.bam > SAM/${SAMPLE}_sort.sam #--We need sorted sam files - to run DEMIC in step4
done
```



``` 
