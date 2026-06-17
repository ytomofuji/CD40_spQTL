#!/bin/bash
#SBATCH --job-name=RNA         # job name
#SBATCH --output=RNA_%A.out       # std out name
#SBATCH --error=RNA_%A.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=4         # num cpu cores / task
#SBATCH --mem=64G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short
## ---------------------------
##
## Script name: bulk_rnaseq_process_like_craft.sh
##
## Purpose of script: Wrapper script to process RNA modality of CRAFT-seq data
##
##
## Date Created: 2026-03-20
##
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

<<SKIP
mkdir -p /path/to/output_directory
cd /path/to/output_directory

dir=/path/to/output_directory
fastq1="/path/to/20260316_MET107_72h_YT13905_S82_L006_R1_001.fastq.gz"
fastq2="/path/to/20260316_MET107_72h_YT13905_S82_L006_R2_001.fastq.gz"
name=MET107_72h
SKIP


dir=$1 #dir=/path/to/output_directory
fastq1=$2
fastq2=$3
name=$4 #name=cDNA1

#load software/conda
export PATH="$PATH:/path/to/software"
eval "$(conda shell.bash hook)"

#load files
STAR=/path/to/output_directory
INDEX="/path/to/output_directory"
ADPT1="AGATCGGAAGAG"
ADPT2="AAGCAGTGGTATCAACGCAGAGTAC"


#1. merge fastq files -> skip for bulk data
mkdir -p ${dir}/${name}
cd ${dir}/${name}


#2. run fastqc for quality check
conda activate craft_seq
export JAVA_OPTS="-Xmx4g"
fastqc ${fastq1} --outdir=./
fastqc ${fastq2} --outdir=./



#3. Cutadapt
# remove adaptor
# remove polyA
# remove Q<10 bases
# remove reads < 30bp
conda activate cutadapt
cutadapt -q 10 -a ${ADPT1} -A ${ADPT1} -a ${ADPT2} -A ${ADPT2} -m 30 --poly-a \
-o trimmed.${name}.R1.fastq.gz -p trimmed.${name}.R2.fastq.gz ${fastq1} ${fastq2}



#4. run fastqc for quality check (after trimming)
conda activate craft_seq
fastqc trimmed.${name}.R1.fastq.gz --outdir=./
fastqc trimmed.${name}.R2.fastq.gz --outdir=./



#5. STAR (not solo) mapping
${STAR} --genomeDir ${INDEX} --runThreadN 4 --readFilesIn trimmed.${name}.R1.fastq.gz trimmed.${name}.R2.fastq.gz --readFilesCommand zcat \
--outFileNamePrefix ${dir}/${name}/STARout --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
--quantMode GeneCounts