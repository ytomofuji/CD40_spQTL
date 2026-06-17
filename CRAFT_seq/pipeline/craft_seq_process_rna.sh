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
## Script name: craft_seq_process_rna.sh
##
## Purpose of script: Wrapper script to process RNA modality of CRAFT-seq data
##
##
## Date Created: 2025-05-06
##
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------



dir=$1 #dir=/path/to/RNA/STAR
layout=$2 # 1. R1 files 2. R2 files: layout=/path/to/RNA/layout/cDNA1.layout.tsv
name=$3 #name=cDNA1

#load software/conda
eval "$(~/miniconda3/bin/conda shell.bash hook)"

#load files
STAR=/path/to/STAR
WHITELIST="/path/to/384RNABarcodes.txt"
INDEX="/path/to/STAR_reference/"
ADPT="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


#1. merge fastq files
mkdir -p ${dir}/${name}
cd ${dir}/${name}

R1_files=`cat ${layout} | awk -F"\t" '{print $1}'`
R2_files=`cat ${layout} | awk -F"\t" '{print $2}'`

for file in ${R1_files}; do cat ${file}; done > ${name}.R1.fastq.gz
for file in ${R2_files}; do cat ${file}; done > ${name}.R2.fastq.gz



#2. run fastqc for quality check
conda activate craft_seq
export JAVA_OPTS="-Xmx4g"
fastqc ${name}.R1.fastq.gz --outdir=./
fastqc ${name}.R2.fastq.gz --outdir=./



#3. Cutadapt
# remove adaptor
# remove polyA
# remove Q<10 bases
# remove reads < 30bp
conda activate cutadapt
cutadapt -q 10 -a ${ADPT} -m 30 --poly-a -o trimmed.${name}.R2.fastq.gz ${name}.R2.fastq.gz

#basically, we assume that read name format is like
# R1: @lh00170:545:22NFCHLT4:4:1101:2294:20197 1:N:0:TCCGCGTTCA+CTAGCGCT
# R2: @lh00170:545:22NFCHLT4:4:1101:2294:20197 2:N:0:TCCGCGTTCA+CTAGCGCT
seqkit seq -n trimmed.${name}.R2.fastq.gz | sed 's/2:N:0/1:N:0/g' > ${name}.kept_ids.txt
seqkit grep -f ${name}.kept_ids.txt -n -o trimmed.${name}.R1.fastq.gz ${name}.R1.fastq.gz

rm ${name}.R1.fastq.gz
rm ${name}.R2.fastq.gz
rm ${name}.kept_ids.txt


#4. run fastqc for quality check (after trimming)
conda activate craft_seq
fastqc trimmed.${name}.R1.fastq.gz --outdir=./
fastqc trimmed.${name}.R2.fastq.gz --outdir=./



#5. STAR solo mapping
${STAR} --genomeDir ${INDEX} --runThreadN 4 --readFilesIn trimmed.${name}.R2.fastq.gz trimmed.${name}.R1.fastq.gz --readFilesCommand zcat \
--outFileNamePrefix ${dir}/${name}/STARout --outSAMtype None --outSAMmode None \
--soloType CB_UMI_Simple --soloCBwhitelist ${WHITELIST} --soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 12 \
--soloCBmatchWLtype 1MM_multi_pseudocounts --soloFeatures GeneFull