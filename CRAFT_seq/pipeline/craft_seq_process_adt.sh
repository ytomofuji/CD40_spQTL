#!/bin/bash
#SBATCH --job-name=ADT         # job name
#SBATCH --output=ADT_%A.out       # std out name
#SBATCH --error=ADT_%A.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=4         # num cpu cores / task
#SBATCH --mem=16G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short
## ---------------------------
##
## Script name: craft_seq_process_adt.sh
##
## Purpose of script: Wrapper script to process ADT modality of CRAFT-seq data
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



dir=$1 #dir=/path/to/ADT/KB
layout=$2 # 1. R1 files 2. R2 files: layout=/path/to/ADT/layout/ADT1.layout.tsv
name=$3 #name=ADT1

#load software/conda
eval "$(~/miniconda3/bin/conda shell.bash hook)"


#load files
WHITELIST="/path/to/384RNABarcodes.txt"
T2G="/path/to/FeaturesMismatch.t2g"
INDEX="/path/to/FeaturesMismatch.idx"


#1. merge fastq files
mkdir -p ${dir}/${name}
cd ${dir}/${name}

R1_files=`cat ${layout} | awk -F"\t" '{print $1}'`
R2_files=`cat ${layout} | awk -F"\t" '{print $2}'`

for file in ${R1_files}; do cat ${file}; done > ${name}.R1.fastq.gz
for file in ${R2_files}; do cat ${file}; done > ${name}.R2.fastq.gz

conda activate craft_seq
export JAVA_OPTS="-Xmx4g"
fastqc ${name}.R1.fastq.gz --outdir=./
fastqc ${name}.R2.fastq.gz --outdir=./

#2. Run Kallisto / bustools
conda activate kb_env2
OUT_DIR=${dir}/${name}/KBout
echo "kb count -i ${INDEX} -x 0,0,14:0,14,26:1,0,0 -t 4 -g ${T2G} -w ${WHITELIST} -o ${OUT_DIR} ${name}.R1.fastq.gz ${name}.R2.fastq.gz" >> ${dir}/${name}/Kallisto_cmd.txt
kb count -i ${INDEX} -x 0,0,14:0,14,26:1,0,0 -t 4 -g ${T2G} -w ${WHITELIST} -o ${OUT_DIR} ${name}.R1.fastq.gz ${name}.R2.fastq.gz


rm ${name}.R1.fastq.gz
rm ${name}.R2.fastq.gz