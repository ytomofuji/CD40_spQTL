#!/bin/bash
#SBATCH --job-name=gDNA         # job name
#SBATCH --output=gDNA_%A.out       # std out name
#SBATCH --error=gDNA_%A.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=4         # num cpu cores / task
#SBATCH --mem=16G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short
## ---------------------------
##
## Script name: craft_seq_process_gDNA.sh
##
## Purpose of script: Wrapper script to process gDNA modality of CRAFT-seq data
##
##
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------
# Script from https://github.com/immunogenomics/craft-seq/blob/main/preprocessing/DNAStep1_Demultiplex.sh


dir=$1 #dir=/path/to/gDNA/CRISPRESSO
layout=$2 # 1. R1 files 2. R2 files: layout=/path/to/gDNA/layout/gDNA1.layout.tsv
name=$3 #name=gDNA1
SEQ=$4 #SEQ=tcccgtctggtcggggcaaaaacaactcacagcggtcagcaagcagccccagaggacgcactgcagaggcagacgaaccatggcgaggtgagaccaggcggcaggaccactgggcgcccgagcgaggcctctgctgactcccctgcccca

samples="/path/to/CRISPR_DNABarcodes.txt"
BATCH_BASE="/path/to/CRISPRessoBatch_R2.txt"


#1. merge fastq files
mkdir -p ${dir}/${name}
cd ${dir}/${name}

R1_files=`cat ${layout} | awk -F"\t" '{print $1}'`
R2_files=`cat ${layout} | awk -F"\t" '{print $2}'`

for file in ${R1_files}; do cat ${file}; done > ${name}.R1.fastq.gz
for file in ${R2_files}; do cat ${file}; done > ${name}.R2.fastq.gz

eval "$(~/miniconda3/bin/conda shell.bash hook)"
conda activate craft_seq
export JAVA_OPTS="-Xmx4g"
fastqc ${name}.R1.fastq.gz --outdir=./
fastqc ${name}.R2.fastq.gz --outdir=./

#2. de-multiplexing
mkdir -p Demultiplex

paste <(zcat ${name}.R1.fastq.gz | paste - - - - ) <(zcat ${name}.R2.fastq.gz | paste - - - - ) | \
awk -v FS="\t" -v OFS="\n" -v n="${name}" '
FNR==NR {samples[$2]=$2; next} 
{
    barcode = substr($2,1,10); 
    if (samples[barcode]) {
        print $1,$2,$3,$4 >> "./Demultiplex/" samples[barcode] "_R1.fastq"; 
        print $5,$6,$7,$8 >> "./Demultiplex/" samples[barcode] "_R2.fastq"
    }
}' $samples -

#gzip them all back up to save space
gzip ./Demultiplex/*.fastq




#3. Run CRISPRESSO
# activate crispresso env
# export LD_LIBRARY_PATH=/path/to/miniconda3/lib:$LD_LIBRARY_PATH

PLATE_BATCH="${dir}/${name}/batch.batch"
head -1 ${BATCH_BASE} > ${PLATE_BATCH}

### Subset batch files to include only the files present in directory
ls ${dir}/${name}/Demultiplex/*_R2.fastq.gz >> dna_filepaths.txt
sed '2,$s/.fastq/.fastq.gz/g' ${BATCH_BASE} | awk -v FASTQ_DIR="${dir}/${name}/Demultiplex/" -v OFS="\t" '$2=FASTQ_DIR$2' | awk 'NR!=1 {print}' | grep -f dna_filepaths.txt >> ${PLATE_BATCH}
rm dna_filepaths.txt

CRISPRessoBatch -o ${dir}/${name} --batch_settings ${PLATE_BATCH}  --skip_failed -p 4 -a ${SEQ} -bo ${dir}/${name} --suppress_plots

# Unzip all Alleles_frequency_tables for easy reading into R
find . -type f -name "Alleles_frequency_table.zip" -exec sh -c 'unzip -qq -d "${1%/*}" "$1"' _ {} \;



#4. Summarize CRISPRESSO result into a single text file (matrix file)
conda activate R
Rscript $(dirname "$0")/RCODE.process.crispresso.r ${dir}/${name} ${name}



#5. remove/compress intermediate files
tar -czf CRISPRessoBatch_on_batch.tar.gz CRISPRessoBatch_on_batch && rm -r CRISPRessoBatch_on_batch
tar -czf Demultiplex.tar.gz Demultiplex && rm -r Demultiplex
rm ${name}.R1.fastq.gz
rm ${name}.R2.fastq.gz