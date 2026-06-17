#ATAC-seq analysis @ 2 cpu

## ---------------------------
##
## Script name: ATAC_seq_mapping.sh
##
## Original file: BA1_mapping.GRCh37.sh
##
## Purpose of script: ATAC-seq alignment and peak calling pipeline
##
##
## Date Created: 2019-01-01
##
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

DIR=$1
FASTQ1=$2
FASTQ2=$3
NAME=$4
REF=$5
peakQC_REF=$6


eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate atac

mkdir -p ${DIR}
cd ${DIR}


#1. fastqc before QC
mkdir -p fastqc_preQC
fastqc --nogroup -t 2 -o ./fastqc_preQC ${FASTQ1}
fastqc --nogroup -t 2 -o ./fastqc_preQC ${FASTQ2}


#2. QC by fastp
#adaptor trimming + quality filtering + length filtering
mkdir -p ${DIR}/fastp
cd ${DIR}/fastp

mkdir -p fastp_html

fastp -i ${FASTQ1} -I ${FASTQ2} \
-o ${NAME}_R1.QCed.fastq.gz -O ${NAME}_R2.QCed.fastq.gz \
-3 -5 -W 4 -M 20 \
-w 2 \
-q 30 \
-n 5 \
-l 25 \
--html fastp_html/${NAME}.fastp.html


#3. fastqc after QC
mkdir -p ${DIR}/fastqc_postQC
cd ${DIR}
fastqc --nogroup -t 2 -o ./fastqc_postQC ./fastp/${NAME}_R1.QCed.fastq.gz
fastqc --nogroup -t 2 -o ./fastqc_postQC ./fastp/${NAME}_R2.QCed.fastq.gz


#4. mapping
mkdir -p ${DIR}/bam
cd ${DIR}/bam

bowtie2 -p 2 -x ${REF} \
-1 ${DIR}/fastp/${NAME}_R1.QCed.fastq.gz -2 ${DIR}/fastp/${NAME}_R2.QCed.fastq.gz -S ${NAME}.sam --maxins 2000

samtools view -@ 2 -h -Sb -F 1804 -f 2 -q 30 ${NAME}.sam | samtools sort -@ 2 > ${NAME}.tmp.bam
samtools index ${NAME}.tmp.bam
samtools view -@ 2 -h -b ${NAME}.tmp.bam {1..22} {X,Y} | samtools sort -@ 2 > ${NAME}.bam
samtools index ${NAME}.bam
rm ${NAME}.sam
rm ${NAME}.tmp.bam
rm ${NAME}.tmp.bam.bai


#5. duplicate removal
picard MarkDuplicates I=${NAME}.bam \
O=rm_dup_${NAME}.bam \
M=rm_dup_${NAME}.metrics.txt \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT

samtools index rm_dup_${NAME}.bam
rm ${NAME}.bam
rm ${NAME}.bam.bai



#7. Peak call with macs2
mkdir -p ${DIR}/peak
cd ${DIR}/peak

macs2 callpeak -t ${DIR}/bam/rm_dup_${NAME}.bam -f BAM -g 2.7e+9 --nomodel --nolambda --keep-dup all --call-summits \
--outdir ./ -n ${NAME} 2> ${NAME}.log

bedtools merge -i ${NAME}_peaks.narrowPeak | \
awk '($3-$2)<=3000 {print $0}' > ${NAME}.non.dup.bed


