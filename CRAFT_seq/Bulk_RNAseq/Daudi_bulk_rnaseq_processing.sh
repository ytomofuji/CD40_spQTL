## ---------------------------
##
## Script name: Daudi_bulk_rnaseq_processing.sh
##
## Original file: 260421.bulk.Daudi.rna.seq.log.sh
##
## Purpose of script: Process Daudi cell bulk RNA-seq data
##
##
## Date Created: 2026-04-21
##
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


#1. make layout file for analysis
cd /path/to/output_directory

for ID in `seq 1 4`; do
for sample in stim nonstim; do
label=${sample}${ID}
R1_file=`ls *_${label}_*R1*.fastq.gz`
R2_file=`ls *_${label}_*R2*.fastq.gz`
echo -e "${label}\t/path/to/input_data/${R1_file}\t/path/to/input_data/${R2_file}"
done
done | sed '1i name\tR1_file\tR2_file' > layout.file.tsv



#2. run STAR mapping
mkdir -p /path/to/output_directory
cd /path/to/output_directory

awk 'NR>1 {print $0}' /path/to/layout.file.tsv | while read line; do
name=`echo -e ${line} | awk '{print $1}'`
R1_file=`echo -e ${line} | awk '{print $2}'`
R2_file=`echo -e ${line} | awk '{print $3}'`

sbatch /path/to/bulk_rnaseq_process_like_craft.sh \
/path/to/output_directory \
${R1_file} ${R2_file} ${name}
done




#3. merge into text file
mkdir -p /path/to/output_directory
cd /path/to/output_directory


Rscript -e '
library(tidyverse)
gene_file="/path/to/geneInfo.tab"
layout_file="/path/to/layout.file.tsv"

layout_df<-read_tsv(layout_file)
for(i in 1:nrow(layout_df)){
	sample_name<-layout_df$name[i]
	file_name<-str_c("/path/to/output_directory",sample_name,"/STARoutReadsPerGene.out.tab")
	tmp<-read_tsv(file_name,col_names=F)%>%dplyr::select(X1,X2)
	colnames(tmp)<-c("gene_id",sample_name)
	if(i==1){
		result<-tmp
	}else{
		result<-full_join(result,tmp,by="gene_id")
	}
}

gene_df<-read_tsv(gene_file,skip=1,col_names=F)
colnames(gene_df)<-c("gene_id","Gene","category")

gene_result<-inner_join(gene_df,result,by="gene_id")

write_tsv(gene_result,"Bulk.CD40L.Daudi.RNA.seq.count.annotation.tsv.gz")
write_tsv(result,"Bulk.CD40L.Daudi.RNA.seq.count.tsv.gz")
'




#4. Summarize fastQC
conda activate multiqc-env
cd /path/to/output_directory

BASE_DIR="/path/to/output_directory"
OUT_DIR="${BASE_DIR}/summary/raw.multiQC_results"

for i in 1 2; do
  TMP_DIR="/path/to/output_directory"
  mkdir -p ${TMP_DIR}

  echo "Finding FastQC files for R${i}..."
  find ${BASE_DIR} -type f -name "20260415_*_R${i}_*_fastqc.zip" -exec cp {} ${TMP_DIR}/ \;

  echo "Running MultiQC for R${i}..."
  multiqc ${TMP_DIR} \
    -n "multiqc_report_pilot_R${i}" \
    --outdir ${OUT_DIR} \
    --title "Pilot_Sequencing_R${i}"

  rm -r ${TMP_DIR}
done


OUT_DIR="${BASE_DIR}/summary/trimmed.multiQC_results"
for i in 1 2; do
  TMP_DIR="/path/to/output_directory"
  mkdir -p ${TMP_DIR}

  echo "Finding FastQC files for R${i}..."
  find ${BASE_DIR} -type f -name "trimmed*.R${i}_*fastqc.zip" -exec cp {} ${TMP_DIR}/ \;

  echo "Running MultiQC for R${i}..."
  multiqc ${TMP_DIR} \
    -n "multiqc_report_pilot_R${i}" \
    --outdir ${OUT_DIR} \
    --title "Pilot_Sequencing_R${i}"

  rm -r ${TMP_DIR}
done
