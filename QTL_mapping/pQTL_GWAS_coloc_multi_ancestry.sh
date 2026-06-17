#!/bin/bash
#SBATCH --job-name=coloc         # job name
#SBATCH --output=coloc.out       # std out name
#SBATCH --error=coloc.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=1         # num cpu cores / task
#SBATCH --mem=32G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short

export OMP_NUM_THREADS=1 
## ---------------------------
##
## Script name: pQTL_GWAS_coloc_multi_ancestry.sh
##
## Original file: 260311_multi.anc.pqtl.gwas.coloc.re.run.sh
##
## Purpose of script: Run colocalization analysis between GWAS and pQTL using multi-ancestry summary statistics
##
##
## Date Created: 2026-03-11
##
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

dir=$1 #/path/to/analysis/B_eQTL/AMP_PBMC/250727_tensorqtl/coloc
pheno_name=$2
gwas_sumstats=$3



#Run colocalization analysis (multi-ancestry pQTL)
mkdir -p ${dir}/${pheno_name}
cd ${dir}/${pheno_name}

eval "$(~/miniconda3/bin/conda shell.bash hook)"
conda activate R

export PATH="$PATH:/path/to/software/genome"

awk -F'\t' 'NR>1 {print $1"\t"$35"\t"$2}' /path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/summary/pQTL.meta.analysis.lead.variants.FDR.0.05.tsv | sed 's/:/\t/g' | \
while read line
do
gene=`echo ${line} | awk '{print $1}'`
cell=`echo ${line} | awk '{print $2}'`
chr=`echo ${line} | awk '{print $3}'`
pos=`echo ${line} | awk '{print $4}'`
ref=`echo ${line} | awk '{print $5}'`
alt=`echo ${line} | awk '{print $6}'`

if [ ! -e ${pheno_name}.${gene}.${cell}.all.coloc.txt ]; then

zcat /path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/pQTL/${cell}.METAL.result.processed.pQTL.tsv.gz | \
awk '$1=="'$gene'" {print $2"\t"$7"\t"$8"\t"$9}' | sed 's/:/\t/g' | sed '1i chr\tpos\tref\talt\tbeta\tse\tp' | gzip -c > ${gene}.${cell}.meta.analysis.sumstats.tsv.gz

chr_num=`echo ${chr} | sed 's/chr//g'`
var_name=`awk '$1=="'$chr_num'" && $4=='$pos' && $6=="'$ref'" && $5=="'$alt'" {print $2}' /path/to/1KG_p3v5.${chr_num}`


#export LD info for locus zoom plot
plink2 --threads 1 --memory 8000 \
--bfile /path/to/1KG_p3v5.${chr_num} \
--r2-unphased --ld-window-r2 0.1 --ld-snp ${var_name} --out ${gene}.${cell}.all --keep-allele-order


Rscript /path/to/code/B_cell/RCODE/RCODE.CITE.seq.pqtl.curated.sumstats.coloc.r \
${gwas_sumstats} \
${gene}.${cell}.meta.analysis.sumstats.tsv.gz \
${gene}.${cell}.all.vcor \
${pos} ${chr} ${gene} ${cell} all ${pheno_name}
fi
done

