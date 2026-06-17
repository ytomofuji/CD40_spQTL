#!/bin/bash
#SBATCH --job-name=tensorqtl        # job name
#SBATCH --output=tensorqtl.out       # std out name
#SBATCH --error=tensorqtl.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=1         # num cpu cores / task
#SBATCH --mem=30G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short
## ---------------------------
##
## Script name: run_tensorqtl.sh
##
## Original file: Shell_run_tensorqtl.eris2.chr.noFDRspline.sh
##
## Purpose of script: Run tensorqtl for eqtl mapping
##
##
## Date Created: 2026-03-05
##
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

dir=$1
genotype=$2
phenotype=$3
covariate=$4
outname=$5
maf=$6



module load singularity
export OMP_NUM_THREADS=1 


mkdir -p ${dir}
cd ${dir}

singularity exec --bind /data:/data /path/to/singularity/tensorqtl_latest.sif \
/bin/bash -c "python3 -m tensorqtl ${genotype} ${phenotype} ${outname} \
--covariates ${covariate} --maf_threshold ${maf} \
--mode cis --qvalue_lambda 0"

singularity exec --bind /data:/data /path/to/singularity/tensorqtl_latest.sif \
/bin/bash -c "python3 -m tensorqtl ${genotype} ${phenotype} ${outname} \
--covariates ${covariate} --maf_threshold ${maf} \
--mode cis_nominal"

singularity exec --bind /data:/data /path/to/singularity/tensorqtl_latest.sif \
/bin/bash -c "python3 -m tensorqtl ${genotype} ${phenotype} ${outname} \
--covariates ${covariate} --maf_threshold ${maf} --cis_output ${dir}/${outname}.cis_qtl.txt.gz \
--mode cis_independent --qvalue_lambda 0"



#convert cis-nominal results to tsv
eval "$(~/miniconda3/bin/conda shell.bash hook)"

for i in `seq 1 22`; do
echo -e "${outname}.cis_qtl_pairs.chr${i}" | tr -d '\n' | \
python -c 'import sys; import pandas as pd; [pd.read_parquet(l+".parquet").to_csv(l+".tsv.gz",sep="\t",index=False) for l in sys.stdin]'
rm ${outname}.cis_qtl_pairs.chr${i}.parquet
done

zcat ${outname}.cis_qtl_pairs.chr1.tsv.gz > ${outname}.cis_qtl_pairs.tsv
rm ${outname}.cis_qtl_pairs.chr1.tsv.gz

for i in `seq 2 22`; do
zcat ${outname}.cis_qtl_pairs.chr${i}.tsv.gz | awk 'NR>1 {print $0}' \
>> ${outname}.cis_qtl_pairs.tsv
rm ${outname}.cis_qtl_pairs.chr${i}.tsv.gz
done
gzip -f ${outname}.cis_qtl_pairs.tsv