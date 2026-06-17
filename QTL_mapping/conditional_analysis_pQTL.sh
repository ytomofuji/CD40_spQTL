#!/bin/bash
#SBATCH --job-name=conditional        # job name
#SBATCH --output=conditional.out       # std out name
#SBATCH --error=conditional.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=1         # num cpu cores / task
#SBATCH --mem=30G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short

module load singularity
export OMP_NUM_THREADS=1 

## ---------------------------
##
## Script name: conditional_analysis_pQTL.sh
##
## Original file: 260310.conditional.analysis.pQTL.all.v3.sh
##
## Purpose of script: Run conditional analysis for pQTL
##
##
## Date Created: 2026-02-13
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
target_gene=$2
condition_variant=$3
celltype=$4




eval "$(~/miniconda3/bin/conda shell.bash hook)"

mkdir -p ${dir}
cd ${dir}

export PATH=$PATH:/path/to/software/general/htslib-1.21/bin
export PATH=$PATH:/apps/released/gcc-toolchain/gcc-4.x/singularity/singularity-3.7.0/bin/
export PATH=$PATH:/path/to/software/genome/METAL/build/bin



#Perez-EUR
conda activate R
Rscript /path/to/code/B_cell/RCODE/RCODE.prepare.tensorqtl.input.for.conditional.analysis.r \
/path/to/analysis/B_eQTL/perez_sle/260305_tensorqtl_pPC/input/pPC2.gPC4.cg.${celltype}_CITE_Phenotype_eur.bed.gz \
/path/to/analysis/B_eQTL/perez_sle/260305_tensorqtl_pPC/input/pPC2.gPC4.cg.${celltype}_covariate_file_eur.tsv \
/path/to/data/Perez/Genotype/EUR \
${target_gene} ${condition_variant} ${target_gene}.${celltype}.Perez.EUR

bgzip -f ${target_gene}.${celltype}.Perez.EUR_CITE_Phenotype.bed


genotype=/path/to/data/Perez/Genotype/EUR
phenotype=${target_gene}.${celltype}.Perez.EUR_CITE_Phenotype.bed.gz
outname=${target_gene}.${celltype}.Perez.EUR
covariate=${target_gene}.${celltype}.Perez.EUR_covariate_file.tsv
maf=0.1
outname=${target_gene}.${celltype}.Perez.EUR

singularity exec --bind /data:/data /path/to/singularity/tensorqtl_latest.sif \
/bin/bash -c "python3 -m tensorqtl ${genotype} ${phenotype} ${outname} \
--covariates ${covariate} --maf_threshold ${maf} \
--mode cis_nominal"

eval "$(~/miniconda3/bin/conda shell.bash hook)"

file=`ls ${outname}.cis_qtl_pairs.chr* | sed 'sX.parquetXXg'`
echo -e ${file} | tr -d '\n' | python -c 'import sys; import pandas as pd; [pd.read_parquet(l+".parquet").to_csv(l+".tsv.gz",sep="\t",index=False) for l in sys.stdin]'
rm ${file}.parquet
file=`ls ${outname}.cis_qtl_pairs.*.tsv.gz`
mv ${file} ${outname}.tsv.gz




#Perez-ASI
conda activate R
Rscript /path/to/code/B_cell/RCODE/RCODE.prepare.tensorqtl.input.for.conditional.analysis.r \
/path/to/analysis/B_eQTL/perez_sle/260305_tensorqtl_pPC/input/pPC2.gPC4.cg.${celltype}_CITE_Phenotype_asi.bed.gz \
/path/to/analysis/B_eQTL/perez_sle/260305_tensorqtl_pPC/input/pPC2.gPC4.cg.${celltype}_covariate_file_asi.tsv \
/path/to/data/Perez/Genotype/ASI \
${target_gene} ${condition_variant} ${target_gene}.${celltype}.Perez.ASI

bgzip -f ${target_gene}.${celltype}.Perez.ASI_CITE_Phenotype.bed


genotype=/path/to/data/Perez/Genotype/ASI
phenotype=${target_gene}.${celltype}.Perez.ASI_CITE_Phenotype.bed.gz
outname=${target_gene}.${celltype}.Perez.ASI
covariate=${target_gene}.${celltype}.Perez.ASI_covariate_file.tsv
maf=0.1
outname=${target_gene}.${celltype}.Perez.ASI

singularity exec --bind /data:/data /path/to/singularity/tensorqtl_latest.sif \
/bin/bash -c "python3 -m tensorqtl ${genotype} ${phenotype} ${outname} \
--covariates ${covariate} --maf_threshold ${maf} \
--mode cis_nominal"

eval "$(~/miniconda3/bin/conda shell.bash hook)"

file=`ls ${outname}.cis_qtl_pairs.chr* | sed 'sX.parquetXXg'`
echo -e ${file} | tr -d '\n' | python -c 'import sys; import pandas as pd; [pd.read_parquet(l+".parquet").to_csv(l+".tsv.gz",sep="\t",index=False) for l in sys.stdin]'
rm ${file}.parquet
file=`ls ${outname}.cis_qtl_pairs.*.tsv.gz`
mv ${file} ${outname}.tsv.gz



if [ "$target_gene" = "FCGR2A_B_C" ]; then
    amp_target="FCGR2A"
else
    amp_target="${target_gene}"
fi



#AMP-RA
conda activate R
Rscript /path/to/code/B_cell/RCODE/RCODE.prepare.tensorqtl.input.for.conditional.analysis.r \
/path/to/analysis/B_eQTL/AMP_PBMC/260305_tensorqtl_pPC/input/pPC5.gPC5.cg.${celltype}_CITE_Phenotype.bed.gz \
/path/to/analysis/B_eQTL/AMP_PBMC/260305_tensorqtl_pPC/input/pPC5.gPC5.cg.${celltype}_covariate_file.tsv \
/path/to/data/AMP/Genotype/AMP \
${amp_target} ${condition_variant} ${target_gene}.${celltype}.AMP.RA

bgzip -f ${target_gene}.${celltype}.AMP.RA_CITE_Phenotype.bed


genotype=/path/to/data/AMP/Genotype/AMP
phenotype=${target_gene}.${celltype}.AMP.RA_CITE_Phenotype.bed.gz
outname=${target_gene}.${celltype}.AMP.RA
covariate=${target_gene}.${celltype}.AMP.RA_covariate_file.tsv
maf=0.05
outname=${target_gene}.${celltype}.AMP.RA

singularity exec --bind /data:/data /path/to/singularity/tensorqtl_latest.sif \
/bin/bash -c "python3 -m tensorqtl ${genotype} ${phenotype} ${outname} \
--covariates ${covariate} --maf_threshold ${maf} \
--mode cis_nominal"

eval "$(~/miniconda3/bin/conda shell.bash hook)"

file=`ls ${outname}.cis_qtl_pairs.chr* | sed 'sX.parquetXXg'`
echo -e ${file} | tr -d '\n' | python -c 'import sys; import pandas as pd; [pd.read_parquet(l+".parquet").to_csv(l+".tsv.gz",sep="\t",index=False) for l in sys.stdin]'
rm ${file}.parquet
file=`ls ${outname}.cis_qtl_pairs.*.tsv.gz`
mv ${file} ${outname}.tsv.gz






#AMP-SLE
conda activate R
Rscript /path/to/code/B_cell/RCODE/RCODE.prepare.tensorqtl.input.for.conditional.analysis.r \
/path/to/analysis/B_eQTL/AMP_PBMC_SLE/260305_CITE_pQTL/input/pPC5.gPC5.cg.${celltype}_CITE_Phenotype.bed.gz \
/path/to/analysis/B_eQTL/AMP_PBMC_SLE/260305_CITE_pQTL/input/pPC5.gPC5.cg.${celltype}_covariate_file.tsv \
/path/to/data/AMP/Genotype/AMP \
${amp_target} ${condition_variant} ${target_gene}.${celltype}.AMP.SLE

bgzip -f ${target_gene}.${celltype}.AMP.SLE_CITE_Phenotype.bed


genotype=/path/to/data/AMP/Genotype/AMP
phenotype=${target_gene}.${celltype}.AMP.SLE_CITE_Phenotype.bed.gz
outname=${target_gene}.${celltype}.AMP.SLE
covariate=${target_gene}.${celltype}.AMP.SLE_covariate_file.tsv
maf=0.05
outname=${target_gene}.${celltype}.AMP.SLE

singularity exec --bind /data:/data /path/to/singularity/tensorqtl_latest.sif \
/bin/bash -c "python3 -m tensorqtl ${genotype} ${phenotype} ${outname} \
--covariates ${covariate} --maf_threshold ${maf} \
--mode cis_nominal"

eval "$(~/miniconda3/bin/conda shell.bash hook)"

file=`ls ${outname}.cis_qtl_pairs.chr* | sed 'sX.parquetXXg'`
echo -e ${file} | tr -d '\n' | python -c 'import sys; import pandas as pd; [pd.read_parquet(l+".parquet").to_csv(l+".tsv.gz",sep="\t",index=False) for l in sys.stdin]'
rm ${file}.parquet
file=`ls ${outname}.cis_qtl_pairs.*.tsv.gz`
mv ${file} ${outname}.tsv.gz





#run meta-analysis
outname=${target_gene}.${celltype}
out_name=${outname}

declare -A datasets=(
    ["${outname}.Perez.EUR.tsv.gz"]="eur"
    ["${outname}.Perez.ASI.tsv.gz"]="asi"
)

layout_file="${out_name}.layout.txt"
echo -e "MARKER GENEREPLACESNP\nALLELE REF ALT\nEFFECT BETA\nPVALUE P\nSTDERR SE\n" > $layout_file
echo -e "SCHEME STDERR" >> $layout_file

valid_datasets=""

for file in "${!datasets[@]}"; do
    label=${datasets[$file]}
    
    if [ -f "$file" ]; then
        echo "Processing ${file}..."
        
        zcat "$file" | \
        sed -e 'sXCD45RA:HI100XCD45RAXg' \
            -e 's/CD45-prot/CD45/g' \
            -e 's/CD45RA-prot/CD45RA/g' \
            -e 's/CD45RO-prot/CD45RO/g' | \
        grep -v "INS" | \
        awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
        sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > "${label}.${out_name}.cis.pQTL.tsv"
        
        echo -e "PROCESS ${label}.${out_name}.cis.pQTL.tsv" >> $layout_file
        
        valid_datasets="$valid_datasets $label"
    else
        echo "Warning: ${file} not found. Skipping."
    fi
done





declare -A datasets=(
    ["${outname}.AMP.RA.tsv.gz"]="AMPRA"
    ["${outname}.AMP.SLE.tsv.gz"]="AMPSLE"
)

for file in "${!datasets[@]}"; do
    label=${datasets[$file]}
    
    if [ -f "$file" ]; then
        echo "Processing ${file}..."
        
        zcat "$file" | \
        sed -e 'sXCD45RA:HI100XCD45RAXg' \
            -e 's/FCGR2A/FCGR2A_B_C/g' \
            -e 's/CD45-prot/CD45/g' \
            -e 's/CD45RA-prot/CD45RA/g' \
            -e 's/CD45RO-prot/CD45RO/g' | \
        grep -v "INS" | \
        awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
        sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > "${label}.${out_name}.cis.pQTL.tsv"
        
        echo -e "PROCESS ${label}.${out_name}.cis.pQTL.tsv" >> $layout_file
        
        valid_datasets="$valid_datasets $label"
    else
        echo "Warning: ${file} not found. Skipping."
    fi
done



# METAL
echo -e "OUTFILE ivw.${out_name}. .cis.pQTL.tsv" >> $layout_file
echo -e "ANALYZE HETEROGENEITY" >> $layout_file

metal < $layout_file

if [ -f "ivw.${out_name}.1.cis.pQTL.tsv" ]; then
    sed 'sXREPLACEX\tXg' "ivw.${out_name}.1.cis.pQTL.tsv" | sed 'sXMarkerNameXGENE\tMarkerNameXg' | gzip -c > "ivw.${out_name}.cis.pQTL.tsv.gz"
    rm "ivw.${out_name}.1.cis.pQTL.tsv"
fi

for label in $valid_datasets; do
    sed 'sXREPLACEX\tXg' "${label}.${out_name}.cis.pQTL.tsv" | gzip -c > "${label}.${out_name}.cis.pQTL.tsv.gz"
    rm "${label}.${out_name}.cis.pQTL.tsv"
done

rm "$layout_file"


conda activate R

Rscript /path/to/code/B_cell/RCODE/RCODE.260217.Process.METAL.cond.EUR.ASI.RA.SLE.v2.r \
.${out_name}.cis.pQTL.tsv.gz ${out_name}

