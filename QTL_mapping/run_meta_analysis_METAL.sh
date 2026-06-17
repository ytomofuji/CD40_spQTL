## ---------------------------
##
## Script name: run_meta_analysis_METAL.sh
##
## Original file: 250911.run.meta.analysis.with.METAL.log.sh
##
## Purpose of script: Run meta-analysis using METAL-software
##
##
## Date Created: 2025-09-11
##
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------


mkdir -p /path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/pQTL/log
cd /path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/pQTL/log




<<pQTL.meta.sh
#!/bin/bash
#SBATCH --job-name=QTL.meta        # job name
#SBATCH --output=QTL.meta.out       # std out name
#SBATCH --error=QTL.meta.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=1         # num cpu cores / task
#SBATCH --mem=100G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short

cd /path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/pQTL

export PATH=$PATH:/path/to/software/genome/METAL/build/bin
export OMP_NUM_THREADS=1 
eval "$(~/miniconda3/bin/conda shell.bash hook)"
conda activate R

cell=$1

asi_file="/path/to/analysis/B_eQTL/perez_sle/250802_tensorqtl_pPC/tensorqtl/pPC2.gPC4/pPC2.gPC4.cg.${cell}.asi.cis_qtl_pairs.tsv.gz"
eur_file="/path/to/analysis/B_eQTL/perez_sle/250802_tensorqtl_pPC/tensorqtl/pPC2.gPC4/pPC2.gPC4.cg.${cell}.eur.cis_qtl_pairs.tsv.gz"
ampra_file="/path/to/analysis/B_eQTL/AMP_PBMC/250802_tensorqtl_pPC/tensorqtl/pPC5.gPC5/pPC5.gPC5.cg.${cell}.cis_qtl_pairs.tsv.gz"
ampsle_file="/path/to/analysis/B_eQTL/AMP_PBMC_SLE/250901_CITE_pQTL/tensorqtl/pPC5.gPC5/pPC5.gPC5.cg.${cell}.cis_qtl_pairs.tsv.gz"
out_name=${cell}



#SNP, REF, ALT, BETA, SE, P, freq, mac
zcat ${asi_file} | sed 'sXCD45RA:HI100XCD45RAXg' | grep -v "INS" | awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > asi.${out_name}.cis.pQTL.tsv

zcat ${eur_file} | sed 'sXCD45RA:HI100XCD45RAXg' | grep -v "INS"  | awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > eur.${out_name}.cis.pQTL.tsv

zcat ${ampra_file} | sed -e 's/FCGR2A/FCGR2A_B_C/g' -e 's/CD45-prot/CD45/g' -e 's/CD45RA-prot/CD45RA/g' -e 's/CD45RO-prot/CD45RO/g' | grep -v "INS"  | awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > AMPRA.${out_name}.cis.pQTL.tsv

zcat ${ampsle_file} | sed 'sXFCGR2AXFCGR2A_B_CXg' | grep -v "INS"  | awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > AMPSLE.${out_name}.cis.pQTL.tsv


echo -e "MARKER GENEREPLACESNP\nALLELE REF ALT\nEFFECT BETA\nPVALUE P\nSTDERR SE\n" > ${out_name}.layout.txt
echo -e "SCHEME STDERR" >> ${out_name}.layout.txt
echo -e "PROCESS asi.${out_name}.cis.pQTL.tsv" >> ${out_name}.layout.txt
echo -e "PROCESS eur.${out_name}.cis.pQTL.tsv" >> ${out_name}.layout.txt
echo -e "PROCESS AMPRA.${out_name}.cis.pQTL.tsv" >> ${out_name}.layout.txt
echo -e "PROCESS AMPSLE.${out_name}.cis.pQTL.tsv" >> ${out_name}.layout.txt
echo -e "OUTFILE ivw.${out_name}. .cis.pQTL.tsv" >> ${out_name}.layout.txt
echo -e "ANALYZE HETEROGENEITY" >> ${out_name}.layout.txt

metal < ${out_name}.layout.txt

sed 'sXREPLACEX\tXg' ivw.${out_name}.1.cis.pQTL.tsv | sed 'sXMarkerNameXGENE\tMarkerNameXg' | gzip -c > ivw.${out_name}.cis.pQTL.tsv.gz
for dataset in asi eur AMPRA AMPSLE ; do
sed 'sXREPLACEX\tXg' ${dataset}.${out_name}.cis.pQTL.tsv | gzip -c > ${dataset}.${out_name}.cis.pQTL.tsv.gz
rm ${dataset}.${out_name}.cis.pQTL.tsv
done

rm ivw.${out_name}.1.cis.pQTL.tsv
rm ${out_name}.layout.txt

Rscript /path/to/code/B_cell/RCODE/RCODE.250910.Process.METAL.pQTL.r \
.${out_name}.cis.pQTL.tsv.gz ${out_name}
pQTL.meta.sh


for cell in T4 T8 NK B cM ncM cDC pDC
do
sbatch pQTL.meta.sh ${cell}
done





#also do same analysis for eQTL
mkdir -p /path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/eQTL/log
cd /path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/eQTL/log

for cell in T4 T8 NK B cM ncM cDC pDC
do
sbatch eQTL.meta.sh ${cell}
done


<<eQTL.meta.sh
#!/bin/bash
#SBATCH --job-name=QTL.meta        # job name
#SBATCH --output=QTL.meta.out       # std out name
#SBATCH --error=QTL.meta.err        # std err name
#SBATCH --time="4-23:00:00"           # time (HH:MM:SS) Normal: "4-23:00:00" , bigmem: "7-23:00:00"
#SBATCH --ntasks=1                # num task
#SBATCH --cpus-per-task=1         # num cpu cores / task
#SBATCH --mem=100G                 # mem 
#SBATCH --partition=normal      # normal/bigmem/long/short

cd /path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/eQTL

export PATH=$PATH:/path/to/software/genome/METAL/build/bin
export OMP_NUM_THREADS=1 
eval "$(~/miniconda3/bin/conda shell.bash hook)"
conda activate R

cell=$1


asi_file="/path/to/analysis/B_eQTL/perez_sle/250804_eQTL_tensorqtl/tensorqtl/ePC5.gPC4/ePC5.gPC4.cg.${cell}.asi.cis_qtl_pairs.tsv.gz"
eur_file="/path/to/analysis/B_eQTL/perez_sle/250804_eQTL_tensorqtl/tensorqtl/ePC5.gPC4/ePC5.gPC4.cg.${cell}.eur.cis_qtl_pairs.tsv.gz"
ampra_file="/path/to/analysis/B_eQTL/AMP_PBMC/250803_eQTL/tensorqtl/ePC5.gPC5/ePC5.gPC5.cg.${cell}.cis_qtl_pairs.tsv.gz"
ampsle_file="/path/to/analysis/B_eQTL/AMP_PBMC_SLE/250907_eQTL/tensorqtl/ePC5.gPC5/ePC5.gPC5.cg.${cell}.cis_qtl_pairs.tsv.gz"
onek1k_file="/path/to/analysis/B_eQTL/OneK1K_bulk_eQTL/250808_tensorqtl_eQTL_Perez_annotation/tensorqtl/ePC10.gPC6/ePC10.gPC6.cg.${cell}.cis_qtl_pairs.tsv.gz"
out_name=${cell}



#SNP, REF, ALT, BETA, SE, P, freq, mac
zcat ${asi_file} | grep -v "INS" | awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > asi.${out_name}.cis.eQTL.tsv

zcat ${eur_file} | grep -v "INS"  | awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > eur.${out_name}.cis.eQTL.tsv

zcat ${ampra_file} | grep -v "INS"  | awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > AMPRA.${out_name}.cis.eQTL.tsv

zcat ${ampsle_file} | grep -v "INS"  | awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > AMPSLE.${out_name}.cis.eQTL.tsv

zcat ${onek1k_file} | grep -v "INS"  | awk -F'[:\t]' 'NR>1 {print $1"REPLACE"$2":"$3":"$4":"$5"\t"$4"\t"$5"\t"$11"\t"$12"\t"$10"\t"$7"\t"$9}' | \
sed '1i GENEREPLACESNP\tREF\tALT\tBETA\tSE\tP\tfreq\tmac' > OneK1K.${out_name}.cis.eQTL.tsv


echo -e "MARKER GENEREPLACESNP\nALLELE REF ALT\nEFFECT BETA\nPVALUE P\nSTDERR SE\n" > ${out_name}.layout.txt
echo -e "SCHEME STDERR" >> ${out_name}.layout.txt
echo -e "PROCESS asi.${out_name}.cis.eQTL.tsv" >> ${out_name}.layout.txt
echo -e "PROCESS eur.${out_name}.cis.eQTL.tsv" >> ${out_name}.layout.txt
echo -e "PROCESS AMPRA.${out_name}.cis.eQTL.tsv" >> ${out_name}.layout.txt
echo -e "PROCESS AMPSLE.${out_name}.cis.eQTL.tsv" >> ${out_name}.layout.txt
echo -e "PROCESS OneK1K.${out_name}.cis.eQTL.tsv" >> ${out_name}.layout.txt
echo -e "OUTFILE ivw.${out_name}. .cis.eQTL.tsv" >> ${out_name}.layout.txt
echo -e "ANALYZE HETEROGENEITY" >> ${out_name}.layout.txt

metal < ${out_name}.layout.txt


sed 'sXREPLACEX\tXg' ivw.${out_name}.1.cis.eQTL.tsv | sed 'sXMarkerNameXGENE\tMarkerNameXg' | gzip -c > ivw.${out_name}.cis.eQTL.tsv.gz
for dataset in asi eur AMPRA AMPSLE OneK1K; do
sed 'sXREPLACEX\tXg' ${dataset}.${out_name}.cis.eQTL.tsv | gzip -c > ${dataset}.${out_name}.cis.eQTL.tsv.gz
rm ${dataset}.${out_name}.cis.eQTL.tsv
done

rm ivw.${out_name}.1.cis.eQTL.tsv
rm ${out_name}.layout.txt

Rscript /path/to/code/B_cell/RCODE/RCODE.250911.Process.METAL.eQTL.r \
.${out_name}.cis.eQTL.tsv.gz ${out_name}

eQTL.meta.sh
