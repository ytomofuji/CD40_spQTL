plink2=/path/to/plink2
snp="chr20:44746982:T:C"

#genotype (C allele dosage)
mkdir -p genotype
echo -e "${snp}\tC" > genotype/count_allele.txt
${plink2} --bfile /path/to/ASI --snp ${snp} --export A --export-allele genotype/count_allele.txt --out genotype/Perez.ASI.lead
${plink2} --bfile /path/to/EUR --snp ${snp} --export A --export-allele genotype/count_allele.txt --out genotype/Perez.EUR.lead
${plink2} --bfile /path/to/all.EUR --snp ${snp} --export A --export-allele genotype/count_allele.txt --out genotype/Perez.allEUR.lead
${plink2} --bfile /path/to/AMP --snp ${snp} --export A --export-allele genotype/count_allele.txt --out genotype/AMP.lead
${plink2} --bfile /path/to/onek1k.allele.aligned.QCed_MAF5.all_chr --snp ${snp} --export A --export-allele genotype/count_allele.txt --out genotype/OneK1K.lead

#spQTL effect size (NB)
Rscript RCODE.CD40.spQTL.NB.effect.size.r CD40.B.spQTL.dataset.list.tsv ${snp}_C CD40.B.spQTL

#count generation
Rscript RCODE.CD40.eQTL.count.simulation.r CD40.B.eQTL.dataset.list.tsv ${snp}_C 1000 20261029 CD40.B.eQTL

#regression analysis (6480 = number of CD40 cis-variants in eQTL meta-analysis)
Rscript RCODE.CD40.eQTL.simulation.regression.r CD40.B.eQTL.dataset.list.tsv CD40.B.eQTL.simulated.counts.rds ${snp}_C 6480 CD40.B.eQTL
