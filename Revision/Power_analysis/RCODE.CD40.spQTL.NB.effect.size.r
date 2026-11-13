## ---------------------------
##
## Script name: RCODE.CD40.spQTL.NB.effect.size.r
##
## Purpose of script: Estimate CD40 spQTL effect (ln fold change per allele) with negative binomial model on raw pseudobulk ADT counts
##
##
## Date Created: 2026-09-29
##
## ---------------------------
##
## Notes:
##   covariates = same as tensorQTL pQTL mapping + log10(total ADT UMI)
##   dataset_list: dataset, raw (pseudobulk raw count rds), cov (tensorQTL covariate file), geno (plink2 --export A), feature, id_fix
##
## ---------------------------

library(tidyverse)
library(MASS)
library(metafor)

argv=commandArgs(T)

dataset_list=argv[1]
snp=argv[2]
out_file=argv[3]

if(FALSE){
dataset_list="CD40.B.spQTL.dataset.list.tsv"
snp="chr20:44746982:T:C_C"
out_file="CD40.B.spQTL"
}

fix_id<-function(x,id_fix){
  if(id_fix %in% c("perez","onek1k")){x<-gsub("-","_",x)}
  if(id_fix=="perez"){x<-gsub("HC_","HC-",x)}
  return(x)
}

ds_df<-read_tsv(dataset_list)

res<-tibble()
for(i in 1:nrow(ds_df)){
d<-ds_df[i,]

raw_mx<-readRDS(d$raw)
colnames(raw_mx)<-fix_id(colnames(raw_mx),d$id_fix)
cnt_df<-tibble(ID=colnames(raw_mx),count=as.numeric(raw_mx[d$feature,]),total_UMI=as.numeric(colSums(raw_mx)))

cov_mx<-read.table(d$cov,sep="\t",header=T,row.names=1,check.names=F)
rownames(cov_mx)<-make.names(rownames(cov_mx))
cov_df<-as.data.frame(t(cov_mx))%>%rownames_to_column("ID")%>%as_tibble()

geno_df<-read_tsv(d$geno)%>%dplyr::select(IID,all_of(snp))
colnames(geno_df)<-c("ID","dosage")

df<-cov_df%>%inner_join(cnt_df,by="ID")%>%inner_join(geno_df,by="ID")%>%mutate(log10_total_UMI=log10(total_UMI))

f_nb<-as.formula(str_c("count ~ dosage + ",str_c(c(rownames(cov_mx),"log10_total_UMI"),collapse=" + ")))
fit_nb<-glm.nb(f_nb,data=df)
s_nb<-summary(fit_nb)$coefficients["dosage",]

res<-bind_rows(res,tibble(dataset=d$dataset,N=nrow(df),ALT_freq=mean(df$dosage)/2,
median_count=median(df$count),median_total_UMI=median(df$total_UMI),
BETA=s_nb[1],SE=s_nb[2],P=s_nb[4],theta=fit_nb$theta,FC=exp(s_nb[1])))
}
write_tsv(res,str_c(out_file,".NB.per.dataset.tsv"))

#fixed effect meta-analysis
fit_meta<-rma(yi=res$BETA,sei=res$SE,method="FE")
meta_df<-tibble(k=nrow(res),BETA=as.numeric(fit_meta$beta),SE=fit_meta$se,P=fit_meta$pval,
lower=fit_meta$ci.lb,upper=fit_meta$ci.ub,I2=fit_meta$I2,Q=fit_meta$QE,Qp=fit_meta$QEp)%>%
mutate(FC=exp(BETA),FC_lower=exp(lower),FC_upper=exp(upper))
write_tsv(meta_df,str_c(out_file,".NB.fe.meta.analysis.tsv"))
