## ---------------------------
##
## Script name: RCODE.CD40.eQTL.count.simulation.r
##
## Purpose of script: Generate CD40 mRNA pseudobulk counts with given allelic fold change based on NB model fitted to real data
##
##
## Date Created: 2026-09-30
##
## ---------------------------
##
## Notes:
##   NB model: count ~ dosage + covariates (tensorQTL eQTL mapping) + log10(library size)
##   log(mu) = eta - b_hat*(g-mean(g)) + beta_sim*(g-mean(g)), count ~ NB(mu, theta)
##   library size / covariates / genotype are kept as real data
##   fold change grid is defined per T allele (beta_sim = -log(FC_T), per C allele)
##
## ---------------------------

library(tidyverse)
library(MASS)

argv=commandArgs(T)

dataset_list=argv[1]
snp=argv[2]
n_rep=as.numeric(argv[3])
seed=as.numeric(argv[4])
out_file=argv[5]

if(FALSE){
dataset_list="CD40.B.eQTL.dataset.list.tsv"
snp="chr20:44746982:T:C_C"
n_rep=1000
seed=20261029
out_file="CD40.B.eQTL"
}

set.seed(seed)

FC_T<-unique(round(c(seq(1,0.90,by=-0.0125),seq(0.90,0.70,by=-0.025)),6))
FC_C<-seq(1,1.30,by=0.05)
beta_grid<-sort(unique(round(c(log(FC_C),-log(FC_T)),12)))

fix_id<-function(x,id_fix){
  if(id_fix %in% c("perez","onek1k")){x<-gsub("-","_",x)}
  if(id_fix=="perez"){x<-gsub("HC_","HC-",x)}
  return(x)
}

ds_df<-read_tsv(dataset_list)

D<-list()
param_df<-tibble()
for(i in 1:nrow(ds_df)){
d<-ds_df[i,]

raw_mx<-readRDS(d$raw)
colnames(raw_mx)<-fix_id(colnames(raw_mx),d$id_fix)

cov_mx<-read.table(d$cov,sep="\t",header=T,row.names=1,check.names=F)
ids<-colnames(cov_mx)
C<-t(as.matrix(cov_mx))

geno<-read_tsv(d$geno)
g<-geno[[snp]][match(ids,geno$IID)]

y<-as.numeric(raw_mx[d$feature,ids])
lib<-as.numeric(colSums(raw_mx[,ids]))

X<-cbind(dosage=g,C,log10_lib=log10(lib))
colnames(X)<-make.names(colnames(X),unique=T)
df<-as.data.frame(X)
df$y<-y
fit<-glm.nb(as.formula(str_c("y ~ ",str_c(colnames(X),collapse=" + "))),data=df)
b_hat<-coef(fit)["dosage"]
gc<-g-mean(g)
eta0<-fit$linear.predictors-b_hat*gc

D[[d$dataset]]<-list(ids=ids,g=g,lib=lib,eta0=eta0,gc=gc,theta=fit$theta)
param_df<-bind_rows(param_df,tibble(dataset=d$dataset,N=length(ids),ALT_freq=mean(g)/2,
mean_count=mean(y),var_count=var(y),median_lib=median(lib),theta=fit$theta,
BETA=b_hat,SE=summary(fit)$coefficients["dosage",2],P=summary(fit)$coefficients["dosage",4]))
}
write_tsv(param_df,str_c(out_file,".NB.parameters.tsv"))

#simulated counts: list(beta) -> list(dataset) -> matrix (sample x replicate)
sim<-list()
for(b in beta_grid){
  sim_b<-list()
  for(ds in names(D)){
    x<-D[[ds]]
    mu<-exp(x$eta0+b*x$gc)
    n<-length(mu)
    sim_b[[ds]]<-matrix(rnbinom(n*n_rep,mu=rep(mu,n_rep),size=x$theta),nrow=n,dimnames=list(x$ids,NULL))
  }
  sim[[as.character(b)]]<-sim_b
}
saveRDS(list(beta_grid=beta_grid,lib=lapply(D,`[[`,"lib"),counts=sim),str_c(out_file,".simulated.counts.rds"))
