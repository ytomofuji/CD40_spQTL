## ---------------------------
##
## Script name: RCODE.CD40.eQTL.simulation.regression.r
##
## Purpose of script: eQTL analysis of simulated CD40 counts with the same pipeline as tensorQTL eQTL mapping and power calculation
##
##
## Date Created: 2026-09-30
##
## ---------------------------
##
## Notes:
##   input = output of RCODE.CD40.eQTL.count.simulation.r
##   log1p(count/library size*1e4) -> INT (Perez ASI and EUR together) -> linear regression with covariates (tensorQTL-like)
##   -> fixed effect meta-analysis (IVW)
##   threshold: 0.05/number of cis-variants
##
## ---------------------------

library(tidyverse)

argv=commandArgs(T)

dataset_list=argv[1]
sim_file=argv[2]
snp=argv[3]
n_cis=as.numeric(argv[4])
out_file=argv[5]

if(FALSE){
dataset_list="CD40.B.eQTL.dataset.list.tsv"
sim_file="CD40.B.eQTL.simulated.counts.rds"
snp="chr20:44746982:T:C_C"
n_cis=6480
out_file="CD40.B.eQTL"
}

apply_INT <- function(x) {
  ranks <- rank(x,ties.method = "average")
  qnorm((ranks - 0.5) / length(x))
}

#residualize phenotype and genotype on covariates (same as tensorQTL)
assoc<-function(Y,g,Q,dof){
  rY<-Y-Q%*%crossprod(Q,Y)
  rg<-as.vector(g-Q%*%crossprod(Q,g))
  sgg<-sum(rg^2)
  beta<-as.vector(crossprod(rg,rY))/sgg
  rss<-colSums(rY^2)-beta^2*sgg
  se<-sqrt(rss/dof/sgg)
  p<-2*pt(-abs(beta/se),dof)
  return(list(beta=beta,se=se,p=p))
}

ds_df<-read_tsv(dataset_list)
sim<-readRDS(sim_file)

D<-list()
for(i in 1:nrow(ds_df)){
d<-ds_df[i,]
cov_mx<-read.table(d$cov,sep="\t",header=T,row.names=1,check.names=F)
ids<-colnames(cov_mx)
C<-t(as.matrix(cov_mx))
geno<-read_tsv(d$geno)
g<-geno[[snp]][match(ids,geno$IID)]
qrC<-qr(cbind(1,C))
Q<-qr.Q(qrC)[,seq_len(qrC$rank),drop=F]
D[[d$dataset]]<-list(g=g,Q=Q,dof=length(ids)-2-ncol(C),int_group=d$int_group)
}

res<-tibble()
for(b in sim$beta_grid){
  Y<-sim$counts[[as.character(b)]]
  I<-list()
  for(grp in unique(ds_df$int_group)){
    dss<-filter(ds_df,int_group==grp)%>%pull(dataset)
    L<-do.call(rbind,lapply(dss,function(ds) log1p(Y[[ds]]/sim$lib[[ds]]*1e4)))
    L<-apply(L,2,apply_INT)
    idx<-rep(dss,sapply(dss,function(ds) nrow(Y[[ds]])))
    for(ds in dss){I[[ds]]<-L[idx==ds,,drop=F]}
  }
  a<-lapply(names(D),function(ds) assoc(I[[ds]],D[[ds]]$g,D[[ds]]$Q,D[[ds]]$dof))
  names(a)<-names(D)
  B<-do.call(rbind,lapply(a,`[[`,"beta"))
  S<-do.call(rbind,lapply(a,`[[`,"se"))
  P<-do.call(rbind,lapply(a,`[[`,"p"))
  W<-1/S^2
  meta_beta<-colSums(W*B)/colSums(W)
  meta_se<-sqrt(1/colSums(W))
  meta_p<-2*pnorm(-abs(meta_beta/meta_se))
  res_b<-bind_rows(
    tibble(dataset=rep(names(D),each=ncol(B)),rep=rep(1:ncol(B),length(D)),BETA=as.vector(t(B)),SE=as.vector(t(S)),P=as.vector(t(P))),
    tibble(dataset="Meta.FE",rep=1:ncol(B),BETA=meta_beta,SE=meta_se,P=meta_p))%>%
  mutate(beta_sim=b,FC_C=exp(b),FC_T=exp(-b))
  res<-bind_rows(res,res_b)
}
write_tsv(res,str_c(out_file,".simulation.regression.all.reps.tsv.gz"))

power_df<-res%>%group_by(beta_sim,FC_C,FC_T,dataset)%>%
summarize(n_rep=n(),mean_BETA=mean(BETA),power_p0.05=mean(P<0.05),power_p1e6=mean(P<1e-6),power_p5e8=mean(P<5e-8),
power_bonferroni=mean(P<0.05/n_cis))%>%ungroup()
write_tsv(power_df,str_c(out_file,".simulation.power.summary.tsv"))
