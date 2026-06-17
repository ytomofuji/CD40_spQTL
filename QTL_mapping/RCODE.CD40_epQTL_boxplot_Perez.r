## ---------------------------
##
## Script name: RCODE.CD40_epQTL_boxplot_Perez.r
##
## Original file: RCODE.260126.CD40.var.e.pQTL.perez.data.ASI.EUR.mix.boxplot.log.r
##
## Purpose of script: Evaluate associations between CD40 variant and CITE-seq proteins (Perez dataset, ASI+EUR)
##
##
## Date Created: 2026-01-26
##
## Email:
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


setwd("/path/to/output_directory")

library(tidyverse)
library(ggrastr)
library(metafor)
library(ggrepel)



metafile="/path/to/analysis/B_eQTL/perez_sle/250203_pseudobulk_eQTL/input/log.cg.B.meta.data.tsv.gz"
exp_mx="/path/to/analysis/B_eQTL/perez_sle/250203_pseudobulk_eQTL/input/cg.B.pseudobulk.log.rds"

genotype1="/path/to/analysis/B_eQTL/perez_sle/cd40_variants/lat.RA.CD40.raw"
genotype2="/path/to/analysis/B_eQTL/perez_sle/cd40_variants/omni.RA.CD40.raw"

outname="Perez.B"



meta<-read_tsv(metafile)%>%filter(Count>=5)
geno_df1<-read_tsv(genotype1)%>%mutate(cohort=1)%>%dplyr::select(1,7,8)
geno_df2<-read_tsv(genotype2)%>%mutate(cohort=2)%>%dplyr::select(1,7,8)
colnames(geno_df1)<-c("individual","variant","cohort")
colnames(geno_df2)<-c("individual","variant","cohort")
geno_df<-bind_rows(geno_df1,geno_df2)

df<-inner_join(meta,geno_df,by=c("ind_cov"="individual"))
sample_list<-unique(df$ind_cov)



exp<-readRDS(exp_mx)
colnames(exp)<-gsub("-","_",colnames(exp))
colnames(exp)<-gsub("HC_","HC-",colnames(exp))
exp<-exp[,sample_list]

exp_df<-as.data.frame(t(as.matrix(exp)))%>%rownames_to_column("ind_cov")%>%as_tibble()%>%
pivot_longer(names_to="Gene",values_to="EXP",-ind_cov)%>%
filter(ind_cov %in% df$ind_cov)

gene="CD40"

n_fun<-function(x){
  return(data.frame(y=max(x)*1.1,label=paste0("N = ",length(x))))
}


plt_df<-filter(exp_df, Gene==gene)%>%inner_join(df, by="ind_cov")
plt_df$variant<-factor(plt_df$variant, levels=c(2,1,0))

p<-ggplot(plt_df,aes(x=variant,y=EXP,fill=variant))+
geom_boxplot(outlier.shape=NA,colour="black")+
geom_jitter_rast(alpha=0.2,size=0.1,height=0,width=0.3,colour="#000000",shape=21)+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=2/1,legend.position="none",
	theme(plot.title = element_text(size = 4,  hjust = 0.5)))+
scale_fill_manual(values = c("2" = "#1b9e77","1" = "#d95f02","0" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
xlab("")+ylab("log(CD40-RNA)")

ggsave(str_c(outname,".CD40.mRNA.exp.boxplot.pdf"),p,width=3,height=3)



p<-ggplot(filter(plt_df,pop_cov=="Asian"),aes(x=variant,y=EXP,fill=variant))+
geom_boxplot(outlier.shape=NA,colour="black")+
geom_jitter_rast(alpha=0.2,size=0.1,height=0,width=0.3,colour="#000000",shape=21)+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=2/1,legend.position="none",
	theme(plot.title = element_text(size = 4,  hjust = 0.5)))+
scale_fill_manual(values = c("2" = "#1b9e77","1" = "#d95f02","0" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
xlab("")+ylab("log(CD40-RNA)")

ggsave(str_c(outname,".CD40.mRNA.exp.boxplot.ASI.pdf"),p,width=3,height=3)


p<-ggplot(filter(plt_df,pop_cov=="European"),aes(x=variant,y=EXP,fill=variant))+
geom_boxplot(outlier.shape=NA,colour="black")+
geom_jitter_rast(alpha=0.2,size=0.1,height=0,width=0.3,colour="#000000",shape=21)+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=2/1,legend.position="none",
	theme(plot.title = element_text(size = 4,  hjust = 0.5)))+
scale_fill_manual(values = c("2" = "#1b9e77","1" = "#d95f02","0" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
xlab("")+ylab("log(CD40-RNA)")

ggsave(str_c(outname,".CD40.mRNA.exp.boxplot.EUR.pdf"),p,width=3,height=3)



metafile="/path/to/analysis/B_eQTL/perez_sle/260305_pQTL_clr_mean/pseudobulk_gex/cg.B.pseudobulk.metadata.tsv.gz"
exp_mx="/path/to/analysis/B_eQTL/perez_sle/260305_pQTL_clr_mean/pseudobulk_gex/cg.B.pseudobulk.mean.rds"

genotype1="/path/to/analysis/B_eQTL/perez_sle/cd40_variants/lat.RA.CD40.raw"
genotype2="/path/to/analysis/B_eQTL/perez_sle/cd40_variants/omni.RA.CD40.raw"

outname="Perez.B"



meta<-read_tsv(metafile)%>%filter(Count>=5)
geno_df1<-read_tsv(genotype1)%>%mutate(cohort=1)%>%dplyr::select(1,7,8)
geno_df2<-read_tsv(genotype2)%>%mutate(cohort=2)%>%dplyr::select(1,7,8)
colnames(geno_df1)<-c("individual","variant","cohort")
colnames(geno_df2)<-c("individual","variant","cohort")
geno_df<-bind_rows(geno_df1,geno_df2)

df<-inner_join(meta,geno_df,by=c("ind_cov"="individual"))
sample_list<-unique(df$ind_cov)



exp<-readRDS(exp_mx)
colnames(exp)<-gsub("-","_",colnames(exp))
colnames(exp)<-gsub("HC_","HC-",colnames(exp))
exp<-exp[,sample_list]

exp_df<-as.data.frame(t(as.matrix(exp)))%>%rownames_to_column("ind_cov")%>%as_tibble()%>%
pivot_longer(names_to="Gene",values_to="EXP",-ind_cov)%>%
filter(ind_cov %in% df$ind_cov)

gene="CD40"

n_fun<-function(x){
  return(data.frame(y=max(x)*1.1,label=paste0("N = ",length(x))))
}


plt_df<-filter(exp_df, Gene==gene)%>%inner_join(df, by="ind_cov")
plt_df$variant<-factor(plt_df$variant, levels=c(2,1,0))

p<-ggplot(plt_df,aes(x=variant,y=EXP,fill=variant))+
geom_boxplot(outlier.shape=NA,colour="black")+
geom_jitter_rast(alpha=0.2,size=0.1,height=0,width=0.3,colour="#000000",shape=21)+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=2/1,legend.position="none",
	theme(plot.title = element_text(size = 4,  hjust = 0.5)))+
scale_fill_manual(values = c("2" = "#1b9e77","1" = "#d95f02","0" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
xlab("")+ylab("clr(CD40-ADT)")

ggsave(str_c(outname,".CD40.ADT.exp.boxplot.pdf"),p,width=3,height=3)




