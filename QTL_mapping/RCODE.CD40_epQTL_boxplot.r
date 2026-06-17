## ---------------------------
##
## Script name: RCODE.CD40_epQTL_boxplot.r
##
## Original file: RCODE.260126.CD40.var.e.pQTL.data.boxplot.r
##
## Purpose of script: Evaluate associations between CD40 variant and CITE-seq proteins (AMP-RA dataset)
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




library(tidyverse)
library(ggrastr)
library(metafor)
library(ggrepel)


argv=commandArgs(T)
metafile=argv[1]
exp_mx=argv[2]
genotype=argv[3]
outname=argv[4]
Phenotype=argv[5]
ind_column=argv[6]




meta<-read_tsv(metafile)%>%filter(Count>=5)
geno_df<-read_tsv(genotype)%>%dplyr::select(2,7)
colnames(geno_df)<-c("individual","variant")

sample_list<-pull(meta,ind_column)



exp<-readRDS(exp_mx)
colnames(exp)<-gsub("-","_",colnames(exp))
sample_list<-gsub("-","_",sample_list)
geno_df$individual<-gsub("-","_",geno_df$individual)
sample_list<-sample_list[sample_list%in%geno_df$individual]

exp<-exp[,sample_list]

exp_df<-as.data.frame(t(as.matrix(exp)))%>%rownames_to_column("individual")%>%as_tibble()%>%
pivot_longer(names_to="Gene",values_to="EXP",-individual)%>%
filter(individual %in% sample_list)



n_fun<-function(x){
  return(data.frame(y=max(x)*1.1,label=paste0("N = ",length(x))))
}


plt_df<-filter(exp_df, Gene==Phenotype)%>%inner_join(geno_df, by="individual")
plt_df$variant<-factor(plt_df$variant, levels=c(2,1,0))

p<-ggplot(plt_df,aes(x=variant,y=EXP,fill=variant))+
geom_boxplot(outlier.shape=NA,colour="black")+
geom_jitter_rast(alpha=0.2,size=0.1,height=0,width=0.3,colour="#000000",shape=21)+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=2/1,legend.position="none",
	plot.title = element_text(size = 4,  hjust = 0.5))+
scale_fill_manual(values = c("2" = "#1b9e77","1" = "#d95f02","0" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
xlab("")+ylab("log(CD40-RNA)")

ggsave(str_c(outname,".CD40.exp.boxplot.pdf"),p,width=3,height=3)





