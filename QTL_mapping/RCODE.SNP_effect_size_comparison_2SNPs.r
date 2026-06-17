## ---------------------------
##
## Script name: RCODE.SNP_effect_size_comparison_2SNPs.r
##
## Original file: RCODE.260112.SNP.effect.size.comparison.across.celltypes.2SNPs.for.Figure2.r
##
## Purpose of script: Make forest plot to compare QTL effect sizes across cell types (two SNPs per gene)
##
##
## Date Created: 2026-01-12
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
library(ggrepel)


argv=commandArgs(T)
input_summary=argv[1]
snp1=argv[2]
snp2=argv[3]
gene=argv[4]


SNP_color<-c("#0F8EAE","#D24655")

data<-read_tsv(input_summary)
df<-filter(data,Gene==gene,SNP%in%c(snp1,snp2))
df$SNP<-factor(df$SNP, levels=rev(c(snp1,snp2)))
df$cell_type<-factor(df$cell_type, levels=rev(c("T4","T8","NK","B","PB","cM","ncM","cDC","pDC")))


#make forest plot
p<-ggplot(df,aes(y=cell_type,x=BETA,color=SNP))+
geom_point(position=position_dodge(width = 0.3))+geom_errorbar(aes(xmin=BETA-SE*1.96,xmax=BETA+SE*1.96),width=0,position=position_dodge(width = 0.3))+
geom_text_repel(data=filter(df,P<0.05),aes(label=str_c("P=",signif(P,digits=2))))+
scale_color_manual(values=SNP_color)+
theme_test()+theme(aspect.ratio=2/1, legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

name<-str_c(c(gene,".effect.size.forest.pdf"),collapse=".")
name<-gsub(":","_",name)
ggsave(name,p,width=2,height=4)

p<-ggplot(df,aes(y=cell_type,x=BETA,color=SNP))+
geom_point(position=position_dodge(width = 0.3))+geom_errorbar(aes(xmin=BETA-SE*1.96,xmax=BETA+SE*1.96),width=0,position=position_dodge(width = 0.3))+
scale_color_manual(values=SNP_color)+
theme_test()+theme(aspect.ratio=2/1, legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

name<-str_c(c(gene,".effect.size.forest.noP.pdf"),collapse=".")
ggsave(name,p,width=2,height=4)



tsvname<-str_c(c(gene,".effect.size.forest.ori.tsv"),collapse=".")
tsvname<-gsub(":","_",tsvname)
write_tsv(df,tsvname)