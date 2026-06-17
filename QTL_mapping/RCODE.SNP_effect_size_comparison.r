## ---------------------------
##
## Script name: RCODE.SNP_effect_size_comparison.r
##
## Original file: RCODE.260112.SNP.effect.size.comparison.across.celltypes.for.Figure2.r
##
## Purpose of script: Make forest plot to compare QTL effect sizes across cell types
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
library(RColorBrewer)

argv=commandArgs(T)
input_summary=argv[1]
snp=argv[2]
gene=argv[3]


cell_color<-c("T4"="#FC7A62","T8"="#EC008C","NK"="#FF9C42","B"="#5D4AAB","PB"="#882255",
                            "cM"="#34BECC","ncM"="#CCECCC","cDC"="#409188","pDC"="#0077b6")

data<-read_tsv(input_summary)
df<-filter(data,Gene==gene,SNP==snp)
df$cell_type<-factor(df$cell_type, levels=rev(c("T4","T8","NK","B","PB","cM","ncM","cDC","pDC")))

#make forest plot
p<-ggplot(df,aes(y=cell_type,x=BETA,color=cell_type))+
geom_point()+geom_errorbar(aes(xmin=BETA-SE*1.96,xmax=BETA+SE*1.96),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(P,digits=2))))+
scale_color_manual(values=cell_color)+
theme_test()+theme(aspect.ratio=2/1, legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

name<-str_c(c(snp,gene,".effect.size.forest.pdf"),collapse=".")
name<-gsub(":","_",name)
ggsave(name,p,width=2,height=4)

tsvname<-str_c(c(snp,gene,".effect.size.forest.ori.tsv"),collapse=".")
tsvname<-gsub(":","_",tsvname)
write_tsv(df,tsvname)