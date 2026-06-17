## ---------------------------
##
## Script name: RCODE.RNA_ADT_QC.r
##
## Original file: 260316_RNA_ADT_QC.log.r
##
## Purpose of script: Process cDNA/ADT data from CRAFT-seq 
##
##
## Date Created: 2025-08-06
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
library(Seurat)
library(harmony)
library(ggpointdensity)
library(viridis)


rna_thr=500
adt_thr=30
mt_thr=20

plt_df<-read_tsv("/path/to/raw.metrics.df.tsv.gz")
plt_df$plate<-factor(plt_df$plate, levels=unique(plt_df$plate))

p<-ggplot(plt_df,aes(x=nCount_RNA,y=nFeature_RNA,color=plate))+
geom_point(size=0.1)+theme_test()+theme(aspect.ratio=1/1)+
scale_color_viridis(option="plasma",discrete=TRUE)+scale_x_log10()+scale_y_log10()+geom_hline(yintercept=rna_thr, linetype="dashed")
ggsave("nCount.Feature_RNA.scatter.pdf",p,width=4,height=3)

p<-ggplot(plt_df,aes(x=nFeature_RNA,fill=plate))+
geom_histogram(position = "stack")+theme_bw()+theme(aspect.ratio=1/1)+
scale_fill_viridis(option="plasma",discrete=TRUE)+
scale_x_log10()+geom_vline(xintercept=rna_thr, linetype="dashed")
ggsave("nFeature_RNA.hist.pdf",p,width=4,height=3)


p<-ggplot(plt_df,aes(x=nCount_ADT,y=nFeature_ADT,color=plate))+
geom_point(size=0.1)+theme_test()+theme(aspect.ratio=1/1)+
scale_color_viridis(option="plasma",discrete=TRUE)+scale_x_log10()+scale_y_log10()+geom_hline(yintercept=adt_thr, linetype="dashed")
ggsave("nCount.Feature_ADT.scatter.pdf",p,width=4,height=3)

p<-ggplot(plt_df,aes(x=nFeature_ADT,fill=plate))+
geom_histogram(position = "stack")+theme_bw()+theme(aspect.ratio=1/1)+
scale_fill_viridis(option="plasma",discrete=TRUE)+
scale_x_log10()+geom_vline(xintercept=adt_thr, linetype="dashed")
ggsave("nFeature_ADT.hist.pdf",p,width=4,height=3)


p<-ggplot(plt_df,aes(x=percent.mt,fill=plate))+
geom_histogram(position = "stack")+theme_test()+theme(aspect.ratio=1/1)+geom_vline(xintercept=mt_thr)+
scale_fill_viridis(option="plasma",discrete=TRUE)
ggsave("percent.mt.hist.pdf",p,width=4,height=3)



p<-ggplot(filter(plt_df,percent.mt<mt_thr),aes(x=nCount_RNA,y=nFeature_RNA,color=plate))+
geom_point(size=0.1)+theme_test()+theme(aspect.ratio=1/1)+
scale_color_viridis(option="plasma",discrete=TRUE)+scale_x_log10()+scale_y_log10()+geom_hline(yintercept=rna_thr, linetype="dashed")
ggsave("nCount.Feature_RNA.scatter.mt.filter.pdf",p,width=4,height=3)

p<-ggplot(filter(plt_df,percent.mt<mt_thr),aes(x=nFeature_RNA,fill=plate))+
geom_histogram(position = "stack")+theme_bw()+theme(aspect.ratio=1/1)+
scale_fill_viridis(option="plasma",discrete=TRUE)+
scale_x_log10()+geom_vline(xintercept=rna_thr, linetype="dashed")
ggsave("nFeature_RNA.hist.mt.filter.pdf",p,width=4,height=3)


p<-ggplot(filter(plt_df,percent.mt>=mt_thr),aes(x=nCount_RNA,y=nFeature_RNA,color=plate))+
geom_point(size=0.1)+theme_test()+theme(aspect.ratio=1/1)+
scale_color_viridis(option="plasma",discrete=TRUE)+scale_x_log10()+scale_y_log10()+geom_hline(yintercept=rna_thr, linetype="dashed")
ggsave("nCount.Feature_RNA.scatter.mt.filter.out.pdf",p,width=4,height=3)

p<-ggplot(filter(plt_df,percent.mt>=mt_thr),aes(x=nFeature_RNA,fill=plate))+
geom_histogram(position = "stack")+theme_bw()+theme(aspect.ratio=1/1)+
scale_fill_viridis(option="plasma",discrete=TRUE)+
scale_x_log10()+geom_vline(xintercept=rna_thr, linetype="dashed")
ggsave("nFeature_RNA.hist.mt.filter.out.pdf",p,width=4,height=3)
