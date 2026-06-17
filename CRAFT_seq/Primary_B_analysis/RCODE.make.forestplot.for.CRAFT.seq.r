## ---------------------------
##
## Script name: RCODE.make.forestplot.for.CRAFT.seq.r
##
## Original file: RCODE.260202.make.forestplot.for.CRAFT.seq.r
##
## Purpose of script: Make forest plot from the table
##
##
## Date Created: 2026-02-02
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
library(viridis)
library(patchwork)
library(ggrastr)
library(ggrepel)


argv=commandArgs(T)

input_file=argv[1]
target_phenotype=argv[2]
modality=argv[3]


if(FALSE){
input_file="/path/to/meta.analysis.result.tsv.gz"
target_phenotype="EBI3"
modality="RNA"
}


data<-read_tsv(input_file)

df <- filter(data,Gene==target_phenotype)%>%
  pivot_longer(
    cols = matches("^(PE1|PE2|ABE|ALL)_"),
    names_to = c("experiment", ".value"),
    names_pattern = "^(PE1|PE2|ABE|ALL)_(.*)"
  )%>%
  mutate(BETA=-BETA)


df$experiment<-factor(df$experiment, levels=c("PE1","PE2","ABE","ALL"))



maxlim<-max(c(df$BETA+df$SE*1.96,0))+0.01
minlim<-min(c(df$BETA-df$SE*1.96,0))-0.01


p<-ggplot(df,aes(x=experiment,y=BETA,color=experiment))+
geom_point(aes(shape=experiment,size=experiment))+geom_errorbar(aes(ymin=BETA-1.96*SE,ymax=BETA+1.96*SE),width=0)+
geom_text_repel(data=subset(df,experiment=="ALL"),aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("PE1"="#FF3333","PE2"="#990000","ABE"="#3399FF","ALL"="#FF1493"))+
scale_shape_manual(values=c("PE1"=20,"PE2"=20,"ABE"=20,"ALL"=18))+
scale_size_manual(values=c("PE1"=3,"PE2"=3,"ABE"=3,"ALL"=4))+
ylim(c(minlim,maxlim))+
theme_test()+theme(aspect.ratio=2/1,legend.position="none")+geom_hline(yintercept=0,linetype="dashed")

ggsave(str_c(target_phenotype,".",modality,".forest.pdf"),p,width=2,height=2)
