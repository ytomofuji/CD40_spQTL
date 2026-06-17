## ---------------------------
##
## Script name: RCODE.make_forest_plot_CRAFT_enh.r
##
## Original file: RCODE.260129.make.forest.plot.enh.CRAFT.narrow.one.log.r
##
## Purpose of script: Generate forest plots for CD40 CRAFT-seq association results
##                    across experimental conditions
##
##
## Date Created: 2026-01-29
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
library(ggrepel)

adt_file="/path/to/CD40.assoc.ADT.summary.tsv"
rna_file="/path/to/CD40.assoc.RNA.summary.tsv"
pe_file="/path/to/CRAFT.seq.CD40.enh.PE.lm.tsv"

adt_df<-read_tsv(adt_file)%>%filter(Condition!="KO")%>%mutate(condition=ifelse(Condition=="KOall","KO",Condition))
rna_df<-read_tsv(rna_file)%>%filter(Condition!="KO")%>%mutate(condition=ifelse(Condition=="KOall","KO",Condition))
pe_df<-read_tsv(pe_file)%>%filter(condition!="KO")%>%mutate(condition=ifelse(condition=="KOall","KO",condition))

adt_df$condition<-factor(adt_df$condition, levels=rev(c("KO","E1","E2","E3","pQTL")))
rna_df$condition<-factor(rna_df$condition, levels=rev(c("KO","E1","E2","E3","pQTL")))
pe_df$condition<-factor(pe_df$condition, levels=rev(c("KO","E1","E2","E3","pQTL")))


p<-ggplot(adt_df,aes(y=condition,x=BETA,color=condition))+
geom_point()+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("E1"="#1F78B4","E2"="#33A02C","E3"="#E31A1C","pQTL"="#FF7F00","KO"="#6A3D9A"))+
theme_test()+theme(aspect.ratio=3/2,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("Enh.CRAFT.seq.ADT.forest.short.narrow.pdf",p,width=2,height=3)


p<-ggplot(rna_df,aes(y=condition,x=BETA,color=condition))+
geom_point()+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("E1"="#1F78B4","E2"="#33A02C","E3"="#E31A1C","pQTL"="#FF7F00","KO"="#6A3D9A"))+
theme_test()+theme(aspect.ratio=3/2,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("Enh.CRAFT.seq.RNA.forest.short.narrow.pdf",p,width=2,height=3)


p<-ggplot(pe_df,aes(y=condition,x=BETA,color=condition))+
geom_point()+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(P,digits=2))))+
scale_color_manual(values=c("E1"="#1F78B4","E2"="#33A02C","E3"="#E31A1C","pQTL"="#FF7F00","KO"="#6A3D9A"))+
theme_test()+theme(aspect.ratio=3/2,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("Enh.CRAFT.seq.PE.forest.short.narrow.pdf",p,width=2,height=3)



#also add replication data (1st CRAFT-seq experiment)
adt_file_old="/path/to/replication/CD40.assoc.ADT.summary.tsv"
rna_file_old="/path/to/replication/CD40.assoc.RNA.summary.tsv"
pe_file_old="/path/to/replication/CRAFT.seq.CD40.1st.PE.lm.tsv"

adt_old_df<-read_tsv(adt_file_old)%>%mutate(condition=str_c(Condition,"_rep"))
rna_old_df<-read_tsv(rna_file_old)%>%mutate(condition=str_c(Condition,"_rep"))
pe_old_df<-read_tsv(pe_file_old)%>%mutate(condition=str_c(condition,"_rep"))

adt_old_df$condition<-factor(adt_old_df$condition, levels=rev(c("KO_rep","pQTL_rep")))
rna_old_df$condition<-factor(rna_old_df$condition, levels=rev(c("KO_rep","pQTL_rep")))
pe_old_df$condition<-factor(pe_old_df$condition, levels=rev(c("KO_rep","pQTL_rep")))


p<-ggplot(adt_old_df,aes(y=condition,x=BETA,color=condition))+
geom_point()+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("E1"="#1F78B4","E2"="#33A02C","E3"="#E31A1C","pQTL"="#FF7F00","KO"="#6A3D9A",
                                 "pQTL_rep" = "#FFB366","KO_rep"="#B9A3D3"))+
theme_test()+theme(aspect.ratio=3/2,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("1st.CRAFT.seq.ADT.forest.short.narrow.pdf",p,width=2,height=2)


p<-ggplot(rna_old_df,aes(y=condition,x=BETA,color=condition))+
geom_point()+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("E1"="#1F78B4","E2"="#33A02C","E3"="#E31A1C","pQTL"="#FF7F00","KO"="#6A3D9A",
                                 "pQTL_rep" = "#FFB366","KO_rep"="#B9A3D3"))+
theme_test()+theme(aspect.ratio=3/2,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("1st.CRAFT.seq.RNA.forest.short.narrow.pdf",p,width=2,height=2)


p<-ggplot(pe_old_df,aes(y=condition,x=BETA,color=condition))+
geom_point()+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(P,digits=2))))+
scale_color_manual(values=c("E1"="#1F78B4","E2"="#33A02C","E3"="#E31A1C","pQTL"="#FF7F00","KO"="#6A3D9A",
                                 "pQTL_rep" = "#FFB366","KO_rep"="#B9A3D3"))+
theme_test()+theme(aspect.ratio=3/2,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("1st.CRAFT.seq.PE.forest.short.narrow.pdf",p,width=2,height=2)



#merge and run

adt_all_df<-bind_rows(adt_df,adt_old_df)
adt_all_df$condition<-factor(adt_all_df$condition, levels=rev(c("KO","E1","E2","E3","pQTL","KO_rep","pQTL_rep")))
rna_all_df<-bind_rows(rna_df,rna_old_df)
rna_all_df$condition<-factor(rna_all_df$condition, levels=rev(c("KO","E1","E2","E3","pQTL","KO_rep","pQTL_rep")))
pe_all_df<-bind_rows(pe_df,pe_old_df)
pe_all_df$condition<-factor(pe_all_df$condition, levels=rev(c("KO","E1","E2","E3","pQTL","KO_rep","pQTL_rep")))

p<-ggplot(adt_all_df,aes(y=condition,x=BETA,color=condition))+
geom_point()+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("E1"="#1F78B4","E2"="#33A02C","E3"="#E31A1C","pQTL"="#FF7F00","KO"="#6A3D9A",
                                 "pQTL_rep" = "#FFB366","KO_rep"="#B9A3D3"))+
theme_test()+theme(aspect.ratio=3/2,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("Enh.1st.CRAFT.seq.ADT.forest.short.narrow.pdf",p,width=2,height=3)


p<-ggplot(rna_all_df,aes(y=condition,x=BETA,color=condition))+
geom_point()+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("E1"="#1F78B4","E2"="#33A02C","E3"="#E31A1C","pQTL"="#FF7F00","KO"="#6A3D9A",
                                 "pQTL_rep" = "#FFB366","KO_rep"="#B9A3D3"))+
theme_test()+theme(aspect.ratio=3/2,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("Enh.1st.CRAFT.seq.RNA.forest.short.narrow.pdf",p,width=2,height=3)


p<-ggplot(pe_all_df,aes(y=condition,x=BETA,color=condition))+
geom_point()+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(P,digits=2))))+
scale_color_manual(values=c("E1"="#1F78B4","E2"="#33A02C","E3"="#E31A1C","pQTL"="#FF7F00","KO"="#6A3D9A",
                                 "pQTL_rep" = "#FFB366","KO_rep"="#B9A3D3"))+
theme_test()+theme(aspect.ratio=3/2,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("Enh.1st.CRAFT.seq.PE.forest.short.narrow.pdf",p,width=2,height=3)
