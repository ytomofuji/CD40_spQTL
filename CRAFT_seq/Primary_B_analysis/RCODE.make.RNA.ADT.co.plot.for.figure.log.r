## ---------------------------
##
## Script name: RCODE.make.RNA.ADT.co.plot.for.figure.log.r
##
## Original file: RCODE.make.RNA.ADT.co.plot.for.figure.log.r
##
## Purpose of script: Compare RNA/ADT modalities for CRAFT-seq
##
##
## Date Created: 2026-02-02
##
## Email:
##
## ---------------------------
##
## Notes: To make concordance, we flipped the allele
##
##
## ---------------------------


library(tidyverse)
library(ggrastr)
library(ggrepel)

argv=commandArgs(T)
RNA_file=argv[1]
ADT_file=argv[2]
bridge_file=argv[3]
gene_list=argv[4]

if(FALSE){
RNA_file="/path/to/RNA.genotype.nb.result.tsv.gz"
ADT_file="/path/to/ADT.genotype.nb.result.tsv.gz"
bridge_file="/path/to/CRAFT.to.Gene.Name.tbl.tsv"
gene_list="FCER2,CD86,CD27,CD40,CD58"
}

rna_data<-read_tsv(RNA_file)
rna_data$Gene<-gsub("_","-",rna_data$Gene)
adt_data<-read_tsv(ADT_file)
bridge_df<-read_tsv(bridge_file)
gene_list<-str_split(gene_list,",")%>%unlist()

df<-inner_join(bridge_df,adt_data,by=c("CRAFT"="Gene"))%>%inner_join(rna_data,by="Gene")%>%
mutate(highlight=ifelse(Gene %in% gene_list,"Yes","No"))%>%
mutate(BETA.x=(-1*BETA.x),BETA.y=(-1*BETA.y))%>%filter(Gene!="PTPRC")


#calculate correlation after removing CD40
cor_df<-filter(df,Gene!="CD40")
cor_test<-cor.test(cor_df$BETA.x,cor_df$BETA.y)


p<-ggplot(df,aes(x=BETA.x,y=BETA.y, color=highlight))+
  geom_text_repel(data=subset(df,Gene%in%gene_list),aes(label=Gene),shape = 21)+
  scale_color_manual(values=c("Yes"="#FF3333","No"="black"))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_errorbar(aes(ymin = BETA.y - 1.96*SE.y,
                    ymax = BETA.y + 1.96*SE.y), width = 0) +
  geom_errorbarh(aes(xmin = BETA.x - 1.96*SE.x,
                     xmax = BETA.x + 1.96*SE.x), height = 0) +
  geom_point(size = 2, shape = 16)+
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "none"
  )+xlab("ADT effect size")+ylab("RNA effect size")+
  annotate(	"text", x = -Inf, y = Inf, label=str_c("P=",signif(cor_test$p.value,digits=2),
  "\nPearson's correlation=",signif(cor_test$estimate,digits=2)),
  hjust = -0.1, vjust = 1.1, size = 3)

ggsave("RNA.vs.ADT.primary.B.CRAFT.seq.pdf",p,width=3,height=3)




p<-ggplot(df,aes(x=BETA.x,y=BETA.y, color=highlight))+
  geom_text_repel(data=subset(df,Gene%in%gene_list),aes(label=Gene),shape = 21)+
  scale_color_manual(values=c("Yes"="#FF3333","No"="black"))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_point(size = 2, shape = 16)+
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "none"
  )+xlab("ADT effect size")+ylab("RNA effect size")+
  annotate(	"text", x = -Inf, y = Inf, label=str_c("P=",signif(cor_test$p.value,digits=2),
  "\nPearson's correlation=",signif(cor_test$estimate,digits=2)),
  hjust = -0.1, vjust = 1.1, size = 3)

ggsave("RNA.vs.ADT.primary.B.CRAFT.seq.no.SE.pdf",p,width=3,height=3)
