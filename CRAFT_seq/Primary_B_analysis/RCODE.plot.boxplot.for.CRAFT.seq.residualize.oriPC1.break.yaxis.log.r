## ---------------------------
##
## Script name: RCODE.plot.boxplot.for.CRAFT.seq.residualize.oriPC1.break.yaxis.log.r
##
## Original file: RCODE.260128.plot.boxplot.for.CRAFT.seq.residualize.oriPC1.break.yaxis.log.r
##
## Purpose of script: Prepare boxplot for the primary B cell CRAFT-seq
##
##
## Date Created: 2026-01-28
##
## Email:
##
## ---------------------------
##
## Notes:
##
##
## ---------------------------

argv=commandArgs(T)
seurat_file=argv[1]

if(FALSE){
seurat_file="/path/to/QCed.merge.primaryB.CRAFT.rds"
}

library(ggplot2)
library(dplyr)
library(Seurat)
library(ggrastr)
library(ggforce)
library(tibble)
library(stringr)

library(viridis)
library(patchwork)

library(ggrepel)
library(scales)
library(ggbreak)

library(tidyverse)

#First, plot the clustering UMAP and cluster ratio per sample

data<-readRDS(seurat_file)
Idents(data)<-"seurat_clusters"

df1 <- FetchData(data, vars = c("rna_CD40", "adt_CD40", "PC_1", "PC_2"),slot="count")%>%as.data.frame()%>%
rownames_to_column("barcode")%>%as_tibble()
df2<-data@meta.data%>%as.data.frame()%>%
rownames_to_column("barcode")%>%as_tibble()

df<-inner_join(df1,df2,by="barcode")%>%
mutate(cluster=case_when(
	seurat_clusters==0 ~ "Activated Naive B cells",
	seurat_clusters==1 ~ "Plasma B cells",
	seurat_clusters==2 ~ "Intermediate B cells",
	seurat_clusters==3 ~ "NR4A+ Activated B cells"))
df$cluster<-factor(df$cluster, levels=c("Activated Naive B cells","Intermediate B cells","Plasma B cells","NR4A+ Activated B cells"))
df$experiment<-factor(df$experiment, levels=c("PE1","PE2","ABE"))
df<-mutate(df,plate_exp=str_c(plate,experiment))%>%
mutate(logRNA=log10(nCount_RNA),logADT=log10(nCount_ADT))

df<-mutate(df,Genotype=case_when(C_allele_dosage==0~"T/T",
	C_allele_dosage==1~"C/T",
	C_allele_dosage==2~"C/C"))
df$Genotype<-factor(df$Genotype ,levels=rev(c("T/T","C/T","C/C")))


#make residualized CD40 expression
rna_cd40 <- MASS::glm.nb(formula = "rna_CD40 ~ plate_exp + logRNA + percent.mt + PC_1 + PC_2", data = df)
df$CD40_RNA_resid<-resid(rna_cd40, type="deviance")

adt_cd40 <- MASS::glm.nb(formula = "adt_CD40 ~ plate_exp + logADT + percent.mt + PC_1 + PC_2", data = df)
df$CD40_ADT_resid<-resid(adt_cd40, type="deviance")

facs_cd40 <- lm(formula = "log_FACS_CD40 ~ plate_exp + PC_1 + PC_2", data = df)
df$CD40_PE_resid<-resid(facs_cd40, type="deviance")


#plot boxplot with all dots
n_fun<-function(x){
  return(data.frame(y=max(x)*1.1,label=paste0("N = ",length(x))))
}




plot_box_rna <- function(DF){
#base plot
  base_p <- ggplot(DF, aes(x = Genotype, y = CD40_RNA_resid)) +
    geom_boxplot(outlier.shape = NA, colour = "black", aes(fill = Genotype)) +
    geom_jitter_rast(alpha = 0.1, size = 0.01, height = 0, width = 0.2, shape = 21) +
    scale_fill_manual(values = c("C/C" = "#1b9e77", "C/T" = "#d95f02", "T/T" = "#7570b3")) +
    theme_test() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  p_top <- base_p +
    coord_cartesian(ylim = c(2, max(df$CD40_RNA_resid))) +
	scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("")

  p_middle <- base_p +
    coord_cartesian(ylim = c(min(df$CD40_RNA_resid), 2)) +
	scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    ylab("CD40-RNA residual")


  combined <- p_top + p_middle + plot_layout(ncol=1,heights = c(0.1,1.2))

  return(combined)
}



plot_box_adt <- function(DF){
#base plot
  base_p <- ggplot(DF, aes(x = Genotype, y = CD40_ADT_resid)) +
    geom_boxplot(outlier.shape = NA, colour = "black", aes(fill = Genotype)) +
    geom_jitter_rast(alpha = 0.1, size = 0.01, height = 0, width = 0.2, shape = 21) +
    scale_fill_manual(values = c("C/C" = "#1b9e77", "C/T" = "#d95f02", "T/T" = "#7570b3")) +
    theme_test() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none",
      axis.line.y = element_line(color = "black", linewidth = 0.5)
    )

  p_top <- base_p +
    coord_cartesian(ylim = c(2, max(df$CD40_ADT_resid))) +
	scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ylab("")

  p_middle <- base_p +
    coord_cartesian(ylim = c(-2, 2)) +
	scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    ylab("CD40-ADT residual")

  p_bottom <- base_p +
    coord_cartesian(ylim = c(min(df$CD40_ADT_resid), -2)) +
	scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    ylab("")

  combined <- p_top + p_middle + p_bottom + plot_layout(ncol=1,heights = c(0.1,1.2,0.1))

  return(combined)
}



p <- plot_box_rna(df)
ggsave("Boxplot.for.all.clst.exp.oriPC.out.ybreak.RNA.pdf",p,width=2,height=4)

p <- plot_box_adt(df)
ggsave("Boxplot.for.all.clst.exp.oriPC.out.ybreak.ADT.pdf",p,width=2,height=4)



experiments<-unique(df$experiment)
for(i in experiments){

tmp<-filter(df,experiment==i)

p <- plot_box_rna(tmp)
ggsave(str_c("Boxplot.for.experiment.",i,".exp.oriPC.out.ybreak.RNA.pdf"),p,width=2,height=4)

p <- plot_box_adt(tmp)
ggsave(str_c("Boxplot.for.experiment.",i,".exp.oriPC.out.ybreak.ADT.pdf"),p,width=2,height=4)

}




#also prepare num-allele/cell bar plot
#manually change the color for the ABE experiment
allele_df<-df%>%group_by(experiment,allele_dosage)%>%summarize(Count=n())
allele_df$allele_dosage<-factor(allele_df$allele_dosage, levels=rev(c(0,1,2)))
allele_df$experiment<-factor(allele_df$experiment, levels=rev(c("PE1","PE2","ABE")))

p<-ggplot(allele_df,aes(x=Count,y=experiment,fill=allele_dosage))+
geom_bar(stat="identity", position="stack")+
geom_text(aes(label = Count),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black",
    fontface = "bold") +
scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
scale_fill_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
theme_test()+theme(aspect.ratio=1/2)

ggsave("Barplot.for.genotype.count.pdf",p,width=4,height=4)


allele_cnt_df<-df%>%group_by(experiment)%>%summarize(Edited=sum(allele_dosage),Non_edited=sum(2-allele_dosage))%>%
pivot_longer(names_to="allele",values_to="Count",-experiment)
allele_cnt_df$allele<-factor(allele_cnt_df$allele, levels=c("Edited","Non_edited"))
allele_cnt_df$experiment<-factor(allele_cnt_df$experiment, levels=rev(c("PE1","PE2","ABE")))

p<-ggplot(allele_cnt_df,aes(x=Count,y=experiment,fill=allele))+
geom_bar(stat="identity", position="stack")+
geom_text(aes(label = Count),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black",
    fontface = "bold") +
scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
theme_test()+theme(aspect.ratio=1/2)

ggsave("Barplot.for.allele.count.pdf",p,width=4,height=4)
