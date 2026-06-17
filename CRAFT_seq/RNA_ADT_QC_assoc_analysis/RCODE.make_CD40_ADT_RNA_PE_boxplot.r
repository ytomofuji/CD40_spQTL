## ---------------------------
##
## Script name: RCODE.make_CD40_ADT_RNA_PE_boxplot.r
##
## Original file: RCODE.260129.make.CD40.ADT.RNA.PE.boxplot.narrow.log.r
##
## Purpose of script: Create boxplots of CD40 RNA/ADT/PE levels by allele dosage
##                    for all CRAFT-seq experimental conditions
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
library(Seurat)
library(ggpointdensity)
library(viridis)
library(lme4)
library(patchwork)
library(ggrastr)
library(ggrepel)
library(furrr)
library(future)
n_fun<-function(x){
  return(data.frame(y=max(x)*1.1,label=paste0("N = ",length(x))))
}


file="/path/to/E1.genotype_kmeans.after.norm.w.CD40.tsv.gz"
E1<-read_tsv(file)%>%filter(condition!="NTC")

file="/path/to/E2.genotype_kmeans.after.norm.w.CD40.tsv.gz"
E2<-read_tsv(file)%>%filter(condition!="NTC")

file="/path/to/E3.genotype_kmeans.after.norm.w.CD40.tsv.gz"
E3<-read_tsv(file)%>%filter(condition!="NTC")

file="/path/to/pQTL_genotype_kmeans.after.norm.w.CD40.tsv.gz"
pQTL<-read_tsv(file)%>%filter(condition!="NTC")

file="/path/to/KO_genotype_kmeans.after.norm.w.CD40.tsv.gz"
KO<-read_tsv(file)%>%filter(condition!="NTC")

file="/path/to/KOall_genotype_kmeans.after.norm.w.CD40.tsv.gz"
KOall<-read_tsv(file)%>%filter(condition!="NTC")

formula="logCD40 ~ num_allele_dosage + DNA_Barcode"
seurat_file="/path/to/raw.RNA.ADT.object.rds"
data<-readRDS(seurat_file)


CD40_RNA<-subset(data, nFeature_RNA > 500 & percent.mt <20)%>%FetchData(vars = "CD40")%>%
  rownames_to_column(var = "cell_barcode") %>%
  as_tibble()%>%rename(CD40_RNA=CD40)

DefaultAssay(data)<-"ADT"

CD40_ADT<-subset(data, nFeature_ADT > 50)%>%FetchData(vars = "CD40")%>%
  rownames_to_column(var = "cell_barcode") %>%
  as_tibble()%>%rename(CD40_ADT=CD40)

CD40_df<-dplyr::select(data@meta.data,cell,plate,name)%>%
mutate(cell_name=str_c(plate,"_",name))%>%as_tibble()%>%
left_join(CD40_ADT,by=c("cell"="cell_barcode"))%>%
left_join(CD40_RNA,by=c("cell"="cell_barcode"))

RNA_max<-1.1*CD40_df$CD40_RNA[!is.na(CD40_df$CD40_RNA)]%>%max()
ADT_max<-1.1*CD40_df$CD40_ADT[!is.na(CD40_df$CD40_ADT)]%>%max()


lm_res<-tibble()
for(dataset in c("KO","KOall","E1","E2","E3","pQTL")){
tmp_data<-eval(parse(text = dataset))
tmp_data$num_allele_dosage<-as.numeric(tmp_data$allele_dosage)
tmp_data$allele_dosage<-factor(as.character(tmp_data$allele_dosage), levels=c("0","1","2"))
tmp_data<-left_join(tmp_data,CD40_df,by="cell_name")%>%mutate(logCD40=log10(CD40))


p_rna_tmp<-ggplot(tmp_data,aes(x=allele_dosage,y=CD40_RNA))+
geom_boxplot(outlier.shape=NA,colour="black",aes(fill=allele_dosage))+
geom_jitter_rast(alpha=0.5,size=0.1,height=0,width=0.3,colour="#000000",shape=21)+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=3/2,legend.position="none",
	plot.title = element_text(size = 4,  hjust = 0.5)
)+
scale_fill_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
xlab(dataset)+ylab("log(CD40-RNA)")+ylim(c(0,RNA_max))


p_adt_tmp<-ggplot(tmp_data,aes(x=allele_dosage,y=CD40_ADT))+
geom_boxplot(outlier.shape=NA,colour="black",aes(fill=allele_dosage))+
geom_jitter_rast(alpha=0.5,size=0.1,height=0,width=0.3,colour="#000000",shape=21)+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=3/2,legend.position="none",
	plot.title = element_text(size = 4,  hjust = 0.5)
)+
scale_fill_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
xlab(dataset)+ylab("clr(CD40-ADT)")+ylim(c(0,ADT_max))


p_pe_tmp<-ggplot(tmp_data,aes(x=allele_dosage,y=logCD40))+
geom_boxplot(outlier.shape=NA,colour="black",aes(fill=allele_dosage))+
geom_jitter_rast(alpha=0.5,size=0.1,height=0,width=0.3,colour="#000000",shape=21)+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=3/2,legend.position="none",
	plot.title = element_text(size = 4,  hjust = 0.5)
)+
scale_fill_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
xlab(dataset)+ylab("log(CD40-PE)")

linear_reg<-lm(formula,tmp_data)
tmp_res<-summary(linear_reg)$coefficients["num_allele_dosage",]%>%as.data.frame()%>%t()%>%as_tibble()
colnames(tmp_res)<-c("BETA","SE","t","P")
tmp_res<-mutate(tmp_res,condition=dataset)
lm_res<-bind_rows(lm_res,tmp_res)

if(dataset=="KO"){
	p_rna<-p_rna_tmp
	p_adt<-p_adt_tmp
	p_pe<-p_pe_tmp
}else{
	p_rna<-p_rna+p_rna_tmp
	p_adt<-p_adt+p_adt_tmp
	p_pe<-p_pe+p_pe_tmp
}
}

p_rna<-p_rna+plot_layout(ncol=6)
p_adt<-p_adt+plot_layout(ncol=6)
p_pe<-p_pe+plot_layout(ncol=6)


ggsave("CRAFT.seq.CD40.enh.RNA.boxplot.narrow.pdf",p_rna,width=12,height=4)
ggsave("CRAFT.seq.CD40.enh.ADT.boxplot.narrow.pdf",p_adt,width=12,height=4)
ggsave("CRAFT.seq.CD40.enh.PE.boxplot.narrow.pdf",p_pe,width=12,height=4)

write_tsv(lm_res, "CRAFT.seq.CD40.enh.PE.narrow.lm.tsv")



#also prepare plot showing the number of genotype

rna_allele_df<-tibble()
adt_allele_df<-tibble()
pe_allele_df<-tibble()


for(dataset in c("KOall","E1","E2","E3","pQTL")){
tmp_data<-eval(parse(text = dataset))
tmp_data$num_allele_dosage<-as.numeric(tmp_data$allele_dosage)
tmp_data$allele_dosage<-factor(as.character(tmp_data$allele_dosage), levels=c("0","1","2"))
tmp_data<-left_join(tmp_data,CD40_df,by="cell_name")%>%mutate(logCD40=log10(CD40))

tmp <- filter(tmp_data,!is.na(CD40_RNA)) %>% group_by(allele_dosage) %>% summarize(Count=n()) %>% mutate(experiment=dataset)
rna_allele_df<-bind_rows(rna_allele_df,tmp)

tmp <- filter(tmp_data,!is.na(CD40_ADT)) %>% group_by(allele_dosage) %>% summarize(Count=n()) %>% mutate(experiment=dataset)
adt_allele_df<-bind_rows(adt_allele_df,tmp)

tmp <- filter(tmp_data,!is.na(CD40)) %>% group_by(allele_dosage) %>% summarize(Count=n()) %>% mutate(experiment=dataset)
pe_allele_df<-bind_rows(pe_allele_df,tmp)
}


rna_allele_df$allele_dosage<-factor(rna_allele_df$allele_dosage, levels=rev(c(0,1,2)))
rna_allele_df$experiment<-factor(rna_allele_df$experiment, levels=rev(c("KOall","E1","E2","E3","pQTL")))

p<-ggplot(rna_allele_df,aes(x=Count,y=experiment,fill=allele_dosage))+
geom_bar(stat="identity", position="stack")+
geom_text(aes(label = Count), 
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black",
    fontface = "bold") +
scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
scale_fill_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
theme_test()+theme(aspect.ratio=1/1)

ggsave("Barplot.for.genotype.count.CRAFT.seq.enh.RNA.pdf",p,width=4,height=4)





adt_allele_df$allele_dosage<-factor(adt_allele_df$allele_dosage, levels=rev(c(0,1,2)))
adt_allele_df$experiment<-factor(adt_allele_df$experiment, levels=rev(c("KOall","E1","E2","E3","pQTL")))

p<-ggplot(adt_allele_df,aes(x=Count,y=experiment,fill=allele_dosage))+
geom_bar(stat="identity", position="stack")+
geom_text(aes(label = Count), 
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black",
    fontface = "bold") +
scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
scale_fill_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
theme_test()+theme(aspect.ratio=1/1)

ggsave("Barplot.for.genotype.count.CRAFT.seq.enh.ADT.pdf",p,width=4,height=4)






pe_allele_df$allele_dosage<-factor(pe_allele_df$allele_dosage, levels=rev(c(0,1,2)))
pe_allele_df$experiment<-factor(pe_allele_df$experiment, levels=rev(c("KOall","E1","E2","E3","pQTL")))

p<-ggplot(pe_allele_df,aes(x=Count,y=experiment,fill=allele_dosage))+
geom_bar(stat="identity", position="stack")+
geom_text(aes(label = Count), 
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black",
    fontface = "bold") +
scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
scale_fill_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
theme_test()+theme(aspect.ratio=1/1)

ggsave("Barplot.for.genotype.count.CRAFT.seq.enh.PE.pdf",p,width=4,height=4)