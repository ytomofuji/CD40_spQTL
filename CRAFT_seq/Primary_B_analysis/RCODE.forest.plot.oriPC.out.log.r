## ---------------------------
##
## Script name: RCODE.forest.plot.oriPC.out.log.r
##
## Original file: RCODE.260318.forest.plot.oriPC.out.log.r
##
## Purpose of script: For RNA/ADT, check association to genotype
##
##
## Date Created: 2026-03-18
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

library(tidyverse)
library(Seurat)
library(ggrepel)
n_fun<-function(x){
  return(data.frame(y=max(x)*1.1,label=paste0("N = ",length(x))))
}



# as a statistical model, try nb model
glm_test <- function(df, full_form_right, null_form_right, variable, gene) {

full_form <- as.formula(str_c(gene, " ~ ", full_form_right))
null_form <- as.formula(str_c(gene, " ~ ", null_form_right))

full_model <- MASS::glm.nb(formula = full_form, data = df)
null_model <- MASS::glm.nb(formula = null_form, data = df)

anova_res <- anova(full_model, null_model, test = "Chisq")
lr_stat <- anova_res$"LR stat."[2]
df_diff <- anova_res$"   df"[2]
log_p_natural <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE, log.p = TRUE)
minus_log10p <- - (log_p_natural / log(10))
LRT_P_val <- 10**(-minus_log10p)


  out <- summary(full_model)$coefficients
  colnames(out) <- c("BETA", "SE", "Z", "P")

  out <- bind_cols(
    tibble(Gene = gene),
    as_tibble(t(out[variable, ]))
  ) %>%
  mutate(LRT_P = LRT_P_val)

  return(out)
}



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




rna_cd40_full <- "C_allele_dosage + plate_exp + logRNA + percent.mt + PC_1 + PC_2"
rna_cd40_null <- "plate_exp + logRNA + percent.mt + PC_1 + PC_2"

adt_cd40_full <- "C_allele_dosage + plate_exp + logADT + percent.mt + PC_1 + PC_2"
adt_cd40_null <- "plate_exp + logADT + percent.mt + PC_1 + PC_2"

test_df<-df
rm(df)

rna_res<-tibble()
rna_res<-bind_rows(rna_res,glm_test(test_df,rna_cd40_full,rna_cd40_null,"C_allele_dosage","rna_CD40")%>%mutate(experiment="All"))
adt_res<-tibble()
adt_res<-bind_rows(adt_res,glm_test(test_df,adt_cd40_full,adt_cd40_null,"C_allele_dosage","adt_CD40")%>%mutate(experiment="All"))

experiments<-unique(test_df$experiment)
for(i in experiments){
tmp<-filter(test_df,experiment==i)
rna_res<-bind_rows(rna_res,glm_test(tmp,rna_cd40_full,rna_cd40_null,"C_allele_dosage","rna_CD40")%>%mutate(experiment=i))
adt_res<-bind_rows(adt_res,glm_test(tmp,adt_cd40_full,adt_cd40_null,"C_allele_dosage","adt_CD40")%>%mutate(experiment=i))
}


write_tsv(rna_res,"all.clst.cis.analysis.oriPC.cov.RNA.tsv")
write_tsv(adt_res,"all.clst.cis.analysis.oriPC.cov.ADT.tsv")



#make forest plot
rna_res$experiment<-factor(rna_res$experiment, levels=rev(c("PE1","PE2","ABE","All")))
adt_res$experiment<-factor(adt_res$experiment, levels=rev(c("PE1","PE2","ABE","All")))

maxlim<-max(c(-rna_res$BETA+rna_res$SE*1.96,-adt_res$BETA+adt_res$SE*1.96))*1.1
minlim<-min(c(-rna_res$BETA-rna_res$SE*1.96,-adt_res$BETA-adt_res$SE*1.96))*1.1


#flip effect so that C allele will be alt
p<-ggplot(rna_res,aes(y=experiment,x=(-BETA),color=experiment))+
geom_point(aes(shape=experiment,size=experiment))+geom_errorbar(aes(xmin=(-BETA)-1.96*SE,xmax=(-BETA)+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("PE1"="#FF3333","PE2"="#990000","ABE"="#3399FF","All"="#FF1493"))+
scale_shape_manual(values=c("PE1"=20,"PE2"=20,"ABE"=20,"All"=18))+
scale_size_manual(values=c("PE1"=2,"PE2"=2,"ABE"=2,"All"=3))+
xlim(c(minlim,maxlim))+
theme_test()+theme(aspect.ratio=2/1,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")

ggsave("all.clst.cis.analysis.oriPC.cov.RNA.forest.pdf",p,width=3,height=3)


p<-ggplot(adt_res,aes(y=experiment,x=(-BETA),color=experiment))+
geom_point(aes(shape=experiment,size=experiment))+geom_errorbar(aes(xmin=(-BETA)-1.96*SE,xmax=(-BETA)+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("PE1"="#FF3333","PE2"="#990000","ABE"="#3399FF","All"="#FF1493"))+
scale_shape_manual(values=c("PE1"=20,"PE2"=20,"ABE"=20,"All"=18))+
scale_size_manual(values=c("PE1"=2,"PE2"=2,"ABE"=2,"All"=3))+
xlim(c(minlim,maxlim))+
theme_test()+theme(aspect.ratio=2/1,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")


ggsave("all.clst.cis.analysis.oriPC.cov.ADT.forest.pdf",p,width=3,height=3)
