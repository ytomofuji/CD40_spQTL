## ---------------------------
##
## Script name: RCODE.make.boxplot.for.primary.B.cell.CRAFT.seq.per.cluster.log.r
##
## Original file: RCODE.260318.make.boxplot.for.primary.B.cell.CRAFT.seq.per.cluster.log.r
##
## Purpose of script: Generate boxplot for primary B cell CRAFT-seq data
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
library(patchwork)
library(ggrastr)
library(ggrepel)


#plot boxplot with all dots
n_fun<-function(x){
  return(data.frame(y=max(x)*1.1,label=paste0("N = ",length(x))))
}



#First, plot the clustering UMAP and cluster ratio per sample

data<-readRDS(seurat_file)
Idents(data)<-"seurat_clusters"

df1 <- FetchData(data, vars = c("rna_CD40", "adt_CD40"),slot="count")%>%as.data.frame()%>%
rownames_to_column("barcode")%>%as_tibble()
df2<-data@meta.data%>%as.data.frame()%>%
rownames_to_column("barcode")%>%as_tibble()
df3 <- FetchData(data, vars = c("rna_CD40", "adt_CD40"),slot="data")%>%as.data.frame()%>%
rownames_to_column("barcode")%>%as_tibble()
colnames(df3)<-c("barcode", "norm_rna_CD40", "norm_adt_CD40")

df<-inner_join(df1,df2,by="barcode")%>%inner_join(df3,by="barcode")%>%
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

#make boxplot per dataset
plot_box<-function(DF){
p_rna<-ggplot(DF,aes(x=Genotype,y=norm_rna_CD40))+
geom_boxplot(outlier.shape=NA,colour="black",aes(fill=Genotype))+
geom_jitter_rast(alpha=0.1,size=0.01,height=0,width=0.2,shape=21,color="black")+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=2/1,legend.position="none",
	theme(plot.title = element_text(size = 4,  hjust = 0.5))
)+
scale_fill_manual(values = c("C/C" = "#1b9e77","C/T" = "#d95f02","T/T" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
ylab("CD40 RNA expression")+facet_wrap(~experiment,ncol=3,scales="free_y")

p_adt<-ggplot(DF,aes(x=Genotype,y=norm_adt_CD40))+
geom_boxplot(outlier.shape=NA,colour="black",aes(fill=Genotype))+
geom_jitter_rast(alpha=0.1,size=0.01,height=0,width=0.2,shape=21,color="black")+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=2/1,legend.position="none",
	theme(plot.title = element_text(size = 4,  hjust = 0.5))
)+
scale_fill_manual(values = c("C/C" = "#1b9e77","C/T" = "#d95f02","T/T" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
ylab("CD40 ADT expression")+facet_wrap(~experiment,ncol=3,scales="free_y")

p_pe<-ggplot(DF,aes(x=Genotype,y=log_FACS_CD40))+
geom_boxplot(outlier.shape=NA,colour="black",aes(fill=Genotype))+
geom_jitter_rast(alpha=0.1,size=0.01,height=0,width=0.2,shape=21,color="black")+
theme_test()+
theme(axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),aspect.ratio=2/1,legend.position="none",
	theme(plot.title = element_text(size = 4,  hjust = 0.5))
)+
scale_fill_manual(values = c("C/C" = "#1b9e77","C/T" = "#d95f02","T/T" = "#7570b3"))+
stat_summary(fun.data=n_fun,geom="text",colour="black",size=3,aes(group=1))+
ylab("CD40-PE")+facet_wrap(~experiment,ncol=3,scales="free_y")

p<-p_rna+p_adt+p_pe+plot_layout(ncol=1)
return(p)
}


plot<-plot_box(df)
ggsave("CD40.primary.B.cell.CRAFT.seq.boxplot.per.dataset.pdf",plot,width=10,height=10)

for(clst_num in seq(0,3)){
	tmp_df<-filter(df,seurat_clusters==clst_num)
	plot<-plot_box(tmp_df)
	ggsave(str_c("clsuter",clst_num,".CD40.primary.B.cell.CRAFT.seq.boxplot.per.dataset.pdf"),plot,width=10,height=10)
}


extract_stats<-function(lm_res,null_res){
	x<-summary(lm_res)$coefficients["C_allele_dosage",]%>%t()%>%as_tibble()
	colnames(x)<-c("BETA","SE","STAT","Wald_P")
	if(length(null_res)>0){
	anova_res <- anova(null_res, lm_res, test = "Chisq")
	lr_stat <- anova_res$"LR stat."[2]
	df_diff <- anova_res$"   df"[2]
	log_p_natural <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE, log.p = TRUE)
	minus_log10p <- - (log_p_natural / log(10))
	P <- 10**(-minus_log10p)
	x<-mutate(x,minus_log10LRT_P=minus_log10p,LRT_P=P)
	}
	return(x)
}




#Run regression analysis for all experiment / per experiment
result<-tibble()


for(targ_experiment in c("All","PE1","PE2","ABE")){

if(targ_experiment %in% c("PE1","PE2","ABE")){
	tmp_df2<-filter(df,experiment==targ_experiment)
}else{
	tmp_df2<-df
}

rna_cd40_reg <- MASS::glm.nb(formula = "rna_CD40 ~ C_allele_dosage + plate_exp + logRNA + percent.mt", data = tmp_df2)
rna_cd40_null <- MASS::glm.nb(formula = "rna_CD40 ~ plate_exp + logRNA + percent.mt", data = tmp_df2)

adt_cd40_reg <- MASS::glm.nb(formula = "adt_CD40 ~ C_allele_dosage + plate_exp + logADT + percent.mt", data = tmp_df2)
adt_cd40_null <- MASS::glm.nb(formula = "adt_CD40 ~ plate_exp + logADT + percent.mt", data = tmp_df2)

facs_cd40_reg <- lm(formula = "log_FACS_CD40 ~ C_allele_dosage + plate_exp", data = tmp_df2)

tmp <-bind_rows(extract_stats(rna_cd40_reg,rna_cd40_null)%>%mutate(Modality="RNA"),
extract_stats(adt_cd40_reg,adt_cd40_null)%>%mutate(Modality="ADT"),
extract_stats(facs_cd40_reg,c())%>%mutate(Modality="PE"))%>%mutate(cluster="all",experiment=targ_experiment)

result<-bind_rows(result,tmp)

for(clst_num in seq(0,3)){
	tmp_df<-filter(tmp_df2,seurat_clusters==clst_num)

if(nrow(tmp_df)>=10){
	rna_cd40_reg <- MASS::glm.nb(formula = "rna_CD40 ~ C_allele_dosage + plate_exp + logRNA + percent.mt", data = tmp_df)
	rna_cd40_null <- MASS::glm.nb(formula = "rna_CD40 ~ plate_exp + logRNA + percent.mt", data = tmp_df)

	adt_cd40_reg <- MASS::glm.nb(formula = "adt_CD40 ~ C_allele_dosage + plate_exp + logADT + percent.mt", data = tmp_df)
	adt_cd40_null <- MASS::glm.nb(formula = "adt_CD40 ~ plate_exp + logADT + percent.mt", data = tmp_df)

	facs_cd40_reg <- lm(formula = "log_FACS_CD40 ~ C_allele_dosage + plate_exp", data = tmp_df)

	tmp <-bind_rows(extract_stats(rna_cd40_reg,rna_cd40_null)%>%mutate(Modality="RNA"),
	extract_stats(adt_cd40_reg,adt_cd40_null)%>%mutate(Modality="ADT"),
	extract_stats(facs_cd40_reg,c())%>%mutate(Modality="PE"))%>%mutate(cluster=str_c("cluster",clst_num),experiment=targ_experiment)

	result<-bind_rows(result,tmp)}
}
}


write_tsv(result,"CD40.primary.B.cell.CRAFT.seq.regression.tsv")
result$Modality<-factor(result$Modality, levels=c("RNA","ADT","PE"))
result$cluster<-factor(result$cluster, levels=unique(result$cluster))
result$experiment<-factor(result$experiment, levels=unique(result$experiment))

#make forest plot for all-combined analysis
plt_df<-filter(result,cluster=="all")
plt_df$experiment<-factor(plt_df$experiment, levels=rev(c("PE1","PE2","ABE","All")))
#rescue for PE analysis
plt_df<-mutate(plt_df,LRT_P=ifelse(is.na(LRT_P),Wald_P,LRT_P))

p<-ggplot(plt_df,aes(y=experiment,x=BETA,color=experiment))+
geom_point(aes(shape=experiment,size=experiment))+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
scale_color_manual(values=c("PE1"="#FF3333","PE2"="#990000","ABE"="#3399FF","All"="#FF1493"))+
scale_shape_manual(values=c("PE1"=20,"PE2"=20,"ABE"=20,"All"=18))+
scale_size_manual(values=c("PE1"=2,"PE2"=2,"ABE"=2,"All"=3))+
theme_test()+theme(aspect.ratio=2/1,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")+
facet_wrap(~Modality,ncol=3,scales="free_x")

ggsave("All.clusters.CD40.primary.B.cell.CRAFT.seq.reg.forest.pdf",p,width=8,height=3)


#per cluster analysis
for(clst_num in seq(0,3)){

	#make forest plot for all-combined analysis
	plt_df<-filter(result,cluster==str_c("cluster",clst_num))
	plt_df$experiment<-factor(plt_df$experiment, levels=rev(c("PE1","PE2","ABE","All")))
	#rescue for PE analysis
	plt_df<-mutate(plt_df,LRT_P=ifelse(is.na(LRT_P),Wald_P,LRT_P))

	p<-ggplot(plt_df,aes(y=experiment,x=BETA,color=experiment))+
	geom_point(aes(shape=experiment,size=experiment))+geom_errorbar(aes(xmin=BETA-1.96*SE,xmax=BETA+1.96*SE),width=0)+
	geom_text_repel(aes(label=str_c("P=",signif(LRT_P,digits=2))))+
	scale_color_manual(values=c("PE1"="#FF3333","PE2"="#990000","ABE"="#3399FF","All"="#FF1493"))+
	scale_shape_manual(values=c("PE1"=20,"PE2"=20,"ABE"=20,"All"=18))+
	scale_size_manual(values=c("PE1"=2,"PE2"=2,"ABE"=2,"All"=3))+
	theme_test()+theme(aspect.ratio=1/1,legend.position="none")+geom_vline(xintercept=0,linetype="dashed")+
	facet_wrap(~Modality,ncol=1,scales="free_x")

	ggsave(str_c("cluster",clst_num,".CD40.primary.B.cell.CRAFT.seq.reg.forest.pdf"),p,width=2,height=8)
}



#make all forest plot (wide)

plt_df<-filter(result,cluster!="all")
plt_df$experiment<-factor(plt_df$experiment, levels=c("PE1","PE2","ABE","All"))
plt_df$cluster<-factor(plt_df$cluster, levels=c("cluster0","cluster2","cluster1","cluster3"))
plt_df<-mutate(plt_df,LRT_P=ifelse(is.na(LRT_P),Wald_P,LRT_P))


pd <- position_dodge(width = 0.6)

p <- ggplot(plt_df, aes(x = cluster, y = BETA, color = experiment)) +
geom_errorbar(aes(ymin = BETA - 1.96 * SE, ymax = BETA + 1.96 * SE),
			width = 0, position = pd) +
geom_point(aes(shape = experiment, size = experiment), position = pd) +
geom_text_repel(data=filter(plt_df,LRT_P<0.05),aes(label = str_c("P=", signif(LRT_P, digits = 2))),
				position = pd, size = 3, show.legend = FALSE) +
scale_color_manual(values = c("PE1" = "#FF3333", "PE2" = "#990000",
							"ABE" = "#3399FF", "All" = "#FF1493")) +
scale_shape_manual(values = c("PE1" = 20, "PE2" = 20, "ABE" = 20, "All" = 18)) +
scale_size_manual(values = c("PE1" = 2, "PE2" = 2, "ABE" = 2, "All" = 3)) +
theme_test() +
theme(aspect.ratio = 1/2, legend.position = "none") +
geom_hline(yintercept = 0, linetype = "dashed") +
facet_wrap(~Modality, ncol = 1, scales = "free_y")


ggsave("CD40.primary.B.cell.CRAFT.seq.reg.forest.per.clusters.pdf",p,width=4,height=5)

p <- ggplot(plt_df, aes(x = cluster, y = BETA, color = experiment)) +
geom_errorbar(aes(ymin = BETA - 1.96 * SE, ymax = BETA + 1.96 * SE),
			width = 0, position = pd) +
geom_point(aes(shape = experiment, size = experiment), position = pd) +
scale_color_manual(values = c("PE1" = "#FF3333", "PE2" = "#990000",
							"ABE" = "#3399FF", "All" = "#FF1493")) +
scale_shape_manual(values = c("PE1" = 20, "PE2" = 20, "ABE" = 20, "All" = 18)) +
scale_size_manual(values = c("PE1" = 2, "PE2" = 2, "ABE" = 2, "All" = 3)) +
theme_test() +
theme(aspect.ratio = 1/2, legend.position = "none") +
geom_hline(yintercept = 0, linetype = "dashed") +
facet_wrap(~Modality, ncol = 1, scales = "free_y")


ggsave("CD40.primary.B.cell.CRAFT.seq.reg.forest.per.clusters.noP.pdf",p,width=4,height=5)
