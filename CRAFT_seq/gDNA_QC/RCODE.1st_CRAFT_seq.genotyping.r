## ---------------------------
##
## Script name: RCODE.1st_CRAFT_seq.genotyping.r
##
## Original file: 250920.genotyping.from.gDNA.log.r
##
## Purpose of script: Do genotyping from gDNA sequencing data from the validation CRAFT-seq experiment (KO/pQTL)
##
## Date Created: 2025-09-20
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
library(tidyplots)
library(patchwork)


argv=commandArgs(T)

input_file=argv[1]

if(FALSE){
input_file="/path/to/collapsed.processed.gDNA.data.tsv"
}

data<-read_tsv(input_file)


#first plot the number of alleles per plate
allele_plt<-dplyr::select(data,Total_allele_count,DNA_Barcode,cell_name)%>%
  dplyr::distinct()
p<-ggplot(allele_plt,aes(x=log10(Total_allele_count)))+
  geom_histogram()+theme_bw()+theme(aspect.ratio=1/1)+geom_vline(xintercept=4)+geom_vline(xintercept=log10(5000))+
  facet_wrap(~DNA_Barcode,ncol=4)

ggsave("Number.of.alleles.per.plates.pdf",p,width=10,height=10)

filter(data,Total_allele_count>5000)%>%pull(cell_name)%>%unique()%>%length()
filter(data,Total_allele_count>10000)%>%pull(cell_name)%>%unique()%>%length()
data%>%pull(cell_name)%>%unique()%>%length()


#check pQTL genotyping
pQTL_df<-data%>%dplyr::select(DNA_Barcode,Well_ID,Total_allele_count,cell_name,Annotation,Ratio,condition)%>%
  filter(condition%in%c("NTC","Variant"))%>%
  pivot_wider(names_from="Annotation",values_from="Ratio")%>%
  mutate(unexpected=KO.mut+VAR.mut+Other+KO)

p<-ggplot(pQTL_df,aes(x=log10(Total_allele_count),y=unexpected,color=condition))+
  geom_point()+geom_vline(xintercept=4)+geom_vline(xintercept=log10(10000))+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("NTC"="#1c39bb","Variant"="#9b111e","KO"="#e1ad01"))+
  facet_wrap(~DNA_Barcode,ncol=4)

ggsave("pQTL.condition.unexpected.alleles.pdf",p,width=10,height=10)

p<-ggplot(pQTL_df,aes(x=WT,y=VAR,color=condition))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("NTC"="#1c39bb","Variant"="#9b111e","KO"="#e1ad01"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("pQTL.condition.WT.vs.pQTL.pdf",p,width=10,height=10)

p<-ggplot(filter(pQTL_df,Total_allele_count>10000),aes(x=WT,y=VAR,color=condition))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("NTC"="#1c39bb","Variant"="#9b111e","KO"="#e1ad01"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("pQTL.condition.WT.vs.pQTL.10000.pdf",p,width=10,height=10)

#check KO genotyping
KO_df<-data%>%dplyr::select(DNA_Barcode,Well_ID,Total_allele_count,cell_name,Annotation,Ratio,condition)%>%
  filter(condition%in%c("NTC","KO"))%>%
  pivot_wider(names_from="Annotation",values_from="Ratio")%>%
  mutate(unexpected=KO.mut+VAR.mut+Other+VAR)

p<-ggplot(KO_df,aes(x=log10(Total_allele_count),y=unexpected,color=condition))+
  geom_point()+geom_vline(xintercept=4)+geom_vline(xintercept=log10(10000))+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("NTC"="#1c39bb","Variant"="#9b111e","KO"="#e1ad01"))+
  facet_wrap(~DNA_Barcode,ncol=4)

ggsave("KO.condition.unexpected.alleles.pdf",p,width=10,height=10)

p<-ggplot(KO_df,aes(x=log10(Total_allele_count),y=unexpected,color=condition))+
  geom_point()+geom_vline(xintercept=4)+geom_vline(xintercept=log10(10000))+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("NTC"="#1c39bb","Variant"="#9b111e","KO"="#e1ad01"))+
  facet_wrap(~DNA_Barcode,ncol=4)

ggsave("KO.condition.WT.vs.KO.pdf",p,width=10,height=10)


p<-ggplot(filter(KO_df,Total_allele_count>10000),aes(x=WT,y=KO,color=condition))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("NTC"="#1c39bb","Variant"="#9b111e","KO"="#e1ad01"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("KO.condition.WT.vs.KO.10000.pdf",p,width=10,height=10)


#Use 10000 as a threshold
data2<-filter(data,Total_allele_count>10000)
pQTL_df2<-filter(pQTL_df,Total_allele_count>10000)
KO_df2<-filter(KO_df,Total_allele_count>10000)




#remove unexpected allele high cells
IQR_QC<-function(x){
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  IQR <- Q3 - Q1
  lower_bound <- Q1 - 3 * IQR #we do not need lower bound
  upper_bound <- Q3 + 3 * IQR
  outliers <- (x > upper_bound)
  return(outliers)
}

pQTL_df2 <- pQTL_df2 %>%
  group_by(DNA_Barcode) %>%
  mutate(unexpected_QC = IQR_QC(unexpected)) %>%
  ungroup()

p<-ggplot(pQTL_df2,aes(x=WT,y=VAR,color=unexpected_QC))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("FALSE"="#01016F","TRUE"="#D8031C"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("pQTL.condition.unexpected.allele.QC.pdf",p,width=10,height=10)



KO_df2 <- KO_df2 %>%
  group_by(DNA_Barcode) %>%
  mutate(unexpected_QC = IQR_QC(unexpected)) %>%
  ungroup()

p<-ggplot(KO_df2,aes(x=WT,y=KO,color=unexpected_QC))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("FALSE"="#01016F","TRUE"="#D8031C"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("KO.condition.unexpected.allele.QC.pdf",p,width=10,height=10)



pQTL_df3<-filter(pQTL_df2,unexpected_QC==FALSE)
KO_df3<-filter(KO_df2,unexpected_QC==FALSE)



#Do clustering using remained cells
#Do clustering at a same time for all plates after normalization
pQTL_genotype_kmeans_norm <- pQTL_df3 %>%
  group_by(DNA_Barcode) %>%
  mutate(min_VAR = min(VAR),max_VAR = max(VAR),min_WT = min(WT),max_WT = max(WT)) %>%
  mutate(VAR_norm = (VAR - min_VAR) / (max_VAR - min_VAR),
  WT_norm = (WT - min_WT) / (max_WT - min_WT)) %>%
  ungroup() %>%
  select(-min_VAR, -max_VAR,-min_WT,-max_WT)

clst<-kmeans(pQTL_genotype_kmeans_norm%>%dplyr::select(VAR_norm,WT_norm),centers=3)
pQTL_genotype_kmeans_norm$clst<-clst$cluster

pQTL_allele_table<-pQTL_genotype_kmeans_norm%>%group_by(clst)%>%summarize(WT_mean=mean(WT_norm),VAR_mean=mean(VAR_norm))%>%
  arrange(VAR_mean)%>%mutate(allele_dosage=row_number()-1)%>%dplyr::select(clst,allele_dosage)%>%
  mutate(allele_dosage=as.factor(allele_dosage))
pQTL_genotype_kmeans_norm<-left_join(pQTL_genotype_kmeans_norm,pQTL_allele_table,by="clst")

p<-ggplot(pQTL_genotype_kmeans_norm,aes(x=WT,y=VAR,color=allele_dosage))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("pQTL.condition.kmeans.after.min.max.pdf",p,width=10,height=10)



KO_genotype_kmeans_norm <- KO_df3 %>%
  group_by(DNA_Barcode) %>%
  mutate(min_KO = min(KO),max_KO = max(KO),min_WT = min(WT),max_WT = max(WT)) %>%
  mutate(KO_norm = (KO - min_KO) / (max_KO - min_KO),
  WT_norm = (WT - min_WT) / (max_WT - min_WT)) %>%
  ungroup() %>%
  select(-min_KO, -max_KO,-min_WT,-max_WT)

clst<-kmeans(KO_genotype_kmeans_norm%>%dplyr::select(KO_norm,WT_norm),centers=3)
KO_genotype_kmeans_norm$clst<-clst$cluster


KO_allele_table<-KO_genotype_kmeans_norm%>%group_by(clst)%>%summarize(WT_mean=mean(WT_norm),KO_mean=mean(KO_norm))%>%
  arrange(KO_mean)%>%mutate(allele_dosage=row_number()-1)%>%dplyr::select(clst,allele_dosage)%>%
  mutate(allele_dosage=as.factor(allele_dosage))
KO_genotype_kmeans_norm<-left_join(KO_genotype_kmeans_norm,KO_allele_table,by="clst")

p<-ggplot(KO_genotype_kmeans_norm,aes(x=WT,y=KO,color=allele_dosage))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("KO.condition.kmeans.after.min.max.pdf",p,width=10,height=10)


pQTL_genotype_kmeans_norm%>%write_tsv("pQTL_genotype_kmeans.after.norm.tsv.gz")
KO_genotype_kmeans_norm%>%write_tsv("KO_genotype_kmeans.after.norm.tsv.gz")


#quickly check the association to CD40-PE
CD40_df<-dplyr::select(data,DNA_Barcode,Well_ID,CD40)%>%dplyr::distinct()
inner_join(pQTL_genotype_kmeans_norm,CD40_df,by=c("DNA_Barcode","Well_ID"))%>%write_tsv("pQTL_genotype_kmeans.after.norm.w.CD40.tsv.gz")
inner_join(KO_genotype_kmeans_norm,CD40_df,by=c("DNA_Barcode","Well_ID"))%>%write_tsv("KO_genotype_kmeans.after.norm.w.CD40.tsv.gz")

CD40_df<-dplyr::select(data,cell_name,CD40)%>%dplyr::distinct()
plot_CD40<-function(df,name){

	pQTL_test<-inner_join(df,CD40_df,by="cell_name")%>%
	mutate(logCD40=log10(CD40))%>%
	mutate(allele_category=ifelse(condition=="NTC","NTC",as.character(allele_dosage)))
	pQTL_test$allele_category<-factor(pQTL_test$allele_category, levels=c("NTC","0","1","2"))

	p<-pQTL_test |>
	tidyplot(x = allele_category, y = logCD40, color = allele_category) |>
	add_boxplot() |>
	add_test_pvalue(ref.group = "0")

	ggsave(str_c(name,".CD40.tidyplot.pdf"),p,width=4,height=4)



	p<-ggplot(pQTL_test,aes(x=allele_category,y=logCD40))+
	geom_boxplot(outlier.shape=NA,colour="black",aes(fill=allele_category))+
	geom_jitter(alpha=0.5,size=0.1,height=0,width=0.3,colour="#000000",shape=21)+
	theme_bw()+
	theme(axis.text.x = element_text(size = 10),
			axis.text.y = element_text(size = 10),aspect.ratio=2/3,legend.position="none")+
	scale_fill_manual(values = c("NTC" = "#BCE1DF","0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
	xlab("")+ylab("log10(CD40)")

	ggsave(str_c(name,".CD40.boxplot.pdf"),p,width=4,height=4)

	p<-p+facet_wrap(~DNA_Barcode,ncol=4,scales="free")
	ggsave(str_c(name,".CD40.boxplot.per.plate.pdf"),p,width=10,height=10)



	p <- ggplot(pQTL_test, aes(x = allele_category, y = logCD40)) +
	geom_violin(aes(fill = allele_category), colour = "black", trim = FALSE) +
	geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", colour = "black") +
	theme_bw() +
	theme(
		axis.text.x = element_text(size = 10),
		axis.text.y = element_text(size = 10),
		aspect.ratio = 2/3,
		legend.position = "none") +
	scale_fill_manual(values = c("NTC" = "#BCE1DF","0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
	xlab("") + ylab("log10(CD40)")

	ggsave(str_c(name,".CD40.violinplot.pdf"),p,width=4,height=4)

	p<-p+facet_wrap(~DNA_Barcode,ncol=4,scales="free")
	ggsave(str_c(name,".CD40.violinplot.per.plate.pdf"),p,width=10,height=10)

	pQTL_test%>%group_by(condition,allele_dosage)%>%summarize(Count=n())%>%write_tsv(str_c(name,".num.cell.per.genotypes.tsv"))

	lm_df<-pQTL_test%>%mutate(is_edited=ifelse(condition=="NTC",0,1))
	model <- lm("logCD40 ~ is_edited + DNA_Barcode", data = lm_df)
	edited<-summary(model)$coefficients["is_edited",]

	lm_df<-filter(pQTL_test, condition!="NTC")%>%mutate(allele_num=as.numeric(allele_dosage))
	model <- lm("logCD40 ~ allele_num + DNA_Barcode", data = lm_df)
	allele_dosage_012<-summary(model)$coefficients["allele_num",]

	lm_df<-filter(pQTL_test, condition!="NTC", allele_dosage %in% c("0","1"))%>%mutate(allele_num=ifelse(allele_dosage=="1",1,0))
	model <- lm("logCD40 ~ allele_num + DNA_Barcode", data = lm_df)
	allele_dosage_01<-summary(model)$coefficients["allele_num",]

	lm_df<-filter(pQTL_test, condition!="NTC", allele_dosage %in% c("0","2"))%>%mutate(allele_num=ifelse(allele_dosage=="2",1,0))
	model <- lm("logCD40 ~ allele_num + DNA_Barcode", data = lm_df)
	allele_dosage_02<-summary(model)$coefficients["allele_num",]

	lm_df<-filter(pQTL_test, condition!="NTC", allele_dosage %in% c("1","2"))%>%mutate(allele_num=ifelse(allele_dosage=="2",1,0))
	model <- lm("logCD40 ~ allele_num + DNA_Barcode", data = lm_df)
	allele_dosage_12<-summary(model)$coefficients["allele_num",]


	summary<-bind_cols(tibble(test=c("is_edited","012","01","02","12")),bind_rows(edited,allele_dosage_012,allele_dosage_01,allele_dosage_02,allele_dosage_12))
	write_tsv(summary,str_c(name,".linear.regression.tsv"))

}


#summarize normalized clustering results
plot_CD40(pQTL_genotype_kmeans_norm, "pQTL_genotype_kmeans_norm")
plot_CD40(KO_genotype_kmeans_norm, "KO_genotype_kmeans_norm")
