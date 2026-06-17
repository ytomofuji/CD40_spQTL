## ---------------------------
##
## Script name: RCODE.stim_Daudi.genotyping.r
##
## Original file: 251019.genotyping.from.gDNA.log.r
##
## Purpose of script: Do genotyping from gDNA sequencing data from the stimulated Daudi CRAFT-seq experiment
##
## Date Created: 2025-10-19
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
  geom_histogram()+theme_bw()+theme(aspect.ratio=1/1)+geom_vline(xintercept=3)+geom_vline(xintercept=4)+
  geom_vline(xintercept=log10(5000))+
  facet_wrap(~DNA_Barcode,ncol=3)

ggsave("Number.of.alleles.per.plates.pdf",p,width=10,height=10)




#check pQTL genotyping
pQTL_df<-data%>%dplyr::select(DNA_Barcode,Well_ID,Total_allele_count,cell_name,Annotation,Ratio,stimulation,condition)%>%
  filter(condition%in%c("NTC","pQTL"))%>%
  pivot_wider(names_from="Annotation",values_from="Ratio")%>%
  mutate(unexpected=KO.mut+VAR.mut+Other+KO)

p<-ggplot(pQTL_df,aes(x=log10(Total_allele_count),y=unexpected,color=condition))+
  geom_point()+geom_vline(xintercept=4)+geom_vline(xintercept=log10(1000))+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("KO"="#4DAF4A","pQTL"="#FF7F00","NTC"="#984EA3","Unassigned"="grey"))+
  facet_wrap(~DNA_Barcode,ncol=4)

ggsave("pQTL.condition.unexpected.alleles.pdf",p,width=10,height=10)

p<-ggplot(pQTL_df,aes(x=WT,y=VAR,color=condition))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("KO"="#4DAF4A","pQTL"="#FF7F00","NTC"="#984EA3","Unassigned"="grey"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("pQTL.condition.WT.vs.pQTL.pdf",p,width=10,height=10)

p<-ggplot(filter(pQTL_df,Total_allele_count>1000),aes(x=WT,y=VAR,color=condition))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("KO"="#4DAF4A","pQTL"="#FF7F00","NTC"="#984EA3","Unassigned"="grey"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("pQTL.condition.WT.vs.pQTL.1000.pdf",p,width=10,height=10)

#check KO genotyping
KO_df<-data%>%dplyr::select(DNA_Barcode,Well_ID,Total_allele_count,cell_name,Annotation,Ratio,stimulation,condition)%>%
  filter(condition%in%c("NTC","KO"))%>%
  pivot_wider(names_from="Annotation",values_from="Ratio")%>%
  mutate(unexpected=KO.mut+VAR.mut+Other+VAR)

p<-ggplot(KO_df,aes(x=log10(Total_allele_count),y=unexpected,color=condition))+
  geom_point()+geom_vline(xintercept=4)+geom_vline(xintercept=log10(1000))+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("KO"="#4DAF4A","pQTL"="#FF7F00","NTC"="#984EA3","Unassigned"="grey"))+
  facet_wrap(~DNA_Barcode,ncol=4)

ggsave("KO.condition.unexpected.alleles.pdf",p,width=10,height=10)

p<-ggplot(KO_df,aes(x=log10(Total_allele_count),y=unexpected,color=condition))+
  geom_point()+geom_vline(xintercept=4)+geom_vline(xintercept=log10(1000))+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("KO"="#4DAF4A","pQTL"="#FF7F00","NTC"="#984EA3","Unassigned"="grey"))+
  facet_wrap(~DNA_Barcode,ncol=4)

ggsave("KO.condition.WT.vs.KO.pdf",p,width=10,height=10)


p<-ggplot(filter(KO_df,Total_allele_count>1000),aes(x=WT,y=KO,color=condition))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("KO"="#4DAF4A","pQTL"="#FF7F00","NTC"="#984EA3","Unassigned"="grey"))+
  facet_wrap(~DNA_Barcode,ncol=4,scales="free")

ggsave("KO.condition.WT.vs.KO.1000.pdf",p,width=10,height=10)


#Use 1000 as a threshold
data2<-filter(data,Total_allele_count>1000)
pQTL_df2<-filter(pQTL_df,Total_allele_count>1000)
KO_df2<-filter(KO_df,Total_allele_count>1000)




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
