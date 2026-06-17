## ---------------------------
##
## Script name: RCODE.2nd_CRAFT_seq.E123.genotyping.r
##
## Original file: 250806.genotyping.from.gDNA.E123.log.r
##
## Purpose of script: Do genotyping from gDNA sequencing data from the 2nd CRAFT-seq experiment
##                    (finemapping: E1/E2/E3 enhancer conditions)
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


library(tidyverse)
library(tidyplots)
library(patchwork)


argv=commandArgs(T)

input_file=argv[1]
name=argv[2]

if(FALSE){
input_file="/path/to/E2.collapsed.processed.gDNA.data.tsv"
name="E2"
}

data<-read_tsv(input_file)


#first plot the number of alleles per plate
allele_plt<-dplyr::select(data,Total_allele_count,DNA_Barcode,cell_name)%>%
  dplyr::distinct()
p<-ggplot(allele_plt,aes(x=log10(Total_allele_count)))+
  geom_histogram()+theme_bw()+theme(aspect.ratio=1/1)+geom_vline(xintercept=3)+geom_vline(xintercept=4)+
  geom_vline(xintercept=log10(5000))+
  facet_wrap(~DNA_Barcode,ncol=4)

ggsave(str_c(name,".Number.of.alleles.per.plates.pdf"),p,width=10,height=10)



#check variant genotyping
variant_df<-data%>%dplyr::select(DNA_Barcode,Well_ID,Total_allele_count,cell_name,Annotation,Ratio,condition)%>%
  filter(condition!="Unassigned")%>%
  pivot_wider(names_from="Annotation",values_from="Ratio")%>%
  mutate(unexpected=VAR.mut+Other)

p<-ggplot(variant_df,aes(x=log10(Total_allele_count),y=unexpected,color=condition))+
  geom_point()+geom_vline(xintercept=4)+geom_vline(xintercept=log10(1000))+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("Enh1"="#4DAF4A","Enh2"="#4DAF4A","Enh3"="#4DAF4A","NTC"="#984EA3","Unassigned"="grey"))+
  facet_wrap(~DNA_Barcode,ncol=3)

ggsave(str_c(name,".condition.unexpected.alleles.pdf"),p,width=10,height=10)

p<-ggplot(variant_df,aes(x=WT,y=VAR,color=condition))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("Enh1"="#4DAF4A","Enh2"="#4DAF4A","Enh3"="#4DAF4A","NTC"="#984EA3","Unassigned"="grey"))+
  facet_wrap(~DNA_Barcode,ncol=3,scales="free")

ggsave(str_c(name,".condition.WT.vs.variant.pdf"),p,width=10,height=10)

p<-ggplot(filter(variant_df,Total_allele_count>1000),aes(x=WT,y=VAR,color=condition))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("Enh1"="#4DAF4A","Enh2"="#4DAF4A","Enh3"="#4DAF4A","NTC"="#984EA3","Unassigned"="grey"))+
  facet_wrap(~DNA_Barcode,ncol=3,scales="free")

ggsave(str_c(name,".condition.WT.vs.variant.1000.pdf"),p,width=10,height=10)


variant_df2<-filter(variant_df,Total_allele_count>1000)




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

variant_df2 <- variant_df2 %>%
  group_by(DNA_Barcode) %>%
  mutate(unexpected_QC = IQR_QC(unexpected)) %>%
  ungroup()

p<-ggplot(variant_df2,aes(x=WT,y=VAR,color=unexpected_QC))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values=c("FALSE"="#01016F","TRUE"="#D8031C"))+
  facet_wrap(~DNA_Barcode,ncol=3,scales="free")

ggsave(str_c(name,".condition.unexpected.allele.QC.pdf"),p,width=10,height=10)




variant_df3<-filter(variant_df2,unexpected_QC==FALSE)



#Do clustering at a same time for all plates after normalization
variant_genotype_kmeans_norm <- variant_df3 %>%
  group_by(DNA_Barcode) %>%
  mutate(min_VAR = min(VAR),max_VAR = max(VAR),min_WT = min(WT),max_WT = max(WT)) %>%
  mutate(VAR_norm = (VAR - min_VAR) / (max_VAR - min_VAR),
  WT_norm = (WT - min_WT) / (max_WT - min_WT)) %>%
  ungroup() %>%
  select(-min_VAR, -max_VAR,-min_WT,-max_WT)

clst<-kmeans(variant_genotype_kmeans_norm%>%dplyr::select(VAR_norm,WT_norm),centers=3)
variant_genotype_kmeans_norm$clst<-clst$cluster

variant_allele_table<-variant_genotype_kmeans_norm%>%group_by(clst)%>%summarize(WT_mean=mean(WT_norm),VAR_mean=mean(VAR_norm))%>%
  arrange(VAR_mean)%>%mutate(allele_dosage=row_number()-1)%>%dplyr::select(clst,allele_dosage)%>%
  mutate(allele_dosage=as.factor(allele_dosage))
variant_genotype_kmeans_norm<-left_join(variant_genotype_kmeans_norm,variant_allele_table,by="clst")

p<-ggplot(variant_genotype_kmeans_norm,aes(x=WT,y=VAR,color=allele_dosage))+
  geom_point(size=0.2)+theme_bw()+theme(aspect.ratio=1/1)+
  scale_color_manual(values = c("0" = "#1b9e77","1" = "#d95f02","2" = "#7570b3"))+
  facet_wrap(~DNA_Barcode,ncol=3,scales="free")

ggsave(str_c(name,".condition.kmeans.after.min.max.pdf"),p,width=10,height=10)



variant_genotype_kmeans_norm%>%write_tsv(str_c(name,".genotype_kmeans.after.norm.tsv.gz"))


cd40_df<-dplyr::select(data,DNA_Barcode,Well_ID,CD40)%>%dplyr::distinct()
inner_join(variant_genotype_kmeans_norm,cd40_df,by=c("DNA_Barcode","Well_ID"))%>%write_tsv(str_c(name,".genotype_kmeans.after.norm.w.CD40.tsv.gz"))
