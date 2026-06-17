## ---------------------------
##
## Script name: RCODE.merge.primary.B.CRAFT.seq.r
##
## Original file: 251218.merge.primary.B.CRAFT.seq.log.r
##
## Purpose of script: Merge 3 CRAFT-seq runs (PE1, PE2, ABE),
##                    perform QC filtering, dimensionality reduction,
##                    and clustering
##
## Date Created: 2025-12-18
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
seurat_PE1_file=argv[1]
seurat_PE2_file=argv[2]
seurat_ABE_file=argv[3]
geno_PE1_file=argv[4]
geno_PE2_file=argv[5]
geno_ABE_file=argv[6]

if(FALSE){
seurat_PE1_file="/path/to/PE1.raw.RNA.ADT.object.rds"
seurat_PE2_file="/path/to/PE2.raw.RNA.ADT.object.rds"
seurat_ABE_file="/path/to/ABE.raw.RNA.ADT.object.rds"
geno_PE1_file="/path/to/PE1.pQTL_genotype_kmeans.after.norm.w.CD40.tsv.gz"
geno_PE2_file="/path/to/PE2.pQTL_genotype_kmeans.after.norm.w.CD40.tsv.gz"
geno_ABE_file="/path/to/ABE.pQTL_genotype_kmeans.after.norm.w.CD40.tsv.gz"
}

library(tidyverse)
library(Seurat)
library(harmony)
library(ggpointdensity)
library(viridis)
library(patchwork)

rna_thr=1000
adt_thr=50
mt_thr=20


file_list<-c(seurat_PE1_file, seurat_PE2_file, seurat_ABE_file)

obj_list <- lapply(file_list, readRDS)

seurat_merged <- merge(
  x = obj_list[[1]],
  y = obj_list[2:3],
  add.cell.ids = c("PE1", "PE2", "ABE")
)

seurat_merged <- JoinLayers(seurat_merged, assay = "RNA")
seurat_merged <- JoinLayers(seurat_merged, assay = "ADT")

#add experiment batch information to the object
meta<-seurat_merged@meta.data
meta$cell<-rownames(meta)
meta<-mutate(meta,experiment=case_when(
	str_detect(cell,"^PE1")~"PE1",
	str_detect(cell,"^PE2")~"PE2",
	str_detect(cell,"^ABE")~"ABE"
	))%>%mutate(plate_exp=str_c(plate,experiment))

seurat_merged@meta.data<-meta
saveRDS(seurat_merged,"raw.merge.primaryB.CRAFT.rds")


#also add genotype information
geno1<-read_tsv(geno_PE1_file)%>%
mutate(experiment="PE1",C_allele_dosage=2-allele_dosage)
geno2<-read_tsv(geno_PE2_file)%>%
mutate(experiment="PE2",C_allele_dosage=2-allele_dosage)
geno3<-read_tsv(geno_ABE_file)%>%
mutate(experiment="ABE",C_allele_dosage=allele_dosage)

geno_df<-bind_rows(geno1,geno2,geno3)%>%mutate(cell_name=str_c(experiment,cell_name,sep="_"))%>%mutate(T_allele_dosage=2-C_allele_dosage)
cell_names<-filter(geno_df,condition=="pQTL")%>%pull(cell_name)

cells_with_geno<-seurat_merged@meta.data %>% mutate(cell_name=str_c(experiment,plate,name,sep="_"))%>%
filter(cell_name %in% cell_names)%>%rownames()

data<-subset(seurat_merged,cells=cells_with_geno)

meta<-data@meta.data
meta<-meta %>% mutate(cell_name=str_c(experiment,plate,name,sep="_"))%>%left_join(geno_df,by=c("cell_name","experiment"))%>%
rename(FACS_CD40=CD40,FACS_CD27=CD27)
rownames(meta)<-meta$cell
data@meta.data<-meta

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data<- NormalizeData(data, normalization.method = "CLR", margin = 2, assay = "ADT")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

saveRDS(data,"raw.merge.primaryB.CRAFT.w.geno.rds")




#make basic plots
#make basic plot (UMI/nGene/%MT) per experiment
plt_df<-as_tibble(data@meta.data)
plt_df$experiment<-factor(plt_df$experiment, levels=c("PE1","PE2","ABE"))
write_tsv(plt_df,"raw.metrics.df.tsv.gz")


p1<-ggplot(plt_df,aes(x=nCount_RNA))+
geom_histogram(position = "stack")+theme_bw()+theme(aspect.ratio=1/1)+
facet_wrap(~experiment,ncol=3)+scale_x_log10()

p2<-ggplot(plt_df,aes(x=nFeature_RNA))+
geom_histogram(position = "stack")+theme_bw()+theme(aspect.ratio=1/1)+
facet_wrap(~experiment,ncol=3)+scale_x_log10()+geom_vline(xintercept=rna_thr)

p3<-ggplot(plt_df,aes(x=percent.mt))+
geom_histogram(position = "stack")+theme_bw()+theme(aspect.ratio=1/1)+geom_vline(xintercept=mt_thr)+
facet_wrap(~experiment,ncol=3)

p4<-ggplot(plt_df,aes(x=nCount_ADT))+
geom_histogram(position = "stack")+theme_bw()+theme(aspect.ratio=1/1)+
facet_wrap(~experiment,ncol=3)+scale_x_log10()

p5<-ggplot(plt_df,aes(x=nFeature_ADT))+
geom_histogram(position = "stack")+theme_bw()+theme(aspect.ratio=1/1)+
facet_wrap(~experiment,ncol=3)+scale_x_log10()+geom_vline(xintercept=adt_thr)

p<-p1+p2+p4+p5+p3+plot_layout(ncol=2)
ggsave("merged.object.metrics.pdf",p,width=10,height=6)



p<-ggplot(plt_df,aes(x=nFeature_RNA,y=nFeature_ADT))+
geom_pointdensity()+theme_bw()+theme(aspect.ratio=1/1)+scale_x_log10()+
scale_y_log10()+scale_color_viridis()+
facet_wrap(~experiment,ncol=3)+geom_vline(xintercept=rna_thr)+geom_hline(yintercept=adt_thr)

ggsave("nFeature_RNA.vs.nFeature_ADT.pdf",p,width=9,height=3)



#filter out low quality cells
cell_ids <- plt_df %>% filter(nFeature_RNA > rna_thr & nFeature_ADT > adt_thr & percent.mt <mt_thr) %>% pull(cell)
data2 <- data[, cell_ids]


data2 <- FindVariableFeatures(data2, selection.method = "vst", nfeatures = 1000)
hvg <- VariableFeatures(data2)
var_regex = '^IG[HJKL]' # remove immunoglobulin genes
hvg = grep(var_regex, hvg, invert=T, value=T)

data2 <- ScaleData(data2, features = rownames(data2))
data2 <- RunPCA(data2, features = hvg)

sdev <- Stdev(data2, reduction = "pca")
explained_variance_ratio <- sdev^2 / sum(sdev^2)
cumulative_explained <- cumsum(explained_variance_ratio)

variance_df <- tibble(
  PC = paste0("PC", seq_along(sdev)),
  StdDev = sdev,
  ExplainedVariance = explained_variance_ratio,
  CumulativeVariance = cumulative_explained
)

variance_df$PC<-factor(variance_df$PC, levels=variance_df$PC)

p<-ggplot(head(variance_df,20),aes(x=PC,y=ExplainedVariance))+
geom_point()+theme_bw()+theme(aspect.ratio=2/3)
ggsave("pca.variance.pdf",p)


#run harmony
data2 <- RunHarmony(data2,
				group.by.vars = c("plate_exp", "experiment"),
				reduction = "pca", assay.use = "RNA", reduction.save = "harmony")

#run clustering
data2<- data2 %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = 0.2)

all_markers <- FindAllMarkers(data2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(all_markers, "find.marker.result.rds")

saveRDS(data2,"QCed.merge.primaryB.CRAFT.rds")
