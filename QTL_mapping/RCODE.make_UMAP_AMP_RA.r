## ---------------------------
##
## Script name: RCODE.make_UMAP_AMP_RA.r
##
## Original file: RCODE.260312.make.UMAP.for.AMP.RA.log.r
##
## Purpose of script: Generate UMAP visualization of protein/RNA expression in AMP-RA dataset
##
##
## Date Created: 2026-03-12
##
## Email:
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


library(Seurat)
library(tidyverse)
library(ggrastr)
library(viridis)
library(patchwork)

argv=commandArgs(T)


protein_name=argv[1]
rna_name=argv[2]
name=argv[3]

#check expression in AMP-RA data set
cell_type_info="/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250909_UMAP_summary/AMP.RA.UMAP.cellbarcode.tsv.gz"
adt_cnt_file="/path/to/adt_raw.rds"
rna_cnt_file="/path/to/GEX_raw.rds"
protein_qc_file="/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/260304_KLD_protein_QC/AMP.RA.keep.protein.list.tsv"


adt_cnt_obj<-readRDS(adt_cnt_file)
rna_cnt_obj<-readRDS(rna_cnt_file)

qc_df<-read_tsv(protein_qc_file)
adt_cnt_obj<-adt_cnt_obj[qc_df$Gene,]

adt_data <- CreateSeuratObject(counts = adt_cnt_obj, project = "ADT", min.cells = 0, min.features = 0)
rna_data <- CreateSeuratObject(counts = rna_cnt_obj, project = "RNA", min.cells = 0, min.features = 0)

meta<-read_tsv(cell_type_info)

cells<-intersect(meta$cell_barcode,colnames(adt_data))
target_cells<-intersect(cells,colnames(rna_data))

adt_data <- NormalizeData(subset(adt_data,cells=target_cells), normalization.method = "CLR", margin = 2)
rna_data <- NormalizeData(subset(rna_data,cells=target_cells), normalization.method = "LogNormalize", scale.factor = 10000)

adt_df <- FetchData(adt_data, vars = c(protein_name)) %>% as_tibble(rownames = "cell_barcode")
rna_df <- FetchData(rna_data, vars = c(rna_name)) %>% as_tibble(rownames = "cell_barcode")
colnames(adt_df)<-c("cell_barcode","ADT_EXP")
colnames(rna_df)<-c("cell_barcode","RNA_EXP")

df<-inner_join(meta,adt_df,by="cell_barcode")%>%inner_join(rna_df,by="cell_barcode")

#make UMAP
p_rna<-ggplot(df, aes(x = RNA_UMAP1, y = RNA_UMAP2, color = RNA_EXP)) +
geom_point_rast(size = 0.01) +
scale_color_viridis(option = "viridis") + 
theme_void() + 
theme(aspect.ratio = 1/1) +
labs(title = "RNA") +
theme(plot.title = element_text(hjust = 0.5, size = 8))

p_adt<-ggplot(df, aes(x = RNA_UMAP1, y = RNA_UMAP2, color = ADT_EXP)) +
geom_point_rast(size = 0.01) +
scale_color_viridis(option = "plasma") + 
theme_void() + 
theme(aspect.ratio = 1/1) +
labs(title = "ADT") +
theme(plot.title = element_text(hjust = 0.5, size = 8))

p<-p_rna+p_adt+plot_layout(ncol=2)
ggsave(str_c(name,".RNA.ADT.UMAP.pdf"),p,width=8,height=3)


p_rna<-ggplot(df, aes(x = RNA_UMAP1, y = RNA_UMAP2, z= RNA_EXP)) +
stat_summary_2d(fun = mean)+
scale_fill_viridis(option = "viridis") + 
theme_void() + 
theme(aspect.ratio = 1/1) +
labs(title = "RNA") +
theme(plot.title = element_text(hjust = 0.5, size = 8))

p_adt<-ggplot(df, aes(x = RNA_UMAP1, y = RNA_UMAP2, z= ADT_EXP)) +
stat_summary_2d(fun = mean)+
scale_fill_viridis(option = "plasma") + 
theme_void() + 
theme(aspect.ratio = 1/1) +
labs(title = "ADT") +
theme(plot.title = element_text(hjust = 0.5, size = 8))

p<-p_rna+p_adt+plot_layout(ncol=2)
ggsave(str_c(name,".RNA.ADT.UMAP.2d.pdf"),p,width=8,height=3)
