## ---------------------------
##
## Script name: RCODE.Xenium.ARBITRATE.TNK.subclustering.log.r
##
## Purpose of script: Re-clustering T/NK cells and annotation of CD4T, CD8T and NK cells
##
##
## Date Created: 2026-11-13
##
## ---------------------------
##
## Notes:
##   vst HVG 2000 (T/NK cells only) -> PCA (30 PCs) -> Harmony (sid) -> uwot UMAP + fuzzy graph -> Louvain (res 0.8)
##   CD4T: cl0, 2, 3, 4, 5, 6, 8, 10 / CD8T: cl1, 9 / NK: cl7
##   output: final cell type table used for the spatial analyses
##
## ---------------------------

#cd /path/to/work_dir
#sbatch: --mem=160G --cpus-per-task=4 --time=08:00:00

library(Seurat)
library(tidyverse)
library(harmony)
library(uwot)
library(patchwork)

argv=commandArgs(T)

work_dir="/path/to/work_dir"
if(length(argv)>0){work_dir=argv[1]}

name="ARBITRATE.pre"
tnk_obj=str_c(work_dir,"/",name,".T_NK.rds")
meta_file=str_c(work_dir,"/",name,".lineage.metadata.tsv.gz")

setwd(work_dir)


#UMAP with uwot (keep fuzzy graph for clustering)
Run_uwot_umap <- function(SeuratObj, reduction = 'harmony', min_dist = 0.01, spread = 1){
    HU <- uwot::umap(SeuratObj@reductions[[reduction]]@cell.embeddings, min_dist = min_dist,
                 spread = spread, ret_extra = 'fgraph', fast_sgd = FALSE)
    colnames(HU$embedding) = c('HUMAP1', 'HUMAP2')
    rownames(HU$fgraph) = colnames(HU$fgraph) = Cells(SeuratObj)
    SeuratObj[['humap']] <- Seurat::CreateDimReducObject(embeddings = HU$embedding, assay = 'RNA', key = 'HUMAP_', global = TRUE)
    HU_graph <- Seurat::as.Graph(HU$fgraph)
    DefaultAssay(HU_graph) <- DefaultAssay(SeuratObj)
    SeuratObj[['humap_fgraph']] <- HU_graph
    return(SeuratObj)
}


data<-readRDS(tnk_obj)

#run general pipeline for clustering
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)

set.seed(0)
data <- RunPCA(data, npcs = 30)

#Batch correction
set.seed(0)
data <- RunHarmony(data, "sid", plot_convergence = FALSE, max_iter = 10, early_stop = FALSE)

#UMAP and clustering on the fuzzy graph
set.seed(0)
data <- Run_uwot_umap(data, spread = 0.8, min_dist = 0.3)
data <- FindClusters(data, graph.name = "humap_fgraph", resolution = 0.8)
print(table(data$humap_fgraph_res.0.8))

#annotation of subclusters
tnk_annotations <- c(
  "0" = "CD4T", "2" = "CD4T", "3" = "CD4T", "4" = "CD4T", "5" = "CD4T", "6" = "CD4T", "8" = "CD4T", "10" = "CD4T",
  "1" = "CD8T", "9" = "CD8T",
  "7" = "NK")
data$TNK_subcluster<-data$humap_fgraph_res.0.8
data$TNK_celltype<-unname(tnk_annotations[as.character(data$humap_fgraph_res.0.8)])
print(table(data$TNK_celltype))

#marker genes
marker_genes<-c("CD3E","CD3G","CD2","CD4","CD40LG","IL7R","CD8A","CD8B","GZMK","NKG7","GNLY","KLRF1","KLRD1","NCAM1","FCGR3A",
"FOXP3","CXCL13","PDCD1","MKI67","TRDC")
marker_genes<-intersect(marker_genes,rownames(data))

Idents(data)<-factor(data$humap_fgraph_res.0.8,levels=sort(as.integer(unique(as.character(data$humap_fgraph_res.0.8)))))
p1<-DotPlot(data, features = marker_genes, dot.scale = 6) + scale_color_distiller(palette = "Reds",direction=1) + theme_test() + RotatedAxis() + ggtitle("subclusters")
Idents(data)<-factor(data$TNK_celltype,levels=c("CD4T","CD8T","NK"))
p2<-DotPlot(data, features = marker_genes, dot.scale = 6) + scale_color_distiller(palette = "Reds",direction=1) + theme_test() + RotatedAxis() + ggtitle("cell types")
p<-p1/p2+plot_layout(heights=c(3,1))
ggsave(str_c(name,".TNK.marker.genes.dot.plot.pdf"),p,width=12,height=10)

p1<-DimPlot(data, reduction = "humap", group.by = "humap_fgraph_res.0.8", label = TRUE, raster = TRUE)+theme(aspect.ratio=1/1)
p2<-DimPlot(data, reduction = "humap", group.by = "TNK_celltype", raster = TRUE)+theme(aspect.ratio=1/1)
ggsave(str_c(name,".TNK.UMAP.pdf"),p1+p2,width=14,height=6)

#save seurat object
data[["RNA"]]$scale.data <- NULL
saveRDS(data,str_c(name,".TNK.subclst.rds"))


#final cell type table
tnk_df<-tibble(cell_id=colnames(data),TNK_subcluster=as.character(data$TNK_subcluster),TNK_celltype=data$TNK_celltype)
meta<-read_tsv(meta_file)%>%left_join(tnk_df,by="cell_id")%>%
mutate(celltype=ifelse(lineage=="T/NK",TNK_celltype,celltype))%>%
dplyr::select(-TNK_celltype)
stopifnot(!any(is.na(meta$celltype)))
write_tsv(meta,str_c(name,".celltype.tsv.gz"))

meta%>%count(celltype)%>%arrange(desc(n))%>%print(n=Inf)
meta%>%count(sid,celltype)%>%pivot_wider(names_from=celltype,values_from=n,values_fill=0)%>%
write_tsv(str_c(name,".celltype.per.sample.tsv"))
