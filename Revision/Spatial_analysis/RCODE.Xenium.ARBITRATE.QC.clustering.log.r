## ---------------------------
##
## Script name: RCODE.Xenium.ARBITRATE.QC.clustering.log.r
##
## Purpose of script: QC and clustering of Xenium data (ARBITRATE, pre-treatment samples)
##
##
## Date Created: 2026-11-13
##
## ---------------------------
##
## Notes:
##   cell QC: nCount_RNA >= 70, nFeature_RNA >= 50, avg_assignment_confidence (Baysor) >= 0.5
##   sample QC: >= 500 cells
##   LogNormalize (scale factor = median nCount_RNA) -> vst HVG 3000 -> PCA (50 PCs) -> Harmony (sid)
##   -> uwot UMAP + fuzzy graph -> Louvain (res 1.5)
##   marker genes of each cluster are compared with AMP phase 2 clusters (top 50 markers)
##
## ---------------------------

#mkdir -p /path/to/work_dir
#cd /path/to/work_dir
#sbatch: --mem=450G --cpus-per-task=16 --time=24:00:00

library(Seurat)
library(tidyverse)
library(Matrix)
library(harmony)
library(uwot)
library(presto)
library(viridis)

argv=commandArgs(T)

work_dir="/path/to/work_dir"
if(length(argv)>0){work_dir=argv[1]}

xenium_obj="/path/to/allcells.rds"
ampp2_obj="/path/to/AMPp2_seuratObj_allcells.202511.rds"
name="ARBITRATE.pre"

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


#load data and select ARBITRATE pre-treatment samples
data<-readRDS(xenium_obj)
data<-subset(data,cells=colnames(data)[data$cohort=="ARBITRATE" & grepl("pre_treatment",data$sid)])
gc(); gc();

#cell QC
data<-subset(data, subset = nCount_RNA >= 70 & nFeature_RNA >= 50 & avg_assignment_confidence >= 0.5)

#sample QC
keep_sid<-data@meta.data%>%count(sid)%>%filter(n>=500)%>%pull(sid)
data<-subset(data, subset = sid %in% keep_sid)
print(ncol(data))
print(table(data$sid))


counts<-LayerData(data, assay = "RNA", layer = "counts")
genes<-rownames(counts)
hit<-genes %in% names(gene_aliases)
genes[hit]<-gene_aliases[genes[hit]]
rownames(counts)<-genes
data<-CreateSeuratObject(counts = counts, meta.data = data@meta.data)
rm(counts); gc(); gc();


#run general pipeline for clustering
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = median(data$nCount_RNA))
print(median(data$nCount_RNA))
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 3000)
data <- ScaleData(data)

set.seed(0)
data <- RunPCA(data)

#Batch correction
set.seed(0)
data <- RunHarmony(data, "sid", plot_convergence = FALSE, max_iter = 10, early_stop = FALSE)

#UMAP and clustering on the fuzzy graph
set.seed(0)
data <- Run_uwot_umap(data, spread = 0.8, min_dist = 0.3)
data <- FindClusters(data, graph.name = "humap_fgraph", resolution = 1.5)
data[["RNA"]]$scale.data <- NULL
print(table(data$seurat_clusters))

#save seurat object
saveRDS(data,str_c(name,".QC.harmony.clst.rds"))
writeLines(VariableFeatures(data),str_c(name,".HVG3000.txt"))


#canonical marker genes per cluster
marker_list <- list(
  t = c("CD3E", "CD8A", "CD4", "IL7R", "NKG7", "NCAM1"),
  b.plasma = c("MS4A1", "CD19", "CD79A", "PRDM1", "MZB1", "XBP1"),
  myeloid = c("CD14", "CD68", "FCGR3A", "LYVE1", "C1QA", "MARCO"),
  mast = c("HDC", "KIT"),
  endothelial = c("PECAM1", "VWF", "CDH5", "FLT1"),
  stromal = c("PRG4", "CD55", "COL6A1", "CXCL12", "PDPN", "THY1", "PDGFRB", "ACTA2"),
  adipocytes = c("ADIPOQ", "PLIN1", "PLIN4"),
  proliferating = c("MKI67", "TOP2A"),
  muscles = c("MYOD1", "CD82"),
  neuronal = c("MPZ", "SOX10", "NGFR"),
  osteoclast = c("ACP5", "CTSK", "MMP9")
)
marker_genes<-intersect(unlist(marker_list),rownames(data))

Idents(data)<-data$seurat_clusters
p<-DotPlot(data, features = marker_genes, dot.scale = 6) +
scale_color_distiller(palette = "Reds",direction=1) + theme_test() + RotatedAxis()
ggsave(str_c(name,".canonical.marker.dot.plot.by.cluster.pdf"),p,width=20,height=10)

umap_df<-bind_cols(data@meta.data%>%dplyr::select(sid,seurat_clusters),as_tibble(Embeddings(data,"humap")))
p<-ggplot(umap_df,aes(x=HUMAP_1,y=HUMAP_2,color=seurat_clusters))+
ggrastr::rasterise(geom_point(size=0.05),dpi=300)+
theme_classic()+theme(aspect.ratio=1/1)+guides(color=guide_legend(override.aes=list(size=3)))
ggsave(str_c(name,".UMAP.by.cluster.pdf"),p,width=8,height=7)
