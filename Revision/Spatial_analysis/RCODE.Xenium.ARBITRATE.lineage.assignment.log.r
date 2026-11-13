## ---------------------------
##
## Script name: RCODE.Xenium.ARBITRATE.lineage.assignment.log.r
##
## Purpose of script: Assign Xenium clusters to lineages and remove low quality cells
##
##
## Date Created: 2026-11-13
##
## ---------------------------
##
## Notes:
##   clusters were assigned based on canonical marker genes and overlap of top 50 marker genes with AMPp2 clusters
##   cells with <= 80% of fuzzy graph neighbors (edge-weighted) from the same lineage -> Low quality
##   B cells = cluster 16, plasma cells = cluster 7
##
## ---------------------------

#cd /path/to/work_dir
#sbatch: --mem=160G --time=02:00:00

library(Seurat)
library(tidyverse)
library(Matrix)

argv=commandArgs(T)

work_dir="/path/to/work_dir"
if(length(argv)>0){work_dir=argv[1]}

name="ARBITRATE.pre"
seurat_file=str_c(work_dir,"/",name,".QC.harmony.clst.rds")
if(length(argv)>1){seurat_file=argv[2]}
frac_thr=0.8

setwd(work_dir)


cluster_annotations <- c(
  "0"  = "Stromal",        # PRG4+ lining
  "1"  = "Stromal",        # CXCL12+ sublining
  "2"  = "Myeloid",
  "3"  = "Mural",
  "4"  = "T/NK",
  "5"  = "Stromal",        # PRG4+ lining
  "6"  = "Stromal",        # sublining
  "7"  = "B/plasma",       # plasma cells
  "8"  = "Myeloid",
  "9"  = "Endothelial",
  "10" = "Myeloid",
  "11" = "Endothelial",
  "12" = "Stromal",
  "13" = "Stromal",        # lining, low nCount
  "14" = "Myeloid",
  "15" = "Myeloid",
  "16" = "B/plasma",       # B cells
  "17" = "Adipocyte",
  "18" = "Mural",
  "19" = "Proliferating",
  "20" = "Neuronal/glial", # MPZ/SOX10/NGFR, low nCount
  "21" = "Mast",
  "22" = "Myeloid",        # DC
  "23" = "Stromal",        # IFN response + stromal markers
  "24" = "Stromal",
  "25" = "Stromal"
)


data<-readRDS(seurat_file)
data$raw_lineage<-unname(cluster_annotations[as.character(data$seurat_clusters)])

#neighbor agreement on the fuzzy graph
graph<-data[["humap_fgraph"]]
lineage<-factor(data$raw_lineage)
one_hot<-sparse.model.matrix(~ lineage - 1, data = data.frame(lineage = lineage))
colnames(one_hot)<-levels(lineage)
nb<-as.matrix(graph %*% one_hot)
same<-nb[cbind(seq_len(nrow(nb)),match(data$raw_lineage,colnames(nb)))]
data$lineage_frac_same<-same/rowSums(nb)
data$lineage<-ifelse(!is.na(data$lineage_frac_same) & data$lineage_frac_same > frac_thr, data$raw_lineage, "Low quality")
rm(graph,nb,one_hot); gc(); gc();

print(sum(data$lineage!="Low quality"))
print(table(data$lineage))

qc_df<-tibble(seurat_clusters=data$seurat_clusters,raw_lineage=data$raw_lineage,pass=data$lineage!="Low quality")%>%
group_by(seurat_clusters,raw_lineage)%>%summarise(n_cells=n(),n_pass=sum(pass),pct_pass=100*mean(pass),.groups="drop")
write_tsv(qc_df,str_c(name,".lineage.QC.by.cluster.tsv"))

#B and plasma cells
data$celltype<-case_when(
  data$lineage=="B/plasma" & data$seurat_clusters==16 ~ "B",
  data$lineage=="B/plasma" & data$seurat_clusters==7 ~ "Plasma",
  TRUE ~ data$lineage)

#output metadata file
data@meta.data%>%rownames_to_column("cell_id")%>%
dplyr::select(cell_id,sid,nCount_RNA,nFeature_RNA,seurat_clusters,raw_lineage,lineage,lineage_frac_same,celltype)%>%
write_tsv(str_c(name,".lineage.metadata.tsv.gz"))

#save B/plasma and T/NK objects for downstream analyses
data[["humap_fgraph"]]<-NULL
saveRDS(subset(data,cells=colnames(data)[data$lineage=="B/plasma"]),str_c(name,".B_plasma.rds"))
saveRDS(subset(data,cells=colnames(data)[data$lineage=="T/NK"]),str_c(name,".T_NK.rds"))
