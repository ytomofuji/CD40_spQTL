## ---------------------------
##
## Script name: RCODE.plot.UMAP.heatmap.for.CRAFT.seq.log.r
##
## Original file: RCODE.260114.plot.UMAP.heatmap.for.CRAFT.seq.log.r
##
## Purpose of script: Prepare figure for showing the overview of the CRAFT-seq dataset
##
##
## Date Created: 2026-01-14
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
library(ggpointdensity)
library(viridis)
library(patchwork)
library(ggrastr)
library(ggrepel)
library(scales)


#First, plot the clustering UMAP and cluster ratio per sample

data<-readRDS(seurat_file)
Idents(data)<-"seurat_clusters"
data <- RunUMAP(data, reduction = "harmony", dims = 1:10, reduction.name = "harmony.umap", reduction.key = "harmonyUMAP_")


df <- FetchData(data, vars = c("experiment","harmonyUMAP_1","harmonyUMAP_2","seurat_clusters","rna_CD40", "adt_CD40","C_allele_dosage"))%>%as.data.frame()%>%
rownames_to_column("barcode")%>%as_tibble()%>%
mutate(cluster=case_when(
	seurat_clusters==0 ~ "Activated Naive B cells",
	seurat_clusters==1 ~ "Plasma B cells",
	seurat_clusters==2 ~ "Intermediate B cells",
	seurat_clusters==3 ~ "NR4A+ Activated B cells"))
df$cluster<-factor(df$cluster, levels=c("Activated Naive B cells","Intermediate B cells","Plasma B cells","NR4A+ Activated B cells"))
df$experiment<-factor(df$experiment, levels=c("PE1","PE2","ABE"))

df<-mutate(df,Genotype=case_when(C_allele_dosage==0~"T/T",
	C_allele_dosage==1~"C/T",
	C_allele_dosage==2~"C/C"))
df$Genotype<-factor(df$Genotype ,levels=c("T/T","C/T","C/C"))

p<-ggplot(df,aes(x=harmonyUMAP_1,y=harmonyUMAP_2,color=cluster))+
geom_point_rast(size=0.5)+
scale_color_manual(values=c(
"Activated Naive B cells"="#D93644",
"Intermediate B cells"="#FFB84C",
"Plasma B cells"="#93C9D9",
"NR4A+ Activated B cells"="#8720D6"
))+
theme_void()+theme(aspect.ratio=1/1)

ggsave("Cluster.UMAP.pdf",p,width=8,height=5)
write_tsv(df, "original.table.tsv.gz")



p<-ggplot(df,aes(x=harmonyUMAP_1,y=harmonyUMAP_2,fill=cluster))+
geom_point_rast(size=1,color="black",shape=21,stroke=0.2)+
scale_fill_manual(values=c(
"Activated Naive B cells"="#D93644",
"Intermediate B cells"="#FFB84C",
"Plasma B cells"="#93C9D9",
"NR4A+ Activated B cells"="#8720D6"
))+
theme_void()+theme(aspect.ratio=1/1)

ggsave("Cluster.bold.UMAP.pdf",p,width=8,height=5)



#plot cluster content ratio
plt_df <- df %>% group_by(experiment,cluster) %>% summarize(Count=n())

p<-ggplot(plt_df,aes(x=experiment,y=Count,fill=cluster))+
geom_bar(stat="identity", position="fill")+
scale_y_continuous(labels = percent, expand = c(0, 0))+
scale_fill_manual(values=c(
"Activated Naive B cells"="#D93644",
"Intermediate B cells"="#FFB84C",
"Plasma B cells"="#93C9D9",
"NR4A+ Activated B cells"="#8720D6"
))+theme_classic()+theme(aspect.ratio=3/2)

ggsave("Cluster.per.experiment.bar.pdf",p,width=6,height=4)


plt_df$experiment<-factor(plt_df$experiment, levels=rev(c("PE1","PE2","ABE")))
plt_df$cluster<-factor(plt_df$cluster, levels=rev(c("Activated Naive B cells","Intermediate B cells","Plasma B cells","NR4A+ Activated B cells")))


p<-ggplot(plt_df,aes(y=experiment,x=Count,fill=cluster))+
geom_bar(stat="identity", position="stack")+
scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
geom_text(aes(label = Count),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black",
    fontface = "bold") +
scale_fill_manual(values=c(
"Activated Naive B cells"="#D93644",
"Intermediate B cells"="#FFB84C",
"Plasma B cells"="#93C9D9",
"NR4A+ Activated B cells"="#8720D6"
))+theme_test()+theme(aspect.ratio=1/2)

ggsave("Cluster.per.experiment.bar.count.pdf",p,width=6,height=4)




#plot Genotype per cluster
plt_df <- df %>% group_by(Genotype,cluster) %>% summarize(Count=n())

p<-ggplot(plt_df,aes(x=cluster,y=Count,fill=Genotype))+
geom_bar(stat="identity", position="fill")+
scale_y_continuous(labels = percent, expand = c(0, 0))+
scale_fill_manual(values = c("C/C" = "#1b9e77","C/T" = "#d95f02","T/T" = "#7570b3"))+
theme_classic()+theme(aspect.ratio=3/2,
axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))

ggsave("Genotype.per.Cluster.bar.pdf",p,width=6,height=4)



#make plot for CD40 expression (mRNA / ADT; on UMAP)
p<-ggplot(arrange(df,rna_CD40),aes(x=harmonyUMAP_1,y=harmonyUMAP_2,color=rna_CD40))+
geom_point_rast(size=0.5)+
scale_color_distiller(palette = "Reds",direction=1)+
theme_void()+theme(aspect.ratio=1/1)

ggsave("RNA.CD40.UMAP.pdf",p,width=8,height=5)


p<-ggplot(arrange(df,adt_CD40),aes(x=harmonyUMAP_1,y=harmonyUMAP_2,color=adt_CD40))+
geom_point_rast(size=0.5)+
scale_color_distiller(palette = "Reds",direction=1)+
theme_void()+theme(aspect.ratio=1/1)

ggsave("ADT.CD40.UMAP.pdf",p,width=8,height=5)



#make heatmap for marker genes
# Naive-like marker genes: TCL1A, FCER2, IL4R, IGHM, IGHD
# GC-like marker genes: CD80, CD86, TRAF1, EBI3, MIR155HG, NFKB1, RELB
# Plasma-like genes: MZB1, XBP1, PRDM1
# Activated genes: NR4A1, NR4A3, TNFRSF18
# CD40
markers <- c("TCL1A", "FCER2", "IL4R", "MS4A1",
             "CD80", "CD86", "TRAF1", "EBI3", "MIR155HG", "NFKB1",
             "CD27", "CD38", "SDC1", "MZB1", "XBP1", "PRDM1", "IRF4",
			 "NR4A1", "NR4A3", "TNFRSF18","ITGAX","TBX21","ZEB2",
             "CD40")


gex_data <- AverageExpression(
  object = data,
  assays = "RNA",
  slot = "data",
  group.by = "seurat_clusters")

plt_df<-apply(gex_data$RNA[markers,],1,scale)%>%as.data.frame()%>%t()%>%as.data.frame()%>%
rownames_to_column("Gene")%>%as_tibble()
colnames(plt_df)<-c("Gene","Activated Naive B cells","Plasma B cells","Intermediate B cells","NR4A+ Activated B cells")
plt_df<-pivot_longer(plt_df,names_to="cell_type",values_to="Z_exp",-Gene)
plt_df$cell_type<-factor(plt_df$cell_type, levels=rev(c("Activated Naive B cells","Intermediate B cells","Plasma B cells","NR4A+ Activated B cells")))
plt_df$Gene<-factor(plt_df$Gene, levels=markers)


range=max(abs(plt_df$Z_exp))*1.01

p1<-ggplot(plt_df,aes(x=Gene,y=cell_type,fill=Z_exp))+
geom_tile(colour="#666666")+theme_test()+
scale_fill_distiller(palette = 'RdYlBu',limits = c(-range,range)) +
scale_y_discrete(drop=FALSE)+
theme(
axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
axis.text.y = element_text(size = 10))+coord_equal()




#Also make ADT heatmap
adt_markers <- c("CD19", "CD20", "CD23", "CD86",
             "CD27", "CD38",
             "CD11c","CD40")


adt_data <- AverageExpression(
  object = data,
  assays = "ADT",
  slot = "data",
  group.by = "seurat_clusters")

adt_plt_df<-apply(adt_data$ADT[adt_markers,],1,scale)%>%as.data.frame()%>%t()%>%as.data.frame()%>%
rownames_to_column("Gene")%>%as_tibble()
colnames(adt_plt_df)<-c("Gene","Activated Naive B cells","Plasma B cells","Intermediate B cells","NR4A+ Activated B cells")
adt_plt_df<-pivot_longer(adt_plt_df,names_to="cell_type",values_to="Z_exp",-Gene)
adt_plt_df$cell_type<-factor(adt_plt_df$cell_type, levels=rev(c("Activated Naive B cells","Intermediate B cells","Plasma B cells","NR4A+ Activated B cells")))
adt_plt_df$Gene<-factor(adt_plt_df$Gene, levels=adt_markers)


adt_range=max(abs(adt_plt_df$Z_exp))*1.01

p2<-ggplot(adt_plt_df,aes(x=Gene,y=cell_type,fill=Z_exp))+
geom_tile(colour="#666666")+theme_test()+
scale_fill_distiller(palette = 'RdYlBu',limits = c(-adt_range,adt_range)) +
scale_y_discrete(drop=FALSE)+
theme(
axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
axis.text.y = element_text(size = 10))+coord_equal()

p<-p1+p2+plot_layout(ncol=2)

ggsave("Heatmap.marker.mRNA.ADT.for.fig.pdf",p,width=10,height=5)


saveRDS(data,"CRAFT.seq.with.new.UMAP.rds")
