## ---------------------------
##
## Script name: RCODE.Daudi_bulk_DEseq2_volcano_GSEA.r
##
## Original file: RCODE.260422.CD40L.Daudi.DEseq2.volcano.GSEA.log.r
##
## Purpose of script: Differential expression and GSEA analysis for
##                    CD40L-stimulated vs unstimulated Daudi cells (bulk RNA-seq)
##
##
## Date Created: 2026-04-22
##
## Email:
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

dir.create("/path/to/output_directory")
setwd("/path/to/output_directory")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(fgsea)
  library(msigdbr)
  library(tidyverse)
  library(ggrastr)
  library(ggrepel)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

raw <- read.table(
  gzfile("/path/to/Bulk.CD40L.Daudi.RNA.seq.count.annotation.tsv.gz"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)

stim_cols    <- c("stim1",    "stim2",    "stim3",    "stim4")
nonstim_cols <- c("nonstim1", "nonstim2", "nonstim3", "nonstim4")
count_cols   <- c(stim_cols, nonstim_cols)

count_mat <- as.matrix(raw[, count_cols])
rownames(count_mat) <- raw$gene_id

lib_sizes   <- colSums(count_mat)
cpm_mat_all <- sweep(count_mat, 2, lib_sizes, "/") * 1e6
keep        <- rowSums(cpm_mat_all >= 1) >= 4

count_mat_filt <- count_mat[keep, ]
storage.mode(count_mat_filt) <- "integer"

ann_filt <- raw[keep, c("gene_id", "Gene", "category")]

col_data <- data.frame(
  condition = factor(
    c(rep("stimulated", 4), rep("unstimulated", 4)),
    levels = c("unstimulated", "stimulated")
  ),
  row.names = count_cols
)

dds <- DESeqDataSetFromMatrix(
  countData = count_mat_filt,
  colData   = col_data,
  design    = ~ condition
)

dds <- DESeq(dds)
res_deseq2 <- results(dds,
                      contrast       = c("condition", "stimulated", "unstimulated"),
                      pAdjustMethod  = "BH")

deseq2_df <- as.data.frame(res_deseq2) %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(ann_filt, by = "gene_id") %>%
  arrange(padj)

deseq2_df<-filter(deseq2_df,!is.na(padj))%>%rename(FDR=padj)

write_tsv(deseq2_df, "DESeq2_results.tsv.gz")

make_volcano<-function(df,p_val_col,gene_list,FDR_THR){
	tmp<-dplyr::select(df,all_of(c("Gene","log2FoldChange",p_val_col)))
	colnames(tmp)<-c("Gene","BETA","pvalue")
	tmp$FDR<-p.adjust(tmp$pvalue, method="BH")
	val1<-filter(tmp,FDR<FDR_THR)%>%pull(pvalue)%>%max()
	val2<-filter(tmp,FDR>=FDR_THR)%>%pull(pvalue)%>%min()
	fdr_thr<-(log10(val1)+log10(val2))*(-0.5)

	tmp<-tmp%>%mutate(highlight=ifelse(Gene %in% gene_list,TRUE,FALSE))%>%
	arrange(pvalue)%>%mutate(category=case_when(
		FDR<FDR_THR&BETA>0 ~ "pos",
		FDR<FDR_THR&BETA<0 ~ "neg",
		TRUE~"others"
	))

	top_genes<-head(tmp,10)%>%pull(Gene)

	pos_cnt<-filter(tmp,category=="pos")%>%nrow()
	neg_cnt<-filter(tmp,category=="neg")%>%nrow()

	p<-ggplot(tmp,aes(x=BETA,
	y=ifelse(-log10(pvalue)>20,20,-log10(pvalue))))+
	geom_point_rast(data=subset(tmp,highlight==FALSE),alpha=0.8,size=1,aes(color=category))+
	geom_point(data=subset(tmp,highlight==TRUE),alpha=0.8,size=2,color="red")+
	geom_text_repel(data=subset(tmp,Gene %in% top_genes & highlight==FALSE),aes(label=Gene,color=category),alpha=0.8,size=2,color="black")+
	geom_text_repel(data=subset(tmp,highlight==TRUE),aes(label=Gene),alpha=0.8,size=2,color="#FF7F0E")+
	theme_test()+
	theme(
	axis.text.x = element_text(size = 10),
	axis.text.y = element_text(size = 10),
	aspect.ratio=1/1)+
	xlab("Effect size")+ylab("-log10(P)")+
	geom_hline(yintercept=fdr_thr,color="black")+
	geom_vline(xintercept=0,color="black",linetype="dashed")+
	scale_color_manual(values=c("pos"="#D62728","neg"="#1F77B4","others"="#555555"))+
	  geom_text(x = -Inf, y = Inf,
              label = neg_cnt,
              hjust = -0.1, vjust = 1.5, size = 3, color = "#1F77B4")+
    geom_text(x = Inf, y = Inf,
              label = pos_cnt,
              hjust = 1.1, vjust = 1.5, size = 3, color = "#D62728")+ylim(c(0,20.1))

	return(list(plot = p, df = tmp))

}

out_name="CD40L.Daudi.bulk"
obj1<-make_volcano(deseq2_df,"pvalue","",0.05)
ggsave(str_c(out_name,".pvalue.volcano.FDR.0.05.pdf"),obj1[["plot"]],width=4,height=4)

unique_genes<-deseq2_df%>%group_by(Gene)%>%summarize(Count=n())%>%filter(Count==1)%>%pull(Gene)%>%unique()
eqtl_df<-filter(deseq2_df,Gene %in% unique_genes)
C2_df <- msigdbr(species = "human", category = "C2")

gene_tbl <- bitr(eqtl_df$Gene, fromType = "SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% as_tibble()
gene_tbl<-left_join(gene_tbl,eqtl_df,by=c("SYMBOL"="Gene"))

C2_pathways = split(x = C2_df$entrez_gene, f = C2_df$gs_name)

ranks<-gene_tbl$stat
names(ranks)<-gene_tbl$ENTREZID
C2_fgseaRes <- fgsea(pathways=C2_pathways, ranks,
                  minSize  = 15,
                  maxSize  = 500)

CD40_related_pathways<-names(C2_pathways)[str_detect(names(C2_pathways),"CD40")]
as_tibble(C2_fgseaRes)%>%filter(pathway %in% CD40_related_pathways)%>%write_tsv("CD40L.Daudi.bulk.gsea.C2.CD40.related.results.tsv")
pathway_list<-as_tibble(C2_fgseaRes)%>%filter(pathway %in% CD40_related_pathways)%>%pull(pathway)
as_tibble(C2_fgseaRes)%>%arrange(padj)%>%write_tsv("CD40L.Daudi.bulk.gsea.C2.results.tsv.gz")

plt_df<-as_tibble(C2_fgseaRes)
p<-ggplot(plt_df,aes(x=NES,y=-log10(padj)))+
geom_point_rast(color="grey",size=1,shape=19)+
geom_hline(yintercept=-log10(0.05))+
geom_vline(xintercept=0,linetype="dashed")+
geom_point(data=subset(plt_df,pathway %in% c("BASSO_CD40_SIGNALING_UP","BASSO_CD40_SIGNALING_DN")),color="red",size=2,shape=19)+
geom_text_repel(data=subset(plt_df,pathway %in% c("BASSO_CD40_SIGNALING_UP","BASSO_CD40_SIGNALING_DN")),aes(label=pathway))+
coord_cartesian(ylim = c(0, NA)) +
theme_test()+theme(aspect.ratio=1/1)+ scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)))

ggsave("C2.volcano.pdf",p,width=3,height=3)

for(target_pathway in c("BASSO_CD40_SIGNALING_UP","BASSO_CD40_SIGNALING_DN")){
	padj=signif(filter(plt_df,pathway==target_pathway)%>%pull(padj),digits=2)
	NES=signif(filter(plt_df,pathway==target_pathway)%>%pull(NES),digits=2)

	pd <- plotEnrichmentData(pathway = C2_pathways[[target_pathway]], stats = ranks )
	p <- ggplot() +
	geom_line(data = pd$curve, aes(x = rank, y = ES), color = "black", linewidth = 1) +
	geom_segment(data = pd$ticks, aes(x = rank, xend = rank, y = -0.05, yend = 0.05), linewidth = 0.2) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
	theme_test() + labs( title = target_pathway, x = "Ranked List Metric", y = "Enrichment Score (ES)" ) + 
	theme( plot.title = element_text(hjust = 0.5, face = "bold"), axis.title = element_text(size = 12) ) +
	annotate(
	"text", x = Inf, y = Inf, 
	label = str_c("Padj=",padj,", NES=",NES),
	hjust = 1.1, vjust = 1.1, 
	size = 4, fontface = "italic", color = "red"
	)

	ggsave(str_c(target_pathway,".CD40L.Daudi.bulk.gsea.plot.pdf"),p,width=3.6,height=2.4)
}

df<-rename(deseq2_df,BETA=log2FoldChange,SE=lfcSE,Z=stat,LRT_P=pvalue)
write_tsv(df,"DESeq2_results.formatted.for.batch.scripts.tsv.gz")