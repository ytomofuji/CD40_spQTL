## ---------------------------
##
## Script name: RCODE.make.gsea.plot.for.figure.log.r
##
## Original file: RCODE.make.gsea.plot.for.figure.log.r
##
## Purpose of script: Generate GSEA mountain plot for figure
##
##
## Date Created: 2026-02-02
##
## Email:
##
## ---------------------------
##
## Notes: To make concordance, we flipped the allele
##
##
## ---------------------------



argv=commandArgs(T)

sumstats_file=argv[1]

if(FALSE){
	sumstats_file="/path/to/RNA.genotype.nb.result.tsv.gz"
}


#Run pathway enrichment analysis
library(msigdbr)
library(fgsea)
library(tidyverse)
library(ggrastr)
library(ggrepel)
library(org.Hs.eg.db)
library(clusterProfiler)

eqtl_df<-read_tsv(sumstats_file)%>%filter(!is.na(Z))
eqtl_df$Gene<-gsub("_","-",eqtl_df$Gene)
unique_genes<-eqtl_df%>%group_by(Gene)%>%summarize(Count=n())%>%filter(Count==1)%>%pull(Gene)%>%unique()
#all genes were unique


C2_df <- msigdbr(species = "human", category = "C2")
C5_df <- msigdbr(species = "human", category = "C5")

#transform gene symbol to entrez ID
gene_tbl <- bitr(eqtl_df$Gene, fromType = "SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") %>% as_tibble()
gene_tbl<-left_join(gene_tbl,eqtl_df,by=c("SYMBOL"="Gene"))

C2_pathways = split(x = C2_df$entrez_gene, f = C2_df$gs_name)

ranks<-gene_tbl$Z*(-1)
names(ranks)<-gene_tbl$ENTREZID
C2_fgseaRes <- fgsea(pathways=C2_pathways, ranks,
                  minSize  = 15,
                  maxSize  = 500)


CD40_related_pathways<-names(C2_pathways)[str_detect(names(C2_pathways),"CD40")]
as_tibble(C2_fgseaRes)%>%filter(pathway %in% CD40_related_pathways)%>%write_tsv("CRAFT.seq.gsea.C2.CD40.related.results.tsv")
pathway_list<-as_tibble(C2_fgseaRes)%>%filter(pathway %in% CD40_related_pathways)%>%pull(pathway)
as_tibble(C2_fgseaRes)%>%arrange(padj)%>%write_tsv("CRAFT.seq.gsea.C2.results.tsv.gz")


#make volcano plot
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
	geom_line(data = pd$curve, aes(x = rank, y = ES), color = "black", size = 1) +
	geom_segment(data = pd$ticks, aes(x = rank, xend = rank, y = -0.05, yend = 0.05), size = 0.2) +
	geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
	theme_test() + labs( title = target_pathway, x = "Ranked List Metric", y = "Enrichment Score (ES)" ) +
	theme( plot.title = element_text(hjust = 0.5, face = "bold"), axis.title = element_text(size = 12) ) +
	annotate(
	"text", x = Inf, y = Inf,
	label = str_c("Padj=",padj,", NES=",NES),
	hjust = 1.1, vjust = 1.1,
	size = 4, fontface = "italic", color = "red"
	)

	ggsave(str_c(target_pathway,".CRAFT.gsea.plot.pdf"),p,width=3.6,height=2.4)
}
