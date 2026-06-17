## ---------------------------
##
## Script name: RCODE.coloc_KLRG1_T8.r
##
## Original file: RCODE.eqtl.pqtl.METAL.coloc.gene.model.KLRG1.T8.r
##
## Purpose of script: Run colocalization analysis for KLRG1 in CD8+ T cells
##
##
## Date Created: 2026-01-09
## 
## Email:
##
## ---------------------------
##
## Notes: Adjust y-axis-scale for figure
##   
##
## ---------------------------



suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggrepel))


argv=commandArgs(T)

input_file=argv[1]
gencode_file=argv[2]
center=as.numeric(argv[3])
window=as.numeric(argv[4])
outname=argv[5]


df<-read_tsv(input_file)

max_eqtl<-(-log10(min(df$p.eqtl)))*1.1
max_pqtl<-(-log10(min(df$p.pqtl)))*1.1

chromosome<-df$chr[1]
start<-center-window
end<-center+window

p1<-ggplot(df,aes(x=pos,y=-log10(p.pqtl),colour=LD))+
geom_vline(xintercept = center, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point(data=subset(df,LD=="top_variant"),size=4,shape=23,fill="#ea1bd1",color="black")+
theme_classic()+theme(aspect.ratio=1/4)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+ylim(c(0,10))+
scale_x_continuous(limits = c(start,end), expand = c(0, 0))

p2<-ggplot(df,aes(x=pos,y=-log10(p.eqtl),colour=LD))+
geom_vline(xintercept = center, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point(data=subset(df,LD=="top_variant"),size=4,shape=23,fill="#ea1bd1",color="black")+
theme_classic()+theme(aspect.ratio=1/4)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+ylim(c(0,max_eqtl))+
scale_x_continuous(limits = c(start,end), expand = c(0, 0))







#load gencode data
gencode_df<-read_tsv(gencode_file,comment="#",col_names=F)

extract_gene<-function(x){
	y<-str_split(x,"; ")%>%unlist()
	gene<-y[str_detect(y,pattern="gene_name")]
	z<-str_split(gene,pattern="\\\"")%>%unlist()
	return(z[2])
}

region<-filter(gencode_df,X1==chromosome&((X4>start&X4<end) | (X5>start&X5<end)))
region$Gene<-lapply(region$X9,extract_gene)%>%unlist()
gene_list<-unique(region$Gene)

plt<-filter(region,str_detect(X9,pattern="protein_coding") | str_detect(X9,pattern="lncRNA"))%>%
mutate(X4=ifelse(X4<start,start,X4),
X5=ifelse(X5>end,end,X5))


p_gene <- ggplot(plt) +
geom_vline(xintercept = center, color = "grey", linetype="dashed") +
geom_segment(data = subset(plt, X3 == "gene"),
				aes(x = X4, xend = X5, y = Gene, yend = Gene, colour = X7),
				linewidth = 0.5) +
geom_segment(data = subset(plt, X3 == "exon"),
				aes(x = X4, xend = X5, y = Gene, yend = Gene, colour = X7),
				linewidth = 1.6) +
geom_text(data=subset(plt,X3=="gene"),aes(label=Gene,x=(X4+X5)/2,y=Gene),size=2,nudge_y = 0.5)+
scale_y_discrete(
	breaks = NULL,
	labels = NULL,
	name = NULL
) +
scale_color_manual(values = c("+" = "#D8031C", "-" = "#01016F")) +
scale_x_continuous(limits = c(start,end), expand = c(0, 0))+
theme_test() +
theme(
	aspect.ratio = 1 / 6,
	panel.grid.major.y = element_blank(), 
	panel.grid.minor.y = element_blank()
)



p<-p2+p1+p_gene+plot_layout(ncol=1)
ggsave(str_c(outname,".ld.locus_plot.pdf"),p,width=6,height=5)
