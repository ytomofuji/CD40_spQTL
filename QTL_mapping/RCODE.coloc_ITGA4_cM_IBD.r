## ---------------------------
##
## Script name: RCODE.coloc_ITGA4_cM_IBD.r
##
## Original file: RCODE.eqtl.pqtl.METAL.coloc.gene.model.ITGA4.cM.IBD.r
##
## Purpose of script: Run colocalization analysis for ITGA4 in monocytes against IBD GWAS
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
pqtl_eqtl_file=argv[1]
pqtl_ibd_file=argv[2]
pqtl_wbc_file=argv[3]
center=as.numeric(argv[4])
window=as.numeric(argv[5])
outname=argv[6]

if(FALSE){
pqtl_eqtl_file="/path/to/ITGA4.cM.pqtl.eqtl.coloc.input.tsv"
pqtl_ibd_file="/path/to/ITGA4.IBD.coloc.input.tsv"
pqtl_wbc_file="/path/to/ITGA4.WBC.coloc.input.tsv"
center=150518448
window=500000
outname="ITGA4.cM"
}

df_eqtl<-read_tsv(pqtl_eqtl_file)
df_ibd<-read_tsv(pqtl_ibd_file)
df_wbc<-read_tsv(pqtl_wbc_file)

max_eqtl<-(-log10(min(df_eqtl$p.eqtl)))*1.1
max_pqtl<-(-log10(min(df_eqtl$p.pqtl)))*1.1
max_ibd<-(-log10(min(df_ibd$p.gwas)))*1.1
max_wbc<-(-log10(min(df_wbc$p.gwas)))*1.1

#merge all data frane
df_ibd_1<-dplyr::select(df_ibd,chr,pos,ref,alt,beta.gwas,p.gwas)
colnames(df_ibd_1)<-c("chr","pos","ref","alt","beta.ibd","p.ibd")

df_wbc_1<-dplyr::select(df_wbc,chr,pos,ref,alt,beta.gwas,p.gwas)
colnames(df_wbc_1)<-c("chr","pos","ref","alt","beta.wbc","p.wbc")

df<-left_join(df_eqtl,df_ibd_1,by=c("chr","pos","ref","alt"))%>%
left_join(df_wbc_1,by=c("chr","pos","ref","alt"))

chromosome<-df$chr[1]
start<-center-window
end<-center+window

p_pqtl<-ggplot(df,aes(x=pos,y=-log10(p.pqtl),colour=LD))+
geom_vline(xintercept = center, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD=="top_variant"),size=4,shape=23,fill="#ea1bd1",color="black")+
theme_classic()+theme(aspect.ratio=1/4)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+ylim(c(0,max_pqtl))+
scale_x_continuous(limits = c(start,end), expand = c(0, 0))

p_eqtl<-ggplot(df,aes(x=pos,y=-log10(p.eqtl),colour=LD))+
geom_vline(xintercept = center, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD=="top_variant"),size=4,shape=23,fill="#ea1bd1",color="black")+
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



p_ibd<-ggplot(df,aes(x=pos,y=-log10(p.ibd),colour=LD))+
geom_vline(xintercept = center, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD=="top_variant"),size=4,shape=23,fill="#ea1bd1",color="black")+
theme_classic()+theme(aspect.ratio=1/4)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+ylim(c(0,max_ibd))+
scale_x_continuous(limits = c(start,end), expand = c(0, 0))

p_wbc<-ggplot(df,aes(x=pos,y=-log10(p.wbc),colour=LD))+
geom_vline(xintercept = center, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD=="top_variant"),size=4,shape=23,fill="#ea1bd1",color="black")+
theme_classic()+theme(aspect.ratio=1/4)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+ylim(c(0,max_wbc))+
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



p<-p_pqtl+p_eqtl+p_ibd+p_wbc+p_gene+plot_layout(ncol=1)
ggsave(str_c(outname,".ld.locus_plot.pdf"),p,width=6,height=9)
