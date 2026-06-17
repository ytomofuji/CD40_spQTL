## ---------------------------
##
## Script name: RCODE.coloc_PDCD1_T48_Asthma.r
##
## Original file: RCODE.260126.eqtl.pqtl.METAL.coloc.gene.model.PDCD1.T48.Asthma.all.in.one.r
##
## Purpose of script: Run colocalization analysis for PDCD1 in T4/T8 cells against Asthma GWAS
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
pqtl_eqtl_file2=argv[2]
pqtl_gwas_file=argv[3]
pos1=as.numeric(argv[4])
pos2=as.numeric(argv[5])
window=as.numeric(argv[6])
outname=argv[7]

if(FALSE){
pqtl_eqtl_file="/path/to/PDCD1.T4.pqtl.eqtl.coloc.input.tsv"
pqtl_eqtl_file2="/path/to/PDCD1.T8.pqtl.eqtl.coloc.input.tsv"
pqtl_gwas_file="/path/to/PDCD1.Asthma.coloc.input.tsv"
pos1=241771790
pos2=241785788
window=500000
outname="PDCD1.T48"
}

df_eqtl<-read_tsv(pqtl_eqtl_file)
df_eqtl2<-read_tsv(pqtl_eqtl_file2)%>%dplyr::select(-LD)
colnames(df_eqtl2)[5:(ncol(df_eqtl2)-1)]<-str_c("T8.",colnames(df_eqtl2)[5:(ncol(df_eqtl2)-1)])
df_gwas<-read_tsv(pqtl_gwas_file)

max_eqtl<-(-log10(min(df_eqtl$p.eqtl)))*1.1
max_pqtl<-(-log10(min(df_eqtl$p.pqtl)))*1.1
max_pqtl2<-(-log10(min(df_eqtl2$T8.p.pqtl)))*1.1
max_gwas<-(-log10(min(df_gwas$p.gwas)))*1.1

#merge all data frane
df<-left_join(df_eqtl,df_gwas,by=c("chr","pos","ref","alt","LD"))%>%
left_join(df_eqtl2,by=c("chr","pos","ref","alt"))

chromosome<-df$chr[1]
center<-(pos1+pos2)/2
start<-center-window
end<-center+window

p_pqtl<-ggplot(df,aes(x=pos,y=-log10(p.pqtl),colour=LD))+
geom_vline(xintercept = pos1, color = "grey", linetype="dashed") +
geom_vline(xintercept = pos2, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,pos==pos1),size=4,shape=23,fill="#D24655",color="black")+
geom_point_rast(data=subset(df,pos==pos2),size=4,shape=23,fill="#0F8EAE",color="black")+
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

p_pqtl2<-ggplot(df,aes(x=pos,y=-log10(T8.p.pqtl),colour=LD))+
geom_vline(xintercept = pos1, color = "grey", linetype="dashed") +
geom_vline(xintercept = pos2, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,pos==pos1),size=4,shape=23,fill="#D24655",color="black")+
geom_point_rast(data=subset(df,pos==pos2),size=4,shape=23,fill="#0F8EAE",color="black")+
theme_classic()+theme(aspect.ratio=1/4)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+ylim(c(0,max_pqtl2))+
scale_x_continuous(limits = c(start,end), expand = c(0, 0))

p_eqtl<-ggplot(df,aes(x=pos,y=-log10(p.eqtl),colour=LD))+
geom_vline(xintercept = pos1, color = "grey", linetype="dashed") +
geom_vline(xintercept = pos2, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,pos==pos1),size=4,shape=23,fill="#D24655",color="black")+
geom_point_rast(data=subset(df,pos==pos2),size=4,shape=23,fill="#0F8EAE",color="black")+
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



p_gwas<-ggplot(df,aes(x=pos,y=-log10(p.gwas),colour=LD))+
geom_vline(xintercept = pos1, color = "grey", linetype="dashed") +
geom_vline(xintercept = pos2, color = "grey", linetype="dashed") +
geom_point_rast(data=subset(df,LD=="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,LD!="[0,0.2]"),size=2)+
geom_point_rast(data=subset(df,pos==pos1),size=4,shape=23,fill="#D24655",color="black")+
geom_point_rast(data=subset(df,pos==pos2),size=4,shape=23,fill="#0F8EAE",color="black")+
theme_classic()+theme(aspect.ratio=1/4)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+ylim(c(0,max_gwas))+
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
geom_vline(xintercept = pos1, color = "grey", linetype="dashed") +
geom_vline(xintercept = pos2, color = "grey", linetype="dashed") +
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



p<-p_pqtl+p_eqtl+p_gwas+p_pqtl2+p_gene+plot_layout(ncol=1)
ggsave(str_c(outname,".ld.locus_plot.for.fig.all.in.one.pdf"),p,width=6,height=9)





