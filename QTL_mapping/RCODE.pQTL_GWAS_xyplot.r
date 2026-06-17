## ---------------------------
##
## Script name: RCODE.pQTL_GWAS_xyplot.r
##
## Original file: RCODE.pQTL.GWAS.xyplot.r
##
## Purpose of script: Make xy plot for colocalization analysis
##
##
## Date Created: 2026-05-19
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

B_file=argv[1]
cM_file=argv[2]
outname=argv[3]
pos1=as.numeric(argv[4])
pos2=as.numeric(argv[5])
window=as.numeric(argv[6])



df_B<-read_tsv(B_file)
df_cM<-read_tsv(cM_file)

df<-inner_join(df_B,df_cM,by=c("chr","pos","ref","alt"))

if(!is.na(window)){
center<-(pos1+pos2)/2
start<-center-window
end<-center+window
df<-filter(df,pos>=start&pos<=end)
}

pqtl_max<-max(-log10(c(df$p.qtl.x,df$p.qtl.y)))*1.1
pqtl_cM_max<-max(-log10(c(df$p.qtl.y)))*1.1
gwas_max<-max(-log10(c(df$p.gwas.x,df$p.gwas.y)))*1.1

#.x = for B, .y = for cM
#B cell
B_plt<-ggplot(df,aes(x=(log10(p.qtl.x))*sign(beta.qtl.x),y=(log10(p.gwas.x))*sign(beta.gwas.x),colour=LD.x))+
geom_hline(yintercept=0,linetype="dashed")+
geom_vline(xintercept=0,linetype="dashed")+
geom_point_rast(data=subset(df,LD.x=="[0,0.2]"),size=1.5)+
geom_point_rast(data=subset(df,LD.x!="[0,0.2]"),size=1.5)+
geom_point(data=subset(df,pos==pos1),size=3,shape=23,fill="#D24655",color="black")+
geom_point(data=subset(df,pos==pos2),size=3,shape=23,fill="#0F8EAE",color="black")+
theme_test()+theme(aspect.ratio=1/1)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+xlim(c(-pqtl_max,pqtl_max))+ylim(c(-gwas_max,gwas_max))


#cM cell
cM_plt<-ggplot(df,aes(x=(log10(p.qtl.y))*sign(beta.qtl.y),y=(log10(p.gwas.y))*sign(beta.gwas.y),colour=LD.x))+
geom_hline(yintercept=0,linetype="dashed")+
geom_vline(xintercept=0,linetype="dashed")+
geom_point_rast(data=subset(df,LD.x=="[0,0.2]"),size=1.5)+
geom_point_rast(data=subset(df,LD.x!="[0,0.2]"),size=1.5)+
geom_point(data=subset(df,pos==pos1),size=3,shape=23,fill="#D24655",color="black")+
geom_point(data=subset(df,pos==pos2),size=3,shape=23,fill="#0F8EAE",color="black")+
theme_test()+theme(aspect.ratio=1/1)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+xlim(c(-pqtl_max,pqtl_max))+ylim(c(-gwas_max,gwas_max))



#cM cell with adjusted scale
cM_plt_low<-ggplot(df,aes(x=(log10(p.qtl.y))*sign(beta.qtl.y),y=(log10(p.gwas.y))*sign(beta.gwas.y),colour=LD.x))+
geom_hline(yintercept=0,linetype="dashed")+
geom_vline(xintercept=0,linetype="dashed")+
geom_point_rast(data=subset(df,LD.x=="[0,0.2]"),size=1.5)+
geom_point_rast(data=subset(df,LD.x!="[0,0.2]"),size=1.5)+
geom_point(data=subset(df,pos==pos1),size=3,shape=23,fill="#D24655",color="black")+
geom_point(data=subset(df,pos==pos2),size=3,shape=23,fill="#0F8EAE",color="black")+
theme_test()+theme(aspect.ratio=1/1)+
scale_color_manual(values=c(
	"top_variant"="#ea1bd1",
	"(0.8,1]"="#da3026",
	"(0.6,0.8]"="#f39a28",
	"(0.4,0.6]"="#56c345",
	"(0.2,0.4]"="#8abce0",
	"[0,0.2]"="#2f1977"
))+xlim(c(-pqtl_cM_max,pqtl_cM_max))+ylim(c(-gwas_max,gwas_max))




p<-B_plt+cM_plt+cM_plt_low+plot_layout(ncol=3)
ggsave(str_c(outname,".xyplot.pdf"),p,width=12,height=3)


filter(df,LD.x!="[0,0.2]"&LD.y!="[0,0.2]")%>%write_tsv(str_c(outname,".problematic.SNPs.tsv"))






#also make 2 LD color ploy
df2<-mutate(df,LD.x=str_c(LD.x,"_1"),LD.y=str_c(LD.y,"_2"))
B_plt<-ggplot(df2,aes(x=(log10(p.qtl.x))*sign(beta.qtl.x),y=(log10(p.gwas.x))*sign(beta.gwas.x),colour=LD.x))+
geom_hline(yintercept=0,linetype="dashed")+
geom_vline(xintercept=0,linetype="dashed")+
geom_point_rast(data=subset(df2,LD.x=="[0,0.2]_1"&LD.y=="[0,0.2]_2"),size=2,color="grey")+
geom_point_rast(data=subset(df2,LD.x!="[0,0.2]_1"),aes(color=LD.x),size=2)+
geom_point_rast(data=subset(df2,LD.y!="[0,0.2]_2"),aes(color=LD.y),size=2)+
geom_point(data=subset(df2,pos==pos1),size=4,shape=23,fill="#D24655",color="black")+
geom_point(data=subset(df2,pos==pos2),size=4,shape=23,fill="#0F8EAE",color="black")+
scale_color_manual(values=c(
  "(0.8,1]_1"   = "#D24655", 
  "(0.6,0.8]_1" = "#E17B85",
  "(0.4,0.6]_1" = "#EEAFB5",
  "(0.2,0.4]_1" = "#F8D7DA", 
  "(0.8,1]_2"   = "#0F8EAE", 
  "(0.6,0.8]_2" = "#57ACC3",
  "(0.4,0.6]_2" = "#91C9D8",
  "(0.2,0.4]_2" = "#C7E6ED"
))+
theme_test()+theme(aspect.ratio=1/1)+
xlim(c(-pqtl_max,pqtl_max))+ylim(c(-gwas_max,gwas_max))




#cM cell
cM_plt<-ggplot(df2,aes(x=(log10(p.qtl.y))*sign(beta.qtl.y),y=(log10(p.gwas.y))*sign(beta.gwas.y),colour=LD.x))+
geom_hline(yintercept=0,linetype="dashed")+
geom_vline(xintercept=0,linetype="dashed")+
geom_point_rast(data=subset(df2,LD.x=="[0,0.2]_1"&LD.y=="[0,0.2]_2"),size=2,color="grey")+
geom_point_rast(data=subset(df2,LD.x!="[0,0.2]_1"),aes(color=LD.x),size=2)+
geom_point_rast(data=subset(df2,LD.y!="[0,0.2]_2"),aes(color=LD.y),size=2)+
geom_point(data=subset(df2,pos==pos1),size=4,shape=23,fill="#D24655",color="black")+
geom_point(data=subset(df2,pos==pos2),size=4,shape=23,fill="#0F8EAE",color="black")+
scale_color_manual(values=c(
  "(0.8,1]_1"   = "#D24655", 
  "(0.6,0.8]_1" = "#E17B85",
  "(0.4,0.6]_1" = "#EEAFB5",
  "(0.2,0.4]_1" = "#F8D7DA", 
  "(0.8,1]_2"   = "#0F8EAE", 
  "(0.6,0.8]_2" = "#57ACC3",
  "(0.4,0.6]_2" = "#91C9D8",
  "(0.2,0.4]_2" = "#C7E6ED"
))+
theme_test()+theme(aspect.ratio=1/1)+
xlim(c(-pqtl_max,pqtl_max))+ylim(c(-gwas_max,gwas_max))



#cM cell with adjusted scale
cM_plt_low<-ggplot(df2,aes(x=(log10(p.qtl.y))*sign(beta.qtl.y),y=(log10(p.gwas.y))*sign(beta.gwas.y),colour=LD.x))+
geom_hline(yintercept=0,linetype="dashed")+
geom_vline(xintercept=0,linetype="dashed")+
geom_point_rast(data=subset(df2,LD.x=="[0,0.2]_1"&LD.y=="[0,0.2]_2"),size=2,color="grey")+
geom_point_rast(data=subset(df2,LD.x!="[0,0.2]_1"),aes(color=LD.x),size=2)+
geom_point_rast(data=subset(df2,LD.y!="[0,0.2]_2"),aes(color=LD.y),size=2)+
geom_point(data=subset(df2,pos==pos1),size=4,shape=23,fill="#D24655",color="black")+
geom_point(data=subset(df2,pos==pos2),size=4,shape=23,fill="#0F8EAE",color="black")+
scale_color_manual(values=c(
  "(0.8,1]_1"   = "#D24655", 
  "(0.6,0.8]_1" = "#E17B85",
  "(0.4,0.6]_1" = "#EEAFB5",
  "(0.2,0.4]_1" = "#F8D7DA", 
  "(0.8,1]_2"   = "#0F8EAE", 
  "(0.6,0.8]_2" = "#57ACC3",
  "(0.4,0.6]_2" = "#91C9D8",
  "(0.2,0.4]_2" = "#C7E6ED"
))+
theme_test()+theme(aspect.ratio=1/1)+
xlim(c(-pqtl_cM_max,pqtl_cM_max))+ylim(c(-gwas_max,gwas_max))


p<-B_plt+cM_plt+cM_plt_low+plot_layout(ncol=3)
ggsave(str_c(outname,".xyplot.bicolor.pdf"),p,width=12,height=3)