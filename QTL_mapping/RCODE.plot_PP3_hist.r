## ---------------------------
##
## Script name: RCODE.plot_PP3_hist.r
##
## Original file: RCODE.260218.plot.PP3.hist.for.figure.log.r
##
## Purpose of script: Plot max PP3 for each protein and conditional analysis results
##
##
## Date Created: 2026-02-16
##
## Email:
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


setwd("/path/to/output_directory")

library(tidyverse)
library(ggrastr)

ori_pqtl="/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/summary/pQTL.meta.analysis.lead.variants.FDR.0.05.tsv"
coloc_file="/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/260213_celltype_pQTL/LD/coloc.across.celltypes.tsv"
coloc_df<-read_tsv(coloc_file)

plt_df <- coloc_df%>%group_by(Gene)%>%
arrange(-PP.H3.abf)%>%slice(1)%>%ungroup()

n_high_pp <- sum(plt_df$PP.H3.abf > 0.5, na.rm = TRUE)
n_middle_pp <- sum(plt_df$PP.H3.abf <= 0.5 & plt_df$PP.H3.abf > 0.25, na.rm = TRUE)
n_low_pp <- sum(plt_df$PP.H3.abf <= 0.25, na.rm = TRUE)


p<-ggplot(plt_df, aes(x=PP.H3.abf))+
geom_histogram(binwidth=0.1,
boundary = 0, closed = "left",fill="#FF4B4B",color="black")+
theme_test()+theme(aspect.ratio=1/1)+
geom_vline(xintercept=0.5,linetype="dashed")+
geom_vline(xintercept=0.25,linetype="dashed")+
annotate("text", 
           x = Inf, y = Inf, 
           label = paste0(n_high_pp, "proteins"), 
           hjust = 1.1, 
           vjust = 1.5, 
           size = 5, 
           fontface = "bold")+
annotate("text", 
           x = 0.5, y = Inf, 
           label = paste0(n_middle_pp, "proteins"), 
           hjust = 1.1, 
           vjust = 1.5, 
           size = 5, 
           fontface = "bold")+
annotate("text", 
           x = 0, y = Inf, 
           label = paste0(n_low_pp, "proteins"), 
           hjust = 0, 
           vjust = 1.5,  
           size = 5, 
           fontface = "bold")

ggsave("max.PP.H3.abf.per.gene.hist.3part.pdf",p,width=4,height=4)

write_tsv(plt_df,"max.PP.H3.abf.per.gene.ori.tsv")

plt_df<-mutate(plt_df, PP3_group=case_when(
	PP.H3.abf>0.5~"PP.H3>0.50",
	PP.H3.abf>0.25~"PP.H3>0.25",
	PP.H3.abf<=0.25~"PP.H3<=0.25"	
))



p<-ggplot(plt_df, aes(x=PP.H3.abf,y=PP.H4.abf,color=PP3_group))+
geom_point()+
theme_test()+theme(aspect.ratio=1/1)+
geom_vline(xintercept=0.5,linetype="dashed")+
geom_vline(xintercept=0.25,linetype="dashed")+
geom_abline(intercept=1,slope=-1,linetype="dashed")+
scale_color_manual(values=c("PP.H3>0.50"="#D8031C","PP.H3>0.25"="#BC80BD","PP.H3<=0.25"="#01016F"))

ggsave("PP.H3.vs.PP.H4.abf.per.gene.point.pdf",p,width=4,height=4)




#also check conditional analysis results
cond_file="/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/260213_celltype_pQTL/conditional_analysis/summary/summary.for.conditional.analysis.tsv"
cond_df<-read_tsv(cond_file)
cond_df<-dplyr::select(cond_df,-c(PP.H3.abf,PP3))
plot_df<-left_join(cond_df,plt_df,by="Gene")

ori_data<-read_tsv(ori_pqtl)
fdr_thr<-max(ori_data$P)


summary_wide <- plot_df %>%
  select(Name, cond, P, PP3_group) %>%
  pivot_wider(names_from = cond, values_from = P)
  


p <- ggplot(summary_wide, aes(x = ifelse(-log10(before)>20,20,-log10(before)), y = ifelse(-log10(after)>20,20,-log10(after)), color = PP3_group)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(fdr_thr), linetype = "dotted", color = "grey30") +
  geom_vline(xintercept = -log10(fdr_thr), linetype = "dotted", color = "grey30") +
  geom_point_rast(alpha = 0.7, size = 2) +
  scale_color_manual(values=c("PP.H3>0.50"="#D8031C","PP.H3>0.25"="#BC80BD","PP.H3<=0.25"="#01016F")) +
  theme_test() +
  theme(aspect.ratio = 1/1) +
  labs(
    x = "-log10(P) before conditioning",
    y = "-log10(P) after conditioning",
    title = "Conditioning effect on spQTL signals"
  )

ggsave("conditioning_scatter_plot.3group.pdf", p, width = 5, height = 5)


plot_df$cond<-factor(plot_df$cond,levels=c("before","after"))

p<-ggplot(plot_df,aes(x=cond,y=(-log10(P)),group=Name,color=PP3_group))+
geom_point_rast()+geom_line()+geom_hline(yintercept=(-log10(fdr_thr)))+
scale_color_manual(values=c("PP.H3>0.50"="#D8031C","PP.H3>0.25"="#BC80BD","PP.H3<=0.25"="#01016F")) +
theme_test()+theme(aspect.ratio=2/1)+facet_wrap(~PP3_group,ncol=3)

ggsave("pval.comp.before.after.conditioning.3group.strat.pdf",p,width=6,height=3)



plot_df<-mutate(plot_df,sig=ifelse(P<fdr_thr,"Yes","No"))

plot_df%>%
group_by(cond,PP3_group,sig)%>%
summarize(Count=n())%>%
write_tsv("Num.sig.hits.after.conditioning.3group.strat.tsv")

plot_df%>%
write_tsv("summary.for.conditional.analysis.3group.strat.tsv")