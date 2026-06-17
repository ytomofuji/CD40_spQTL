## ---------------------------
##
## Script name: RCODE.pQTL_PPP_effect_size_comparison.r
##
## Original file: RCODE.pQTL.PPP.effect.size.comparison.figure.v2.log.r
##
## Purpose of script: Compare pQTL and UKB-PPP effect sizes
##
##
## Date Created: 2026-03-10
##
## Email:
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


library(tidyverse)
library(viridis)
library(ggrastr)
library(ggrepel)

setwd("/path/to/output_directory")

pqtl_file="/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/summary/pQTL.meta.analysis.lead.variants.FDR.0.05.tsv"
PPP_file="/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/260211_UKBPPP_pQTL/reference_data/UKBPPP.sumstats.for.CITE.pQTL.tsv.gz"
annov_file="/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/annovar/annovar.5pop.0.8.nonsyn.pqtl.tag.variants.tsv"
coloc_file="/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/260211_UKBPPP_pQTL/summary/PPP.pQTL.coloc.summary.tsv"

n_fun<-function(x){
  return(data.frame(y=max(x)*1.1,label=paste0("N = ",length(x))))
}

pqtl_df<-read_tsv(pqtl_file)%>%dplyr::select(all_of(c("Gene","SNP","CHR","POS","REF","ALT","cell_type","BETA","SE","P")))
ppp_df<-read_tsv(PPP_file)%>%dplyr::select(all_of(c("Gene","SNP","CHR","POS","REF","ALT","BETA","SE","LOG10P")))
colnames(pqtl_df)<-c("Gene","SNP","CHR","POS","REF","ALT","cell_type","BETA_pQTL","SE_pQTL","P_pQTL")
colnames(ppp_df)<-c("Gene","SNP","CHR","POS","REF","ALT","BETA_PPP","SE_PPP","log10P_PPP")
merge_df<-left_join(pqtl_df,ppp_df,by=c("Gene","SNP","CHR","POS","REF","ALT"))


#FCGR2A_B_C does not map to single-gene in reference file. Remove it for this analysis
df<-filter(merge_df,!is.na(BETA_PPP) & Gene != "FCGR2A_B_C")

#also add colocalization / non-syn information
#for non-syn information, prioritize: "startloss" -> "stopgain" -> "frameshift substitution" -> "nonsynonymous SNV"
nonsyn_df<-read_tsv(annov_file)%>%mutate(merge_key_gene=ifelse(Gene.refGene %in% c("FCGR2A","FCGR2B","FCGR2C"),"FCGR2A_B_C",Gene.refGene))
nonsyn_anno<-filter(nonsyn_df,X7==merge_key_gene)%>%dplyr::select(X8,merge_key_gene,Gene.refGene,ExonicFunc.refGene)%>%dplyr::distinct()
df<-left_join(df,nonsyn_anno,by=c("Gene"="merge_key_gene","cell_type"="X8"))

coloc_df<-read_tsv(coloc_file)

df <- left_join(df,coloc_df,by=c("Gene"="protein","cell_type"="cell"))
df$ExonicFunc.refGene[is.na(df$ExonicFunc.refGene)]<-"None"
df<-mutate(df,plot_cat=ifelse(ExonicFunc.refGene=="None","Safe","Nonsyn"))
df<-mutate(df,coloc_cat=case_when(
	PP.H4.abf>0.75~"PP4>0.75",
	PP.H4.abf>0.50~"PP4>0.50",
	PP.H4.abf>0.25~"PP4>0.25",
	PP.H4.abf<=0.25~"PP4<=0.25"
	))%>%
	mutate(coloc_nopQTL_cat=case_when(
	PP.H0.abf+PP.H1.abf>0.75~"no_pQTL>0.75",
	PP.H0.abf+PP.H1.abf>0.50~"no_pQTL>0.50",
	PP.H0.abf+PP.H1.abf>0.25~"no_pQTL>0.25",
	PP.H0.abf+PP.H1.abf<=0.25~"no_pQTL<=0.25"
	))%>%
	mutate(coloc_nocoloc_cat=case_when(
	PP.H2.abf+PP.H3.abf>0.75~"nocoloc>0.75",
	PP.H2.abf+PP.H3.abf>0.50~"nocoloc>0.50",
	PP.H2.abf+PP.H3.abf>0.25~"nocoloc>0.25",
	PP.H2.abf+PP.H3.abf<=0.25~"nocoloc<=0.25"
	))



#make plots
celltype_colors <- c(
  "T4"="#FC7A62","T8"="#EC008C","NK"="#FF9C42","B"="#5D4AAB",
  "cM"="#34BECC","ncM"="#CCECCC","cDC"="#409188","pDC"="#0077b6"
)

exonic_colors <- c(
  "None" = "#bdbdbd",
  "frameshift substitution" = "#e41a1c",
  "nonsynonymous SNV" = "#377eb8"
)



#effect size comparison
p <- ggplot(df, aes(x = BETA_pQTL, y = BETA_PPP,
                    color = cell_type, fill = cell_type)) +
  geom_errorbar(aes(ymin = BETA_PPP,
                    ymax = BETA_PPP), width = 0) +
  geom_errorbarh(aes(xmin = BETA_pQTL,
                     xmax = BETA_pQTL), height = 0) +
  geom_point(size = 1, shape = 21) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(values = celltype_colors) +
  scale_fill_manual(values = celltype_colors) +
  xlab("pQTL effect sizes") +
  ylab("PPP effect sizes") 

ggsave("PPP.vs.pQTL.simple.pdf",p,width=3,height=3)







df<-mutate(df,rep_category=case_when(
	BETA_PPP*BETA_pQTL>0 & 10**(-log10P_PPP) < 0.05/nrow(df) ~ str_c("P < 0.05/",nrow(df)),
	BETA_PPP*BETA_pQTL>0 & 10**(-log10P_PPP) < 0.05 ~ "P < 0.05",
	BETA_PPP*BETA_pQTL>0  ~ "Same effect direction",
	BETA_PPP*BETA_pQTL<=0  ~ "Others"	
))

tbl1<-df%>% group_by(rep_category) %>% summarize(Count=n())
tbl1$rep_category<-factor(tbl1$rep_category,
levels=rev(c(str_c("P < 0.05/",nrow(df)),"P < 0.05","Same effect direction","Others")))


plot_data <- tbl1 %>% 
    mutate(Total = sum(Count)) %>%
    mutate(
        Percentage = Count / Total,
        LabelText = paste0(Count, " (", scales::percent(Percentage, accuracy = 0.1), ")")
    ) %>%
    arrange(desc(rep_category)) %>% 
    mutate(
        y_end = cumsum(Percentage),
        y_start = lag(y_end, default = 0),
        LabelPosition = (y_start + y_end) / 2
    ) %>%
    arrange(rep_category)

p <- ggplot(plot_data, aes(x = 2, y = Percentage, fill = rep_category)) +
    geom_bar(
        stat = "identity",
        color = "white", 
        linewidth = 0.8
    ) +
    geom_text(
        aes(y = LabelPosition, label = LabelText),
        color = "black", 
        size = 3.5
    ) +
    coord_polar(theta = "y", start = 0) +
    labs(
        title = "Distribution of Categories",
        fill = "Category",
        x = NULL, y = NULL 
    ) +
    theme_void() + 
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
	scale_fill_manual(values=rev(c("#D93644","#FFB84C","#93C9D9","#E5E5E5")))

ggsave("pQTL.to.PPP.pval.simple.pie.pdf",p,width=6,height=3)




tbl2<-df%>% group_by(coloc_cat) %>% summarize(Count=n())
tbl2$coloc_cat<-factor(tbl2$coloc_cat,
levels=c("PP4<=0.25","PP4>0.25","PP4>0.50","PP4>0.75"))



plot_data2 <- tbl2 %>% 
    mutate(Total = sum(Count)) %>%
    mutate(
        Percentage = Count / Total,
        LabelText = paste0(Count, " (", scales::percent(Percentage, accuracy = 0.1), ")")
    ) %>%
    arrange(desc(coloc_cat)) %>% 
    mutate(
        y_end = cumsum(Percentage),
        y_start = lag(y_end, default = 0),
        LabelPosition = (y_start + y_end) / 2
    ) %>%
    arrange(coloc_cat)

p <- ggplot(plot_data2, aes(x = 2, y = Percentage, fill = coloc_cat)) +
    geom_bar(
        stat = "identity",
        color = "white", 
        linewidth = 0.8
    ) +
    geom_text(
        aes(y = LabelPosition, label = LabelText),
        color = "black", 
        size = 3.5
    ) +
    coord_polar(theta = "y", start = 0) +
    labs(
        title = "Distribution of Categories",
        fill = "Category",
        x = NULL, y = NULL 
    ) +
    theme_void() + 
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
	scale_fill_manual(values=rev(c("#D93644","#FFB84C","#93C9D9","#E5E5E5")))
ggsave("pQTL.to.PPP.coloc.simple.pie.pdf",p,width=6,height=3)








#effect size comparison with exonic outcome
df$plot_cat<-factor(df$plot_cat,levels=c("Safe","Nonsyn"))


p <- ggplot(df, aes(x = BETA_pQTL, y = BETA_PPP,
                    color = cell_type, fill = cell_type)) +
  geom_errorbar(aes(ymin = BETA_PPP,
                    ymax = BETA_PPP), width = 0) +
  geom_errorbarh(aes(xmin = BETA_pQTL,
                     xmax = BETA_pQTL), height = 0) +
  geom_point(size = 1, shape = 21) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(values = celltype_colors) +
  scale_fill_manual(values = celltype_colors) +
  xlab("pQTL effect sizes") +
  ylab("PPP effect sizes") + facet_wrap(~plot_cat,ncol=2)

ggsave("PPP.vs.pQTL.simple.exonic.pdf",p,width=8,height=3)







df<-mutate(df,rep_category=case_when(
	BETA_PPP*BETA_pQTL>0 & 10**(-log10P_PPP) < 0.05/nrow(df) ~ str_c("P < 0.05/",nrow(df)),
	BETA_PPP*BETA_pQTL>0 & 10**(-log10P_PPP) < 0.05 ~ "P < 0.05",
	BETA_PPP*BETA_pQTL>0  ~ "Same effect direction",
	BETA_PPP*BETA_pQTL<=0  ~ "Others"	
))

tbl1<-df%>% group_by(rep_category,plot_cat) %>% summarize(Count=n())
tbl1$rep_category<-factor(tbl1$rep_category,
levels=rev(c(str_c("P < 0.05/",nrow(df)),"P < 0.05","Same effect direction","Others")))
tbl1$plot_cat<-factor(tbl1$plot_cat,levels=c("Safe","Nonsyn"))


plot_data <- tbl1 %>% group_by(plot_cat) %>%
    mutate(Total = sum(Count)) %>%
    mutate(
        Percentage = Count / Total,
        LabelText = paste0(Count, " (", scales::percent(Percentage, accuracy = 0.1), ")")
    ) %>%
    arrange(desc(rep_category)) %>% 
    mutate(
        y_end = cumsum(Percentage),
        y_start = lag(y_end, default = 0),
        LabelPosition = (y_start + y_end) / 2
    ) %>%
    arrange(rep_category)

p <- ggplot(plot_data, aes(x = 2, y = Percentage, fill = rep_category)) +
    geom_bar(
        stat = "identity",
        color = "white", 
        linewidth = 0.8
    ) +
    geom_text(
        aes(y = LabelPosition, label = LabelText),
        color = "black", 
        size = 3.5
    ) +
    coord_polar(theta = "y", start = 0) +
    labs(
        title = "Distribution of Categories",
        fill = "Category",
        x = NULL, y = NULL 
    ) +
    theme_void() + 
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
	scale_fill_manual(values=rev(c("#D93644","#FFB84C","#93C9D9","#E5E5E5")))+facet_wrap(~plot_cat,ncol=2)

ggsave("pQTL.to.PPP.pval.simple.pie.exonic.pdf",p,width=8,height=3)




tbl2<-df%>% group_by(coloc_cat,plot_cat) %>% summarize(Count=n())
tbl2$coloc_cat<-factor(tbl2$coloc_cat,
levels=c("PP4<=0.25","PP4>0.25","PP4>0.50","PP4>0.75"))
tbl2$plot_cat<-factor(tbl2$plot_cat,levels=c("Safe","Nonsyn"))


plot_data2 <- tbl2 %>% group_by(plot_cat) %>%
    mutate(Total = sum(Count)) %>%
    mutate(
        Percentage = Count / Total,
        LabelText = paste0(Count, " (", scales::percent(Percentage, accuracy = 0.1), ")")
    ) %>%
    arrange(desc(coloc_cat)) %>% 
    mutate(
        y_end = cumsum(Percentage),
        y_start = lag(y_end, default = 0),
        LabelPosition = (y_start + y_end) / 2
    ) %>%
    arrange(coloc_cat)

p <- ggplot(plot_data2, aes(x = 2, y = Percentage, fill = coloc_cat)) +
    geom_bar(
        stat = "identity",
        color = "white", 
        linewidth = 0.8
    ) +
    geom_text(
        aes(y = LabelPosition, label = LabelText),
        color = "black", 
        size = 3.5
    ) +
    coord_polar(theta = "y", start = 0) +
    labs(
        title = "Distribution of Categories",
        fill = "Category",
        x = NULL, y = NULL 
    ) +
    theme_void() + 
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
	scale_fill_manual(values=rev(c("#D93644","#FFB84C","#93C9D9","#E5E5E5")))+facet_wrap(~plot_cat,ncol=2)

ggsave("pQTL.to.PPP.coloc.simple.pie.exonic.pdf",p,width=8,height=3)



