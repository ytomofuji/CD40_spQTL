## ---------------------------
##
## Script name: RCODE.pQTL_eQTL_effect_size_comparison.r
##
## Original file: RCODE.pQTL.eQTL.effect.size.comparison.for.simple.fig.v2.log.r
##
## Purpose of script: Compare eQTL and pQTL effect sizes
##
##
## Date Created: 2026-01-09
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
df<-read_tsv("pQTL.eQTL.effect.size.comparison.input.tsv")


#make plots
celltype_colors <- c(
  "T4"="#FC7A62","T8"="#EC008C","NK"="#FF9C42","B"="#5D4AAB",
  "cM"="#34BECC","ncM"="#CCECCC","cDC"="#409188","pDC"="#0077b6"
)

df$plot_cat<-factor(df$plot_cat, levels=c("Safe","Nonsyn"))


p <- ggplot(df, aes(x = BETA_pQTL, y = BETA_eQTL,
                    color = cell_type)) +
  rasterise(geom_errorbar(aes(ymin = BETA_eQTL - SE_eQTL,
                    ymax = BETA_eQTL + SE_eQTL), width = 0), dpi = 300) +
  rasterise(geom_errorbarh(aes(xmin = BETA_pQTL - SE_pQTL,
                     xmax = BETA_pQTL + SE_pQTL), height = 0), dpi = 300) +
  geom_point_rast()+    
  geom_point(data = filter(df,
  (cell_type=="B"&Gene=="CD40")|
  (cell_type=="cM"&Gene=="ITGA4")|  
  (cell_type=="T4"&Gene=="PDCD1")), 
             aes(x = BETA_pQTL, y = BETA_eQTL), 
             color="black",size=2) + 
  geom_errorbar(data = filter(df,
  (cell_type=="B"&Gene=="CD40")|
  (cell_type=="cM"&Gene=="ITGA4")|  
  (cell_type=="T4"&Gene=="PDCD1")),aes(ymin = BETA_eQTL - SE_eQTL*1.96,
                    ymax = BETA_eQTL + SE_eQTL*1.96), width = 0, color="black",linewidth=0.75) +
  geom_errorbarh(data = filter(df,
  (cell_type=="B"&Gene=="CD40")|
  (cell_type=="cM"&Gene=="ITGA4")|  
  (cell_type=="T4"&Gene=="PDCD1")),aes(xmin = BETA_pQTL - SE_pQTL*1.96,
                    xmax = BETA_pQTL + SE_pQTL*1.96), height = 0, color="black",linewidth=0.75) +
  geom_label_repel(data = filter(df,
  (cell_type=="B"&Gene=="CD40")|
  (cell_type=="cM"&Gene=="ITGA4")|  
  (cell_type=="T4"&Gene=="PDCD1")),
                   aes(x = BETA_pQTL, y = BETA_eQTL, label = Gene),
                   fill = "white",       
                   box.padding = 0.5,    
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  scale_color_manual(values = celltype_colors) +
  xlab("pQTL effect sizes") +
  ylab("eQTL effect sizes") + facet_wrap(~plot_cat,ncol=2)

ggsave("eQTL.vs.pQTL.CI.highlight.genes.syn.pdf",p,width=8,height=3)



p <- ggplot(df, aes(x = BETA_pQTL, y = BETA_eQTL,
                    color = cell_type)) +
  geom_point_rast()+    
  geom_point(data = filter(df,
  (cell_type=="B"&Gene=="CD40")|
  (cell_type=="cM"&Gene=="ITGA4")|  
  (cell_type=="T4"&Gene=="PDCD1")), 
             aes(x = BETA_pQTL, y = BETA_eQTL), 
             color="black",size=2) + 
  geom_label_repel(data = filter(df,
  (cell_type=="B"&Gene=="CD40")|
  (cell_type=="cM"&Gene=="ITGA4")|  
  (cell_type=="T4"&Gene=="PDCD1")),
                   aes(x = BETA_pQTL, y = BETA_eQTL, label = Gene),
                   fill = "white",       
                   box.padding = 0.5,    
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  scale_color_manual(values = celltype_colors) +
  xlab("pQTL effect sizes") +
  ylab("eQTL effect sizes") + facet_wrap(~plot_cat,ncol=2)

ggsave("eQTL.vs.pQTL.highlight.genes.syn.pdf",p,width=8,height=3)





p <- ggplot(df, aes(x = BETA_pQTL, y = BETA_eQTL,
                    color = cell_type)) +
  geom_point_rast()+    
  geom_point(data = filter(df,
  (cell_type=="B"&Gene=="CD40")|
  (cell_type=="cM"&Gene=="ITGA4")|  
  (cell_type=="T4"&Gene=="PDCD1")), 
             aes(x = BETA_pQTL, y = BETA_eQTL), 
             color="black",size=2) + 
  geom_label_repel(data = filter(df,
  (cell_type=="B"&Gene=="CD40")|
  (cell_type=="cM"&Gene=="ITGA4")|  
  (cell_type=="T4"&Gene=="PDCD1")),
                   aes(x = BETA_pQTL, y = BETA_eQTL, label = Gene),
                   fill = "white",       
                   box.padding = 0.5,    
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  scale_color_manual(values = celltype_colors) +
  xlab("pQTL effect sizes") +
  ylab("eQTL effect sizes") 

ggsave("eQTL.vs.pQTL.highlight.genes.pdf",p,width=5,height=3)





df<-mutate(df,rep_category=case_when(
	BETA_eQTL*BETA_pQTL>0 & P_eQTL < 0.05/nrow(df) ~ str_c("P < 0.05/",nrow(df)),
	BETA_eQTL*BETA_pQTL>0 & P_eQTL < 0.05 ~ "P < 0.05",
	BETA_eQTL*BETA_pQTL>0  ~ "Same effect direction",
	BETA_eQTL*BETA_pQTL<=0  ~ "Others"	
))

tbl1<-df%>% group_by(rep_category,plot_cat) %>% summarize(Count=n())
tbl1$rep_category<-factor(tbl1$rep_category,
levels=rev(c(str_c("P < 0.05/",nrow(df)),"P < 0.05","Same effect direction","Others")))

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
	scale_fill_manual(values=rev(c("#D93644","#FFB84C","#93C9D9","#E5E5E5"))) + facet_wrap(~plot_cat,ncol=2)

ggsave("pQTL.to.eQTL.pval.category.pie.syn.pdf",p,width=8,height=3)




tbl2<-df%>% group_by(coloc_cat,plot_cat) %>% summarize(Count=n())
tbl2$coloc_cat<-factor(tbl2$coloc_cat,
levels=c("PP4<=0.25","PP4>0.25","PP4>0.50","PP4>0.75"))




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
	scale_fill_manual(values=rev(c("#D93644","#FFB84C","#93C9D9","#E5E5E5"))) + facet_wrap(~plot_cat,ncol=2)

ggsave("pQTL.to.eQTL.coloc.category.pie.syn.pdf",p,width=8,height=3)



