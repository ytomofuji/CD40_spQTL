## ---------------------------
##
## Script name: RCODE.check_effect_consistency.r
##
## Original file: RCODE.260314.check.consistency.of.effects.log.r
##
## Purpose of script: Check consistency of pQTL/eQTL effect sizes across datasets
##
##
## Date Created: 2026-03-14
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
library(patchwork)
library(ggrepel)

df<-read_tsv("../pQTL.meta.analysis.lead.variants.FDR.0.05.tsv")


celltype_colors <- c(
  "T4"="#FC7A62","T8"="#EC008C","NK"="#FF9C42","B"="#5D4AAB",
  "cM"="#34BECC","ncM"="#CCECCC","cDC"="#409188","pDC"="#0077b6"
)


# define pairs
pairs <- list(
  c("Perez.ASI", "Perez.EUR"),
  c("Perez.ASI", "AMP.RA"),
  c("Perez.EUR", "AMP.RA"),
  c("Perez.ASI", "AMP.SLE"),
  c("Perez.EUR", "AMP.SLE"),
  c("AMP.RA", "AMP.SLE")
)


# make pairwise plot
make_plot <- function(x, y) {
  x_beta <- paste0(x, "_BETA")
  y_beta <- paste0(y, "_BETA")
  x_se <- paste0(x, "_SE")
  y_se <- paste0(y, "_SE")

  valid_idx <- complete.cases(df[[x_beta]], df[[y_beta]])
  n_valid <- sum(valid_idx)
  
  ggplot(df, aes(x = .data[[x_beta]], y = .data[[y_beta]], color = cell_type)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    theme_test(base_size = 10) +
    theme(aspect.ratio = 1,
          legend.position = "none",
          panel.grid = element_blank()) +
    scale_color_manual(values = celltype_colors) +
    xlab(paste0(x, " BETA")) +
    ylab(paste0(y, " BETA")) +
    ggtitle(sprintf("%s vs %s  (%d signals)", x, y, n_valid))
}

plots <- lapply(pairs, function(p) make_plot(p[1], p[2]))

p <- (
  plots[[1]] + plot_spacer() + plot_spacer() +
  plots[[2]] + plots[[3]] + plot_spacer() +
  plots[[4]] + plots[[5]] + plots[[6]]
) + plot_layout(ncol = 3, nrow = 3, guides = 'collect')

ggsave("effect.size.comparisons.across.datasets.pdf",p,width=9,height=9)




# Pearson correlation for each plot
# make pairwise plot
annotate_correlation <- function(cor_result, xpos = -Inf, ypos = Inf, hjust = -0.1, vjust = 1.2) {
  annotate("text",
           x = xpos, y = ypos,
           label = sprintf("Spearman cor. = %.2f\np = %.1e", cor_result$estimate, cor_result$p.value),
           hjust = hjust, vjust = vjust, size = 3.5)
}

make_plot_wcor <- function(x, y) {
  x_beta <- paste0(x, "_BETA")
  y_beta <- paste0(y, "_BETA")
  x_se <- paste0(x, "_SE")
  y_se <- paste0(y, "_SE")

  valid_idx <- complete.cases(df[[x_beta]], df[[y_beta]])
  n_valid <- sum(valid_idx)
  cor_result <- suppressWarnings(cor.test(df[[x_beta]], df[[y_beta]], method = "pearson", use = "complete.obs"))
  
  ggplot(df, aes(x = .data[[x_beta]], y = .data[[y_beta]], color = cell_type)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
	annotate_correlation(cor_result) +
    theme_test(base_size = 10) +
    theme(aspect.ratio = 1,
          legend.position = "none",
          panel.grid = element_blank()) +
    scale_color_manual(values = celltype_colors) +
    xlab(paste0(x, " BETA")) +
    ylab(paste0(y, " BETA")) +
    ggtitle(sprintf("%s vs %s  (%d signals)", x, y, n_valid))
}



plots <- lapply(pairs, function(p) make_plot_wcor(p[1], p[2]))

p <- (
  plots[[1]] + plot_spacer() + plot_spacer() +
  plots[[2]] + plots[[3]] + plot_spacer() +
  plots[[4]] + plots[[5]] + plots[[6]]
) + plot_layout(ncol = 3, nrow = 3, guides = 'collect')

ggsave("effect.size.comparisons.across.datasets.w.pearson.pdf",p,width=9,height=9)



#also make venn diagram for trait overlap
library(ggVennDiagram)



#plot dataset
plt_df<-df%>%mutate(name=str_c(Gene,cell_type,SNP))

Perez_ASI<-filter(plt_df,!is.na(Perez.ASI_BETA))%>%pull(name)
Perez_EUR<-filter(plt_df,!is.na(Perez.EUR_BETA))%>%pull(name)
AMP_RA<-filter(plt_df,!is.na(AMP.RA_BETA))%>%pull(name)
AMP_SLE<-filter(plt_df,!is.na(AMP.SLE_BETA))%>%pull(name)


venn_list <- list(
  Perez_ASI = Perez_ASI,
  Perez_EUR = Perez_EUR,
  AMP_RA = AMP_RA,
  AMP_SLE = AMP_SLE
)

p <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F6E6FA", high = "#9B59B6") +
  theme_void() +
  theme(
    text = element_text(family = "Helvetica", size = 10),
    legend.position = "none"
  )

ggsave("Sig.hit.overlap.venn.pdf",p)



#also make plot for the number of original traits
all_df<-read_tsv("../pQTL.meta.analysis.full.result.tsv.gz")

Perez_ASI<-filter(all_df,!is.na(Perez.ASI_BETA))%>%pull(Gene)%>%unique()
Perez_EUR<-filter(all_df,!is.na(Perez.EUR_BETA))%>%pull(Gene)%>%unique()
AMP_RA<-filter(all_df,!is.na(AMP.RA_BETA))%>%pull(Gene)%>%unique()
AMP_SLE<-filter(all_df,!is.na(AMP.SLE_BETA))%>%pull(Gene)%>%unique()


venn_list <- list(
  Perez_ASI = Perez_ASI,
  Perez_EUR = Perez_EUR,
  AMP_RA = AMP_RA,
  AMP_SLE = AMP_SLE
)

p <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F6E6FA", high = "#9B59B6") +
  theme_void() +
  theme(
    text = element_text(family = "Helvetica", size = 10),
    legend.position = "none"
  )



ggsave("Trait.overlap.venn.pdf",p)