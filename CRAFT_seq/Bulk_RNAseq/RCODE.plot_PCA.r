## ---------------------------
##
## Script name: RCODE.plot_PCA.r
##
## Original file: RCODE.260502.plot.PCA.log.r
##
## Purpose of script: 
##
##
## Date Created: 2026-05-02
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

library(tidyverse)
library(patchwork) # For combining plots
library(ggfortify) # For autoplot(prcomp)
library(broom)
library(ggrastr)
library(ggrepel)
library(viridis)



plot_data<-read_tsv("/path/to/PCA.before.adjustment.tsv.gz")
plot_data$TimePoint<-factor(plot_data$TimePoint, levels=unique(plot_data$TimePoint))

plt_scores <- ggplot(plot_data, aes(x = PC1, y = PC2, color = TimePoint, shape=Genotype)) +
geom_point(size = 1.5, alpha = 0.8) +
scale_color_viridis(discrete=TRUE,option="plasma") +
scale_shape_manual(values=c("C"=16,"T"=17)) +
theme_test() +
theme(
  aspect.ratio = 1/1,
  strip.background = element_blank(),
  strip.text = element_text(face = "bold")
)

ggsave("PCA.PC1.2.for.figure.pdf",plt_scores,width=5,height=3)
