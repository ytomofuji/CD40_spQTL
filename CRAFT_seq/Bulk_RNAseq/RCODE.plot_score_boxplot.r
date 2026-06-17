## ---------------------------
##
## Script name: RCODE.plot_score_boxplot.r
##
## Original file: RCODE.260502.plot.score.boxplot.log.r
##
## Purpose of script: 
##
##
## Date Created: 2024-05-02
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



plot_data<-read_tsv("/path/to/score.original.df.tsv.gz")




genotype_colors <- c("C" = "#00a0e9", "T" = "#e60012")

plot_data<-mutate(plot_data,Value=ifelse(Score_Type=="pos_weighted",-Value,Value))
type_list<-unique(plot_data$Score_Type)
plot_data$TimePoint<-factor(plot_data$TimePoint, levels=unique(plot_data$TimePoint))


for(score_type in type_list){
	tmp<-filter(plot_data,Score_Type==score_type)
	plt_scores <- ggplot(tmp, aes(x = TimePoint, y = Value, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 0.5, alpha = 0.8,  shape=21) +
  facet_wrap(~ Score_Type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = genotype_colors) +
  theme_test() +
  theme(
    aspect.ratio = 2/3,
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(title = score_type,
       x = "Time Point",
       y = "Calculated Score")

ggsave(str_c(score_type,".plot.scores.per.genotype.across.timepoints.pdf"),plt_scores,width=3,height=2)
}