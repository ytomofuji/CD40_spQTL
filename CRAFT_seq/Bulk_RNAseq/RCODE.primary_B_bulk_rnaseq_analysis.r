## ---------------------------
##
## Script name: RCODE.primary_B_bulk_rnaseq_analysis.r
##
## Original file: RCODE.260502.analysis.bulk.B.cell.RNA.seq.log.r
##
## Purpose of script: Bulk RNA-seq analysis of primary B cells (CC vs TT genotype)
##                    with CRAFT-seq effect size comparison
##
##
## Date Created: 2026-04-30
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
library(patchwork) # For combining plots
library(ggfortify) # For autoplot(prcomp)
library(broom)
library(ggrastr)
library(ggrepel)


data<-read_tsv("/path/to/Bulk.B.cell.C6.to.T9.Daudi.RNA.seq.count.annotation.tsv.gz")
data$Gene<-gsub("_","-",data$Gene)
qc_df<-read_tsv("/path/to/RNA.QC.information.tsv")
sex_df<-read_tsv("/path/to/Sex_inference_result.tsv")%>%dplyr::select(donor,inferred_sex)
#Since C-allele is the effect allele, flip the beta to align the effect allele to other analyses (T-allele)
craft_df<-read_tsv("/path/to/clst0.Primary.B.pQTL.ADT.LRT.P.volcano.FDR.05.tsv.gz")%>%
mutate(BETA=(-BETA))
craft_df$Gene<-gsub("_","-",craft_df$Gene)

qc_df<-qc_df%>%mutate(tmp=`Sample Name`)%>%separate(tmp,into=c("donor","time"),"_")%>%inner_join(sex_df,by="donor")%>%mutate(is_female=ifelse(inferred_sex=="Female",1,0))
qc_df$`Sample Name`<-gsub("C6_","T9_",qc_df$`Sample Name`)

raw_counts_matrix <- data %>%
  dplyr::select(starts_with("C"), starts_with("T")) %>% 
  dplyr::select(-category) %>% 
  as.matrix()


total_counts <- colSums(raw_counts_matrix)

cpm_matrix <- sweep(raw_counts_matrix, 2, total_counts, "/") * 1e6

prop_cpm_above_1 <- rowMeans(cpm_matrix >= 1)

passing_genes <- prop_cpm_above_1 >= 0.50

plt_hist <- tibble(Proportion = prop_cpm_above_1) %>%
  ggplot(aes(x = Proportion)) +
  geom_histogram(binwidth = 0.05, fill = "darkseagreen", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0.50, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Proportion of Samples with CPM >= 1 per Gene",
       x = "Proportion of Samples",
       y = "Gene Count") +
  theme_minimal()

ggsave("Hist.CPM.detection.ratio.pdf", plt_hist, width=3, height=3)

filtered_raw_counts <- raw_counts_matrix[passing_genes, ]
filtered_gene_metadata <- data[passing_genes, ] %>% select(gene_id, Gene, category)

current_total_counts <- colSums(filtered_raw_counts)
ln_cpm_pre <- log(sweep(filtered_raw_counts, 2, current_total_counts, "/") * 1e6 + 1)

sample_info <- data.frame(Sample_ID = colnames(ln_cpm_pre)) %>%
  mutate(
    Genotype = factor(substr(Sample_ID, 1, 1)), # 'C' or 'T'
    TimePoint = factor(sub(".*_", "", Sample_ID), levels = c("0h", "4h", "24h", "72h")),
    log10_TotalCount = log10(current_total_counts)
  ) %>%
  left_join(qc_df, by = c("Sample_ID" = "Sample Name"))


ln_cpm_post <- apply(ln_cpm_pre, 1, function(y) {
  regression_df <- cbind(sample_info, y = y)
  fit <- lm(y ~ TimePoint + Genotype + RIN + log10_TotalCount + is_female, data = regression_df)
  
  coeffs <- coef(fit)
  adjustment <- (coeffs["RIN"] * (regression_df$RIN - mean(regression_df$RIN, na.rm=T))) + 
                (coeffs["log10_TotalCount"] * (regression_df$log10_TotalCount - mean(regression_df$log10_TotalCount, na.rm=T)))+
				(coeffs["is_female"] * (regression_df$is_female - mean(regression_df$is_female, na.rm=T)))
  
  return(y - adjustment)
})
ln_cpm_post <- t(ln_cpm_post)





run_and_plot_pca <- function(data_matrix, meta_data, title_prefix, pc_x = 1, pc_y = 2) {
  pca_input <- t(data_matrix)
  pca_res <- prcomp(pca_input, scale. = TRUE)
  pca_tibble <- as_tibble(pca_res$x, rownames = "Sample_ID") %>%
    left_join(meta_data, by = c("Sample_ID" = "Sample_ID"))
  base_plot <- function(fill_var, var_label, is_discrete = TRUE) {
    p <- autoplot(pca_res, data = meta_data, 
                  x = pc_x, y = pc_y,
                  colour = fill_var, size = 2) +
      labs(title = paste0(title_prefix, ": Colored by ", var_label),
           subtitle = paste0("PC", pc_x, " vs PC", pc_y)) +
      theme_test() + 
      theme(aspect.ratio = 1/1)
    
    if (is_discrete) {
      p <- p + scale_color_viridis_d(option = "plasma")
    } else {
      p <- p + scale_color_viridis_c(option = "magma")
    }
    return(p)
  }
  
  p1 <- base_plot("TimePoint", "Time Point", TRUE)
  p2 <- base_plot("Genotype", "Genotype", TRUE)
  p3 <- base_plot("RIN", "RIN (Quality)", FALSE)
  p4 <- base_plot("log10_TotalCount", "log10(Read Count)", FALSE)
  
  combined_plot <- (p1 + p2) / (p3 + p4)
  return(list(plot = combined_plot, data = pca_tibble))
}

pca_pre_plots <- run_and_plot_pca(ln_cpm_pre, sample_info, "Pre-Correction")
pca_post_plots <- run_and_plot_pca(ln_cpm_post, sample_info, "Post-Correction")

ggsave("PCA.before.adjustment.pdf",pca_pre_plots$plot,width=8,height=5)
ggsave("PCA.after.adjustment.pdf",pca_post_plots$plot,width=8,height=5)

write_tsv(pca_pre_plots$data,"PCA.before.adjustment.tsv.gz")
write_tsv(pca_post_plots$data,"PCA.after.adjustment.tsv.gz")


# Also for PC3 and 4
run_and_plot_pca_pc34 <- function(data_matrix, meta_data, title_prefix) {
  pca_input <- t(data_matrix)
  pca_res <- prcomp(pca_input, scale. = TRUE)
  base_plot <- function(fill_var, var_label, is_discrete = TRUE) {
    p <- autoplot(pca_res, data = meta_data, x = 3, y = 4, colour = fill_var, size = 2) +
      labs(title = paste0(title_prefix, ": Colored by ", var_label)) +
      theme_test()+theme(aspect.ratio=1/1)
    
    if (is_discrete) {
      p <- p + scale_color_viridis_d(option = "plasma")
    } else {
      p <- p + scale_color_viridis_c(option = "magma")
    }
    return(p)
  }
  
  p1 <- base_plot("TimePoint", "Time Point", TRUE)
  p2 <- base_plot("Genotype", "Genotype", TRUE)
  p3 <- base_plot("RIN", "RIN (Quality)", FALSE)
  p4 <- base_plot("log10_TotalCount", "log10(Read Count)", FALSE)
  (p1 + p2) / (p3 + p4)
}

pca_pre_plots_pc34 <- run_and_plot_pca_pc34(ln_cpm_pre, sample_info, "Pre-Correction")
pca_post_plots_pc34 <- run_and_plot_pca_pc34(ln_cpm_post, sample_info, "Post-Correction")

ggsave("PCA_pc34.before.adjustment.pdf",pca_pre_plots_pc34,width=8,height=5)
ggsave("PCA_pc34.after.adjustment.pdf",pca_post_plots_pc34,width=8,height=5)






rownames(ln_cpm_post)<-data$Gene[passing_genes]
rownames(ln_cpm_pre)<-data$Gene[passing_genes]

valid_degs <- filter(craft_df,FDR < 0.05) %>%
  filter(Gene %in% rownames(ln_cpm_post))

pos_genes <- valid_degs %>% filter(category == "pos")
neg_genes <- valid_degs %>% filter(category == "neg")

#num filter passed genes
write_tsv(tibble(Pos_Gene=nrow(pos_genes),Neg_Gene=nrow(neg_genes)),"Number.of.genes.used.for.score.tsv")

calculate_gene_set_score <- function(expr_matrix, target_genes_df, score_name) {
  genes <- target_genes_df$Gene
  sub_matrix <- expr_matrix[genes, , drop = FALSE]
  simple_sum <- colSums(sub_matrix)
  weights <- target_genes_df$BETA
  weighted_matrix <- sweep(sub_matrix, 1, weights, "*")
  weighted_sum <- colSums(weighted_matrix)
  res <- tibble(
    Sample_ID = colnames(expr_matrix),
    !!paste0(score_name, "_simple") := simple_sum,
    !!paste0(score_name, "_weighted") := weighted_sum
  )
  
  return(res)
}

#pos = C-allele-associated-score
#neg = T-allele-associated-score
pos_scores <- calculate_gene_set_score(ln_cpm_post, pos_genes, "pos")
neg_scores <- calculate_gene_set_score(ln_cpm_post, neg_genes, "neg")

final_score_df <- pca_post_plots$data %>%
  left_join(pos_scores, by = "Sample_ID") %>%
  left_join(neg_scores, by = "Sample_ID") %>%
  mutate(combined_weighted = neg_weighted + pos_weighted)


score_vars <- c("pos_simple", "pos_weighted", "neg_simple", "neg_weighted", "combined_weighted")

test_results <- final_score_df %>%
  pivot_longer(cols = all_of(score_vars), names_to = "Score_Type", values_to = "Value") %>%
  group_by(Score_Type, TimePoint) %>%
  summarise(
    t_test_p = t.test(Value ~ Genotype)$p.value,
    wilcox_p = wilcox.test(Value ~ Genotype)$p.value,
    .groups = "drop"
  )

plot_data <- final_score_df %>%
  pivot_longer(cols = all_of(score_vars), names_to = "Score_Type", values_to = "Value") %>%
  mutate(
    Score_Type = factor(Score_Type, levels = score_vars),
    Genotype = factor(Genotype, levels = c("C", "T")) # C on left, T on right
  )

genotype_colors <- c("C" = "#00a0e9", "T" = "#e60012")

plt_scores <- ggplot(plot_data, aes(x = TimePoint, y = Value, fill = Genotype)) +
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
  labs(title = "Gene Set Scores across TimePoints",
       x = "Time Point",
       y = "Calculated Score")

ggsave("Plot.scores.per.genotype.across.timepoints.pdf",plt_scores,width=8,height=4)

write_tsv(test_results,"score.tests.results.tsv")
write_tsv(plot_data,"score.original.df.tsv.gz")



analysis_ready_df <- as_tibble(ln_cpm_pre, rownames = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample_ID", values_to = "Expression") %>%
  left_join(sample_info, by = "Sample_ID") %>%
  mutate(Genotype = factor(Genotype, levels = c("C", "T"))) %>% # Set 'C' as reference
  filter(Gene %in% craft_df$Gene)

write_tsv(analysis_ready_df,"lm.input.df.tsv.gz")

de_results <- analysis_ready_df %>%
  group_by(Gene, TimePoint) %>%
  nest() %>%
  mutate(
    mean_C   = map_dbl(data, ~ mean(.x$Expression[.x$Genotype == "C"], na.rm = TRUE)),
    mean_T   = map_dbl(data, ~ mean(.x$Expression[.x$Genotype == "T"], na.rm = TRUE)),
    median_C = map_dbl(data, ~ median(.x$Expression[.x$Genotype == "C"], na.rm = TRUE)),
    median_T = map_dbl(data, ~ median(.x$Expression[.x$Genotype == "T"], na.rm = TRUE)),
    lm_res = map(data, ~ tidy(lm(Expression ~ Genotype + is_female + RIN + log10_TotalCount, data = .x))),
    ttest_res = map(data, ~ tidy(t.test(Expression ~ Genotype, data = .x))),
    wilcox_res = map(data, ~ tidy(wilcox.test(Expression ~ Genotype, data = .x)))
  )

final_de_table <- de_results %>%
  mutate(
    lm_target = map(lm_res, ~ filter(.x, term == "GenotypeT")),
    lm_beta   = map_dbl(lm_target, ~ .x$estimate),
    lm_se     = map_dbl(lm_target, ~ .x$std.error),
    lm_p      = map_dbl(lm_target, ~ .x$p.value),
    ttest_diff = map_dbl(ttest_res, ~ .x$estimate2 - .x$estimate1),
    ttest_p    = map_dbl(ttest_res, ~ .x$p.value),
    wilcox_p = map_dbl(wilcox_res, ~ .x$p.value)
  ) %>%
  select(Gene, TimePoint, mean_C, mean_T, median_C, median_T, lm_beta, lm_se, lm_p, ttest_diff, ttest_p, wilcox_p) %>%
  ungroup()

final_de_table <- final_de_table %>%
  group_by(TimePoint) %>%
  mutate(
    lm_fdr     = p.adjust(lm_p, method = "fdr"),
    ttest_fdr  = p.adjust(ttest_p, method = "fdr"),
    wilcox_fdr = p.adjust(wilcox_p, method = "fdr")
  ) %>%
  ungroup()


write_tsv(final_de_table,"genotype.deg.analysis.results.across.time.points.tsv.gz")

colnames(craft_df)<-c("Gene", "CRAFT_BETA_T_allele","CRAFT_pvalue","CRAFT_FDR","highlight","category")
all_df<-inner_join(craft_df,final_de_table,by="Gene")

write_tsv(all_df,"genotype.deg.analysis.results.across.time.points.with.CRAFT.seq.results.tsv.gz")



all_df<-mutate(all_df, sign_lm_p=ifelse(lm_beta*CRAFT_BETA_T_allele>0,lm_p/2,1-(lm_p/2)))%>%
mutate(FDR_pass=ifelse(CRAFT_FDR<0.05,"Yes","No"))


#show 2 plots (signed-QQ and effect size scatter plot)
#for effect size scatter, calculate correlation and sign-test
#make volcano plot for lm results

#signed-QQ
#make QQ plot with stratification based on FDR 5%
make_qq<-function(DF,pval.col,FDR_column){
df<-dplyr::select(DF,all_of(c("Gene",pval.col,FDR_column)))
colnames(df)<-c("Gene","meta.P","is_CRAFTseq")
df<-filter(df,!is.na(meta.P))

plt_df<-group_by(df,is_CRAFTseq)%>%
arrange(meta.P)%>%filter(Gene!="CD40")%>%
mutate(RANK=row_number())%>%mutate(EXP=(RANK-0.5)/max(RANK))%>%ungroup()

a<-filter(plt_df,is_CRAFTseq=="Yes")%>%pull(meta.P)
b<-filter(plt_df,is_CRAFTseq=="No")%>%pull(meta.P)
test<-wilcox.test(a,b)
P<-test$p.value

p<-ggplot(plt_df,aes(x=-log10(EXP),
y=ifelse(-log10(meta.P)>10,10,-log10(meta.P)),color=is_CRAFTseq))+
geom_point_rast(alpha=0.8,size=1)+
geom_point(data=subset(plt_df,Gene %in% c("CD40")),alpha=0.8,size=1,color="#00CD00")+
geom_text_repel(data=subset(plt_df,Gene %in% c("CD40")),aes(label=Gene))+
theme_test()+
theme(
axis.text.x = element_text(size = 10),
axis.text.y = element_text(size = 10),
aspect.ratio=1/1)+
geom_hline(yintercept=-log10(0.05/nrow(df)),colour="black",linetype="dashed")+
geom_abline(intercept=0,slope=1,colour="black",linetype="dashed")+
scale_color_manual(values=c("Yes"="#D8031C","No"="#01016F"))+
xlab("Expected(P)")+ylab("Observed(P)")+
ggtitle(str_c("Wilcox.test P = ",signif(P,digits=3),"\nLamnda: Yes ",
signif(median(a),digits=3),", No ",signif(median(b),digits=3)))

return(p)
}

for(timepoint in c("0h","4h","24h","72h")){
p<-make_qq(filter(all_df,TimePoint==timepoint),"sign_lm_p","FDR_pass")
ggsave(str_c(timepoint,".conc.to.CRAFT.seq.signed.QQ.qqplot.pdf"),p,width=5,height=3)
}






#Effect size scatter

for(timepoint in c("0h","4h","24h","72h")){
df<-filter(all_df,TimePoint==timepoint)%>%filter(FDR_pass=="Yes")%>%
mutate(consistent=ifelse(lm_beta*CRAFT_BETA_T_allele>0,"Yes","No"))%>%filter(!is.na(lm_beta),!is.na(CRAFT_BETA_T_allele))%>%
filter(Gene!="CD40")

cor_df<-df
consistent_cases <- sum(sign(cor_df$lm_beta) == sign(cor_df$CRAFT_BETA_T_allele), na.rm = TRUE)
total_cases <- sum(!is.na(cor_df$lm_beta) & !is.na(cor_df$CRAFT_BETA_T_allele))
sign_test_res <- binom.test(
  x = consistent_cases, 
  n = total_cases, 
  p = 0.5, 
  alternative = "greater")

sink(str_c(timepoint,".sign_test_results.txt"))
print(sign_test_res)
sink()

sp_test_res<-cor.test(cor_df$lm_beta,cor_df$CRAFT_BETA_T_allele,method="spearman")
sink(str_c(timepoint,".spearman_test_results.txt"))
print(sp_test_res)
sink()


ylim_val<-c(min(df$lm_beta)*1.1,max(df$lm_beta)*1.1)
xlim_val<-c(min(df$CRAFT_BETA_T_allele)*1.1,max(df$CRAFT_BETA_T_allele)*1.1)

q1 <- sum(df$CRAFT_BETA_T_allele > 0 & df$lm_beta > 0, na.rm = TRUE) # RU
q2 <- sum(df$CRAFT_BETA_T_allele < 0 & df$lm_beta > 0, na.rm = TRUE) # LU
q3 <- sum(df$CRAFT_BETA_T_allele < 0 & df$lm_beta < 0, na.rm = TRUE) # LD
q4 <- sum(df$CRAFT_BETA_T_allele > 0 & df$lm_beta < 0, na.rm = TRUE) # RD


p<-ggplot(df,aes(x=CRAFT_BETA_T_allele,y=lm_beta, color=consistent))+
  scale_color_manual(values=c("Yes"="#D62728","No"="#1F77B4"))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_point(size = 1, shape = 16)+
  theme_test() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "none"
  )+xlab("CRAFT-seq effect")+ylab("QTL effect")+xlim(xlim_val)+ylim(ylim_val)+
  annotate(	"text", x = -Inf, y = Inf, label=str_c("Psign=",signif(sign_test_res$p.value,digits=2),
  "\nconc ratio=",signif(sign_test_res$estimate,digits=2)),
    hjust = -0.1, vjust = 1.1, size = 3)+
    annotate("text", x = xlim_val[2], y = ylim_val[2], label = q1, hjust = 1.2, vjust = 1.5, size = 4, color = "#D62728") + # Q1
  annotate("text", x = xlim_val[1], y = ylim_val[2], label = q2, hjust = -0.2, vjust = 1.5, size = 4, color = "#1F77B4") + # Q2
  annotate("text", x = xlim_val[1], y = ylim_val[1], label = q3, hjust = -0.2, vjust = -0.5, size = 4, color = "#D62728") + # Q3
  annotate("text", x = xlim_val[2], y = ylim_val[1], label = q4, hjust = 1.2, vjust = -0.5, size = 4, color = "#1F77B4")# Q4

ggsave(str_c(timepoint,".craft.bulk.B.CTcomp.coplot.pdf"),p,width=4,height=3)

}





#make volcano plot for lm results
make_volcano<-function(df,p_val_col,gene_list,FDR_THR){
	tmp<-dplyr::select(df,all_of(c("Gene","lm_beta",p_val_col)))
	colnames(tmp)<-c("Gene","lm_beta","pvalue")
	tmp$FDR<-p.adjust(tmp$pvalue, method="BH")
	val1<-filter(tmp,FDR<FDR_THR)%>%pull(pvalue)%>%max()
	val2<-filter(tmp,FDR>=FDR_THR)%>%pull(pvalue)%>%min()
	fdr_thr<-(log10(val1)+log10(val2))*(-0.5)

	tmp<-tmp%>%mutate(highlight=ifelse(Gene %in% gene_list,TRUE,FALSE))%>%
	arrange(pvalue)%>%mutate(category=case_when(
		FDR<FDR_THR&lm_beta>0 ~ "pos",
		FDR<FDR_THR&lm_beta<0 ~ "neg",
		TRUE~"others"
	))

	#also annotate top 10 genes
	top_genes<-head(tmp,10)%>%pull(Gene)

	pos_cnt<-filter(tmp,category=="pos")%>%nrow()
	neg_cnt<-filter(tmp,category=="neg")%>%nrow()

	p<-ggplot(tmp,aes(x=lm_beta,
	y=ifelse(-log10(pvalue)>20,20,-log10(pvalue))))+
	geom_point_rast(data=subset(tmp,highlight==FALSE),alpha=0.8,size=1,aes(color=category))+
	geom_point(data=subset(tmp,highlight==TRUE),alpha=0.8,size=2,color="red")+
	geom_text_repel(data=subset(tmp,Gene %in% top_genes & highlight==FALSE),aes(label=Gene,color=category),alpha=0.8,size=2,color="black")+
	geom_text_repel(data=subset(tmp,highlight==TRUE),aes(label=Gene),alpha=0.8,size=2,color="#FF7F0E")+
	theme_test()+
	theme(
	axis.text.x = element_text(size = 10),
	axis.text.y = element_text(size = 10),
	aspect.ratio=1/1)+
	xlab("Effect size")+ylab("-log10(P)")+
	geom_hline(yintercept=fdr_thr,color="black")+
	geom_vline(xintercept=0,color="black",linetype="dashed")+
	scale_color_manual(values=c("pos"="#D62728","neg"="#1F77B4","others"="#555555"))+
	  geom_text(x = -Inf, y = Inf,
              label = neg_cnt,
              hjust = -0.1, vjust = 1.5, size = 3, color = "#1F77B4")+
    geom_text(x = Inf, y = Inf,
              label = pos_cnt,
              hjust = 1.1, vjust = 1.5, size = 3, color = "#D62728")+ylim(c(0,20.1))

	return(list(plot = p, df = tmp))

}



for(timepoint in c("0h","4h","24h","72h")){
	obj<-make_volcano(filter(all_df,TimePoint==timepoint),"lm_p","CD40",0.05)
	ggsave(str_c(timepoint,".CC.vs.TT.lm.P.volcano.FDR.0.05.pdf"),obj[["plot"]],width=4,height=4)
}


