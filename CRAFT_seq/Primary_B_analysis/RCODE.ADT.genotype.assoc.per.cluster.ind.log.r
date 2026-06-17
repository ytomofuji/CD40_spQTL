## ---------------------------
##
## Script name: RCODE.ADT.genotype.assoc.per.cluster.ind.log.r
##
## Original file: 260109.RCODE.ADT.genotype.assoc.per.cluster.ind.log.r
##
## Purpose of script: For ADT, check association to genotype
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




argv=commandArgs(T)
seurat_file=argv[1]
name=argv[2]
is_perm=as.logical(argv[3])
Ratio_THR=as.numeric(argv[4])
UMI_THR=as.numeric(argv[5])
target_clst=argv[6]


if(FALSE){
seurat_file="/path/to/QCed.merge.primaryB.CRAFT.rds"
name="test"
is_perm=TRUE
Ratio_THR=0.5
UMI_THR=1
target_clst="0"
}


library(tidyverse)
library(Seurat)
library(ggpointdensity)
library(viridis)
library(lme4)
library(patchwork)
library(ggrastr)
library(ggrepel)
library(furrr)
library(future)
n_fun<-function(x){
  return(data.frame(y=max(x)*1.1,label=paste0("N = ",length(x))))
}



# as a statistical model, try nb model
glm_test <- function(df, full_form_right, null_form_right, variable, gene) {

  full_form <- as.formula(str_c(gene, " ~ ", full_form_right))
  null_form <- as.formula(str_c(gene, " ~ ", null_form_right))

  # Return NA if either model fails
  tryCatch({
    full_model <- MASS::glm.nb(formula = full_form, data = df)

    # Fit the null model and calculate LRT p-value.
    # If the null model fails, LRT_P will be NA.
    LRT_P_val <- NA # Initialize LRT_P_val to NA
    null_model <- tryCatch(
      MASS::glm.nb(formula = null_form, data = df),
      error = function(e) {
        warning(paste("Null model failed for gene", gene, ":", e$message))
        return(NULL)
      }
    )

    if (!is.null(null_model)) {
          anova_res <- anova(full_model, null_model, test = "Chisq")
          lr_stat <- anova_res$"LR stat."[2]
          df_diff <- anova_res$"   df"[2]
          log_p_natural <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE, log.p = TRUE)
          minus_log10p <- - (log_p_natural / log(10))
          LRT_P_val <- 10**(-minus_log10p)
    }

    out <- summary(full_model)$coefficients
    colnames(out) <- c("BETA", "SE", "Z", "P")

    out <- bind_cols(
      tibble(Gene = gene),
      as_tibble(t(out[variable, ]))
    ) %>%
    mutate(LRT_P = LRT_P_val)

    return(out)

  }, error = function(e) {
    warning(paste("Full model failed for gene", gene, ":", e$message))
    return(tibble(Gene = gene, BETA = NA, SE = NA, Z = NA, P = NA, LRT_P = NA))
  })
}





data<-readRDS(seurat_file)
target_clst<-str_split(target_clst,",")%>%unlist()
target_clst<-as.numeric(target_clst)
#select specific clusters
data<-subset(data,seurat_clusters%in%target_clst)

ADT_mx <- GetAssayData(data, slot = "counts", assay="ADT")
num_cell<-ncol(ADT_mx)
ADT_qc<-tibble(Gene=rownames(ADT_mx),UMI=(rowSums(ADT_mx)/num_cell),
Ratio=(rowSums(ADT_mx>0)/num_cell))

p<-ggplot(ADT_qc,aes(x=UMI,y=Ratio))+scale_x_log10()+
geom_point_rast(size=0.2)+geom_point(data=subset(ADT_qc,Gene=="CD40"),size=1,color="red")+
geom_point(data=subset(ADT_qc,Gene=="CD27"),size=1,color="red")+
geom_text_repel(data=subset(ADT_qc,Gene=="CD40"),aes(label=Gene),size=3,color="red")+
geom_text_repel(data=subset(ADT_qc,Gene=="CD27"),aes(label=Gene),size=3,color="red")+
theme_bw()+theme(aspect.ratio=1/1)+
geom_vline(xintercept=UMI_THR,color="red")+geom_hline(yintercept=Ratio_THR,color="red")

ggsave(str_c(name,".ADT.qc.pdf"),p,width=5,height=5)
write_tsv(ADT_qc,str_c(name,".ADT.qc.tsv.gz"))


qc_passed_ADTs<-filter(ADT_qc,UMI>UMI_THR,Ratio>Ratio_THR)%>%pull(Gene)

ADT_mx <- GetAssayData(data, slot = "counts", assay="ADT")
ADT<-ADT_mx[qc_passed_ADTs,]%>%as.data.frame()%>%t()%>%as.data.frame()%>%rownames_to_column("cell")%>%as_tibble()

colnames(ADT)<-gsub("-","_",colnames(ADT))
qc_passed_ADTs<-gsub("-","_",qc_passed_ADTs)

df<-inner_join(ADT,as_tibble(data@meta.data),by="cell")%>%mutate(logADT=log10(nCount_ADT))%>%mutate(plate_exp=str_c(plate,experiment))


if(is_perm==TRUE){
	for(x in unique(df$plate_exp)){
   df$C_allele_dosage[df$plate_exp == x] <- sample(df$C_allele_dosage[df$plate_exp == x])
}
}



original_df<-df
original_name<-name

for(EXP in c(unique(original_df$experiment),"all")){

if(EXP=="all"){
df<-original_df
name<-original_name
}else{
df<-filter(original_df,experiment==EXP)
name<-str_c(EXP,".",original_name)
}



plan(multisession, workers = 8)

ADT.result <- future_map_dfr(
  qc_passed_ADTs,
  ~ glm_test(
    df = df,
    full_form_right = "C_allele_dosage + plate_exp + logADT + percent.mt",
    null_form_right = "plate_exp + logADT + percent.mt",
    variable = "C_allele_dosage",
    gene = .x
  ))



#save all results into tsv files
ADT.result %>%write_tsv(str_c(name,".ADT.genotype.nb.result.tsv.gz"))


