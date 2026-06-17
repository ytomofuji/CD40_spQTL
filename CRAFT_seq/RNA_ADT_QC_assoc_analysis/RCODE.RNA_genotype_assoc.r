## ---------------------------
##
## Script name: RCODE.RNA_genotype_assoc.r
##
## Original file: 260505.RCODE.RNA.only.genotype.assoc.QC.modified.error.resolve.log.r
##
## Purpose of script: Test associations between RNA expression and genotypes,
##                    with robust handling of non-convergent null models
##
##
## Date Created: 2026-05-05
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
genotype_file=argv[2]
name=argv[3]
RNA_THR=as.numeric(argv[4])
MT_THR=as.numeric(argv[5])
is_perm=as.logical(argv[6])
Ratio_THR=as.numeric(argv[7])
UMI_THR=as.numeric(argv[8])




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



# as a statistical model, try nb model with fallback for init.theta
glm_test <- function(df, full_form_right, null_form_right, variable, gene) {

  full_form <- as.formula(str_c(gene, " ~ ", full_form_right))
  null_form <- as.formula(str_c(gene, " ~ ", null_form_right))
  
  # Helper function to fit glm.nb with a fallback to init.theta = 1
  fit_nb_safe <- function(form, data) {
    tryCatch({
      # 1st attempt: Default initialization
      MASS::glm.nb(formula = form, data = data)
    }, error = function(e) {
      # 2nd attempt: Retry with manual init.theta
      tryCatch({
        MASS::glm.nb(formula = form, data = data, init.theta = 1)
      }, error = function(e2) {
        # Both attempts failed
        return(NULL)
      })
    })
  }

  full_model <- fit_nb_safe(full_form, df)
  
  # If the full model fails even after retry, return NA tibble
  if (is.null(full_model)) {
    warning(paste("Full model failed for gene", gene, "after fallback"))
    return(tibble(Gene = gene, BETA = NA, SE = NA, Z = NA, P = NA, LRT_P = NA))
  }

  null_model <- fit_nb_safe(null_form, df)
  
  LRT_P_val <- NA
  
  if (!is.null(null_model)) {
    tryCatch({
      anova_res <- anova(full_model, null_model)
      # MASS::glm.nb's anova uses "LR stat." and "   df"
      lr_stat <- anova_res$"LR stat."[2]
      df_diff <- anova_res$"   df"[2]
      
      log_p_natural <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE, log.p = TRUE)
      minus_log10p <- - (log_p_natural / log(10))
      LRT_P_val <- 10**(-minus_log10p)
    }, error = function(e) {
      warning(paste("LRT failed for gene", gene, ":", e$message))
    })
  } else {
    warning(paste("Null model failed for gene", gene, "after fallback"))
  }

  tryCatch({
    out_mat <- summary(full_model)$coefficients
    colnames(out_mat) <- c("BETA", "SE", "Z", "P")
    
    # Ensure the variable exists in the coefficients
    if (!variable %in% rownames(out_mat)) {
      return(tibble(Gene = gene, BETA = NA, SE = NA, Z = NA, P = NA, LRT_P = LRT_P_val))
    }

    out <- bind_cols(
      tibble(Gene = gene),
      as_tibble(t(out_mat[variable, ]))
    ) %>%
    mutate(LRT_P = LRT_P_val)
    
    return(out)
    
  }, error = function(e) {
    warning(paste("Coefficient extraction failed for gene", gene, ":", e$message))
    return(tibble(Gene = gene, BETA = NA, SE = NA, Z = NA, P = NA, LRT_P = LRT_P_val))
  })
}




data<-readRDS(seurat_file)
geno<-read_tsv(genotype_file)


#merge genotype info into seurat object
new_meta<-data@meta.data
new_meta<-new_meta%>%mutate(cell_name=str_c(plate,"_",name))%>%left_join(geno,by="cell_name")%>%
mutate(allele_cat=case_when(allele_dosage==0~"Ref/Ref",allele_dosage==1~"Ref/Alt",allele_dosage==2~"Alt/Alt"))
new_meta$allele_cat<-factor(new_meta$allele_cat, levels=c("Ref/Ref","Ref/Alt","Alt/Alt"))
rownames(new_meta)<-new_meta$cell
data@meta.data<-new_meta
cell_list<-filter(new_meta,!is.na(allele_dosage),nFeature_RNA>RNA_THR,percent.mt<MT_THR)%>%rownames()

data<-subset(data,cells=cell_list)

#make raw count data
#protein-wise analysis (RNA)
edited_cell_list<-filter(new_meta,!is.na(allele_dosage),condition!="NTC",nFeature_RNA>RNA_THR,percent.mt<MT_THR)%>%rownames()

RNA_mx <- GetAssayData(data, slot = "counts", assay="RNA")
num_cell<-ncol(RNA_mx[,edited_cell_list])
RNA_qc<-tibble(Gene=rownames(RNA_mx),UMI=(rowSums(RNA_mx[,edited_cell_list])/num_cell),
Ratio=(rowSums(RNA_mx[,edited_cell_list]>0)/num_cell))

p<-ggplot(RNA_qc,aes(x=UMI,y=Ratio))+scale_x_log10()+
geom_point_rast(size=0.2)+geom_point(data=subset(RNA_qc,Gene=="CD40"),size=1,color="red")+
geom_point(data=subset(RNA_qc,Gene=="CD27"),size=1,color="red")+
geom_text_repel(data=subset(RNA_qc,Gene=="CD40"),aes(label=Gene),size=3,color="red")+
geom_text_repel(data=subset(RNA_qc,Gene=="CD27"),aes(label=Gene),size=3,color="red")+
theme_bw()+theme(aspect.ratio=1/1)+
geom_vline(xintercept=UMI_THR,color="red")+geom_hline(yintercept=Ratio_THR,color="red")

ggsave(str_c(name,".RNA.qc.pdf"),p,width=5,height=5)
write_tsv(RNA_qc,str_c(name,".RNA.qc.tsv.gz"))


qc_passed_RNAs<-filter(RNA_qc,UMI>UMI_THR,Ratio>Ratio_THR)%>%pull(Gene)

RNA_mx <- GetAssayData(data, slot = "counts", assay="RNA")
RNA<-RNA_mx[qc_passed_RNAs,]%>%as.data.frame()%>%t()%>%as.data.frame()%>%rownames_to_column("cell")%>%as_tibble()

colnames(RNA)<-gsub("-","_",colnames(RNA))
qc_passed_RNAs<-gsub("-","_",qc_passed_RNAs)

df<-inner_join(RNA,new_meta,by="cell")%>%filter(condition!="NTC")%>%mutate(logRNA=log10(nCount_RNA))%>%
filter(!is.na(allele_dosage),!is.na(percent.mt))


if(is_perm==TRUE){
df[["allele_dosage"]] <- sample(df[["allele_dosage"]])
}




plan(multisession, workers = 8)

RNA.result <- future_map_dfr(
  qc_passed_RNAs,
  ~ glm_test(
    df = df,
    full_form_right = "allele_dosage + plate + logRNA + percent.mt",
    null_form_right = "plate + logRNA + percent.mt",
    variable = "allele_dosage",
    gene = .x
  ))



#save all results into tsv files
RNA.result %>%write_tsv(str_c(name,".RNA.genotype.nb.result.tsv.gz"))


