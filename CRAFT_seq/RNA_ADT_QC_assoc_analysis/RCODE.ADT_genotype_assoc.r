## ---------------------------
##
## Script name: RCODE.ADT_genotype_assoc.r
##
## Original file: 260505.RCODE.ADT.only.genotype.assoc.QC.modified.error.resolve.log.r
##
## Purpose of script: Test associations between ADT expression and genotypes,
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
ADT_THR=as.numeric(argv[4])
is_perm=as.logical(argv[5])
Ratio_THR=as.numeric(argv[6])
UMI_THR=as.numeric(argv[7])



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




#also try linear model because ADT is not as sparse as RNA data
lm_test<-function(df,full_form_right,null_form_right,variable,gene){

	full_form<-str_c(gene,"~",full_form_right)
	null_form<-str_c(gene,"~",null_form_right)

	tryCatch({
		full_model <- lm(formula = full_form, data = df)
		null_model <- lm(formula = null_form, data = df)

		model_lrt <- as.data.frame(anova(null_model, full_model))

		out <- summary(full_model)$coefficients
		colnames(out)<-c("BETA","SE","T","P")
		out<-bind_cols(tibble(Gene=gene),as_tibble(t(out[variable,])))%>%mutate(LRT_P=model_lrt[2,"Pr(>F)"])
		return(out) }, error = function(e) {
		return(tibble(Gene = gene, BETA = NA, SE = NA, Z = NA, P = NA, LRT_P = NA))
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
cell_list<-filter(new_meta,!is.na(allele_dosage),nFeature_ADT>ADT_THR)%>%rownames()

data<-subset(data,cells=cell_list)

#make raw count data
#protein-wise analysis (ADT)
edited_cell_list<-filter(new_meta,!is.na(allele_dosage),condition!="NTC",nFeature_ADT>ADT_THR)%>%rownames()

adt_mx <- GetAssayData(data, slot = "counts", assay="ADT")
num_cell<-ncol(adt_mx[,edited_cell_list])
ADT_qc<-tibble(Gene=rownames(adt_mx),UMI=(rowSums(adt_mx[,edited_cell_list])/num_cell),
Ratio=(rowSums(adt_mx[,edited_cell_list]>0)/num_cell))

p<-ggplot(ADT_qc,aes(x=UMI,y=Ratio))+scale_x_log10()+
geom_point_rast(size=0.2)+geom_point(data=subset(ADT_qc,Gene=="CD40"),size=1,color="red")+
geom_point(data=subset(ADT_qc,str_detect(Gene,pattern="-IgG")),size=1,color="blue")+
geom_text_repel(data=subset(ADT_qc,Gene=="CD40"),aes(label=Gene),size=3,color="red")+
geom_text_repel(data=subset(ADT_qc,str_detect(Gene,pattern="-IgG")),aes(label=Gene),size=3,color="blue")+
theme_bw()+theme(aspect.ratio=1/1)+
geom_vline(xintercept=UMI_THR,color="red")+geom_hline(yintercept=Ratio_THR,color="red")

ggsave(str_c(name,".adt.qc.pdf"),p,width=5,height=5)
write_tsv(ADT_qc,str_c(name,".adt.qc.tsv.gz"))

qc_passed_adts<-filter(ADT_qc,UMI>UMI_THR,Ratio>Ratio_THR)%>%pull(Gene)

adt_mx <- GetAssayData(data, slot = "counts", assay="ADT")
adt_clr_mx <- GetAssayData(data, slot = "data", assay="ADT")
adt<-adt_mx[qc_passed_adts,]%>%as.data.frame()%>%t()%>%as.data.frame()%>%rownames_to_column("cell")%>%as_tibble()
adt_clr<-adt_clr_mx[qc_passed_adts,]%>%as.data.frame()%>%t()%>%as.data.frame()%>%rownames_to_column("cell")%>%as_tibble()

colnames(adt)<-gsub("-","_",colnames(adt))
colnames(adt_clr)<-gsub("-","_",colnames(adt_clr))
qc_passed_adts<-gsub("-","_",qc_passed_adts)

df<-inner_join(adt,new_meta,by="cell")%>%filter(condition!="NTC")%>%mutate(logRNA=log10(nCount_RNA),logADT=log10(nCount_ADT))%>%
filter(!is.na(allele_dosage),!is.na(percent.mt))
df_clr<-inner_join(adt_clr,new_meta,by="cell")%>%filter(condition!="NTC")%>%mutate(logRNA=log10(nCount_RNA),logADT=log10(nCount_ADT))%>%
filter(!is.na(allele_dosage),!is.na(percent.mt))


if(is_perm==TRUE){
df[["allele_dosage"]] <- sample(df[["allele_dosage"]])
df_clr[["allele_dosage"]] <- sample(df_clr[["allele_dosage"]])
}




plan(multisession, workers = 8)

ADT.result <- future_map_dfr(
  qc_passed_adts,
  ~ glm_test(
    df = df,
    full_form_right = "allele_dosage + plate + logADT + percent.mt",
    null_form_right = "plate + logADT + percent.mt",
    variable = "allele_dosage",
    gene = .x
  ))



ADT.lm.result <- future_map_dfr(
  qc_passed_adts,
  ~ lm_test(
    df = df_clr,
    full_form_right = "allele_dosage + plate + logADT + percent.mt",
    null_form_right = "plate + logADT + percent.mt",
    variable = "allele_dosage",
    gene = .x
  ))


#save all results into tsv files
ADT.result %>%write_tsv(str_c(name,".ADT.genotype.nb.result.tsv.gz"))
ADT.lm.result %>%write_tsv(str_c(name,".ADT.genotype.lm.result.tsv.gz"))



