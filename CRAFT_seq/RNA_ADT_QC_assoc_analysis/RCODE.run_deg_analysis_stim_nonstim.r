## ---------------------------
##
## Script name: RCODE.run_deg_analysis_stim_nonstim.r
##
## Original file: RCODE.260205.run.deg.analysis.stim.nonstim.r
##
## Purpose of script: Analyze differential expression/abundance between stimulated
##                    and non-stimulated cells by allele dosage
##
##
## Date Created: 2026-02-05
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
stim_status=argv[2]
cond_status=argv[3]
downsample=as.logical(argv[4])
modality=argv[5]
UMI_THR=as.numeric(argv[6])
Ratio_THR=as.numeric(argv[7])
name=argv[8]




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

set.seed(9)

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

meta<-data@meta.data

#select condition
if(cond_status=="pQTL"){
cells<-filter(meta,condition==cond_status)%>%
mutate(allele_dosage=pQTL_allele_dosage)%>%filter(!is.na(allele_dosage))
}else if (cond_status=="KO") {
cells<-filter(meta,condition==cond_status)%>%
mutate(allele_dosage=KO_allele_dosage)%>%filter(!is.na(allele_dosage))
}


if(downsample){
	num_cells<-filter(cells,stimulation!=stim_status)%>%
	group_by(allele_dosage)%>%summarize(Count=n())

	df<-tibble()

	for(i in 1:nrow(num_cells)){
		dosage<-num_cells$allele_dosage[i]
		cnt<-num_cells$Count[i]
		tmp<-filter(cells,stimulation==stim_status)%>% 
  		slice_sample(n = cnt)

		df<-bind_rows(df,tmp)
	}
}else{
	df<-filter(cells,stimulation==stim_status)
}


#extract expression and apply QC
mx <- GetAssayData(subset(data,cells=df$cell), slot = "counts", assay=modality)
num_cell<-ncol(mx)
qc<-tibble(Gene=rownames(mx),UMI=(rowSums(mx)/num_cell),
Ratio=(rowSums(mx>0)/num_cell))


p<-ggplot(qc,aes(x=UMI,y=Ratio))+scale_x_log10()+
geom_point_rast(size=0.2)+geom_point(data=subset(qc,Gene=="CD40"),size=1,color="red")+
geom_point(data=subset(qc,Gene=="CD27"),size=1,color="red")+
geom_text_repel(data=subset(qc,Gene=="CD40"),aes(label=Gene),size=3,color="red")+
geom_text_repel(data=subset(qc,Gene=="CD27"),aes(label=Gene),size=3,color="red")+
theme_bw()+theme(aspect.ratio=1/1)+
geom_vline(xintercept=UMI_THR,color="red")+geom_hline(yintercept=Ratio_THR,color="red")

ggsave(str_c(name,".qc.pdf"),p,width=5,height=5)
write_tsv(qc,str_c(name,".qc.tsv.gz"))


qc_passed<-filter(qc,UMI>UMI_THR,Ratio>Ratio_THR)%>%pull(Gene)
name_conv<-tibble(ori_name=qc_passed, conv_name=gsub("-","_",qc_passed))
rownames(mx)<-gsub("-","_",rownames(mx))

test_df<-mx%>%as.data.frame()%>%t()%>%as.data.frame()%>%rownames_to_column("cell")%>%as_tibble()%>%
inner_join(df,by="cell")%>%as_tibble()%>%rename(plate=plate.x)%>%mutate(logRNA=log10(nCount_RNA),logADT=log10(nCount_ADT))



plan(multisession, workers = 4)

if(modality=="RNA"){
RNA.result <- future_map_dfr(
  name_conv$conv_name,
  ~ glm_test(
    df = test_df,
    full_form_right = "allele_dosage + plate + logRNA + percent.mt",
    null_form_right = "plate + logRNA + percent.mt",
    variable = "allele_dosage",
    gene = .x
  )
)

RNA.result<-left_join(name_conv,RNA.result,by=c("conv_name"="Gene"))%>%dplyr::select(-conv_name)%>%rename(Gene=ori_name)
RNA.result %>%write_tsv(str_c(name,".genotype.nb.result.tsv.gz"))
}


if(modality=="ADT"){
ADT.result <- future_map_dfr(
  name_conv$conv_name,
  ~ glm_test(
    df = test_df,
    full_form_right = "allele_dosage + plate + logADT + percent.mt",
    null_form_right = "plate + logADT + percent.mt",
    variable = "allele_dosage",
    gene = .x
  )
)

ADT.result<-left_join(name_conv,ADT.result,by=c("conv_name"="Gene"))%>%dplyr::select(-conv_name)%>%rename(Gene=ori_name)

ADT.result %>%write_tsv(str_c(name,".genotype.nb.result.tsv.gz"))
}


write_tsv(test_df,str_c(name,".genotype.nb.ori.df.tsv.gz"))

test_df%>%
group_by(allele_dosage)%>%summarize(Count=n())%>%
write_tsv(str_c(name,".num.cells.tsv"))