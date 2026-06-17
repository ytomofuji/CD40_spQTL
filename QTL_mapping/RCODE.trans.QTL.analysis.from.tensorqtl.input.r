## ---------------------------
##
## Script name: RCODE.trans.QTL.analysis.from.tensorqtl.input.r
##
## Original file: RCODE.trans.QTL.analysis.from.tensorqtl.input.r
##
## Purpose of script: Run trans-QTL analysis using tensorqtl input
##
##
## Date Created: 2025-09-13
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




argv=commandArgs(T)
exp_file=argv[1]
cov_file=argv[2]
plink_file=argv[3]
target_variant=argv[4]
name=argv[5]



PLINK2 <- "/path/to/software/genome/plink2"

exp_df<-read_table(exp_file)
exp_mx<-dplyr::select(exp_df,-c(1,2,3,4))%>%t()
colnames(exp_mx)<-str_c("feature",seq(1,length(exp_df$gene_id)))
feature_tbl<-tibble(Feature=exp_df$gene_id, Gene=str_c("feature",seq(1,length(exp_df$gene_id))))

exp_mx<-rownames_to_column(as.data.frame(exp_mx),"ID")%>%as_tibble()
cov_df<-read.table(cov_file)%>%t()%>%as.data.frame()%>%rownames_to_column("ID")%>%as_tibble()
colnames(cov_df)<-gsub("\\-","_",colnames(cov_df))
colnames(cov_df)<-gsub("\\.","_",colnames(cov_df))

exp_trait_list<-exp_df$gene_id
cov_list<-colnames(cov_df)[2:ncol(cov_df)]
#call plink2 function for preparing genotype file
variant_info<-str_split(target_variant,":")%>%unlist()
tibble(name=target_variant, allele=variant_info[4])%>%write_tsv(str_c(name,".variant.export.info.txt"), col_names=F)

command <- paste(
  PLINK2,
  "--export 'A'",
  "--threads 1",
  "--memory 4000",
  "--bfile", plink_file,
  "--export-allele", str_c(name,".variant.export.info.txt"),
  "--extract", str_c(name,".variant.export.info.txt"),
  "--out ", name,
  "--keep-allele-order"
)

system(command)
file.remove(str_c(name,".variant.export.info.txt"))

geno_df<-read_table(str_c(name,".raw"))%>%dplyr::select(2,7)
colnames(geno_df)<-c("ID","allele_dosage")
#remove all characters which may cause problem in merge
exp_mx$ID<-gsub("[X\\.\\-\\_\\-]","",exp_mx$ID)
cov_df$ID<-gsub("[X\\.\\-\\_\\-]","",cov_df$ID)
geno_df$ID<-gsub("[X\\.\\-\\_\\-]","",geno_df$ID)

input_df<-inner_join(cov_df, exp_mx, by="ID")%>%inner_join(geno_df,by="ID")


#run linear regression across traits
# (exp_trait_list) ~ allele_dosage + cov_list
lm_test<-function(df,full_form_right,null_form_right,variable,gene){

	full_form<-str_c(gene,"~",full_form_right)
	null_form<-str_c(gene,"~",null_form_right)

	tryCatch({
		full_model <- lm(formula = full_form, data = df)
		null_model <- lm(formula = null_form, data = df)

		model_lrt <- as.data.frame(anova(null_model, full_model))

		out <- summary(full_model)$coefficients
		colnames(out)<-c("BETA","SE","T","P")
		out<-bind_cols(tibble(Gene=gene),as_tibble(t(out[variable,])))%>%
		mutate(LRT_P=model_lrt[2,"Pr(>F)"],N=length(summary(full_model)$residuals))
		return(out) }, error = function(e) {
		return(tibble(Gene = gene, BETA = NA, SE = NA, Z = NA, P = NA, LRT_P = NA, N = NA))
	})
}

exp_trait_list<-str_c("feature",seq(1,length(exp_df$gene_id)))
cov_list<-colnames(cov_df)[2:ncol(cov_df)]

lm.result <- map_dfr(
  exp_trait_list,
  ~ lm_test(
    df = input_df,
    full_form_right = str_c("allele_dosage+",str_c(cov_list,collapse="+")),
    null_form_right = str_c(cov_list,collapse="+"),
    variable = "allele_dosage",
    gene = .x
  ))

lm.result<-inner_join(feature_tbl,lm.result, by="Gene")%>%dplyr::select(-Gene)%>%rename("Gene"="Feature")

write_tsv(lm.result, str_c(name,".trans.analysis.tsv.gz"))