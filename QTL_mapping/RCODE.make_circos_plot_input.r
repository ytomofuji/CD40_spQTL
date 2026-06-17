## ---------------------------
##
## Script name: RCODE.make_circos_plot_input.r
##
## Original file: RCODE.251204.make.circos.plot.input.log.r
##
## Purpose of script: Generate circos/fuji plot input for summarizing pQTL results
##
##
## Date Created: 2025-12-04
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
library(scales)
# the input format looks like following,
#We may also want to make a column showing significance (Bonferroni or FDR)
#LOCUS_ID	CATEGORY	TRAIT	CHR	BP	MARKER	GENE
#1	Blood pressure	PP	1	10796866	rs880315	CASZ1
#1	Blood pressure	SBP	1	10796866	rs880315	CASZ1
#1	Blood pressure	MAP	1	10796866	rs880315	CASZ1
#1	Blood pressure	DBP	1	10796866	rs880315	CASZ1
#2	Kidney-related	sCr	1	16316119	rs848302	ZBTB17,C1orf64
#2	Kidney-related	eGFR	1	16316119	rs848302	ZBTB17,C1orf64
#3	Kidney-related	eGFR	1	23692229	rs4525087	ZNF436
#3	Kidney-related	sCr	1	23702384	rs9887774	ZNF436-AS1,TCEA3
#4	Protein	AG	1	25286151	rs4648888	RUNX3
#4	Protein	NAP	1	25286151	rs4648888	RUNX3


fdr_file<-"/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/summary/pQTL.meta.analysis.lead.variants.FDR.0.05.tsv"
bf_file<-"/path/to/analysis/B_eQTL/merge_datasets/250908_Perez_AMP_RA_SLE/250910_METAL/summary/pQTL.meta.analysis.lead.variants.Bonferroni.0.05.tsv"

fdr_data<-read_tsv(fdr_file)
bf_data<-read_tsv(bf_file)

fdr_plt<-dplyr::select(fdr_data,cell_type,Gene,SNP,CHR,POS)
bf_plt<-dplyr::select(bf_data,cell_type,Gene,SNP,CHR,POS)%>%mutate(Category="BF")

plt<-left_join(fdr_plt,bf_plt,by=c("cell_type","Gene","SNP","CHR","POS"))
plt$cell_type<-factor(plt$cell_type, levels=c("T4","T8","NK","B","cM","ncM","cDC","pDC"))
plt$Category[is.na(plt$Category)]<-"FDR"
plt<-arrange(plt,CHR,POS)

#Gene = locus ID
gene_list<-unique(plt$Gene)
gene_tbl<-tibble(Gene=gene_list,LOCUS_ID=seq(1,length(gene_list)))
plt<-left_join(plt,gene_tbl,by="Gene")%>%mutate(cell_type2=cell_type)

out<-dplyr::select(plt,LOCUS_ID,cell_type,cell_type2,CHR,POS,SNP,Gene,Category)
colnames(out)<-c("LOCUS_ID","CATEGORY","TRAIT","CHR","BP","MARKER","GENE","Significance")
write_tsv(out,"pQTL.fuji.input.txt")


#also output traitlist
trait_df<-tibble(
	CATEGORY=c("T4","T8","NK","B","cM","ncM","cDC","pDC"),
	TRAIT=c("T4","T8","NK","B","cM","ncM","cDC","pDC"),
	COLOR=c("#FC7A62","#EC008C","#FF9C42","#5D4AAB",
	"#34BECC","#CCECCC","#409188","#0077b6")
)

write_tsv(trait_df,"pQTL.fuji.traitlist.txt")


#make barplot for summarizeing the number of the association
num_assoc<-out %>% group_by(CATEGORY,Significance) %>% summarize(Count=n())
num_assoc$Significance<-factor(num_assoc$Significance, levels=rev(c("BF","FDR")))
num_assoc$CATEGORY<-factor(num_assoc$CATEGORY, levels=rev(c("T4","T8","NK","B","cM","ncM","cDC","pDC")))

COLOR_BF <- c("#FC7A62","#EC008C","#FF9C42","#5D4AAB","#34BECC","#CCECCC","#409188","#0077b6")
names(COLOR_BF) <- unique(num_assoc$CATEGORY)
COLOR_FDR <- c("#FCCABF", "#E866B0", "#FFC597", "#988FEA","#93E0E6", "#E3F5E6", "#9AC0BC", "#66B0E0")
names(COLOR_FDR) <- unique(num_assoc$CATEGORY)

fill_mapping <- setNames(c(COLOR_BF,COLOR_FDR),
    c(paste(names(COLOR_BF), "BF", sep = "_"),paste(names(COLOR_FDR), "FDR", sep = "_")))


plot_data <- num_assoc %>%
    mutate(fill_key = paste(CATEGORY, Significance, sep = "_"))

plot_data$fill_key<-factor(plot_data$fill_key, levels=rev(plot_data$fill_key))



p<-ggplot(plot_data, aes(y = CATEGORY, x = Count, fill = fill_key)) +
    geom_bar(
        stat = "identity",
        position = "stack",
        color = "black",  
        linewidth = 0.5   
    ) +
	geom_text(
        aes(label = Count),  
        position = position_stack(vjust = 0.5), 
        color = "black",     
        size = 3.5           
    ) +
    scale_fill_manual(
        values = fill_mapping ,
        labels = c("FDR" = "FDR","BF" = "BF"),
        breaks = c("T4_FDR","T4_BF"),
        name = "Significance"
    ) +
    labs(
        y = "Cell Type (Category)",
        x = "Number of Associations (Count)",
        title = "Stacked Bar Chart of Associations by Cell Type and Significance"
    ) +
    theme_test() +
    guides(fill = guide_legend(override.aes = list(fill = c(COLOR_FDR[1],COLOR_BF[1]),
                                                   color = "black")))


ggsave("Number.of.assocaitions.bar.pdf",p,width=4,height=4)