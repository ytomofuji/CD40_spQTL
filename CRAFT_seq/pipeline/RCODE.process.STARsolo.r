## ---------------------------
##
## Script name: RCODE.process.STARsolo.r
##
## Purpose of script: Process/merge STARsolo output
##
##
## Date Created: 2025-05-08
##
## Email:
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


source("misc_funcs.R")
source("utils.R")
library(stringr)
library(Matrix)

argv=commandArgs(T)

layout_file=argv[1]
out_name=argv[2]

if(FALSE){
layout_file="/path/to/processed/cDNA/layout.file.tsv"
out_name="CD40_Daudi"
}


layout<-read_tsv(layout_file, col_names=F)

exprs = do.call(cbind, lapply(seq(1,nrow(layout)), function(i){
    name=layout$X1[i]
    run=layout$X2[i]
    exprs = read10x_mtx(run, name)
})
)

meta = data.frame(cell = colnames(exprs),
                  plate = exprs %>% colnames %>% str_split('_') %>% map(2) %>% unlist)
rownames(meta) = colnames(exprs)

saveRDS(exprs, str_c(out_name,'_RNA_all.rds'))
write_tsv(meta, str_c(out_name,'_RNA_all.meta.tsv'))