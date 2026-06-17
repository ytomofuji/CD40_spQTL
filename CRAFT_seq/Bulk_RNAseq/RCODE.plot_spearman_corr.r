## ---------------------------
##
## Script name: RCODE.plot_spearman_corr.r
##
## Original file: RCODE.260502.plot.spearman.corr.log.r
##
## Purpose of script: Plot spearman correlation across time points
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





setwd("/path/to/output_directory")

library(tidyverse)

extract_spearman_results <- function(input_text) {
  s_stat <- as.numeric(str_match(input_text, "S = ([0-9.]+),") [2])
  p_value <- as.numeric(str_match(input_text, "p-value [=<]+ ([0-9.e-]+)") [2])
  rho <- as.numeric(str_match(input_text, "rho\\s+(-?[0-9.e-]+)") [2])
  
  res <- tibble(
    S_Statistic = s_stat,
    P_Value     = p_value,
    Rho         = rho,
    Correlation_Type = "Spearman"
  )
  
  return(res)
}

files <- str_c("/path/to/output_directory",c(0,4,24,72),"h.spearman_test_results.txt")
all_results <- map_dfr(files, ~ {
  txt <- read_file(.x)
  extract_spearman_results(txt) %>% mutate(File = .x)
})

all_results$timepoints<-factor(c("0h","4h","24h","72h"),levels=c("0h","4h","24h","72h"))

p<-ggplot(all_results,aes(x=timepoints,y=Rho))+
geom_point(shape=21,size=2)+theme_test()+theme(aspect.ratio=2/3)+
geom_hline(yintercept=0,linetype="dashed")+ylim(c(-0.2,0.6))

ggsave("spearman.summary.plot.pdf",p,width=2,height=2)

all_results<-mutate(all_results,signP=ifelse(Rho>0,P_Value/2,1-(P_Value/2)))
write_tsv(all_results,"spearman.summary.tsv")