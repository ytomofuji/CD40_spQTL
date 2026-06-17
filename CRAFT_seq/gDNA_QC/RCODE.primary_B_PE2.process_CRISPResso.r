## ---------------------------
##
## Script name: RCODE.primary_B_PE2.process_CRISPResso.r
##
## Original file: 251130.proccess.collapse.gDNA.log.r
##
## Purpose of script: Process CRISPResso2 output from the PE2 CRAFT-seq experiment
##                    (primary B cells), classify alleles using matchPattern,
##                    and summarize per-cell allele counts
##
## Date Created: 2025-11-30
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
input_file=argv[1]
index_sorting=argv[2]

if(FALSE){
input_file="/path/to/gDNA_all.tsv.gz"
index_sorting="/path/to/index.sorting.tsv"
}

library(tidyverse)
library(Biostrings)
library(patchwork)
library(tidyplots)



#extract gRNA region
query_sequence <- "ctgcagaggcagacgaaccatggcgaggtgagaccaggc"
VAR_sequence <- "ctgcagaggcagacgaaccatagcgaggtgagaccaggc"

max_mismatch <- 4

data<-read_tsv(input_file)

# Apply matchPattern
df_matched <- data %>%
  mutate(
    match_result = map(Aligned_Sequence, ~{
      subject_seq <- DNAString(.x)
      match_positions <- matchPattern(query_sequence, subject_seq, max.mismatch = max_mismatch, with.indels = TRUE)

      if (length(match_positions) == 1) {
        tibble(
          Start = start(match_positions),
          End = end(match_positions),
          Matched_Sequence = as.character(subseq(subject_seq, start(match_positions), end(match_positions)))
        )
      } else {
        tibble(Start = NA, End = NA, Matched_Sequence = NA)
      }
    })
  ) %>%
  unnest(match_result)


#collapse sequences into classes like WT allele or pQTL allele
df_matched<-df_matched%>%
  mutate(VAR_pos=substr(Matched_Sequence, nchar(Matched_Sequence) - 17, nchar(Matched_Sequence) - 17))%>%
  mutate(KO_pos=substr(Matched_Sequence, nchar(Matched_Sequence) - 18, nchar(Matched_Sequence) - 18))%>%
  mutate(Annotation = case_when(
    Matched_Sequence == toupper(query_sequence) ~ "WT",
    Matched_Sequence == toupper(VAR_sequence) ~ "VAR",
    VAR_pos == "A"  ~ "VAR.mut",
    TRUE ~ "Other"
  ))

df_matched$Matched_Sequence[is.na(df_matched$Matched_Sequence)]<-"Diff"

allele_count_total<-df_matched%>%group_by(Matched_Sequence,Annotation)%>%
  summarize(Count=sum(`#Reads`))%>%ungroup()%>%arrange(-Count)
allele_count_total$Matched_Sequence<-factor(allele_count_total$Matched_Sequence, levels=allele_count_total$Matched_Sequence)

allele_count_per_plate<-df_matched%>%group_by(Matched_Sequence,Annotation,DNA_Barcode)%>%
  summarize(Count=sum(`#Reads`))%>%ungroup()%>%arrange(-Count)
allele_count_per_plate$Matched_Sequence<-factor(allele_count_per_plate$Matched_Sequence, levels=allele_count_total$Matched_Sequence)

p1<-ggplot(dplyr::slice(allele_count_total,1:20),aes(x=Matched_Sequence,y=Count,fill=Annotation))+
  geom_bar(stat="identity")+theme_test()+theme(aspect.ratio=1/1,axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_y_log10()

top20<-dplyr::slice(allele_count_total,1:20)%>%pull(Matched_Sequence)
p2<-ggplot(filter(allele_count_per_plate,Matched_Sequence%in%top20),aes(x=Matched_Sequence,y=Count,fill=Annotation))+
  geom_bar(stat="identity")+theme_test()+theme(aspect.ratio=1/1,axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_y_log10()+facet_wrap(~DNA_Barcode,ncol=3)

p<-p1+p2+plot_layout(ncol=2)
ggsave("Number.of.alleles.top20.pdf",p,width=20,height=8)

annotation_count_total<-df_matched%>%group_by(Annotation)%>%
  summarize(Count=sum(`#Reads`))%>%arrange(-Count)

annotation_count_per_plate<-df_matched%>%group_by(Annotation,DNA_Barcode)%>%
  summarize(Count=sum(`#Reads`))%>%arrange(-Count)

annotation_count_per_plate<-bind_rows(annotation_count_per_plate,mutate(annotation_count_total,DNA_Barcode="All"))
annotation_count_per_plate$DNA_Barcode<-factor(annotation_count_per_plate$DNA_Barcode, levels=c("All",str_c("Plate",c(seq(1,12)))))
annotation_count_per_plate$Annotation<-factor(annotation_count_per_plate$Annotation, levels=c("VAR","WT","VAR.mut","Other"))
p<-ggplot(annotation_count_per_plate,aes(x=DNA_Barcode,y=Count,fill=Annotation))+
  geom_bar(stat = "identity", position = "fill")+theme_test()+theme(aspect.ratio=1/1,axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_hline(yintercept=0.9, linetype="dashed")+
  geom_hline(yintercept=0.95, linetype="dashed")+
  geom_hline(yintercept=0.85, linetype="dashed")

ggsave("Ratio.of.allele.annotation.pdf",p,width=8,height=8)


#plot 1. ratio of the allele-annotation / total number of reads for each plate
total_allele_num<-df_matched%>%group_by(DNA_Barcode,Well_ID)%>%summarize(Total_allele_count=sum(`#Reads`))%>%ungroup()
genotype_df<-df_matched%>%group_by(DNA_Barcode,Well_ID,Annotation)%>%
  summarize(Count=sum(`#Reads`))%>%left_join(total_allele_num,by=c("DNA_Barcode","Well_ID"))%>%ungroup()%>%
  mutate(Ratio=Count/Total_allele_count)%>%
  mutate(cell_name=str_c(DNA_Barcode,"_",Well_ID))

tmp_df<-genotype_df%>%dplyr::select(cell_name,Annotation,Ratio)%>%
  pivot_wider(names_from="Annotation",values_from="Ratio")%>%
  pivot_longer(names_to="Annotation",values_to="Ratio",-cell_name)
tmp_df$Ratio[is.na(tmp_df$Ratio)]<-0

genotype_df<-genotype_df%>%dplyr::select(-c(Ratio,Annotation,Count,Ratio))%>%dplyr::distinct()%>%
  inner_join(tmp_df,by="cell_name")


facs<-read_tsv(index_sorting)%>%
  mutate(cell_name=str_c("Plate",cell_name))%>%
  mutate(Row = str_extract(well, "^[A-P]"),
         Col = str_extract(well, "[0-9]+"))


true_genotype_df<-inner_join(genotype_df,facs,by="cell_name")

p<-ggplot(filter(true_genotype_df, Annotation=="WT"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point(size=0.1)+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=4)+geom_vline(xintercept=log10(5000))+
  geom_vline(xintercept=log10(10000))
ggsave("WT.raw.genotype.pdf",p,width=10,height=10)

p<-ggplot(filter(true_genotype_df, Annotation=="VAR"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point(size=0.1)+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=4)+geom_vline(xintercept=log10(5000))+
  geom_vline(xintercept=log10(10000))
ggsave("VAR.raw.genotype.pdf",p,width=10,height=10)



p<-ggplot(filter(true_genotype_df, Annotation=="Other"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point()+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=3)+
  geom_vline(xintercept=log10(50))+geom_hline(yintercept=0.05)
ggsave("Other.corrected.genotype.pdf",p,width=10,height=5)

p<-ggplot(filter(true_genotype_df, Annotation=="VAR.mut"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point()+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=3)+
  geom_vline(xintercept=log10(50))+geom_hline(yintercept=0.05)
ggsave("VAR.mut.corrected.genotype.pdf",p,width=10,height=5)

write_tsv(true_genotype_df, "collapsed.processed.gDNA.data.tsv")
