## ---------------------------
##
## Script name: RCODE.2nd_CRAFT_seq.process_CRISPResso.r
##
## Original file: 250806.proccess.collapse.gDNA.log.r
##
## Purpose of script: Process CRISPResso2 output from the 2nd CRAFT-seq experiment (finemapping:
##                    E1/E2/E3 enhancers and pQTL/KO), classify alleles using matchPattern,
##                    and summarize per-cell allele counts
##
## Date Created: 2025-08-06
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
library(Biostrings)
library(patchwork)
library(tidyplots)

argv=commandArgs(T)
input_file=argv[1]
index_sorting=argv[2]

if(FALSE){
input_file="/path/to/gDNA_all.tsv.gz"
index_sorting="/path/to/index.sorting.tsv"
}

all_data<-read_tsv(input_file)


data<-filter(all_data, DNA_Barcode %in% c("Plate10","Plate11","Plate12"))
#start from pQTL genotypes (plate 10-12)
#extract gRNA region
query_sequence <- "aactcacagcggtcagcaagcagccccagaggacgcactgcagaggcagacgaaccatggcgaggtgagaccaggc"
KO_sequence <- "aactcgcagcggtcagcaagcagccccagaggacgcactgcagaggcagacgaaccatggcgaggtgagaccaggc"
KO2_sequence <- "aactcgcggcggtcagcaagcagccccagaggacgcactgcagaggcagacgaaccatggcgaggtgagaccaggc"
KOoff_sequence <- "aactcacggcggtcagcaagcagccccagaggacgcactgcagaggcagacgaaccatggcgaggtgagaccaggc"
VAR_sequence <- "aactcacagcggtcagcaagcagccccagaggacgcactgcagaggcagacgaaccatagcgaggtgagaccaggc"
max_mismatch <- 4


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


#collapse sequences into classes like WT allele or KO allele or pQTL allele
df_matched<-df_matched%>%
  mutate(VAR_pos=substr(Matched_Sequence, nchar(Matched_Sequence) - 17, nchar(Matched_Sequence) - 17))%>%
  mutate(KO_pos=substr(Matched_Sequence, 6, 6))%>%
  mutate(Annotation = case_when(
    Matched_Sequence == toupper(query_sequence) ~ "WT",
    Matched_Sequence == toupper(KO_sequence) ~ "KO",
    Matched_Sequence == toupper(KO2_sequence) ~ "KO2",
    Matched_Sequence == toupper(KOoff_sequence) ~ "KOoff",
    Matched_Sequence == toupper(VAR_sequence) ~ "VAR",
    KO_pos == "G"  ~ "KO.mut",
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
  geom_bar(stat="identity")+theme_test()+theme(aspect.ratio=1/1,axis.text.x = element_text(size=2,angle = 90, hjust = 1, vjust = 0.5))+
  scale_y_log10()

top20<-dplyr::slice(allele_count_total,1:20)%>%pull(Matched_Sequence)
p2<-ggplot(filter(allele_count_per_plate,Matched_Sequence%in%top20),aes(x=Matched_Sequence,y=Count,fill=Annotation))+
  geom_bar(stat="identity")+theme_test()+theme(aspect.ratio=1/1,axis.text.x = element_text(size=2,angle = 90, hjust = 1, vjust = 0.5))+
  scale_y_log10()+facet_wrap(~DNA_Barcode,ncol=3)

p<-p1+p2+plot_layout(ncol=2)
ggsave("pQTL.KO.Number.of.alleles.top20.pdf",p,width=20,height=8)

annotation_count_total<-df_matched%>%group_by(Annotation)%>%
  summarize(Count=sum(`#Reads`))%>%arrange(-Count)

annotation_count_per_plate<-df_matched%>%group_by(Annotation,DNA_Barcode)%>%
  summarize(Count=sum(`#Reads`))%>%arrange(-Count)

annotation_count_per_plate<-bind_rows(annotation_count_per_plate,mutate(annotation_count_total,DNA_Barcode="All"))
annotation_count_per_plate$Annotation<-factor(annotation_count_per_plate$Annotation, levels=c("VAR","KO","KO2","KOoff","WT","KO.mut","VAR.mut","Other"))
p<-ggplot(annotation_count_per_plate,aes(x=DNA_Barcode,y=Count,fill=Annotation))+
  geom_bar(stat = "identity", position = "fill")+theme_test()+theme(aspect.ratio=1/1,axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_hline(yintercept=0.9, linetype="dashed")+
  geom_hline(yintercept=0.95, linetype="dashed")+
  geom_hline(yintercept=0.85, linetype="dashed")

ggsave("pQTL.KO.Ratio.of.allele.annotation.pdf",p,width=8,height=8)


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

p<-ggplot(filter(true_genotype_df, Annotation=="KO"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point(size=0.1)+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=4)+geom_vline(xintercept=log10(5000))+
  geom_vline(xintercept=log10(10000))
ggsave("pQTL.KO.KO.raw.genotype.pdf",p,width=10,height=10)

p<-ggplot(filter(true_genotype_df, Annotation=="KO2"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point(size=0.1)+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=4)+geom_vline(xintercept=log10(5000))+
  geom_vline(xintercept=log10(10000))
ggsave("pQTL.KO.KO2.raw.genotype.pdf",p,width=10,height=10)

p<-ggplot(filter(true_genotype_df, Annotation=="KOoff"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point(size=0.1)+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=4)+geom_vline(xintercept=log10(5000))+
  geom_vline(xintercept=log10(10000))
ggsave("pQTL.KO.KOoff.raw.genotype.pdf",p,width=10,height=10)

p<-ggplot(filter(true_genotype_df, Annotation=="WT"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point(size=0.1)+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=4)+geom_vline(xintercept=log10(5000))+
  geom_vline(xintercept=log10(10000))
ggsave("pQTL.KO.WT.raw.genotype.pdf",p,width=10,height=10)

p<-ggplot(filter(true_genotype_df, Annotation=="VAR"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point(size=0.1)+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=4)+geom_vline(xintercept=log10(5000))+
  geom_vline(xintercept=log10(10000))
ggsave("pQTL.KO.VAR.raw.genotype.pdf",p,width=10,height=10)



p<-ggplot(filter(true_genotype_df, Annotation=="Other"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point()+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=3)+
  geom_vline(xintercept=log10(50))+geom_hline(yintercept=0.05)
ggsave("pQTL.KO.Other.corrected.genotype.pdf",p,width=10,height=5)

p<-ggplot(filter(true_genotype_df, Annotation=="KO.mut"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point()+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=3)+
  geom_vline(xintercept=log10(50))+geom_hline(yintercept=0.05)
ggsave("pQTL.KO.KO.mut.corrected.genotype.pdf",p,width=10,height=5)

p<-ggplot(filter(true_genotype_df, Annotation=="VAR.mut"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point()+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=3)+
  geom_vline(xintercept=log10(50))+geom_hline(yintercept=0.05)
ggsave("pQTL.KO.VAR.mut.corrected.genotype.pdf",p,width=10,height=5)

write_tsv(true_genotype_df, "pQTL.KO.collapsed.processed.gDNA.data.tsv")


#also do genotyping for enhancers 1-3

for(dataset in c("E1","E2","E3")){

if(dataset=="E1"){

max_mismatch <- 3
data<-filter(all_data, DNA_Barcode %in% c("Plate1","Plate2","Plate3"))

query_sequence <- "atgtcatgaggacactcagc"
VAR_sequence <- "atgtcgtgaggacactcagc"


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
  mutate(VAR_pos=substr(Matched_Sequence, 6, 6))%>%
  mutate(Annotation = case_when(
    Matched_Sequence == toupper(query_sequence) ~ "WT",
    Matched_Sequence == toupper(VAR_sequence) ~ "VAR",
    VAR_pos == "G"  ~ "VAR.mut",
    TRUE ~ "Other"
  ))

}else if(dataset=="E2"){

max_mismatch <- 4
data<-filter(all_data, DNA_Barcode %in% c("Plate4","Plate5","Plate6"))

query_sequence <- "ttcctgcccaggcttcaggccggaatctacagccctctaaag"
VAR_sequence <- "ttcctgcccaggcttcaggcaggaatctacagccctctaaag"

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
  mutate(VAR_pos=substr(Matched_Sequence, 21, 21))%>%
  mutate(Annotation = case_when(
    Matched_Sequence == toupper(query_sequence) ~ "WT",
    Matched_Sequence == toupper(VAR_sequence) ~ "VAR",
    VAR_pos == "A"  ~ "VAR.mut",
    TRUE ~ "Other"
  ))

}else if(dataset=="E3"){

max_mismatch <- 3
data<-filter(all_data, DNA_Barcode %in% c("Plate7","Plate8","Plate9"))

query_sequence <- "ggccctcccaggcctgctcttggttt"
VAR_sequence <- "ggccctcccaggcctgctcttgattt"

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
  mutate(VAR_pos=substr(Matched_Sequence, 23, 23))%>%
  mutate(Annotation = case_when(
    Matched_Sequence == toupper(query_sequence) ~ "WT",
    Matched_Sequence == toupper(VAR_sequence) ~ "VAR",
    VAR_pos == "A"  ~ "VAR.mut",
    TRUE ~ "Other"
  ))

}


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
ggsave(str_c(dataset,".Number.of.alleles.top20.pdf"),p,width=20,height=8)

annotation_count_total<-df_matched%>%group_by(Annotation)%>%
  summarize(Count=sum(`#Reads`))%>%arrange(-Count)

annotation_count_per_plate<-df_matched%>%group_by(Annotation,DNA_Barcode)%>%
  summarize(Count=sum(`#Reads`))%>%arrange(-Count)

annotation_count_per_plate<-bind_rows(annotation_count_per_plate,mutate(annotation_count_total,DNA_Barcode="All"))
annotation_count_per_plate$Annotation<-factor(annotation_count_per_plate$Annotation, levels=c("VAR","WT","VAR.mut","Other"))
p<-ggplot(annotation_count_per_plate,aes(x=DNA_Barcode,y=Count,fill=Annotation))+
  geom_bar(stat = "identity", position = "fill")+theme_test()+theme(aspect.ratio=1/1,axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  geom_hline(yintercept=0.9, linetype="dashed")+
  geom_hline(yintercept=0.95, linetype="dashed")+
  geom_hline(yintercept=0.85, linetype="dashed")

ggsave(str_c(dataset,".Ratio.of.allele.annotation.pdf"),p,width=8,height=8)


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
ggsave(str_c(dataset,".WT.raw.genotype.pdf"),p,width=10,height=10)

p<-ggplot(filter(true_genotype_df, Annotation=="VAR"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point(size=0.1)+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=4)+geom_vline(xintercept=log10(5000))+
  geom_vline(xintercept=log10(10000))
ggsave(str_c(dataset,".VAR.raw.genotype.pdf"),p,width=10,height=10)



p<-ggplot(filter(true_genotype_df, Annotation=="Other"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point()+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=3)+
  geom_vline(xintercept=log10(50))+geom_hline(yintercept=0.05)
ggsave(str_c(dataset,".Other.corrected.genotype.pdf"),p,width=10,height=5)


p<-ggplot(filter(true_genotype_df, Annotation=="VAR.mut"),aes(x=log10(Total_allele_count),y=Ratio,color=condition))+
  geom_point()+theme_bw()+theme(aspect.ratio=1/1)+facet_wrap(~DNA_Barcode,ncol=3)+
  geom_vline(xintercept=log10(50))+geom_hline(yintercept=0.05)
ggsave(str_c(dataset,".VAR.mut.corrected.genotype.pdf"),p,width=10,height=5)

write_tsv(true_genotype_df, str_c(dataset,".collapsed.processed.gDNA.data.tsv"))
}
