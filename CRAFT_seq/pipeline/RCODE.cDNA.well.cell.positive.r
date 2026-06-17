## ---------------------------
##
## Script name: RCODE.cDNA.well.cell.positive.r
##
## Purpose of script: Check the position of well which were positive for cells
##
##
## Date Created: 2024-05-09
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
library(viridis)

argv=commandArgs(T)

full_input_file=argv[1]
input_file=argv[2]
barcode_file=argv[3]

if(FALSE){
full_input_file="/path/to/processed/cDNA/CD40_Daudi_RNA_all.meta.tsv"
input_file="/path/to/processed/cDNA/filtered.CD40_Daudi_RNA_all.meta.tsv"
barcode_file="/path/to/well_DNA_RNA_barcode.tsv"
}



data<-read_tsv(input_file)
raw_data<-read_tsv(full_input_file)
barcode<-read_tsv(barcode_file)

plate_df<-data %>% dplyr::select(cell) %>%
separate(cell,into=c("RNA_barcode","plate")) %>%
inner_join(barcode,by="RNA_barcode")%>%
  extract(name, into = c("row", "column"), regex = "([A-Za-z]+)([0-9]+)", convert = TRUE)%>%
  mutate(is_cell=TRUE)


raw_plate_df<-raw_data %>% dplyr::select(cell) %>%
separate(cell,into=c("RNA_barcode","plate")) %>%
inner_join(barcode,by="RNA_barcode")%>%
  extract(name, into = c("row", "column"), regex = "([A-Za-z]+)([0-9]+)", convert = TRUE)%>%
  dplyr::select(plate,row,column)

plate_list <- unique(raw_plate_df$plate)
row_list <- unique(raw_plate_df$row)
column_list <- unique(raw_plate_df$column)

plt_df <- expand_grid(
  plate = plate_list,
  row = row_list,
  column = column_list
)%>%left_join(plate_df,by=c("plate","row","column"))


plt_df$is_cell[is.na(plt_df$is_cell)]<-FALSE
plt_df$plate<-factor(plt_df$plate, levels=unique(plt_df$plate))
plt_df$row<-factor(plt_df$row, levels=sort(unique(plt_df$row)))
plt_df$column<-factor(plt_df$column, levels=sort(unique(plt_df$column)))

p<-ggplot(plt_df, aes(x=column,y=row,fill=is_cell))+
  geom_tile()+theme_bw()+scale_fill_manual(values=c("TRUE"="#D8031C","FALSE"="#01016F"))+coord_fixed()+
  facet_wrap(~plate,ncol=3)

ggsave("cell.pos.neg.per.plate.pdf",p,width=12,height=12)

plt_df2<-plt_df%>%group_by(row,column)%>%summarize(Count=sum(is_cell))
p<-ggplot(plt_df2, aes(x=column,y=row,fill=Count))+coord_fixed()+
  geom_tile()+theme_bw()+scale_fill_viridis(option="inferno")

ggsave("cell.pos.neg.all.plate.pdf",p,width=5,height=5)