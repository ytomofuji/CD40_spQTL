## ---------------------------
##
## Script name: RCODE.extract.Ig.junction.reads.per.cell.log.r
##
## Purpose of script: Scan CRAFT-seq BAMs for reads spanning the CH1 splice
##   junction of each IGH constant gene, classify as productive (J-CH1) or
##   sterile (I-exon-CH1) junctions, and collapse to UMI-deduplicated counts
##   per cell barcode
##
##
## Date Created: 2026-11-13
##
## ---------------------------
##
## Notes:
##   BAMs must be STAR-aligned (STARsolo CB_UMI_Simple) with CB/UB/NH tags
##   retained (--outSAMattributes NH ... CB UB ...), coordinate-sorted and
##   indexed. One subdirectory per plate under bam_root, named
##   <experiment>_<plate>/STARoutAligned.sortedByCoord.out.bam
##
## ---------------------------

argv=commandArgs(T)
bam_root=argv[1]
out_file=argv[2]

if(FALSE){
	bam_root="/path/to/STAR_BAM"
	out_file="Ig.junction.UMI.per.cell.tsv.gz"
}

# mkdir -p /path/to/Ig_class
# cd /path/to/Ig_class

library(tidyverse)

module_cmd<-"module load SAMtools/1.22.1-GCC-13.3.0 2>/dev/null;"

experiments<-c("PE1","PE2","ABE")
plates<-paste0("Plate",1:12)
PLATES<-as.vector(outer(experiments,plates,paste,sep="_"))

#IGH constant gene CH1 exon boundary + sterile I-exon donor site (GRCh38, 10x GRCh38-2020-A genes.gtf)
genes<-tribble(
	~gene, ~ch1_end, ~i_exon_donor,
	"IGHA2", 105588395, 105591672,
	"IGHE",  105601728, 105605036,
	"IGHG4", 105626066, 105629001,
	"IGHG2", 105644790, 105647807,
	"IGHA1", 105708665, 105712643,
	"IGHG1", 105743071, 105747325,
	"IGHG3", 105771405, 105775026,
	"IGHD",  105845678, NA,
	"IGHM",  105856218, 105856714)

#functional IGHJ1-6 start coordinates (GRCh38, 10x GRCh38-2020-A genes.gtf)
IGHJ_POS<-c(105863198,105863814,105864215,105864587,105865199,105865407)

get_splice_edges<-function(pos,cigar){
	ops<-str_extract_all(cigar,"[0-9]+[MIDNSHP=X]")[[1]]
	cur<-pos
	edges<-list()
	for(op in ops){
		len<-as.numeric(str_sub(op,1,-2))
		type<-str_sub(op,-1,-1)
		if(type%in%c("M","D","=","X")){
			cur<-cur+len
		}else if(type=="N"){
			edges[[length(edges)+1]]<-c(low=cur,high=cur+len-1)
			cur<-cur+len
		}
	}
	edges
}

classify_junction<-function(pos,cigar,ch1_end,i_exon_donor){
	if(!str_detect(cigar,"N")) return(NA_character_)
	edges<-get_splice_edges(pos,cigar)
	hit<-NA_character_
	for(e in edges){
		if(abs(e["low"]-1-ch1_end)<=3){
			if(any(abs(e["high"]-IGHJ_POS)<=3)) hit<-"productive"
			else if(!is.na(i_exon_donor)&&abs(e["high"]-i_exon_donor)<=3) hit<-"sterile"
		}
	}
	hit
}

all_reads<-map_dfr(PLATES,function(plate){
	bam<-file.path(bam_root,plate,"STARoutAligned.sortedByCoord.out.bam")
	if(!file.exists(bam)) return(tibble())
	map_dfr(1:nrow(genes),function(i){
		g<-genes$gene[i]; ch1<-genes$ch1_end[i]; donor<-genes$i_exon_donor[i]
		region<-sprintf("chr14:%d-%d",ch1-5,ch1+5)
		cmd<-sprintf("%s samtools view %s '%s'",module_cmd,bam,region)
		lines<-system(cmd,intern=TRUE)
		if(length(lines)==0) return(tibble())
		fields<-str_split_fixed(lines,"\t",11)
		df<-tibble(
			pos=as.numeric(fields[,4]),
			cigar=fields[,6],
			cb=str_remove(str_extract(lines,"CB:Z:[ACGT]+"),"CB:Z:"),
			ub=str_remove(str_extract(lines,"UB:Z:[ACGTN]+"),"UB:Z:"),
			nh=as.numeric(str_remove(str_extract(lines,"NH:i:[0-9]+"),"NH:i:")))
		df<-df%>%filter(!is.na(cb),nh==1)
		df$junction<-map2_chr(df$pos,df$cigar,classify_junction,ch1_end=ch1,i_exon_donor=donor)
		df<-df%>%filter(!is.na(junction))
		if(nrow(df)==0) return(tibble())
		df%>%mutate(gene=g,plate=plate)
	})
})

all_reads<-all_reads%>%
	mutate(experiment=str_extract(plate,"^[A-Z0-9]+"),
		plate_name=str_extract(plate,"Plate[0-9]+"),
		barcode=str_c(experiment,"_",cb,"_",plate_name))

#UMI dedup: collapse to unique (barcode, gene, junction type, UMI)
umi_counts<-all_reads%>%filter(!is.na(ub),ub!="")%>%
	distinct(barcode,gene,junction,ub)%>%
	count(barcode,gene,junction,name="n_umis")

write_tsv(umi_counts,out_file)
