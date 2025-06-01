# This is the way used by OncoHarmony
# Rscript file_exists/他给的生成方法/generate_hg38_bed.R 
library(dplyr)
library(data.table)
library(stringr)
# gtf <- rtracklayer::import('file_exists/ref_annot.gtf.gz')
# gtf_df <- as.data.frame(gtf)

# mamba install -c bioconda -c conda-forge singularity
# https://github.com/NBISweden/AGAT
# conda package is also available
# mamba install -c bioconda -c conda-forge agat
# Another GFF Analysis Toolkit (AGAT) - Version: v1.4.0     
# mamba install -c bioconda -c conda-forge bedops

# Add introns to gtf, convert to gff3
# run under temp dir
#
# agat_sp_add_introns.pl -g ~/fusion/file_exists/ref_annot.gtf.gz -o genes.INTRONS.gff3
# gff2bed < genes.INTRONS.gff3 > genes.INTRONS.agat.bed

# $ cut -f 8 genes.INTRONS.agat.bed | sort | uniq -c
#  876858 CDS
# 1628304 exon
#   61583 gene
# 1377973 intron
#       4 RNA
#     128 Selenocysteine
#   96798 start_codon
#   90604 stop_codon
#  250327 transcript
#  375027 UTR

# 实际上这个文件是包含所有的
total.introns.bed <- fread(file="file_exists/temp/genes.INTRONS.agat.bed", header = FALSE, stringsAsFactors = F, sep="\t", na.strings = "",data.table = F)
colnames(total.introns.bed) <- c("chr","start","end","gene_id","tmp","strand","gene_biotype","type","V9","description")
total.introns.bed$transcript_id <- gsub("\\;.*","",str_split_fixed(total.introns.bed$description,"transcript_id=",n=2)[,2])
total.introns.bed$gene_name <-gsub("\\;.*","",str_split_fixed(total.introns.bed$description,"gene_name=",n=2)[,2])

transcript_ids <- unique(total.introns.bed$transcript_id)
NA %in% transcript_ids
"" %in% transcript_ids
head(transcript_ids)
nchar(transcript_ids[1])
head(total.introns.bed$gene_id)
head(total.introns.bed)

total.introns.bed$gene_id = substr(total.introns.bed$gene_id, 1, 15)
total.introns.bed$transcript_id = substr(total.introns.bed$transcript_id, 1, 15)
transcript_ids <- setdiff(unique(total.introns.bed$transcript_id), "")
transcript_ids[1]
all(nchar(transcript_ids) == 15)
transcript_ids[nchar(transcript_ids) != 15 | !startsWith(transcript_ids, "ENST")]

total.introns.bed[total.introns.bed$transcript_id == "ENSG00000260596", ]
transcript_ids = setdiff(transcript_ids, transcript_ids[nchar(transcript_ids) != 15 | !startsWith(transcript_ids, "ENST")])
transcript_ids[nchar(transcript_ids) != 15 | !startsWith(transcript_ids, "ENST")]

file.to_write <- "file_exists/temp/cleaned_metafusion_hg38_gene.bed"

if(file.exists(file.to_write) ) {file.remove(file.to_write)}

#START CLOCK: THE INDEXING TAKES A LONG TIME, LIKE 5 HOURS
ptm <- proc.time()

# Index each transcript feature, incrementing when an intron is passed
## metafusion expects exon count 0 to (N(exons)-1)
## Forward strand: Exon 0 == Exon 1
### Reverse strand: Exon 0 == LAST EXON IN TRANSCRIPT
for (i in 1:length(transcript_ids)){
  id = transcript_ids[i]
  message("handling ", id, ", handled ", i, "/", length(transcript_ids))
  transcript <- total.introns.bed[total.introns.bed$transcript_id == id,]
  # Remove exons if coding gene, since "exon" and "CDS" are duplicates of one another
  if ("CDS" %in% transcript$type){
    transcript <- transcript[!transcript$type == "exon",]
  }
  # Order features by increasing bp 
  transcript <- transcript[order(transcript$start, decreasing = FALSE),]
  # Index features
  idx <- 0
  for (i in 1:nrow(transcript)){
    transcript$idx [i]<- idx
    if (transcript$type[i] == "intron"){
      idx <- idx + 1
    }
  }
  # REFORMAT TRANSCRIPT
  #Change strand info (+ --> f, - --> r)
  if (unique(transcript$strand) == "+"){
    transcript$strand <- 'f'
  } else if  (unique(transcript$strand) == "-"){
    transcript$strand <- 'r'
  } else {
    errorCondition("Strand info for this transcript is inconsistent")
  }
  #Add "chr" prefix to chromosomes
  #transcript$chr <- sapply("chr", paste0,  transcript$chr)
  #Change CDS --> cds ### IF A TRANSCRIPT LACKS "CDS" THIS LINE WILL DO NOTHING, Changing exon values to UTRs later 
  if ("CDS" %in% unique(transcript$type)){transcript[transcript$type == "CDS",]$type <- "cds"}
  
  ## DETERMING UTR3 and UTR5
  
  
  ### INSTEAD OF START AND STOP, USE CDS LOCATIONS AND STRAND INFORMATION.....
  if ("UTR" %in% unique(transcript$type)){
    if( unique(transcript$strand) == "f"){
      #Forward strand 
      start_coding <- min(transcript[transcript$type == "cds","start"])
      stop_coding <-  max(transcript[transcript$type == "cds","end"])
      transcript$type[transcript$end <= start_coding &  transcript$type == "UTR"] <- "utr5"
      transcript$type[transcript$start >= stop_coding & transcript$type == "UTR"] <- "utr3"
    }else {
      start_coding <- max(transcript[transcript$type == "cds","end"])
      stop_coding <- min(transcript[transcript$type == "cds","start"])
      transcript$type[transcript$end <= start_coding &  transcript$type == "UTR"] <- "utr3"
      transcript$type[transcript$start >= stop_coding & transcript$type == "UTR"] <- "utr5"
    }
    
  }
  
  transcript <- transcript[,c("chr", "start", "end", "transcript_id", "type", "idx", "strand", "gene_name", "gene_id" )]
  write.table(transcript, file.to_write, append=TRUE, sep="\t", quote=F,  row.names=F, col.names=F)
}

time <- proc.time() - ptm
print(time)
# 
# user    system   elapsed 
# 16657.116    32.227 16741.382 

new.bed <- fread(file.to_write,data.table = F)
colnames(new.bed) <- c("chr","start","end","transcript_id","type","idx","strand","gene_name","gene_id")

#### Any exon that remains after the cds change, is likely and untranslated region. change below

# Basically, subfeatures which are "exon" need to be changed (i.e. exon --> utr3/utr5)
#Forward strand
new.bed$type[new.bed$strand == "f" &  new.bed$type == "exon" ] <- "utr5"
#Reverse strand
new.bed$type[new.bed$strand == "r" &  new.bed$type == "exon"]<- "utr3"
expected_types <- c("cds","intron","utr3","utr5")
new.bed.ready <- new.bed[new.bed$type %in% c(expected_types),]

write.table(new.bed.ready, "file_exists/hg38_gene.bed",  sep="\t", quote=F,  row.names=F, col.names=F)
