library(sankeyD3)
library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(tidyr)
#supfig1c
metafusion <- fread("C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv" , sep = "\t") |>
  separate_rows(gene1, gene1_type, sep = ",") |>
  separate_rows(gene2, gene2_type, sep = ",") 

coding_gene <- c( "protein_coding")
lncRNA_gene <- c("lncRNA")
other_noncoding_RNA_gene <- c("misc_RNA",  "rRNA", "TEC","snoRNA", "ncRNA","snRNA","miRNA" )
immune_gene <- c("IG_D_gene","TR_J_gene", "TR_V_gene", "TR_C_gene")
pseudogene <- c("transcribed_processed_pseudogene" , "transcribed_unprocessed_pseudogene",
                "unprocessed_pseudogene", "processed_pseudogene","pseudogene", "transcribed_unitary_pseudogene",
                "unitary_pseudogene","polymorphic_pseudogene","translated_unprocessed_pseudogene",
                "IG_V_pseudogene",  "rRNA_pseudogene","TR_V_pseudogene" )

DB_CANCER_HITS <- c(
  "Mitelman", "chimerdb_omim", "chimerdb_pubmed", "Cosmic",
  "YOSHIHARA_TCGA", "ChimerKB", "ChimerPub", "ChimerSeq",
  "Klijn_CellLines", "Larsson_TCGA", "CCLE_StarF2019",
  "HaasMedCancer", "TCGA_StarF2019", "TumorFusionsNAR2018",
  "DEEPEST2019", "GUO2018CR_TCGA", "DepMap2023"
)
fusion_distribution_well <-metafusion |>
  mutate(
    sample_type = case_when(
      sample_type == "Blood" ~ "Tumor",  
      TRUE ~ sample_type),
    annotate = case_when(
      cancer_db_hits == "" ~ "Not Reported",
      sapply(strsplit(cancer_db_hits, ","), function(x) any(trimws(x) %in% DB_CANCER_HITS)) ~ "Reported Tumor Fusion",
      TRUE ~ "Reported Normal Fusion"
    ),
    gene1_type_sum = case_when(
      gene1_type %in% coding_gene ~ "5' CG",
      gene1_type %in% lncRNA_gene ~ "5' LG",
      gene1_type %in% other_noncoding_RNA_gene ~ "5' NG",
      gene1_type %in% immune_gene ~ "5' IG",
      gene1_type %in% pseudogene ~ "5' PG",
      TRUE ~ "Unknown"
    ),
    gene2_type_sum = case_when(
      gene2_type %in% coding_gene ~ "3' CG",
      gene2_type %in% lncRNA_gene ~ "3' LG",
      gene2_type %in% other_noncoding_RNA_gene ~ "3' NG",
      gene2_type %in% immune_gene ~ "3' IG",
      gene2_type %in% pseudogene ~ "3' PG",
      TRUE ~ "Unknown"
    )
  )

data_Fusion_type <- fusion_distribution_well |> 
  count(sample_type, inferred_fusion_type) |>  
  arrange(n)
colnames(data_Fusion_type) <- c("source", "target", "value")

data_Gene_score <- fusion_distribution_well |> 
  count(inferred_fusion_type, Score) |>  
  arrange(n)
colnames(data_Gene_score) <- c("source", "target", "value")

data_annotate <- fusion_distribution_well |> 
  count(Score, annotate) |>  
  arrange(n)
colnames(data_annotate) <- c("source", "target", "value")

data_5 <- fusion_distribution_well |> 
  count(annotate, gene1_type_sum) |>  
  arrange(n)
colnames(data_5) <- c("source", "target", "value")

data_3 <- fusion_distribution_well |> 
  count(gene1_type_sum, gene2_type_sum) |>  
  arrange(n)
colnames(data_3) <- c("source", "target", "value")

links  <- rbind(data_Fusion_type,
  data_Gene_score,
  data_annotate, 
  data_5,
  data_3) |>
  as.data.table()
nodes <- data.frame(name = unique(c(links$source, links$target)))
links$source <- match(links$source, nodes$name)-1
links$target <- match(links$target, nodes$name)-1
nodes$color<- c("#bf1d2d", "#293890", 
                "#374E55FF","#DF8F44FF","#00A1D5FF","#B24745FF","#6A6599FF",
                "#80CDC1","#A6611A","#DFC27D","#018571",
                "#eeca40", "#fd763f", "#23bac5",
                "#631879FF","#008280FF", "#3B4992FF","#008B45FF", "#EE0000FF",
                "#008280FF","#631879FF", "#EE0000FF","#008B45FF", "#3B4992FF" 
)
links$color <- c()
p <- sankeyD3::sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
                             Target = "target", Value = "value", NodeID = "name",
                             height=1000,width=1500, NodeColor = "color",
                             linkColor = "#A0A0A0",
                             #height=6000, width=6000, 
                             units = 'TWh', numberFormat=".0f", 
                             fontSize = 18,  
                             nodeWidth = 60,  
                             nodePadding = 18  )
p
saveNetwork(p,"C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/supfig1/sankey.html")