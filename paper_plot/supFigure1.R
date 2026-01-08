### Figure2a
library(sankeyD3)
library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(tidyr)
library(stringr)
library(sankeyD3)
metafusion <- fread("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv" , sep = "\t") 
cohort1 <- fread("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/cohort_meta_check.csv")
set.seed(2024)
coding_gene <- c( "protein_coding")
lncRNA_gene <- c("lncRNA")
other_noncoding_RNA_gene <- c("misc_RNA",  "rRNA", "TEC","snoRNA", "ncRNA","snRNA","miRNA" )
immune_gene <- c("IG_D_gene","TR_J_gene", "TR_V_gene", "TR_C_gene")
pseudogene <- c("transcribed_processed_pseudogene" , "transcribed_unprocessed_pseudogene",
                "unprocessed_pseudogene", "processed_pseudogene","pseudogene", "transcribed_unitary_pseudogene",
                "unitary_pseudogene","polymorphic_pseudogene","translated_unprocessed_pseudogene",
                "IG_V_pseudogene",  "rRNA_pseudogene","TR_V_pseudogene" )
#CG : protein_coding
#LG : lncRNA
#NG : misc_RNA, TEC, rRNA, snoRNA, ncRNA, snRNA, miRNA
#IG : IG_D_gene, TR_J_gene, TR_V_gene, TR_C_gene
#PG : pseudogene
tumor_markers <- c("Mitelman", "chimerdb_omim", "chimerdb_pubmed", "Cosmic", 
                   "YOSHIHARA_TCGA", "ChimerKB", "ChimerPub", "ChimerSeq", 
                   "Klijn_CellLines", "Larsson_TCGA", "CCLE_StarF2019", 
                   "HaasMedCancer", "TCGA_StarF2019", "TumorFusionsNAR2018", 
                   "DEEPEST2019", "GUO2018CR_TCGA", "DepMap2023")
normal_markers <- c("GTEx_recurrent_StarF2019", "DGD_PARALOGS", "HGNC_GENEFAM", "BodyMap", "Greger_Normal", "ConjoinG", "Babiceanu_Normal")
tumor_pattern <- paste0("\\b", paste(tumor_markers, collapse = "\\b|\\b"), "\\b")
normal_pattern <- paste0("\\b", paste(normal_markers, collapse = "\\b|\\b"), "\\b")
side_class <- function(x) {
  vapply(strsplit(as.character(x), ",", fixed = TRUE), function(t) {
    t <- trimws(t)
    if (length(t) == 0 || all(is.na(t) | t == "")) return(NA_character_)
    if (any(t == "protein_coding"))          return("CG")
    if (any(t == "lncRNA"))                  return("LG")
    if (any(grepl("^(IG|TR)_", t)))          return("IG")
    if (any(grepl("pseudogene", t)))         return("PG")
    "NG"
  }, character(1))
}

hit_any <- function(x, pats) {
  x <- as.character(x); x[is.na(x)] <- ""
  pats <- unique(pats[!is.na(pats) & nzchar(pats)])
  if (!length(pats)) return(rep(FALSE, length(x)))
  Reduce(`|`, lapply(pats, function(p) grepl(p, x, ignore.case = TRUE)))
}
fusion_distribution_well <- metafusion |>
  mutate(
    side5 = side_class(gene1_type),
    side3 = side_class(gene2_type),
    gene1_type_sum = paste0("5' ", side5),
    gene2_type_sum = paste0("3' ", side3),
    hit_t = hit_any(cancer_db_hits, tumor_pattern),
    hit_n = hit_any(cancer_db_hits, normal_pattern),
    reported = case_when(
      hit_t & !hit_n                              ~ "Tumor only",
      !hit_t & hit_n                              ~ "Control only",
      hit_t & hit_n                               ~ "Both",
      is.na(cancer_db_hits) | cancer_db_hits == "" ~ "Not reported",
      TRUE                                        ~ "Other"
    )
  )

## ------------------------------------------------------------------
## 2) 逐列计数（相邻两列构成一段 flow）
## ------------------------------------------------------------------
sample_order <- c("Tumor", "Tumor_Blood", "Tumor_Cell",
                  "Normal", "Normal_Adjacent", "Normal_Blood")

d1 <- fusion_distribution_well |>
  count(sample_type, inferred_fusion_type, name = "value") |>
  arrange(match(sample_type, sample_order), -value) |>
  select(source = sample_type, target = inferred_fusion_type, value)

d2 <- fusion_distribution_well |>
  count(inferred_fusion_type, Score, name = "value") |>
  arrange(-value) |>
  select(source = inferred_fusion_type, target = Score, value)

d3 <- fusion_distribution_well |>
  count(Score, reported, name = "value") |>
  arrange(-value) |>
  select(source = Score, target = reported, value)

d4 <- fusion_distribution_well |>
  count(reported, gene1_type_sum, name = "value") |>
  arrange(-value) |>
  select(source = reported, target = gene1_type_sum, value)

d5 <- fusion_distribution_well |>
  count(gene1_type_sum, gene2_type_sum, name = "value") |>
  arrange(-value) |>
  select(source = gene1_type_sum, target = gene2_type_sum, value)

## 各段统一为 character/整数，避免 rbind 类型冲突
fix <- function(d) d |> mutate(source = as.character(source),
                               target = as.character(target),
                               value  = as.integer(value))
links <- rbind(fix(d1), fix(d2), fix(d3), fix(d4), fix(d5)) |> as.data.frame()

## ------------------------------------------------------------------
## 3) 节点表 + 按名字取色（不再按位置硬编码，避免错位）
## ------------------------------------------------------------------
nodes <- data.frame(
  name = unique(c(links$source, links$target)),
  stringsAsFactors = FALSE
)

pal <- c(
  "Not reported" = "#707b89",
  "Tumor only"   = "#bf1d2d",
  "Control only" = "#2E7D32",
  "Both"         = "#293890",
  "Tumor" = "#7C000C", "Tumor_Blood" = "#E60D0C", "Cell Model" = "#FACCBD",
  "Normal" = "#43A047FF", "Normal_Adjacent" = "#1B5E20FF", "Normal_Blood" = "#CDDC39FF",
  "CodingFusion" = "#6A6599FF", "TruncatedCoding" = "#B24745FF",
  "TruncatedNoncoding" = "#00A1D5FF", "NoHeadGene" = "#374E55FF", "SameGene" = "#707b89",
  "0" = "#D0E8F2", "1" = "#ADD0E1", "2" = "#5F9BBC", "3" = "#19506F",
  "5' CG" = "#FFB703", "5' LG" = "#CD5582", "5' IG" = "#FF7F00",
  "5' NG" = "#898CF3", "5' PG" = "#5D2B84",
  "3' CG" = "#FFB703", "3' LG" = "#CD5582", "3' IG" = "#FF7F00",
  "3' NG" = "#898CF3", "3' PG" = "#5D2B84"
)
nodes$color <- unname(pal[nodes$name])
nodes$color[is.na(nodes$color)] <- "#CCCCCC"

## 索引转 0-based（sankeyD3 要求）
links$source <- match(links$source, nodes$name) - 1L
links$target <- match(links$target, nodes$name) - 1L

## ------------------------------------------------------------------
## 4) 出图
## ------------------------------------------------------------------
p <- sankeyD3::sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
                             Target = "target", Value = "value", NodeID = "name",
                             height = 1000, width = 1500,
                             NodeColor = "color",
                             linkColor = "#A0A0A0",
                             numberFormat = ".0f",
                             fontSize = 18,
                             nodeWidth = 60,
                             nodePadding = 30)
p
saveNetwork(p,"C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图1/sankey.html")
