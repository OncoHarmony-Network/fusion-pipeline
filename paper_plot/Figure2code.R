library(ComplexHeatmap)
library(ggplot2)
library(survival)
library(survminer)
library(data.table)
library(dplyr)
library(grid)
library(tidyr)
library(reticulate)
library(plotly)
library(ggpubr)
library(ggthemes)
library(ggprism)
library(ggsignif)
library(stringr)
library(patchwork)
library(reticulate)
library(htmlwidgets)
#Figure2A
#use_virtualenv("r-reticulate", required = TRUE)
use_python("C:/Program Files/Python312/python.exe", required = TRUE)
py_install(c("plotly", "kaleido"))
py_config()
output_dir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure2/新版/Figure2A/" 
file_meta <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
folder_path <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo"
file_names <- list.files(folder_path, full.names = FALSE, pattern = "_info.tsv$")
Credibility_gene <- fread(file_meta, sep = "\t") 
cohort_sum_all <- sub("_info.tsv$", "", file_names) 
cohort1 <- fread("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/cohort_meta_check.csv")

OS_all <- data.table()
for(cohort_id in cohort_sum_all ) {
  Need <- fread(paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/", cohort_id, "_info.tsv")) 
  Need1 <- Need
  columns_to_select <- c("Patient_ID", "Run", "Sample_type")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    select(all_of(existing_columns))
  
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  if(length(missing_columns) > 0) {
    for(col in missing_columns) {
      Need2[[col]] <- NA_character_
    }
  }
  
  Need2$cohort <- cohort_id
  OS_all <- rbind(Need2, OS_all)
}

OS_all_filter <- OS_all |>
  distinct() |>
  left_join(
    cohort1 %>% select(id, cancer_type),
    by = c("cohort" = "id")
  ) |>
  rename(cancer = cancer_type) |>
  mutate(cancer = case_when(
    # 1. 正常样本
    cohort %in% c("TCGA-NORMAL", "CPTAC-NORMAL") ~ "Normal Adjacent",
    cohort == "TARGET-NORMAL" & Sample_type == "Bone Marrow Normal" ~ "Normal",
    cohort == "TARGET-NORMAL" & Sample_type == "Blood Derived Normal" ~ "Normal Blood",
    cohort == "TARGET-LAML" & Sample_type == "Bone Marrow Normal" ~ "Normal",
    cohort == "TARGET-LAML" & Sample_type == "Blood Derived Normal" ~ "Normal Blood",
    
    # 2. 细胞系
    cohort == "TARGET-CELL" ~ "Cell Model",
    
    # 5. 其余 TCGA/CPTAC/TARGET 肿瘤样本
    grepl("^CPTAC-|^TARGET-|^TCGA-", cohort) ~ sub("^.*?-", "", cohort),

    TRUE ~ cancer
  )) 


coding_gene <- c( "protein_coding")
lncRNA_gene <- c("lncRNA")
other_noncoding_RNA_gene <- c("misc_RNA",  "rRNA", "TEC","snoRNA", "ncRNA","snRNA","miRNA" )
immune_gene <- c("IG_D_gene","TR_J_gene", "TR_V_gene", "TR_C_gene")
pseudogene <- c("transcribed_processed_pseudogene" , "transcribed_unprocessed_pseudogene",
                "unprocessed_pseudogene", "processed_pseudogene","pseudogene", "transcribed_unitary_pseudogene",
                "unitary_pseudogene","polymorphic_pseudogene","translated_unprocessed_pseudogene",
                "IG_V_pseudogene",  "rRNA_pseudogene","TR_V_pseudogene" )

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

Tumor_ICB_OS <- Credibility_gene |>
  as.data.table() |>
  mutate(    side5 = side_class(gene1_type),
             side3 = side_class(gene2_type),
             gene1_type_sum = paste0("5' ", side5),
             gene2_type_sum = paste0("3' ", side3)
  ) |>
  select(Sample_ID, gene1, gene2, gene1_type_sum, gene2_type_sum, cohort) |>
  mutate(gene_type = paste0(gene1_type_sum,"::", gene2_type_sum, " rate")) |>
  dplyr::rename(Run = Sample_ID) |>
  distinct()
all_samples <- OS_all_filter %>%
  select(Run, cohort, cancer) %>%   # 请根据实际列名调整（可能为 cancer 或 cohort）
  distinct()

# 2. Fusion⁺ 样本的样本-基因类型对应关系（去重）
fusion_type_info <- Tumor_ICB_OS %>%
  select(Run, cohort, gene_type) %>%
  distinct() %>%
  mutate(value = 1)

# 3. 所有基因类型组合
all_gene_types <- unique(fusion_type_info$gene_type)

# 4. 转为宽格式（每个样本一行，基因类型为列）
sample_type_wide <- fusion_type_info %>%
  pivot_wider(
    id_cols = c(Run, cohort),
    names_from = gene_type,
    values_from = value,
    values_fill = 0
  )

# 5. 合并到所有样本，缺失列补 0
all_samples <- all_samples %>%
  left_join(sample_type_wide, by = c("Run", "cohort")) %>%
  mutate(across(all_of(all_gene_types), ~ replace_na(.x, 0)))

all_samples <- all_samples %>%
  mutate(Fusion_flag = if_else(rowSums(across(all_of(all_gene_types))) > 0, 1, 0))

# 7. 按癌种分组，计算总样本数、阳性样本数、阳性率
result_all <- all_samples %>%
  group_by(cancer) %>%
  summarise(
    total_runs = n(),
    detected_runs = sum(Fusion_flag, na.rm = TRUE),
    Fusion_rate = mean(Fusion_flag, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(Fusion_rate = sprintf("%.2f%%", Fusion_rate * 100)) |>
  dplyr::rename(cancer_type = cancer)

cancer_type <- result_all$cancer_type
for(cancer_type1 in cancer_type){
  ids <- c("Root1", "Root", "nodetected", "detected")
  parents <- c("", "Root1", "Root", "Root")
  colors <- c("#ffffff", "#ffffff","#808080", "#32037d")
  data <- result_all |> filter(cancer_type == cancer_type1)
  values <- c(data$total_runs, data$total_runs, data$total_runs - data$detected_runs, data$detected_runs)
  labels <- c("","","","")
  p_figure2a_part1 <- plot_ly(ids = ids, 
                              labels = labels, 
                              parents = parents, 
                              type = 'sunburst',
                              # textfont = list(color = "#ffffff",size = 166),
                              # insidetextorientation = "radial",
                              values = values, 
                              sort = FALSE,
                              branchvalues = 'total',
                              rotation = 90,
                              marker = list(colors = colors))
  html_file <- paste0(output_dir, "/", cancer_type1, "_Figure2a.html")
  saveWidget(p_figure2a_part1, html_file, selfcontained = TRUE)
  
  
}
edge_path <- "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"

html_files <- list.files(output_dir, pattern = "\\.html$", full.names = TRUE)

for(html_file in html_files){
  png_file <- sub("\\.html$", ".png", html_file)
  # 核心命令与 Chrome 完全一致[reference:4][reference:5]
  cmd <- sprintf('"%s" --headless --disable-gpu --screenshot="%s" --window-size=1000,1000 --hide-scrollbars "file:///%s"',
                 edge_path, png_file, html_file)
  system(cmd)
  cat("Converted:", basename(png_file), "\n")
}


### Figure2B
data_figure2_pre <- fread("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv" , sep = "\t") |>
  as.data.frame()   |>
  select(Fusion_ID, Sample_ID, cohort) |>
  left_join(
    OS_all_filter |>
      select(Run, cancer) |>
      distinct(),          # 防止一个 Run 对应多个 cancer 导致行膨胀
    by = c("Sample_ID" = "Run")
  )|>
  # 4) 每个 cancer 的独特 fusion 数
  group_by(cancer) |>
  summarise(
    n_fusion = n_distinct(Fusion_ID),
    .groups = "drop"
  ) |>
  arrange(desc(n_fusion))

p1 <- ggplot(data= data_figure2_pre,aes(x = reorder(cancer, n_fusion , FUN = sum), y = n_fusion )) +
  geom_bar(stat="identity",position = "dodge", fill = "#1c3c63")+
  theme_classic() +
  theme(  text = element_text(size = 14),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
          axis.text.y = element_text(size = 14),  
          axis.title.y = element_text(size = 24),
          legend.title = element_blank()
  )+
  guides(fill = "none") +
  labs(x=NULL,y= "Number of chimeric transcripts") 
p1

ggsave("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure2/新版/Figure2B.pdf", 
       plot = p1,         
       width = 10,                
       height = 7,              
       dpi = 600,                
       units = "in",              
       device = "pdf",           
       bg = "white")  
ggsave("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure2/新版/Figure2B.png", 
       plot = p1,         
       width = 10,                
       height = 7,                  
       dpi = 600,                
       units = "in",              
       device = "png",           
       bg = "white")  


### Figure2c
Tumor_data <- fread("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv" , sep = "\t") |>
  filter(sample_type == "Tumor" | sample_type == "Tumor_Blood")  |> 
  select(c(Fusion_ID,cancer_db_hits)) |>
  distinct() 

group_chimer      <- c("ChimerKB", "ChimerPub", "ChimerSeq", "chimerdb_omim", "chimerdb_pubmed")
group_tcga        <- c("DEEPEST2019", "YOSHIHARA_TCGA", "Larsson_TCGA", "TCGA_StarF2019", "TumorFusionsNAR2018", "GUO2018CR_TCGA")
group_ccle        <- c("Klijn_CellLines", "CCLE_StarF2019", "DepMap2023")
group_karyotype   <- "Mitelman"
group_reviews     <- c("Cosmic", "HaasMedCancer")
tumor_markers <- c("Mitelman", "chimerdb_omim", "chimerdb_pubmed", "Cosmic", "YOSHIHARA_TCGA", 
                   "ChimerKB", "ChimerPub", "ChimerSeq", "Klijn_CellLines", "Larsson_TCGA", 
                   "CCLE_StarF2019", "HaasMedCancer", "TCGA_StarF2019", "TumorFusionsNAR2018", 
                   "DEEPEST2019", "GUO2018CR_TCGA", "DepMap2023")
normal_markers <- c("GTEx_recurrent_StarF2019", "DGD_PARALOGS", "HGNC_GENEFAM", "BodyMap", "Greger_Normal", "ConjoinG", "Babiceanu_Normal")
gene_up_tumor <- Tumor_data |>
  distinct() |> 
  mutate(
    Group_TCGA = ifelse(rowSums(sapply(group_tcga, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_CCLE = ifelse(rowSums(sapply(group_ccle, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Karyotype = ifelse(rowSums(sapply(group_karyotype, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Chimer = ifelse(rowSums(sapply(group_chimer, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Reviews = ifelse(rowSums(sapply(group_reviews, grepl, x = cancer_db_hits)) > 0, 1, 0),
    Group_Not_Reported = as.integer(trimws(cancer_db_hits) == ""),
    Group_Normal = ifelse(rowSums(sapply(normal_markers, grepl, x = cancer_db_hits)) > 0, 1, 0)
  )  
TCGA <- gene_up_tumor |> 
  select(Fusion_ID, Group_TCGA) |>
  filter(Group_TCGA == 1) |>
  rename(Group = Group_TCGA)|>
  mutate(Group = "Group_TCGA") 
CCLE <- gene_up_tumor |> 
  select(Fusion_ID, Group_CCLE)|>
  filter(Group_CCLE == 1) |>
  rename(Group = Group_CCLE) |>
  mutate(Group = "Group_CCLE") 
Karyotype <- gene_up_tumor |> 
  select(Fusion_ID, Group_Karyotype)|>
  filter(Group_Karyotype == 1) |>
  rename(Group = Group_Karyotype) |>
  mutate(Group = "Group_Karyotype") 
Chimer <- gene_up_tumor |> 
  select(Fusion_ID, Group_Chimer )|>
  filter(Group_Chimer == 1) |>
  rename(Group = Group_Chimer) |>
  mutate(Group = "Group_Chimer") 
Reviews <- gene_up_tumor |> 
  select(Fusion_ID, Group_Reviews)|>
  filter(Group_Reviews == 1) |>
  rename(Group = Group_Reviews) |>
  mutate(Group = "Group_Reviews") 
Not_Reported <- gene_up_tumor |> 
  select(Fusion_ID, Group_Not_Reported) |>
  filter(Group_Not_Reported == 1) |>
  rename(Group = Group_Not_Reported) |>
  mutate(Group = "Group_Not_Reported") 
Normal <- gene_up_tumor |> 
  select(Fusion_ID, Group_Normal) |>
  filter(Group_Normal == 1) |>
  rename(Group = Group_Normal) |>
  mutate(Group = "Group_Not_Tumour") 

Group_All <- rbind(Normal,TCGA, CCLE, Karyotype, Chimer, Reviews,  Not_Reported)
labels <- c("Group_Karyotype","Group_Not_Tumour","Group_Reviews","Group_CCLE","Group_Chimer","Group_TCGA", "Group_Not_Reported")
Group_All <- Group_All %>%
  mutate(Group = factor(Group, levels = labels)) %>%
  arrange(Group)
lists <- lapply(split(Group_All, Group_All$Group), function(x) x$Fusion_ID)
m <- make_comb_mat(lists)
cs <- comb_size(m) 
colors <- c("#f3a59a", "#d7c9e2", "#a6ddea", "#80d0c3", "#9eaac4", "#c2c8da", "#d8cdc1")
outdir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure2/新版"

make_upset <- function() {
  p <- draw(ComplexHeatmap::UpSet(
    m,
    pt_size         = unit(6, "mm"),
    lwd             = 2.4,
    bg_col          = rev(colors),
    set_order       = labels,
    row_names_gp    = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 14),
    top_annotation  = upset_top_annotation(
      m,
      height               = unit(6, "cm"),            # y 轴变长
      ylim                 = c(0, max(cs) * 1.6),      # 留白放数字
      bar_width            = 0.5,                      # x 方向间隙变大
      gp                   = gpar(fill = "black", col = NA),
      axis_param           = list(gp = gpar(fontsize = 14), labels_rot = 0),
      annotation_name_gp   = gpar(fontsize = 14),
      annotation_name_rot  = 90,
      annotation_name_side = "left"
    )
  ))
  od <- column_order(p)
  decorate_annotation("intersection_size", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
              default.units = "native", just = "left", rot = 90,
              gp = gpar(fontsize = 14))
  })
}

pdf(file.path(outdir, "Figure2c.pdf"), width = 12, height = 7); make_upset(); dev.off()
png(file.path(outdir, "Figure2c.png"), width = 12, height = 7, units = "in", res = 600); make_upset(); dev.off()