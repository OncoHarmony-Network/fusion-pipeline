library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tidyr)
library(stringr)
library(openxlsx)
  library(writexl)
output_dir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure2/新版/Figure2A/" 
file_meta <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
folder_path <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/"
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
    cohort == "TARGET-CELL" ~ "Cell",
    
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
  mutate(
    side5 = side_class(gene1_type),
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
    sample_count = n(),                               # 样本总数 (代替原来的 total_runs)
    detected_runs = sum(Fusion_flag, na.rm = TRUE),   # 检出融合的样本数
    Fusion_rate = mean(Fusion_flag, na.rm = TRUE),    # 总体融合率
    
    # 批量计算所有基因类型列的均值（即发生率）
    across(all_of(all_gene_types), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  # 将所有的比例数值转换为 "xx.xx%" 格式的字符串
  # 修改列名，把生成的 "CG::CG rate" 改名为 "CG::CG_rate" 以符合你的要求
  rename_with(~ gsub(" rate", "_rate", .x), all_of(all_gene_types)) %>%
  # 重命名 cancer 列
  dplyr::rename(cancer_type = cancer) %>%
  # 如果不需要 detected_runs 这一列（你的示例结果里没有），可以把它去掉
  select(-detected_runs)

write_xlsx(
  result_all ,
  path = "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/附表/TableS2.xlsx"
)
