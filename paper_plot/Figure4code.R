library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tidyr)
library(stringr)
library(readxl)
library(patchwork)
library(doParallel)
library(tidyverse)
library(rstatix)
library(writexl)
out_dir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4"

file_meta <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta) |>
  mutate(
    Fusion_Class = case_when(
      inferred_fusion_type %in% c("CodingFusion") ~ "Coding",
      inferred_fusion_type %in% c("TruncatedCoding", "TruncatedNoncoding") ~ "Truncated",
      inferred_fusion_type %in% c("ReadThrough", "SameGene", "Intergenic", "NoHeadGene") ~ "Passenger",
      TRUE ~ ""
    )
  )
cohort1 <- fread("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/cohort_meta_check.csv")

plot_os_by_group <- function(data, 
                             group_levels = c( "Fusion-","Fusion+"),
                             choose = "OS",
                             title = "OS", 
                             strata_var = "cohort",
                             filename = NULL) {
  
  data$event <- factor(data$event, levels = group_levels)   # ← 参照组放前面
  group_counts <- table(data$event)
  ln <- group_counts[group_levels[1]]   
  hn <- group_counts[group_levels[2]]    
  
  if(choose == "OS"){
    fit <- survfit(Surv(OS_Time, OS_Status) ~ event, data = data)
    formula_cox <- as.formula(paste0("Surv(OS_Time, OS_Status) ~ event + strata(", strata_var, ")"))
  }else{
    fit <- survfit(Surv(PFS_Time, PFS_Status) ~ event, data = data)
    formula_cox <- as.formula(paste0("Surv(PFS_Time, PFS_Status) ~ event + strata(", strata_var, ")"))
  }
  
  # 🌟 2. 核心修改：使用分层 Cox 模型统一计算 HR、CI 和 P 值
  cox_model <- coxph(formula_cox, data = data)
  cox_summary <- summary(cox_model)
  
  # 提取 P 值
  p.val1 <- cox_summary$coefficients[1, "Pr(>|z|)"] 
  p.val <- ifelse(p.val1 < 0.001, "\nP< 0.001", paste0("\nP=", round(p.val1, 3)))
  
  # 提取 HR 和 95% CI (注意：这里提取的是指数化后的值)
  HR <- round(cox_summary$coefficients[1, "exp(coef)"], 2)
  low95 <- round(cox_summary$conf.int[1, "lower .95"], 2)
  up95 <- round(cox_summary$conf.int[1, "upper .95"], 2)
  
  HR_str <- paste("Hazard Ratio =", HR)
  CI_str <- paste("95% CI:", paste(low95, up95, sep = " - "))
  
  # 构建图注
  caption1 <- paste0(
    paste0(group_levels[2], " : ", hn),
    "\n", paste0(group_levels[1], " : ", ln),
    p.val,
    "\n", HR_str,
    "\n", CI_str
  )
  data$event <- factor(data$event, levels = rev(group_levels) )
  
  if(choose == "OS"){
    fit <- survfit(Surv(OS_Time, OS_Status) ~ event, data = data)
  }else{
    fit <- survfit(Surv(PFS_Time, PFS_Status) ~ event, data = data)
  }
  p1 <- ggsurvplot(fit,
                   data = data,
                   linetype = "strata",
                   risk.table = TRUE,
                   risk.table.col = "strata",
                   legend.title = "",
                   pval = FALSE,
                   title = title,
                   font.y = c(20, "plain", "black"),
                   font.legend = c(15, "plain", "black"),
                   font.tickslab = c(15, "plain", "black"),
                   risk.table.fontsize = 4,
                   palette = c("#ed0000ff","#00468bff" ),
                   legend.labs = rev(group_levels),
                   xlab = "Time(Days)"
  )
  
  # 添加注释并调整主题
  p1$plot <- p1$plot +
    annotate("text", x = Inf, y = Inf, label = caption1, 
             vjust = 1, hjust = 1, size = 4) +
    theme(
      text = element_text(size = 14),
      legend.position = "none",
      legend.location = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      title = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(margin = margin(r = -8)), 
      plot.margin = margin(t = 5, r = 5, b = 0, l = 4)
    )
  
  # 调整 risk.table
  p1$table <- p1$table + 
    theme(
      text = element_text(size = 14),
      legend.position = "none",
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
    )
  p1$table$labels$title <- ""
  
  # 合并 plot 和 table
  p2 <- (p1$plot / p1$table) +
    plot_layout(heights = c(4, 1)) +
    plot_annotation(
      theme = theme(
        text = element_text(size = 14),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
      )
    )
  
  # 保存（如果提供了文件名）
  if (!is.null(filename)) {
    ggsave(filename = filename, plot = p2,  
           width = 6, 
           height = 5,
           dpi = 600)
    
  }
  
}

#Figure4A
cohort_OS <- cohort1 |> filter(OS == TRUE) 
cohort_sum_OS <- (cohort_OS$id) |> unique()
OS_all <- data.table()
for(cohort_id in cohort_sum_OS){
  Need <- fread(paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/",cohort_id,"_info.tsv"))
  if(filter(cohort_OS,id == cohort_id)$time == "Pre/On" ){
    if(cohort_id == "PRJNA356761"){
      Need1 <- Need |>
        filter(grepl("Pre", Sample_ID))
    }else if(cohort_id == "EGAD00001006282"){
      Need1 <- Need |>
        filter(grepl("_SCREEN", Sample_ID))
    }else{
      Need1 <- Need |>
        filter(grepl("Pre", Response2))
    }}else if(cohort_id  %in% c("PRJNA482620", "PRJNA744780")){
      Need1 <- Need |>
        filter(grepl("Pre", Response2))
    }else if(cohort_id == "PHS001464"){
      Need1 <- Need |>
        mutate(Sample_type = case_when(
          TUMOR_COLLECTION_TIMEPOINT ==  "Acquired Resistance" ~ "On",
          TUMOR_COLLECTION_TIMEPOINT ==  "Pre-Immunotherapy"~ "Pre",
          TUMOR_COLLECTION_TIMEPOINT ==  "Off Therapy Recurrence"~ "On",
          TUMOR_COLLECTION_TIMEPOINT ==  "Diagnosis"~ "Pre",
          TUMOR_COLLECTION_TIMEPOINT ==  "Mixed Response (Site of Primary Resistance)"~ "On",
          TRUE ~ NA
        ) )
    }else{
      Need1 <- Need
    }
  
  columns_to_select <- c("Patient_ID", "OS_Time", "OS_Status", "Run", "Treatment", "Age", "Sex")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    select(all_of(existing_columns)) |>
    filter(!is.na(OS_Time) & !is.na(OS_Status))
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  if (length(missing_columns) > 0) {
    for (col in missing_columns) {
      Need2[[col]] <- NA
    }
  }
  
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort <- cohort_id
  if(cohort_id == "EGAD00001003977"){
    keep <- c("Run", "Sample_ID")
    cols <- setdiff(names(Need2), keep)
    i_a <- which(Need2$Patient_ID == "10285" & Need2$Run == "EGAR00001652956")
    i_b <- which(Need2$Patient_ID == "10285" & Need2$Run == "EGAR00001653203")
    stopifnot(length(i_a) == 1, length(i_b) == 1)
    Need3 <- as.data.frame(Need2)
    Need3[i_a, cols] <- Need3[i_b, cols]
  }else{
    Need3 <- Need2
  }
  OS_all <- rbind(Need3, OS_all) |> 
    filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
}
Tumor_ICB_OS <- fusion_all |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_sum_OS) |> 
  select(Sample_ID, gene1, gene2,cohort, Score)|>
  rename(Run = Sample_ID) |>
  distinct()
cancer_type1 = "ALL"
cohort_ALL<- cohort1 |> filter(OS == TRUE )
cohort_select <- (cohort_ALL$id) |> unique()
fusion_select <- Tumor_ICB_OS |>
  filter(cohort %in% cohort_select) |> 
  as.data.table()
cohort_select <- fusion_select$cohort |> unique()

fusion_data_subset <- OS_all |>
  filter(cohort %in% cohort_select) |>
  mutate(event = case_when(
    Run %in% fusion_select$Run  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                           TRUE ~ event))  |>
  select(Patient_ID, event, OS_Time, OS_Status, cohort) |>
  distinct()

plot_os_by_group(data = fusion_data_subset, 
                 title = paste0("ICB-OS"), 
                 choose = "OS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4/Figure4A.pdf")
)
plot_os_by_group(data = fusion_data_subset, 
                 title = paste0("ICB-OS"),
                 choose = "OS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4/Figure4A.png")
)

#Figure4B
cohort_meta <- cohort1 |>
  distinct(cohort = id, cancer_type) |>
  group_by(cohort) |>
  summarise(cancer = paste(unique(na.omit(cancer_type)), collapse = "/"),
            .groups = "drop")

source_cox <- function(data, meta_tbl) {
  data |>
    mutate(event = factor(event, levels = c("Fusion-","Fusion+"))) |>
    group_by(cohort) |>
    group_modify(~{
      df <- .x
      n_fp <- sum(df$event == "Fusion+", na.rm = TRUE)
      n_fm <- sum(df$event == "Fusion-", na.rm = TRUE)
      if (n_fp < 2 || n_fm < 2 || length(unique(na.omit(df$OS_Status))) < 2)
        return(tibble(HR = NA_real_, lower = NA_real_, upper = NA_real_,
                      p = NA_real_, n = nrow(df),
                      n_event = sum(df$OS_Status, na.rm = TRUE)))
      s <- summary(coxph(Surv(OS_Time, OS_Status) ~ event, data = df))
      tibble(HR = s$conf.int["eventFusion+","exp(coef)"],
             lower = s$conf.int["eventFusion+","lower .95"],
             upper = s$conf.int["eventFusion+","upper .95"],
             p = s$coefficients["eventFusion+","Pr(>|z|)"],
             n = nrow(df), n_event = sum(df$OS_Status, na.rm = TRUE))
    }) |>
    ungroup() |>
    left_join(meta_tbl, by = "cohort") |>
    mutate(label = paste(cancer, cohort, sep = "_"))
}
cox_summary <- source_cox(fusion_data_subset, cohort_meta)
cohort_PFS <- cohort1 |> filter(PFS == TRUE) 
cohort_sum_PFS <- cohort_PFS$id |> unique()
PFS_all <- data.table()
for(cohort_id in cohort_sum_PFS){
  Need <- fread(paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/",cohort_id,"_info.tsv"))
  if(filter(cohort_PFS,id == cohort_id)$time == "Pre/On" ){
    if(cohort_id == "PRJNA356761"){
      Need1 <- Need |>
        filter(grepl("Pre", Sample_ID))
    }else if(cohort_id == "EGAD00001006282"){
      Need1 <- Need |>
        filter(grepl("_SCREEN", Sample_ID))
    }else if(cohort_id == "PHS001464"){
      Need2 <- Need2 |>
        mutate(Sample_type = case_when(
          TUMOR_COLLECTION_TIMEPOINT ==  "Acquired Resistance" ~ "On",
          TUMOR_COLLECTION_TIMEPOINT ==  "Pre-Immunotherapy"~ "Pre",
          TUMOR_COLLECTION_TIMEPOINT ==  "Off Therapy Recurrence"~ "On",
          TUMOR_COLLECTION_TIMEPOINT ==  "Diagnosis"~ "Pre",
          TUMOR_COLLECTION_TIMEPOINT ==  "Mixed Response (Site of Primary Resistance)"~ "On",
          TRUE ~ NA
        ) )
    }else{
      Need1 <- Need |>
        filter(grepl("Pre", Response2))
    }}else if(cohort_id %in% c("PRJNA482620", "PRJNA744780")){
      Need1 <- Need |>
        filter(grepl("Pre", Response2))
    }else{
      Need1 <- Need
    }
  columns_to_select <- c("Patient_ID", "PFS_Time", "PFS_Status", "Run", "Treatment", "Age", "Sex")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    dplyr::select(all_of(existing_columns)) |>
    filter(!is.na(PFS_Time) & !is.na(PFS_Status))
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  if (length(missing_columns) > 0) {
    for (col in missing_columns) {
      Need2[[col]] <- NA
    }
  }
  
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort <- cohort_id
  PFS_all <- rbind(Need2,PFS_all) |> filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
}
Tumor_ICB_PFS <- fusion_all |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_sum_PFS) |> 
  select(Sample_ID, gene1, gene2,cohort, Score)|>
  rename(Run = Sample_ID) |>
  distinct()
cancer_type1 = "ALL"
cohort_ALL<- cohort1 |> filter(PFS == TRUE )
cohort_select <- (cohort_ALL$id) |> unique()
fusion_select <- Tumor_ICB_PFS |>
  filter(cohort %in% cohort_select) |> 
  as.data.table()
cohort_select <- fusion_select$cohort |> unique()

#Figure4B
fusion_data_subset <- PFS_all |>
  filter(cohort %in% cohort_select) |>
  mutate(event = case_when(
    Run %in% fusion_select$Run  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                           TRUE ~ event))  |>
  select(Patient_ID, event, PFS_Time, PFS_Status, cohort)  |>
  distinct()
plot_os_by_group(data = fusion_data_subset, 
                 title = paste0("ICB-PFS"), 
                 choose = "PFS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4/Figure4B.pdf")
)
plot_os_by_group(data = fusion_data_subset, 
                 title = paste0("ICB-PFS"), 
                 choose = "PFS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4/Figure4B.png")
)


#Figure4C
#Figure4C
library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(forestploter)
library(grid)
library(stringr)
file_meta <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
metafusion <- fread(file_meta, sep = "\t")

Credibility_gene <- metafusion |>
  filter(sample_type %in% c("Tumor", "Tumor_Blood")) |>
  as.data.table()

# ---------------- 2. 读入临床，拼成一张表 ----------------
folder_path <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/"
clin_files  <- list.files(folder_path,
                          pattern = "^(CPTAC-|TARGET-|TCGA-).*_info\\.tsv$",
                          full.names = TRUE)

cohort_ids <- gsub("_", "-", sub("_info\\.tsv$", "", basename(clin_files)))
keep       <- !cohort_ids %in% c("CPTAC-NORMAL", "TARGET-NORMAL",
                                 "TCGA-NORMAL",  "TARGET-CELL")
clin_files <- clin_files[keep]
cohort_ids <- cohort_ids[keep]

all_clin <- rbindlist(lapply(seq_along(clin_files), function(i){
  fread(clin_files[i]) |>
    select(Patient_ID, Run, OS_Time, OS_Status) |>
    mutate(cohort = cohort_ids[i])
}), fill = TRUE) |>
  filter(!is.na(OS_Time), !is.na(OS_Status)) |>
  filter(!(grepl("^CPTAC-", cohort) & str_count(Run, "-") %in% c(5, 8, 11, 14))) |>
  distinct()

# 患者级 OS（一患者多样本取首条；如需其他规则改这里）
pat_os <- all_clin[, .(OS_Time   = OS_Time[1],
                       OS_Status = OS_Status[1]),
                   by = .(cohort, Patient_ID)]

# 样本 -> 患者映射
smp2pat <- all_clin[, .(Run, Patient_ID, cohort)] |> unique()
stopifnot(!anyDuplicated(smp2pat$Run))

# ---------------- 3. 患者层面 Fusion status ----------------
pat_fusion_pos <- Credibility_gene |>
  dplyr::distinct(Sample_ID) |>
  inner_join(smp2pat, by = c("Sample_ID" = "Run")) |>
  dplyr::distinct(cohort, Patient_ID) |>
  as.data.table()

stopifnot(!anyDuplicated(pat_fusion_pos, by = c("cohort", "Patient_ID")))

message("Fusion+ patients: ", nrow(pat_fusion_pos),
        " across ", uniqueN(pat_fusion_pos$cohort), " cohorts")

# ---------------- 4. 逐 cohort Cox ----------------
min_pos <- 3
min_neg <- 3

run_cohort <- function(ch){
  dat <- copy(pat_os[cohort == ch])
  if(nrow(dat) == 0) return(NULL)
  
  pos_ids <- pat_fusion_pos[cohort == ch, Patient_ID]
  dat[, event := ifelse(Patient_ID %in% pos_ids, "Fusion+", "Fusion-")]
  
  n_pos <- sum(dat$event == "Fusion+")
  n_neg <- sum(dat$event == "Fusion-")
  if(n_pos < min_pos || n_neg < min_neg) return(NULL)
  
  dat[, event := factor(event, levels = c("Fusion-", "Fusion+"))]
  
  fit <- tryCatch(
    suppressWarnings(coxph(Surv(OS_Time, OS_Status) ~ event, data = dat)),
    error = function(e) NULL)
  if(is.null(fit)) return(NULL)
  
  s <- tryCatch(summary(fit), error = function(e) NULL)
  if(is.null(s) || nrow(s$coefficients) < 1) return(NULL)
  
  data.table(cohort     = ch,
             n_pos      = n_pos,
             n_neg      = n_neg,
             HR         = s$coefficients[1, "exp(coef)"],
             low        = s$conf.int[1, "lower .95"],
             up         = s$conf.int[1, "upper .95"],
             p          = s$coefficients[1, "Pr(>|z|)"],
             n_strata   = 1L,
             is_summary = FALSE)
}

res_cohort <- rbindlist(lapply(unique(pat_os$cohort), run_cohort), fill = TRUE)
if(nrow(res_cohort) == 0)
  stop("No cohort passed filters — check Fusion/patient mapping and n thresholds.")

# ---------------- 5. TCGA / TARGET / CPTAC 分层 Cox 汇总 ----------------
run_strata_summary <- function(cohort_set, label,
                               min_pos = 3, min_neg = 3, min_strata = 2){
  if(length(cohort_set) < min_strata) return(NULL)
  dat <- copy(pat_os[cohort %in% cohort_set])
  if(nrow(dat) == 0) return(NULL)
  
  pos_ids <- pat_fusion_pos[cohort %in% cohort_set]
  key_pos <- paste(pos_ids$cohort, pos_ids$Patient_ID)
  dat[, event := ifelse(paste(cohort, Patient_ID) %in% key_pos,
                        "Fusion+", "Fusion-")]
  
  # 只保留两类都有的 cohort（否则该层对 event 系数无信息）
  dat <- dat[, if(uniqueN(event) == 2) .SD, by = cohort]
  
  n_pos    <- sum(dat$event == "Fusion+")
  n_neg    <- sum(dat$event == "Fusion-")
  n_strata <- uniqueN(dat$cohort)
  if(n_pos < min_pos || n_neg < min_neg || n_strata < min_strata) return(NULL)
  
  dat[, event := factor(event, levels = c("Fusion-", "Fusion+"))]
  fit <- tryCatch(
    suppressWarnings(coxph(Surv(OS_Time, OS_Status) ~ event + strata(cohort),
                           data = dat)),
    error = function(e) NULL)
  if(is.null(fit)) return(NULL)
  
  s <- tryCatch(summary(fit), error = function(e) NULL)
  if(is.null(s) || nrow(s$coefficients) < 1) return(NULL)
  
  data.table(cohort     = label,
             n_pos      = n_pos,
             n_neg      = n_neg,
             HR         = s$coefficients[1, "exp(coef)"],
             low        = s$conf.int[1, "lower .95"],
             up         = s$conf.int[1, "upper .95"],
             p          = s$coefficients[1, "Pr(>|z|)"],
             n_strata   = n_strata,
             is_summary = TRUE)
}

cohorts_tcga   <- unique(pat_os$cohort[grepl("^TCGA-",   pat_os$cohort)])
cohorts_target <- unique(pat_os$cohort[grepl("^TARGET-", pat_os$cohort)])
cohorts_cptac  <- unique(pat_os$cohort[grepl("^CPTAC-",  pat_os$cohort)])

res_summary <- rbindlist(list(
  run_strata_summary(cohorts_tcga,   "TCGA-ALL (stratified)"),
  run_strata_summary(cohorts_target, "TARGET-ALL (stratified)"),
  run_strata_summary(cohorts_cptac,  "CPTAC-ALL (stratified)")
), fill = TRUE)

# ---------------- 6. 合并、拆分、FDR (分家族 BH) ----------------
res <- rbindlist(list(res_cohort, res_summary), fill = TRUE)
res[, is_summary := fifelse(is.na(is_summary), FALSE, is_summary)]
res[, finite := is.finite(HR) & is.finite(low) & is.finite(up) & !is.na(p)]

res_note <- res[finite == FALSE]
res_ok   <- res[finite == TRUE]

# 4 个家族：TCGA 单队列 / TARGET 单队列 / CPTAC 单队列 / SUMMARY
res_ok[, source := fifelse(is_summary, "SUMMARY",
                           fifelse(grepl("^TCGA",   cohort), "TCGA",
                                   fifelse(grepl("^TARGET", cohort), "TARGET",
                                           fifelse(grepl("^CPTAC",  cohort), "CPTAC", NA_character_))))]

stopifnot(!anyNA(res_ok$source))   # 若有行 source 为 NA，说明 cohort 命名不符预期

res_ok[, q   := p.adjust(p, method = "BH"), by = source]
res_ok[, sig := q < 0.05]
res_ok[, sig_nominal := p < 0.05]

# ---------------- 9. 显著 / 不显著 ----------------
summary_rows  <- res_ok[is_summary == TRUE]
cohort_sig    <- res_ok[is_summary == FALSE & sig == TRUE]
cohort_nonsig <- res_ok[is_summary == FALSE & sig == FALSE]

message("Cohorts sig(q<0.05): ",
        res_ok[is_summary == FALSE & sig == TRUE, .N],
        " | summary sig: ",
        res_ok[is_summary == TRUE & sig == TRUE, .N],
        " | summary rows present: ",
        res_ok[is_summary == TRUE, .N],
        " | Inf/non-estimable: ", nrow(res_note))
print(res_ok[is_summary == TRUE, .(cohort, n_pos, n_neg, HR, low, up, p, q, sig)])

# ---------------- 8. make_forest：footnote 换行 + 高度 ----------------
tm <- forest_theme(
  base_size = 14,
  colhead  = list(fg_params = list(hjust = 0)),
  core     = list(padding = unit(c(8, 10), "mm")),
  ci_pch = 15, ci_col = "black", ci_fill = "black",
  ci_alpha = 0.8, ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.2,
  refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey20"),
  footnote_gp = gpar(cex = 0.9, fontface = "bold", col = "black"),
  vertline_lwd = 1, vertline_lty = "dashed", vertline_col = "grey20"
)
make_forest <- function(dt, title, out_pdf, max_rows = 80,
                        AXIS_MAX = 20, x_trans = "log") {
  
  if (nrow(dt) == 0) { message("empty: ", title); return(invisible(NULL)) }
  if (!"is_summary" %in% names(dt)) dt[, is_summary := FALSE]
  
  dt <- copy(dt)
  dt[, row_order := fifelse(is_summary, 0L, 1L)]
  setorder(dt, row_order)
  dt <- dt[1:min(.N, max_rows)]
  
  x_lo <- 1 / AXIS_MAX
  x_hi <- AXIS_MAX
  
  dt[, `:=`(
    low_p = pmin(pmax(low, x_lo), x_hi),
    up_p  = pmax(pmin(up,  x_hi), x_lo),
    HR_p  = pmin(pmax(HR,  x_lo), x_hi),
    cens  = (!is.na(low) & low < x_lo) | (!is.na(up) & up > x_hi)
  )]
  dt[low_p > up_p, `:=`(low_p = up_p, up_p = low_p)]
  
  fmt_p <- function(x, na = "") {
    out <- ifelse(is.na(x), na, sprintf("%.3f", x))
    out[which(!is.na(x) & x < 0.001)] <- "<0.001"
    out
  }
  
  dt[, `Fusion+/Fusion-` := paste0(n_pos, "/", n_neg)]
  dt[, `HR (95% CI)` := paste0(
    sprintf("%.2f (%.2f-%.2f)", HR, low, up),
    fifelse(cens, "\u2020", "")
  )]
  dt[, p_lab := fmt_p(p)]
  dt[, q_lab := fmt_p(q)]
  dt[, blank := paste(rep(" ", 30), collapse = " ")]
  
  plot_dt <- as.data.frame(dt[, .(cohort, `Fusion+/Fusion-`, blank,
                                  `HR (95% CI)`, p_lab, q_lab)])
  colnames(plot_dt) <- c("Cohort", "Binary fusion status", " ", "HR (95% CI)", "p", "q (BH)")
  
  note_txt <- if (exists("res_note") && nrow(res_note) > 0) {
    paste0("Not shown (HR = Inf / non-estimable): ",
           paste(res_note$cohort, collapse = "; "))
  } else ""
  
  footnote_txt <- paste(
    "\n", "\n", "\n", "\n",
    "GDC OS: Fusion+ vs Fusion-",
    "Stratified summary (~ Fusion + strata(cohort)).",
    "q: BH-FDR.",
    note_txt, sep = "\n"
  )
  
  wrap_lines <- function(txt, width = 110) {
    lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
    lines <- unlist(lapply(lines, function(x) if (nchar(x) == 0) "" else strwrap(x, width = width)))
    paste(lines, collapse = "\n")
  }
  footnote_txt <- wrap_lines(footnote_txt, width = 110)
  n_fn_lines <- stringr::str_count(footnote_txt, "\n") + 1
  h <- max(6, nrow(plot_dt) * 0.45 + 3 + n_fn_lines * 0.22)
  
  p <- forest(plot_dt,
              est = dt$HR_p, lower = dt$low_p, upper = dt$up_p,
              ci_column = 3, ref_line = 1,
              xlim = c(x_lo, x_hi), x_trans = x_trans,
              ticks_at = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20),
              ci_pch = ifelse(dt$is_summary, 18, 15),
              ci_lwd = ifelse(dt$is_summary, 2.0, 1.5),
              title = title, footnote = footnote_txt, theme = tm)
  
  p <- insert_text(p, text = "Fusion+/Fusion-", row = 1, col = 2, just = "left",
                   gp = gpar(cex = 0.7, col = "black", fontface = "italic"))
  
  n_sum <- sum(dt$is_summary)
  if (n_sum > 0 && n_sum < nrow(dt)) {
    p <- add_border(p, part = "body", row = n_sum + 1,
                    gp = gpar(lwd = 0.8, col = "grey40"))
  }
  
  ## ---- 尺寸：按内容实测，别写死 10 in（写死会把左右两列裁掉）--------
  wh <- forestploter::get_wh(p, unit = "in")
  w  <- wh[1] * 1.02
  h  <- wh[2] * 1.02 + n_fn_lines * 0.175 + 0.15   # get_wh 不含脚注行高
  
  ext <- tolower(tools::file_ext(out_pdf))
  if (ext == "png") {
    png(out_pdf, width = w, height = h, units = "in", res = 300)
  } else {
    pdf(out_pdf, width = w, height = h, onefile = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  
  print(p)
  invisible(out_pdf)
}
make_forest(rbind(summary_rows, cohort_sig, fill = TRUE)[order(is_summary, q)],
            "GDC OS pooled estimate (individual cohorts shown; q < 0.05)",
            file.path(out_dir, "Figure4C.pdf"))


#Figure4d
cohort_ids <- gsub("_", "-", sub("_info\\.tsv$", "", basename(clin_files)))
keep       <- !cohort_ids %in% c("CPTAC-NORMAL", "TARGET-NORMAL",
                                 "TCGA-NORMAL",  "TARGET-CELL",
                                 "TARGET-NBL", "TARGET-GNB",
                                 "TARGET-LAML", "TCGA-LAML")
clin_files <- clin_files[keep]
cohort_ids <- cohort_ids[keep]

read_one_pfs <- function(f, ch){
  dt <- fread(f)
  data.table(
    Patient_ID = dt$Patient_ID,
    Run        = dt$Run,
    PFS_Time   = if("PFS_Time"   %in% names(dt)) as.numeric(dt$PFS_Time)   else NA_real_,
    PFS_Status = if("PFS_Status" %in% names(dt)) as.integer(dt$PFS_Status) else NA_integer_,
    cohort     = ch
  )
}

all_clin_pfs <- rbindlist(lapply(seq_along(clin_files), function(i){
  read_one_pfs(clin_files[i], cohort_ids[i])
}), fill = TRUE) |>
  filter(!is.na(PFS_Time), !is.na(PFS_Status)) |>
  filter(!(grepl("^CPTAC-", cohort) & str_count(Run, "-") %in% c(5, 8, 11, 14))) |>
  distinct() 

message("PFS usable cohorts: ", uniqueN(all_clin_pfs$cohort),
        " | patients: ", nrow(all_clin_pfs))

# 患者级 PFS（一患者多样本取首条）
pat_pfs <- all_clin_pfs[, .(PFS_Time   = PFS_Time[1],
                            PFS_Status = PFS_Status[1]),
                        by = .(cohort, Patient_ID)]

# 样本 -> 患者映射（PFS 版，与 OS 可能不同队列集）
smp2pat_pfs <- all_clin_pfs[, .(Run, Patient_ID, cohort)] |> unique()
stopifnot(!anyDuplicated(smp2pat_pfs$Run))

# ---------------- 3. 患者层面 Fusion status（复用同一 meta） ----------------
pat_fusion_pos <- Credibility_gene |>
  dplyr::distinct(Sample_ID) |>
  inner_join(smp2pat_pfs, by = c("Sample_ID" = "Run")) |>
  dplyr::distinct(cohort, Patient_ID) |>
  as.data.table()

stopifnot(!anyDuplicated(pat_fusion_pos, by = c("cohort", "Patient_ID")))

message("Fusion+ patients with PFS: ", nrow(pat_fusion_pos),
        " across ", uniqueN(pat_fusion_pos$cohort), " cohorts")

# ---------------- 4. 逐 cohort Cox (PFS) ----------------
min_pos <- 3
min_neg <- 3

run_cohort_pfs <- function(ch){
  dat <- copy(pat_pfs[cohort == ch])
  if(nrow(dat) == 0) return(NULL)
  
  pos_ids <- pat_fusion_pos[cohort == ch, Patient_ID]
  dat[, event := ifelse(Patient_ID %in% pos_ids, "Fusion+", "Fusion-")]
  
  n_pos <- sum(dat$event == "Fusion+")
  n_neg <- sum(dat$event == "Fusion-")
  if(n_pos < min_pos || n_neg < min_neg) return(NULL)
  
  dat[, event := factor(event, levels = c("Fusion-", "Fusion+"))]
  fit <- tryCatch(
    suppressWarnings(coxph(Surv(PFS_Time, PFS_Status) ~ event, data = dat)),
    error = function(e) NULL)
  if(is.null(fit)) return(NULL)
  
  s <- tryCatch(summary(fit), error = function(e) NULL)
  if(is.null(s) || nrow(s$coefficients) < 1) return(NULL)
  
  data.table(cohort     = ch,
             n_pos      = n_pos,
             n_neg      = n_neg,
             HR         = s$coefficients[1, "exp(coef)"],
             low        = s$conf.int[1, "lower .95"],
             up         = s$conf.int[1, "upper .95"],
             p          = s$coefficients[1, "Pr(>|z|)"],
             n_strata   = 1L,
             is_summary = FALSE)
}

res_cohort_pfs <- rbindlist(lapply(unique(pat_pfs$cohort), run_cohort_pfs), fill = TRUE)
if(nrow(res_cohort_pfs) == 0)
  stop("No cohort passed PFS filters — check PFS columns and n thresholds.")

# ---------------- 5. TCGA / TARGET / CPTAC 分层 Cox 汇总 (PFS) ----------------
run_strata_summary_pfs <- function(cohort_set, label,
                                   min_pos = 3, min_neg = 3, min_strata = 2){
  if(length(cohort_set) < min_strata) return(NULL)
  dat <- copy(pat_pfs[cohort %in% cohort_set])
  if(nrow(dat) == 0) return(NULL)
  
  pos_ids <- pat_fusion_pos[cohort %in% cohort_set]
  key_pos <- paste(pos_ids$cohort, pos_ids$Patient_ID)
  dat[, event := ifelse(paste(cohort, Patient_ID) %in% key_pos,
                        "Fusion+", "Fusion-")]
  
  dat <- dat[, if(uniqueN(event) == 2) .SD, by = cohort]
  
  n_pos    <- sum(dat$event == "Fusion+")
  n_neg    <- sum(dat$event == "Fusion-")
  n_strata <- uniqueN(dat$cohort)
  if(n_pos < min_pos || n_neg < min_neg || n_strata < min_strata) return(NULL)
  
  dat[, event := factor(event, levels = c("Fusion-", "Fusion+"))]
  fit <- tryCatch(
    suppressWarnings(coxph(Surv(PFS_Time, PFS_Status) ~ event + strata(cohort),
                           data = dat)),
    error = function(e) NULL)
  if(is.null(fit)) return(NULL)
  
  s <- tryCatch(summary(fit), error = function(e) NULL)
  if(is.null(s) || nrow(s$coefficients) < 1) return(NULL)
  
  data.table(cohort     = label,
             n_pos      = n_pos,
             n_neg      = n_neg,
             HR         = s$coefficients[1, "exp(coef)"],
             low        = s$conf.int[1, "lower .95"],
             up         = s$conf.int[1, "upper .95"],
             p          = s$coefficients[1, "Pr(>|z|)"],
             n_strata   = n_strata,
             is_summary = TRUE)
}

cohorts_tcga   <- unique(pat_pfs$cohort[grepl("^TCGA-",   pat_pfs$cohort)])
cohorts_target <- unique(pat_pfs$cohort[grepl("^TARGET-", pat_pfs$cohort)])
cohorts_cptac  <- unique(pat_pfs$cohort[grepl("^CPTAC-",  pat_pfs$cohort)])

res_summary_pfs <- rbindlist(list(
  run_strata_summary_pfs(cohorts_tcga,   "TCGA-ALL (stratified)"),
  run_strata_summary_pfs(cohorts_target, "TARGET-ALL (stratified)"),
  run_strata_summary_pfs(cohorts_cptac,  "CPTAC-ALL (stratified)")
), fill = TRUE)

# ---------------- 6. 合并、拆分、4 家族 BH-FDR ----------------
res <- rbindlist(list(res_cohort_pfs, res_summary_pfs), fill = TRUE)
res[, is_summary := fifelse(is.na(is_summary), FALSE, is_summary)]
res[, finite := is.finite(HR) & is.finite(low) & is.finite(up) & !is.na(p)]

res_note <- res[finite == FALSE]
res_ok   <- res[finite == TRUE]

res_ok[, source := fifelse(is_summary, "SUMMARY",
                           fifelse(grepl("^TCGA",   cohort), "TCGA",
                                   fifelse(grepl("^TARGET", cohort), "TARGET",
                                           fifelse(grepl("^CPTAC",  cohort), "CPTAC", NA_character_))))]
stopifnot(!anyNA(res_ok$source))

res_ok[, q   := p.adjust(p, method = "BH"), by = source]
res_ok[, sig := q < 0.05]
res_ok[, sig_nominal := p < 0.05]

message("PFS cohorts sig(q<0.05): ",
        res_ok[is_summary == FALSE & sig == TRUE, .N],
        " | summary sig: ",
        res_ok[is_summary == TRUE & sig == TRUE, .N],
        " | summary rows present: ",
        res_ok[is_summary == TRUE, .N],
        " | Inf/non-estimable: ", nrow(res_note))
print(res_ok[is_summary == TRUE, .(cohort, n_pos, n_neg, HR, low, up, p, q, sig)])

# ---------------- 9. 显著 / 不显著两张 PFS 图 ----------------
summary_rows  <- res_ok[is_summary == TRUE]
cohort_sig    <- res_ok[is_summary == FALSE & sig == TRUE]
cohort_nonsig <- res_ok[is_summary == FALSE & sig == FALSE]

# ---------------- 8. make_forest：footnote 换行 + 高度 ----------------
tm <- forest_theme(
  base_size = 14,
  colhead  = list(fg_params = list(hjust = 0)),
  core     = list(padding = unit(c(8, 10), "mm")),
  ci_pch = 15, ci_col = "black", ci_fill = "black",
  ci_alpha = 0.8, ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.2,
  refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey20"),
  footnote_gp = gpar(cex = 0.9, fontface = "bold", col = "black"),
  vertline_lwd = 1, vertline_lty = "dashed", vertline_col = "grey20"
)
make_forest <- function(dt, title, out_pdf, max_rows = 80,
                        AXIS_MAX = 20, x_trans = "log") {
  
  if (nrow(dt) == 0) { message("empty: ", title); return(invisible(NULL)) }
  if (!"is_summary" %in% names(dt)) dt[, is_summary := FALSE]
  
  dt <- copy(dt)
  dt[, row_order := fifelse(is_summary, 0L, 1L)]
  setorder(dt, row_order)
  dt <- dt[1:min(.N, max_rows)]
  
  x_lo <- 1 / AXIS_MAX
  x_hi <- AXIS_MAX
  
  dt[, `:=`(
    low_p = pmin(pmax(low, x_lo), x_hi),
    up_p  = pmax(pmin(up,  x_hi), x_lo),
    HR_p  = pmin(pmax(HR,  x_lo), x_hi),
    cens  = (!is.na(low) & low < x_lo) | (!is.na(up) & up > x_hi)
  )]
  dt[low_p > up_p, `:=`(low_p = up_p, up_p = low_p)]
  
  fmt_p <- function(x, na = "") {
    out <- ifelse(is.na(x), na, sprintf("%.3f", x))
    out[which(!is.na(x) & x < 0.001)] <- "<0.001"
    out
  }
  
  dt[, `Fusion+/Fusion-` := paste0(n_pos, "/", n_neg)]
  dt[, `HR (95% CI)` := paste0(
    sprintf("%.2f (%.2f-%.2f)", HR, low, up),
    fifelse(cens, "\u2020", "")
  )]
  dt[, p_lab := fmt_p(p)]
  dt[, q_lab := fmt_p(q)]
  dt[, blank := paste(rep(" ", 30), collapse = " ")]
  
  plot_dt <- as.data.frame(dt[, .(cohort, `Fusion+/Fusion-`, blank,
                                  `HR (95% CI)`, p_lab, q_lab)])
  colnames(plot_dt) <- c("Cohort", "Binary Fusion Status", " ", "HR (95% CI)", "p", "q (BH)")
  
  note_txt <- if (exists("res_note") && nrow(res_note) > 0) {
    paste0("Not shown (HR = Inf / non-estimable): ",
           paste(res_note$cohort, collapse = "; "))
  } else ""
  
  footnote_txt <- paste(
    "\n", "\n", "\n", "\n",
    "GDC PFS: Fusion+ vs Fusion-.",
    "Stratified summary (~ Fusion + strata(cohort)).",
    "q: BH-FDR.",
    note_txt, sep = "\n"
  )
  
  wrap_lines <- function(txt, width = 110) {
    lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
    lines <- unlist(lapply(lines, function(x) if (nchar(x) == 0) "" else strwrap(x, width = width)))
    paste(lines, collapse = "\n")
  }
  footnote_txt <- wrap_lines(footnote_txt, width = 110)
  n_fn_lines <- stringr::str_count(footnote_txt, "\n") + 1
  h <- max(6, nrow(plot_dt) * 0.45 + 3 + n_fn_lines * 0.22)
  
  p <- forest(plot_dt,
              est = dt$HR_p, lower = dt$low_p, upper = dt$up_p,
              ci_column = 3, ref_line = 1,
              xlim = c(x_lo, x_hi), x_trans = x_trans,
              ticks_at = c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10),
              ci_pch = ifelse(dt$is_summary, 18, 15),
              ci_lwd = ifelse(dt$is_summary, 2.0, 1.5),
              title = title, footnote = footnote_txt, theme = tm)
  
  p <- insert_text(p, text = "Fusion+/Fusion-", row = 1, col = 2, just = "left",
                   gp = gpar(cex = 0.7, col = "black", fontface = "italic"))
  
  n_sum <- sum(dt$is_summary)
  if (n_sum > 0 && n_sum < nrow(dt)) {
    p <- add_border(p, part = "body", row = n_sum + 1,
                    gp = gpar(lwd = 0.8, col = "grey40"))
  }
  
  ## ---- 尺寸：按内容实测，别写死 10 in（写死会把左右两列裁掉）--------
  wh <- forestploter::get_wh(p, unit = "in")
  w  <- wh[1] * 1.02
  h  <- wh[2] * 1.02 + n_fn_lines * 0.175 + 0.15   # get_wh 不含脚注行高
  
  ext <- tolower(tools::file_ext(out_pdf))
  if (ext == "png") {
    png(out_pdf, width = w, height = h, units = "in", res = 300)
  } else {
    pdf(out_pdf, width = w, height = h, onefile = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  
  print(p)
  invisible(out_pdf)
}
make_forest(rbind(summary_rows, cohort_sig, fill = TRUE)[order(is_summary, q)],
            "GDC PFS pooled estimate (individual cohorts shown; q < 0.05)",
            file.path(out_dir, "Figure4D.pdf"))
message("Done. PFS outputs in: ", out_dir)


#Figure4E
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(sandwich); library(lmtest)
})

out_dir  <- Sys.getenv("FIG4E_OUT",
                       unset = "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
data_dir  <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据"
folder_path <- file.path(data_dir, "Clininfo/")

## ---- 1. 融合分类（原来靠 figure4EF.R 里的 fusion_all，这里内置）-----
fusion_all <- fread(file.path(data_dir, "metafusion_calls.star.tsv")) |>
  mutate(Fusion_Class = case_when(
    inferred_fusion_type %in% c("CodingFusion") ~ "Coding",
    inferred_fusion_type %in% c("TruncatedCoding", "TruncatedNoncoding") ~ "Truncated",
    inferred_fusion_type %in% c("ReadThrough", "SameGene", "Intergenic", "NoHeadGene") ~ "Passenger",
    TRUE ~ ""))

file_names <- list.files(folder_path, full.names = FALSE,
                         pattern = "^(CPTAC-|TARGET|TCGA-).*_info\\.tsv$")
cohort_sum      <- sub("_info.tsv$", "", file_names) |> setdiff("TCGA-NORMAL")
cohort_sum_tcga <- cohort_sum[grepl("^TCGA", cohort_sum)]

## ---- 2. 临床：ploidy / purity / WGD --------------------------------
TCGA_all <- rbindlist(lapply(cohort_sum_tcga, function(cohort_id) {
  fread(file.path(folder_path, paste0(gsub("_", "-", cohort_id), "_info.tsv"))) |>
    select(ploidy, purity, Genome_doublings, Sample_ID) |>
    mutate(Sample_ID = as.character(Sample_ID), cohort = cohort_id,
           ploidy = as.numeric(ploidy), purity = as.numeric(purity),
           Genome_doublings = as.numeric(Genome_doublings)) |>
    group_by(Sample_ID, cohort) |>
    summarise(ploidy = if (all(is.na(ploidy))) NA_real_ else max(ploidy, na.rm = TRUE),
              purity = if (all(is.na(purity))) NA_real_ else max(purity, na.rm = TRUE),
              Genome_doublings = if (all(is.na(Genome_doublings))) NA_real_ else max(Genome_doublings, na.rm = TRUE),
              .groups = "drop") |>
    filter(!is.na(purity) & !is.na(ploidy))
}))

## ---- 3. TMB / HRD --------------------------------------------------
mc3 <- fread(file.path(data_dir, "mc3.v0.2.8.PUBLIC.xena"))
silent <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
            "IGR", "RNA", "Splice_Region")
tmb_all <- mc3 |>
  filter(!effect %in% silent) |>
  rename(Sample_ID = sample) |>
  count(Sample_ID, name = "n_mutations") |>
  mutate(TMB = n_mutations / 38, Sample_ID = substr(Sample_ID, 1, 15)) |>
  group_by(Sample_ID) |> slice_max(TMB, n = 1, with_ties = FALSE) |> ungroup()

hrd_raw <- fread(file.path(data_dir, "TCGA.HRD_withSampleID.txt"))
loh_df  <- as.data.frame(t(hrd_raw[, -1]), stringsAsFactors = FALSE)
colnames(loh_df) <- hrd_raw[[1]]
loh_df$Sample_ID <- substr(rownames(loh_df), 1, 15)
loh_df <- loh_df |>
  transmute(Sample_ID, HRD = as.numeric(HRD)) |>
  group_by(Sample_ID) |> slice_max(HRD, n = 1, with_ties = FALSE) |> ungroup()

## ---- 4. 合并 + 分组 ------------------------------------------------
fusion_select <- fusion_all |> filter(cohort %in% cohort_sum_tcga)

plot_df <- TCGA_all |>
  left_join(tmb_all, by = "Sample_ID") |>
  left_join(loh_df, by = "Sample_ID") |>
  filter(!is.na(TMB) & !is.na(HRD)) |>
  left_join(fusion_select |> select(Sample_ID, Fusion_Class), by = "Sample_ID") |>
  group_by(Sample_ID) |>
  mutate(event = case_when(all(is.na(Fusion_Class)) ~ "Fusion-", TRUE ~ "Fusion+")) |>
  select(Sample_ID, event, HRD, TMB, ploidy, purity, cohort) |>
  distinct() |>
  ungroup() |>
  mutate(event = factor(event, levels = c("Fusion-", "Fusion+")))

cat("分析集 n =", nrow(plot_df), "  Fusion+ =", sum(plot_df$event == "Fusion+"),
    "  Fusion- =", sum(plot_df$event == "Fusion-"), "  癌种 =", uniqueN(plot_df$cohort), "\n")

## ---- 5. 队列固定效应 + HC3 ------------------------------------------
## scale = "auto_log": 全为正 -> 取 log（比值可直接解释成倍数）；
##                       有 0 -> 用 log1p（比值解释成 (值+1) 的倍数）
## scale = "identity": 原尺度（有界变量，如 purity）
adj_effect <- function(d, y, scale = "auto_log") {
  ok <- !is.na(d[[y]]) & !is.na(d$event) & !is.na(d$cohort)
  v  <- d[[y]][ok]
  if (scale == "auto_log") {
    if (all(v > 0)) { v2 <- log(v);   kind <- "ratio"; used <- "log" }
    else            { v2 <- log1p(v); kind <- "ratio"; used <- "log1p" }
  } else { v2 <- v; kind <- "diff"; used <- "identity" }
  
  dd <- data.frame(y = v2, fusion = d$event[ok], cohort = d$cohort[ok])
  m  <- lm(y ~ fusion + cohort, data = dd)
  ct <- coeftest(m, vcov. = vcovHC(m, type = "HC3"))
  b  <- ct["fusionFusion+", "Estimate"]
  se <- ct["fusionFusion+", "Std. Error"]
  p  <- ct["fusionFusion+", "Pr(>|t|)"]
  
  if (kind == "ratio")
    list(kind = kind, scale_used = used, est = exp(b), lo = exp(b - 1.96 * se),
         hi = exp(b + 1.96 * se), p = p, n = sum(ok))
  else
    list(kind = kind, scale_used = used, est = b, lo = b - 1.96 * se,
         hi = b + 1.96 * se, p = p, n = sum(ok))
}

METRICS <- data.table(
  var   = c("TMB", "HRD", "ploidy", "purity"),
  label = c("TMB (mut/Mb)", "HRD Score", "Ploidy", "Tumor Purity"),
  scale = c("auto_log", "auto_log", "auto_log", "identity"))

res <- rbindlist(lapply(seq_len(nrow(METRICS)), function(i) {
  a <- adj_effect(plot_df, METRICS$var[i], METRICS$scale[i])
  med <- tapply(plot_df[[METRICS$var[i]]], plot_df$event, median, na.rm = TRUE)
  data.table(var = METRICS$var[i], label = METRICS$label[i], kind = a$kind,
             scale_used = a$scale_used, n = a$n,
             median_fp = med[["Fusion+"]], median_fn = med[["Fusion-"]],
             est = a$est, lo = a$lo, hi = a$hi, p = a$p)
}))

## 图上/表里的那句话
res[, effect_lab := ifelse(kind == "ratio",
                           "Adj. GMR (95% CI)",
                           "Adj. mean difference (95% CI)")]
res[, value_lab := ifelse(kind == "ratio",
                          sprintf("%.2f (%.2f-%.2f)", est, lo, hi),
                          sprintf("%+.3f (%.3f to %+.3f)", est, lo, hi))]

res[, p_txt := ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p))]

res[, strip := paste0(label, "\n",
                      effect_lab, "\n",     # ← 这里多一个 \n 就是空行
                      value_lab, "\n",
                      p_txt)]

cat("\n== 队列固定效应校正后的效应量 ==\n")
print(res[, .(label, scale_used, n, median_fp = round(median_fp, 2), median_fn = round(median_fn, 2),
              est = signif(est, 3), lo = signif(lo, 3), hi = signif(hi, 3), p = signif(p, 3))])

## ---- 6. 图 ---------------------------------------------------------
plot_long <- plot_df |>
  pivot_longer(cols = METRICS$var, names_to = "var", values_to = "value") |>
  filter(!is.na(value)) |>
  left_join(res[, .(var, strip)], by = "var") |>
  mutate(strip = factor(strip, levels = res$strip))

pal <- c("Fusion-" = "#00468bff", "Fusion+" = "#ed0000ff")

p <- ggplot(plot_long, aes(x = event, y = value, fill = event)) +
  geom_boxplot(width = 0.6, outlier.size = 0.4, alpha = 0.85, colour = "grey20") +
  facet_wrap(~ strip, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = NULL,
       title = "TCGA Pan-cancer: Fusion status vs genomic features") +
  theme_bw(base_size = 15) +
  theme(legend.position  = "none",
        strip.text       = element_text(size = 14, face = "bold",  lineheight = 0.95),
        axis.text.x      = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 14),
        plot.title       = element_text(size = 25, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())
p
ggsave(file.path(out_dir, "Figure4E.pdf"), p, width = 15, height = 7, dpi = 600)
ggsave(file.path(out_dir, "Figure4E.png"), p, width = 15, height = 7, dpi = 600)


#Figure4E
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(sandwich); library(lmtest)
})

out_dir  <- Sys.getenv("FIG4E_OUT",
                       unset = "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
data_dir  <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据"
folder_path <- file.path(data_dir, "Clininfo/")

## ---- 1. 融合分类（原来靠 figure4EF.R 里的 fusion_all，这里内置）-----
fusion_all <- fread(file.path(data_dir, "metafusion_calls.star.tsv")) |>
  mutate(Fusion_Class = case_when(
    inferred_fusion_type %in% c("CodingFusion") ~ "Coding",
    inferred_fusion_type %in% c("TruncatedCoding", "TruncatedNoncoding") ~ "Truncated",
    inferred_fusion_type %in% c("ReadThrough", "SameGene", "Intergenic", "NoHeadGene") ~ "Passenger",
    TRUE ~ ""))

file_names <- list.files(folder_path, full.names = FALSE,
                         pattern = "^(CPTAC-|TARGET|TCGA-).*_info\\.tsv$")
cohort_sum      <- sub("_info.tsv$", "", file_names) |> setdiff("TCGA-NORMAL")
cohort_sum_tcga <- cohort_sum[grepl("^TCGA", cohort_sum)]

## ---- 2. 临床：ploidy / purity / WGD --------------------------------
TCGA_all <- rbindlist(lapply(cohort_sum_tcga, function(cohort_id) {
  fread(file.path(folder_path, paste0(gsub("_", "-", cohort_id), "_info.tsv"))) |>
    select(ploidy, purity, Genome_doublings, Sample_ID) |>
    mutate(Sample_ID = as.character(Sample_ID), cohort = cohort_id,
           ploidy = as.numeric(ploidy), purity = as.numeric(purity),
           Genome_doublings = as.numeric(Genome_doublings)) |>
    group_by(Sample_ID, cohort) |>
    summarise(ploidy = if (all(is.na(ploidy))) NA_real_ else max(ploidy, na.rm = TRUE),
              purity = if (all(is.na(purity))) NA_real_ else max(purity, na.rm = TRUE),
              Genome_doublings = if (all(is.na(Genome_doublings))) NA_real_ else max(Genome_doublings, na.rm = TRUE),
              .groups = "drop") |>
    filter(!is.na(purity) & !is.na(ploidy))
}))

## ---- 3. TMB / HRD --------------------------------------------------
mc3 <- fread(file.path(data_dir, "mc3.v0.2.8.PUBLIC.xena"))
silent <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
            "IGR", "RNA", "Splice_Region")
tmb_all <- mc3 |>
  filter(!effect %in% silent) |>
  rename(Sample_ID = sample) |>
  count(Sample_ID, name = "n_mutations") |>
  mutate(TMB = n_mutations / 38, Sample_ID = substr(Sample_ID, 1, 15)) |>
  group_by(Sample_ID) |> slice_max(TMB, n = 1, with_ties = FALSE) |> ungroup()

hrd_raw <- fread(file.path(data_dir, "TCGA.HRD_withSampleID.txt"))
loh_df  <- as.data.frame(t(hrd_raw[, -1]), stringsAsFactors = FALSE)
colnames(loh_df) <- hrd_raw[[1]]
loh_df$Sample_ID <- substr(rownames(loh_df), 1, 15)
loh_df <- loh_df |>
  transmute(Sample_ID, HRD = as.numeric(HRD)) |>
  group_by(Sample_ID) |> slice_max(HRD, n = 1, with_ties = FALSE) |> ungroup()

## ---- 4. 合并 + 分组 ------------------------------------------------
fusion_select <- fusion_all |> filter(cohort %in% cohort_sum_tcga)

plot_df <- TCGA_all |>
  left_join(tmb_all, by = "Sample_ID") |>
  left_join(loh_df, by = "Sample_ID") |>
  filter(!is.na(TMB) & !is.na(HRD)) |>
  left_join(fusion_select |> select(Sample_ID, Fusion_Class), by = "Sample_ID") |>
  group_by(Sample_ID) |>
  mutate(event = case_when(all(is.na(Fusion_Class)) ~ "Fusion-", TRUE ~ "Fusion+")) |>
  select(Sample_ID, event, HRD, TMB, ploidy, purity, cohort) |>
  distinct() |>
  ungroup() |>
  mutate(event = factor(event, levels = c("Fusion-", "Fusion+")))

cat("分析集 n =", nrow(plot_df), "  Fusion+ =", sum(plot_df$event == "Fusion+"),
    "  Fusion- =", sum(plot_df$event == "Fusion-"), "  癌种 =", uniqueN(plot_df$cohort), "\n")

## ---- 5. 队列固定效应 + HC3 ------------------------------------------
## scale = "auto_log": 全为正 -> 取 log（比值可直接解释成倍数）；
##                       有 0 -> 用 log1p（比值解释成 (值+1) 的倍数）
## scale = "identity": 原尺度（有界变量，如 purity）
adj_effect <- function(d, y, scale = "auto_log") {
  ok <- !is.na(d[[y]]) & !is.na(d$event) & !is.na(d$cohort)
  v  <- d[[y]][ok]
  if (scale == "auto_log") {
    if (all(v > 0)) { v2 <- log(v);   kind <- "ratio"; used <- "log" }
    else            { v2 <- log1p(v); kind <- "ratio"; used <- "log1p" }
  } else { v2 <- v; kind <- "diff"; used <- "identity" }
  
  dd <- data.frame(y = v2, fusion = d$event[ok], cohort = d$cohort[ok])
  m  <- lm(y ~ fusion + cohort, data = dd)
  ct <- coeftest(m, vcov. = vcovHC(m, type = "HC3"))
  b  <- ct["fusionFusion+", "Estimate"]
  se <- ct["fusionFusion+", "Std. Error"]
  p  <- ct["fusionFusion+", "Pr(>|t|)"]
  
  if (kind == "ratio")
    list(kind = kind, scale_used = used, est = exp(b), lo = exp(b - 1.96 * se),
         hi = exp(b + 1.96 * se), p = p, n = sum(ok))
  else
    list(kind = kind, scale_used = used, est = b, lo = b - 1.96 * se,
         hi = b + 1.96 * se, p = p, n = sum(ok))
}

METRICS <- data.table(
  var   = c("TMB", "HRD", "ploidy", "purity"),
  label = c("TMB (mut/Mb)", "HRD Score", "Ploidy", "Tumor Purity"),
  scale = c("auto_log", "auto_log", "auto_log", "identity"))

res <- rbindlist(lapply(seq_len(nrow(METRICS)), function(i) {
  a <- adj_effect(plot_df, METRICS$var[i], METRICS$scale[i])
  med <- tapply(plot_df[[METRICS$var[i]]], plot_df$event, median, na.rm = TRUE)
  data.table(var = METRICS$var[i], label = METRICS$label[i], kind = a$kind,
             scale_used = a$scale_used, n = a$n,
             median_fp = med[["Fusion+"]], median_fn = med[["Fusion-"]],
             est = a$est, lo = a$lo, hi = a$hi, p = a$p)
}))

## 图上/表里的那句话
res[, effect_lab := ifelse(kind == "ratio",
                           "Adj. GMR (95% CI)",
                           "Adj. mean difference (95% CI)")]
res[, value_lab := ifelse(kind == "ratio",
                          sprintf("%.2f (%.2f-%.2f)", est, lo, hi),
                          sprintf("%+.3f (%.3f to %+.3f)", est, lo, hi))]

res[, p_txt := ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p))]

res[, strip := paste0(label, "\n",
                      effect_lab, "\n",     # ← 这里多一个 \n 就是空行
                      value_lab, "\n",
                      p_txt)]

cat("\n== 队列固定效应校正后的效应量 ==\n")
print(res[, .(label, scale_used, n, median_fp = round(median_fp, 2), median_fn = round(median_fn, 2),
              est = signif(est, 3), lo = signif(lo, 3), hi = signif(hi, 3), p = signif(p, 3))])

## ---- 6. 图 ---------------------------------------------------------
plot_long <- plot_df |>
  pivot_longer(cols = METRICS$var, names_to = "var", values_to = "value") |>
  filter(!is.na(value)) |>
  left_join(res[, .(var, strip)], by = "var") |>
  mutate(strip = factor(strip, levels = res$strip))

pal <- c("Fusion-" = "#00468bff", "Fusion+" = "#ed0000ff")

p <- ggplot(plot_long, aes(x = event, y = value, fill = event)) +
  geom_boxplot(width = 0.6, outlier.size = 0.4, alpha = 0.85, colour = "grey20") +
  facet_wrap(~ strip, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal) +
  labs(x = NULL, y = NULL,
       title = "TCGA Pan-cancer: Fusion status vs genomic features") +
  theme_bw(base_size = 15) +
  theme(legend.position  = "none",
        strip.text       = element_text(size = 14, face = "bold",  lineheight = 0.95),
        axis.text.x      = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y      = element_text(size = 14),
        plot.title       = element_text(size = 25, face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank())
p
ggsave(file.path(out_dir, "Figure4E.pdf"), p, width = 15, height = 7, dpi = 600)
ggsave(file.path(out_dir, "Figure4E.png"), p, width = 15, height = 7, dpi = 600)


#Figure4F/G
suppressPackageStartupMessages({
  library(data.table); library(survival); library(ggplot2); library(dplyr)
})
standard_sex <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("M", "MALE")]   <- "Male"
  x[x %in% c("F", "FEMALE")] <- "Female"
  x[!x %in% c("Male", "Female")] <- NA_character_
  x
}
out_dir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4"
data_dir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据"
save_tab <- function(x, file, dir = out_dir)
  fwrite(x, file.path(dir, paste0(file, ".tsv")), sep = "\t")
## ---- 0. 参数：要改的东西全在这里 ------------------------------------
REQ_EXTRA <- c("Stage_ajcc", "Sex", "Age", "purity", "ploidy")

SEQ <- data.table(
  model = c("Unadjusted",
            "+ TMB", "+ Tumor Purity", "+ Ploidy","+ HRD score", 
            "+ Sex", "+ Age (>=60)", "+ Stage"),
  grp = factor(c("None",
                 rep("Genomic", 4), rep("Clinical", 3)),
               levels = c("None", "Genomic", "Clinical")),
  covs = list(character(0),
              "TMB", "purity", "ploidy", "HRD", 
              "Sex", "Age60", "Stage")
)
GRP_COL <- c(None = "grey35", Genomic = "#2E75B6",
             Clinical = "#B22222")

## ---- 1. 读入融合 ----------------------------------------------------
mf  <- fread(file.path(data_dir, "metafusion_calls.star.tsv"))
fus <- mf[sample_type %in% c("Tumor", "Tumor_Blood"),
          .(n_calls = uniqueN(Fusion_ID)), by = Sample_ID]
fus[, `:=`(Sample_ID = substr(Sample_ID, 1, 15), fusion = "Fusion+")]
fus <- fus[!duplicated(Sample_ID)]

## ---- 2. 读入临床 ----------------------------------------------------
files <- list.files(file.path(data_dir, "Clininfo"),
                    pattern = "^TCGA-.*_info\\.tsv$", full.names = TRUE)
clin <- rbindlist(lapply(files, function(f) {
  x    <- fread(f)
  keep <- intersect(c("Sample_ID", "OS_Time", "OS_Status", "PFS_Time", "PFS_Status",
                      "Stage_ajcc",  "Sex", 
                      "Age", "purity", "ploidy", "cohort"), names(x))
  x[, ..keep]
}), fill = TRUE)
clin <- clin[!grepl("NORMAL", cohort)]
clin[, Sample_ID := substr(as.character(Sample_ID), 1, 15)]
clin <- clin[!duplicated(Sample_ID)]

## ---- 3. TMB / HRD ---------------------------------------------------
mc3    <- fread(file.path(data_dir, "mc3.v0.2.8.PUBLIC.xena"))
silent <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
            "IGR", "RNA", "Splice_Region")
tmb <- mc3[!effect %in% silent, .(n_mut = .N), by = .(sample)]
setnames(tmb, "sample", "Sample_ID")
tmb[, `:=`(Sample_ID = substr(Sample_ID, 1, 15), TMB = n_mut / 38)]
tmb <- tmb[!duplicated(Sample_ID), .(Sample_ID, TMB)]

hrd_raw <- fread(file.path(data_dir, "TCGA.HRD_withSampleID.txt"))
hrd     <- as.data.frame(t(hrd_raw[, -1]), stringsAsFactors = FALSE)
colnames(hrd) <- hrd_raw[[1]]
hrd$Sample_ID <- substr(rownames(hrd), 1, 15)
hrd <- as.data.table(hrd)[, .(Sample_ID, HRD = as.numeric(HRD))]
hrd <- hrd[!is.na(HRD)][!duplicated(Sample_ID)]

## ---- 4. 合并 + 清洗（清洗必须在筛行之前） ---------------------------
d <- Reduce(function(a, b) merge(a, b, by = "Sample_ID", all.x = TRUE),
            list(clin, fus, tmb, hrd))
d[is.na(fusion), fusion := "Fusion-"]
d[, fusion := factor(fusion, levels = c("Fusion-", "Fusion+"))]

norm_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "NOT REPORTED", "NOT APPLICABLE", "STAGE X", "UNKNOWN")] <- NA
  x <- sub("^STAGE\\s*", "", x)
  x <- sub("[ABC]$", "", x)
  x[!x %in% c("I", "II", "III", "IV")] <- NA
  x
}

d[, `:=`(
  Stage_chr = norm_stage(Stage_ajcc),
  Sex_c     = standard_sex(Sex),          # 见下方函数
  Age       = suppressWarnings(as.numeric(Age)),
  purity    = suppressWarnings(as.numeric(purity)),
  ploidy    = suppressWarnings(as.numeric(ploidy)),
  OS_Time   = suppressWarnings(as.numeric(OS_Time)),
  PFS_Time  = suppressWarnings(as.numeric(PFS_Time)),
  OS_Status = suppressWarnings(as.integer(OS_Status)),
  PFS_Status= suppressWarnings(as.integer(PFS_Status))
)]
d[, `:=`(Age60 = as.integer(Age >= 60))]

## ---- 5. 分析集：OS / PFS 各自固定 -----------------------------------
req_os  <- c("OS_Time",  "OS_Status",  "Stage_chr",  "Sex_c")
req_pfs <- c("PFS_Time", "PFS_Status", "Stage_chr",  "Sex_c")
if (any(c("Age", "purity", "ploidy") %in% REQ_EXTRA))
  req_os <- c(req_os, "Age", "purity", "ploidy")[c(TRUE, TRUE, TRUE, TRUE, 
                                                   c("Age","purity","ploidy") %in% REQ_EXTRA)]

d_os  <- d[complete.cases(d[, ..req_os])] |>
  filter(!is.na(OS_Time), !is.na(OS_Status),
         !is.na(PFS_Time), !is.na(PFS_Status))
d_pfs <- d[complete.cases(d[, ..req_pfs])] |>
  filter(!is.na(OS_Time), !is.na(OS_Status),
         !is.na(PFS_Time), !is.na(PFS_Status))

## 因子水平：必须在筛行之后设，否则会留下空水平
mk_fac <- function(dt) {
  dt[, `:=`(
    Stage = factor(Stage_chr, levels = c("I", "II", "III", "IV")),
    Sex   = factor(Sex_c, levels = c("Female", "Male"))
  )][]
}
d_os  <- mk_fac(d_os);  d_pfs <- mk_fac(d_pfs)

## 逐步校正的共同集：所有 SEQ 里出现过的协变量都要齐
## 注意：必须取 SEQ$covs 这一列；unlist(SEQ) 会把 model/grp 两个因子列
##       一并摊平（因子变成整数码），得到 "1","2",... 这种列名而报错。
seq_vars  <- unique(unlist(SEQ$covs))
stopifnot(length(seq_vars) > 0,
          all(seq_vars %in% names(d_os)),
          all(seq_vars %in% names(d_pfs)))
d_seq_os  <- d_os[ complete.cases(d_os[,  ..seq_vars])]
d_seq_pfs <- d_pfs[complete.cases(d_pfs[, ..seq_vars])]

## ---- 6. 逐步校正模型（F / G） ---------------------------------------
run_seq <- function(dat, timevar, statvar, title, file, drop_cohort = NULL) {
  x0 <- copy(dat)
  if (!is.null(drop_cohort)) x0 <- x0[!cohort %in% drop_cohort]
  
  rows <- lapply(seq_len(nrow(SEQ)), function(i) {
    covs <- SEQ$covs[[i]]
    rhs  <- paste(c("fusion", covs, "strata(cohort)"), collapse = " + ")
    f    <- as.formula(paste0("Surv(", timevar, ", ", statvar, ") ~ ", rhs))
    fit  <- try(coxph(f, data = x0), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    cf <- coef(fit)
    if (!"fusionFusion+" %in% names(cf)) return(NULL)
    b  <- cf[["fusionFusion+"]]
    v  <- vcov(fit)["fusionFusion+", "fusionFusion+"]
    z  <- b / sqrt(v)                                   # ← 补上
    data.table(model = SEQ$model[i], grp = SEQ$grp[i],
               n = fit$n, events = fit$nevent,
               HR = exp(b),
               lo = exp(b - 1.96 * sqrt(v)),
               hi = exp(b + 1.96 * sqrt(v)),
               P  = 2 * pnorm(-abs(z)))                 # ← 只留一个 )
  })
  S <- rbindlist(rows)
  fp <- function(x) ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  S[, lab := paste0("P = ", fp(P))]
  
  S[, model := factor(model, levels = rev(as.character(SEQ$model)))]
  save_tab(S, file)
  
  p <- ggplot(S, aes(HR, model, colour = grp)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = .18, linewidth = .7) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = GRP_COL, name = "Adjusted for", drop = FALSE) +
    labs(x = "Fusion+ vs Fusion- hazard ratio", y = NULL, title = title) +
    theme_bw(base_size = 20) +
    geom_text(aes(x = hi, label = lab), hjust = -0.12, vjust = 0.5,
              size = 5, colour = "grey20", show.legend = FALSE) +
    ## 对称轴：0.05 – 20（以 HR = 1 为中心）
    scale_x_log10(breaks = c(0.1, 0.5, 1, 2),
                  labels = scales::label_number(accuracy = 0.01, drop0trailing = TRUE),
                  expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(xlim = c(1/2, 2)) +
    theme(legend.position = "right",
          legend.key = element_blank(),
          axis.text.x      = element_text(size = 20, hjust = 1),
          axis.text.y      = element_text(size = 20),
          plot.title       = element_text(size = 25, face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, paste0(file, ".pdf")), p, width = 12, height = 6, dpi = 600)
  ggsave(file.path(out_dir, paste0(file, ".png")), p, width = 12, height = 6, dpi = 600)
  S
}
mods <- c(TMB = "TMB", HRD = "HRD", ploidy = "ploidy", purity = "purity",
          Stage = "Stage", Age60 = "Age60")
mods <- mods[sapply(mods, function(m) sum(!is.na(d_os[[m]])) > 300)]

make_mod <- function(dt, m) {
  x <- dt[[m]]
  switch(m,
         "Stage"    = factor(x, levels = c("I","II","III","IV")),
         "Age60"    = factor(x, levels = c(0,1), labels = c("<60", ">=60")),
         { q <- quantile(x, c(1/3, 2/3), na.rm = TRUE)
         cut(x, c(-Inf, q, Inf), labels = c("low","mid","high")) })
}

S_os  <- run_seq(d_seq_os,  "OS_Time",  "OS_Status",
                 "TCGA: Fusion status and OS after sequential adjustment",                     
                 "Figure4F")
S_pfs <- run_seq(d_seq_pfs, "PFS_Time", "PFS_Status",
                 "TCGA: Fusion status and PFS after sequential adjustment",
                 "Figure4G")

## Figure 4H / 4I  —— fusion × modifier 交互（LRT）与层内效应
suppressPackageStartupMessages({
  library(data.table); library(survival); library(ggplot2)
})

data_dir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据"
out_dir  <- Sys.getenv("FIG4HI_OUT",
                       unset = "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. 融合 --------------------------------------------------------
mf  <- fread(file.path(data_dir, "metafusion_calls.star.tsv"))
fus <- mf[sample_type %in% c("Tumor", "Tumor_Blood"),
          .(n_calls = uniqueN(Fusion_ID)), by = Sample_ID]
fus[, `:=`(Sample_ID = substr(Sample_ID, 1, 15), fusion = "Fusion+")]
fus <- fus[!duplicated(Sample_ID)]

## ---- 2. 临床 --------------------------------------------------------
files <- list.files(file.path(data_dir, "Clininfo"),
                    pattern = "^TCGA-.*_info\\.tsv$", full.names = TRUE)
clin <- rbindlist(lapply(files, function(f) {
  x    <- fread(f)
  keep <- intersect(c("Sample_ID", "OS_Time", "OS_Status", "PFS_Time", "PFS_Status",
                      "Stage_ajcc", "Sex", "Age", "purity", "ploidy", "cohort"), names(x))
  x[, ..keep]
}), fill = TRUE)
clin <- clin[!grepl("NORMAL", cohort)]
clin[, Sample_ID := substr(as.character(Sample_ID), 1, 15)]
clin <- clin[!duplicated(Sample_ID)] |>
  filter(!is.na(OS_Time), !is.na(OS_Status),
         !is.na(PFS_Time), !is.na(PFS_Status))

## ---- 3. TMB / HRD ---------------------------------------------------
mc3    <- fread(file.path(data_dir, "mc3.v0.2.8.PUBLIC.xena"))
silent <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
            "IGR", "RNA", "Splice_Region")
tmb <- mc3[!effect %in% silent, .(n_mut = .N), by = .(sample)]
setnames(tmb, "sample", "Sample_ID")
tmb[, `:=`(Sample_ID = substr(Sample_ID, 1, 15), TMB = n_mut / 38)]
tmb <- tmb[!duplicated(Sample_ID), .(Sample_ID, TMB)]

hrd_raw <- fread(file.path(data_dir, "TCGA.HRD_withSampleID.txt"))
hrd     <- as.data.frame(t(hrd_raw[, -1]), stringsAsFactors = FALSE)
colnames(hrd) <- hrd_raw[[1]]
hrd$Sample_ID <- substr(rownames(hrd), 1, 15)
hrd <- as.data.table(hrd)[, .(Sample_ID, HRD = as.numeric(HRD))]
hrd <- hrd[!is.na(HRD)][!duplicated(Sample_ID)]

## ---- 4. 合并 --------------------------------------------------------
d <- Reduce(function(a, b) merge(a, b, by = "Sample_ID", all.x = TRUE),
            list(clin, fus, tmb, hrd))
d[is.na(fusion), fusion := "Fusion-"]

d[, fusion := factor(fusion, levels = c("Fusion-", "Fusion+"))]

## [FIX 2] 原始 Stage_ajcc 有大小写/空格差异，先规范化再抽数字
stage_num <- function(x) {
  x   <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(x))
  out[grepl("\\bIV\\b",  x)] <- 4L
  out[grepl("\\bIII\\b", x)] <- 3L
  out[grepl("\\bII\\b",  x)] <- 2L
  out[is.na(out) & grepl("\\bI\\b", x)] <- 1L
  out
}

## [FIX 3] purity / ploidy 显式转数值，否则一旦读成字符，quantile() 直接报错
standard_sex <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("M", "MALE")]   <- "Male"
  x[x %in% c("F", "FEMALE")] <- "Female"
  x[!x %in% c("Male", "Female")] <- NA_character_
  x
}
d[, `:=`(Age      = suppressWarnings(as.numeric(Age)),
         purity   = suppressWarnings(as.numeric(purity)),
         ploidy   = suppressWarnings(as.numeric(ploidy)),
         Sex_c    = standard_sex(Sex),
         Stage    = stage_num(Stage_ajcc),
         OS_Time  = suppressWarnings(as.numeric(OS_Time)),
         OS_Status  = suppressWarnings(as.integer(OS_Status)),
         PFS_Time   = suppressWarnings(as.numeric(PFS_Time)),
         PFS_Status = suppressWarnings(as.integer(PFS_Status)))]
d[, `:=`(Age60 = as.integer(Age >= 60))]


UNIFY <- TRUE
if (UNIFY) {
  sv <- c("TMB", "purity", "ploidy", "HRD", "Sex_c", "Age60", "Stage")
  d  <- d[complete.cases(d[, ..sv])]
}

## ---- 5. 修饰变量：按指定顺序 7 个，分面顺序即此顺序 ---------------
mods <- c(Stage = "Stage", Sex = "Sex_c", Age = "Age60",
          purity = "purity", ploidy = "ploidy", TMB = "TMB", HRD = "HRD")
mods <- mods[sapply(mods, function(m) sum(!is.na(d[[m]])) > 300)]

## [FIX 4] tertile 断点可能并列（purity/ploidy 大量重复值），
##         cut() 遇到非唯一断点会直接 error。这里先去重，退化时用中位数二分。
tertile <- function(x) {
  q <- unique(as.numeric(quantile(x, c(1/3, 2/3), na.rm = TRUE)))
  if (length(q) == 2) return(cut(x, c(-Inf, q, Inf), labels = c("low", "mid", "high")))
  q <- unique(as.numeric(quantile(x, 1/2, na.rm = TRUE)))
  if (length(q) == 1) return(cut(x, c(-Inf, q, Inf), labels = c("low", "high")))
  stop("modifier 无法分组")
}

make_mod <- function(x, m) {
  switch(m,
         "Stage"    = factor(x, levels = 1:4, labels = c("I", "II", "III", "IV")),
         "Age60"    = factor(x, levels = c(0, 1), labels = c("<60", ">=60")),
         "Sex_c"    = factor(x, levels = c("Female", "Male")),
         tertile(x))
}

## ---- 6. 交互模型 ----------------------------------------------------
fit_interaction <- function(dd, timevar, statvar, mvar = "Mf") {
  f0 <- as.formula(paste0("Surv(", timevar, ", ", statvar, ") ~ fusion + ",  mvar, " + strata(cohort)"))
  f1 <- as.formula(paste0("Surv(", timevar, ", ", statvar, ") ~ fusion * ", mvar, " + strata(cohort)"))
  m0 <- coxph(f0, data = dd)
  m1 <- coxph(f1, data = dd)
  a  <- anova(m0, m1)
  P  <- if ("Pr(>|Chi|)" %in% names(a)) a[["Pr(>|Chi|)"]][nrow(a)] else NA_real_
  Df <- a[["Df"]][nrow(a)]
  b  <- coef(m1); V <- vcov(m1)
  base <- "fusionFusion+"
  ## [FIX 5] 系数没对上就别静默返回 NULL，直接报出来
  if (!base %in% names(b))
    stop("找不到系数 ", base, "；现有：", paste(names(b), collapse = ", "))
  lv   <- levels(dd[[mvar]])
  rows <- lapply(seq_along(lv), function(i) {
    idx <- if (i == 1L) base else c(base, paste0(base, ":", mvar, lv[i]))
    if (!all(idx %in% names(b))) return(NULL)
    est <- sum(b[idx])
    se  <- sqrt(sum(V[idx, idx, drop = FALSE]))
    sub <- dd[dd[[mvar]] == lv[i], ]
    data.table(lvl = lv[i], n = nrow(sub),
               events = sum(sub[[statvar]] == 1, na.rm = TRUE),
               HR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
  })
  tab <- rbindlist(rows)
  tab[, `:=`(P_int = P, Df_int = Df, n_model = m1$n)]
  tab
}

run_panel <- function(dat, timevar, statvar, title, file, drop_cohort = NULL,
                      w = 15, h = 10) {
  x0 <- copy(dat)
  if (!is.null(drop_cohort)) x0 <- x0[!cohort %in% drop_cohort]
  
  ## 只去掉时间/状态缺失与 cohort 缺失，保留 time = 0（coxph 接受 0，不报错）
  x0 <- x0[!is.na(get(timevar)) & !is.na(get(statvar)) &
             get(timevar) >= 0 & !is.na(cohort)]
  
  flow <- rbindlist(lapply(mods, function(m) {
    dd <- x0[!is.na(get(m))]
    dd[, Mf := make_mod(get(m), m)]
    dd <- dd[!is.na(Mf)]
    tab <- fit_interaction(dd, timevar, statvar)
    tab[, mod := m]
    tab[, n_used := nrow(dd)]
    tab
  }), fill = TRUE)
  tab <- flow
  tab[, dropped_mod_na := nrow(x0) - n_used, by = mod]
  qtab <- unique(tab[, .(mod, P_int)])
  qtab[, q_int := p.adjust(P_int, method = "BH", n = sum(!is.na(P_int)))]
  tab[qtab, q_int := i.q_int, on = "mod"]
  lab_map <- c(Stage = "Stage", Sex_c = "Sex",
               Age60 = "Age", purity = "Tumor Purity", ploidy = "Ploidy",
               TMB = "TMB", HRD = "HRD score")
  tab[, modlab := factor(lab_map[mod], levels = unname(lab_map[mods]))]
  tab[, ev_model := sum(events), by = mod]
  tab[, lvl := factor(as.character(lvl),
                      levels = c("low", "mid", "high",
                                 "Female", "Male", "<60", ">=60",
                                 "I", "II", "III", "IV"))]
  
  tab[, p_lab := ifelse(is.na(P_int), "= NA",
                        ifelse(P_int < 0.001, "< 0.001", sprintf("= %.2g", P_int)))]
  tab[, q_lab := ifelse(is.na(q_int), "= NA",
                        ifelse(q_int < 0.001, "< 0.001", sprintf("= %.2g", q_int)))]
  grp_lab <- c(Stage = "4 levels", Sex_c = "2 groups", Age60 = "2 groups",
               purity = "tertiles", ploidy = "tertiles",
               TMB = "tertiles", HRD = "tertiles")
  tab[, strip := sprintf("%s\nP(int) %s\nq(int) %s (%s)",
                         modlab, p_lab, q_lab, grp_lab[mod])]
  tab[, strip := factor(strip, levels = unique(strip[order(modlab)]))]
  
  ncol_f <- min(4, length(unique(tab$mod)))
  nm <- unique(tab$n_model); ev <- unique(tab$ev_model)
  pH <- ggplot(tab, aes(x = HR, y = lvl)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                  width = 0.18, linewidth = 0.7) +
    geom_point(size = 2.8, colour = "grey15") +
    facet_wrap(~ strip, ncol = ncol_f, scales = "free_y") +
    ## 对称轴：0.05 – 20（以 HR = 1 为中心）
    scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5),
                  labels = scales::label_number(accuracy = 0.01, drop0trailing = TRUE),
                  expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(xlim = c(1/5, 5)) +
    labs(x = "Fusion+ vs Fusion- hazard ratio within stratum", y = NULL, title = title) +
    theme_bw(base_size = 20) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text   = element_text(size = 18, face = "bold"),
          axis.title.x  = element_text(size = 20, face = "bold"),
          axis.title.y  = element_text(size = 20, face = "bold"),
          plot.title   = element_text(size = 25, face = "bold"),
          plot.caption = element_text(size = 15, colour = "grey35", hjust = 0))
  ggsave(file.path(out_dir, paste0(file, ".pdf")), pH, width = w, height = h, dpi = 600)
  ggsave(file.path(out_dir, paste0(file, ".png")), pH, width = w, height = h, dpi = 600)
  fwrite(tab, file.path(out_dir, paste0(file, ".tsv")), sep = "\t")
  fwrite(unique(tab[, .(mod, modlab, n_model, ev_model, P_int, Df_int)]),
         file.path(out_dir, paste0(file, "_interaction.tsv")), sep = "\t")
  
  cat("\n==", file, "==  分析集 n =", nrow(x0), "\n")
  print(tab[, .(mod, lvl, n, events, HR = round(HR, 2), lo = round(lo, 2),
                hi = round(hi, 2), P_int = signif(P_int, 3), Df_int, n_model)])
  invisible(pH)
}

pH_os  <- run_panel(d, "OS_Time",  "OS_Status",
                    "P(interaction): likelihood-ratio test of the fusion status × modifier interaction in OS",
                    "Figure4H")
pH_pfs <- run_panel(d, "PFS_Time", "PFS_Status",
                    "P(interaction): likelihood-ratio test of the fusion status × modifier interaction in PFS",
                    "Figure4I")
