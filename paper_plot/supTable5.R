# ============================================================
# 0. 依赖
# ============================================================
library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(patchwork)
library(future)
library(future.apply)
library(tidyverse)
set.seed(2024)
plan(multisession)
# ------------------------------------------------------------
# 1. 读取 + 预处理 fusion
# ------------------------------------------------------------
file_meta  <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
metafusion <- fread(file_meta, sep = "\t")

Credibility_gene <- metafusion |>
  dplyr::filter(sample_type %in% c("Tumor", "Tumor_Blood")) 
cohort1 <- fread("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/cohort_meta_check.csv")

# ------------------------------------------------------------
# 2. cohort → cancer 映射
# ------------------------------------------------------------
cohort2cancer <- unique(cohort1[, c("id", "cancer_type")])
setnames(cohort2cancer, c("id", "cancer_type"), c("cohort", "cancer"))

# ------------------------------------------------------------
# 3. OS 数据预处理
# ------------------------------------------------------------
cohort_sum_OS <- cohort1 |> filter(OS == TRUE) |> pull(id) |> unique()

OS_all <- data.table()
for (cohort_id in cohort_sum_OS) {
  Need <- fread(paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/",
                       cohort_id, "_info.tsv")) |>
    filter(!is.na(OS_Time) & !is.na(OS_Status))
  
  ICB_pre_on_OS <- cohort1 |> filter(OS == TRUE & time == "Pre/On") |> pull(id)
  
  if (cohort_id %in% ICB_pre_on_OS) {
    if (cohort_id == "PRJNA356761") {
      Need1 <- Need |> filter(grepl("Pre", Sample_ID))
    } else if (cohort_id == "EGAD00001006282") {
      Need1 <- Need |> filter(grepl("_SCREEN", Run))
    } else {
      Need1 <- Need |> filter(grepl("Pre", Response2))
    }
  } else if (cohort_id %in% c("PRJNA482620", "PRJNA744780")) {
    Need1 <- Need |> filter(grepl("Pre", Response2))
  } else {
    Need1 <- Need
  }
  
  columns_to_select <- c("Patient_ID","Run","Age","Sex","Treatment",
                         "Therapy_Type","OS_Time","OS_Status")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |> select(all_of(existing_columns))
  
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  for (col in missing_columns) Need2[[col]] <- NA
  
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort     <- cohort_id
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
    filter(!Treatment %in% c("Docetaxel", "Observation", "Sunitinib"))
}

# ------------------------------------------------------------
# 4. PFS 数据预处理
# ------------------------------------------------------------
cohort_sum_PFS <- cohort1 |> filter(PFS == TRUE) |> pull(id) |> unique()

PFS_all <- data.table()
for (cohort_id in cohort_sum_PFS) {
  Need <- fread(paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/",
                       cohort_id, "_info.tsv"))
  
  ICB_pre_on_PFS <- cohort1 |> filter(PFS == TRUE & time == "Pre/On") |> pull(id)
  
  if (cohort_id %in% ICB_pre_on_PFS) {
    if (cohort_id == "PRJNA356761") {
      Need1 <- Need |> filter(grepl("Pre", Sample_ID))
    } else if (cohort_id == "EGAD00001006282") {
      Need1 <- Need |> filter(grepl("_SCREEN", Run))
    } else {
      Need1 <- Need |> filter(grepl("Pre", Response2))
    }
  } else if (cohort_id %in% c("PRJNA482620", "PRJNA744780")) {
    Need1 <- Need |> filter(grepl("Pre", Response2))
  } else {
    Need1 <- Need
  }
  
  columns_to_select <- c("Patient_ID","PFS_Time","PFS_Status",
                         "Run","Treatment","Age","Sex")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    dplyr::select(all_of(existing_columns)) |>
    filter(!is.na(PFS_Time) & !is.na(PFS_Status))
  
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  for (col in missing_columns) Need2[[col]] <- NA
  
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort     <- cohort_id
  
  PFS_all <- rbind(Need2, PFS_all) |>
    filter(!Treatment %in% c("Docetaxel", "Observation", "Sunitinib"))
}

# ------------------------------------------------------------
# 5. 构建 Tumor_ICB_OS / PFS
# ------------------------------------------------------------
Tumor_ICB_OS <- Credibility_gene |>
  as.data.table() |>
  dplyr::filter(sample_type != "Normal") |>
  dplyr::filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  dplyr::filter(cohort %in% cohort_sum_OS) |>
  select(Sample_ID, gene1, gene2, cohort) |>
  rename(Run = Sample_ID) |>
  distinct()

Tumor_ICB_PFS <- Credibility_gene |>
  as.data.table() |>
  dplyr::filter(sample_type != "Normal") |>
  dplyr::filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  dplyr::filter(cohort %in% cohort_sum_PFS) |>
  select(Sample_ID, gene1, gene2, cohort) |>
  rename(Run = Sample_ID) |>
  distinct()

OS_all        <- merge(OS_all,        cohort2cancer, by = "cohort", all.x = TRUE)
PFS_all       <- merge(PFS_all,       cohort2cancer, by = "cohort", all.x = TRUE)
Tumor_ICB_OS  <- merge(Tumor_ICB_OS,  cohort2cancer, by = "cohort", all.x = TRUE)
Tumor_ICB_PFS <- merge(Tumor_ICB_PFS, cohort2cancer, by = "cohort", all.x = TRUE)

# ------------------------------------------------------------
# 6. 单基因分析函数（★ 全程患者级）
# ------------------------------------------------------------
analyze_gene_group <- function(gene, side, fusion_sub, surv_sub_base,
                               time_col, status_col,
                               out_dir, suffix,
                               analysis_label, min_pos = 3,
                               cohort_filter = NULL) {
  library(data.table); library(dplyr); library(survival)
  library(survminer); library(ggplot2); library(tidyr)
  library(stringr); library(patchwork)
  
  gene_vec  <- str_split(gene, "\\s*,\\s*")[[1]]
  gene_name <- gsub(",\\s*", "_", gene)
  
  if (!is.null(cohort_filter)) {
    fusion_sub    <- fusion_sub    |> dplyr::filter(cohort %in% cohort_filter)
    surv_sub_base <- surv_sub_base |> dplyr::filter(cohort %in% cohort_filter)
    if (nrow(fusion_sub) == 0 || nrow(surv_sub_base) == 0) return(NULL)
  }
  
  # ★ 改动 1：只看对应端。5' = gene1，3' = gene2
  hit_rows <- if (side == "5'") fusion_sub[gene1 %in% gene_vec]
  else              fusion_sub[gene2 %in% gene_vec]
  
  fusion_select <- hit_rows |> distinct(Patient_ID, cohort, cancer)
  if (nrow(fusion_select) == 0) return(NULL)
  
  keep_cancers <- unique(fusion_select$cancer)
  surv_use <- surv_sub_base |> filter(cancer %in% keep_cancers)
  if (nrow(surv_use) == 0) return(NULL)
  
  surv_use <- surv_use |>
    mutate(Patient_ID = as.character(Patient_ID)) |>
    mutate(is_pos = Patient_ID %in% fusion_select$Patient_ID) |>
    group_by(Patient_ID) |>
    mutate(event = if (any(is_pos)) "Fusion+" else "Fusion-") |>
    ungroup() |>
    select(Patient_ID, cohort, cancer, event,
           all_of(time_col), all_of(status_col)) |>
    distinct()
  
  hn <- length(unique(surv_use$Patient_ID[surv_use$event == "Fusion+"]))
  ln <- length(unique(surv_use$Patient_ID[surv_use$event == "Fusion-"]))
  if (hn < min_pos) return(NULL)
  
  surv_use$event <- factor(surv_use$event, levels = c("Fusion-", "Fusion+"))
  
  use_strata <- length(unique(surv_use$cohort)) >= 2
  rhs <- if (use_strata) "event + strata(cohort)" else "event"
  fml <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", rhs))
  fit_cox <- tryCatch(coxph(fml, data = surv_use), error = function(e) NULL)
  if (is.null(fit_cox)) return(NULL)
  
  cf <- summary(fit_cox)$coefficients
  if (!"eventFusion+" %in% rownames(cf)) return(NULL)
  beta <- cf["eventFusion+", "coef"]
  se   <- cf["eventFusion+", "se(coef)"]
  pval <- cf["eventFusion+", "Pr(>|z|)"]
  HR   <- exp(beta)
  low  <- exp(beta - qnorm(0.975) * se)
  up   <- exp(beta + qnorm(0.975) * se)
  
  data.table(
    `Fusion Gene` = gene_name,
    Side          = side,            # ★ 改动 2：多一列，不写进基因名
    Cancer        = analysis_label,
    Cohort        = paste(sort(unique(surv_use$cohort)), collapse = ", "),
    `Fusion+`     = hn,
    `Fusion-`     = ln,
    p_raw         = pval,
    HR            = round(HR, 2),
    CI            = paste0(round(low, 2), " - ", round(up, 2)),
    Status        = suffix
  )
}

# ------------------------------------------------------------
# 7. 主循环（★ 患者级 cohort 门槛）
# ------------------------------------------------------------
gene_long <- function(dt) {
  out <- rbind(
    dt[, .(gene = gene1, Side = "5'", Patient_ID, cohort, cancer)],
    dt[, .(gene = gene2, Side = "3'", Patient_ID, cohort, cancer)]
  )
  out[!is.na(gene) & nzchar(gene)]
}

run_one_level <- function(fusion_all, surv_all,
                          time_col, status_col, suffix, out_base,
                          min_pos = 3) {
  library(future.apply); library(data.table); library(dplyr)
  
  clean_df <- function(dt) {
    dt <- as.data.table(dt)
    dt <- dt[, !duplicated(names(dt)), with = FALSE]
    for (drop in c("cancer.x", "cancer.y")) if (drop %in% names(dt)) dt[, (drop) := NULL]
    dt
  }
  fusion_all <- clean_df(fusion_all)
  surv_all   <- clean_df(surv_all)
  stopifnot("cancer" %in% names(fusion_all), "cancer" %in% names(surv_all))
  
  run2patient <- unique(surv_all[, .(Run, Patient_ID)])
  run2patient[, `:=`(Run = as.character(Run), Patient_ID = as.character(Patient_ID))]
  stopifnot(uniqueN(run2patient$Run) == nrow(run2patient))
  fusion_all[, Run := as.character(Run)]
  fusion_all <- merge(fusion_all, run2patient, by = "Run", all = FALSE)
  if (nrow(fusion_all) == 0) {
    message(sprintf("[%s] 无任何 fusion Run 有生存数据，跳过", suffix))
    return(data.table())
  }
  
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
  res_all <- list()
  for (cs in unique(fusion_all$cancer)) {
    if (is.na(cs)) next
    fusion_c <- fusion_all[cancer == cs]
    surv_c   <- surv_all[cancer == cs]
    if (nrow(fusion_c) == 0 || nrow(surv_c) == 0) next
    
    gl_c <- gene_long(fusion_c)
    # ★ 门槛按 (gene, Side, cohort) 的患者数算
    gene_cohort_n_c <- gl_c[, .(n = uniqueN(Patient_ID)),
                            by = .(gene, Side, cohort)]
    gene_ok <- unique(gene_cohort_n_c[n >= min_pos, .(gene, Side)])
    if (!nrow(gene_ok)) {
      message(sprintf("[%s|%s] 无 gene–side 达到 patient-level cohort n >= %d，跳过",
                      suffix, cs, min_pos))
      next
    }
    message(sprintf("[%s|%s] %d gene–side combos (patient-level cohort n >= %d)",
                    suffix, cs, nrow(gene_ok), min_pos))
    
    lst <- future_lapply(
      seq_len(nrow(gene_ok)),
      function(i) {
        g   <- gene_ok$gene[i]
        sd  <- gene_ok$Side[i]
        ok_cohorts <- gene_cohort_n_c[gene == g & Side == sd & n >= min_pos,
                                      unique(cohort)]
        if (length(ok_cohorts) == 0) return(NULL)
        analyze_gene_group(
          g, sd, fusion_c, surv_c, time_col, status_col,
          out_base, paste0(suffix, "_", cs),
          analysis_label = cs, min_pos = min_pos,
          cohort_filter = ok_cohorts)
      },
      future.seed = TRUE
    )
    res_all <- c(res_all, lst[!sapply(lst, is.null)])
  }
  
  # ---------- 维度 2：Pan-Cancer ----------
  gl_all <- gene_long(fusion_all)
  gene_cohort_n <- gl_all[, .(n_pos = uniqueN(Patient_ID)),
                          by = .(gene, Side, cohort, cancer)]
  gene_ok <- unique(gene_cohort_n[n_pos >= min_pos, .(gene, Side)])
  
  if (nrow(gene_ok)) {
    message(sprintf("[%s|Pan] %d gene–side combos (patient-level cohort n >= %d)",
                    suffix, nrow(gene_ok), min_pos))
    
    lst <- future_lapply(
      seq_len(nrow(gene_ok)),
      function(i) {
        g  <- gene_ok$gene[i]
        sd <- gene_ok$Side[i]
        gc <- gene_cohort_n[gene == g & Side == sd]
        ok_cohorts <- gc[n_pos >= min_pos, unique(cohort)]
        if (length(ok_cohorts) == 0) return(NULL)
        cs_keep <- gc[cohort %in% ok_cohorts, unique(cancer)]
        if (length(cs_keep) < 2) return(NULL)
        # ★ 总阳性数也按端数
        total_pos <- gl_all[cohort %in% ok_cohorts & gene == g & Side == sd,
                            uniqueN(Patient_ID)]
        if (total_pos < min_pos) return(NULL)
        label <- paste0("Pan-Cancer(",
                        paste(sort(cs_keep), collapse = ","), ")")
        analyze_gene_group(
          g, sd, fusion_all, surv_all, time_col, status_col,
          out_base, paste0(suffix, "_PanCancer"),
          analysis_label = label, min_pos = min_pos,
          cohort_filter = ok_cohorts)
      },
      future.seed = TRUE
    )
    res_all <- c(res_all, lst[!sapply(lst, is.null)])
  } else {
    message(sprintf("[%s|Pan] 无 gene–side 达到 patient-level cohort n >= %d，跳过",
                    suffix, min_pos))
  }
  
  rbindlist(res_all, use.names = TRUE, fill = TRUE)
}

# ------------------------------------------------------------
# 8. 执行 OS / PFS（天然分开，各自独立 FDR）
# ------------------------------------------------------------
base_dir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/附表"

OS_all_cancers  <- run_one_level(Tumor_ICB_OS,  OS_all,
                                 "OS_Time",  "OS_Status", "OS",
                                 base_dir, min_pos = 3)
PFS_all_cancers <- run_one_level(Tumor_ICB_PFS, PFS_all,
                                 "PFS_Time", "PFS_Status", "PFS",
                                 base_dir, min_pos = 3)

# ------------------------------------------------------------
# 9. FDR：Per-Cancer 各癌种内合并；Pan-Cancer 合并后 FDR
FDR_BY_SIDE <- FALSE   # TRUE = 5'/3' 各自一族；FALSE = 两端合并一族（与 Table S7 一致）

apply_fdr <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  dt <- as.data.table(dt)
  dt[, is_pan := grepl("^Pan-Cancer", Cancer)]
  
  if (FDR_BY_SIDE) {
    dt[is_pan == FALSE, `q (BH)` := p.adjust(p_raw, "BH"), by = .(Status, Cancer, Side)]
    dt[is_pan == TRUE,  `q (BH)` := p.adjust(p_raw, "BH"), by = .(Status, Side)]

    dt[is_pan == FALSE, n_tests_family := sum(!is.na(p_raw)), by = .(Status, Cancer, Side)]
    dt[is_pan == TRUE,  n_tests_family := sum(!is.na(p_raw)), by = .(Status, Side)]
  } else {
    dt[is_pan == FALSE, `q (BH)` := p.adjust(p_raw, "BH"), by = .(Status, Cancer)]
    dt[is_pan == TRUE,  `q (BH)` := p.adjust(p_raw, "BH"), by = .(Status)]
    dt[is_pan == FALSE, n_tests_family := sum(!is.na(p_raw)), by = .(Status, Cancer)]
    dt[is_pan == TRUE,  n_tests_family := sum(!is.na(p_raw)), by = .(Status)]
  }
  
  setorder(dt, Status, is_pan, Cancer, Side, p_raw)
  dt[, p                := signif(p_raw, 3)]
  dt[, `q (BH)`        := signif(`q (BH)`, 3)]
  dt[, c("is_pan", "p_raw") := NULL]
  
  setcolorder(dt, c("Fusion Gene", "Side", "Cancer", "Cohort",
                    "Fusion+", "Fusion-",
                    "p", "q (BH)",  "n_tests_family",
                    "HR", "CI", "Status"))
  dt[]
}
OS_all_cancers  <- apply_fdr(OS_all_cancers)
PFS_all_cancers <- apply_fdr(PFS_all_cancers)
# ------------------------------------------------------------
# 10. 合并输出
# ------------------------------------------------------------
All_cancers <- rbindlist(
  list(OS_all_cancers, PFS_all_cancers),
  use.names = TRUE, fill = TRUE
) 

write_xlsx(
  All_cancers,
  path = "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/附表/TableS5.xlsx"
)

plan(sequential)
message("All done!")
