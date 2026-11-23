library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tidyr)
library(stringr)
library(patchwork)
library(forestploter)
library(grid)
set.seed(2024)
out_dir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图2/"
cohort <- fread("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/cohort_meta_check.csv")
km_os_by_group <- function(data, 
                           group_levels = c( "Fusion-","Fusion+"),
                           choose = "OS",
                           strata_var = "cohort") {
  
  data$event <- factor(data$event, levels = group_levels)   # ← 参照组放前面
  group_counts <- table(data$event)
  ln <- group_counts[group_levels[1]]   
  hn <- group_counts[group_levels[2]]    
  if(hn < 2) next
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
  p.val1 <- round(cox_summary$coefficients[1, "Pr(>|z|)"], 3)
  
  # 提取 HR 和 95% CI (注意：这里提取的是指数化后的值)
  HR <- round(cox_summary$coefficients[1, "exp(coef)"], 2)
  low95 <- round(cox_summary$conf.int[1, "lower .95"], 2)
  up95 <- round(cox_summary$conf.int[1, "upper .95"], 2)
  cancer_type   <- cohort$cancer_type[cohort$id == cohort_choose]
  
  results <- data.frame(
    hn = hn,
    ln = ln,
    HR    = HR,
    low95 = low95,
    up95  = up95 ,
    p     = p.val1 
  )
  return(results)
}
tm <- forest_theme(
  base_size = 10,
  colhead  = list(fg_params = list(hjust = 0)),
  core     = list(padding = unit(c(5, 5), "mm")),
  ci_pch = 15, ci_col = "black", ci_fill = "black",
  ci_alpha = 0.8, ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.2,
  refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey20"),
  footnote_gp = gpar(cex = 0.8, fontface = "plain", col = "black"),
  vertline_lwd = 1, vertline_lty = "dashed", vertline_col = "grey20"
)
file_meta <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
metafusion <- fread(file_meta, sep = "\t")

Credibility_gene <- metafusion |>
  filter(sample_type %in% c("Tumor", "Tumor_Blood")) |>
  as.data.table()

#supFigure2A
cohort_OS <- cohort |> filter(OS == TRUE) 
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
  cancer_type   <- cohort$cancer_type[cohort$id == cohort_id]
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort <- cohort_id
  Need2$cancer <- cancer_type
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
  
  OS_all <- rbind(Need3,OS_all) |> filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
}
Tumor_ICB_OS <- Credibility_gene |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_sum_OS) |> 
  select(Sample_ID, gene1, gene2, cohort)|>
  rename(Run = Sample_ID) |>
  distinct()
cancer_type1 = "ALL"
cohort_ALL<- cohort |> filter(OS == TRUE )
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
  select(Patient_ID, event, OS_Time, OS_Status, cohort, cancer) |>
  distinct() |>
  filter(!is.na(OS_Status) & !is.na(OS_Time) )

result_all <- data.table()
cohort_cancer <- fusion_data_subset$cancer |> unique()
for (cohort_choose in cohort_cancer) {
  data_choose <- fusion_data_subset |> filter(cancer == cohort_choose)
  result_choose <- km_os_by_group(data = data_choose, 
                                  choose = "OS"
  )
  result_choose$cohort_name <- cohort_choose
  result_all <- rbind(result_all, result_choose)
}
result_all <-  result_all |> 
  filter(!is.infinite(up95) & !is.na(up95) & cohort_name != "BRCA") 
result_all[, `q(BH)` := p.adjust(p, method = "BH")]
out_path <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图2/supFigure2A.pdf"
make_forest_hr <- function(dt, title, out_pdf, max_rows = 80, theme = NULL,
                           x_max = 20, width = NULL) {
  
  if (nrow(dt) == 0) { message("empty: ", title); return(invisible(NULL)) }
  dt <- data.table::copy(data.table::as.data.table(dt))
  
  # 1. 兼容不同版本的输入列名
  if ("cohort_name" %in% names(dt) && !"cohort" %in% names(dt)) data.table::setnames(dt, "cohort_name", "cohort")
  if ("hn" %in% names(dt) && !"n_pos" %in% names(dt)) data.table::setnames(dt, "hn", "n_pos")
  if ("ln" %in% names(dt) && !"n_neg" %in% names(dt)) data.table::setnames(dt, "ln", "n_neg")
  if ("low95" %in% names(dt) && !"low" %in% names(dt)) data.table::setnames(dt, "low95", "low")
  if ("up95" %in% names(dt) && !"up" %in% names(dt)) data.table::setnames(dt, "up95", "up")
  if ("q(BH)" %in% names(dt) && !"q" %in% names(dt)) data.table::setnames(dt, "q(BH)", "q")
  if (!"is_summary" %in% names(dt)) dt[, is_summary := FALSE]
  if (!"q" %in% names(dt)) dt[, q := NA_real_]
  
  # 2. 排序并截断
  dt[, row_order := data.table::fifelse(is_summary, 0L, 1L)]
  data.table::setorder(dt, row_order)
  dt <- dt[1:min(.N, max_rows)]
  
  # 3. 文本：始终显示真实的 HR 与 95% CI（上限不截成 20）
  dt[, `Fusion+/Fusion-` := paste0(n_pos, "/", n_neg)]
  dt[, `HR (95% CI)` := sprintf("%.2f (%.2f\u2013%.2f)", HR, low, up)]
  dt[, p_lab := ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))]
  dt[, q_lab := ifelse(is.na(q), "", ifelse(q < 0.001, "<0.001", sprintf("%.3f", q)))]
  
  # 森林图列宽度：嫌窄/宽就改这里的空格数
  dt[, ` ` := paste(rep(" ", 60), collapse = "")]
  
  data_rows <- data.frame(
    Cohort                 = dt$cohort,
    `Binary fusion status` = dt$`Fusion+/Fusion-`,
    ` `                    = dt[[" "]],
    `HR (95% CI)`          = dt$`HR (95% CI)`,
    p                      = dt$p_lab,
    `q (BH)`               = dt$q_lab,
    check.names            = FALSE, stringsAsFactors = FALSE
  )
  
  # 伪表头行：直接复制 data_rows 的结构，保证列名/列序完全一致
  header_row <- data_rows[1, , drop = FALSE]
  header_row[] <- ""
  header_row[[2]] <- "Fusion+/Fusion-"
  header_row[[3]] <- dt[[" "]][1]
  plot_dt <- rbind(header_row, data_rows)   
  # 4. 轴：以 1 为中心左右对称（对数轴下 1/x_max 与 x_max 关于 1 对称）
  use_log <- all(c(dt$HR, dt$low, dt$up) > 0, na.rm = TRUE)
  X_HI <- x_max
  X_LO <- if (use_log) 1 / X_HI else -X_HI
  xlim <- c(X_LO, X_HI)
  
  ticks <- if (use_log) c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
  if (!is.null(ticks)) {
    ticks <- ticks[ticks >= X_LO & ticks <= X_HI]
    if (length(ticks) < 3) ticks <- NULL
  }
  
  # 关键：不裁剪数据。forestploter 会把超出 xlim 的 CI 画到轴端并自动加箭头
  est_vec   <- c(NA, dt$HR)
  lower_vec <- c(NA, dt$low)
  upper_vec <- c(NA, dt$up)
  
  n_clip <- sum(dt$low < X_LO | dt$up > X_HI, na.rm = TRUE)
  
  foot <- paste0(
    "\n\n\n\n\n\n",                       # 空行：把脚注推到横轴数字与 xlab 之下
    "ICB Cox: Surv ~ Binary fusion status, cohort as strata.\n",
    "BRCA were not analysed because the Fusion- group was too small (n = 1)\n",
    "Not shown (HR = Inf / non−estimable): GBM",
    "q: Benjamini\u2013Hochberg FDR.\n"
  )
  
  wrap_lines <- function(txt, width = 110) {
    lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
    lines <- unlist(lapply(lines, function(x) {
      if (nchar(x) == 0) return("")
      strwrap(x, width = width)
    }))
    paste(lines, collapse = "\n")
  }
  footnote_txt <- wrap_lines(foot, width = 110)
  
  # 5. 主题：斑马纹 + 左对齐 + 黑色方块 + 显示横轴刻度（原来 xaxis_gp col = NA 把刻度藏了）
  if (is.null(theme)) {
    row_fill <- rep(c("white", "#E8F3E8"), length.out = nrow(plot_dt))
    theme <- forestploter::forest_theme(
      base_size    = 11,
      base_family  = "sans",
      ci_pch       = 15,
      ci_col       = "black",
      ci_fill      = "black",
      ci_lwd       = 1.4,
      ci_Theight   = 0.15,
      refline_gp   = grid::gpar(lwd = 1, lty = "dashed", col = "grey25"),
      vertline_lwd = 0.5,
      vertline_lty = "solid",
      vertline_col = "grey90",
      xaxis_gp     = grid::gpar(fontsize = 7.5, col = "grey30"),
      xlab_gp      = grid::gpar(fontsize = 9, col = "grey30"),
      title_just   = "left",
      core = list(
        bg_params = list(fill = row_fill, col = NA, lwd = 0),
        fg_params = list(hjust = 0, x = 0.02, fontsize = 10)
      ),
      colhead = list(
        fg_params = list(hjust = 0, x = 0.02, fontsize = 10, fontface = "bold")
      )
    )
  }
  
  vlines <- if (use_log) c(0.1, 0.5, 2, 5, 10) else NULL
  if (!is.null(vlines)) vlines <- vlines[vlines > xlim[1] & vlines < xlim[2]]
  
  p <- forestploter::forest(
    plot_dt,
    est = est_vec, lower = lower_vec, upper = upper_vec,
    ci_column = 3,
    ref_line = 1,
    vert_line = vlines,
    xlim = xlim,
    x_trans = if (use_log) "log" else "none",
    ticks_at = ticks,
    xlab = if (use_log) "Hazard ratio (log scale)" else "Hazard ratio",
    title = title,
    footnote = footnote_txt,
    theme = theme
  )
  
  p <- forestploter::edit_plot(p, row = 1, col = 2, which = "text",
                               gp = grid::gpar(fontface = "italic", cex = 0.9))
  
  n_sum <- sum(dt$is_summary)
  if (n_sum > 0 && n_sum < nrow(dt)) {
    p <- forestploter::add_border(p, part = "body", row = n_sum + 2,
                                  gp = grid::gpar(lwd = 0.8, col = "grey40"))
  }
  
  # 6. 尺寸：按实际内容量出来，不再用行数估，避免大片留白
  wh <- forestploter::get_wh(p)
  w  <- if (is.null(width)) as.numeric(wh["width"]) + 0.20 else width
  h  <- as.numeric(wh["height"]) + 0.35
  
  grDevices::pdf(out_pdf, width = w, height = h, onefile = FALSE)
  print(p)
  grDevices::dev.off()
  
  invisible(list(plot = p, width = w, height = h, n_clipped = n_clip))
}
make_forest_hr(result_all, 
               title = "ICB OS: Fusion+ vs Fusion-", 
               out_pdf = out_path)



#supFigure2B
cohort_PFS <- cohort |> filter(PFS == TRUE) 
cohort_sum_PFS <- (cohort_PFS$id) |> unique()
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
    }else{
      Need1 <- Need |>
        filter(grepl("Pre", Response2))
    }}else if(cohort_id  %in% c("PRJNA482620", "PRJNA744780")){
      Need1 <- Need |>
        filter(grepl("Pre", Response2))
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
      Need1 <- Need
    }
  
  columns_to_select <- c("Patient_ID", "PFS_Time", "PFS_Status", "Run", "Treatment", "Age", "Sex")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    select(all_of(existing_columns)) |>
    filter(!is.na(PFS_Time) & !is.na(PFS_Status))
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  if (length(missing_columns) > 0) {
    for (col in missing_columns) {
      Need2[[col]] <- NA
    }
  }
  cancer_type   <- cohort$cancer_type[cohort$id == cohort_id]
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort <- cohort_id
  Need2$cancer <- cancer_type
  PFS_all <- rbind(Need2,PFS_all) |> filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
}
Tumor_ICB_PFS <- Credibility_gene |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_sum_PFS) |> 
  select(Sample_ID, gene1, gene2,cohort)|>
  rename(Run = Sample_ID) |>
  distinct()
cancer_type1 = "ALL"
cohort_ALL<- cohort |> filter(PFS == TRUE )
cohort_select <- (cohort_ALL$id) |> unique()
fusion_select <- Tumor_ICB_PFS |>
  filter(cohort %in% cohort_select) |> 
  as.data.table()
cohort_select <- fusion_select$cohort |> unique()

fusion_data_subset <- PFS_all |>
  filter(cohort %in% cohort_select) |>
  mutate(event = case_when(
    Run %in% fusion_select$Run  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                           TRUE ~ event))  |>
  select(Patient_ID, event, PFS_Time, PFS_Status, cohort,cancer) |>
  distinct() |>
  filter(cancer != "SGC")
result_all <- data.table()
cohort_cancer <- fusion_data_subset$cancer |> unique()
for (cohort_choose in cohort_cancer) {
  data_choose <- fusion_data_subset |> filter(cancer == cohort_choose)
  result_choose <- km_os_by_group(data = data_choose, 
                                  choose = "PFS"
  )
  result_choose$cohort_name <- cohort_choose
  result_all <- rbind(result_all, result_choose)
}
result_all <-  result_all |> filter(!is.infinite(up95) & !is.na(up95))
result_all[, `q(BH)` := p.adjust(p, method = "BH")]
out_path <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图2/supFigure2B.pdf"
make_forest_hr <- function(dt, title, out_pdf, max_rows = 80, theme = NULL,
                           x_max = 20, width = NULL) {
  
  if (nrow(dt) == 0) { message("empty: ", title); return(invisible(NULL)) }
  dt <- data.table::copy(data.table::as.data.table(dt))
  
  # 1. 兼容不同版本的输入列名
  if ("cohort_name" %in% names(dt) && !"cohort" %in% names(dt)) data.table::setnames(dt, "cohort_name", "cohort")
  if ("hn" %in% names(dt) && !"n_pos" %in% names(dt)) data.table::setnames(dt, "hn", "n_pos")
  if ("ln" %in% names(dt) && !"n_neg" %in% names(dt)) data.table::setnames(dt, "ln", "n_neg")
  if ("low95" %in% names(dt) && !"low" %in% names(dt)) data.table::setnames(dt, "low95", "low")
  if ("up95" %in% names(dt) && !"up" %in% names(dt)) data.table::setnames(dt, "up95", "up")
  if ("q(BH)" %in% names(dt) && !"q" %in% names(dt)) data.table::setnames(dt, "q(BH)", "q")
  if (!"is_summary" %in% names(dt)) dt[, is_summary := FALSE]
  if (!"q" %in% names(dt)) dt[, q := NA_real_]
  
  # 2. 排序并截断
  dt[, row_order := data.table::fifelse(is_summary, 0L, 1L)]
  data.table::setorder(dt, row_order)
  dt <- dt[1:min(.N, max_rows)]
  
  # 3. 文本：始终显示真实的 HR 与 95% CI（上限不截成 20）
  dt[, `Fusion+/Fusion-` := paste0(n_pos, "/", n_neg)]
  dt[, `HR (95% CI)` := sprintf("%.2f (%.2f\u2013%.2f)", HR, low, up)]
  dt[, p_lab := ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))]
  dt[, q_lab := ifelse(is.na(q), "", ifelse(q < 0.001, "<0.001", sprintf("%.3f", q)))]
  
  # 森林图列宽度：嫌窄/宽就改这里的空格数
  dt[, ` ` := paste(rep(" ", 60), collapse = "")]
  
  data_rows <- data.frame(
    Cohort                 = dt$cohort,
    `Binary fusion status` = dt$`Fusion+/Fusion-`,
    ` `                    = dt[[" "]],
    `HR (95% CI)`          = dt$`HR (95% CI)`,
    p                      = dt$p_lab,
    `q (BH)`               = dt$q_lab,
    check.names            = FALSE, stringsAsFactors = FALSE
  )
  
  # 伪表头行：直接复制 data_rows 的结构，保证列名/列序完全一致
  header_row <- data_rows[1, , drop = FALSE]
  header_row[] <- ""
  header_row[[2]] <- "Fusion+/Fusion-"
  header_row[[3]] <- dt[[" "]][1]
  plot_dt <- rbind(header_row, data_rows)   
  # 4. 轴：以 1 为中心左右对称（对数轴下 1/x_max 与 x_max 关于 1 对称）
  use_log <- all(c(dt$HR, dt$low, dt$up) > 0, na.rm = TRUE)
  X_HI <- x_max
  X_LO <- if (use_log) 1 / X_HI else -X_HI
  xlim <- c(X_LO, X_HI)
  
  ticks <- if (use_log) c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
  if (!is.null(ticks)) {
    ticks <- ticks[ticks >= X_LO & ticks <= X_HI]
    if (length(ticks) < 3) ticks <- NULL
  }
  
  # 关键：不裁剪数据。forestploter 会把超出 xlim 的 CI 画到轴端并自动加箭头
  est_vec   <- c(NA, dt$HR)
  lower_vec <- c(NA, dt$low)
  upper_vec <- c(NA, dt$up)
  
  n_clip <- sum(dt$low < X_LO | dt$up > X_HI, na.rm = TRUE)
  
  foot <- paste0(
    "\n\n\n\n\n\n",                       # 空行：把脚注推到横轴数字与 xlab 之下
    " ICB PFS Cox: Surv ~ Binary fusion status, cohort as strata.\n",
    "q: Benjamini\u2013Hochberg FDR.\n"
  )
  
  wrap_lines <- function(txt, width = 110) {
    lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
    lines <- unlist(lapply(lines, function(x) {
      if (nchar(x) == 0) return("")
      strwrap(x, width = width)
    }))
    paste(lines, collapse = "\n")
  }
  footnote_txt <- wrap_lines(foot, width = 110)
  
  # 5. 主题：斑马纹 + 左对齐 + 黑色方块 + 显示横轴刻度（原来 xaxis_gp col = NA 把刻度藏了）
  if (is.null(theme)) {
    row_fill <- rep(c("white", "#E8F3E8"), length.out = nrow(plot_dt))
    theme <- forestploter::forest_theme(
      base_size    = 11,
      base_family  = "sans",
      ci_pch       = 15,
      ci_col       = "black",
      ci_fill      = "black",
      ci_lwd       = 1.4,
      ci_Theight   = 0.15,
      refline_gp   = grid::gpar(lwd = 1, lty = "dashed", col = "grey25"),
      vertline_lwd = 0.5,
      vertline_lty = "solid",
      vertline_col = "grey90",
      xaxis_gp     = grid::gpar(fontsize = 7.5, col = "grey30"),
      xlab_gp      = grid::gpar(fontsize = 9, col = "grey30"),
      title_just   = "left",
      core = list(
        bg_params = list(fill = row_fill, col = NA, lwd = 0),
        fg_params = list(hjust = 0, x = 0.02, fontsize = 10)
      ),
      colhead = list(
        fg_params = list(hjust = 0, x = 0.02, fontsize = 10, fontface = "bold")
      )
    )
  }
  
  vlines <- if (use_log) c(0.1, 0.5, 2, 5, 10) else NULL
  if (!is.null(vlines)) vlines <- vlines[vlines > xlim[1] & vlines < xlim[2]]
  
  p <- forestploter::forest(
    plot_dt,
    est = est_vec, lower = lower_vec, upper = upper_vec,
    ci_column = 3,
    ref_line = 1,
    vert_line = vlines,
    xlim = xlim,
    x_trans = if (use_log) "log" else "none",
    ticks_at = ticks,
    xlab = if (use_log) "Hazard ratio (log scale)" else "Hazard ratio",
    title = title,
    footnote = footnote_txt,
    theme = theme
  )
  
  p <- forestploter::edit_plot(p, row = 1, col = 2, which = "text",
                               gp = grid::gpar(fontface = "italic", cex = 0.9))
  
  n_sum <- sum(dt$is_summary)
  if (n_sum > 0 && n_sum < nrow(dt)) {
    p <- forestploter::add_border(p, part = "body", row = n_sum + 2,
                                  gp = grid::gpar(lwd = 0.8, col = "grey40"))
  }
  
  # 6. 尺寸：按实际内容量出来，不再用行数估，避免大片留白
  wh <- forestploter::get_wh(p)
  w  <- if (is.null(width)) as.numeric(wh["width"]) + 0.20 else width
  h  <- as.numeric(wh["height"]) + 0.35
  
  grDevices::pdf(out_pdf, width = w, height = h, onefile = FALSE)
  print(p)
  grDevices::dev.off()
  
  invisible(list(plot = p, width = w, height = h, n_clipped = n_clip))
}
make_forest_hr(result_all, 
               title = "ICB PFS: Fusion+ vs Fusion-", 
               out_pdf = out_path)


#supFigure2C
cohort_Response <- cohort |> filter(Response == TRUE)
cohort_sum_Response <- unique(cohort_Response$id)
fusion_subset  <- fread(file_meta) |>
  as.data.table() |>
  filter(sample_type %in% c("Tumor", "Tumor_Blood")) |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_sum_Response) |>
  select(Sample_ID, gene1, gene2, cohort) |>
  rename(Run = Sample_ID) |>
  distinct()
result_all <- data.table()
Response_all <- data.table()
for(cohort_id in cohort_sum_Response){
  path <- paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/", cohort_id, "_info.tsv")
  if(!file.exists(path)) next
  Need <- fread(path)
  if(filter(cohort_Response, id == cohort_id)$time == "Pre/On" ){
    if(cohort_id == "PRJNA356761"){
      Need1 <- Need |> filter(grepl("Pre", Sample_ID))
    }else if(cohort_id == "EGAD00001006282"){
      Need1 <- Need |> filter(grepl("_SCREEN", Sample_ID))
    }else{
      Need1 <- Need |> filter(grepl("Pre", Response2))
    }
  }else if(cohort_id %in% c("PRJNA482620", "PRJNA744780")){
    Need1 <- Need |> filter(grepl("Pre", Response2))
  }else if(cohort_id == "PHS001464"){
    Need1 <- Need |> mutate(
      Sample_type = case_when(
        TUMOR_COLLECTION_TIMEPOINT ==  "Acquired Resistance" ~ "On",
        TUMOR_COLLECTION_TIMEPOINT ==  "Pre-Immunotherapy"~ "Pre",
        TUMOR_COLLECTION_TIMEPOINT ==  "Off Therapy Recurrence"~ "On",
        TUMOR_COLLECTION_TIMEPOINT ==  "Diagnosis"~ "Pre",
        TUMOR_COLLECTION_TIMEPOINT ==  "Mixed Response (Site of Primary Resistance)"~ "On",
        TRUE ~ NA_character_
      ) )
  }else{
    Need1 <- Need
  }
  columns_to_select <- c("Patient_ID", "Response2","Treatment", "Run")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |> select(all_of(existing_columns))
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  if (length(missing_columns) > 0) {
    for (col in missing_columns) { Need2[[col]] <- NA }
  }
  cancer_type   <- cohort$cancer_type[cohort$id == cohort_id]
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort <- cohort_id
  Need2$cancer <- cancer_type
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
  Response_all <- rbind(Need3, Response_all) |> filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib"))
}

fusion_data_subset<- Response_all  |>
  mutate(event = case_when(
    Run %in% fusion_subset$Run  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                           TRUE ~ event))  |>
  group_by(cohort) |>
  filter(n_distinct(event) == 2) |>
  ungroup()  |>
  select(Patient_ID, event, Response2, cohort,cancer) |>
  distinct() |>
  rename(Response = Response2) |> 
  filter(Response != "" & !is.na(Response ) & Response != "not collected" ) |> 
  as.data.table()
fusion_data_subset$Response  <- gsub("Pre-", "", fusion_data_subset$Response )
fusion_data_subset$Response <- factor(fusion_data_subset$Response, levels = c("NR", "R"))
cancer_Response <- fusion_data_subset$cancer |> unique()
result_all <- data.table()
for (cohort_choose in cancer_Response) {
  data_choose <- fusion_data_subset|> filter(cancer == cohort_choose)
  if (length(unique(data_choose$event)) < 2 ||
      length(unique(data_choose$Response)) < 2) {
    warning(paste("跳过", cohort_choose, "：event 或 Response 只有一个类别"))
    next
  }
  contingency_table <- table(data_choose$event, data_choose$Response)
  data_choose$event    <- factor(data_choose$event,    levels = c("Fusion+", "Fusion-"))
  tab <- table(data_choose$event, data_choose$Response)
  stopifnot(all(dim(tab) == c(2, 2)))
  
  a  <- tab["Fusion+", "R"];  b  <- tab["Fusion+", "NR"]
  cc <- tab["Fusion-", "R"];  dd <- tab["Fusion-", "NR"]
  
  ## ---- 2. Haldane-Anscombe OR + 95% CI + log2 OR ------------------------
  fisher_or <- function(a, b, c, d) {
    a <- a + .5; b <- b + .5; c <- c + .5; d <- d + .5        # 四格各 +0.5
    or <- (a * d) / (b * c)
    se <- sqrt(1/a + 1/b + 1/c + 1/d)
    c(OR = or, lo = exp(log(or) - 1.96 * se), hi = exp(log(or) + 1.96 * se),
      log2OR = log2(or))
  }
  
  o     <- fisher_or(a, b, cc, dd)
  p_val <- fisher.test(tab)$p.value
  
  result_choose <- data.frame(
    cohort = cohort_choose,
    hn = tab["Fusion+", "R"] + tab["Fusion+", "NR"],
    ln = tab["Fusion-", "R"] + tab["Fusion-", "NR"],
    OR    = o[["OR"]],  
    low95 = o[["lo"]],
    up95  = o[["hi"]],
    p     = p_val
  )
  result_all <- rbind(result_all, result_choose)
}
result_all <- result_all |> 
  filter(!is.infinite(OR) & low95 != 0 & up95 != 0)|>
  filter(hn >= 3 & ln >= 3) |>
  mutate(q = p.adjust(p, method = "BH"))|>
  arrange(p)
  
make_forest_or <- function(dt, title, out_pdf, max_rows = 80, theme = NULL,
                           x_max = 20, width = NULL) {
  
  if (nrow(dt) == 0) { message("empty: ", title); return(invisible(NULL)) }
  dt <- data.table::copy(data.table::as.data.table(dt))
  
  # 1. 兼容不同版本的输入列名
  if ("cohort_name" %in% names(dt) && !"cohort" %in% names(dt)) data.table::setnames(dt, "cohort_name", "cohort")
  if ("hn" %in% names(dt) && !"n_pos" %in% names(dt)) data.table::setnames(dt, "hn", "n_pos")
  if ("ln" %in% names(dt) && !"n_neg" %in% names(dt)) data.table::setnames(dt, "ln", "n_neg")
  if ("low95" %in% names(dt) && !"low" %in% names(dt)) data.table::setnames(dt, "low95", "low")
  if ("up95" %in% names(dt) && !"up" %in% names(dt)) data.table::setnames(dt, "up95", "up")
  if ("q(BH)" %in% names(dt) && !"q" %in% names(dt)) data.table::setnames(dt, "q(BH)", "q")
  if (!"is_summary" %in% names(dt)) dt[, is_summary := FALSE]
  if (!"q" %in% names(dt)) dt[, q := NA_real_]
  
  # 2. 排序并截断
  dt[, row_order := data.table::fifelse(is_summary, 0L, 1L)]
  data.table::setorder(dt, row_order)
  dt <- dt[1:min(.N, max_rows)]
  
  # 3. 格式化显示的文本（OR 列保持原样，写真实 CI）
  dt[, `Fusion+/Fusion-` := paste0(n_pos, "/", n_neg)]
  dt[, `OR (95% CI)` := sprintf("%.2f (%.2f\u2013%.2f)", OR, low, up)]
  dt[, p_lab := ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))]
  dt[, q_lab := ifelse(is.na(q), "", ifelse(q < 0.001, "<0.001", sprintf("%.3f", q)))]
  
  # 森林图列宽度：嫌窄/宽就改这里的空格数
  dt[, ` ` := paste(rep(" ", 50), collapse = "")]
  
  data_rows <- data.frame(
    Cohort        = dt$cohort,
    `Binary fusion status`        = dt$`Fusion+/Fusion-`,
    ` `           = dt[[" "]],
    `OR (95% CI)` = dt$`OR (95% CI)`,
    p             = dt$p_lab,
    `q (BH)`      = dt$q_lab,
    check.names   = FALSE, stringsAsFactors = FALSE
  )
  
  # 伪表头行：复制 data_rows 的结构，列名/列序必然一致
  header_row <- data_rows[1, , drop = FALSE]
  header_row[] <- ""
  header_row[[2]] <- "Fusion+/Fusion-"
  header_row[[3]] <- dt[[" "]][1]
  
  plot_dt <- rbind(header_row, data_rows)
  
  # 4. 轴：以 1 为中心左右对称（对数轴下 1/x_max 与 x_max 关于 1 对称）
  use_log <- all(c(dt$OR, dt$low, dt$up) > 0, na.rm = TRUE)
  X_HI <- x_max
  X_LO <- if (use_log) 1 / X_HI else -X_HI
  xlim <- c(X_LO, X_HI)
  
  ticks <- if (use_log) c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
  if (!is.null(ticks)) {
    ticks <- ticks[ticks >= X_LO & ticks <= X_HI]
    if (length(ticks) < 3) ticks <- NULL
  }
  
  # 关键：不裁剪数据。forestploter 会把超出 xlim 的 CI 画到轴端并自动加箭头
  est_vec   <- c(NA, dt$OR)
  lower_vec <- c(NA, dt$low)
  upper_vec <- c(NA, dt$up)
  
  n_clip <- sum(dt$low < X_LO | dt$up > X_HI, na.rm = TRUE)
  
  # 脚注文字与原来一致，只是多了前面的空行（显示横轴后必须下移，否则压在刻度上）
  footnote_txt <- paste(
    "\n\n\n\n\n\n",
    "ICB: Binary fusion status ~ Response.\n",
    "SGC and GBM were not analysed because the fusion-negative group was too small (n = 2 in each)\n",
    "q: BH-FDR.\n"
  )
  
  wrap_lines <- function(txt, width = 110) {
    lines <- strsplit(txt, "\n", fixed = TRUE)[[1]]
    lines <- unlist(lapply(lines, function(x) {
      if (nchar(x) == 0) return("")
      strwrap(x, width = width)
    }))
    paste(lines, collapse = "\n")
  }
  footnote_txt <- wrap_lines(footnote_txt, width = 110)
  
  # 5. 主题：斑马纹 + 左对齐 + 黑色方块 + 显示横轴刻度
  if (is.null(theme)) {
    row_fill <- rep(c("white", "#E8F3E8"), length.out = nrow(plot_dt))
    theme <- forestploter::forest_theme(
      base_size    = 11,
      base_family  = "sans",
      ci_pch       = 15,
      ci_col       = "black",
      ci_fill      = "black",
      ci_lwd       = 1.4,
      ci_Theight   = 0.15,
      refline_gp   = grid::gpar(lwd = 1, lty = "dashed", col = "grey25"),
      vertline_lwd = 0.5,
      vertline_lty = "solid",
      vertline_col = "grey90",
      xaxis_gp     = grid::gpar(fontsize = 7.5, col = "grey30"),
      xlab_gp      = grid::gpar(fontsize = 9, col = "grey30"),
      title_just   = "left",
      core = list(
        bg_params = list(fill = row_fill, col = NA, lwd = 0),
        fg_params = list(hjust = 0, x = 0.02, fontsize = 10)
      ),
      colhead = list(
        fg_params = list(hjust = 0, x = 0.02, fontsize = 10, fontface = "bold")
      )
    )
  }
  
  vlines <- if (use_log) c(0.1, 0.5, 2, 5, 10) else NULL
  if (!is.null(vlines)) vlines <- vlines[vlines > xlim[1] & vlines < xlim[2]]
  
  p <- forestploter::forest(
    plot_dt,
    est = est_vec, lower = lower_vec, upper = upper_vec,
    ci_column = 3,
    ref_line = 1,
    vert_line = vlines,
    xlim = xlim,
    x_trans = if (use_log) "log" else "none",
    ticks_at = ticks,
    xlab = if (use_log) "Odds ratio (log scale)" else "Odds ratio",
    title = title,
    footnote = footnote_txt,
    theme = theme
  )
  
  p <- forestploter::edit_plot(p, row = 1, col = 2, which = "text",
                               gp = grid::gpar(fontface = "italic", cex = 0.9))
  
  n_sum <- sum(dt$is_summary)
  if (n_sum > 0 && n_sum < nrow(dt)) {
    p <- forestploter::add_border(p, part = "body", row = n_sum + 2,
                                  gp = grid::gpar(lwd = 0.8, col = "grey40"))
  }
  
  # 6. 尺寸：按实际内容量出来，避免大片留白
  wh <- forestploter::get_wh(p)
  w  <- if (is.null(width)) as.numeric(wh["width"]) + 0.20 else width
  h  <- as.numeric(wh["height"]) + 0.35
  
  grDevices::pdf(out_pdf, width = w, height = h, onefile = FALSE)
  print(p)
  grDevices::dev.off()
  
  invisible(list(plot = p, width = w, height = h, n_clipped = n_clip))
}
make_forest_or(result_all, 
               title = "ICB Response: Fusion+ vs Fusion-", 
               out_pdf = "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图2/supFigure2C.pdf")

#supFigure2D
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


message("Cohorts sig(q<0.05): ",
        res_ok[is_summary == FALSE & sig == TRUE, .N],
        " | summary sig: ",
        res_ok[is_summary == TRUE & sig == TRUE, .N],
        " | summary rows present: ",
        res_ok[is_summary == TRUE, .N],
        " | Inf/non-estimable: ", nrow(res_note))
print(res_ok[is_summary == TRUE, .(cohort, n_pos, n_neg, HR, low, up, p, q, sig)])

# ---------------- 8. make_forest：footnote 换行 + 高度 ----------------

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
  
  ext <- tolower(tools::file_ext(out_pdf))
  if (ext == "png") {
    png(out_pdf, width = 10, height = h, units = "in", res = 300)
  } else {
    pdf(out_pdf, width = 10, height = h, onefile = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  
  print(p)
  invisible(out_pdf)
}
data <- res_ok[order(q)] 
make_forest(data,
            "GDC OS",
            paste0(out_dir, "supFigure2D.pdf"))


#supFigure2E
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

data <- res_ok[order(q)]
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
  
  ext <- tolower(tools::file_ext(out_pdf))
  if (ext == "png") {
    png(out_pdf, width = 10, height = h, units = "in", res = 300)
  } else {
    pdf(out_pdf, width = 10, height = h, onefile = FALSE)
  }
  on.exit(dev.off(), add = TRUE)
  
  print(p)
  invisible(out_pdf)
}
make_forest(data,
            "GDC PFS",
            paste0(out_dir, "supFigure2E.pdf"))
