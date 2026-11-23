## ---- 0. 依赖 ----------------------------------------------------------------
need_pkgs <- c("data.table", "dplyr", "tidyr", "purrr", "tibble",
               "ggplot2", "scales", "stringr", "survival", "cowplot")
for (p in need_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(need_pkgs, library, character.only = TRUE))
set.seed(2024)

## ---- 1. 路径与开关 ----------------------------------------------------------
file_meta <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
meta_file <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/cohort_meta_check.csv"
clin_dir  <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/"
out_dir   <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图3/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
## ---- 2. 通用统计函数 --------------------------------------------------------
fisher_or <- function(a, b, c, d) {
  a_c <- a + .5; b_c <- b + .5; c_c <- c + .5; d_c <- d + .5        
  or <- (a_c * d_c) / (b_c * c_c)
  se <- sqrt(1/a_c + 1/b_c + 1/c_c + 1/d_c)
  c(OR = or, lo = exp(log(or) - 1.96 * se), hi = exp(log(or) + 1.96 * se))
}

## ---- 3. 读入融合表，展开成 (Run, fusion) 长表 --------------------------------
fusion_raw <- fread(file_meta, sep = "\t")
cohort1 <- fread(meta_file)
CANCERS_SEL    <- c("NSCLC", "SKCM", "BLCA", "KIRC", "STAD", "SGC")  
BOX_FROM_FOREST <- TRUE              
EXCLUDE_CANCER <- c("BRCA")          
MIN_N          <- 3              
LOG_Y_BOX      <- FALSE              
USE_FDR_COLOUR <- TRUE               
FDR_METHOD     <- "BH"               

COL_R  <- "#9dd3af"                  
COL_NR <- "#f2f2f2"                  
COL_HI <- "#c0504d"

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))
}

pick_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x)) x[1] else NA_real_
}

SURV_UNIT   <- "day"                 
SURV_DIV    <- if (SURV_UNIT == "day") 30.4375 else 1
SURV_COLS   <- c("OS_Time", "OS_Status", "PFS_Time", "PFS_Status")

if ("OS" %in% names(cohort1))  cohort1[, OS_flag  := as.logical(OS)]  else cohort1[, OS_flag  := TRUE]
if ("PFS" %in% names(cohort1)) cohort1[, PFS_flag := as.logical(PFS)] else cohort1[, PFS_flag := TRUE]
cohort_flag <- cohort1 |> select(cohort = id, OS_flag, PFS_flag) |> distinct()

cox_by_cancer <- function(d, endpoint) {
  purrr::map_dfr(sort(unique(d$cancer)), function(cc) {
    dd <- d[d$cancer == cc, ] |> as.data.frame()
    dd$cohort <- factor(dd$cohort)
    n_ch <- length(unique(dd$cohort))
    n  <- nrow(dd)
    dat <- as.data.frame(dd)
    dat$cohort <- factor(dat$cohort)                       # ← 因子化，droplevels 才能用
    if (nlevels(droplevels(dat$cohort)) > 1){
        if(endpoint == "OS"){
          fit <- survfit(Surv(OS_Time, OS_Status) ~ TFB_group, data = dat)
          formula_cox <- as.formula(paste0("Surv(OS_Time, OS_Status) ~ TFB_group + strata(cohort)"))
        }else{
          fit <- survfit(Surv(PFS_Time, PFS_Status) ~ TFB_group, data = dat)
          formula_cox <- as.formula(paste0("Surv(PFS_Time, PFS_Status) ~ TFB_group + strata(cohort)"))
        }
      }else{
        if(endpoint == "OS"){
          fit <- survfit(Surv(OS_Time, OS_Status) ~ TFB_group, data = dat)
          formula_cox <- as.formula(paste0("Surv(OS_Time, OS_Status) ~ TFB_group"))
        }else{
          fit <- survfit(Surv(PFS_Time, PFS_Status) ~ TFB_group, data = dat)
          formula_cox <- as.formula(paste0("Surv(PFS_Time, PFS_Status) ~ TFB_group"))
        }
      }
    cox_model <- coxph(formula_cox, data = dd)
    cox_summary <- summary(cox_model)
    s1 <- cox_summary$conf.int

    
    tibble(cancer = cc, endpoint = endpoint, n = n, 
           n_cohort = n_ch,
           p =  cox_summary$coefficients[1, "Pr(>|z|)"],
           HR    = s1[1, "exp(coef)"],
           low95 = s1[1, "lower .95"],
           up95  = s1[1, "upper .95"])
  })
}


#OS/PFS生存森林图
cohort_OS <- cohort1 |> filter(OS == TRUE | PFS == TRUE) 
cohort_sum_OS <- (cohort_OS$id) |> unique()

Tumor_ICB <- fusion_raw |>
  filter(sample_type == "Tumor") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_sum_OS) |>
  select(Sample_ID, Fusion_ID, cohort) |>
  rename(Run = Sample_ID) |>
  distinct()
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
  c_cancer <- cohort1 |> filter(id == cohort_id ) |> pull(cancer_type )
  columns_to_select <- c("Patient_ID", "OS_Time", "OS_Status", "PFS_Time", "PFS_Status","Run", "Treatment", "Age", "Sex")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    select(all_of(existing_columns))
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  if (length(missing_columns) > 0) {
    for (col in missing_columns) {
      Need2[[col]] <- NA
    }
  }
  
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort <- cohort_id
  Need2$cancer <- c_cancer
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
  
  OS_all <- rbind(Need3,OS_all) |>   filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
}
data_all <- OS_all|>
  left_join(Tumor_ICB, by = c("cohort", "Run")) |>
  select(cohort,  Patient_ID, Fusion_ID, cancer) |>
  distinct()|>
  group_by(cohort, Patient_ID) |>
  summarise(TFB = n_distinct(Fusion_ID[!is.na(Fusion_ID)])) |>
  ungroup()

surv_info <- OS_all |>
            select(cancer,
                   Patient_ID, cohort,
                   OS_Time, OS_Status,
                   PFS_Time,PFS_Status) |>
           distinct()

plot_forest_surv <- function(tab, endpoint) {
  d <- tab |>
    filter(endpoint == !!endpoint, !is.na(HR)) |>
    arrange(HR) |>
    mutate(y     = row_number(),
           label = sprintf("%s (n=%d, %d cohort%s)", cancer, n,
                           n_cohort, ifelse(n_cohort > 1, "s", "")),
           right = sprintf("%.2f (%.2f-%.2f)\np =%s, q (BH) =%s",
                           HR, low95, up95, fmt_p(p), fmt_p(p_adj)),
           hi_col = if (USE_FDR_COLOUR) !is.na(p_adj) & p_adj < 0.05 else p < 0.05)
  if (!nrow(d)) return(NULL)
  
  ## ---- 轴范围：右侧固定在 20，左侧取 1/20，log 尺度下以 HR = 1 为中心左右对称 ----
  x_hi <- 20
  x_lo <- 1 / x_hi            # = 0.05
  
  d <- d |>
    mutate(lo   = pmax(low95, x_lo),
           hi   = pmin(up95,  x_hi),
           c_lo = low95 < x_lo,
           c_hi = up95  > x_hi,
           HR_p = pmin(pmax(HR, x_lo), x_hi),   # 点在轴外时贴边，避免画到留白里
           label = factor(label, levels = label))
  
  ## 对称轴下的整齐刻度（只取落在 0.05-20 内的）
  brk <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
  brk <- brk[brk >= x_lo & brk <= x_hi]
  
  ggplot(d, aes(y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(xmin = lo, xmax = hi, colour = hi_col),
                  orientation = "y", width = .18, linewidth = .55) +
    geom_segment(data = dplyr::filter(d, c_hi),
                 aes(x = x_hi, xend = x_hi * 1.35, y = label, yend = label,
                     colour = hi_col),
                 arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
                 linewidth = .55, show.legend = FALSE) +
    geom_segment(data = dplyr::filter(d, c_lo),
                 aes(x = x_lo, xend = x_lo * 0.72, y = label, yend = label,
                     colour = hi_col),
                 arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
                 linewidth = .55, show.legend = FALSE) +
    geom_point(aes(x = HR_p, colour = hi_col, size = n), shape = 15) +
    geom_text(aes(x = x_hi * 1.55, label = right), hjust = 0, size = 2.9,
              lineheight = .9, colour = "grey15") +
    scale_x_log10(breaks = brk,
                  labels = scales::label_number(accuracy = 0.01,
                                                drop0trailing = TRUE)) +
    scale_colour_manual(values = c(`TRUE` = COL_HI, `FALSE` = "grey35"),
                        guide = "none") +
    scale_size_continuous(range = c(1.4, 4.2), guide = "none") +
    annotation_logticks(sides = "b", colour = "grey55", linewidth = .3) +
    labs(x = sprintf("%s: hazard ratio, High vs Low TFB (log scale)", endpoint),
         y = NULL,
         title = sprintf("ICB-%s-TFB", endpoint),
         caption = sprintf(paste0("Cox model with cohort as stratum; HR < 1 = better %s.\n",
                                  "Bars: 95%% CI; arrows: CI beyond axis. q (%s), %d cancers."),
                           endpoint, FDR_METHOD, nrow(d))) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 12, margin = margin(r = 20)),
      plot.title = element_text(size = 12, face = "bold", hjust = 0),
      plot.caption = element_text(size = 11, colour = "grey40"),
      plot.margin = margin(8, 235, 8, 8)) +
    coord_cartesian(xlim = c(x_lo, x_hi), clip = "off")
}

for (ep in c("OS", "PFS")) {
  if(ep == "OS"){
    
    data_tab  <- data_all |>
      left_join(surv_info, by = c("cohort", "Patient_ID")) |> 
      filter(!is.na(OS_Time), !is.na(OS_Status)) |>
      filter(!cancer %in% EXCLUDE_CANCER)
    
    thr_tab <- data_tab |> 
      group_by(cancer) |>
      summarise(n = n(), thr = floor(median(TFB)), .groups = "drop")
    
    patient_tfb <- data_tab  |>
      left_join(thr_tab |> select(cancer, thr), by = "cancer") |>
      mutate(TFB_group = factor(ifelse(TFB > thr, "High TFB", "Low TFB"), levels = c("Low TFB", "High TFB")))|>
      group_by(cancer)
    surv_tab  <- cox_by_cancer(patient_tfb, ep) |>
                 group_by(endpoint) |>
                 mutate(p_adj    = p.adjust(p))
    figure_name = "supFigure3A"
    p_ep <- plot_forest_surv(surv_tab, ep)
    }else{
      data_tab  <- data_all |>
        left_join(surv_info, by = c("cohort", "Patient_ID")) |> 
        filter(!is.na(PFS_Time), !is.na(PFS_Status)) |>
        filter(!cancer %in% EXCLUDE_CANCER)
      
      thr_tab <- data_tab |> 
        group_by(cancer) |>
        summarise(n = n(), thr = floor(median(TFB)), .groups = "drop")
      
      patient_tfb <- data_tab  |>
        left_join(thr_tab |> select(cancer, thr), by = "cancer") |>
        mutate(TFB_group = factor(ifelse(TFB > thr, "High TFB", "Low TFB"), levels = c("Low TFB", "High TFB")))|>
        group_by(cancer)
      surv_tab  <- cox_by_cancer(patient_tfb, ep) |>
        group_by(endpoint) |>
        mutate(p_adj    = p.adjust(p))
      
    figure_name = "supFigure3B"
    p_ep <- plot_forest_surv(surv_tab, ep) 
    }
  
  if (is.null(p_ep)) { message("无 ", ep, " 数据，跳过森林图"); next }
  nrow_ep <- sum(surv_tab$endpoint == ep & !is.na(surv_tab$HR))
  suppressWarnings({
    ggsave(file.path(out_dir, paste0(figure_name, ".pdf")), p_ep,
           width = 10, height = 0.42 * nrow_ep + 2.2)
    ggsave(file.path(out_dir,paste0(figure_name,  ".png")), p_ep,
           width = 10, height = 0.42 * nrow_ep + 2.2, dpi = 600)
  })
}


#supFigure3C,D
if (!"Response2" %in% names(cohort1) && "Response" %in% names(cohort1))
  cohort1[, Response2 := Response]
if (!"cancer" %in% names(cohort1) && "cancer_type" %in% names(cohort1))
  cohort1[, cancer := cancer_type]
if (!"time" %in% names(cohort1)) cohort1[, time := NA_character_]

cohort_resp <- cohort1 |> filter(Response2 == TRUE)
cohort_ids  <- unique(cohort_resp$id)

Tumor_ICB <- fusion_raw |>
  filter(sample_type == "Tumor") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_ids) |>
  select(Sample_ID, Fusion_ID, cohort) |>
  select(Run = Sample_ID, cohort, Fusion_ID) |>
  distinct()

## ---- 4. 汇总临床信息 --------------------------------------------------------
Response_all <- data.table()
for (cohort_id in cohort_ids) {
  need_file <- paste0(clin_dir, cohort_id, "_info.tsv")
  if (!file.exists(need_file)) next
  Need <- fread(need_file)
  
  c_time   <- cohort1[id == cohort_id, time] |> unique()
  c_cancer <- cohort1[id == cohort_id, cancer] |> unique() |> first()
  
  if ("Pre/On" %in% c_time) {
    Need1 <- if (cohort_id == "PRJNA356761" && "Sample_ID" %in% names(Need)) {
      Need |> filter(grepl("Pre", Sample_ID))
    } else if (cohort_id == "EGAD00001006282" && "Sample_ID" %in% names(Need)) {
      Need |> filter(grepl("_SCREEN", Sample_ID))
    } else if ("Response2" %in% names(Need)) {
      Need |> filter(grepl("Pre", Response2))
    } else Need
  } else if (cohort_id %in% c("PRJNA482620", "PRJNA744780") && "Response2" %in% names(Need)) {
    Need1 <- Need |> filter(grepl("Pre", Response2))
  } else if (cohort_id == "PHS001464" && "TUMOR_COLLECTION_TIMEPOINT" %in% names(Need)) {
    Need1 <- Need |>
      mutate(Sample_type = case_when(
        TUMOR_COLLECTION_TIMEPOINT == "Acquired Resistance"                       ~ "On",
        TUMOR_COLLECTION_TIMEPOINT == "Pre-Immunotherapy"                         ~ "Pre",
        TUMOR_COLLECTION_TIMEPOINT == "Off Therapy Recurrence"                    ~ "On",
        TUMOR_COLLECTION_TIMEPOINT == "Diagnosis"                                 ~ "Pre",
        TUMOR_COLLECTION_TIMEPOINT == "Mixed Response (Site of Primary Resistance)" ~ "On",
        TRUE ~ NA_character_)) |>
      filter(Sample_type == "Pre")
  } else {
    Need1 <- Need
  }
  
  cols_keep <- c("Patient_ID", "Run", "Response2", "Treatment", intersect(SURV_COLS, names(Need1)))
  cols_keep <- cols_keep[cols_keep %in% names(Need1)]
  Need2 <- Need1 |> select(all_of(cols_keep))
  for (col in setdiff(c("Patient_ID", "Run", "Response2", "Treatment", SURV_COLS), names(Need2)))
    Need2[[col]] <- NA
  
  Need2 <- Need2 |>
    mutate(Patient_ID = as.character(Patient_ID), cancer = c_cancer, cohort = cohort_id) |>
    filter(!Treatment %in% c("Docetaxel", "Observation", "Sunitinib"))
  
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
  
  Response_all <- rbind(Need3, Response_all, fill = TRUE)
}

## ---- 5. 患者层 TFB 与分组 ---------------------------------------------------
Response_all1 <- Response_all |>
  filter(!is.na(Response2), Response2 != "", !grepl("On", Response2)) |>
  mutate(Response = gsub("^Pre-", "", Response2), Response = trimws(Response)) |>
  filter(Response %in% c("R", "NR")) |>
  left_join(Tumor_ICB, by = c("cohort" = "cohort", "Run" = "Run")) |>
  group_by(cohort, Patient_ID) |>
  summarise(cancer = first(cancer), Response = first(Response),
            TFB = n_distinct(Fusion_ID[!is.na(Fusion_ID)]), .groups = "drop")

patient_tfb <- Response_all1 |> ungroup() |> filter(!cancer %in% EXCLUDE_CANCER)

thr_tab <- patient_tfb |> group_by(cancer) |>
  summarise(n = n(), thr = floor(median(TFB)), .groups = "drop")

patient_tfb <- patient_tfb |>
  left_join(thr_tab |> select(cancer, thr), by = "cancer") |>
  mutate(TFB_group = factor(ifelse(TFB > thr, "High TFB", "Low TFB"), levels = c("Low TFB", "High TFB")),
         Response  = factor(Response, levels = c("R", "NR")))

## ---- 5b. 计算森林图统计量 --------------------
forest_tab <- patient_tfb |>
  filter(!is.na(Response), !is.na(TFB_group)) |>
  group_by(cancer) |>
  group_modify(~ {
    tab <- table(factor(.x$TFB_group, levels = c("Low TFB", "High TFB")),
                 factor(.x$Response,  levels = c("NR", "R")))
    n <- nrow(.x)
    if (!all(dim(tab) == 2) || any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
      return(tibble(n = n, n_low = NA_integer_, n_high = NA_integer_,
                    R_low = NA_integer_, R_high = NA_integer_,
                    OR = NA_real_, low95 = NA_real_, up95 = NA_real_, p = NA_real_))
    }
    
    a <- tab["High TFB", "R"]; b <- tab["High TFB", "NR"]
    c <- tab["Low TFB", "R"];  d <- tab["Low TFB", "NR"]
    
    o <- fisher_or(a, b, c, d)
    p_val <- fisher.test(tab)$p.value
    
    tibble(n = n, n_low = c + d, n_high = a + b, R_low = c, R_high = a,
           OR = o[["OR"]], low95 = o[["lo"]], up95 = o[["hi"]], p = p_val)
  }) |>
  ungroup() |>
  filter(!is.na(OR), n >= MIN_N, is.finite(low95), is.finite(up95), low95 > 0, up95 > 0) |>
  mutate(p_adj = p.adjust(p, method = FDR_METHOD), sig = p < 0.05, sig_fdr = p_adj < 0.05)


## --- supFigure4C,D  --------
plot_forest_surv <- function(tab, endpoint, hr_lab = "High vs Low TFB",
                             title_lab = NULL, axis_max = 20,
                             arrow_k = 0.05, text_k = 0.073) {
  ## 局部向量化 p 格式化：不依赖全局 fmt_p 是否写成 if/else 标量版
  fp <- function(x) {
    out <- ifelse(is.na(x), "NA", sprintf("%.3f", x))
    out[which(!is.na(x) & x < 0.001)] <- "< 0.001"
    out
  }
  
  d <- tab |> filter(endpoint == !!endpoint, !is.na(HR)) |> arrange(HR) |>
    mutate(y = row_number(),
           label = sprintf("%s (n=%d, %d cohort%s)", cancer, n, n_cohort,
                           ifelse(n_cohort > 1, "s", "")),
           right = sprintf("%.2f (%.2f-%.2f)\np =%s, q (BH) =%s",
                           HR, low95, up95, fp(p), fp(p_adj)),   # ← 这里换掉 fmt_p
           hi_col = if (USE_FDR_COLOUR) !is.na(p_adj) & p_adj < 0.05 else p < 0.05)
  if (!nrow(d)) return(NULL)
  
  ## axis_max = NULL：按数据自动取以 1 为中心、左右对称的范围（上限仍 ≤ 20）
  if (is.null(axis_max)) {
    need <- suppressWarnings(max(c(1 / d$low95, d$up95, 1), na.rm = TRUE)) * 1.05
    if (!is.finite(need)) need <- 20
    ladder <- c(1.1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20)
    axis_max <- ladder[which(ladder >= need)[1]]
    if (is.na(axis_max)) axis_max <- 20
    axis_max <- min(axis_max, 20)
  }
  x_hi <- axis_max; x_lo <- 1 / x_hi; ratio <- x_hi / x_lo
  
  d <- d |> mutate(lo = pmax(low95, x_lo), hi = pmin(up95, x_hi),
                   c_lo = low95 < x_lo, c_hi = up95 > x_hi,
                   HR_p = pmin(pmax(HR, x_lo), x_hi),
                   label = factor(label, levels = label))
  
  if (abs(x_hi - 20) < 1e-9) {
    brk <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20)
  } else {
    fam <- sort(unique(as.vector(outer(c(1, 2, 5), 10^(-2:1)))))
    fam <- fam[fam <= x_hi * 1.001]
    brk <- sort(unique(c(fam, 1 / fam)))          # 1-2-5 阶梯自身关于 1 对称
  }
  brk <- brk[brk >= x_lo * 0.999 & brk <= x_hi * 1.001]
  if (length(brk) < 3) brk <- c(x_lo, 1, x_hi)
  
  x_txt <- x_hi * ratio^text_k; x_arr <- x_hi * ratio^arrow_k
  x_arr_l <- x_lo / ratio^arrow_k
  
  ggplot(d, aes(y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(xmin = lo, xmax = hi, colour = hi_col),
                  orientation = "y", width = .18, linewidth = .55) +
    geom_segment(data = dplyr::filter(d, c_hi),
                 aes(x = x_hi, xend = x_arr, y = label, yend = label, colour = hi_col),
                 arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
                 linewidth = .55, show.legend = FALSE) +
    geom_segment(data = dplyr::filter(d, c_lo),
                 aes(x = x_lo, xend = x_arr_l, y = label, yend = label, colour = hi_col),
                 arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
                 linewidth = .55, show.legend = FALSE) +
    geom_point(aes(x = HR_p, colour = hi_col, size = n), shape = 15) +
    geom_text(aes(x = x_txt, label = right), hjust = 0, size = 2.9,
              lineheight = .9, colour = "grey15") +
    scale_x_log10(breaks = brk,
                  labels = scales::label_number(accuracy = 0.01, drop0trailing = TRUE)) +
    scale_colour_manual(values = c(`TRUE` = COL_HI, `FALSE` = "grey35"), guide = "none") +
    scale_size_continuous(range = c(1.4, 4.2), guide = "none") +
    annotation_logticks(sides = "b", colour = "grey55", linewidth = .3) +
    labs(x = sprintf("%s: hazard ratio, %s (log scale)", endpoint, hr_lab),
         y = NULL,
         title = if (is.null(title_lab)) sprintf("ICB-%s-TFB", endpoint) else title_lab,
         caption = sprintf(paste0("Cox model with cohort as stratum; HR < 1 = better %s.\n",
                                  "Bars: 95%% CI; arrows: CI beyond axis. q (%s), %d cancers."),
                           endpoint, FDR_METHOD, nrow(d))) +
    theme_classic(base_size = 10) +
    theme(axis.text.y = element_text(size = 12, margin = margin(r = 20)),
          plot.title = element_text(size = 12, face = "bold", hjust = 0),
          plot.caption = element_text(size = 11, colour = "grey40"),
          plot.margin = margin(8, 235, 8, 8)) +
    coord_cartesian(xlim = c(x_lo, x_hi), clip = "off")
}
cox_by_cancer_cont <- function(d, endpoint, unit = 1) {
  cols <- if (endpoint == "OS") c("OS_Time", "OS_Status") else c("PFS_Time", "PFS_Status")
  purrr::map_dfr(sort(unique(d$cancer)), function(cc) {
    dd <- as.data.frame(d[d$cancer == cc, c("cohort", "TFB", cols), drop = FALSE])
    dd <- dd[stats::complete.cases(dd), ]
    if (nrow(dd) < 3) return(NULL)                                    # 样本太少
    dd$cohort <- factor(dd$cohort)
    n_ch <- nlevels(droplevels(dd$cohort))
    if (length(unique(dd$TFB)) < 2) return(NULL)                      # TFB 无变异
    if (sum(dd[[cols[2]]] == 1, na.rm = TRUE) < 3) return(NULL)       # 事件太少
    dd$TFB_u <- dd$TFB / unit
    rhs <- sprintf("Surv(%s, %s) ~ TFB_u", cols[1], cols[2])
    if (n_ch > 1) rhs <- paste0(rhs, " + strata(cohort)")
    fit <- tryCatch(coxph(stats::as.formula(rhs), data = dd),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    s1 <- summary(fit)$conf.int
    tibble(cancer = cc, endpoint = endpoint, n = nrow(dd), n_cohort = n_ch,
           p     = summary(fit)$coefficients[1, "Pr(>|z|)"],
           HR    = s1[1, "exp(coef)"],
           low95 = s1[1, "lower .95"],
           up95  = s1[1, "upper .95"])
  })
}

## ---- 出图 ---------------------------------------------------------------
CONT_UNIT     <- 1            # 每 1 个融合；改成 sd(...) 即「每 1 个 SD」
CONT_ONLY_POS <- FALSE        # TRUE: 仅在 TFB > 0 的患者中做（去零膨胀）
CONT_AXIS_MAX <- NULL         # NULL = 自适应对称轴；20 = 固定 0.05-20
CONT_FIG      <- c(OS = "supFigure3C", PFS = "supFigure3D")

for (ep in c("OS", "PFS")) {
  d_tab <- data_all |>
    select(cohort, Patient_ID, TFB) |>
    left_join(surv_info, by = c("cohort", "Patient_ID")) |>
    filter(if (ep == "OS") !is.na(OS_Time) & !is.na(OS_Status)
           else !is.na(PFS_Time) & !is.na(PFS_Status)) |>
    filter(!is.na(TFB), TFB >= 0, !cancer %in% EXCLUDE_CANCER)
  if (CONT_ONLY_POS) d_tab <- filter(d_tab, TFB > 0)
  if (!nrow(d_tab)) { message("无 ", ep, " 数据，跳过"); next }
  
  surv_cont <- cox_by_cancer_cont(d_tab, ep, unit = CONT_UNIT) |>
    group_by(endpoint) |>
    mutate(p_adj = p.adjust(p, method = FDR_METHOD)) |>
    ungroup()
  
  p_ep <- plot_forest_surv(surv_cont, ep,
                           hr_lab = if (isTRUE(CONT_UNIT == 1))
                             "per 1 fusion (continuous)"
                           else sprintf("per %.2f fusions (continuous)", CONT_UNIT),
                           title_lab = sprintf("ICB-%s-TFB (continuous)", ep),
                           axis_max = CONT_AXIS_MAX)
  if (is.null(p_ep)) next
  
  nrow_ep <- sum(surv_cont$endpoint == ep & !is.na(surv_cont$HR))
  h <- 0.42 * nrow_ep + 2.2
  suppressWarnings({
    ggsave(file.path(out_dir, paste0(CONT_FIG[[ep]], ".pdf")), p_ep, width = 10, height = h)
    ggsave(file.path(out_dir, paste0(CONT_FIG[[ep]], ".png")), p_ep, width = 10, height = h, dpi = 600)
  })
}




## ---- supFigure3E-------------------------
if (BOX_FROM_FOREST) CANCERS_SEL <- forest_tab$cancer[order(forest_tab$OR)]

box_df <- patient_tfb |>
  filter(cancer %in% CANCERS_SEL, !is.na(Response))
box_df$cancer <- factor(box_df$cancer, levels = intersect(CANCERS_SEL, unique(box_df$cancer)))

p_lab <- purrr::map_dfr(levels(droplevels(box_df$cancer)), function(cc) {
  d <- box_df[box_df$cancer == cc, ]
  p <- NA_real_
  if (nrow(d) >= 3 && length(unique(d$Response)) == 2) {
    p <- suppressWarnings(wilcox.test(d$TFB ~ d$Response, exact = FALSE)$p.value)
  }
  data.frame(cancer = cc, n = nrow(d), p = p)
})

p_lab <- p_lab |>
  mutate(p_adj = p.adjust(p, method = FDR_METHOD),
         cancer = factor(cancer, levels = levels(box_df$cancer)),
         lab = sprintf('atop(italic(p) == "%s", italic(q (BH)) == "%s")', fmt_p(p), fmt_p(p_adj)),
         bold = if (USE_FDR_COLOUR) !is.na(p_adj) & p_adj < 0.05 else !is.na(p) & p < 0.05)

strip_lab <- box_df |> count(cancer) |> mutate(lab = sprintf("%s (n=%d)", cancer, n)) |> pull(lab, name = cancer)
y_expand <- if (LOG_Y_BOX) expansion(mult = c(0.02, 0.25)) else expansion(mult = c(0.05, 0.16))
n_box <- nlevels(droplevels(box_df$cancer))
box_ncol <- min(3, n_box)
box_nrow <- ceiling(n_box / box_ncol)

p_box <- ggplot(box_df, aes(x = Response, y = TFB + if (LOG_Y_BOX) 1 else 0)) +
  geom_boxplot(aes(fill = Response), width = .62, outlier.shape = NA, colour = "grey25", alpha = .95) +
  geom_jitter(aes(colour = Response), width = .14, size = 1, alpha = .45) +
  geom_text(data = p_lab, aes(x = 1.5, y = Inf, label = lab), parse = TRUE, 
            vjust = 1.55, size = 3.3, lineheight = .95, fontface = ifelse(p_lab$bold, "bold", "plain"), inherit.aes = FALSE) +
  facet_wrap(~ cancer, scales = "free_y", nrow = box_nrow, ncol = box_ncol, labeller = as_labeller(strip_lab)) +
  scale_fill_manual(values = c(R = COL_R, NR = COL_NR)) +
  scale_colour_manual(values = c(R = "#4a7c59", NR = "grey55")) +
  scale_y_continuous(expand = y_expand) +
  labs(x = NULL, y = "TFB (fusion count)", title = "TFB and ICB response") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey95", colour = NA),
        plot.caption = element_text(size = 12, colour = "grey40"),
        strip.text = element_text(face = "bold"), 
        plot.title = element_text(size = 15,hjust = .5, face = "bold"))

if (LOG_Y_BOX) p_box <- p_box + scale_y_log10()

ggsave(file.path(out_dir, "supFigure3E.pdf"), p_box, width = 10, height = 10.6)
ggsave(file.path(out_dir, "supFigure3E.png"), p_box, width = 10, height = 10.6, dpi = 600)

## ---- 7. Stacked bar 数据 --------------------
bar_stats <- patient_tfb |>
  filter(!is.na(Response), !is.na(TFB_group)) |>
  count(cancer, TFB_group, Response, name = "n") |>
  group_by(cancer, TFB_group) |>
  mutate(total = sum(n), pct = 100 * n / total) |> ungroup()

bar_p <- patient_tfb |>
  filter(!is.na(Response), !is.na(TFB_group)) |>
  group_by(cancer) |>
  group_modify(~ {
    tab <- table(.x$TFB_group, .x$Response)
    if (all(dim(tab) == 2) && all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
      a <- tab["High TFB", "R"]; b <- tab["High TFB", "NR"]
      c <- tab["Low TFB", "R"];  d <- tab["Low TFB", "NR"]
      o <- fisher_or(a, b, c, d)
      p_val <- fisher.test(tab)$p.value
      tibble(p = p_val, OR = o[["OR"]], low95 = o[["lo"]], up95 = o[["hi"]], n = nrow(.x))
    } else {
      tibble(p = NA_real_, OR = NA_real_, low95 = NA_real_, up95 = NA_real_, n = nrow(.x))
    }
  }) |>
  ungroup() |>
  mutate(p_adj = p.adjust(p, method = FDR_METHOD))

## ---- supFigure3F --------------------------------

mh_est <- tryCatch({
  lv  <- forest_tab$cancer
  arr <- array(NA_real_, dim = c(2, 2, length(lv)), dimnames = list(c("Low TFB", "High TFB"), c("NR", "R"), lv))
  ok  <- logical(length(lv))
  for (i in seq_along(lv)) {
    d <- patient_tfb[patient_tfb$cancer == lv[i] & !is.na(patient_tfb$Response), ]
    tb <- table(factor(d$TFB_group, levels = c("Low TFB", "High TFB")), factor(d$Response, levels = c("NR", "R")))
    if (all(dim(tb) == 2) && all(rowSums(tb) > 0) && all(colSums(tb) > 0)) {
      arr[, , i] <- as.matrix(tb); ok[i] <- TRUE
    }
  }
  if (sum(ok) < 1) stop("no valid stratum")
  mh <- mantelhaen.test(arr[, , ok, drop = FALSE])
  data.frame(OR = unname(mh$estimate), low95 = mh$conf.int[1], up95 = mh$conf.int[2], p = mh$p.value)
}, error = function(e) { message("M-H 合并失败：", conditionMessage(e)); NULL })

n_coh_tab <- patient_tfb |> distinct(cancer, cohort) |> count(cancer, name = "n_coh")

forest_plot_df <- forest_tab |>
  left_join(n_coh_tab, by = "cancer") |>
  arrange(OR) |>
  mutate(y = row_number(),
         label = sprintf("%s (n=%d, %d cohort%s)", cancer, n, n_coh, ifelse(n_coh > 1, "s", "")),
         right = sprintf("%.2f (%.2f-%.2f)\np =%s, q (BH) =%s", OR, low95, up95, fmt_p(p), fmt_p(p_adj)),
         hi_col = if (USE_FDR_COLOUR) sig_fdr else sig)

forest_plot_df$label <- factor(forest_plot_df$label, levels = forest_plot_df$label[order(forest_plot_df$y)])

X_LO <- 0.02; X_HI <- 100
forest_plot_df <- forest_plot_df |> mutate(lo = pmax(low95, X_LO), hi = pmin(up95, X_HI), c_lo = low95 < X_LO, c_hi = up95 > X_HI)

p_forest <- ggplot(forest_plot_df, aes(y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = hi_col), orientation = "y", width = .18, linewidth = .55) +
  geom_segment(data = dplyr::filter(forest_plot_df, c_hi), aes(x = hi, xend = X_HI * 1.35, y = label, yend = label, colour = hi_col), arrow = arrow(length = unit(2.2, "mm"), type = "closed"), linewidth = .55, show.legend = FALSE) +
  geom_segment(data = dplyr::filter(forest_plot_df, c_lo), aes(x = lo, xend = X_LO * 0.72, y = label, yend = label, colour = hi_col), arrow = arrow(length = unit(2.2, "mm"), type = "closed"), linewidth = .55, show.legend = FALSE) +
  geom_point(aes(x = OR, colour = hi_col, size = n), shape = 15) +
  geom_text(aes(x = X_HI * 1.55, label = right), hjust = 0, size = 2.9, lineheight = .9, colour = "grey15") +
  scale_x_log10(breaks = c(0.02, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10, 50), labels = c("0.02", "0.05", "0.1", "0.25", "0.5", "1", "2", "5", "10", "50")) +
  scale_colour_manual(values = c(`TRUE` = "#c0504d", `FALSE` = "grey35"), guide = "none") +
  scale_size_continuous(range = c(1.4, 4.2), guide = "none") +
  annotation_logticks(sides = "b", colour = "grey55", linewidth = .3) +
  labs(x = "Odds ratio for response: High vs Low TFB (log scale)", y = NULL,
       title = "Pan-cancer association between TFB and ICB response",
       caption = sprintf(paste0("OR < 1 = high TFB, lower response. Bars: \n95%% CI (Fisher's exact with Haldane correction); ",
                                "arrows: CI beyond axis.\nq (%s) across %d cancers."), FDR_METHOD, nrow(forest_plot_df))) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 11), 
        plot.title = element_text(size = 12,face = "bold", hjust = 0),
        plot.caption = element_text(size = 11, colour = "grey40", hjust = 0), plot.margin = margin(8, 235, 8, 8)) +
  coord_cartesian(xlim = c(X_LO, X_HI), clip = "off")

suppressWarnings({
  ggsave(file.path(out_dir, "supFigure3F.pdf"), p_forest, width = 10, height = 0.42 * nrow(forest_plot_df) + 2.2)
  ggsave(file.path(out_dir, "supFigure3F.png"), p_forest, width = 10, height = 0.42 * nrow(forest_plot_df) + 2.2, dpi = 600)
})

