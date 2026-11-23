suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(stringr); library(writexl)
})
library(writexl)
set.seed(2024)
root      <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3"
file_meta <- file.path(root, "数据", "metafusion_calls.star.tsv")
Credibility_gene <-fread(file_meta) |>
  separate_rows(gene1, sep = ",") |>
  separate_rows(gene2, sep = ",") |>
  dplyr::filter(sample_type %in% c("Tumor", "Tumor_Blood")) |>
  mutate(fusion = paste0(gene1, "::", gene2)) |>   # 保留 fusion 列，供其他用途
  setDT()

## ---- 0. 估计量：不变 ------------------------------------------------------
fisher_or <- function(a, b, c, d) {
  a <- a + .5; b <- b + .5; c <- c + .5; d <- d + .5
  or <- (a * d) / (b * c)
  se <- sqrt(1/a + 1/b + 1/c + 1/d)
  c(OR = or, lo = exp(log(or) - 1.96 * se), hi = exp(log(or) + 1.96 * se),
    log2OR = log2(or))
}
fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

MIN_N_FUSION <- 3
time_scenarios <- c("Pre", "On")
cohort_sum_irAE <- c("PHS003284", "PRJNA795330", "PHS003316", "PRJNA482620")

## ---- 1. prep_cohort：完全不变 --------------------------------------------
prep_cohort <- function(cohort_id) {
  dt <- fread(paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/", cohort_id, "_info.tsv"))
  
  if (cohort_id == "PRJNA795330") {
    dt[, time_point := fifelse(grepl("Pre", Sampling), "Pre",
                               fifelse(grepl("On",  Sampling), "On", NA_character_))]
    ae_cols <- grep("^irAE_days_", names(dt), value = TRUE)
    for (cl in ae_cols) set(dt, j = cl, value = !is.na(dt[[cl]]))
    dt[, irAE := Reduce(`|`, lapply(ae_cols, function(cl) dt[[cl]]))]
    setnames(dt, ae_cols, gsub("^irAE_days_", "ae_", ae_cols))
    type_cols <- gsub("^irAE_days_", "ae_", ae_cols)
    
  } else if (cohort_id == "PHS003316") {
    dt[, time_point := fifelse(grepl("Pre", Response2), "Pre",
                               fifelse(grepl("On",  Response2), "On", NA_character_))]
    dt[, irAE := !is.na(Autoimmune.toxicities) & Autoimmune.toxicities != "None"]
    type_cols <- character(0)
    
  } else if (cohort_id == "PHS003284") {
    dt[, time_point := "Pre"]
    parts <- c("1st irAE", "2nd irAE", "Grade >= 3 irAE")
    parts <- parts[parts %in% names(dt)]
    dt[, all_irAE_str := do.call(paste, c(lapply(parts, function(p) dt[[p]]), sep = ";"))]
    lvl1 <- c("Pneumonitis", "Hypothyroid", "Rash", "Infusion reaction",
              "Uveitis", "ALT Increased", "Bilirubin Increased")
    for (ae in lvl1)
      set(dt, j = paste0("ae_", gsub(" ", "_", ae)),
          value = !is.na(dt$all_irAE_str) & grepl(ae, dt$all_irAE_str, ignore.case = TRUE))
    type_cols <- paste0("ae_", gsub(" ", "_", lvl1))
    dt[, irAE := Reduce(`|`, lapply(type_cols, function(cl) dt[[cl]]))]
    
  } else if (cohort_id == "PRJNA482620") {
    dt[, time_point := "Pre"]
    map <- c("Endocrine toxicities (yes/No)" = "ae_Endocrine",
             "Rash (yes/no)"                 = "ae_Rash",
             "Liver Disfunction"             = "ae_Liver")
    for (old in names(map))
      set(dt, j = map[[old]],
          value = !is.na(dt[[old]]) & tolower(trimws(dt[[old]])) != "no")
    type_cols <- unname(map)
    dt[, irAE := Reduce(`|`, lapply(type_cols, function(cl) dt[[cl]]))]
  }
  
  dt <- dt[!is.na(time_point)]
  list(dt = dt, type_cols = type_cols)
}

## ---- 2. 主循环：gene1 / gene2 分开做 --------------------------------------
result_all <- list()

for (cohort_id in cohort_sum_irAE) {
  prep   <- prep_cohort(cohort_id)
  Need1  <- prep$dt
  types  <- c("ALL" = "irAE", setNames(prep$type_cols, sub("^ae_", "", prep$type_cols)))
  
  cohort_gene_all <- Credibility_gene |> dplyr::filter(cohort == cohort_id)
  
  # 5' 端(gene1) 与 3' 端(gene2) 各自单独成池，不再合并去重
  gene_sides <- list(
    "5'" = list(col = "gene1", pool = unique(cohort_gene_all$gene1)),
    "3'" = list(col = "gene2", pool = unique(cohort_gene_all$gene2))
  )
  
  for (time_pt in time_scenarios) {
    Need2 <- Need1[time_point == time_pt]
    if (nrow(Need2) == 0) next
    
    for (ae_name in names(types)) {
      ae_col <- types[[ae_name]]
      
      for (side in names(gene_sides)) {                 # 先分 5' / 3'
        gcol      <- gene_sides[[side]]$col
        gene_pool <- gene_sides[[side]]$pool
        gene_pool <- gene_pool[!is.na(gene_pool) & gene_pool != ""]
        
        for (g in gene_pool) {
          # 只在对应端匹配——5' 只看 gene1，3' 只看 gene2
          hit <- unique(cohort_gene_all[cohort_gene_all[[gcol]] == g, Sample_ID])
          if (length(hit) == 0) next
          
          d <- Need2[, .(pos = any(Run %in% hit),
                         ae  = any(get(ae_col))),
                     by = Patient_ID][!is.na(ae)]
          
          if (sum(d$pos) < MIN_N_FUSION) next
          
          a  <- sum( d$pos &  d$ae); b  <- sum( d$pos & !d$ae)
          cc <- sum(!d$pos &  d$ae); dd <- sum(!d$pos & !d$ae)
          if (a + b < MIN_N_FUSION || cc + dd == 0) next
          
          o <- fisher_or(a, b, cc, dd)
          p <- tryCatch(fisher.test(matrix(c(a, b, cc, dd), 2))$p.value,
                        error = function(e) NA_real_)
          
          result_all[[length(result_all) + 1]] <- data.table(
            Gene = g, Side = side,            # 基因名不加前缀，用 Side 区分
            cohort = cohort_id, irAE = ae_name, Time = time_pt,
            n_pos = a + b, n_neg = cc + dd,
            pos_ae = a, pos_noae = b, neg_ae = cc, neg_noae = dd,
            OR = o[["OR"]], lo = o[["lo"]], hi = o[["hi"]], log2OR = o[["log2OR"]], p = p)
        }
      }
    }
  }
}

res <- rbindlist(result_all)

## ---- 3. FDR / 排序 / 输出：不变 -------------------------------------------
res[, FDR     := p.adjust(p, method = "BH"), by = .(Time, irAE, Side)]
res[, FDR_all := p.adjust(p, method = "BH")]

setorder(res, p, na.last = TRUE)
write_xlsx(res, path = "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/附表/TableS7.xlsx")
