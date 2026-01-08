suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(writexl)
})
set.seed(2024)
## ---------- 0. paths ----------
root      <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3"
file_meta <- file.path(root, "数据", "metafusion_calls.star.tsv")
tme_dir   <- file.path(root, "数据")
out_xlsx  <- file.path(root, "附表", "TableS6.xlsx")

gene_name     <- "LSAMP"
cohort_choose <- c("PHS000452_PD1", "PHS000452_CTLA4")

Need1 <- fread(file.path(root, "数据/Clininfo/PHS000452_PD1_info.tsv"))
Need2 <- fread(file.path(root, "数据/Clininfo/PHS000452_CTLA4_info.tsv"))
sample_all <- c(Need1$Run, Need2$Run) |> unique()
On_sample <- Need1 |> filter(grepl("On", Response2)) |> pull(Run)

## ---------- 1. fusion status ----------
fusion_select <- fread(file_meta, sep = "\t") |>
  filter(sample_type != "Normal",
         cohort %in% cohort_choose,
         gene2 == gene_name) |>
  distinct(Run = Sample_ID)

## ---------- 2. read TME ----------
tme <- lapply(list.files(tme_dir, pattern = "_TME\\.tsv\\.gz$", full.names = TRUE), function(f) {
  df <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  df$cohort <- sub("_TME\\.tsv\\.gz$", "", basename(f))
  df
}) |>
  bind_rows() |>
  filter(cohort %in% cohort_choose) |>
  filter(!ID %in% On_sample)

meta_cols <- c("ID", "cohort")
drop_cols <- c("P-value_CIBERSORT", "Correlation_CIBERSORT", "RMSE_CIBERSORT",
               "P-value_CIBERSORT-ABS", "Correlation_CIBERSORT-ABS", "RMSE_CIBERSORT-ABS")
value_cols <- setdiff(names(tme), c(meta_cols, drop_cols))

tme <- tme |>
  mutate(purity_raw = as.numeric(TumorPurity_estimate),
         purity_lgt = qlogis(pmin(pmax(purity_raw, 0.01), 0.99)),
         fusion = factor(if_else(ID %in% fusion_select$Run, "pos", "neg"),
                         levels = c("neg", "pos")))

## ---------- 3. CLR —— 只对三个真正的成分数据方法 ----------
## CIBERSORT / EPIC / quanTIseq : 每样本细胞类型和 = 1  -> clr
## CIBERSORT-ABS / TIMER / MCPcounter / xCell / ESTIMATE / IPS -> 非成分，禁止 clr
clr_group <- function(d, cols) {
  if (!length(cols)) return(d)
  m <- as.matrix(as.data.frame(d)[, cols, drop = FALSE])
  m[!is.finite(m)] <- NA
  for (i in seq_len(nrow(m))) {
    v <- m[i, ]
    if (all(is.na(v))) next
    pos <- v[!is.na(v) & v > 0]
    if (!length(pos)) { m[i, ] <- 0; next }
    delta <- min(pos) / 2                 # 零值替换；严谨可用 cmultRepl
    v[!is.na(v) & v == 0] <- delta
    lv <- log(v)
    m[i, ] <- lv - mean(lv, na.rm = TRUE)
  }
  d[cols] <- as.data.frame(m)
  d
}

## 每个方法必须作为独立成分空间分别 clr（和都为 1，但互不混合）
comp_groups <- list(
  CIBERSORT = grep("CIBERSORT$", value_cols, ignore.case = TRUE, value = TRUE),
  EPIC      = grep("EPIC$",      value_cols, ignore.case = TRUE, value = TRUE),
  quanTIseq = grep("quanTIseq$", value_cols, ignore.case = TRUE, value = TRUE)
)
comp_groups <- comp_groups[lengths(comp_groups) > 0]
stopifnot(length(comp_groups) == 3)

for (g in names(comp_groups)) {
  cols <- comp_groups[[g]]
  s <- rowSums(as.data.frame(tme)[, cols, drop = FALSE], na.rm = TRUE)
  message(sprintf("%-10s : %2d cols | 行和 %.3f – %.3f",
                  g, length(cols), min(s, na.rm = TRUE), max(s, na.rm = TRUE)))
  tme <- clr_group(tme, cols)
}

## ---------- 4. long + within-cohort z ----------
long <- tme |>
  pivot_longer(all_of(value_cols), names_to = "variable", values_to = "value") |>
  filter(!is.na(value)) |>
  group_by(variable, cohort) |>
  mutate(value_z = as.numeric(scale(value))) |>
  ungroup()

## ---------- 5. main model (one BH family) ----------
fit_one <- function(d, adjust = FALSE) {
  if (n_distinct(d$fusion) < 2) return(NULL)
  f <- if (adjust){value_z ~ fusion + cohort + purity_lgt} else {value_z ~ fusion + cohort}
  m  <- lm(f, data = d)
  co <- summary(m)$coefficients
  if (!"fusionpos" %in% rownames(co)) return(NULL)
  ci <- confint(m)["fusionpos", ]
  data.table(variable = d$variable[1],
             beta = co["fusionpos", 1], se = co["fusionpos", 2],
             p = co["fusionpos", 4], low = ci[1], high = ci[2],
             n_pos = sum(d$fusion == "pos"), n_neg = sum(d$fusion == "neg"))
}

res      <- rbindlist(lapply(split(long, long$variable), fit_one))
res[, q  := p.adjust(p, method = "BH")]            # 全部 157 个变量同一族
res[, direction := ifelse(beta > 0, "higher in LSAMP+", "lower in LSAMP+")]

res_adj  <- rbindlist(lapply(split(long, long$variable), fit_one, adjust = TRUE))
res_adj <- res_adj[, .(variable,
                       beta_adj = beta, p_adj = p,
                       low_adj  = low,  high_adj = high,
                       q_adj    = p.adjust(p, method = "BH"))]

## ---------- 6. per-cohort (unadjusted, transparency only) ----------
per_cohort <- rbindlist(lapply(split(long, list(long$variable, long$cohort), drop = TRUE), function(d) {
  if (n_distinct(d$fusion) < 2) return(NULL)
  m  <- lm(value_z ~ fusion, data = d); co <- summary(m)$coefficients
  if (!"fusionpos" %in% rownames(co)) return(NULL)
  data.table(variable = d$variable[1], cohort = d$cohort[1],
             beta_cohort = co["fusionpos", 1], p_cohort = co["fusionpos", 4],
             n_pos_cohort = sum(d$fusion == "pos"),
             n_neg_cohort = sum(d$fusion == "neg"))
}))

dir_tab <- per_cohort[, .(dir_consistent = .N == 2 && all(sign(beta_cohort) == sign(beta_cohort[1]))),
                      by = variable]

## ---------- 7. assemble + write ----------
tab <- res |>
  left_join(res_adj, by = "variable") |>
  left_join(per_cohort[cohort == "PHS000452_PD1",
                       .(variable, beta_PD1 = beta_cohort, p_PD1 = p_cohort)], by = "variable") |>
  left_join(per_cohort[cohort == "PHS000452_CTLA4",
                       .(variable, beta_CTLA4 = beta_cohort, p_CTLA4 = p_cohort)], by = "variable") |>
  left_join(dir_tab, by = "variable") |>
  arrange(p) |>
  mutate(gene = "LSAMP",
         side = "3'") |>
  select(gene, side,variable,  beta, se, p, low, high, n_pos, n_neg, q, 
         beta_adj, p_adj, q_adj,
         beta_PD1, p_PD1, beta_CTLA4, p_CTLA4, 
         beta_adj, p_adj, q_adj) 
  

dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)
ok <- tryCatch({ write_xlsx(list(TableS6 = as.data.frame(tab)), out_xlsx); TRUE },
               error = function(e) { message("write failed: ", conditionMessage(e)); FALSE })
if (!ok) {
  out_xlsx <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/附表/TableS6.xlsx"
  write_xlsx(list(TableS6 = as.data.frame(tab)), out_xlsx)
  message("saved to: ", out_xlsx)
}