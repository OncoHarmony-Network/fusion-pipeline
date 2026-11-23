suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(swimplot); library(ggnewscale); library(scales); library(forcats)
})

## ---- 0. 路径 --------------------------------------------------------------
file_meta   <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
clin_dir    <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/"
meta_file   <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/cohort_meta_check.csv"
out_dir     <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图5/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

MIN_N_FUSION <- 3          # Fusion+ 少于此数不出图/不做检验
WINDOW_B     <- "Pre"      
TAG_B <- switch(WINDOW_B,
                Any = "any window (Pre or On)",
                Pre = "pre-treatment samples",
                On  = "on-treatment samples")
IRAE_TYPES   <- c("Skin", "Flu-like", "Gastrointestinal", "Musculoskeletal",
                  "Neurologic", "Pulmonary", "Thyroid")
CANCER_COL   <- c(Pre = "#DD5F60", On = "#3FB8AF")

## ---- 1. 工具函数 ----------------------------------------------------------
# Fisher OR（Haldane-Anscombe 校正，处理 0 格）+ 95% CI
fisher_or <- function(a, b, c, d) {
  a <- a + .5; b <- b + .5; c <- c + .5; d <- d + .5
  or <- (a * d) / (b * c)
  se <- sqrt(1/a + 1/b + 1/c + 1/d)
  c(OR = or, lo = exp(log(or) - 1.96 * se), hi = exp(log(or) + 1.96 * se))
}

# Clopper-Pearson 区间
cp_ci <- function(k, n) {
  if (n == 0) return(c(NA, NA))
  ci <- suppressWarnings(binom.test(k, n)$conf.int)
  c(ci[1] * 100, ci[2] * 100)
}

fmt_p <- function(p) ifelse(is.na(p), "NA",
                            ifelse(p < 0.001, "<0.001", sprintf("%.4f", p)))

## ---- 2. 数据 --------------------------------------------------------------
fusion_all <- fread(file_meta) |>
  separate_rows(gene1, gene1_type,sep = ",") |>
  separate_rows(gene2, gene2_type,sep = ",") 

cohort1 <- fread(meta_file)

# --- 2a. PRJNA795330（唯一有大量 irAE + 大量融合的队列） --------------------
CID <- "PRJNA795330"
Need <- fread(paste0(clin_dir, CID, "_info.tsv")) |>
  mutate(Sample_type = if_else(Sampling == "Pre", "Pre", "On"),
         irAE_any    = if_else(rowSums(across(starts_with("irAE_days_"), ~ !is.na(.))) > 0, "yes", "no"))

ae_safe <- setNames(gsub("-", "_", IRAE_TYPES), IRAE_TYPES)
for (ae in IRAE_TYPES) {
  col <- paste0("irAE_days_", ae)
  Need[[paste0("ae_", ae_safe[[ae]])]] <- if_else(!is.na(Need[[col]]), "yes", "no")
}

# 患者层面的窗口暴露：Pre / On 
fus_795 <- fusion_all |> filter(cohort == CID) |> select(Sample_ID, gene1, gene2) |> distinct()

exposure <- function(genes, window = "Any", side = c("any", "5p", "3p")) {
  side <- match.arg(side)
  hit <- switch(side,
                any = fus_795 |> filter(gene1 %in% genes | gene2 %in% genes),
                `5p` = fus_795 |> filter(gene1 %in% genes),
                `3p` = fus_795 |> filter(gene2 %in% genes)) |>
    pull(Sample_ID) |> unique()
  d <- Need
  if (window %in% c("Pre", "On")) d <- d |> filter(Sample_type == window)
  d |>
    mutate(pos = Run %in% hit) |>
    group_by(Patient_ID) |>
    summarise(pos = any(pos), irAE_any = first(irAE_any),
              across(starts_with("ae_"), first), .groups = "drop")
}

## ---- 3. Panel A：游泳图 ---------------------------------------------------
dA <- Need |>
  filter(Sample_type == "Pre") |>
  rowwise() |>
  mutate(irAE_time = suppressWarnings(max(c_across(starts_with("irAE_days_")), na.rm = TRUE))) |>
  ungroup() |>
  filter(irAE_any == "yes", is.finite(irAE_time)) |>
  arrange(desc(irAE_time)) |>
  as.data.frame()

gA_meta <- dA |> mutate(has_fusion = if_else(Run %in% fus_795$Sample_ID, "Fusion+", "Fusion-")) |>
  select(Patient_ID, has_fusion) |>
  mutate(facet = "Fusion") |>
  rename(event = has_fusion)

timeline <- dA |>
  select(Patient_ID, starts_with("irAE_days_")) |>
  pivot_longer(-Patient_ID, names_to = "ae", values_to = "time") |>
  filter(!is.na(time)) |>
  mutate(ae = sub("irAE_days_", "", ae))

dA$fus_status <- if_else(dA$Run %in% fus_795$Sample_ID, "Fusion+", "Fusion-")

# 手工画游泳图（不用 swimplot，避免 coord_flip + 离散轴的坑）
lv <- dA$Patient_ID[order(dA$irAE_time)]
dA$Patient_ID <- factor(dA$Patient_ID, levels = lv)
timeline_f <- timeline |> mutate(Patient_ID = factor(Patient_ID, levels = lv))

pA <- ggplot() +
  geom_segment(data = dA, aes(x = 0, xend = irAE_time, y = Patient_ID, yend = Patient_ID,
                              colour = fus_status),
               linewidth = 1.5, lineend = "butt") +
  scale_colour_manual(name = "Fusion status",
                      values = c(`Fusion+` = "#4a7ebb", `Fusion-` = "#e8e8e8")) +
  ggnewscale::new_scale_colour() +
  geom_point(data = timeline_f, aes(x = time, y = Patient_ID, shape = ae, colour = ae),
             size = 1.5, stroke = .5) +
  scale_colour_manual(name = "irAE event",
                      values = c(Skin = "#B092B6", `Flu-like` = "#4f9ec4",
                                 Gastrointestinal = "#3f9e8e", Musculoskeletal = "black",
                                 Neurologic = "#C74D26", Pulmonary = "#d98b25",
                                 Thyroid = "#7aa84f")) +
  scale_shape_manual(name = "irAE event",
                     values = c(Skin = 17, `Flu-like` = 16, Gastrointestinal = 18,
                                Musculoskeletal = 3, Neurologic = 5,
                                Pulmonary = 15, Thyroid = 14)) +
  scale_x_continuous(expand = expansion(mult = c(0, .02))) +
  labs(x = "Days since treatment", y = NULL,
       title = "Timing of irAE relative to ICB start (PRJNA795330, pre-treatment samples)",
       subtitle = sprintf("%d patients with at least one documented irAE; bars end at the time of last irAE", nrow(dA))) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.key.size = unit(.42, "cm"), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5),
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 6.4, colour = "grey30"))

ggsave(file.path(out_dir, "FigureS5_A.pdf"), pA, width = 7.2, height = 8.6)
ggsave(file.path(out_dir, "FigureS5_A.png"), pA, width = 7.2, height = 8.6, dpi = 600)

## ---- 4. Panel B：gene x irAE 效应量热图（替代原来的单基因堆叠柱） --------
TABLE_S7 <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/附表/TableS7.xlsx"
library(readxl)
COHORT_B     <- "PRJNA795330"
WINDOW_B     <- "Pre"
WINDOW_C     <- WINDOW_B          # 与 B 同一时间窗
TAG_C        <- "pre-treatment"
MIN_N_FUSION <- 3                 # 仅用于 caption 说明（表里的行已按此筛过）
## 统一 p 的显示：两位有效数字
fmt_p2 <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "g", digits = 2))
## ---- 1. 读表 ---------------------------------------------------------------
tab <- read_excel(TABLE_S7, sheet = 1) |> as.data.frame()
need_cols <- c("Gene", "Side", "cohort", "irAE", "Time", "n_pos", "n_neg",
               "pos_ae", "pos_noae", "neg_ae", "neg_noae",
               "OR", "lo", "hi", "p", "FDR")
miss <- setdiff(need_cols, names(tab))
if (length(miss)) stop("TableS7 缺列: ", paste(miss, collapse = ", "))
tab <- tab |>
  mutate(Side = as.character(Side),
         Time = as.character(Time),
         cohort = as.character(cohort))
#supFigure5B
dB0 <- tab |>
  filter(Time == WINDOW_B, cohort == COHORT_B) |>
  mutate(
    side = factor(Side, levels = c("5'", "3'")),
    gene = Gene,
    l2or = log2(OR),
    lab  = ifelse(!is.na(p) & p < 0.05, fmt_p2(p), "")
  ) |>
  filter(!Gene %in% DROP_GENES )
DROP_GENES <- c(
  "AL031432.2", "AL669831.1", "AC123912.3",   # pseudogene
  "HBA1", "HBA2", "HBB", "HBG1",              # hemoglobin
  "ACTB", "RPS27",                            # housekeeping
  "AL645608.2", "AL645608.6"                  # lncRNA
)

IRAE_ORDER <- c("ALL", setdiff(sort(unique(dB0$irAE)), "ALL"))
dB0 <- dB0 |> mutate(irAE = factor(irAE, levels = IRAE_ORDER))
stopifnot(!anyNA(dB0$irAE), !anyNA(dB0$side), !anyNA(dB0$l2or))
dB0 <- dB0 |>
  mutate(
    row_id = paste0(gene, "\u2016", side),
    ordv   = -log10(p)
  )
ord <- dB0 |>
  group_by(side, row_id) |>
  summarise(v = suppressWarnings(max(ordv, na.rm = TRUE)), .groups = "drop") |>
  arrange(side, v, row_id)
dB0$row_id <- factor(dB0$row_id, levels = ord$row_id)
strip_suffix <- function(x) sub("\u2016.*$", "", x)
MIN_FDR   <- suppressWarnings(min(dB0$FDR, na.rm = TRUE))
N_FAM     <- paste0(sprintf("%d", min(dB0$n_pos, na.rm = TRUE)), "-",
                    sprintf("%d", max(dB0$n_pos, na.rm = TRUE)))
N_G_5     <- n_distinct(dB0$gene[dB0$side == "5'"])
N_G_3     <- n_distinct(dB0$gene[dB0$side == "3'"])
N_CELL    <- nrow(dB0)
pB <- ggplot(dB0, aes(irAE, row_id, fill = l2or)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = lab), size = 2.0, colour = "grey10") +
  facet_wrap(~ side, ncol = 2, scales = "free_y") +
  scale_y_discrete(labels = strip_suffix) +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#c0392b",
                       midpoint = 0, limits = c(-3, 3), oob = squish,
                       name = "log2 OR") +
  labs(x = NULL, y = NULL,
       title = paste0("Fusion partner genes versus irAE in PRJNA795330\n",
                      "(pre-treatment samples)"),
       caption = paste0("Ten genes of four classes were excluded a priori; source data in sTable 7."))  +
  theme_minimal(base_size = 7.5) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6.6),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 6.2, colour = "grey30"),
        plot.caption = element_text(size = 5.8, colour = "grey35", hjust = 0))
pB 
ggsave(file.path(out_dir, "FigureS5_B.pdf"), pB, width = 7.2, height = 4.6)
ggsave(file.path(out_dir, "FigureS5_B.png"), pB, width = 7.2, height = 4.6, dpi = 600)

#supFigure5C
cC <- tab |>
  filter(Gene == "JAZF1", Side == "5'", irAE == "Flu-like",
         Time == WINDOW_C, cohort == COHORT_B)
if (nrow(cC) != 1) stop("TableS7 里 JAZF1/5'/Flu-like 的行数为 ", nrow(cC), "（应为 1）")
n_pos_C <- cC$n_pos
dC_long <- tibble::tibble(
  grp = factor(c("Fusion-", "Fusion-", "Fusion+", "Fusion+"),
               levels = c("Fusion-", "Fusion+")),
  flu = factor(c("Flu-like yes", "Flu-like no", "Flu-like yes", "Flu-like no"),
               levels = c("Flu-like yes", "Flu-like no")),
  n   = c(cC$neg_ae, cC$neg_noae, cC$pos_ae, cC$pos_noae)
) |>
  group_by(grp) |>
  mutate(pct = n / sum(n) * 100, N = sum(n)) |>
  ungroup()
grp_lab <- dC_long |> distinct(grp, N) |>
  mutate(lab = sprintf("%s\n(n=%d)", as.character(grp), N))
p_val <- cC$p
OR    <- cC$OR
low95 <- cC$lo
up95  <- cC$hi
pC <- ggplot(dC_long, aes(grp, pct, fill = flu)) +
  geom_col(width = .62, colour = "black", linewidth = .3) +
  geom_text(aes(label = ifelse(pct >= 6, sprintf("%d (%.1f%%)", n, pct), "")),
            position = position_stack(vjust = .5), size = 2.4) +
  scale_fill_manual(values = c(`Flu-like yes` = "#c0392b", `Flu-like no` = "#4a7ebb")) +
  scale_x_discrete(labels = setNames(grp_lab$lab, grp_lab$grp)) +
  labs(x = NULL, y = "Patients (%)", fill = NULL,
       title = "JAZF1-5' chimeric transcript and \nflu-like irAE",
       subtitle = sprintf("Fisher's exact P = %s\nOR %s (95%% CI %s-%s)",
                          fmt_p2(p_val),
                          formatC(OR, format = "f", digits = 2),
                          formatC(low95, format = "f", digits = 2),
                          formatC(up95,  format = "f", digits = 2)),
       caption = sprintf("%s, %s samples (Table S7).", COHORT_B, TAG_C)) +
  theme_classic(base_size = 8) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 6, colour = "grey30"),
        plot.caption = element_text(size = 6, colour = "grey35", hjust = 0))
pC
ggsave(file.path(out_dir, "FigureS5_C.pdf"), pC, width = 2.6, height = 3.0)
ggsave(file.path(out_dir, "FigureS5_C.png"), pC, width = 2.6, height = 3.0, dpi = 600)


## ---- 6. Panel D：检测率 by sample type x irAE（加 CI 与检验） -------------
build_D <- function(ae_col, label) {
  Need |>
    mutate(has_fusion = Run %in% fus_795$Sample_ID) |>
    group_by(Sample_type, status = .data[[ae_col]]) |>
    summarise(k = n_distinct(Patient_ID[has_fusion]),
              n = n_distinct(Patient_ID), .groups = "drop") |>
    rowwise() |>
    mutate(pct = k / n * 100,
           lo = cp_ci(k, n)[1], hi = cp_ci(k, n)[2],
           irAE_type = label) |>
    ungroup()
}

dD <- bind_rows(lapply(c(list(list("irAE_any", "Any")),
                         setNames(lapply(IRAE_TYPES, function(x) list(paste0("ae_", ae_safe[[x]]), x)), IRAE_TYPES)),
                       function(z) build_D(z[[1]], z[[2]]))) |>
  mutate(y1 = 60)

# 每个 panel × sample type 内 irAE+ vs irAE- 的 Fisher
dD_p <- dD |>
  group_by(irAE_type, Sample_type) |>
  summarise(k = list(k), n = list(n), .groups = "drop") |>
  rowwise() |>
  mutate(p = {
    yy <- unlist(k); nn <- unlist(n)
    if (length(yy) == 2) fisher.test(matrix(c(yy, nn - yy), nrow = 2))$p.value else NA_real_
  }) |>
  ungroup() |>
  mutate(lab = sprintf("P = %s", fmt_p(p)),
         x = Sample_type) 

dD$irAE_type <- factor(dD$irAE_type, levels = c("Any", IRAE_TYPES))
Y_TOP <- 68     # 轴上限
Y_P   <- 65      # P 值标签基线（vjust = 1，向下排）
Y_KN  <- 60      # k/n 标签基线（vjust = 0，向上排）

pD <- ggplot(dD, aes(Sample_type, pct, fill = status)) +
  geom_col(position = position_dodge(.78), width = .72,
           colour = "black", linewidth = .25) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(.78), width = .18,
                linewidth = .28, colour = "grey25") +
  geom_text(aes(label = sprintf("%d/%d", k, n), y = Y_KN),
            position = position_dodge(.78), vjust = 0, size = 1.8,
            show.legend = FALSE) +
  geom_text(data = dD_p, aes(x = x, y = Y_P, label = lab),
            inherit.aes = FALSE, vjust = 1, size = 1.7) +
  facet_wrap(~ irAE_type, nrow = 2) +
  scale_fill_manual(values = c(yes = "#DD5F60", no = "#8fb8d8"),
                    labels = c(yes = "irAE+", no = "irAE-"), name = NULL) +
  coord_cartesian(ylim = c(0, Y_TOP), expand = FALSE) +
  scale_y_continuous(breaks = seq(0, 60, by = 20)) +
  labs(x = NULL, y = "Fusion detection rate (%)",
       title = "Fusion detection rate by sample type and irAE status (PRJNA795330)",
       subtitle = "Bars: proportion of patients with >=1 fusion; whiskers: Clopper-Pearson 95% CI; P: Fisher exact, irAE+ vs irAE- within sample type") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7, face = "bold"),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 6, colour = "grey30"))
pD
ggsave(file.path(out_dir, "FigureS5_D.pdf"), pD, width = 7.2, height = 3)
ggsave(file.path(out_dir, "FigureS5_D.png"), pD, width = 7.2, height = 3, dpi = 600)

## ---- 7. Panel E：Pre vs On 检测率（重算）+ 配对 McNemar -------------------
both_cohorts <- cohort1 |> filter(grepl("Pre", time) & grepl("On", time)) |> pull(id) |> unique()

sample_type_of <- function(d, cid) {
  if (cid == "PRJNA356761") return(if_else(grepl("Pre", d$Sample_ID), "Pre", "On"))
  if (cid == "PRJNA795330") return(if_else(d$Sampling == "Pre", "Pre", "On"))
  if (cid == "EGAD00001006282") return(if_else(grepl("_SCREEN", d$Sample_ID), "Pre", "On"))
  if (cid == "PHS001464") return(case_when(
    d$TUMOR_COLLECTION_TIMEPOINT == "Acquired Resistance" ~ "On",
    d$TUMOR_COLLECTION_TIMEPOINT == "Pre-Immunotherapy" ~ "Pre",
    d$TUMOR_COLLECTION_TIMEPOINT == "Off Therapy Recurrence" ~ "On",
    d$TUMOR_COLLECTION_TIMEPOINT == "Diagnosis" ~ "Pre",
    TRUE ~ NA_character_))
  case_when(grepl("Pre", d$Response2) | grepl("Pre", d$Response) ~ "Pre",
            grepl("On",  d$Response2) | grepl("On",  d$Response) ~ "On",
            TRUE ~ NA_character_)
}

paired_list <- list()
for (cid in setdiff(both_cohorts, "PHS001464")) {
  f <- file.path(clin_dir, paste0(cid, "_info.tsv"))
  if (!file.exists(f)) next
  d <- fread(f)
  d <- d |> mutate(Sample_type = sample_type_of(d, cid))
  if (!"Patient_ID" %in% names(d)) next
  paired_list[[cid]] <- d |> filter(!is.na(Sample_type)) |>
    select(Patient_ID, Run, Sample_type) |>
    mutate(cohort = cid)
}
paired <- rbindlist(paired_list, fill = TRUE) |>
  group_by(cohort, Patient_ID) |>
  filter(all(c("Pre", "On") %in% Sample_type)) |>
  ungroup()

fus_pair <- fusion_all |> filter(cohort %in% unique(paired$cohort)) |>
  select(cohort, Sample_ID, gene1, gene2) |> distinct() |>
  mutate(pos = TRUE)

## 癌种一律从 cohort1 查（id -> cancer_type），不硬编码
cancer_map <- cohort1 |> select(id, cancer_type) |> distinct()
stopifnot(!anyDuplicated(cancer_map$id))   # id 必须唯一，否则 join 会炸行

paired <- paired |>
  left_join(fus_pair, by = c("cohort", "Run" = "Sample_ID")) |>
  mutate(pos = !is.na(pos)) |>
  group_by(cohort, Patient_ID, Sample_type) |>
  summarise(pos = any(pos), .groups = "drop") |>
  pivot_wider(names_from = Sample_type, values_from = pos, values_fill = FALSE) |>
  left_join(cancer_map, by = c("cohort" = "id"))

## 查不到就按队列编号单独成组并报警，绝不静默归到某个癌种
unmapped <- unique(paired$cohort[is.na(paired$cancer_type)])
if (length(unmapped))
  warning("cohort 不在 cohort1$id 中，已按编号单独成组：", paste(unmapped, collapse = ", "),
          call. = FALSE)

paired <- paired |>
  mutate(cancer = coalesce(cancer_type, cohort)) |>
  select(-cancer_type)

summarise_E <- function(d, label) {
  if (nrow(d) == 0) return(NULL)
  b  <- sum(d$Pre & !d$On)      # Pre+ / On-
  c_ <- sum(!d$Pre & d$On)      # Pre- / On+
  ## 精确 McNemar = 不一致对上的双侧二项检验；b+c = 0 时无定义
  p_mc <- if (b + c_ == 0) NA_real_ else binom.test(b, b + c_, 0.5)$p.value
  data.frame(cancer = label, n = nrow(d),
             Pre = sum(d$Pre) / nrow(d) * 100,
             On  = sum(d$On)  / nrow(d) * 100,
             p_mcnemar = p_mc,
             discordant = b + c_)
}

dE <- bind_rows(
  lapply(sort(unique(stats::na.omit(paired$cancer))),
         function(cc) summarise_E(paired[which(paired$cancer == cc), , drop = FALSE], cc)),
  summarise_E(paired, "All"))
cat("paired cohorts:", paste(unique(paired$cohort), collapse = ", "), "\n")
cat("paired cancers:", paste(unique(paired$cancer), collapse = ", "), "\n")

## 自检：各组 n 之和必须等于 All
stopifnot(sum(dE$n[-nrow(dE)]) == dE$n[nrow(dE)])

dE$cancer <- factor(dE$cancer, levels = dE$cancer[order(dE$n)])
dE_long <- dE |> select(cancer, n, Pre, On, p_mcnemar, discordant) |>
  pivot_longer(c(Pre, On), names_to = "Sample_type", values_to = "pct")


p_lab <- function(p) ifelse(is.na(p), "n/a", fmt_p(p))
pE <- ggplot(dE_long, aes(cancer, pct, fill = Sample_type)) +
  geom_col(position = position_dodge(.75), width = .7, colour = "black", linewidth = .25) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(.75),
            vjust = -0.45, size = 2.1) +
  geom_text(data = dE, aes(cancer, Inf, label = sprintf("P = %s", p_lab(p_mcnemar))),
            inherit.aes = FALSE, vjust = 1.3, size = 2.1) +
  scale_fill_manual(values = CANCER_COL, name = NULL) +
  scale_x_discrete(labels = setNames(sprintf("%s (n=%d)", dE$cancer, dE$n), as.character(dE$cancer))) +
  scale_y_continuous(expand = expansion(mult = c(0, .18)), limits = c(0, 108)) +
  labs(x = NULL, y = "Fusion detection rate (%)",
       title = "Fusion detection in paired pre- and on-treatment samples",
       subtitle = "Patients with both samples; P: McNemar exact test on the paired 2x2 table\nn/a: no discordant pairs, test not evaluable.") +
  theme_classic(base_size = 7.5) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 6, colour = "grey30"))

ggsave(file.path(out_dir, "FigureS5_E.pdf"), pE, width = 4.4, height = 3.2)
ggsave(file.path(out_dir, "FigureS5_E.png"), pE, width = 4.4, height = 3.2, dpi = 600)

## ---- 8. 预览拼图 ----------------------------------------------------------
top  <- pB + pC + plot_layout(widths = c(3, 1.4))
prev <- top / pD  / pE
prev
ggsave(file.path(out_dir, "FigureS5_BCD_E.pdf"), prev,
       width = 8, height = 10)
