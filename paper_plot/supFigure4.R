library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tidyr)
library(stringr)
library(patchwork)
library(rstatix)
library(IOBR)
library(forestploter	)
library(grid)
set.seed(2024)

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
    yalb <- "Overall survival %"
  }else{
    fit <- survfit(Surv(PFS_Time, PFS_Status) ~ event, data = data)
    yalb <- "Progression Free Survival %"
  }
  p1 <- ggsurvplot(fit,
                   data = data,
                   linetype = "strata",
                   risk.table = TRUE,
                   risk.table.col = "strata",
                   legend.title = "",
                   pval = FALSE,
                   title = title,
                   ylab = yalb,
                   font.y = c(20, "plain", "black"),
                   font.legend = c(15, "plain", "black"),
                   font.title = c(20, "bold", "black"),
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

file_meta <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta , sep = "\t") 
cohort1 <- fread("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/cohort_meta_check.csv")
#sup4A
gene_name <- "LSAMP"
cohort_choose <- c("PHS000452_PD1", "PHS000452_CTLA4")
fusion_select <- fusion_all |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_choose ) |> 
  select(Sample_ID, gene1, gene2,cohort)|>
  rename(Run = Sample_ID) |>
  filter( gene2 %in% gene_name) |>
  as.data.table() |>
  distinct()
OS_all <- data.table()

for(cohort_id in cohort_choose){
  Need <- fread(paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/",cohort_id,"_info.tsv"))
  if(filter(cohort1 ,id == cohort_id)$time == "Pre/On" ){
    if(cohort_id == "PRJNA356761"){
      Need1 <- Need |>
        filter(grepl("Pre", Sample_ID))
    }else if(cohort_id == "EGAD00001006282"){
      Need1 <- Need |>
        filter(grepl("_SCREEN", Sample_ID))
    }else{
      Need1 <- Need |>
        filter(grepl("Pre", Response2))
    }}else if(cohort_id == "PRJNA744780"){
      Need1 <- Need |>
        filter(grepl("Pre", Response2))
    }else{
      Need1 <- Need
    }
  columns_to_select <- c("Patient_ID", "OS_Time", "OS_Status",  
                         "PFS_Time", "PFS_Status", 
                         "Run", "Treatment","Response2")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    select(all_of(existing_columns)) #|>   filter(!is.na(OS_Time) & !is.na(OS_Status))
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
  
  OS_all <- rbind(Need3,OS_all) |> 
    filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) |> 
    filter(!grepl("On", Response2)) 
}
fusion_data_subset <- OS_all |>
  filter(cohort == "PHS000452_CTLA4")|>
  mutate(event = case_when(
    Run %in% fusion_select$Run  ~  "LSAMP Fusion+",
    TRUE ~ "LSAMP Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "LSAMP Fusion+") ~ "LSAMP Fusion+",
                           TRUE ~ event))  |>
  group_by(cohort) |>
  filter(n_distinct(event) == 2) |>
  ungroup()  |>
  select(Patient_ID, event, OS_Time, OS_Status, PFS_Time, PFS_Status, Response2,cohort) |>
  distinct() 
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("PHS000452_CTLA4-OS\nLSAMP-3'"), 
                 choose = "OS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4A.pdf")
)
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("PHS000452_CTLA4-OS\nLSAMP-3'"), 
                 choose = "OS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4A.png")
)
#sup4B
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c("LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("PHS000452_CTLA4-LSAMP-3'"), 
                 choose = "PFS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4B.pdf")
)
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("PHS000452_CTLA4-LSAMP-3'"), 
                 choose = "PFS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4B.png")
)
#sup4C
fusion_subset <- fusion_data_subset |>
  select(Patient_ID, Response2, event) |>
  rename(Response = Response2) |> 
  filter(Response != "" & !is.na(Response ) & Response != "not collected" ) |> 
  distinct() |>
  as.data.table()
fusion_subset$Response  <- gsub("Pre-", "", fusion_subset$Response )
fusion_subset$Response <- factor(fusion_subset$Response, levels = c("NR", "R"))
fusion_subset$event    <- factor(fusion_subset$event,    levels = c("LSAMP Fusion+", "LSAMP Fusion-"))
tab <- table(fusion_subset$event, fusion_subset$Response)
stopifnot(all(dim(tab) == c(2, 2)))

a  <- tab["LSAMP Fusion+", "R"];  b  <- tab["LSAMP Fusion+", "NR"]
cc <- tab["LSAMP Fusion-", "R"];  dd <- tab["LSAMP Fusion-", "NR"]

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

fmt_p <- function(p) if (p < 0.001) "< 0.001" else sprintf("%.3f", p)

p_text <- bquote(bold(
  bolditalic(P) == .(fmt_p(p_val)) ~
    "; OR " == .(sprintf("%.2f", o[["OR"]])) ~
    "(95% CI " == .(sprintf("%.2f", o[["lo"]])) * "-" *
    .(sprintf("%.2f", o[["hi"]])) * ")"
))

p <- ggstatsplot::ggbarstats(
  data            = fusion_subset,
  x               = Response,      
  y               = event,      
  paired          = FALSE,      
  results.subtitle = FALSE,  
  label           = "both", 
  label.args       = list(color = NA, alpha = 0),  
  title           = paste0("PHS000452_CTLA4-LSAMP-3'"),
  ggtheme = theme_classic()
) + scale_fill_manual(values = c(
  "R" = "#9dd3af",
  "NR" = "#f2f2f2"
)) + 
  labs(
    subtitle = p_text, x = NULL, y = NULL) +  theme(
      text = element_text( size = 14),
      title = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold")
    )
p


filename <- paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4C.png")
ggsave(filename = filename, 
       plot = p,  
       width = 6, 
       height = 6,
       dpi = 600)

filename <- paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4C.pdf")
ggsave(filename = filename, 
       plot = p,  
       width = 6, 
       height = 6,
       dpi = 600)

#sup4D
fusion_data_subset <- OS_all |>
  filter(cohort == "PHS000452_PD1")|>
  mutate(event = case_when(
    Run %in% fusion_select$Run  ~  "LSAMP Fusion+",
    TRUE ~ "LSAMP Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "LSAMP Fusion+") ~ "LSAMP Fusion+",
                           TRUE ~ event))  |>
  group_by(cohort) |>
  filter(n_distinct(event) == 2) |>
  ungroup()  |>
  select(Patient_ID, event, OS_Time, OS_Status, PFS_Time, PFS_Status, Response2,cohort) |>
  distinct() 
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("PHS000452_PD1-LSAMP-3'"), 
                 choose = "OS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4D.pdf")
)
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("PHS000452_PD1-LSAMP-3'"), 
                 choose = "OS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4D.png")
)
#sup4E
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c("LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("PHS000452_PD1-LSAMP-3'"), 
                 choose = "PFS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4E.pdf")
)
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("PHS000452_PD1-LSAMP-3'"), 
                 choose = "PFS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4E.png")
)

#sup4F
fusion_subset <- fusion_data_subset |>
  select(Patient_ID, Response2, event) |>
  rename(Response = Response2) |> 
  filter(Response != "" & !is.na(Response ) & Response != "not collected" ) |> 
  distinct() |>
  as.data.table()
fusion_subset$Response  <- gsub("Pre-", "", fusion_subset$Response )
fusion_subset$Response <- factor(fusion_subset$Response, levels = c("NR", "R"))
fusion_subset$event    <- factor(fusion_subset$event,    levels = c("LSAMP Fusion+", "LSAMP Fusion-"))
tab <- table(fusion_subset$event, fusion_subset$Response)
stopifnot(all(dim(tab) == c(2, 2)))

a  <- tab["LSAMP Fusion+", "R"];  b  <- tab["LSAMP Fusion+", "NR"]
cc <- tab["LSAMP Fusion-", "R"];  dd <- tab["LSAMP Fusion-", "NR"]

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

fmt_p <- function(p) if (p < 0.001) "< 0.001" else sprintf("%.3f", p)

p_text <- bquote(bold(
  bolditalic(P) == .(fmt_p(p_val)) ~
    "; OR " == .(sprintf("%.2f", o[["OR"]])) ~
    "(95% CI " == .(sprintf("%.2f", o[["lo"]])) * "-" *
    .(sprintf("%.2f", o[["hi"]])) * ")"
))

p <- ggstatsplot::ggbarstats(
  data            = fusion_subset,
  x               = Response,      
  y               = event,      
  paired          = FALSE,      
  results.subtitle = FALSE,  
  label           = "both", 
  label.args       = list(color = NA, alpha = 0),  
  title           = paste0("PHS000452_PD1-LSAMP-3'"),
  ggtheme = theme_classic()
) + scale_fill_manual(values = c(
  "R" = "#9dd3af",
  "NR" = "#f2f2f2"
)) + 
  labs(
    subtitle = p_text, x = NULL, y = NULL) +  theme(
      text = element_text( size = 14),
      title = element_text( size = 20, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold")
    )



filename <- paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4F.png")
ggsave(filename = filename, 
       plot = p,  
       width = 6, 
       height = 6,
       dpi = 600)

filename <- paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4F.pdf")
ggsave(filename = filename, 
       plot = p,  
       width = 6, 
       height = 6,
       dpi = 600)

#sup4G
sym_log_axis <- function(low, upper, ref = 1, pad = 1.03, max_ticks = 6) {
  lo   <- suppressWarnings(min(c(low,   ref), na.rm = TRUE))
  hi   <- suppressWarnings(max(c(upper, ref), na.rm = TRUE))
  need <- max(ref / lo, hi / ref)                     # 需要的对数半宽
  ladder <- c(1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 50)
  m <- ladder[which(ladder >= need * pad)[1]]         # 向外取整齐的边界
  if (is.na(m)) m <- need * pad * 1.1
  xlim <- c(ref / m, ref * m)                         # ★ log 空间严格对称 → 1 居中
  pool <- c(0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 4, 5, 10, 20)
  ticks <- pool[pool >= xlim[1] * 0.999 & pool <= xlim[2] * 1.001]
  if (length(ticks) > max_ticks) {                    # 刻度太密就换稀疏档
    coarse <- c(0.05, 0.1, 0.5, 1, 2, 10, 20)
    ticks <- pool[pool %in% coarse & pool >= xlim[1] * 0.999 & pool <= xlim[2] * 1.001]
  }
  list(xlim = xlim, ticks_at = sort(unique(ticks)))
}
out_dir <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure4"
file_meta <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
metafusion <- fread(file_meta, sep = "\t")

# ★ 只保留 LSAMP 参与的融合
Credibility_gene <- metafusion |>
  filter(sample_type %in% c("Tumor", "Tumor_Blood")) |>
  filter( gene2 == "LSAMP") |>
  as.data.table()

message("LSAMP fusion records: ", nrow(Credibility_gene),
        " | unique samples: ", uniqueN(Credibility_gene$Sample_ID))

# ---------------- 2. 临床 ----------------
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

pat_os <- all_clin[, .(OS_Time   = OS_Time[1],
                       OS_Status = OS_Status[1]),
                   by = .(cohort, Patient_ID)]

smp2pat <- all_clin[, .(Run, Patient_ID, cohort)] |> unique()
stopifnot(!anyDuplicated(smp2pat$Run))

# ---------------- 3. 患者层面 LSAMP+ ----------------
pat_fusion_pos <- Credibility_gene |>
  dplyr::distinct(Sample_ID) |>
  inner_join(smp2pat, by = c("Sample_ID" = "Run")) |>
  dplyr::distinct(cohort, Patient_ID) |>
  as.data.table()

stopifnot(!anyDuplicated(pat_fusion_pos, by = c("cohort", "Patient_ID")))

message("LSAMP+ patients with OS: ", nrow(pat_fusion_pos),
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
             n_strata   = 1L)
}

res_cohort <- rbindlist(lapply(unique(pat_os$cohort), run_cohort), fill = TRUE)
if(nrow(res_cohort) == 0)
  stop("No cohort passed LSAMP filters — check mapping and n thresholds.")

# ---------------- 5. FDR（仅 per-cohort，按 TCGA / TARGET / CPTAC 三族 BH） ----------------
res <- copy(res_cohort)
res[, finite := is.finite(HR) & is.finite(low) & is.finite(up) & !is.na(p)]

res_note <- res[finite == FALSE]
res_ok   <- res[finite == TRUE]

res_ok[, source := fifelse(grepl("^TCGA",   cohort), "TCGA",
                           fifelse(grepl("^TARGET", cohort), "TARGET",
                                   fifelse(grepl("^CPTAC",  cohort), "CPTAC", NA_character_)))]
stopifnot(!anyNA(res_ok$source))

res_ok[, q   := p.adjust(p, method = "BH"), by = source]
res_ok[, sig := q < 0.05]

message("LSAMP OS — cohorts passed: ", nrow(res_ok),
        " | sig(q<0.05): ", res_ok[sig == TRUE, .N],
        " | Inf/non-estimable: ", nrow(res_note))
print(res_ok[order(source, q),
             .(cohort, n_pos, n_neg, HR, low, up, p, q, sig)])


# ---------------- 6. 森林图 ----------------
tm <- forest_theme(
  base_size = 10,
  colhead  = list(fg_params = list(hjust = 0)),
  core     = list(padding = unit(c(5, 5), "mm")),
  ci_pch = 15, ci_col = "black", ci_fill = "black",
  ci_alpha = 0.8, ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.2,
  refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey20"),
  footnote_gp = gpar(cex = 0.6, fontface = "plain", col = "black"),
  vertline_lwd = 1, vertline_lty = "dashed", vertline_col = "grey20"
)

make_forest <- function(dt, title, out_pdf){
  if(nrow(dt) == 0){ message("empty: ", title); return(invisible(NULL)) }
  
  dt <- copy(dt)
  setorder(dt, source, q)
  
  dt[, `Fusion+/Fusion-` := paste0(n_pos, "/", n_neg)]
  dt[, `HR (95% CI)`     := sprintf("%.2f (%.2f-%.2f)", HR, low, up)]
  dt[, p_lab             := ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))]
  dt[, q_lab             := ifelse(is.na(q), "",
                                   ifelse(q < 0.001, "<0.001", sprintf("%.3f", q)))]
  dt[, blank             := paste(rep(" ", 30), collapse = " ")]
  
  plot_dt <- as.data.frame(dt[, .(cohort, `Fusion+/Fusion-`, blank,
                                  `HR (95% CI)`, p_lab, q_lab)])
  colnames(plot_dt) <- c("Cohort", "Fusion", " ",
                         "HR (95% CI)", "p", "q (BH)")
  
  ax <- sym_log_axis(dt$low, dt$up)
  
  note_txt <- if(nrow(res_note) > 0)
    paste0("Not shown (HR = Inf / non-estimable): ",
           paste(res_note$cohort, collapse = "; ")) else ""
  
  footnote_txt <- paste(
    "\n",
    "\n",
    "\n",
    "\n",
    "Per-cohort Cox:\nSurv ~ LSAMP Fusion status.\n",
    "q: BH-FDR within source family (TCGA / TARGET / CPTAC).\n",
    note_txt
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
  n_fn_lines <- str_count(footnote_txt, "\n") + 1
  h <- max(6, nrow(plot_dt) * 0.45 + 3 + n_fn_lines * 0.22)
  
  pdf(out_pdf, width = 10, height = h, onefile = FALSE)
  p <- forest(plot_dt,
              est = dt$HR, lower = dt$low, upper = dt$up,
              ci_column = 3,
              ref_line = 1,
              xlim     = ax$xlim,      # ← 对称边界
              ticks_at = ax$ticks_at,  # ← 整齐刻度
              x_trans  = "log",        # ★ 关键：对数轴
              ci_pch = 15, ci_lwd = 1.5,
              title = title, footnote = footnote_txt, theme = tm)
  p <- insert_text(p, text = "Fusion+/Fusion-", row = 1, col = 2,
                   just = "left",
                   gp = gpar(cex = 0.7, col = "black", fontface = "italic"))
  print(p); dev.off()
}

make_forest(res_ok,
            "GDC-OS-LSAMP-3'",
            "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4G.pdf")
#sup4h
folder_path <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/Clininfo/"
clin_files  <- list.files(folder_path,
                          pattern = "^(CPTAC-|TARGET-|TCGA-).*_info\\.tsv$",
                          full.names = TRUE)

cohort_ids <- gsub("_", "-", sub("_info\\.tsv$", "", basename(clin_files)))

# ★ PFS 版额外排除无 PFS 数据的 4 个队列
keep <- !cohort_ids %in% c("CPTAC-NORMAL", "TARGET-NORMAL",
                           "TCGA-NORMAL",  "TARGET-CELL",
                           "TARGET-NBL",   "TARGET-GNB",
                           "TARGET-LAML",  "TCGA-LAML")
clin_files <- clin_files[keep]
cohort_ids <- cohort_ids[keep]

read_one_pfs <- function(f, ch){
  dt <- fread(f)
  data.table(
    Patient_ID = dt$Patient_ID,
    Run        = dt$Run,
    PFS_Time   = if ("PFS_Time"   %in% names(dt)) as.numeric(dt$PFS_Time)   else NA_real_,
    PFS_Status = if ("PFS_Status" %in% names(dt)) as.integer(dt$PFS_Status) else NA_integer_,
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

# 样本 -> 患者映射
smp2pat_pfs <- all_clin_pfs[, .(Run, Patient_ID, cohort)] |> unique()
stopifnot(!anyDuplicated(smp2pat_pfs$Run))

# ---------------- 3. 患者层面 LSAMP+（基于 PFS 队列） ----------------
pat_fusion_pos <- Credibility_gene |>
  dplyr::distinct(Sample_ID) |>
  inner_join(smp2pat_pfs, by = c("Sample_ID" = "Run")) |>
  dplyr::distinct(cohort, Patient_ID) |>
  as.data.table()

stopifnot(!anyDuplicated(pat_fusion_pos, by = c("cohort", "Patient_ID")))

message("LSAMP+ patients with PFS: ", nrow(pat_fusion_pos),
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
  
  data.table(cohort = ch,
             n_pos  = n_pos,
             n_neg  = n_neg,
             HR     = s$coefficients[1, "exp(coef)"],
             low    = s$conf.int[1, "lower .95"],
             up     = s$conf.int[1, "upper .95"],
             p      = s$coefficients[1, "Pr(>|z|)"])
}

res_cohort_pfs <- rbindlist(lapply(unique(pat_pfs$cohort), run_cohort_pfs), fill = TRUE)
if(nrow(res_cohort_pfs) == 0)
  stop("No cohort passed LSAMP PFS filters — check mapping and n thresholds.")

# ---------------- 5. FDR（仅 per-cohort，按 TCGA / TARGET / CPTAC 三族 BH） ----------------
res <- copy(res_cohort_pfs)
res[, finite := is.finite(HR) & is.finite(low) & is.finite(up) & !is.na(p)]

res_note <- res[finite == FALSE]
res_ok   <- res[finite == TRUE]

res_ok[, source := fifelse(grepl("^TCGA",   cohort), "TCGA",
                           fifelse(grepl("^TARGET", cohort), "TARGET",
                                   fifelse(grepl("^CPTAC",  cohort), "CPTAC", NA_character_)))]
stopifnot(!anyNA(res_ok$source))

res_ok[, q   := p.adjust(p, method = "BH"), by = source]
res_ok[, sig := q < 0.05]

message("LSAMP PFS — cohorts passed: ", nrow(res_ok),
        " | sig(q<0.05): ", res_ok[sig == TRUE, .N],
        " | Inf/non-estimable: ", nrow(res_note))
print(res_ok[order(source, q),
             .(cohort, n_pos, n_neg, HR, low, up, p, q, sig)])

fwrite(res_ok, file.path(out_dir, "Figure4F_LSAMP_PFS_cox.tsv"), sep = "\t")

# ---------------- 6. 森林图 ----------------
tm <- forest_theme(
  base_size = 10,
  colhead  = list(fg_params = list(hjust = 0)),
  core     = list(padding = unit(c(5, 5), "mm")),
  ci_pch = 15, ci_col = "black", ci_fill = "black",
  ci_alpha = 0.8, ci_lty = 1, ci_lwd = 1.5, ci_Theight = 0.2,
  refline_gp  = gpar(lwd = 1, lty = "dashed", col = "grey20"),
  footnote_gp = gpar(cex = 0.6, fontface = "plain", col = "black"),
  vertline_lwd = 1, vertline_lty = "dashed", vertline_col = "grey20"
)

make_forest <- function(dt, title, out_pdf){
  if(nrow(dt) == 0){ message("empty: ", title); return(invisible(NULL)) }
  
  dt <- copy(dt)
  setorder(dt, source, q)
  
  dt[, `Fusion+/Fusion-` := paste0(n_pos, "/", n_neg)]
  dt[, `HR (95% CI)`     := sprintf("%.2f (%.2f-%.2f)", HR, low, up)]
  dt[, p_lab             := ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))]
  dt[, q_lab             := ifelse(is.na(q), "",
                                   ifelse(q < 0.001, "<0.001", sprintf("%.3f", q)))]
  dt[, blank             := paste(rep(" ", 30), collapse = " ")]
  
  plot_dt <- as.data.frame(dt[, .(cohort, `Fusion+/Fusion-`, blank,
                                  `HR (95% CI)`, p_lab, q_lab)])
  colnames(plot_dt) <- c("Cohort", "Fusion", " ",
                         "HR (95% CI)", "p", "q (BH)")
  
  ax <- sym_log_axis(dt$low, dt$up)
  
  note_txt <- if(nrow(res_note) > 0)
    paste0("Not shown (HR = Inf / non-estimable): ",
           paste(res_note$cohort, collapse = "; ")) else ""
  
  footnote_txt <- paste(
    "\n",
    "\n",
    "\n",
    "\n",
    "Per-cohort Cox:\nSurv ~ LSAMP Fusion status.\n",
    "q: BH-FDR within source family (TCGA / TARGET / CPTAC).\n",
    note_txt
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
  n_fn_lines <- str_count(footnote_txt, "\n") + 1
  h <- max(6, nrow(plot_dt) * 0.45 + 3 + n_fn_lines * 0.22)
  
  pdf(out_pdf, width = 10, height = h, onefile = FALSE)
  p <- forest(plot_dt,
              est = dt$HR, lower = dt$low, upper = dt$up,
              ci_column = 3,
              ref_line = 1,
              xlim     = ax$xlim,      # ← 对称边界
              ticks_at = ax$ticks_at,  # ← 整齐刻度
              x_trans  = "log",        # ★ 关键：对数轴
              ci_pch = 15, ci_lwd = 1.5,
              title = title, footnote = footnote_txt, theme = tm)
  p <- insert_text(p, text = "Fusion+/Fusion-", row = 1, col = 2,
                   just = "left",
                   gp = gpar(cex = 0.7, col = "black", fontface = "italic"))
  print(p); dev.off()
}

make_forest(res_ok,
            "GDC-PFS-LSAMP-3'",
            "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/附图4/supFigure4H.pdf")



