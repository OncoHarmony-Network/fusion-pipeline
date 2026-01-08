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
library(IOBR)
library(ggrepel)
library(scales)
library(adjustedCurves)
file_meta <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/数据/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta , sep = "\t") 
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

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(ggplot2); library(scales); library(stringr)
})
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_repel) suppressPackageStartupMessages(library(ggrepel))

## ---------- paths ----------
root    <- "C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5"
in_xlsx <- file.path("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/附表/TableS5.xlsx")
out_png <- file.path(root, "Figure5A.png")
out_pdf <- file.path(root, "Figure5A.pdf")

## ---------- read + tidy ----------
d <- read_excel(in_xlsx) %>%
  rename(gene = `Fusion Gene`, cancer = Cancer, n_pos = `Fusion+`, side = Side,
         n_neg = `Fusion-`, p = p, fdr = `q (BH)`,
         hr = HR, ci = CI, status = Status) %>%
  mutate(
    # ★ 1. side 显式转因子，并给图上要显示的标签
    side = factor(side,
                  levels = c("5'", "3'"),
                  labels = c("5'", "3'")),
    endpoint = ifelse(grepl("^OS", status), "OS", "PFS"),
    log2hr   = log2(hr),
    log2lo   = log2(as.numeric(str_extract(ci, "^[0-9.]+"))),
    log2hi   = log2(as.numeric(str_extract(ci, "[0-9.]+$"))),
    nlog10fdr = -log10(fdr),
    cls = case_when(
      fdr < 0.05 ~ "FDR < 0.05",
      p   < 0.05 ~ "P < 0.05 only",
      TRUE       ~ "not significant"
    ),
    cls = factor(cls, levels = c("FDR < 0.05", "P < 0.05 only", "not significant"))
  )

## ★ 2. 标签按 (gene, side) 分组，否则 3′ 面板会标上 5′ 的点
lab <- d %>%
  filter(fdr < 0.05) %>%
  group_by(gene, side) %>%
  slice_max(order_by = nlog10fdr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(fdr)

## ★ 3. 统计注释改成每个分面各一条
ann <- d %>%
  group_by(side) %>%
  summarise(n_tests  = n(),
            n_genes  = n_distinct(gene),
            n_fdr    = sum(fdr < 0.05),
            n_gene_f = n_distinct(gene[fdr < 0.05]),
            n_nom    = sum(p < 0.05),
            .groups  = "drop") %>%
  mutate(txt = sprintf("%d comparisons | %d genes\n%d with FDR < 0.05 (%d genes)\n%d with P < 0.05",
                       n_tests, n_genes, n_fdr, n_gene_f, n_nom))

## ---------- plotting ----------
p <- ggplot(d, aes(log2hr, nlog10fdr)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_point(data = filter(d, cls == "not significant"),
             aes(size = n_pos), shape = 21, fill = "grey80", colour = "grey55",
             stroke = 0.2, alpha = 0.75) +
  geom_point(data = filter(d, cls == "P < 0.05 only"),
             aes(size = n_pos, shape = endpoint), fill = "#EF8A62", colour = "grey30",
             stroke = 0.2, alpha = 0.9) +
  geom_point(data = filter(d, cls == "FDR < 0.05"),
             aes(size = n_pos, shape = endpoint), fill = "#B2182B", colour = "black",
             stroke = 0.3) +
  geom_label(data = ann, aes(x = -Inf, y = Inf, label = txt),
             inherit.aes = FALSE, hjust = -0.06, vjust = 1.35,
             size = 4, lineheight = 1.05, colour = "grey20",
             fill = "white", linewidth = 0, alpha = 0.9,
             label.padding = unit(0.18, "lines")) +
  scale_shape_manual(values = c(OS = 21, PFS = 24), name = "Endpoint") +
  scale_size_continuous(range = c(0.8, 4.2), name = "Fusion-positive\npatients",
                        breaks = c(3, 5, 10)) +
  scale_x_continuous(breaks = log2(c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32)),
                     labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8", "16", "32")) +
  labs(x = "Hazard ratio (fusion-positive vs fusion-negative, log scale)",
       y = expression(-log[10]~"(FDR)")) +
  facet_wrap(~ side, ncol = 2) +          # ★ 4. 分面，这就是 5′/3′ 的呈现方式
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.25),
        legend.position = "right",
        panel.spacing.x = unit(3, "lines"),
        panel.spacing.y = unit(2, "lines"),
        legend.key.height = unit(1.1, "lines"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(size = 20, face = "bold"),
        text = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14),
        plot.margin = margin(14, 20, 14, 14))

if (nrow(lab) > 0) {
  if (has_repel) {
    set.seed(1)
    p <- p + geom_text_repel(data = lab, aes(label = gene), size = 6,
                             fontface = "italic", min.segment.length = 0,
                             segment.size = 0.5, segment.colour = "grey40",
                             box.padding = 0.45, point.padding = 0.3,
                             max.overlaps = Inf, seed = 1)
  } else {
    p <- p + geom_text(data = lab, aes(label = gene), size = 6,
                       fontface = "italic", vjust = -0.9, check_overlap = TRUE)
  }
}
p
ggsave(out_png, p, width = 12, height = 6, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 12, height = 6, device = cairo_pdf, bg = "white")


#Figure5B
gene_name <- "FGFR3"
cohort_choose <- c("EGAD00001003977", "EGAD00001007575")
fusion_select <- fusion_all |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_choose ) |> 
  select(Sample_ID, gene1, gene2,cohort)|>
  rename(Run = Sample_ID) |>
  filter(gene1 %in% gene_name ) |>
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
                         "Sex", "Run", "Treatment","Response2")
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

  OS_all <- rbind(Need3,OS_all) |> 
            filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
  
}
fusion_data_subset <- OS_all |>
  mutate(event = case_when(
    Run %in% fusion_select$Run  ~  "FGFR3 Fusion+",
    TRUE ~ "FGFR3 Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "FGFR3 Fusion+") ~ "FGFR3 Fusion+",
                           TRUE ~ event))  |>
  ungroup()|>
  select(Patient_ID, event, OS_Time, OS_Status,Sex,  Response2,cohort) |>
  distinct() 

data1 <- fusion_data_subset |>
  group_by(event) |>
  summarise(count_event = n()) |>
  ungroup() |> 
  mutate(total_count = sum(count_event), 
         percentage1 = (count_event / total_count) * 100)

data2 <- fusion_data_subset |>
  group_by(event, Sex) |>
  summarise(count_event = n()) |>
  ungroup() |> 
  mutate(total_count = sum(count_event), 
         percentage2 = (count_event / total_count) * 100) 
data2$Sex2 <- c("F1","M1","F2", "M2")
data2$Sex2 <-  factor(data2$Sex2, levels = c("F2", "M2","F1","M1"))

colors <- c(
  "LSAMP Fusion+" = "#f28680",
  "LSAMP Fusion-" = "#5861ac",
  "F" = "#ff7f00",
  "M" = "#ffe0c1"
)
p <- ggplot() + 
  geom_col(data = data1,
           aes(x = 3,
               y= percentage1, 
               fill= event),
           width= 2, 
           color= 'white') +
  geom_col(data = data2,
           aes(x = 5, 
               y= percentage2, 
               fill= Sex2),
           width= 1, 
           color= 'white') +
  coord_polar(theta = "y") +
  scale_fill_manual(
    name = "Sex",
    values = c(
      "FGFR3 Fusion+" = "#f28680",
      "FGFR3 Fusion-" = "#5861ac",
      "F1" = "#ff7f00",
      "F2" = "#ff7f00",
      "M1" = "#ffe0c1",
      "M2" = "#ffe0c1"),
    labels = c(
      "FGFR3 Fusion+" = "FGFR3 Fusion+",
      "FGFR3 Fusion-" = "FGFR3 Fusion-",
      "F1" = "F",
      "F2" = "F",
      "M1" = "M",
      "M2" = "M")
  ) +
  coord_polar(theta = "y") +
  xlim(c(0, 5.5)) +
  theme_void()+ 
  theme(
    legend.position = "right",
    strip.text.x = element_text(size = 14), 
    legend.title = element_text(size = 15), 
    legend.text = element_text(size = 14)
  )
p 


ggsave(filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5B.png"),
       plot = p,         
       width = 7,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "png",           
       bg = "white")
ggsave(filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5B.pdf"),
       plot = p,         
       width = 7,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "pdf",           
       bg = "white")

  #Figure5C
  plot_os_by_group(data = fusion_data_subset, 
                   group_levels = c( "FGFR3 Fusion-","FGFR3 Fusion+"),
                   title = paste0("ICB-BLCA-OS-FGFR3-5'"), 
                   choose = "OS",
                   filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5C.pdf")
  )
  plot_os_by_group(data = fusion_data_subset, 
                   group_levels = c( "FGFR3 Fusion-","FGFR3 Fusion+"),
                   title = paste0("ICB-BLCA-FGFR3-5'"), 
                   choose = "OS",
                   filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5C.png")
  )

  #Figure5D
  fusion_subset <- fusion_data_subset |>
    select(Patient_ID, Response2, event) |>
    rename(Response = Response2) |> 
    filter(Response != "" & !is.na(Response ) & Response != "not collected" ) |> 
    distinct() |>
    as.data.table()
  fusion_subset$Response  <- gsub("Pre-", "", fusion_subset$Response )
  fusion_subset$Response <- factor(fusion_subset$Response, levels = c("NR", "R"))
  fusion_subset$event    <- factor(fusion_subset$event,    levels = c("FGFR3 Fusion+", "FGFR3 Fusion-"))
  tab <- table(fusion_subset$event, fusion_subset$Response)
  stopifnot(all(dim(tab) == c(2, 2)))
  
  a  <- tab["FGFR3 Fusion+", "R"];  b  <- tab["FGFR3 Fusion+", "NR"]
  cc <- tab["FGFR3 Fusion-", "R"];  dd <- tab["FGFR3 Fusion-", "NR"]
  
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
    title           = paste0("ICB-BLCA-FGFR3-5'"),
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
  
  filename <- paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5D.png")
  ggsave(filename = filename, 
         plot = p,  
         width = 5, 
         height = 6,
         dpi = 600)
  
  filename <- paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5D.pdf")
  ggsave(filename = filename, 
         plot = p,  
         width = 5, 
         height = 6,
         dpi = 600)
  
#Figure5E
gene_name <- "LSAMP"
cohort_choose <- c("PHS000452_PD1","PHS000452_CTLA4")
fusion_select <- fusion_all |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_choose ) |> 
  select(Sample_ID, gene1, gene2,cohort)|>
  rename(Run = Sample_ID) |>
  filter(gene2 %in% gene_name) |>
  as.data.table() |>
  distinct()


OS_all <- data.table()
for(cohort_id in cohort_choose){
  Need <- fread(paste0("C:/Users/97245/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv"))
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
  if(cohort_id == "PHS000452_PD1"){
    Need1$primary <- Need1$Primary_Type
    Need1$LDH <- Need1$LDH_Elevated
  }
  if(cohort_id == "PHS000452_CTLA4"){
    Need1$M_Stage <- Need1$M
  }
  columns_to_select <- c("Patient_ID",
                         "OS_Time", "OS_Status",  
                         "PFS_Time", "PFS_Status",
                         "LDH", "Run", "Treatment", 
                         "primary","M_Stage",
                         "Age", "Sex", "Response2")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    select(all_of(existing_columns)) |>  
    filter(!is.na(OS_Time) & !is.na(OS_Status))|>  
    filter(!is.na(PFS_Time) & !is.na(PFS_Status))
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  if (length(missing_columns) > 0) {
    for (col in missing_columns) {
      Need2[[col]] <- NA
    }
  }
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort <- cohort_id
  OS_all <- rbind(Need2,OS_all) |> 
    filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) |> 
    filter(!grepl("On", Response2)) 
}

  
fusion_data_subset <- OS_all |>
  mutate(event = case_when(
    Run %in% fusion_select$Run  ~  "LSAMP Fusion+",
    TRUE ~ "LSAMP Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = if (any(event == "LSAMP Fusion+")) {
    "LSAMP Fusion+"
  } else {
    first(event)          
  }) |>
  group_by(cohort) |>
  filter(n_distinct(event) == 2) |>
  ungroup()  |>
  select(Patient_ID, event, OS_Time, OS_Status, PFS_Time, PFS_Status, Sex, Response2,cohort) |>
  distinct() 

data1 <- fusion_data_subset |>
  group_by(event) |>
  summarise(count_event = n()) |>
  ungroup() |> 
  mutate(total_count = sum(count_event), 
         percentage1 = (count_event / total_count) * 100)

data2 <- fusion_data_subset |>
  group_by(event, Sex) |>
  summarise(count_event = n()) |>
  ungroup() |> 
  mutate(total_count = sum(count_event), 
         percentage2 = (count_event / total_count) * 100) 
data2$Sex2 <- c("F1","M1","F2", "M2")
data2$Sex2 <-  factor(data2$Sex2, levels = c("F2", "M2","F1","M1"))

colors <- c(
  "LSAMP Fusion+" = "#f28680",
  "LSAMP Fusion-" = "#5861ac",
  "F" = "#ff7f00",
  "M" = "#ffe0c1"
)

p <- ggplot() + 
  geom_col(data = data1,
           aes(x = 3,
               y= percentage1, 
               fill= event),
           width= 2, 
           color= 'white') +
  geom_col(data = data2,
           aes(x = 5, 
               y= percentage2, 
               fill= Sex2),
           width= 1, 
           color= 'white') +
  coord_polar(theta = "y") +
  scale_fill_manual(
    name = "Sex",
    values = c(
      "LSAMP Fusion+" = "#f28680",
      "LSAMP Fusion-" = "#5861ac",
      "F1" = "#ff7f00",
      "M1" = "#ffe0c1",
      "F2" = "#ff7f00",
      "M2" = "#ffe0c1"),
    labels = c(
      "LSAMP Fusion+" = "LSAMP Fusion+",
      "LSAMP Fusion-" = "LSAMP Fusion-",
      "F1" = "F",
      "F2" = "F",
      "M1" = "M",
      "M2" = "M")
  ) +
  coord_polar(theta = "y") +
  xlim(c(0, 5.5)) +
  theme_void()+ 
  theme(
    legend.position = "right",
    strip.text.x = element_text(size = 14), 
    legend.title = element_text(size = 15), 
    legend.text = element_text(size = 14)
  )
p 
ggsave(filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5E.png"),
       plot = p,         
       width = 7,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "png",           
       bg = "white")
ggsave(filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5E.pdf"),
       plot = p,         
       width = 7,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "pdf",           
       bg = "white")


#Figure5F
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("ICB-SKCM-LSAMP-3'"), 
                 choose = "OS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5F.pdf")
)
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("ICB-SKCM-LSAMP-3'"), 
                 choose = "OS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5F.png")
)

#Figure5G
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("ICB-SKCM-LSAMP-3'"), 
                 choose = "PFS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5G.pdf")
)
plot_os_by_group(data = fusion_data_subset, 
                 group_levels = c( "LSAMP Fusion-","LSAMP Fusion+"),
                 title = paste0("ICB-SKCM-LSAMP-3'"), 
                 choose = "PFS",
                 filename = paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5G.png")
)

#Figure5H
fusion_data_subset$event <- factor(fusion_data_subset$event, levels = c("LSAMP Fusion+", "LSAMP Fusion-"))
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
  title           = paste0("ICB-SKCM-LSAMP-3'"),
  ggtheme = theme_classic()
) + scale_fill_manual(values = c(
  "R" = "#9dd3af",
  "NR" = "#f2f2f2"
)) + 
  labs(
    subtitle = p_text, x = NULL, y = NULL) +  
  theme(
    text = element_text( size = 14),
    title = element_text(hjust = 0.5, size = 20, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold")
  )


filename <- paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5H.png")
ggsave(filename = filename, 
       plot = p,  
       width = 5, 
       height = 6,
       dpi = 600)

filename <- paste0("C:/Users/97245/Nutstore/1/ShinyIO协作组/IO-GeneFusion/BMC Medicine-R3/图和附图/Figure5/Figure5H.pdf")
ggsave(filename = filename, 
       plot = p,  
       width = 5, 
       height = 6,
       dpi = 600)
