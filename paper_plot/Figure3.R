library(ggplot2)
library(survival)
library(survminer)
library(data.table)
library(dplyr)
library(grid)
library(tidyr)
#Sys.setenv(RETICULATE_PYTHON = "/c/Users/Administrator/AppData/Local/Programs/Python/Python313/python")
library(reticulate)
library(plotly)
library(ggpubr)
library(ggthemes)
library(ggprism)
library(ggsignif)
library(stringr)
library(patchwork)
library(forestploter)
library(reticulate)
library(plotly)

file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta) |>
  tidyr::separate_rows(gene1, gene1_type, sep = ",") |>
  tidyr::separate_rows(gene2, gene2_type, sep = ",") 
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv")
#Figure3A
cohort_OS <- cohort1 |> filter(OS == TRUE) 
cohort_sum_OS <- (cohort_OS$id) |> unique()
OS_all <- data.table()
for(cohort_id in cohort_sum_OS){
  Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv"))
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
    }}else if(cohort_id == "PRJNA744780"){
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
  
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort <- cohort_id
  OS_all <- rbind(Need2,OS_all) |> filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
}
Tumor_ICB_OS <- fusion_all |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_sum_OS) |> 
  select(Sample_ID, gene1, gene2,cohort)|>
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
    select(Patient_ID, event, OS_Time, OS_Status) |>
    distinct() 
  hn = nrow(filter(fusion_data_subset, event == "Fusion+"))
  ln = nrow(fusion_data_subset) - hn 
  fusion_data_subset$event <- factor(fusion_data_subset$event, level = c("Fusion-", "Fusion+"))
  data.survdiff <- survdiff(Surv(OS_Time, OS_Status) ~ event, data = fusion_data_subset)
  fit <- survfit(Surv(OS_Time, OS_Status) ~ event, data = fusion_data_subset)
  p.val1 <- round(1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1), 3)
  p.val <- ifelse(p.val1 < 0.001, "\np< 0.001", paste0("\np:",p.val1))
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")
  
  caption1 <- paste0(
    paste0("Fusion+ :"), hn,
    paste0("\n", "Fusion- :"), ln,
    p.val,
    "\n", HR, "\n", CI
  )
  p1 <- ggsurvplot(fit,
                   linetype = "strata",
                   risk.table = TRUE,
                   risk.table.col = "strata",
                   legend.title = "",
                   pval = FALSE,      
                   font.y = c(20, "plain", "black"),    
                   font.legend = c(15, "plain", "black"), 
                   font.tickslab = c(15, "plain", "black"), 
                   risk.table.fontsize  = 4, 
                   palette = c( "#00468bff", "#ed0000ff"),
                   title = paste0("ICB ",cancer_type1, " OS"),
                   legend.labs = c( "Fusion-", "Fusion+"),
                   xlab = "Time(Days)"
  ) 
  p1$plot <- p1$plot + 
    annotate("text", Inf, Inf, label = caption1,  vjust = 1, hjust = 1,size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 4)
    ) 
  
  p1$table$labels$title <- ""
  p1$table$theme$legend.position <- "none"
  
  p1$table <- p1$table +
    theme(
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
    )
  p2 <- (p1$plot / p1$table) +
    plot_layout(heights = c(4, 1)) +
    plot_annotation(
      theme = theme(
        panel.spacing = unit(0, "lines")  ,
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
      )
    )
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/Figure2A.png"),
   plot = p2, dpi = 600, width = 6, height = 5) 

#Figure3B
cohort_PFS <- cohort1 |> filter(PFS == TRUE) 
cohort_sum_PFS <- cohort_PFS$id |> unique()
PFS_all <- data.table()
for(cohort_id in cohort_sum_PFS){
  Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv"))
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
    }}else if(cohort_id == "PRJNA744780"){
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
  dplyr::select(Sample_ID, gene1, gene2,cohort)|>
  dplyr::rename(Run = Sample_ID) |>
  distinct()
cancer_pfs <-  filter(cohort1, PFS == TRUE )$cancer |> unique()
cancer_type1 = "ALL"
cohort_ALL<- cohort1 |> filter(PFS == TRUE ) 
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
    select(Patient_ID, event, PFS_Time, PFS_Status) |>
    distinct() 
  hn = nrow(filter(fusion_data_subset, event == "Fusion+"))
  ln = nrow(fusion_data_subset) - hn 
  fusion_data_subset$event <- factor(fusion_data_subset$event, level = c("Fusion-", "Fusion+"))
  data.survdiff <- survdiff(Surv(PFS_Time, PFS_Status) ~ event, data = fusion_data_subset)
  fit <- survfit(Surv(PFS_Time, PFS_Status) ~ event, data = fusion_data_subset)
  p.val1 <- round(1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1), 3)
  p.val <- ifelse(p.val1 < 0.001, "\np< 0.001", paste0("\np:",p.val1))
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")
  
  caption1 <- paste0(
    paste0("Fusion+ :"), hn,
    paste0("\n", "Fusion- :"), ln,
    p.val,
    "\n", HR, "\n", CI
  )
  p1 <- ggsurvplot(fit,
                   linetype = "strata",
                   risk.table = TRUE,
                   risk.table.col = "strata",
                   legend.title = "",
                   pval = FALSE,      
                   font.y = c(20, "plain", "black"),    
                   font.legend = c(15, "plain", "black"), 
                   font.tickslab = c(15, "plain", "black"), 
                   risk.table.fontsize  = 4, 
                   palette = c( "#00468bff", "#ed0000ff"),
                   title = paste0("ICB ",cancer_type1, " PFS"),
                   legend.labs = c( "Fusion-", "Fusion+"),
                   xlab = "Time(Days)"
  ) 
  p1$plot <- p1$plot + 
    annotate("text", Inf, Inf, label = caption1,  vjust = 1, hjust = 1,size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 4)
    ) 
  
  p1$table$labels$title <- ""
  p1$table$theme$legend.position <- "none"
  
  p1$table <- p1$table +
    theme(
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
    )
  p2 <- (p1$plot / p1$table) +
    plot_layout(heights = c(4, 1)) +
    plot_annotation(
      theme = theme(
        panel.spacing = unit(0, "lines")  ,
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
      )
    )
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/Figure2B.png"), plot = p2, dpi = 600, width = 6, height = 5) 

#Figure3F,G
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
metafusion <- fread(file_meta , sep = "\t")
Credibility_gene <- metafusion |>
  separate_rows(gene1, sep = ",") |>
  separate_rows(gene2, sep = ",")  |> 
  dplyr::filter(sample_type == "Tumor") |>
  mutate(pair = paste0(gene1,"::",gene2))
folder_path <- "C:/Users/Administrator/Desktop/fsdownload/Clininfo_GDC/"
file_names <- list.files(folder_path, full.names = FALSE, pattern = "_info.tsv$")


#Figure3C gdc pfs
cohort_sum_PFS <- file_names[!grepl("^TARGET", file_names)]
cohort_sum_PFS <- sub("_info.tsv$", "", cohort_sum_PFS ) |>
  setdiff(c("CPTAC-NORMAL","TARGET-NORMAL", "TCGA-NORMAL", "TCGA-LAML"))
cohort_sum_tcga   <- cohort_sum_PFS[grepl("^TCGA", cohort_sum_PFS)]
cohort_sum_cptac  <- cohort_sum_PFS[grepl("^CPTAC", cohort_sum_PFS)] 
cohort_sum_PFS1 <- c("TCGA-ALL", "CPTAC-ALL", cohort_sum_PFS)
data_all_pfs <- data.table()
for(cohort_id in cohort_sum_PFS1){
cohort_id2 <- gsub("_", "-", cohort_id)
if(cohort_id == "TCGA-ALL"){
  Need_all <- data.table()
  for(cohort_id3 in cohort_sum_tcga){
    Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id3,"_info.tsv"))
    Need_all <- rbind(Need_all, Need)
  }
  fusion_select <- Credibility_gene |> 
    dplyr::filter(cohort %in% cohort_sum_tcga)|> 
    as.data.table()
}else if(cohort_id == "CPTAC-ALL"){
  Need_all <- data.table()
  for(cohort_id3 in cohort_sum_cptac){
    Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id3,"_info.tsv"))
    Need_all <- rbind(Need_all, Need)
  }
  fusion_select <- Credibility_gene |> 
    dplyr::filter(cohort %in% cohort_sum_cptac)|> 
    as.data.table()
}else{
  Need_all <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id2,"_info.tsv"))
  fusion_select <- Credibility_gene |> 
    dplyr::filter(cohort == cohort_id2)|> 
    as.data.table()
}

fusion_data_subset <-   Need_all |> 
  mutate(event = case_when(
    Run %in% fusion_select$Sample_ID  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                           TRUE ~ event))  |>
  select(Patient_ID, event, PFS_Time, PFS_Status) |>
  distinct() 
hn = nrow(filter(fusion_data_subset, event == "Fusion+"))
ln = nrow(fusion_data_subset) - hn 
fusion_data_subset$event <- factor(fusion_data_subset$event, level = c("Fusion-", "Fusion+"))
data.survdiff <- survdiff(Surv(PFS_Time, PFS_Status) ~ event, data = fusion_data_subset)
fit <- survfit(Surv(PFS_Time, PFS_Status) ~ event, data = fusion_data_subset)
p.val1 <- round(1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1), 3)
p.val <- ifelse(p.val1 < 0.001, "p < 0.001", p.val1)
HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
data <- data.table(cohort = cohort_id, `Fusion+` = hn, 
                   `Fusion-` = ln, p = p.val, HR = HR, CIlow = low95, CIup = up95)
data_all_pfs  <- rbind(data_all_pfs, data)
}
fwrite(data_all_pfs, "C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3C_data.tsv")

tm <- forest_theme(
  base_size = 10,  
  colhead = list(fg_params = list(hjust = 0)),  
  core = list(padding = unit(c(5,5), "mm")),  
  
  ci_pch = 15,  
  ci_col = "black",  
  ci_fill = "black",  
  ci_alpha = 0.8, 
  ci_lty = 1, 
  ci_lwd = 1.5, 
  ci_Theight = 0.2,  
  
  refline_gp = gpar(
    lwd = 1,  
    lty = "dashed",  
    col = "grey20"  
  ),
  footnote_gp = gpar(
    cex = 0.6, 
    fontface = "plain", 
    col = "black"  
  ),
  
  vertline_lwd = 1,  
  vertline_lty = "dashed",  
  vertline_col = "grey20",  
  
)
data_all_pfs <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3C_data.tsv")
data_all_pfs$HR <- round(data_all_pfs$HR,2) |> as.numeric()
data_all_pfs$CIlow <- round(data_all_pfs$CIlow,2) |> as.numeric()
data_all_pfs$CIup <- round(data_all_pfs$CIup,2) |> as.numeric()
data_all_pfs1 <-  data_all_pfs |> 
  filter(!is.na(HR)) |> 
  as.data.frame() |>
  mutate(Fusion = paste0(`Fusion+`, "/",`Fusion-`),
         `HR(95%CI)` = paste0(HR, "(",CIlow," - ", CIup,")")) 
data_all_pfs1$` ` <- paste(rep(" ", 20), collapse = " ")
pdf("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3C.pdf", 
width = 7, height = 14)

p <- forest(
  data_all_pfs1[,c(1,8,10,4,9)],
  title = "GDC PFS",
  est = data_all_pfs1$HR,
  lower = data_all_pfs1$CIlow,
  upper = data_all_pfs1$CIup,
  ci_column = 3,
  ref_line = 1,
  xlim = c(0, 3),
  ticks_at = c(0, 0.5, 1, 2, 3),
  footnote = "\n\n\nCPTAC-BRCA (100/33),\nCPTAC-COAD (100/330),\nCPTAC-OV (87/14) p=1, not shown.\nIncluded in CPTAC-ALL.",
  theme = tm
)
p <- insert_text(
  p,
  text = "Fusion+/Fusion-",
  row = 1,
  col = 2,
  just = "left",

  gp = gpar(cex = 0.7, col = "black", fontface = "italic")
)
plot(p)
dev.off()

#figure3D gdc os
cohort_sum_OS <- sub("_info.tsv$", "", file_names) |> setdiff(c("CPTAC-NORMAL","TARGET-NORMAL", "TARGET-CELL","TCGA-NORMAL"))
cohort_sum_tcga   <- cohort_sum_OS[grepl("^TCGA", cohort_sum_OS)]
cohort_sum_cptac  <- cohort_sum_OS[grepl("^CPTAC", cohort_sum_OS)]
cohort_sum_target  <- cohort_sum_OS[grepl("^TARGET", cohort_sum_OS)]
cohort_sum_OS1 <- c("TCGA-ALL", "CPTAC-ALL", "TARGET-ALL",cohort_sum_OS)
data_all_os <- data.table()
for(cohort_id in cohort_sum_OS1){
  cohort_id2 <- gsub("_", "-", cohort_id)
  if(cohort_id == "TCGA-ALL"){
    Need_all <- data.table()
    for(cohort_id3 in cohort_sum_tcga){
      Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id3,"_info.tsv"))
      Need_all <- rbind(Need_all, Need)
    }
    fusion_select <- Credibility_gene |> 
      dplyr::filter(cohort %in% cohort_sum_tcga)|> 
      as.data.table()
  }else if(cohort_id == "CPTAC-ALL"){
    Need_all <- data.table()
    for(cohort_id3 in cohort_sum_cptac){
      Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id3,"_info.tsv"))
      Need_all <- rbind(Need_all, Need)
    }
    fusion_select <- Credibility_gene |> 
      dplyr::filter(cohort %in% cohort_sum_cptac)|> 
      as.data.table()
  }else if(cohort_id == "TARGET-ALL"){
    Need_all <- data.table()
    for(cohort_id3 in cohort_sum_target){
      Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id3,"_info.tsv"))
      Need_all <- rbind(Need_all, Need)
    }
    fusion_select <- Credibility_gene |> 
      dplyr::filter(cohort %in% cohort_sum_target)|> 
      as.data.table()
  }else{
    Need_all <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id2,"_info.tsv"))
    fusion_select <- Credibility_gene |> 
      dplyr::filter(cohort == cohort_id2)|> 
      as.data.table()
  }
  
  fusion_data_subset <-   Need_all |> 
    mutate(event = case_when(
      Run %in% fusion_select$Sample_ID  ~  "Fusion+",
      TRUE ~ "Fusion-"
    )) |>
    group_by(Patient_ID) |>
    mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                             TRUE ~ event))  |>
    select(Patient_ID, event, OS_Time, OS_Status) |>
    filter(!is.na(OS_Time) & !is.na(OS_Status)) |>
    distinct() 
  hn = nrow(filter(fusion_data_subset, event == "Fusion+"))
  ln = nrow(fusion_data_subset) - hn 
  fusion_data_subset$event <- factor(fusion_data_subset$event, level = c("Fusion-", "Fusion+"))
  data.survdiff <- survdiff(Surv(OS_Time, OS_Status) ~ event, data = fusion_data_subset)
  fit <- survfit(Surv(OS_Time, OS_Status) ~ event, data = fusion_data_subset)
  p.val1 <- round(1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1), 3)
  p.val <- ifelse(p.val1 < 0.001, "p < 0.001", p.val1)
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  data <- data.table(cohort = cohort_id, `Fusion+` = hn, 
                     `Fusion-` = ln, p = p.val, HR = HR, CIlow = low95, CIup = up95)
  data_all_os  <- rbind(data_all_os, data)
}
fwrite(data_all_os, "C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3D_data.tsv")

data_all_os <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3D_data.tsv")
data_all_os$HR <- round(data_all_os$HR,2) |> as.numeric()
data_all_os$CIlow <- round(data_all_os$CIlow,2) |> as.numeric()
data_all_os$CIup <- round(data_all_os$CIup,2) |> as.numeric()
data_all_os1 <-  data_all_os |> 
  filter(HR != Inf) |> 
  as.data.frame() |>
  mutate(Fusion = paste0(`Fusion+`, "/",`Fusion-`),
         `HR(95%CI)` = paste0(HR, "(",CIlow," - ", CIup,")")) 
data_all_os1$` ` <- paste(rep(" ", 20), collapse = " ")
pdf("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3C.pdf", 
    width = 7, height = 16)
#    width = 7, height = 14, units = 'in', res = 600)
p <- forest(
  data_all_os1[,c(1,8,10,4,9)],
  title = "GDC OS",
  est = data_all_os1$HR,
  lower = data_all_os1$CIlow,
  upper = data_all_os1$CIup,
  ci_column = 3,
  ref_line = 1,
  xlim = c(0, 3),
  ticks_at = c(0, 0.5, 1, 2, 3),
  footnote = "\n\n\n\nTCGA-PRAD (471/26),\nCPTAC-BRCA (100/33),\nTCGA-TGCT (92/58),\nCPTAC-OV (87/14) HR=Inf, not shown.\nIncluded in TCGA-ALL/CPTAC-ALL.",
  theme = tm
)
p <- insert_text(
  p,
  text = "Fusion+/Fusion-",
  row = 1,
  col = 2,
  just = "left",
  gp = gpar(cex = 0.7, col = "black", fontface = "italic")
)
plot(p)
dev.off()


#Figure3G
NSCLC_OS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_NSCLC.OS.tsv")  
SKCM_OS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_SKCM.OS.tsv") 
GBM_OS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_GBM.OS.tsv") 
STAD_OS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_STAD.OS.tsv") 
BLCA_OS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_BLCA.OS.tsv") 
HNSC_OS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_HNSC.OS.tsv") 

OS_all <- rbind(NSCLC_OS,SKCM_OS,GBM_OS,STAD_OS,BLCA_OS,HNSC_OS) |> 
  mutate(rate = `Fusion+` / (`Fusion+` + `Fusion-`) ) |> 
  filter(rate > 0.01) |> 
  select(-rate)
OS_all$HR <- gsub("Hazard Ratio = ", "", OS_all$HR)
OS_all$CI <- gsub("95% CI: ", "", OS_all$CI)
OS_all$p <- gsub("p:", "", OS_all$p)
OS_all$p <- gsub("p< 0.001", "0.001", OS_all$p)
OS_all$p <- as.numeric(OS_all$p)
OS_all1 <- OS_all |>
  group_by(cancer) |>
  summarise(Fusion_count = n_distinct(Fusion),
            signf = sum(p < 0.05),
            non.signf = Fusion_count - signf) |>
  ungroup() |>
  pivot_longer(cols = c(signf, non.signf),
               names_to = "Fusion_type",
               values_to = "values") |>
  filter(!values == 0) |>
  mutate(
    Status = "OS")

NSCLC_PFS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_NSCLC.PFS.tsv") 
SKCM_PFS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_SKCM.PFS.tsv") 
GBM_PFS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_GBM.PFS.tsv") 
STAD_PFS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_STAD.PFS.tsv") 
KIRC_PFS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_KIRC.PFS.tsv") 
SGC_PFS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/ICB_KM/icb_SGC.PFS.tsv") 

PFS_all <- rbind(NSCLC_PFS,SKCM_PFS,GBM_PFS,STAD_PFS,KIRC_PFS,SGC_PFS)|> 
  mutate(rate = `Fusion+` / (`Fusion+` + `Fusion-`) ) |> 
  filter(rate > 0.01) |> 
  select(-rate)
PFS_all$HR <- gsub("Hazard Ratio = ", "", PFS_all$HR)
PFS_all$CI <- gsub("95% CI: ", "", PFS_all$CI)
PFS_all$p <- gsub("p:", "", PFS_all$p)
PFS_all$p <- gsub("p< 0.001", "0.001", PFS_all$p)
PFS_all$p <- as.numeric(PFS_all$p)
PFS_all1 <- PFS_all |>
  group_by(cancer) |>
  summarise(Fusion_count = n_distinct(Fusion),
            signf = sum(p < 0.05),
            non.signf = Fusion_count - signf) |>
  ungroup() |>
  pivot_longer(cols = c(signf, non.signf),
               names_to = "Fusion_type",
               values_to = "values") |>
  filter(!values == 0) |>
  mutate(Status = "PFS" )
result_all <- rbind(PFS_all1, OS_all1) |> 
  mutate(rate =  (values / Fusion_count) * 100, 
         type = paste0(cancer,"_",Status),
         Fusion_type = factor(Fusion_type, levels = c("signf", "non.signf"))
  ) |> 
  arrange(desc(rate))
result_all <- result_all[c(13:21,4:12,1,2,3),]
result_all$type <-  factor(result_all$type, levels = unique(result_all$type) )

result_all <-   result_all  |> 
  mutate(type_num = row_number()) 
colors <- c(signf = "#e1abbc", non.signf = "#6a73cf")
p <- ggplot(result_all, aes(x = type, y = rate, fill = Fusion_type)) + 
  geom_col(width = 0.55, position = position_stack(reverse = TRUE), show.legend = FALSE) + 
  geom_line(aes(colour = Fusion_type)) +
  scale_x_discrete( limits =   c("SKCM_OS" , "NSCLC_PFS",
                                 "STAD_PFS","BLCA_OS",   "KIRC_PFS",  "SKCM_PFS", 
                                 "GBM_OS",  "NSCLC_OS",    "STAD_OS", 
                                 "GBM_PFS","SGC_PFS",   "HNSC_OS" ), 
                    labels = c("SKCM_OS(142/833)" , "NSCLC_PFS(10/65)", 
                               "STAD_PFS(4/26)", "BLCA_OS(8/55)",  "KIRC_PFS(4/30)",  "SKCM_PFS(34/288)", 
                               "GBM_OS(3/29)",  "NSCLC_OS(4/65)",    "STAD_OS(1/26)", 
                               "GBM_PFS(0/27)",  "SGC_PFS(0/4)",   "HNSC_OS(0/1)" )) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20)) + 
  scale_fill_manual(values = colors) +  
  scale_colour_manual(values = colors) +
  theme_classic(base_size = 14) + 
  theme(axis.line = element_line(size = 0.3), 
        axis.ticks = element_line(size = 0.3), 
        axis.ticks.length = unit(0.2, "cm"), 
        panel.grid = element_blank(), 
        axis.title.y = element_text(size = 10), 
        panel.border = element_blank(), 
        axis.text.x = element_text(colour = "black",angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top") + 
  labs(subtitle = "", 
       colour = "", 
       fill = "", 
       x = "", 
       y = "")
p

ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3E.png"), 
       plot = p, dpi = 600, width = 8, height = 5)
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3E.pdf"), 
       plot = p, dpi = 600, width = 8, height = 5)

#Figure3E
GDC_OS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/GDC_KM/GDC_OS.tsv")|> 
  mutate(Status = "OS")  |> 
  mutate(rate = `Fusion+` / (`Fusion+` + `Fusion-`) ) |> 
  filter(rate > 0.01) |> 
  select(-rate)

GDC_ALL <- rbind(GDC_OS) |>
  group_by(cohort) |>
  summarise(
    Fusion_count = n_distinct(Fusion),
    signf = sum(p < 0.05),
    non.signf = Fusion_count - signf,
    .groups = "drop" 
  ) |>
  pivot_longer(cols = c(signf, non.signf),
               names_to = "Fusion_type",
               values_to = "values") |>
  mutate(rate = (values / Fusion_count) * 100, 
         Fusion_type = factor(Fusion_type, levels = c("signf", "non.signf")),
         type = paste0(cohort, "_OS")) |>
  filter(values != 0) |>
  arrange(desc(rate))
GDC_ALL <- GDC_ALL[c(47:77,16:46,2,14,13,8,12,11,3,15,9,6,4,10,1,5,7),]
GDC_ALL$type <-  factor(GDC_ALL$type, levels = unique(GDC_ALL$type) )

GDC_ALL1 <- GDC_ALL %>%
  mutate(type_num = row_number()) %>%
  mutate(
    labels = case_when(
      Fusion_type == "non.signf" & rate == 100  ~  paste0(type,"(","0/",Fusion_count,")"),
      Fusion_type == "signf"  ~ paste0(type,"(",values,"/",Fusion_count,")"),
      TRUE ~  NA 
    ) ) %>%
  group_by(type) %>%
  mutate(
    type_num = case_when(
      Fusion_type == "non.signf" & any(Fusion_type == "signf") ~ type_num[Fusion_type == "signf"][1],
      TRUE ~ type_num ),
    labels = case_when(
      Fusion_type == "non.signf" & rate != 100  ~ first(labels[Fusion_type == "signf"]),
      TRUE ~  labels
    )
  ) %>%
  ungroup()
colors <- c(signf = "#e1abbc", non.signf = "#6a73cf")
p <- ggplot(GDC_ALL1, aes(x = type, y = rate, fill = Fusion_type)) + 
  geom_col(width = 0.55, position = position_stack(reverse = TRUE), show.legend = FALSE) + 
  geom_line(aes(colour = Fusion_type)) +
  scale_x_discrete(labels = c(unique(GDC_ALL1$labels))) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20)) + 
  scale_fill_manual(values = colors) +  
  scale_colour_manual(values = colors) +
  theme_classic(base_size = 14) + 
  theme(plot.subtitle = element_text(face = "bold", size = 18, hjust = -0.2, vjust = -3), 
        axis.line = element_line(size = 0.3), 
        axis.ticks = element_line(size = 0.3), 
        axis.ticks.length = unit(0.2, "cm"), 
        axis.title.y = element_text(size = 10), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text.x = element_text(colour = "black",angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top") + 
  labs(subtitle = "", 
       colour = "", 
       fill = "", 
       x = "", 
       y = "")
p

ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3F.png"), 
       plot = p, dpi = 600, width = 8, height = 5)
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3F.pdf"), 
       plot = p, dpi = 600, width = 8, height = 5)


#Figure3F
GDC_PFS <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/计算结果/GDC_KM/GDC_PFS.tsv")|> 
  mutate(Status = "PFS")  |> 
  mutate(rate = `Fusion+` / (`Fusion+` + `Fusion-`) ) |> 
  filter(rate > 0.01) |> 
  select(-rate)

GDC_ALL <- GDC_PFS |>
  group_by(cohort) |>
  summarise(
    Fusion_count = n_distinct(Fusion),
    signf = sum(p < 0.05),
    non.signf = Fusion_count - signf,
    .groups = "drop" 
  ) |>
  pivot_longer(cols = c(signf, non.signf),
               names_to = "Fusion_type",
               values_to = "values") |>
  mutate(rate = (values / Fusion_count) * 100, 
         Fusion_type = factor(Fusion_type, levels = c("signf", "non.signf")),
         type = paste0(cohort, "_PFS" )) |>
  filter(values != 0) |>
  arrange(desc(rate))
GDC_ALL <- GDC_ALL[c(42:66,17:41,12,14,11,13,10,7,4,15,16,8,3,6,1,5,9,2),]
GDC_ALL$type <-  factor(GDC_ALL$type, levels = unique(GDC_ALL$type) )

GDC_ALL1 <- GDC_ALL %>%
  mutate(type_num = row_number()) %>%
  mutate(
    labels = case_when(
      Fusion_type == "non.signf" & rate == 100  ~  paste0(type,"(","0/",Fusion_count,")"),
      Fusion_type == "signf"  ~ paste0(type,"(",values,"/",Fusion_count,")"),
      TRUE ~  NA 
    ) ) %>%
  group_by(type) %>%
  mutate(
    type_num = case_when(
      Fusion_type == "non.signf" & any(Fusion_type == "signf") ~ type_num[Fusion_type == "signf"][1],
      TRUE ~ type_num ),
    labels = case_when(
      Fusion_type == "non.signf" & rate != 100  ~ first(labels[Fusion_type == "signf"]),
      TRUE ~  labels
    )
  ) %>%
  ungroup()
colors <- c(signf = "#e1abbc", non.signf = "#6a73cf")
p <- ggplot(GDC_ALL1, aes(x = type, y = rate, fill = Fusion_type)) + 
  geom_col(width = 0.55, position = position_stack(reverse = TRUE), show.legend = FALSE) + 
  geom_line(aes(colour = Fusion_type)) +
  scale_x_discrete(labels = c(unique(GDC_ALL1$labels))) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20)) + 
  scale_fill_manual(values = colors) +  
  scale_colour_manual(values = colors) +
  theme_classic(base_size = 14) + 
  theme(plot.subtitle = element_text(face = "bold", size = 18, hjust = -0.2, vjust = -3), 
        axis.line = element_line(size = 0.3), 
        axis.ticks = element_line(size = 0.3), 
        axis.ticks.length = unit(0.2, "cm"), 
        axis.title.y = element_text(size = 10), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text.x = element_text(colour = "black",angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "top") + 
  labs(subtitle = "", 
       colour = "", 
       fill = "", 
       x = "", 
       y = "")
p

ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3G.png"), 
       plot = p, dpi = 600, width = 8, height = 5)
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure3G.pdf"), 
       plot = p, dpi = 600, width = 8, height = 5)