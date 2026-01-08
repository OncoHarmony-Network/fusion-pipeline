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
#Figure3a
#use_condaenv("myenv", required = TRUE)
use_condaenv("C:/Users/Administrator/miniconda3/envs/myenv", required = TRUE)
##use_python("C:/Users/Administrator/miniconda3/envs/myenv/python.exe", required = TRUE)
#use_python("C:/Users/Administrator/AppData/Local/Programs/Python/Python313/python.exe", required = TRUE)
metafusion_high <- fread("C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv" , sep = "\t") |> 
  filter( sample_type %in% c("Tumor","Normal","Blood"))
tumor <-  metafusion_high |>
  filter(sample_type == "Tumor" | sample_type == "Blood") |>
  select(c(Fusion_ID, Sample_ID, cohort))
normal <- metafusion_high |>
  filter(sample_type == "Normal")  |>
  select(c(Fusion_ID, Sample_ID, cohort))  

all_data_Run_pro <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/sample_id.tsv") |>
  rename(Sample_ID = Run , 
         cohort = "cohort_id")

### normal_high
all_normal <- all_data_Run_pro |>
  filter(cancer_type == "NORMAL") |>
  mutate(
    gene_name = case_when(
      Sample_ID %in% normal$Sample_ID & cohort %in% normal$cohort ~ "Fusion+",
      TRUE ~ "Fusion-"
    )
  )

result_normal <- all_normal |>
  group_by(cancer_type) |>
  summarise(
    detected_runs = sum(gene_name == "Fusion+"),
    total_runs = n(),
    detection_rate = round(detected_runs / total_runs, 3)
  ) 
### tumor_high
all_tumor <- all_data_Run_pro |>
  filter(cancer_type != "NORMAL") |>
  mutate(
    gene_name = case_when(
      Sample_ID %in% tumor$Sample_ID & cohort %in% tumor$cohort ~ "Fusion+",
      TRUE ~ "Fusion-"
    )
  )

result_tumor <- all_tumor |>
  group_by(cancer_type) |>
  summarise(
    detected_runs = sum(gene_name == "Fusion+"),
    total_runs = n(),
    detection_rate = round(detected_runs / total_runs, 3)
  ) 
result_all <- rbind(result_normal, result_tumor)
cancer_type <- result_all$cancer_type
for(cancer_type1 in cancer_type){
  ids <- c("Root1", "Root", "nodetected", "detected")
  parents <- c("", "Root1", "Root", "Root")
  colors <- c("#ffffff", "#ffffff","#808080", "#32037d")
  data <- result_all |> filter(cancer_type == cancer_type1)
  values <- c(data$total_runs, data$total_runs, data$total_runs - data$detected_runs, data$detected_runs)
  labels <- c("","","","")
  p_figure2a_part1 <- plot_ly(ids = ids, 
                              labels = labels, 
                              parents = parents, 
                              type = 'sunburst',
                              # textfont = list(color = "#ffffff",size = 166),
                              # insidetextorientation = "radial",
                              values = values, 
                              sort = FALSE,
                              branchvalues = 'total',
                              rotation = 90,
                              marker = list(colors = colors))
  save_image(p_figure2a_part1, file = paste0("C:/Users/Administrator/Desktop/Figure3a_data/", cancer_type1,"_Figure2a.pdf"), height = 1000, width = 1000)
  save_image(p_figure2a_part1, file = paste0("C:/Users/Administrator/Desktop/Figure3a_data/", cancer_type1, "_Figure2a.png"), height = 1000, width = 1000)
}
#Figure3b,c code遗失了
#Figure3b
Donut_Chart_data <- readRDS("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Donut_Chart_data.rds")
p_figure3b<-plot_ly(Donut_Chart_data ,
                    ids = ~ids,
                    labels = ~labels,
                    parents= ~parents,
                    type = 'sunburst',
                    marker = list(colors = Donut_Chart_data$colors))
save_image(p_figure3b, file ="C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/dount_icb.png", 
           height = 1000, width = 1000)
#Figure3c
Donut_Chart_data_gdc <- readRDS("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Donut_Chart_data_gdc.rds")
p_figure3c <-plot_ly(Donut_Chart_data_gdc ,
                    ids = ~ids,
                    labels = ~labels,
                    parents= ~parents,
                    type = 'sunburst',
                    marker = list(colors = Donut_Chart_data_gdc$colors))
p_figure3c
save_image(p_figure3c, file ="C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/dount_gdc.png", 
           height = 1000, width = 1000)

#sup2
#OS

file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta) |>
  tidyr::separate_rows(gene1, gene1_type, sep = ",") |>
  tidyr::separate_rows(gene2, gene2_type, sep = ",") 
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv")
#icb OS Figure3c
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


#icb OS Figure3c
Tumor_ICB_OS <- fusion_all |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_sum_OS) |> 
  select(Sample_ID, gene1, gene2,cohort)|>
  rename(Run = Sample_ID) |>
  distinct()
cancer_os <-  filter(cohort1, OS == TRUE )$cancer |> unique()
cancer_os <- c(cancer_os,"ALL")
for(cancer_type1 in cancer_os){
  if(cancer_type1 == "ALL"){
    cohort_ALL<- cohort1 |> filter(OS == TRUE ) 
  }else{
    cohort_ALL<- cohort1 |> filter(OS == TRUE & cancer_type == cancer_type1) 
  }
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
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup2/",cancer_type1,"_OS",".png"), plot = p2, dpi = 600, width = 6, height = 5) 
}
cohort_os <-  filter(cohort1, OS == TRUE )$id |> unique()
for(cohort_type1 in cohort_os){
  
  cohort_select <- cohort_type1
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
  if(ln == 0 | hn ==0){
    next
  }else{
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
                     title = paste0("ICB ",cohort_type1, " OS"),
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
    ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup2/cohort_survival/",cohort_type1,"_OS",".png"), plot = p2, dpi = 600, width = 6, height = 5) 
  }
}

#PFS icb figure3d
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
cancer_pfs <- c(cancer_pfs,"ALL")
for(cancer_type1 in cancer_pfs){
  if(cancer_type1 == "ALL"){
    cohort_ALL<- cohort1 |> filter(PFS == TRUE ) 
  }else{
    cohort_ALL<- cohort1 |> filter(PFS == TRUE & cancer_type == cancer_type1) 
  }
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
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup2/",cancer_type1,"_PFS",".png"), plot = p2, dpi = 600, width = 6, height = 5) 
}
cohort_pfs <-  filter(cohort1, PFS == TRUE )$id |> unique()
for(cohort_type1 in cohort_pfs){
  cohort_select <- cohort_type1
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
  if(hn == 0 | ln == 0){
    next
  }else{
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
                     title = paste0("ICB ",cohort_type1, " PFS"),
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
    ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup2/cohort_survival/",cohort_type1,"_PFS",".png"), plot = p2, dpi = 600, width = 6, height = 5) 
  }
}


#fig2m
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta , sep = "\t") |>
  separate_rows(gene1, sep = ",") |>
  separate_rows(gene2, sep = ",")
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv") |> 
  filter(Response2 == TRUE)

cohort_select2 <-  cohort1$id |> unique()
Response_all <- data.table()
for(cohort_id in cohort_select2){
  Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv"))
  if(filter(cohort1,id == cohort_id)$time == "Pre/On" ){
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
  
  columns_to_select <- c("Patient_ID", "Run", "Response2", "Treatment")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    select(all_of(existing_columns)) 
  missing_columns <- columns_to_select[!columns_to_select %in% colnames(Need2)]
  if (length(missing_columns) > 0) {
    for (col in missing_columns) {
      Need2[[col]] <- NA
    }
  }
  matched_indices <- match(cohort_id , cohort1$id)
  matched_cancer <- cohort1$cancer[matched_indices]
  Need2$cancer <- matched_cancer
  Need2$Patient_ID <- as.character(Need2$Patient_ID)
  Need2$cohort <- cohort_id
  Response_all <- rbind(Need2,Response_all) |> 
    filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
}
Tumor_ICB <- fusion_all |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_select2) |> 
  select(Sample_ID, gene1, gene2,cohort)|>
  mutate(fusion = paste0(gene1, "::", gene2)) |>
  select(-c(gene1, gene2)) |> 
  rename(Run = Sample_ID) |>
  distinct()
Response_all1 <- Response_all |>
  filter(Response2 != "") |>
  filter(!grepl("On", Response2)) |>
  mutate(Fusion = case_when(
    Run %in% Tumor_ICB$Run  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(Fusion = case_when(any(Fusion == "Fusion+") ~ "Fusion+",
                            TRUE ~ Fusion)) |>
  select(Patient_ID, Response2, Fusion, cancer, cohort) |>
  distinct()  
cancer_all <- Response_all1$cancer |> unique() |> setdiff("BRCA")
cancer_all <- c(cancer_all, "ALL")
for (cancer1 in cancer_all){
  if(cancer1 == "ALL"){
    fusion_subset <- Response_all1 |> 
      as.data.table()
  }else{
    fusion_subset <- Response_all1 |> 
      filter(cancer == cancer1) |> 
      as.data.table() |>
      rename(Response = Response2)
  }
  fusion_subset$Response <- gsub("Pre-", "", fusion_subset$Response )
  fusion_subset[, `:=`(
    Fusion = factor(Fusion, levels = c("Fusion+", "Fusion-")),
    Response  = factor(Response,  levels = c("NR", "R"))
  )]
  
  p <- ggstatsplot::ggbarstats(
    data            = fusion_subset,
    x               = Response,      
    y               =  Fusion,     
    paired          = FALSE,      
    results.subtitle = TRUE,  
    label           = "both",     
    title           = "Fusion status vs irAE"
  )
  subtitle <- p$labels$subtitle
  p_value <- regmatches(subtitle, regexec('italic\\(p\\) == "([0-9.]+)"', subtitle))[[3]][3]
  p <- ggstatsplot::ggbarstats(
    data            = fusion_subset,
    x               = Response,     
    y               = Fusion,      
    paired          = FALSE,      
    results.subtitle = FALSE,  
    label           = "both",     
    title           = cancer1
  )
  p <- p + labs(
    subtitle = paste0("p: ",  p_value)
  )+
    ggthemes::theme_base() +
    scale_fill_manual(
      values = c(
        "R" = "#9dd3af",
        "NR" = "#f2f2f2")
    ) 
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup2/fusion and response/",cancer1,"_fusion_response.pdf"), 
         plot = p, dpi = 600, width = 5, height = 7)
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup2/fusion and response/",cancer1,"_fusion_response.png"), 
         plot = p, dpi = 600, width = 5, height = 7)
}
cohort_all <- Response_all1$cohort |> unique() |> setdiff(c("PRJNA312948", "PHS003316",
                                                            "PHS001919","PHS001427", "PHS001038"))
for (cohort_name in cohort_all){
  
  fusion_subset <- Response_all1 |> 
    filter(cohort == cohort_name) |> 
    as.data.table()
  
  fusion_subset[, `:=`(
    Fusion = factor(Fusion, levels = c("Fusion+", "Fusion-")),
    Response2  = factor(Response2,  levels = c("Pre-NR", "Pre-R"))
  )]
  
  p <- ggstatsplot::ggbarstats(
    data            = fusion_subset,
    x               = Response2,      
    y               =  Fusion,     
    paired          = FALSE,      
    results.subtitle = TRUE,  
    label           = "both",     
    title           = "Fusion status vs irAE"
  )
  subtitle <- p$labels$subtitle
  p_value <- regmatches(subtitle, regexec('italic\\(p\\) == "([0-9.]+)"', subtitle))[[3]][3]
  p <- ggstatsplot::ggbarstats(
    data            = fusion_subset,
    x               = Response2,     
    y               = Fusion,      
    paired          = FALSE,      
    results.subtitle = FALSE,  
    label           = "both",     
    title           = "Fusion event vs Response"
  )
  p <- p + labs(
    subtitle = paste0("p: ",  p_value)
  )+
    ggthemes::theme_base() +
    scale_fill_manual(
      values = c(
        "Pre-R" = "#9dd3af",
        "Pre-NR" = "#f2f2f2")
    ) 
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup2/fusion cohort and response/",cohort_name,"_fusion_response.pdf"), 
         plot = p, dpi = 600, width = 5, height = 6)
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup2/fusion cohort and response/",cohort_name,"_fusion_response.png"), 
         plot = p, dpi = 600, width = 5, height = 6)
}





