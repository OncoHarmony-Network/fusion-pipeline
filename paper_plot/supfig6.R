
# Fusion- AND Fusion+ 
# Fusion+ AND low Fusion
#sup4 fusion counts and survival
library(ggplot2)
library(survival)
library(survminer)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)


#icb os
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta , sep = "\t") |>
  separate_rows(gene1, sep = ",") |>
  separate_rows(gene2, sep = ",")
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv") |> 
  filter(Response2 == TRUE)

cohort_select2 <-  cohort1$id |> unique()
OS_all <- data.table()
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
  
  columns_to_select <- c("Patient_ID", "OS_Time", "OS_Status", "Run", "Response2", "Treatment",
                         "OS_Time", "OS_Status", "PFS_Time", "PFS_Status")
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
  OS_all <- rbind(Need2,OS_all) |> 
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
OS_all1 <- OS_all |>
  left_join(Tumor_ICB, by = c("cohort" = "cohort", "Run" = "Run")) |>
  filter(!grepl("On", Response2)) |>
  filter(Response2 != "") |>
  group_by(Patient_ID) |>
  mutate(all_fusion_na = all(is.na(fusion))) |>
  filter(all_fusion_na | (!all_fusion_na & !is.na(fusion))) |>
  ungroup() |>
  filter(!is.na(OS_Time) &  !is.na(OS_Status) ) |>
  select(Patient_ID, OS_Time, OS_Status, fusion, cancer, cohort) |>
  distinct()
fusion_data_long <- OS_all1 %>%
  group_by(cancer, Patient_ID) %>%
  summarise(
    unique_fusions = ifelse(all(is.na(fusion)), 0, length(unique(na.omit(fusion)))),
    .groups = 'drop',
    OS_Time = first(OS_Time),
    OS_Status = first(OS_Status)
  ) %>%
  ungroup() 
cancer_type_count <- fusion_data_long$cancer |> unique() |> setdiff("BRCA")
OS_result <- data.table()
for(cancer_type1 in cancer_type_count){
  # "BLCA"  "GBM"   "HNSC"  "NSCLC" "SKCM"  "STAD" 
  
  fusion_data_long_cancer <- fusion_data_long   |>
    filter(cancer == cancer_type1) 
  
  if(cancer_type1 == "BLCA"){
    trans_result = 3
  }else if(cancer_type1 == "KIRC"){
    trans_result = 2
  }else if(cancer_type1 == "STAD"){
    trans_result = 4
  }else if(cancer_type1 == "HNSC"){
    trans_result = 2
  }else if(cancer_type1 == "SGC"){
    trans_result = 5
  }else if(cancer_type1 == "SKCM"){
    trans_result = 3
  }else if(cancer_type1 == "GBM"){
    trans_result = 8
  }else if(cancer_type1 == "NSCLC"){
    trans_result = 2
  }

   fusion_data_long1 <- fusion_data_long_cancer |> 
    mutate(event = case_when(
      unique_fusions < trans_result ~ "Low Fusion+",
      TRUE ~ "High Fusion+"
    ))
  fusion_data_long1$event <- factor( fusion_data_long1$event, level = c("Low Fusion+", "High Fusion+"))
  data.survdiff <- survdiff(Surv(OS_Time, OS_Status) ~ event, data =  fusion_data_long1)
  fit <- survfit(Surv(OS_Time, OS_Status) ~ event, data =  fusion_data_long1)
  
  hn = nrow(filter( fusion_data_long1, event == "High Fusion+"))
  ln = nrow( fusion_data_long1) - hn 
   fusion_data_long1$event <- factor( fusion_data_long1$event, level = c("Low Fusion+", "High Fusion+"))
  data.survdiff <- survdiff(Surv(OS_Time, OS_Status) ~ event, data =  fusion_data_long1)
  fit <- survfit(Surv(OS_Time, OS_Status) ~ event, data =  fusion_data_long1)
  p.val1 <- round(1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1), 3)
  p.val <- ifelse(p.val1 < 0.001, "\np< 0.001", paste0("\np:",p.val1))
  HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
  up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")
  
  caption1 <- paste0(
    paste0("High TFB :"), hn,
    paste0("\n", "Low TFB :"), ln,
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
                   title = paste0("ICB-",cancer_type1,"-OS"),
                   legend.labs = c("Low Fusion+", "High Fusion+"),
                   xlab = "Time(Days)"
  ) 
  p1$plot <- p1$plot + 
    annotate("text", Inf, Inf, label = caption1,  vjust = 1, hjust = 1, size = 4) +
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
  file_path <- paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup6/counts_KM/",
                      cancer_type1,"_",trans_result,"_OS.png")
  ggsave(filename = file_path,
         plot = p2, dpi = 600, width = 10, height = 6)
}


#icb pfs

PFS_all1 <- OS_all |>
  left_join(Tumor_ICB, by = c("cohort" = "cohort", "Run" = "Run")) |>
  filter(!grepl("On", Response2)) |>
  filter(Response2 != "") |>
  group_by(Patient_ID) |>
  mutate(all_fusion_na = all(is.na(fusion))) |>
  filter(all_fusion_na | (!all_fusion_na & !is.na(fusion))) |>
  ungroup() |>
  filter(!is.na(PFS_Time) &  !is.na(PFS_Status) ) |>
  select(Patient_ID, PFS_Time, PFS_Status, fusion, cancer, cohort) |>
  distinct()
fusion_data_long <- PFS_all1 %>%
  group_by(cancer, Patient_ID) %>%
  summarise(
    unique_fusions = ifelse(all(is.na(fusion)), 0, length(unique(na.omit(fusion)))),
    .groups = 'drop',
    PFS_Time = first(PFS_Time),
    PFS_Status = first(PFS_Status)
  ) %>%
  ungroup() 
cancer_type_count <- fusion_data_long$cancer |> unique()
PFS_result  <- data.table()
for(cancer_type1 in cancer_type_count){
  # "GBM"   "KIRC"  "NSCLC" "SGC"   "SKCM"  "STAD" 
  fusion_data_long_cancer <- fusion_data_long   |>
    filter(cancer == cancer_type1)

  if(cancer_type1 == "BLCA"){
    trans_result = 3
  }else if(cancer_type1 == "KIRC"){
    trans_result = 2
  }else if(cancer_type1 == "STAD"){
    trans_result = 4
  }else if(cancer_type1 == "HNSC"){
    trans_result = 2
  }else if(cancer_type1 == "SGC"){
    trans_result = 5
  }else if(cancer_type1 == "SKCM"){
    trans_result = 3
  }else if(cancer_type1 == "GBM"){
    trans_result = 8
  }else if(cancer_type1 == "NSCLC"){
    trans_result = 2
  }
    
    fusion_data_long1 <- fusion_data_long_cancer |> 
      mutate(event = case_when(
        unique_fusions < trans_result ~ "Low Fusion+",
        TRUE ~ "High Fusion+"
      ))
    fusion_data_long1$event <- factor( fusion_data_long1$event, level = c("Low Fusion+", "High Fusion+"))
    data.survdiff <- survdiff(Surv(PFS_Time, PFS_Status) ~ event, data =  fusion_data_long1)
    fit <- survfit(Surv(PFS_Time, PFS_Status) ~ event, data =  fusion_data_long1)
    
    hn = nrow(filter( fusion_data_long1, event == "High Fusion+"))
    ln = nrow( fusion_data_long1) - hn 
    fusion_data_long1$event <- factor( fusion_data_long1$event, level = c("Low Fusion+", "High Fusion+"))
    data.survdiff <- survdiff(Surv(PFS_Time, PFS_Status) ~ event, data =  fusion_data_long1)
    fit <- survfit(Surv(PFS_Time, PFS_Status) ~ event, data =  fusion_data_long1)
    p.val1 <- round(1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1), 3)
    p.val <- ifelse(p.val1 < 0.001, "\np< 0.001", paste0("\np:",p.val1))
    HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
    up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
    low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
    HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
    CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")
    
    caption1 <- paste0(
      paste0("High TFB :"), hn,
      paste0("\n", "Low TFB :"), ln,
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
                     legend.labs = c("Low Fusion+", "High Fusion+"),
                     title = paste0("ICB-",cancer_type1,"-PFS"),
                     xlab = "Time(Days)"
    ) 
    
    p1$plot <- p1$plot + 
      annotate("text", Inf, Inf, label = caption1,  vjust = 1, hjust = 1, size = 4) +
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
    file_path <- paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup6/counts_KM/",cancer_type1,"_",trans_result,"_PFS.png")
    ggsave(filename = file_path,
           plot = p2, dpi = 600, width = 10, height = 6)
  }
  
 
  


