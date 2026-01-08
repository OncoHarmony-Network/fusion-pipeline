library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tidyr)
library(stringr)
library(patchwork)
#sup4 gdc os
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
metafusion <- fread(file_meta , sep = "\t")
Credibility_gene <- metafusion |>
  separate_rows(gene1, sep = ",") |>
  separate_rows(gene2, sep = ",")  |> 
  dplyr::filter(sample_type == "Tumor")  |>
  filter(grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(!cohort %in% c("CPTAC-NORMAL", "TCGA-NORMAL","TARGET-CELL", "TARGET-NORMAL")) |>
  filter(gene1 == "LSAMP" | gene2 == "LSAMP")

gene1 <- Credibility_gene |> select(c(cohort, gene1, Sample_ID)) |> rename(gene = gene1)
gene2 <- Credibility_gene |> select(c(cohort, gene2, Sample_ID)) |> rename(gene = gene2)
fusion_select2 <- rbind(gene1,gene2) %>%
  distinct()%>%
  group_by(cohort, gene) %>%
  filter(n_distinct(Sample_ID) >= 2) %>%
  ungroup()

cohort_sum_OS1 <- fusion_select2$cohort |> unique()

for(cohort_id in cohort_sum_OS1){
  cohort_id2 <- gsub("_", "-", cohort_id)
  Need_all <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id2,"_info.tsv"))
  fusion_select <- Credibility_gene |> 
    dplyr::filter(cohort == cohort_id2)|> 
    as.data.table()
  
  fusion_data_subset <-   Need_all |> 
    mutate(event = case_when(
      Run %in% fusion_select$Sample_ID  ~  "Fusion+",
      TRUE ~ "Fusion-"
    )) |>
    group_by(Patient_ID) |>
    mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                             TRUE ~ event))  |>
    filter(!is.na(OS_Time) & !is.na(OS_Status)) |>
    select(Patient_ID, event, OS_Time, OS_Status) |>
    distinct() 
  hn <- nrow(filter(fusion_data_subset, event == "Fusion+")) 
  n_all <- length(fusion_data_subset$Patient_ID)
  if(hn < 3 ){
    next()
  }else{
    ln = nrow(fusion_data_subset) - hn 
    fusion_data_subset$event <- factor(fusion_data_subset$event, level = c("Fusion-", "Fusion+"))
    data.survdiff <- survdiff(Surv(OS_Time, OS_Status) ~ event, data = fusion_data_subset)
    fit <- survfit(Surv
                   (OS_Time, OS_Status) ~ event, data = fusion_data_subset)
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
                     legend.labs = c( "Fusion-", "Fusion+"),
                     title = paste0(cohort_id,"-OS"),
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
    ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup4/",cohort_id,".OS.png"), 
           plot = p2, dpi = 600, width = 6, height = 5)
    
  }
  
} 


#sup4 GDC PFS
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
metafusion <- fread(file_meta , sep = "\t")
Credibility_gene <- metafusion |>
  separate_rows(gene1, sep = ",") |>
  separate_rows(gene2, sep = ",")  |> 
  dplyr::filter(sample_type == "Tumor")  |>
  filter(grepl("TCGA|CPTAC", cohort)) |>
  filter(!cohort %in% c("CPTAC-NORMAL", "TCGA-NORMAL", "TCGA-LAML")) |>
  filter(gene1 == "LSAMP" | gene2 == "LSAMP")

gene1 <- Credibility_gene |> select(c(cohort, gene1, Sample_ID)) |> rename(gene = gene1)
gene2 <- Credibility_gene |> select(c(cohort, gene2, Sample_ID)) |> rename(gene = gene2)
fusion_select2 <- rbind(gene1,gene2) %>%
  distinct()%>%
  group_by(cohort, gene) %>%
  filter(n_distinct(Sample_ID) >= 2) %>%
  ungroup()

cohort_sum_PFS1 <- fusion_select2$cohort |> unique()

for(cohort_id in cohort_sum_PFS1){
  cohort_id2 <- gsub("_", "-", cohort_id)
  Need_all <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id2,"_info.tsv"))
  fusion_select <- Credibility_gene |> 
    dplyr::filter(cohort == cohort_id2)|> 
    as.data.table()
  
  fusion_data_subset <-   Need_all |> 
    mutate(event = case_when(
      Run %in% fusion_select$Sample_ID  ~  "Fusion+",
      TRUE ~ "Fusion-"
    )) |>
    group_by(Patient_ID) |>
    mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                             TRUE ~ event))  |>
    filter(!is.na(PFS_Time) & !is.na(PFS_Status)) |>
    select(Patient_ID, event, PFS_Time, PFS_Status) |>
    distinct() 
  hn <- nrow(filter(fusion_data_subset, event == "Fusion+")) 
  n_all <- length(fusion_data_subset$Patient_ID)
  if(hn < 3 ){
    next()
  }else{
    ln = nrow(fusion_data_subset) - hn 
    fusion_data_subset$event <- factor(fusion_data_subset$event, level = c("Fusion-", "Fusion+"))
    data.survdiff <- survdiff(Surv(PFS_Time, PFS_Status) ~ event, data = fusion_data_subset)
    fit <- survfit(Surv
                   (PFS_Time, PFS_Status) ~ event, data = fusion_data_subset)
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
                     legend.labs = c( "Fusion-", "Fusion+"),
                     title = paste0(cohort_id,"-PFS"),
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
    ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup4/",cohort_id,".PFS.png"), 
           plot = p2, dpi = 600, width = 6, height = 5)
    
  }
  
} 