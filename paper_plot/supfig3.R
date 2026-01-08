library(data.table)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tidyr)
library(stringr)
library(patchwork)

file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta , sep = "\t") |>
  separate_rows(gene1, sep = ",") |>
  separate_rows(gene2, sep = ",")
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv")
cohort_SKCM <- cohort1 |> filter(OS == TRUE & PFS == TRUE & cancer_type == "SKCM") 
cohort_select <- (cohort_SKCM$id) |> unique()
Tumor_ICB_OS <- fusion_all |>
  as.data.table() |>
  filter(sample_type != "Normal") |>
  filter(!grepl("TCGA|CPTAC|TARGET", cohort)) |>
  filter(cohort %in% cohort_select) |> 
  select(Sample_ID, gene1, gene2,cohort)|>
  rename(Run = Sample_ID) |>
  distinct()
gene_vec <- "LSAMP"
fusion_select1 <- Tumor_ICB_OS |>
  filter(cohort %in% cohort_select)  |> 
  filter(gene1 %in% gene_vec | gene2 %in% gene_vec) |>
  as.data.table()
cohort_select2 <- fusion_select1$cohort |> unique() |> setdiff("PRJEB23709")
OS_all <- data.table()
for(cohort_id in cohort_select2){
  Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv"))
  Need1 <- Need |>
           filter(grepl("Pre", Response2)) 
  if(cohort_id == "PHS000452_PD1"){
    Need1$primary <- Need1$Primary_Type
    Need1$LDH <- Need1$LDH_Elevated
  }
  if(cohort_id == "PHS000452_CTLA4"){
    Need1$M_Stage <- Need1$M
  }
  columns_to_select <- c("Patient_ID", "OS_Time", "OS_Status",  
                         "PFS_Time", "PFS_Status", "LDH",
                         "Run", "Treatment", "primary","M_Stage",
                         "Age", "Sex", "Response2")
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
  OS_all <- rbind(Need2,OS_all) |> filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
}

#supA  #supF  
for (cohort_name in cohort_select2) {
  fusion_select <- fusion_select1 |>
                   filter(cohort == cohort_name)
  
  fusion_data_subset <- OS_all |>
    filter(cohort == cohort_name) |>
    mutate(event = case_when(
      Run %in% fusion_select$Run  ~  "Fusion+",
      TRUE ~ "Fusion-"
    )) |>
    group_by(Patient_ID) |>
    mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                             TRUE ~ event))  |>
    select(Patient_ID, event, OS_Time, OS_Status, Sex, Response2) |>
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
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup3/supfig3_",cohort_name,"_OS.png"),
         plot = p2, dpi = 600, width = 6, height = 5)
  
  
  
}
#supB  #supG  
for (cohort_name in cohort_select2) {
  fusion_select <- fusion_select1 |>
    filter(cohort == cohort_name)
  
  fusion_data_subset <- OS_all |>
    filter(cohort == cohort_name) |>
    mutate(event = case_when(
      Run %in% fusion_select$Run  ~  "Fusion+",
      TRUE ~ "Fusion-"
    )) |>
    group_by(Patient_ID) |>
    mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                             TRUE ~ event))  |>
    select(Patient_ID, event, PFS_Time, PFS_Status, Sex, Response2) |>
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
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup3/supfig3_",cohort_name,"_PFS.png"),
         plot = p2, dpi = 600, width = 6, height = 5)
  
  
  
}
#supD  #supH  
for (cohort_name in cohort_select2) {
  fusion_select <- fusion_select1 |>
    filter(cohort == cohort_name)
  
  fusion_data_subset <- OS_all |>
    filter(cohort == cohort_name) |>
    mutate(event = case_when(
      Run %in% fusion_select$Run  ~  "Fusion+",
      TRUE ~ "Fusion-"
    )) |>
    group_by(Patient_ID) |>
    mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                             TRUE ~ event))  |>
    select(Patient_ID, event,  Sex, Response2) |>
    distinct()
  fusion_subset <- fusion_data_subset |>
    select(Patient_ID, Response2, event) |> 
    filter(Response2 != "") |> 
    distinct() |>
    as.data.table() |>
    rename(Response = Response2,
           Fusion = event) 
  fusion_subset$Response  <- gsub("Pre-", "", fusion_subset$Response )
  fusion_subset[, `:=`(
    Fusion = factor(Fusion, levels = c("Fusion-", "Fusion+")),
    Response  = factor(Response,  levels = c("NR", "R"))
  )]
  p <- ggstatsplot::ggbarstats(
    data            = fusion_subset,
    x               = Response,      
    y               = Fusion,       
    paired          = FALSE,      
    results.subtitle = TRUE,  
    label           = "both",     
    title           = "Fusion status vs irAE",
    type = "nonparametric" 
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
    label.args       = list(color = NA, alpha = 0),  
    title           = paste0(cohort_name,"-LSAMP"),
    ggtheme = theme_classic()
  ) + scale_fill_manual(values = c(
    "R" = "#9dd3af",
    "NR" = "#f2f2f2"
  )) + 
    labs(
      subtitle = paste0("p: ",  p_value)
    ) +  theme(
      plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图3/",cohort_name,"_response.png"),
         plot = p, dpi = 600, width = 5, height = 6)
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图3/",cohort_name,"_response.pdf"),
         plot = p, dpi = 600, width = 5, height = 6)
  

  
  
  
}

#sup primary
for (cohort_name in cohort_select2) {
  fusion_select <- fusion_select1 |>
    filter(cohort == cohort_name)
  
  fusion_data_subset <- OS_all |>
    filter(cohort == cohort_name) |>
    mutate(event = case_when(
      Run %in% fusion_select$Run  ~  "Fusion+",
      TRUE ~ "Fusion-"
    )) |>
    group_by(Patient_ID) |>
    mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                             TRUE ~ event))  |>
    select(Patient_ID, event,  primary) |>
    distinct()
  fusion_subset <- fusion_data_subset |>
    select(Patient_ID, primary, event) |> 
    filter(primary != "") |> 
    distinct() |>
    as.data.table()
  row.names(fusion_subset) <-  fusion_subset$Patient_ID
  fusion_subset1 <- fusion_subset |> 
                    select(-Patient_ID)  
  if(cohort_name  == "PHS000452_PD1"){
    fusion_subset1[, `:=`(
      event = factor(event, levels = c("Fusion-", "Fusion+")),
      primary  = factor(primary,  levels = c("occult", "skin", "mucosal", "acral" ))
    )]
  }else{
    fusion_subset1[, `:=`(
      event = factor(event, levels = c("Fusion-", "Fusion+")),
      primary  = factor(primary,  levels = c("occult",  "skin",    "mucosal")) 
    )]
  }
  p <- ggstatsplot::ggbarstats(
    data            = fusion_subset1,
    x               = event,      
    y               =  primary,       
    paired          = FALSE,      
    results.subtitle = TRUE,  
    label           = "both",     
    title           = "Fusion status vs irAE",
    type = "nonparametric" 
  )
  subtitle <- p$labels$subtitle
  p_value <- regmatches(subtitle, regexec('italic\\(p\\) == "([0-9.]+)"', subtitle))[[3]][3]
  
  p <- ggstatsplot::ggbarstats(
    data            = fusion_subset1,
    x               = event,      
    y               =  primary,         
    paired          = FALSE,      
    results.subtitle = FALSE,  
    label           = "both", 
    label.args       = list(color = NA, alpha = 0),  
    title           = paste0(cohort_name,"-LSAMP"),
    ggtheme = theme_classic()
  ) + scale_fill_manual(values = c(
    "occult" = "#f9e7a7", 
    "skin" = "#ef776b", 
    "mucosal" = "#43a3ef",
    "acral" = "#93c8c0",
    "Fusion-" =  "#3c9bc9", 
    "Fusion+" = "#fc757b"
  )) + 
    labs(
      subtitle = paste0("p: ",  p_value)
    ) +  theme(
      plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup3/",
                           cohort_name,"_primary.png"),
         plot = p, dpi = 600, width = 5, height = 6)
  
}
#sup M_Stage
for (cohort_name in cohort_select2) {
  fusion_select <- fusion_select1 |>
    filter(cohort == cohort_name)
  
  fusion_data_subset <- OS_all |>
    filter(cohort == cohort_name) |>
    mutate(event = case_when(
      Run %in% fusion_select$Run  ~  "Fusion+",
      TRUE ~ "Fusion-"
    )) |>
    group_by(Patient_ID) |>
    mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                             TRUE ~ event))  |>
    select(Patient_ID, event,  M_Stage) |>
    distinct()
  fusion_subset <- fusion_data_subset |>
    select(Patient_ID, M_Stage, event) |> 
    filter(M_Stage != "") |> 
    distinct() |>
    as.data.table()
  row.names(fusion_subset) <-  fusion_subset$Patient_ID
  fusion_subset <- fusion_subset |> select(-Patient_ID)
  if(cohort_name  == "PHS000452_PD1"){
  fusion_subset[, `:=`(
      event = factor(event, levels = c("Fusion-", "Fusion+")),
      M_Stage  = factor(M_Stage,  levels = c("M1a", "M1b", "M1c" , "IIIC"))
    )]
  }else{
    fusion_subset[, `:=`(
      event = factor(event, levels = c("Fusion-", "Fusion+")),
      M_Stage  = factor(M_Stage,  levels = c("M1a", "M1b", "M1c" ))
    )]  
  }
  p <- ggstatsplot::ggbarstats(
    data            = fusion_subset,
    x               = event,      
    y               =  M_Stage,       
    paired          = FALSE,      
    results.subtitle = TRUE,  
    label           = "both",     
    title           = "Fusion status vs irAE",
    type = "nonparametric" 
  )
  subtitle <- p$labels$subtitle
  p_value <- regmatches(subtitle, regexec('italic\\(p\\) == "([0-9.]+)"', subtitle))[[3]][3]
  
  p <- ggstatsplot::ggbarstats(
    data            = fusion_subset,
    x               = event,      
    y               =  M_Stage,         
    paired          = FALSE,      
    results.subtitle = FALSE,  
    label           = "both", 
    label.args       = list(color = NA, alpha = 0),  
    title           = paste0(cohort_name,"-LSAMP"),
    ggtheme = theme_classic()
  ) + scale_fill_manual(values = c(
    "Fusion-" =  "#3c9bc9", 
    "Fusion+" = "#fc757b"
  )) + 
    labs(
      subtitle = paste0("p: ",  p_value)
    ) +  theme(
      plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup3/",
                           cohort_name,"_M_Stage.png"),
         plot = p, dpi = 600, width = 5, height = 6)
  
}



