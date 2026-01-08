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



#figure4A braf
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
gene_vec <- "BRAF"
fusion_select1 <- Tumor_ICB_OS |>
  filter(cohort %in% cohort_select)  |> 
  filter(gene1 %in% gene_vec | gene2 %in% gene_vec) |>
  as.data.table()
cohort_select2 <- fusion_select1$cohort |> unique()
OS_all <- data.table()
for(cohort_id in cohort_select2){
  Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv"))
  if(filter(cohort_SKCM ,id == cohort_id)$time == "Pre/On" ){
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
  
  columns_to_select <- c("Patient_ID", "OS_Time", "OS_Status", "PFS_Time", "PFS_Status", 
                         "Run", "Treatment", "Age", "Sex", "Response2")
  existing_columns <- columns_to_select[columns_to_select %in% colnames(Need1)]
  Need2 <- Need1 |>
    select(all_of(existing_columns)) #|>   filter(!is.na(OS_Time) & !is.na(OS_Status)) |>   filter(!is.na(PFS_Time) & !is.na(PFS_Status))
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

fusion_data_subset <- OS_all |>
  mutate(event = case_when(
    Run %in% fusion_select1$Run  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                           TRUE ~ event))  |>
  select(Patient_ID, event, OS_Time, OS_Status,PFS_Time, PFS_Status,
         Sex, Response2,cohort) |>
  group_by(cohort) |>
  filter(n_distinct(event) == 2) |>
  ungroup()  |>
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
data2$Sex2 <- c("F1", "M1",  "F2", "M2")
data2$Sex2 <-  factor(data2$Sex2, levels = c("F2", "M2","F1","M1"))

# sex and fusion
data3 <- fusion_data_subset  |>
  as.data.table()
data3 [, `:=`(
  event = factor(event, levels = c("Fusion-", "Fusion+")),
  Response2  = factor(Sex,  levels = c("M", "F"))
)]
#data3$event与data3$Sex的prop.test P = 1
p <- ggstatsplot::ggbarstats(
  data            = data3,
  x               = Sex,      
  y               = event,       
  paired          = FALSE,      
  results.subtitle = TRUE,  
  label           = "both",     
  title           = "Fusion vs Sex"
)
p
#plot
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
    name = "Fusion",
    values = c("Fusion+" = "#f28680", "Fusion-" = "#5861ac"),
    na.translate = FALSE
  ) +
  scale_fill_manual(
    name = "Sex",
    values = c(
      "Fusion+" = "#f28680",
      "Fusion-" = "#5861ac",
      "F1" = "#ff7f00",
      "M1" = "#ffe0c1",
      "F2" = "#ff7f00",
      "M2" = "#ffe0c1"),
    labels = c(
      "Fusion+" = "Fusion+",
      "Fusion-" = "Fusion-",
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

ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4A.png"),
       plot = p,         
       width = 7,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "png",           
       bg = "white")
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4A.pdf"),
       plot = p,         
       width = 7,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "pdf",           
       bg = "white")

#Figure4B

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
ggsave(filename = "C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4B.png", 
       plot = p2, dpi = 600, width = 6, height = 5)

#Figure4C
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
ggsave(filename = "C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4C.png", 
       plot = p2, dpi = 600, width = 6, height = 5)

#Figure4D
#bug 
fusion_subset1 <- fusion_data_subset |>
  select(Patient_ID, Response2, event) |> 
  filter(Response2 != "") |> 
  distinct() |>
  as.data.table()|>
  rename(Response = Response2) 
fusion_subset1$Response  <- gsub("Pre-", "", fusion_subset1$Response )
fusion_subset1[, `:=`(
  event = factor(event, levels = c("Fusion-", "Fusion+")),
  Response  = factor(Response,  levels = c("NR", "R"))
)]

a2 <-  fusion_subset1 |>
  group_by(event, Response) |>
  summarise(count_event = n()) |>
  mutate(total_count = sum(count_event), 
         percentage1 = (count_event / total_count) * 100)

p <- ggstatsplot::ggbarstats(
  data            = fusion_subset1,
  x               = Response,      
  y               = event,       
  paired          = FALSE,      
  results.subtitle = TRUE,  
  label           = "both",     
  title           = "Fusion status vs irAE"
)
subtitle <- p$labels$subtitle
p_value <- regmatches(subtitle, regexec('italic\\(p\\) == "([0-9.]+)"', subtitle))[[3]][3]

p <- ggstatsplot::ggbarstats(
  data            = fusion_subset1,
  x               = Response,      
  y               = event,      
  paired          = FALSE,      
  results.subtitle = FALSE,  
  label           = "both", 
  label.args       = list(color = NA, alpha = 0),  
  title           = "BRAF",
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
p
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4D.png"),
       plot = p, dpi = 600, width = 5, height = 6)
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4D.pdf"),
       plot = p, dpi = 600, width = 5, height = 6)

#Figure4E
gene_vec <- "LSAMP"
fusion_select1 <- Tumor_ICB_OS |>
  filter(cohort %in% cohort_select)  |> 
  filter(gene1 %in% gene_vec | gene2 %in% gene_vec) |>
  as.data.table()
cohort_select2 <- fusion_select1$cohort |> unique()
OS_all <- data.table()
for(cohort_id in cohort_select2){
  Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv"))
  if(filter(cohort_SKCM ,id == cohort_id)$time == "Pre/On" ){
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
  columns_to_select <- c("Patient_ID", "OS_Time", "OS_Status",  
                         "PFS_Time", "PFS_Status", "LDH",
                         "Run", "Treatment", "primary","M_Stage",
                         "Age", "Sex", "Response2")
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
  OS_all <- rbind(Need2,OS_all) |> filter(!Treatment %in% c("Docetaxel","Observation","Sunitinib")) 
}

fusion_data_subset <- OS_all |>
  mutate(event = case_when(
    Run %in% fusion_select1$Run  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) |>
  group_by(Patient_ID) |>
  mutate(event = case_when(any(event == "Fusion+") ~ "Fusion+",
                           TRUE ~ event))  |>
  group_by(cohort) |>
  filter(n_distinct(event) == 2) |>
  ungroup()  |>
  select(Patient_ID, event, OS_Time, OS_Status, Sex, Response2,cohort) |>
  distinct() 
#Figure4E

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
  "Fusion+" = "#f28680",
  "Fusion-" = "#5861ac",
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
      "Fusion+" = "#f28680",
      "Fusion-" = "#5861ac",
      "F1" = "#ff7f00",
      "M1" = "#ffe0c1",
      "F2" = "#ff7f00",
      "M2" = "#ffe0c1"),
    labels = c(
      "Fusion+" = "Fusion+",
      "Fusion-" = "Fusion-",
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
 
 
ggsave("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4E.png", 
       plot = p,         
       width = 7,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "png",           
       bg = "white") 
ggsave("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4E.pdf", 
       plot = p,         
       width = 7,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "pdf",           
       bg = "white") 


#Figure2F
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
ggsave(filename = "C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4F.png", plot = p2, dpi = 600, width = 6, height = 5)


#Figure2G
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
ggsave(filename = "C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4G.png", plot = p2, dpi = 600, width = 6, height = 5)


#Figure4H
fusion_subset1 <- fusion_data_subset |>
                 select(Patient_ID, Response2, event) |> 
                 filter(Response2 != "") |> 
                 distinct() |>
                 as.data.table() |>
                 rename(Response = Response2) 
fusion_subset1$Response  <- gsub("Pre-", "", fusion_subset1$Response )
fusion_subset1[, `:=`(
  event = factor(event, levels = c("Fusion-", "Fusion+")),
  Response  = factor(Response,  levels = c("NR", "R"))
)]

a2 <-  fusion_subset1 |>
  group_by(event, Response) |>
  summarise(count_event = n()) |>
  mutate(total_count = sum(count_event), 
         percentage1 = (count_event / total_count) * 100)

p <- ggstatsplot::ggbarstats(
  data            = fusion_subset1,
  x               = Response2,      
  y               = event,       
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
  x               = Response,      
  y               = event,      
  paired          = FALSE,      
  results.subtitle = FALSE,  
  label           = "both", 
  label.args       = list(color = NA, alpha = 0),  
  title           = "LSAMP",
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
p
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4H.png"),
       plot = p, dpi = 600, width = 5, height = 6)
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure4/Figure4H.pdf"),
       plot = p, dpi = 600, width = 5, height = 6)







