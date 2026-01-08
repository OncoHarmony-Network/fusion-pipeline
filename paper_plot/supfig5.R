
# Fusion- and Fusion+ 
# Fusion+ and low Fusion
#sup4 fusion counts and survival
library(ggplot2)
library(survival)
library(survminer)
library(data.table)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggthemes)
library(stringr)
library(patchwork)

#tfb AND response 
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
  left_join(Tumor_ICB, by = c("cohort" = "cohort", "Run" = "Run")) |>
  group_by(Patient_ID) |>
  mutate(all_fusion_na = all(is.na(fusion))) |>
  filter(all_fusion_na | (!all_fusion_na & !is.na(fusion))) |>
  ungroup() |>
  select(Patient_ID, Response2, fusion, cancer, cohort) |>
  distinct() |>
  rename(Response = Response2)
Response_all1$Response <- gsub("Pre-", "", Response_all1$Response )
cancer_all <- Response_all1$cancer |> unique() |> setdiff("BRCA")
for (cancer1 in cancer_all){
  fusion_data_long_cancer <- Response_all1 %>%
    filter(cancer == cancer1) %>%
    as.data.table() %>%
    group_by(cancer, Patient_ID) %>%
    summarise(
      unique_fusions = ifelse(all(is.na(fusion)), 0, length(unique(na.omit(fusion)))),
      .groups = 'drop',
     Response = first(Response)
    ) %>%
    ungroup() 
  
  if(cancer1 == "BLCA"){
    trans_result = 3
  }else if(cancer1 == "KIRC"){
    trans_result = 2
  }else if(cancer1 == "STAD"){
    trans_result = 4
  }else if(cancer1 == "HNSC"){
    trans_result = 2
  }else if(cancer1 == "SGC"){
    trans_result = 5
  }else if(cancer1 == "SKCM"){
    trans_result = 3
  }else if(cancer1 == "GBM"){
    trans_result = 8
  }else if(cancer1 == "NSCLC"){
    trans_result = 2
  }
    
    fusion_subset <- fusion_data_long_cancer |> 
      mutate(TFB = case_when(
        unique_fusions < trans_result ~ "Low TFB",
        TRUE ~ "High TFB"
      ))|> 
      as.data.table()
    
    fusion_subset[, `:=`(
      fusion = factor(TFB, levels = c("Low TFB", "High TFB")),
      Response  = factor(Response,  levels = c("NR", "R"))
    )]
    
    p <- ggstatsplot::ggbarstats(
      data            = fusion_subset,
      x               = Response,      
      y               =  fusion,     
      paired          = FALSE,      
      results.subtitle = TRUE,  
      label           = "both",     
      title           = ""
    )
    subtitle <- p$labels$subtitle
    p_value <- regmatches(subtitle, regexec('italic\\(p\\) == "([0-9.]+)"', subtitle))[[3]][3]
    p <- ggstatsplot::ggbarstats(
      data            = fusion_subset,
      x               = Response,     
      y               = fusion,      
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
    
    ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup5/counts_response/",cancer1,"_",trans_result,"_response.png"), 
           plot = p, dpi = 600, width = 5, height = 7)
  }




#Figurea,c,e,g,h,i,k
library(ggplot2)
library(data.table)
library(dplyr)
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
  left_join(Tumor_ICB, by = c("cohort" = "cohort", "Run" = "Run")) |>
  group_by(Patient_ID) |>
  mutate(all_fusion_na = all(is.na(fusion))) |>
  filter(all_fusion_na | (!all_fusion_na & !is.na(fusion))) |>
  ungroup() |>
  select(Patient_ID, Response2, fusion, cancer, cohort) |>
  distinct()


cancer_all <- Response_all1$cancer |> unique() |> setdiff("BRCA")
for (cancer1 in cancer_all){
  fusion_data_long_cancer <- Response_all1 %>%
    filter(cancer == cancer1) %>%
    as.data.table() %>%
    group_by(cancer, Patient_ID) %>%
    summarise(
      unique_fusions = ifelse(all(is.na(fusion)), 0, length(unique(na.omit(fusion)))),
      .groups = 'drop',
      Response2 = first(Response2)
    ) %>%
    ungroup() |>
    mutate(Response = Response2)
  fusion_data_long <- fusion_data_long_cancer
  fusion_data <-  fusion_data_long   |>
    group_by(unique_fusions) |> 
    summarise(count_event = n()) |>
    ungroup() |> 
    mutate(total_count = sum(count_event), 
           percentage = (count_event / total_count) * 100)
  
  result <- fusion_data|>
    mutate(cum_percentage = cumsum(percentage)) |>
    filter(cum_percentage > 50) |>
    slice(1) |>
    select(unique_fusions, cum_percentage)
  
  # 输出结果
  print(paste0(cancer1,result))
  
  fusion_data_long_cancer$Response <- gsub("Pre-", "", fusion_data_long_cancer$Response )
  
  bar_data <- fusion_data_long_cancer %>%
    group_by(Response, unique_fusions) %>%
    summarise(count = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = Response, values_from = count, values_fill = 0) |>
    arrange(unique_fusions) |>
    rename(TFB = unique_fusions)
bar_data$TFB <- as.character(bar_data$TFB)  
bar_data$TFB <- factor(bar_data$TFB, levels = bar_data$TFB)
col <- c("R" = "#9dd3af", "NR" = "#f2f2f2") 

p2<-ggplot(bar_data,aes(y=TFB))+
  geom_bar(aes(x=R),stat="identity",fill="#9dd3af")+
  geom_bar(aes(x=-NR),stat="identity",fill="#f2f2f2")+
  theme_bw()+
  theme(
    panel.border=element_blank(),
    axis.line.x=element_line(),
    axis.line.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.y = element_text(size = 15), 
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 25),
    axis.text.x = element_text(size = 15),
    panel.grid=element_blank(),
    panel.background=element_blank(),
    plot.background=element_blank()
  )+
  labs(x=NULL,y="TFB",title = cancer1)
print(p2)

ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/sup5/",
                         cancer1,"_TFB&Response.pdf"), 
       plot = p2, dpi = 600, width = 6, height = 8)
}








