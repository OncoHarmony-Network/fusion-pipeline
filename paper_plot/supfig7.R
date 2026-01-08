#sup7A
library(ggstatsplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(swimplot)
library(patchwork)
library(matrixStats)
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta) |>
  tidyr::separate_rows(gene1, gene1_type, sep = ",") |>
  tidyr::separate_rows(gene2, gene2_type, sep = ",") 
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv") 
cohort_id <- "PRJNA795330"
gene_name <- "AL645608.2"
fusion_select <- fusion_all |>
  filter(cohort == cohort_id) |> 
  filter(gene1 == gene_name | gene2 == gene_name)
gene_name1 <- "JAZF1"
fusion_select1 <- fusion_all |>
  filter(cohort == cohort_id) |> 
  filter(gene1 == gene_name1| gene2 == gene_name1)

fusion_select_all <- fusion_all |>
  filter(cohort == cohort_id)

Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv")) |>
  filter(Sampling == "Pre") |> 
  mutate(
    irAE = if_else(rowSums(across(starts_with("irAE_days_"), ~ !is.na(.))) > 0, 1, 0),
    Fusions = case_when(
      Run %in% fusion_select$Sample_ID  ~  "AL645608.2 Fusion+",
      Run %in% fusion_select1$Sample_ID  ~  "JAZF1 Fusion+",
      TRUE ~ "Fusion-")
  ) |>
  filter(irAE == 1) |>
  rowwise() %>%
  mutate(
    irAE_time = max(c_across(starts_with("irAE_days_")), na.rm = TRUE)
  ) %>%
  ungroup() |>
  arrange(desc(irAE_time)) |> as.data.frame() 
Need_metadata1 <- Need |> 
  select(Patient_ID, Fusions)  |>
  mutate(
    event = case_when(
      Fusions == "AL645608.2 Fusion+" ~ "Fusion+",
      TRUE ~ "Fusion-"),
    Fusions = "AL645608.2 Fusion"
  )
Need_metadata2 <- Need |> 
  select(Patient_ID, Fusions)  |> 
  mutate(
    event = case_when(
      Fusions == "JAZF1 Fusion+" ~ "Fusion+",
      TRUE ~ "Fusion-"
    ),
    Fusions = "JAZF1 Fusion")


Need_metadata3 <-  fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv")) %>%
  filter(Sampling == "Pre") %>% 
  mutate(
    irAE = if_else(rowSums(across(starts_with("irAE_days_"), ~ !is.na(.))) > 0, 1, 0),
    event = case_when(
      Run %in% fusion_select_all$Sample_ID    ~ "Fusion+",
      TRUE ~ "Fusion-"),
    Fusions = "Fusion"
  )%>%
  filter(irAE == 1) %>%
  rowwise() %>%
  mutate(
    irAE_time = max(c_across(starts_with("irAE_days_")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(desc(irAE_time)) %>% 
  as.data.frame() %>%
  select(Patient_ID, Fusions, event)  

Need_metadata <-rbind(Need_metadata1,Need_metadata2,Need_metadata3)  %>%
  as.data.frame()

response_timeline <- Need %>%
  select(-irAE_time) %>%
  pivot_longer(
    cols = starts_with("irAE_days_"),
    names_to = "Events_types",
    values_to = "time"
  ) %>%
  mutate(
    Events_types = sub("irAE_days_", "", Events_types),
    time = as.numeric(time)
  ) %>%
  filter(!is.na(time))  |>
  select(Patient_ID, Events_types, time) |>as.data.frame()

p <- swimmer_plot(df = Need,
                  id = "Patient_ID",
                  end = "irAE_time",
                  name_fill = "Treatment",
                  id_order = "irAE_time",
                  increasing = FALSE,
                  col = NA,
                  alpha = 1,
                  width = 0.8) +
  swimmer_points(df_points = response_timeline,
                 id = "Patient_ID",
                 time = "time",
                 name_shape = "Events_types",
                 name_col = "Events_types",
                 size = 2) +
  coord_flip(clip = 'off', ylim = c(0, 380)) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 380, 10)) +
  scale_fill_manual(name = "Treatment",
                    values = c("ICI" = "#ffe0c1")) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  scale_color_manual(name="irAE events",
                     values=c(Skin = "#B092B6", `Flu-like` = "#AED2E2", Gastrointestinal = "#93c8c0",  Musculoskeletal = "black", 
                              Neurologic = "#C74D26", Pulmonary = "#e38d26",  Thyroid = "#a4c97c"),
                     breaks=c("Skin", "Flu-like", "Gastrointestinal", "Musculoskeletal",  "Neurologic",      
                              "Pulmonary",        "Thyroid"))+
  scale_shape_manual(name = "irAE events",
                     values = c(Skin = 17, `Flu-like` = 16, Gastrointestinal = 18,  Musculoskeletal = 3, Pulmonary = 15,
                                Neurologic = 5,  Thyroid = 14),
                     breaks = c("Skin", "Flu-like", "Gastrointestinal", "Musculoskeletal",  "Neurologic",      
                                "Pulmonary",        "Thyroid")) +
  guides(
    color = guide_legend(override.aes = list(size = 14)),  
    shape = guide_legend(override.aes = list(size = 14))
  ) + 
  labs(y = "Days since Treatment", x = "Patients enrolled") +
  theme_classic(base_size = 14) +
  theme(axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.key.spacing.y = unit(0.2, "cm"),
        legend.key.size = unit(1, 'cm'), 
        legend.title = element_text(size=16),
        legend.text = element_text(size=14)
  )


print(p)


a <- ggplot() +
  geom_tile(data = Need_metadata |> filter(Fusions == "AL645608.2 Fusion"),
            aes(x = Fusions, y = Patient_ID, fill = event),
            stat = "identity", position = "identity",
            width = 0.7, height = 0.7,
            na.rm = FALSE, color = 'black', linewidth = 0.6) +
  scale_fill_manual(name = "AL645608.2 Fusion",
                    values = c("Fusion-" = "#6AADD6", "Fusion+" = "#203468"),
                    na.value = "white") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = Need_metadata |> filter(Fusions == "JAZF1 Fusion"),
            aes(x = Fusions, y = Patient_ID, fill = event),
            stat = "identity", position = "identity",
            width = 0.7, height = 0.7,
            na.rm = FALSE, color = 'black', linewidth = 0.6) +
  scale_fill_manual(name = "JAZF1 Fusion",
                    values = c("Fusion-" = "#F69173", "Fusion+" ="#98361F"),
                    na.value = "white") +
  guides(fill = guide_legend(nrow = 1)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = Need_metadata |> filter(Fusions == "Fusion"),
            aes(x = Fusions, y = Patient_ID, fill = event),
            stat = "identity", position = "identity",
            width = 0.7, height = 0.7,
            na.rm = FALSE, color = 'black', linewidth = 0.6) +
  scale_fill_manual(name = "Fusion",
                    values = c("Fusion-" = "#B4D88A", "Fusion+" ="#00675E"),
                   na.value = "white") +
#  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "Patients enrolled") +
  theme_classic(base_size = 14) +
  theme(
    plot.background = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks = element_blank(),
    legend.key.spacing.y = unit(0.2, "cm"),
    legend.key.size = unit(1, 'cm'), 
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    legend.position = "right"
  )


print(a)
swim <- a + p + plot_layout(width = c(1,19), guides ="collect")
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图7/sup7A.pdf"), 
       plot = swim, dpi = 600, width = 20, height = 30)


#sup7B
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta) |>
  tidyr::separate_rows(gene1, gene1_type, sep = ",") |>
  tidyr::separate_rows(gene2, gene2_type, sep = ",") 
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv") 
cohort_id <- "PRJNA795330"
gene_name <- "AL645608.2"
Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv")) |>
  filter(Sampling == "Pre") |>
  mutate(irAE = case_when(
    !is.na(`irAE_days_Skin`) ~ "Skin yes",
    TRUE ~ "Skin no" ))
sample_all <- Need$Run
fusion_select <- fusion_all |>
  filter(cohort == cohort_id & Sample_ID %in% sample_all) 

fusion_select1 <- fusion_select |> filter(gene1 == gene_name| gene2 == gene_name)
fusion_subset <- Need |>
  mutate(event = case_when(
    Run %in% fusion_select1$Sample_ID  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) 

fusion_subset[, `:=`(
  event = factor(event, levels = c("Fusion+", "Fusion-")),
  irAE  = factor(irAE,  levels = c("Skin yes", "Skin no"))
)]

p <- ggstatsplot::ggbarstats(
  data            = fusion_subset,
  x               =  irAE,      
  y               =  event,     
  paired          = FALSE,      
  results.subtitle = TRUE,  
  label           = "both",     
  title           = "Fusion status vs irAE"
)
subtitle <- p$labels$subtitle
p_value <- regmatches(subtitle, regexec('italic\\(p\\) == "([0-9.]+)"', subtitle))[[3]][3]
p <- ggstatsplot::ggbarstats(
  data            = fusion_subset,
  x               =  irAE,     
  y               =  event,      
  paired          = FALSE,      
  results.subtitle = FALSE,  
  label           = "both",     
  title           = gene_name
)
p <- p + labs(
  subtitle = paste0("p: ",  p_value)
)
p
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图7/sup7B.pdf"), 
       plot = p, dpi = 600, width = 4, height = 6)

#sup7C
cohort_id <- "PRJNA795330"
gene_name <- "JAZF1"
Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv")) |>
  filter(Sampling == "Pre") |>
  mutate(irAE = case_when(
    !is.na(`irAE_days_Flu-like`) ~ "Flu-like yes",
    TRUE ~ "Flu-like no" ))
sample_all <- Need$Run
fusion_select <- fusion_all |>
  filter(cohort == cohort_id & Sample_ID %in% sample_all) 

fusion_select1 <- fusion_select |> filter(gene1 == gene_name| gene2 == gene_name)
fusion_subset <- Need |>
  mutate(event = case_when(
    Run %in% fusion_select1$Sample_ID  ~  "Fusion+",
    TRUE ~ "Fusion-"
  )) 

fusion_subset[, `:=`(
  event = factor(event, levels = c("Fusion+", "Fusion-")),
  irAE  = factor(irAE,  levels = c("Flu-like yes", "Flu-like no"))
)]

p <- ggstatsplot::ggbarstats(
  data            = fusion_subset,
  x               =  irAE,      
  y               =  event,     
  paired          = FALSE,      
  results.subtitle = TRUE,  
  label           = "both",     
  title           = "Fusion status vs irAE"
)
p_value <- "1.71e-03"
p <- ggstatsplot::ggbarstats(
  data            = fusion_subset,
  x               =  irAE,     
  y               =  event,      
  paired          = FALSE,      
  results.subtitle = FALSE,  
  label           = "both",     
  title           = gene_name
)
p <- p + labs(
  subtitle = paste0("p: ",  p_value)
)
p
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图7/sup7C.pdf"), 
       plot = p, dpi = 600, width = 4, height = 6)

#sup7A
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta) |>
  tidyr::separate_rows(gene1, gene1_type, sep = ",") |>
  tidyr::separate_rows(gene2, gene2_type, sep = ",") |>
  mutate(fusion = paste0(gene1, "::", gene2)) |> 
  filter(cohort == "PRJNA795330") |>
  select(Sample_ID, fusion) |>
  distinct() |> 
  rename(Run = Sample_ID)
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv")
cohort_id <- "PRJNA795330"
Need <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv")) |>
  mutate(Sample_type = case_when(
    Sampling == "Pre" ~ "Pre",
    TRUE ~ "On"
  ) ) |>  
  left_join(fusion_all, by = join_by("Run"))

Need_all <- Need |>
  mutate(cancer = "Pan-Cancer") |>
  mutate(irAE = if_else(rowSums(across(starts_with("irAE_days_"), ~ !is.na(.))) > 0, "yes", "no")) |>
  group_by(Sample_type, irAE) |>
  summarise(
    Sample_type = first(Sample_type),
    Numerator = n_distinct(Patient_ID[!is.na(fusion)]),
    Denominator = n_distinct(Patient_ID),
    Percentage = Numerator / Denominator * 100,
    irAE = first(irAE),
    irAE_type = "ALL"
  ) |>
  ungroup()
irAE_types <- c("Skin", "Flu-like", "Gastrointestinal", "Musculoskeletal", "Pulmonary", "Neurologic", "Thyroid")

results <- list()
for (irAE_type in irAE_types) {
  Need1 <- Need |>
    mutate(cancer = "Pan-Cancer") |>
    mutate(irAE = if_else(!is.na(get(paste0("irAE_days_", irAE_type))), "yes", "no")) |>
    group_by(Sample_type, irAE) |>
    summarise(
      Sample_type = first(Sample_type),
      Numerator = n_distinct(Patient_ID[!is.na(fusion)]),
      Denominator = n_distinct(Patient_ID),
      Percentage = Numerator / Denominator * 100,
      irAE = first(irAE),
      irAE_type = irAE_type
    ) |>
    ungroup()
  
  results[[irAE_type]] <- Need1
}

final_result <- bind_rows(results) |> rbind(Need_all)  

data <- final_result %>%
  mutate(group = paste0(Sample_type,"_",irAE))

colors <- c("Pre_no" = "#93c8c0", "On_no" = "#72c3a3", "Pre_yes" = "#f39da0", "On_yes" = "#e84445")

p <- ggplot(data, aes(x = group, y = Percentage, fill = group)) +
  geom_bar(stat = "identity", position = "dodge",width = 0.8,size = 1.2,alpha=0.7) +
  facet_wrap(~ irAE_type) +
  labs(x = "irAE Status",
       y = "Fusion Detetct Percentage",
       fill = "irAE Status") +
  scale_color_manual(values=colors)+ 
  scale_fill_manual(values=colors)+  
  theme_bw()+ 
  geom_signif(test = t.test,
              test.args = list(var.equal = TRUE, alternative = "two.sided"),
              size = 0.8, parse = TRUE) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(colour='black',size=15),
        axis.title.y = element_text(size = 15)
  ) 
p
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图7/sup7A.png"), 
       plot = p, dpi = 600, width = 10, height = 7)
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图7/sup7A.pdf"), 
       plot = p, dpi = 600, width = 10, height = 7)


#sup7E
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta) |>
  tidyr::separate_rows(gene1, gene1_type, sep = ",") |>
  tidyr::separate_rows(gene2, gene2_type, sep = ",") |>
  mutate(fusion = paste0(gene1, "::", gene2)) 
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv")
cohort_both <- cohort1 |> filter(time == "Pre/On") 
cohort_sum_both  <- (cohort_both$id) |> unique() |> setdiff("PHS001464")
both_all <- data.table()
for(cohort_id in cohort_sum_both){
  Need2 <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv"))
  if(cohort_id == "PRJNA356761"){
    Need2 <- Need2 |>
      mutate(Sample_type = case_when(
        grepl("Pre", Sample_ID)  ~ "Pre",
        grepl("On", Sample_ID)  ~ "On",
        TRUE ~ "On"
      ) )
  }else if(cohort_id == "PRJNA795330"){
    Need2 <- Need2 |>
      mutate(Sample_type = case_when(
        Sampling == "Pre" ~ "Pre",
        TRUE ~ "On"
      ) )
  }else if(cohort_id == "EGAD00001006282"){
    Need2 <- Need2 |>
      mutate(Sample_type = case_when(
        grepl("_SCREEN", Sample_ID) ~ "Pre",
        TRUE ~ "On"
      ) )
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
    Need2<- Need2 |>
      mutate(Sample_type = case_when(
        grepl("Pre", Response2) | grepl("Pre", Response) ~ "Pre",
        grepl("On",  Response2) | grepl("On",  Response) ~ "On",
        TRUE ~ NA
      ) )
  }
  Need1 <- Need2 |> select("Patient_ID", "Run", "Sample_type")
  Need1$Patient_ID <- as.character(Need1$Patient_ID)
  Need1$cohort <- cohort_id
  
  both_all <- rbind(Need1,both_all) 
  
}
both_all1 <-   both_all |> 
  filter(!is.na(Sample_type) ) |>
  group_by(Patient_ID) |>
  filter(all(c("Pre", "On") %in% Sample_type)) |>
  ungroup() 

cohort_select <- both_all1$cohort |> unique()
fusion_select <- fusion_all |>
  filter(cohort %in% cohort_select)|>
  select(cohort,Sample_ID, fusion) |>
  distinct() |> 
  rename(Run = Sample_ID)
both_all2 <-   both_all1 |> 
  left_join(fusion_select, by = join_by("Run","cohort")) |>
  mutate(cancer = case_when(
    cohort == "PRJNA940989" ~ "SGC",
    cohort == "PHS003316" ~ "BRCA",
    cohort == "PRJNA795330" ~ "Pan",
    TRUE ~ "SKCM"
  ))


#ALL  On 46.2%(177/383),  Pre 46.7(179/383)
#BRCA On 100%(3/3),       Pre 76.7(2/3)
#SGC  On 83.3%(14/17),    Pre 94.1%(16/17)
#SKCM On 69%(69/100),     Pre 77%(77/100)      
#Pan  On 34.6%(91/263),   Pre 31.9%(84/263) 
data <- data.frame(
  cancer = rep(c("ALL", "BRCA", "SGC", "SKCM", "Pan"), each = 2),
  Sample_type = rep(c("On", "Pre"), times = 5),
  Numerator = c(177, 179, 3, 2, 14, 16, 69, 77, 91, 84),
  Denominator = c(383, 383, 3, 3, 17, 17, 100, 100, 263, 263),
  Percentage = c(46.2, 46.7, 100, 76.7, 83.3, 94.1, 69, 77, 34.6, 31.9)
) |> arrange(desc(Sample_type))
data$Sample_type <- factor(data$Sample_type, levels = c("Pre", "On"))
colors <- c(Pre= "#DD5F60", On = "#7DDFD7")
p <- ggplot()+
  geom_bar(data= data,mapping=aes(x=cancer,y=Percentage,fill=Sample_type),
           size = 1.2,alpha=0.7,                
           position="dodge", stat="identity",width = 0.8)+    
  scale_color_manual(values=colors)+ 
  scale_fill_manual(values=colors)+  
  theme_bw()+ labs(y="Fusion detect rate", x="")+  
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(colour='black',size=12),
        axis.title.y = element_text(size = 10)
  ) +
  scale_x_discrete(labels = c("ALL(383)", "BRCA(3)","SGC(17)","Pan-Cancer(263)", "SKCM(100)"))
p
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图7/sup7E.png"), 
       plot = p, dpi = 600, width = 8, height = 5)
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图7/sup7E.pdf"), 
       plot = p, dpi = 600, width = 8, height = 5)

#sup7F
file_meta <- "C:/Users/Administrator/Desktop/fsdownload/metafusion_calls.star.tsv"
fusion_all <- fread(file_meta) |>
  tidyr::separate_rows(gene1, gene1_type, sep = ",") |>
  tidyr::separate_rows(gene2, gene2_type, sep = ",") |>
  mutate(fusion = paste0(gene1, "::", gene2)) 
cohort1 <- fread("C:/Users/Administrator/Desktop/fsdownload/cohorts_meta1.csv")
cohort_both <- cohort1 |> filter(time == "Pre/On") 
cohort_sum_both  <- (cohort_both$id) |> unique() |> setdiff("PHS001464")
both_all <- data.table()
for(cohort_id in cohort_sum_both){
  Need2 <- fread(paste0("C:/Users/Administrator/Desktop/fsdownload/Clininfo/",cohort_id,"_info.tsv"))
  if(cohort_id == "PRJNA356761"){
    Need2 <- Need2 |>
      mutate(Sample_type = case_when(
        grepl("Pre", Sample_ID)  ~ "Pre",
        grepl("On", Sample_ID)  ~ "On",
        TRUE ~ "On"
      ) )
  }else if(cohort_id == "PRJNA795330"){
    Need2 <- Need2 |>
      mutate(Sample_type = case_when(
        Sampling == "Pre" ~ "Pre",
        TRUE ~ "On"
      ) )
  }else if(cohort_id == "EGAD00001006282"){
    Need2 <- Need2 |>
      mutate(Sample_type = case_when(
        grepl("_SCREEN", Sample_ID) ~ "Pre",
        TRUE ~ "On"
      ) )
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
    Need2<- Need2 |>
      mutate(Sample_type = case_when(
        grepl("Pre", Response2) | grepl("Pre", Response) ~ "Pre",
        grepl("On",  Response2) | grepl("On",  Response) ~ "On",
        TRUE ~ NA
      ) )
  }
  Need1 <- Need2 |> select("Patient_ID", "Run", "Sample_type")
  Need1$Patient_ID <- as.character(Need1$Patient_ID)
  Need1$cohort <- cohort_id
  
  both_all <- rbind(Need1,both_all) 
  
}
both_all1 <-   both_all |> 
  filter(!is.na(Sample_type) ) |>
  group_by(Patient_ID) |>
  filter(all(c("Pre", "On") %in% Sample_type)) |>
  ungroup() 

cohort_select <- both_all1$cohort |> unique()
fusion_select <- fusion_all |>
  filter(cohort %in% cohort_select)|>
  select(cohort,Sample_ID, fusion) |>
  distinct() |> 
  rename(Run = Sample_ID)
both_all2 <-   both_all1 |> 
  left_join(fusion_select, by = join_by("Run","cohort")) |>
  mutate(cancer = case_when(
    cohort == "PRJNA940989" ~ "SGC",
    cohort == "PHS003316" ~ "BRCA",
    cohort == "PRJNA795330" ~ "Pan",
    TRUE ~ "SKCM"
  )) |> 
  as.data.frame() |> 
  select(Patient_ID, Sample_type, fusion, cancer) |> 
  mutate(fusion = case_when(
    is.na(fusion) ~ "Not detected",
    TRUE ~ fusion
  )) |>
  mutate(Sample_type = factor(Sample_type, levels = c("Pre","On")))

pair_stat <- both_all2 %>% 
  mutate(detected = fusion != "Not detected") %>%   
  group_by(Patient_ID, cancer, fusion) %>% 
  summarise(
    pre_det = any(Sample_type == "Pre" & detected),
    on_det  = any(Sample_type == "On"  & detected),
    .groups = "drop"
  ) %>% 
  mutate(flag = case_when(
    !pre_det & !on_det ~ "Neither",   # 都未有fusion检出
    pre_det & on_det   ~ "Both",      # 都有fusion检出且为同一个
    pre_det            ~ "Pre_only",  # 仅 Pre 检出了这个fusion，或者On没有检出
    TRUE               ~ "On_only"    # 仅 On  检出了这个fusion，或者Pre没有检出
  )) |> 
  select(cancer, flag) |> 
  as.data.table() 

cancer_flag_pct <- pair_stat[, .N, by = .(cancer, flag)
][, pct := round(N/sum(N)*100, 1), by = cancer][]


p <- ggplot(cancer_flag_pct, aes(x = cancer, y = pct, fill = flag)) +
  geom_col(width = 0.7, color = "black") +        # 堆叠柱
  scale_fill_manual(values = c("Both"      = "#1f77b4",
                               "Neither"   = "#7f7f7f",
                               "On_only"   = "#d62728",
                               "Pre_only"  = "#ff7f0e")) +
  labs(x = NULL,
       y = "Percentage (%)",
       title = "Fusion detection rate: Pre vs On-treatment by cancer type") +
  theme_minimal(base_size = 14) +
  theme( panel.grid = element_blank(),
         plot.background = element_blank(),
         axis.line = element_blank(),
         axis.text = element_text(color = "black"),
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
         axis.ticks = element_blank(),
         legend.key.spacing.y = unit(0.2, "cm"),
         legend.key.size = unit(1, 'cm'), 
         legend.title = element_text(size=16),
         legend.text = element_text(size=14),
         legend.position = "right")
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图7/sup7F.png"), 
       plot = p, dpi = 600, width = 8, height = 5)
ggsave(filename = paste0("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/附图/附图7/sup7F.pdf"), 
       plot = p, dpi = 600, width = 8, height = 5)