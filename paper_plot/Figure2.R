### Figure2a
library(sankeyD3)
library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(tidyr)

#Figure2a
Donut_Chart_data <- readRDS("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/Donut_Chart_data.rds")
p_figure2a<-plot_ly(Donut_Chart_data ,
                    ids = ~ids,
                    labels = ~labels,
                    parents= ~parents,
                    type = 'sunburst',
                    marker = list(colors = Donut_Chart_data$colors))
save_image(p_figure2a, file ="C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/Figure2a.png", 
           height = 1000, width = 1000)

#Figure2b
Donut_Chart_data_gdc <- readRDS("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Donut_Chart_data_gdc.rds")
p_figure2b <-plot_ly(Donut_Chart_data_gdc ,
                    ids = ~ids,
                    labels = ~labels,
                    parents= ~parents,
                    type = 'sunburst',
                    marker = list(colors = Donut_Chart_data_gdc$colors))
p_figure2b
save_image(p_figure2b, file ="C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure3/Figure2b.png", 
           height = 1000, width = 1000)


### Figure 2c
all_data_Run_pro <- fread("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/sample_id.tsv")
data_figure2_pre <- fread("C:/Users/Administrator/Desktop/fsdownload/fusion_distribution.tsv" , sep = "\t") |>
                    separate_rows(gene1, sep = ",") |>
                    separate_rows(gene2, sep = ",") |>
                    mutate(fusion = paste0(gene1, "::", gene2)) |>
                    select(c(fusion,Sample_ID,cohort)) |>
                    distinct() 
data_figure2 <- data_figure2_pre |>
  mutate(
    cancer_type = case_when(
      Sample_ID %in% all_data_Run_pro$Run ~ 
      all_data_Run_pro$cancer_type[match(Sample_ID, all_data_Run_pro$Run)],
      TRUE ~ NA_character_    
    )) |> 
  select(c(cancer_type, fusion)) |>
  distinct() |>
  group_by(cancer_type)|>
  summarise(
    n_rows = n(),  
    .groups = 'drop'  
  ) |>
  arrange(desc(n_rows)) 

p1 <- ggplot(data= data_figure2,aes(x = reorder(cancer_type, -n_rows, FUN = sum), y = n_rows)) +
  geom_bar(stat="identity",position = "dodge", fill = "#1c3c63")+
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),  
        axis.title.y = element_text(size = 30),
        legend.title = element_blank()
  )+
  guides(fill = "none") +
  labs(x=NULL,y="Fusions Counts") 
p1

ggsave("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/Figure2c.pdf", 
       plot = p1,         
       width = 14,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "pdf",           
       bg = "white")  
ggsave("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/Figure2c.png", 
       plot = p1,         
       width = 14,                
       height = 7,                
       dpi = 600,                
       units = "in",              
       device = "png",           
       bg = "white")  


### Figure2d
### Tumor data
Tumor_data <- fread("C:/Users/Administrator/Desktop/fsdownload/fusion_distribution.tsv" , sep = "\t") |>
  filter(sample_type == "Tumor" | sample_type == "Blood")  |> 
  separate_rows(gene1, sep = ",") |>
  separate_rows(gene2, sep = ",") |>
  mutate(fusion = paste0(gene1, "::", "gene2")) |>
  select(c(fusion,cancer_db_hits)) |>
  distinct() 

group_chimer      <- c("ChimerKB", "ChimerPub", "ChimerSeq", "chimerdb_omim", "chimerdb_pubmed")
group_tcga        <- c("DEEPEST2019", "YOSHIHARA_TCGA", "Larsson_TCGA", "TCGA_StarF2019", "TumorFusionsNAR2018", "GUO2018CR_TCGA")
group_ccle        <- c("Klijn_CellLines", "CCLE_StarF2019", "DepMap2023")
group_karyotype   <- "Mitelman"
group_reviews     <- c("Cosmic", "HaasMedCancer")
tumor_markers <- c("Mitelman", "chimerdb_omim", "chimerdb_pubmed", "Cosmic", "YOSHIHARA_TCGA", 
                   "ChimerKB", "ChimerPub", "ChimerSeq", "Klijn_CellLines", "Larsson_TCGA", 
                   "CCLE_StarF2019", "HaasMedCancer", "TCGA_StarF2019", "TumorFusionsNAR2018", 
                   "DEEPEST2019", "GUO2018CR_TCGA", "DepMap2023")
gene_up_tumor <- Tumor_data |>
  distinct() |> 
  mutate(
    Group_TCGA = ifelse(rowSums(sapply(group_tcga, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_CCLE = ifelse(rowSums(sapply(group_ccle, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Karyotype = ifelse(rowSums(sapply(group_karyotype, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Chimer = ifelse(rowSums(sapply(group_chimer, grepl, x = cancer_db_hits))> 0, 1, 0),
    Group_Reviews = ifelse(rowSums(sapply(group_reviews, grepl, x = cancer_db_hits)) > 0, 1, 0),
    Group_Not_Reported = ifelse(rowSums(sapply(tumor_markers, grepl, x = cancer_db_hits)) == 0, 1, 0),
  )  
TCGA <- gene_up_tumor |> 
  select(fusion, Group_TCGA) |>
  filter(Group_TCGA == 1) |>
  rename(Group = Group_TCGA)|>
  mutate(Group = "Group_TCGA") 
CCLE <- gene_up_tumor |> 
  select(fusion, Group_CCLE)|>
  filter(Group_CCLE == 1) |>
  rename(Group = Group_CCLE) |>
  mutate(Group = "Group_CCLE") 
Karyotype <- gene_up_tumor |> 
  select(fusion, Group_Karyotype)|>
  filter(Group_Karyotype == 1) |>
  rename(Group = Group_Karyotype) |>
  mutate(Group = "Group_Karyotype") 
Chimer <- gene_up_tumor |> select(fusion, Group_Chimer )|>
  filter(Group_Chimer == 1) |>
  rename(Group = Group_Chimer) |>
  mutate(Group = "Group_Chimer") 
Reviews <- gene_up_tumor |> select(fusion, Group_Reviews)|>
  filter(Group_Reviews == 1) |>
  rename(Group = Group_Reviews) |>
  mutate(Group = "Group_Reviews") 
Not_ReportedT <- gene_up_tumor |> 
  select("fusion", "Group_Not_Reported") |>
  filter(Group_Not_Reported == 1) |>
  rename(Group = Group_Not_Reported) |>
  mutate(Group = "Group_Not_Reported") 

Group_All <- rbind(TCGA, CCLE, Karyotype, Chimer, Reviews, Not_ReportedT)
labels <- c("Group_Karyotype","Group_Reviews","Group_CCLE","Group_Chimer","Group_TCGA","Group_Not_Reported")
Group_All <- Group_All %>%
  mutate(Group = factor(Group, levels = labels)) %>%
  arrange(Group)
lists <- lapply(split(Group_All, Group_All$Group), function(x) x$fusion)
m <- make_comb_mat(lists)
cs <- comb_size(m) 
colors <- c("#f3a59a",  "#a6ddea",  "#80d0c3", "#9eaac4", "#c2c8da", "#d8cdc1")
pdf("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/Figure2d.pdf", width = 8, height = 5)
p_tumor = draw(ComplexHeatmap::UpSet(m,                  
                                     top_annotation = upset_top_annotation(m,  
                                                                           ylim = c(0, max(cs)*1.5),
                                                                           annotation_name_rot = 90,                                                      
                                                                           annotation_name_side = "left",height = unit(4, "cm")),               
                                     bg_col = rev(colors),              
                                     set_order = labels)) 
od = column_order(p_tumor)
decorate_annotation("intersection_size", {  
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2., "pt"), default.units = "native", just = "left", gp = gpar(fontsize = 10),rot = 90)})

dev.off()
png("C:/Users/Administrator/Nutstore/1/IO-GeneFusion/投稿准备/Figure2/Figure2d.png", width = 8, height = 4, units = "in", res = 600)
p_tumor = draw(ComplexHeatmap::UpSet(m,                  
                                     top_annotation = upset_top_annotation(m,  
                                                                           ylim = c(0, max(cs)*1.5),
                                                                           annotation_name_rot = 90,                                                      
                                                                           annotation_name_side = "left",height = unit(4, "cm")),               
                                     bg_col = rev(colors),              
                                     set_order = labels)) 
od = column_order(p_tumor)
decorate_annotation("intersection_size", {  
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2., "pt"), default.units = "native", just = "left", gp = gpar(fontsize = 10),rot = 90)})

dev.off()

#Figure2E
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
  p_Figure2E_part1 <- plot_ly(ids = ids, 
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
  save_image(p_Figure2E_part1, file = paste0("C:/Users/Administrator/Desktop/Figure2E_data/", cancer_type1,"_Figure2E.pdf"), height = 1000, width = 1000)
  save_image(p_Figure2E_part1, file = paste0("C:/Users/Administrator/Desktop/Figure2E_data/", cancer_type1, "_Figure2E.png"), height = 1000, width = 1000)
}