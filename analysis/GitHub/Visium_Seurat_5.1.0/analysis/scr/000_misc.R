# AIM ---------------------------------------------------------------------
# misc task from martina

# LIBRARIES ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
# library(SeuratData)
library(patchwork)
library(harmony)
library(ggrepel)
library(ComplexHeatmap)
library(pals)
library(scales)
library(SeuratWrappers)
library(presto)
library(glmGamPoi)


# provide the number of spots per manual annotation -----------------------
# read in the data
list_brain <- readRDS("../../Visium_brainMS_single_sample/out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")

LUT_sample <- read_csv("../../Visium_brainMS_single_sample/data/sample_classification.csv") %>%
  mutate(dataset = paste0("brain_",sample))

# pull all the metadata in one
df_meta_tot <- lapply(list_brain, function(sl){
  test <- sl@meta.data %>%
    rownames_to_column("barcode")
  return(test)
}) %>%
  bind_rows()

df_meta_tot2 <- df_meta_tot %>%
  # add more info based on the manual segmentation
  # Martina asked to gather the periplaques categories in one.
  mutate(manual_segmentation2 = case_when(manual_segmentation %in% c("periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6") ~ "WM",
                                          TRUE ~ manual_segmentation)) %>%
  
  left_join(LUT_sample,by="dataset") %>%
  mutate(manual_segmentation3 = case_when(manual_segmentation2 == "edge" ~ paste0("edge_",type),
                                          TRUE ~ manual_segmentation2))

# count the number of spots per manual annotation
df_summary_spots <- df_meta_tot2 %>%
  group_by(sample,type,manual_segmentation3) %>%
  summarise(n_spot = n(),.groups = "drop") %>%
  group_by(sample,type) %>%
  mutate(n_spot_tot = sum(n_spot))

# save the table
write_tsv(df_summary_spots,"../out/table/000_df_summary_spots_Visium.tsv")

