# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(hdf5r)
library(limma)
library(future)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
list_spotlight <- readRDS("out/object/list_brain_all_spotlight.rds")

# extract the the metadata of the spatial dataset
list_brain <- pmap(list(list_spotlight,names(list_spotlight)), function(x,name){
  
  # # expression values
  # df_exp <- x@assays$SCT@data["TSPO",] %>% 
  #   data.frame(exp = .) %>% 
  #   rownames_to_column("barcodes") %>% 
  #   mutate(gene = "TSPO") %>% 
  #   mutate(dataset = paste0("brain_",name))
  
  # extract the metadata for the dataset
  df_meta <- x@meta.data %>% 
    rownames_to_column("barcodes") %>% 
    mutate(dataset = paste0("brain_",name))
  
  # return(list(exp = df_exp,
  #             meta = df_meta))
  return(list(meta = df_meta))
})

# build a full matrix for the expression information on TSPO
df_dec <- lapply(list_brain,function(x){
  x$meta
}) %>% 
  bind_rows()

dim(df_dec)

# laod the segmentation from Sofia ----------------------------------------
folder3 <- "data/segmentation/sofia/all/draft_02/"
file_id3 <- dir(folder3) %>% 
  str_subset("sofia_segmentation") %>% 
  str_remove_all(pattern = "_sofia_segmentation_no_vessel.csv")

# build a full matrix for the expression information on TSPO
df_sofi_segmentation <- lapply(file_id3,function(x){
  
  # keep track of the progress of the reading
  print(x)
  
  df_meta <- read_csv(paste0("data/segmentation/sofia/all/draft_02/",x,"_sofia_segmentation_no_vessel.csv")) %>% 
    mutate(dataset = paste0("brain_",x)) %>% 
    # some point do not have an annotation fiter them out
    filter(!is.na(sofia_segmentation))
  
  return(df_meta)
  
}) %>% 
  bind_rows()

dim(df_sofi_segmentation)

# laod the TSPO threshold from Sofia --------------------------------------
folder4 <- "data/segmentation/sofia/all/tspo_pos/"
file_id4 <- dir(folder4) %>% 
  str_remove_all(pattern = "tspo|.csv")

# build a full matrix for the expression information on TSPO
df_sofi_tspo <- lapply(file_id4,function(x){
  
  # keep track of the progress of the reading
  print(x)
  
  df_meta <- read_csv(paste0(folder4,"tspo",x,".csv")) %>% 
    mutate(dataset = paste0("brain_V",x)) %>% 
    # some point do not have an annotation fiter them out
    filter(!is.na(tspo_class))
  
  return(df_meta)
  
}) %>% 
  bind_rows()

dim(df_sofi_tspo)

table(df_sofi_tspo$tspo_class)

# add the annotation of the sample type based on the slide
LUT_sample <- read_csv("data/sample_classification.csv") %>% 
  mutate(dataset = paste0("brain_",sample))

# join the datasets
df_tot <- left_join(df_dec,LUT_sample,by = "dataset") %>% 
  inner_join(df_sofi_segmentation,c("barcodes","dataset")) %>%
  inner_join(df_sofi_tspo,c("barcodes","dataset")) %>%
  # pivot_longer(cols = AST:VAS,names_to = "cell_type",values_to = "prop") %>% 
  # after chatting with sofia she decided to remove the BASAL sample
  filter(type != "BASAL")

table(df_tot$sofia_segmentation)

# fix the meta
df_tot_fix <- df_tot %>% 
  mutate(sofia_segmentation = case_when(sofia_segmentation == "vessels"~"vessel",
                                        sofia_segmentation == "WM_1layer"~"WM_01_layer",
                                        sofia_segmentation == "WM_2layer"~"WM_02_layer",
                                        sofia_segmentation == "WM_3layer"~"WM_03_layer",
                                        sofia_segmentation == "WM_4lyer"~"WM_04_layer",
                                        sofia_segmentation == "WM_4layer"~"WM_04_layer",
                                        sofia_segmentation == "WM_5layer"~"WM_05_layer",
                                        sofia_segmentation == "WM_6layer"~"WM_06_layer",
                                        sofia_segmentation == "WM_7layer"~"WM_07_layer",
                                        sofia_segmentation == "WM_8layer"~"WM_08_layer",
                                        sofia_segmentation == "WM_9layer"~"WM_09_layer",
                                        sofia_segmentation == "WM_10layer"~"WM_10_layer",
                                        T~sofia_segmentation)) %>% 
  mutate(sofia_segmentation_short = case_when(str_detect(sofia_segmentation,pattern = "WM_")~"WM",
                                              T~sofia_segmentation))

# check the metadata
table(df_tot_fix$sofia_segmentation,useNA="ifany")
table(df_tot_fix$sofia_segmentation_short,useNA="ifany")
table(df_tot_fix$tspo_class,useNA="ifany")
sum(table(df_tot_fix$sofia_segmentation_short))

# save the dataset to share
df_tot_fix %>% 
  write_tsv("out/table/df_tot_TSPO_sofia.tsv")

# trimm only the column of interest
df_tot_fix %>% 
  dplyr::select(barcodes,dataset,class_top,type,sofia_segmentation_short,tspo_class) %>% 
  write_tsv("out/table/df_tot_TSPO_sofia_minimal_column.tsv")

# summarise
df_tot_fix %>% 
  dplyr::select(barcodes,dataset,class_top,type,sofia_segmentation_short,tspo_class) %>% 
  group_by(dataset,type,tspo_class) %>% 
  summarise(count_spot = n()) %>% 
  write_tsv("out/table/df_tot_TSPO_sofia_summary_dataset_tspo.tsv")

df_tot_fix %>% 
  dplyr::select(barcodes,dataset,class_top,type,sofia_segmentation_short,tspo_class) %>% 
  group_by(dataset,type,class_top,tspo_class) %>% 
  summarise(count_spot = n()) %>% 
  write_tsv("out/table/df_tot_TSPO_sofia_summary_dataset_cellType_tspo.tsv")

df_tot_fix %>% 
  dplyr::select(barcodes,dataset,class_top,type,sofia_segmentation_short,tspo_class) %>% 
  group_by(dataset,type,class_top,sofia_segmentation_short,tspo_class) %>% 
  summarise(count_spot = n()) %>% 
  write_tsv("out/table/df_tot_TSPO_sofia_summary_dataset_cellType_sofiaSegmentation_tspo.tsv")

df_tot_fix %>% 
  dplyr::select(barcodes,dataset,class_top,type,sofia_segmentation_short,tspo_class) %>% 
  group_by(dataset,type,tspo_class) %>% 
  summarise(count_spot = n()) %>% 
  
  ungroup() %>% 
  summarise(s = sum(count_spot))
