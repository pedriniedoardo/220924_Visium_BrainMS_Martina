# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)

# # where class_top is immune, what is the predominant population at the edge?
# df_meta_subset <- read_tsv("out/table/test_IMM_whole_slide_SPOTLight.tsv")
# 
# df_class_top_IMM <- df_meta_subset %>%
#   dplyr::select(barcodes,MG.Homeo:MIMS.iron) %>%
#   pivot_longer(names_to = "class_top_IMM",values_to = "score",-barcodes) %>%
#   group_by(barcodes) %>%
#   arrange(desc(score)) %>%
#   dplyr::slice(1) %>%
#   # top_n(wt = score,n = 1) %>%
#   dplyr::select(-score)
# 
# meta_fix <- left_join(df_meta_subset,df_class_top_IMM,"barcodes") %>% 
#   # add the new catrgory based on the deconvolution of the IMM class
#   mutate(class_top2 = case_when(class_top == "IMM"~class_top_IMM,
#                                 T~class_top))

# read in data ------------------------------------------------------------
# read in the deconvolution result from each slice
folder2 <- "out/table/"
file_id2 <- dir(folder2) %>% 
  str_subset("SPOTLight_ASTRO") %>% 
  str_remove_all(pattern = "_SPOTLight_ASTRO_WholeSlide.tsv")

# build a full matrix for the expression information on TSPO
list_meta_dec <- lapply(file_id2,function(x){
  
  # keep track of the progress of the reading
  print(x)
  
  df_meta <- read_tsv(paste0("out/table/",x,"_SPOTLight_ASTRO_WholeSlide.tsv")) %>% 
    mutate(dataset = paste0("brain_",x))
  
  return(df_meta)
  
}) %>% 
  setNames(file_id2)

# read in the processed data of the brain slices
list_brain <- readRDS("out/object/list_brain_all.rds")

# add the spotlight annotation --------------------------------------------
list_brain_spotlight_ASTRO <- lapply(file_id2,function(x){
  # check the progress
  print(x)
  
  # read in the seurat object -----------------------------------------------
  test <- list_brain[[x]]
  
  # read in the metadata ----------------------------------------------------
  # wrangling pick the top class per score
  df_class_top_ASTRO <- list_meta_dec[[x]] %>% 
    dplyr::select(barcodes,AIMS:ASTRO.sen) %>%
    pivot_longer(names_to = "class_top_ASTRO",values_to = "score",-barcodes) %>%
    group_by(barcodes) %>%
    arrange(desc(score)) %>%
    dplyr::slice(1) %>%
    # top_n(wt = score,n = 1) %>%
    dplyr::select(-score)
  
  # join the two metadata
  meta_fix <- left_join(list_meta_dec[[x]],df_class_top_ASTRO,"barcodes") %>%
    # add the new catrgory based on the deconvolution of the IMM class
    mutate(class_top3 = case_when(class_top == "AST"~class_top_ASTRO,
                                  T~class_top))
  
  # add the meta from the spotlight
  # class(list_brain$V01@meta.data)
  test@meta.data <- meta_fix %>%
    column_to_rownames("barcodes")
  
  return(test)
}) %>% 
  setNames(file_id2)

# save a list of all the preprocessed brain slices with spotligth annotation
saveRDS(list_brain_spotlight_ASTRO,"out/object/list_brain_spotlight_ASTRO.rds")
list_brain_spotlight_ASTRO <- readRDS("out/object/list_brain_spotlight_ASTRO.rds")

# loop all the plotting ---------------------------------------------------
# x <- file_id2[1]
lapply(file_id2,function(x){
  # check the progress
  print(x)
  
  # read in the seurat object -----------------------------------------------
  test <- list_brain_spotlight_ASTRO[[x]]
  
  # make the different classes as categorical
  # SpatialFeaturePlot(test, features = "IMM") + theme(legend.position = "right")
  pdf(paste0("out/image/10_class_deconvolution_ASTROSubset_",x,".pdf"),width = 15,height = 15)
  
  # reorder the levels
  test$class_top3 <- factor(test$class_top3,levels = c("NEU","OLIGO","OPC","IMM","VAS","LYM","ASTRO.nr","ASTRO.reactive","ASTRO.sen","AIMS"))
  
  Idents(test) <- "class_top3"
  print(SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE, ncol = 3,alpha = 0.5) +
          plot_annotation(title = "class_top3"))
  dev.off()
})

# x <- "V01"
# # try to change the alpha of the plotting
# lapply(file_id2,function(x){
#   # check the progress
#   print(x)
#   
#   # read in the seurat object -----------------------------------------------
#   test <- list_brain_spotlight_ASTRO[[x]]
#   
#   # make the different classes as categorical
#   # SpatialFeaturePlot(test, features = "IMM") + theme(legend.position = "right")
#   pdf(paste0("out/image/10_class_deconvolution_ASTROSubset_",x,"_2.pdf"),width = 15,height = 15)
#   
#   # reorder the levels
#   test$class_top3 <- factor(test$class_top3,levels = c("NEU","OLIGO","OPC","IMM","VAS","LYM","ASTRO.nr","ASTRO.reactive","ASTRO.sen","AIMS"))
#   
#   Idents(test) <- "class_top3"
#   print(SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE, ncol = 3,alpha = 0.8) +
#           plot_annotation(title = "class_top3"))
#   dev.off()
# })

# see correlation with location -------------------------------------------
# load the dataset
list_spotlight <- readRDS("out/object/list_brain_spotlight_ASTRO.rds")

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

# add the annotation of the sample type based on the slide
LUT_sample <- read_csv("data/sample_classification.csv") %>% 
  mutate(dataset = paste0("brain_",sample))

# join the datasets
df_tot <- left_join(df_dec,LUT_sample,by = "dataset") %>% 
  inner_join(df_sofi_segmentation,c("barcodes","dataset")) %>%
  # inner_join(df_sofi_tspo,c("barcodes","dataset")) %>%
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
  write_tsv("out/table/df_tot_SPOTLight_ASTRO_sofia.tsv")

df_tot_fix %>% 
  group_by(dataset,type) %>% 
  summarise(n = n())

# plot the summarisation --------------------------------------------------
# dotplot
df_tot_fix %>% 
  dplyr::select(barcodes,dataset,class_top3,type,sofia_segmentation_short) %>% 
  group_by(dataset,type,class_top3) %>% 
  summarise(count_spot = n()) %>% 
  ungroup() %>% 
  group_by(dataset) %>% 
  mutate(tot = sum(count_spot)) %>% 
  ungroup() %>% 
  mutate(prop = count_spot/tot) %>% 
  mutate(type = factor(type,levels = c("NAWM","CORTEX","CI","A","CA","RM"))) %>% 
  ggplot(aes(x=type,y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),alpha=0.6) +
  facet_wrap(~class_top3,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("out/image/summary_SPOTLight_ASTRO_TypeCellType_dotplot.pdf",width = 10,height = 10)

# the same as above as an heatmap
mat_02 <- df_tot_fix %>% 
  dplyr::select(barcodes,dataset,class_top3,type,sofia_segmentation_short) %>% 
  group_by(type,class_top3) %>% 
  summarise(count_spot = n()) %>% 
  ungroup() %>% 
  group_by(type) %>% 
  mutate(tot = sum(count_spot)) %>% 
  ungroup() %>% 
  mutate(prop = count_spot/tot) %>% 
  ungroup() %>% 
  group_by(class_top3) %>% 
  mutate(norm = (prop-mean(prop))/sd(prop)) %>% 
  dplyr::select(type,class_top3,norm) %>% 
  pivot_wider(names_from = class_top3,values_from = norm) %>% 
  column_to_rownames("type")

# filter the matrix to avoid too many NA
id_02 <- is.na(mat_02) %>%
  colSums() < 2
id2_02 <- is.na(mat_02[,id_02]) %>%
  rowSums() < 2

mat_02_filter <- mat_02[id2_02,id_02]

# build the annotation on the heatmap
# sample_ordered <- data.frame(sample = str_extract(colnames(mat_short_filter),pattern = c("V\\d+"))) %>% 
#   left_join(LUT_sample)
# 
# table(sample_ordered$type)
# 
# sample_ordered_conversion <- sample_ordered %>% pull(type)
# area_ordered <- str_extract(colnames(mat_short_filter),pattern = c("WM|core|cortex|edge|periplaque|active lesion"))
# 
# column_ha <- HeatmapAnnotation(sample = sample_ordered_conversion,
#                                area = area_ordered,
#                                col = list(area = c("WM"="beige",
#                                                    "core"="red",
#                                                    "cortex"="gray",
#                                                    "edge"="orange",
#                                                    "periplaque"="green",
#                                                    "active lesion"="yellow"),
#                                           sample = c("NAWM"="cyan",
#                                                      "BASAL"="gray",
#                                                      "CORTEX"="green",
#                                                      "CA"="red",
#                                                      "A"="orange",
#                                                      "CI"="yellow",
#                                                      "RM"="blue")))
pdf("out/image/summary_SPOTLight_ASTRO_TypeCellType_heatmap.pdf",width = 5,height = 4)
# Heatmap(mat_01_filter,name = "prop",
#         top_annotation = column_ha)
Heatmap(mat_02_filter,name = "Z score")
dev.off()

# -------------------------------------------------------------------------
# dotplot
df_tot_fix %>% 
  dplyr::select(barcodes,dataset,class_top3,type,sofia_segmentation_short) %>% 
  group_by(dataset,sofia_segmentation_short,class_top3) %>% 
  summarise(count_spot = n()) %>% 
  ungroup() %>% 
  group_by(dataset,sofia_segmentation_short) %>% 
  mutate(tot = sum(count_spot)) %>% 
  ungroup() %>% 
  mutate(prop = count_spot/tot) %>% 
  mutate(sofia_segmentation_short = factor(sofia_segmentation_short,levels = c("WM","cortex","periplaque","edge","core","core (RM)","active lesion"))) %>% 
  ggplot(aes(x=sofia_segmentation_short,y=prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1),alpha=0.6) +
  facet_wrap(~class_top3,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("out/image/summary_SPOTLight_ASTRO_LocationCellType_dotplot.pdf",width = 10,height = 10)

# # the same as above as an heatmap
# mat_01 <- df_tot_fix %>% 
#   dplyr::select(barcodes,dataset,class_top2,type,sofia_segmentation_short) %>% 
#   group_by(type,class_top2) %>% 
#   summarise(count_spot = n()) %>% 
#   ungroup() %>% 
#   group_by(type) %>% 
#   mutate(tot = sum(count_spot)) %>% 
#   ungroup() %>% 
#   mutate(prop = count_spot/tot) %>% 
#   dplyr::select(type,class_top2,prop) %>% 
#   pivot_wider(names_from = class_top2,values_from = prop) %>% 
#   column_to_rownames("type")
# 
# # filter the matrix to avoid too many NA
# id_01 <- is.na(mat_01) %>%
#   colSums() == 0
# id2_01 <- is.na(mat_01[,id_01]) %>%
#   rowSums() == 0
# 
# mat_01_filter <- mat_01[id2_01,id_01]
# 
# # build the annotation on the heatmap
# # sample_ordered <- data.frame(sample = str_extract(colnames(mat_short_filter),pattern = c("V\\d+"))) %>% 
# #   left_join(LUT_sample)
# # 
# # table(sample_ordered$type)
# # 
# # sample_ordered_conversion <- sample_ordered %>% pull(type)
# # area_ordered <- str_extract(colnames(mat_short_filter),pattern = c("WM|core|cortex|edge|periplaque|active lesion"))
# # 
# # column_ha <- HeatmapAnnotation(sample = sample_ordered_conversion,
# #                                area = area_ordered,
# #                                col = list(area = c("WM"="beige",
# #                                                    "core"="red",
# #                                                    "cortex"="gray",
# #                                                    "edge"="orange",
# #                                                    "periplaque"="green",
# #                                                    "active lesion"="yellow"),
# #                                           sample = c("NAWM"="cyan",
# #                                                      "BASAL"="gray",
# #                                                      "CORTEX"="green",
# #                                                      "CA"="red",
# #                                                      "A"="orange",
# #                                                      "CI"="yellow",
# #                                                      "RM"="blue")))
# pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_cor_TSPO_pearson_SofiaShort.pdf",width = 10,height = 4)
# # Heatmap(mat_01_filter,name = "prop",
# #         top_annotation = column_ha)
# Heatmap(mat_01_filter,name = "prop")
# dev.off()

# the same as above as an heatmap
mat_12 <- df_tot_fix %>% 
  dplyr::select(barcodes,dataset,class_top3,type,sofia_segmentation_short) %>% 
  group_by(sofia_segmentation_short,class_top3) %>% 
  summarise(count_spot = n()) %>% 
  ungroup() %>% 
  group_by(sofia_segmentation_short) %>% 
  mutate(tot = sum(count_spot)) %>% 
  ungroup() %>% 
  mutate(prop = count_spot/tot) %>% 
  ungroup() %>% 
  group_by(class_top3) %>% 
  mutate(norm = (prop-mean(prop))/sd(prop)) %>% 
  dplyr::select(sofia_segmentation_short,class_top3,norm) %>% 
  pivot_wider(names_from = class_top3,values_from = norm) %>% 
  column_to_rownames("sofia_segmentation_short")

# filter the matrix to avoid too many NA
id_12 <- is.na(mat_12) %>%
  colSums() < 5
id2_12 <- is.na(mat_12[,id_12]) %>%
  rowSums() < 4

mat_12_filter <- mat_12[id2_12,id_12]

# build the annotation on the heatmap
# sample_ordered <- data.frame(sample = str_extract(colnames(mat_short_filter),pattern = c("V\\d+"))) %>% 
#   left_join(LUT_sample)
# 
# table(sample_ordered$type)
# 
# sample_ordered_conversion <- sample_ordered %>% pull(type)
# area_ordered <- str_extract(colnames(mat_short_filter),pattern = c("WM|core|cortex|edge|periplaque|active lesion"))
# 
# column_ha <- HeatmapAnnotation(sample = sample_ordered_conversion,
#                                area = area_ordered,
#                                col = list(area = c("WM"="beige",
#                                                    "core"="red",
#                                                    "cortex"="gray",
#                                                    "edge"="orange",
#                                                    "periplaque"="green",
#                                                    "active lesion"="yellow"),
#                                           sample = c("NAWM"="cyan",
#                                                      "BASAL"="gray",
#                                                      "CORTEX"="green",
#                                                      "CA"="red",
#                                                      "A"="orange",
#                                                      "CI"="yellow",
#                                                      "RM"="blue")))
pdf("out/image/summary_SPOTLight_ASTRO_LocationCellType_heatmap.pdf",width = 5,height = 4)
# Heatmap(mat_01_filter,name = "prop",
#         top_annotation = column_ha)
Heatmap(mat_12_filter,name = "Z score")
dev.off()

saveRDS(df_tot_fix,"out/object/df_tot_fix_ASTRO_dec.rds")
