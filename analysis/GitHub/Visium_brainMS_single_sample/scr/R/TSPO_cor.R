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
library(ggrepel)

# # single sample loading ---------------------------------------------------
# brain <- readRDS("../out_large/brain_V01.rds")
# 
# df_goi <- SpatialFeaturePlot(brain, features = c("TSPO"))
# df_goi$data %>% 
#   rownames_to_column("barcodes") %>% 
#   filter(barcodes%in%c("ACGTTCACTATGCCGC-1"))
# 
# # dim(brain@assays$Spatial@data)
# # 
# # brain@assays$Spatial@data["TSPO",] %>% 
# #   data.frame(exp = .) %>% 
# #   rownames_to_column("barcodes") %>% 
# #   filter(barcodes%in%c("ACGTTCACTATGCCGC-1"))
# 
# dim(brain@assays$SCT@data)
# 
# df_exp <- brain@assays$SCT@data["TSPO",] %>% 
#   data.frame(exp = .) %>% 
#   rownames_to_column("barcodes") %>% 
#   mutate(gene = "TSPO")
# 
# # confirm the data is the same of the one collected from the SpatilaFeaturePlot funciton
# df_exp %>% 
#   filter(barcodes%in%c("ACGTTCACTATGCCGC-1"))
# 
# # meta deconvolution
# df_meta <- read_tsv("out/table/V01_SPOTLight.tsv")
# 
# df_tot <- left_join(df_meta,df_exp,"barcodes") %>% 
#   pivot_longer(cols = AST:VAS,names_to = "cell_type",values_to = "prop")
# 
# # df_tot %>% 
# #   filter(barcodes%in%c("AAACAAGTATCTCCCA-1"))
# # 
# # df_meta %>% 
# #   filter(barcodes%in%c("AAACAAGTATCTCCCA-1")) %>% 
# #   pivot_longer(cols = AST:VAS,names_to = "cell_type",values_to = "prop") %>% 
# #   pull(prop)
# 
# # # AAACACCAATAACTGC-1
# # df_meta %>% 
# #   filter(barcodes%in%c("AAACACCAATAACTGC-1")) %>% 
# #   pivot_longer(cols = AST:VAS,names_to = "cell_type",values_to = "prop") %>% 
# #   pull(prop) %>% sum()
# 
# df_tot %>% 
#   ggplot(aes(x=prop,y=exp))+geom_point()+facet_wrap(~cell_type)+geom_smooth()
# 
# mat <- df_tot %>% 
#   group_by(cell_type) %>% 
#   summarise(cor_TSPO = cor(prop,exp)) %>% 
#   column_to_rownames("cell_type")
# 
# Heatmap(mat)  

# read in the data --------------------------------------------------------
list_spotlight <- readRDS("out/object/list_brain_all_spotlight.rds")

# extract the expression values of TSPO and of the metadata
list_brain <- pmap(list(list_spotlight,names(list_spotlight)), function(x,name){
  
  # expression values
  df_exp <- x@assays$SCT@data["TSPO",] %>% 
    data.frame(exp = .) %>% 
    rownames_to_column("barcodes") %>% 
    mutate(gene = "TSPO") %>% 
    mutate(dataset = paste0("brain_",name))
  
  # extract the metadata for the dataset
  df_meta <- x@meta.data %>% 
    rownames_to_column("barcodes") %>% 
    mutate(dataset = paste0("brain_",name))
  
  return(list(exp = df_exp,
              meta = df_meta))
})

# save the two elements as independent tables
df_exp <- lapply(list_brain,function(x){
  x$exp
}) %>% 
  bind_rows()

dim(df_exp)

# build a full matrix for the expression information on TSPO
df_dec <- lapply(list_brain,function(x){
  x$meta
}) %>% 
  bind_rows()

dim(df_dec)

# add the annotation of the sample type based on the slide
LUT_sample <- read_csv("data/sample_classification.csv") %>% 
  mutate(dataset = paste0("brain_",sample))

# join the two datasets
df_tot <- left_join(df_dec,LUT_sample,by = "dataset") %>% 
  left_join(df_exp,c("barcodes","dataset")) %>% 
  pivot_longer(cols = AST:VAS,names_to = "cell_type",values_to = "prop") %>% 
  # after chatting with sofia she decided to remove the BASAL sample
  filter(type != "BASAL")

# # load all the sample -----------------------------------------------------
# # define all the files available
# folder <- "../out_large/"
# file_id <- dir(folder) %>% 
#   str_subset("^brain_V") %>% 
#   str_remove_all(pattern = ".rds")
# 
# # build a full matrix for the expression information on TSPO
# df_exp <- lapply(file_id,function(x){
#   
#   # keep track of the progress of the reading
#   print(x)
#   
#   # read in the data
#   brain <- readRDS(paste0("../out_large/",x,".rds"))
#   
#   # extract the expression of the feature
#   df_exp <- brain@assays$SCT@data["TSPO",] %>% 
#     data.frame(exp = .) %>% 
#     rownames_to_column("barcodes") %>% 
#     mutate(gene = "TSPO") %>% 
#     mutate(dataset = x)
#   
#   return(df_exp)
# 
# }) %>% 
#   bind_rows()

# dim(df_exp)

# # read in the deconvolution result from each slice
# folder2 <- "out/table/"
# file_id2 <- dir(folder2) %>% 
#   str_subset("SPOTLight") %>% 
#   str_remove_all(pattern = "_SPOTLight.tsv")

# # build a full matrix for the expression information on TSPO
# df_dec <- lapply(file_id2,function(x){
#   
#   # keep track of the progress of the reading
#   print(x)
#   
#   df_meta <- read_tsv(paste0("out/table/",x,"_SPOTLight.tsv")) %>% 
#     mutate(dataset = paste0("brain_",x))
#   
#   return(df_meta)
#   
# }) %>% 
#   bind_rows()
# 
# dim(df_dec)
# 
# # add the annotation of the sample type based on the slide
# LUT_sample <- read_csv("data/sample_classification.csv") %>% 
#   mutate(dataset = paste0("brain_",sample))
# 
# # join the two datasets
# df_tot <- left_join(df_dec,LUT_sample,by = "dataset") %>% 
#   left_join(df_exp,c("barcodes","dataset")) %>% 
#   pivot_longer(cols = AST:VAS,names_to = "cell_type",values_to = "prop") %>% 
#   # after chatting with sofia she decided to remove the BASAL sample
#   filter(type != "BASAL")

# plot total correlation --------------------------------------------------
df_tot %>% 
  ggplot(aes(x=prop,y=exp)) +
  geom_point(alpha=0.2) +
  facet_grid(dataset~cell_type)+geom_smooth()+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("out/image/sofia_TSPO/scatterplot_cor_TSPO.pdf",height = 20,width = 20)

# df_tot %>%
#   dplyr::filter(dataset == "brain_V01") %>% 
#   dplyr::select(prop,exp) %>%
#   mutate(exp_noise = exp+rnorm(n = n())/100) %>% 
#   ggplot(aes(x=prop,y=exp_noise)) +
#   stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   theme(
#     legend.position='none'
#   )

# plot correlation split --------------------------------------------------
mat0 <- df_tot %>% 
  group_by(cell_type,dataset) %>% 
  summarise(cor_TSPO = cor(prop,exp,method = "pearson")) %>%
  pivot_wider(names_from = dataset,values_from = cor_TSPO) %>% 
  column_to_rownames("cell_type")

# build the heatmap annotation
sample_ordered0 <- data.frame(sample = str_extract(colnames(mat0),pattern = c("V\\d+"))) %>% 
  left_join(LUT_sample)

table(sample_ordered0$type)

sample_ordered_conversion0 <- sample_ordered0 %>% pull(type)
column_ha0 <- HeatmapAnnotation(sample = sample_ordered_conversion0,
                               col = list(sample = c("NAWM"="cyan",
                                                   "BASAL"="gray",
                                                   "CORTEX"="green",
                                                   "CA"="red",
                                                   "A"="orange",
                                                   "CI"="yellow",
                                                   "RM"="blue")))

pdf("out/image/sofia_TSPO/heatmap_cor_TSPO_pearson.pdf",width = 6,height = 4)
Heatmap(mat0,name = "cor",
        top_annotation = column_ha0)  
dev.off()

# do the same but compound the slice in sample types
mat02 <- df_tot %>% 
  group_by(cell_type,type) %>% 
  summarise(cor_TSPO = cor(prop,exp,method = "pearson")) %>%
  pivot_wider(names_from = type,values_from = cor_TSPO) %>% 
  column_to_rownames("cell_type")

# build the heatmap annotation
pdf("out/image/sofia_TSPO/heatmap_cor_TSPO_pearson2.pdf",width = 4,height = 4)
Heatmap(mat02,name = "cor")  
dev.off()

# spearman
mat1 <- df_tot %>% 
  group_by(cell_type,dataset) %>% 
  summarise(cor_TSPO = cor(prop,exp,method = "spearman")) %>%
  pivot_wider(names_from = dataset,values_from = cor_TSPO) %>% 
  column_to_rownames("cell_type")

# build the annotation
sample_ordered1 <- data.frame(sample = str_extract(colnames(mat1),pattern = c("V\\d+"))) %>% 
  left_join(LUT_sample)

table(sample_ordered1$type)

sample_ordered_conversion1 <- sample_ordered1 %>% pull(type)
column_ha1 <- HeatmapAnnotation(sample = sample_ordered_conversion1,
                               col = list(sample = c("NAWM"="cyan",
                                                     "BASAL"="gray",
                                                     "CORTEX"="green",
                                                     "CA"="red",
                                                     "A"="orange",
                                                     "CI"="yellow",
                                                     "RM"="blue")))


pdf("out/image/sofia_TSPO/heatmap_cor_TSPO_spearman.pdf",width = 6,height = 4)
Heatmap(mat1,name = "cor",
        top_annotation = column_ha1)
dev.off()

# do the same but compound the slice in sample types
mat12 <- df_tot %>% 
  group_by(cell_type,type) %>% 
  summarise(cor_TSPO = cor(prop,exp,method = "spearman")) %>%
  pivot_wider(names_from = type,values_from = cor_TSPO) %>% 
  column_to_rownames("cell_type")

# build the heatmap annotation
pdf("out/image/sofia_TSPO/heatmap_cor_TSPO_spearman2.pdf",width = 4,height = 4)
Heatmap(mat12,name = "cor")  
dev.off()

# plot total expression per type ------------------------------------------
# is there an overall higher expression of TSPO in the slide overall?
df_tot %>% 
  group_by(sample,type) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  group_by(type) %>% 
  mutate(avg = mean(avg_TSPO)) %>%
  ungroup() %>% 
  mutate(type = fct_reorder(type,avg)) %>% 
  # mutate(type = factor(type,levels = c("NAWM","CORTEX","A","CA","CI","RM","BASAL"))) %>% 
  ggplot(aes(x=type,y=avg_TSPO))+
  # geom_boxplot(outlier.shape = NA,alpha = 0.5)+
  # stat_summary(fun.data = mean_cl_normal,col="red",geom = "errorbar",width=0.05,position=position_dodge(width = 0.75))+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+
  # stat_summary(fun.data = mean_cl_normal,col="red",geom = "point",shape=4,size=2,position=position_dodge(width = 0.75))+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  geom_text_repel(aes(label = sample))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/sofia_TSPO/scatter_avg_exp_TSPO_sample.pdf",width = 5,height = 5)

# also without the labels
df_tot %>% 
  group_by(sample,type) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  group_by(type) %>% 
  mutate(avg = mean(avg_TSPO)) %>%
  ungroup() %>% 
  mutate(type = fct_reorder(type,avg)) %>% 
  # mutate(type = factor(type,levels = c("NAWM","CORTEX","A","CA","CI","RM","BASAL"))) %>% 
  ggplot(aes(x=type,y=avg_TSPO))+
  # geom_boxplot(outlier.shape = NA,alpha = 0.5)+
  # stat_summary(fun.data = mean_cl_normal,col="red",geom = "errorbar",width=0.05,position=position_dodge(width = 0.75))+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+
  # stat_summary(fun.data = mean_cl_normal,col="red",geom = "point",shape=4,size=2,position=position_dodge(width = 0.75))+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  # geom_text_repel(aes(label = sample))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/sofia_TSPO/scatter_avg_exp_TSPO_sample2.pdf",width = 5,height = 5)

# rank the prop of the annotated cell type per spot
# df_top_rank <- df_tot %>%
#   group_by(barcodes,sample) %>%
#   # mutate(rank = rank(dplyr::desc(prop),ties.method = "first")) %>%
#   # filter(rank == 1)
#   top_n(wt = prop,n = 1)
# 
# 
# df_top_rank %>%
#   # average expression per sample, type and cell_type
#   group_by(sample,type,cell_type) %>%
#   summarise(avg_TSPO = mean(exp)) %>%
#   ungroup() %>%
#   # average expression per type fo the sample
#   group_by(type) %>%
#   mutate(avg_type = mean(avg_TSPO)) %>%
#   ungroup() %>%
#   mutate(type = fct_reorder(type,avg_type)) %>%
#   # average expression per cell_type
#   group_by(cell_type) %>%
#   mutate(avg_cell_type = mean(avg_TSPO)) %>%
#   ungroup() %>%
#   mutate(cell_type = fct_reorder(cell_type,dplyr::desc(avg_cell_type))) %>%
#   # ggplot(aes(x=type,y=avg_TSPO,col=cell_type)) +
#   ggplot(aes(x=type,y=avg_TSPO)) +
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
#   facet_wrap(~cell_type)+
#   theme(axis.text.x = element_text(hjust = 1,angle = 90),
#         strip.background = element_blank())+
#   stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)
# ggsave("out/image/sofia_TSPO/scatter_avg_exp_TSPO_sampleCellType_TopAssigned.pdf",width = 9,height = 9)

# plot expression of TSPO top assignment ----------------------------------
df_tot %>% 
  # average expression per sample, type and cell_type
  group_by(sample,type,class_top) %>%
  summarise(avg_TSPO = mean(exp)) %>%
  ungroup() %>%
  # average expression per type fo the sample
  group_by(type) %>%
  mutate(avg_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(type = fct_reorder(type,avg_type)) %>%
  # average expression per cell_type
  group_by(class_top) %>%
  mutate(avg_cell_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(class_top = fct_reorder(class_top,dplyr::desc(avg_cell_type))) %>%
  # ggplot(aes(x=type,y=avg_TSPO,col=cell_type)) +
  ggplot(aes(x=type,y=avg_TSPO)) +
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  facet_wrap(~class_top)+
  theme(axis.text.x = element_text(hjust = 1,angle = 90),
        strip.background = element_blank())+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)
ggsave("out/image/sofia_TSPO/scatter_avg_exp_TSPO_sampleCellType_TopAssigned.pdf",width = 9,height = 9)

# # render as an heatmap
# mat_tot <- df_tot %>%
#   # average expression per sample, type and cell_type
#   group_by(sample,type,class_top) %>% 
#   summarise(avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   mutate(sample_type = paste0(sample,"_",type)) %>% 
#   dplyr::select(-c(sample,type)) %>% 
#   pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
#   column_to_rownames("class_top")
# 
# # build the heatmap annotation
# sample_ordered_tot <- data.frame(sample = str_extract(colnames(mat_tot),pattern = c("V\\d+"))) %>% 
#   left_join(LUT_sample)
# 
# table(sample_ordered_tot$type)
# 
# sample_ordered_tot_conversion <- sample_ordered_tot %>% pull(type)
# column_ha_tot <- HeatmapAnnotation(sample = sample_ordered_tot_conversion,
#                                    col = list(sample = c("NAWM"="cyan",
#                                                          "BASAL"="gray",
#                                                          "CORTEX"="green",
#                                                          "CA"="red",
#                                                          "A"="orange",
#                                                          "CI"="yellow",
#                                                          "RM"="blue")))
# 
# pdf("out/image/sofia_TSPO/heatmap_avg_exp_TSPO_sampleCellType_TopAssigned.pdf",width = 6,height = 4)
# Heatmap(mat_tot,name = "avg exp TSPO",
#         top_annotation = column_ha_tot,
#         col = viridis::viridis(20,option = "inferno"))  
# dev.off()

# render as an heatmap
mat_tot0 <- df_tot %>%
  # average expression per sample, type and cell_type
  group_by(sample,type,class_top) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  mutate(sample_type = paste0(sample,"_",type)) %>% 
  dplyr::select(-c(sample,type)) %>% 
  pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  column_to_rownames("class_top")

# build the heatmap annotation
sample_ordered_tot0 <- data.frame(sample = str_extract(colnames(mat_tot0),pattern = c("V\\d+"))) %>% 
  left_join(LUT_sample)

table(sample_ordered_tot0$type)

sample_ordered_tot_conversion0 <- sample_ordered_tot0 %>% pull(type)
column_ha_tot0 <- HeatmapAnnotation(sample = sample_ordered_tot_conversion0,
                                   col = list(sample = c("NAWM"="cyan",
                                                         "BASAL"="gray",
                                                         "CORTEX"="green",
                                                         "CA"="red",
                                                         "A"="orange",
                                                         "CI"="yellow",
                                                         "RM"="blue")))

pdf("out/image/sofia_TSPO/heatmap_avg_exp_TSPO_sampleCellType_TopAssigned.pdf",width = 6,height = 4)
Heatmap(mat_tot0,name = "avg exp TSPO",
        top_annotation = column_ha_tot0,
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# do the same but compound the slice in sample types
mat_tot02 <- df_tot %>%
  # average expression per sample, type and cell_type
  group_by(type,class_top) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  # mutate(sample_type = paste0(sample,"_",type)) %>% 
  # dplyr::select(-c(sample,type)) %>% 
  # pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  pivot_wider(names_from = type,values_from = avg_TSPO) %>% 
  column_to_rownames("class_top")

# build the heatmap annotation
pdf("out/image/sofia_TSPO/heatmap_avg_exp_TSPO_sampleCellType_TopAssigned2.pdf",width = 4,height = 4)
Heatmap(mat_tot02,name = "avg exp TSPO",
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# plot expression of TSPO top 50 assignment -------------------------------
df_tot %>% 
  # average expression per sample, type and cell_type
  group_by(sample,type,class_50) %>%
  summarise(avg_TSPO = mean(exp)) %>%
  ungroup() %>%
  # average expression per type fo the sample
  group_by(type) %>%
  mutate(avg_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(type = fct_reorder(type,avg_type)) %>%
  # average expression per cell_type
  group_by(class_50) %>%
  mutate(avg_cell_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(class_50 = fct_reorder(class_50,dplyr::desc(avg_cell_type))) %>%
  # ggplot(aes(x=type,y=avg_TSPO,col=cell_type)) +
  ggplot(aes(x=type,y=avg_TSPO)) +
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  facet_wrap(~class_50)+
  theme(axis.text.x = element_text(hjust = 1,angle = 90),
        strip.background = element_blank())+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)
ggsave("out/image/sofia_TSPO/scatter_avg_exp_TSPO_sampleCellType_Top50Assigned.pdf",width = 9,height = 9)

# render as an heatmap
mat_tot50 <- df_tot %>%
  # average expression per sample, type and cell_type
  group_by(sample,type,class_50) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  mutate(sample_type = paste0(sample,"_",type)) %>% 
  dplyr::select(-c(sample,type)) %>% 
  pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  column_to_rownames("class_50")

# build the heatmap annotation
sample_ordered_tot50 <- data.frame(sample = str_extract(colnames(mat_tot50),pattern = c("V\\d+"))) %>% 
  left_join(LUT_sample)

table(sample_ordered_tot50$type)

sample_ordered_tot_conversion50 <- sample_ordered_tot50 %>% pull(type)
column_ha_tot50 <- HeatmapAnnotation(sample = sample_ordered_tot_conversion50,
                                   col = list(sample = c("NAWM"="cyan",
                                                         "BASAL"="gray",
                                                         "CORTEX"="green",
                                                         "CA"="red",
                                                         "A"="orange",
                                                         "CI"="yellow",
                                                         "RM"="blue")))

pdf("out/image/sofia_TSPO/heatmap_avg_exp_TSPO_sampleCellType_Top50Assigned.pdf",width = 6,height = 4)
Heatmap(mat_tot50,name = "avg exp TSPO",
        top_annotation = column_ha_tot50,
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# do the same but compound the slice in sample types
mat_tot502 <- df_tot %>%
  # average expression per sample, type and cell_type
  group_by(type,class_50) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  # mutate(sample_type = paste0(sample,"_",type)) %>% 
  # dplyr::select(-c(sample,type)) %>% 
  # pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  pivot_wider(names_from = type,values_from = avg_TSPO) %>% 
  column_to_rownames("class_50")

# build the heatmap annotation
pdf("out/image/sofia_TSPO/heatmap_avg_exp_TSPO_sampleCellType_Top50Assigned2.pdf",width = 4,height = 4)
Heatmap(mat_tot502,name = "avg exp TSPO",
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# plot expression of TSPO top 70 assignment -------------------------------
df_tot %>% 
  # average expression per sample, type and cell_type
  group_by(sample,type,class_70) %>%
  summarise(avg_TSPO = mean(exp)) %>%
  ungroup() %>%
  # average expression per type fo the sample
  group_by(type) %>%
  mutate(avg_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(type = fct_reorder(type,avg_type)) %>%
  # average expression per cell_type
  group_by(class_70) %>%
  mutate(avg_cell_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(class_70 = fct_reorder(class_70,dplyr::desc(avg_cell_type))) %>%
  # ggplot(aes(x=type,y=avg_TSPO,col=cell_type)) +
  ggplot(aes(x=type,y=avg_TSPO)) +
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  facet_wrap(~class_70)+
  theme(axis.text.x = element_text(hjust = 1,angle = 90),
        strip.background = element_blank())+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)
ggsave("out/image/sofia_TSPO/scatter_avg_exp_TSPO_sampleCellType_Top70Assigned.pdf",width = 9,height = 9)

# render as an heatmap
mat_tot70 <- df_tot %>%
  # average expression per sample, type and cell_type
  group_by(sample,type,class_70) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  mutate(sample_type = paste0(sample,"_",type)) %>% 
  dplyr::select(-c(sample,type)) %>% 
  pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  column_to_rownames("class_70")

# build the heatmap annotation
sample_ordered_tot70 <- data.frame(sample = str_extract(colnames(mat_tot70),pattern = c("V\\d+"))) %>% 
  left_join(LUT_sample)

table(sample_ordered_tot70$type)

sample_ordered_tot_conversion70 <- sample_ordered_tot70 %>% pull(type)
column_ha_tot70 <- HeatmapAnnotation(sample = sample_ordered_tot_conversion70,
                                     col = list(sample = c("NAWM"="cyan",
                                                           "BASAL"="gray",
                                                           "CORTEX"="green",
                                                           "CA"="red",
                                                           "A"="orange",
                                                           "CI"="yellow",
                                                           "RM"="blue")))

pdf("out/image/sofia_TSPO/heatmap_avg_exp_TSPO_sampleCellType_Top70Assigned.pdf",width = 6,height = 4)
Heatmap(mat_tot70,name = "avg exp TSPO",
        top_annotation = column_ha_tot70,
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# do the same but compound the slice in sample types
mat_tot702 <- df_tot %>%
  # average expression per sample, type and cell_type
  group_by(type,class_70) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  # mutate(sample_type = paste0(sample,"_",type)) %>% 
  # dplyr::select(-c(sample,type)) %>% 
  # pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  pivot_wider(names_from = type,values_from = avg_TSPO) %>% 
  column_to_rownames("class_70")

# build the heatmap annotation
pdf("out/image/sofia_TSPO/heatmap_avg_exp_TSPO_sampleCellType_Top70Assigned2.pdf",width = 4,height = 4)
Heatmap(mat_tot702,name = "avg exp TSPO",
        col = viridis::viridis(20,option = "inferno"))  
dev.off()






























# # assign only the spot that have more than 50 % of the specific cell type
# df_top_rank50 <- df_tot %>% 
#   group_by(barcodes,sample) %>% 
#   mutate(rank = rank(desc(prop),ties.method = "first")) %>% 
#   filter(rank == 1) %>% 
#   filter(prop > 0.5)
# 
# df_top_rank50 %>%
#   # average expression per sample, type and cell_type
#   group_by(sample,type,cell_type) %>% 
#   summarise(avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   # average expression per type fo the sample
#   group_by(type) %>% 
#   mutate(avg_type = mean(avg_TSPO)) %>%
#   ungroup() %>% 
#   mutate(type = fct_reorder(type,avg_type)) %>% 
#   # average expression per cell_type
#   group_by(cell_type) %>% 
#   mutate(avg_cell_type = mean(avg_TSPO)) %>%
#   ungroup() %>% 
#   mutate(cell_type = fct_reorder(cell_type,desc(avg_cell_type))) %>%
#   # ggplot(aes(x=type,y=avg_TSPO,col=cell_type)) +
#   ggplot(aes(x=type,y=avg_TSPO)) +
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
#   facet_wrap(~cell_type)+
#   theme(axis.text.x = element_text(hjust = 1,angle = 90),
#         strip.background = element_blank())+
#   stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)
# ggsave("out/image/scatter_avg_exp_TSPO_sampleCellType_Above50.pdf",width = 9,height = 9)
# 
# # render as an heatmap
# mat_50 <- df_top_rank50 %>%
#   # average expression per sample, type and cell_type
#   group_by(sample,type,cell_type) %>% 
#   summarise(avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   mutate(sample_type = paste0(sample,"_",type)) %>% 
#   dplyr::select(-c(sample,type)) %>% 
#   pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
#   column_to_rownames("cell_type")
# 
# # build the heatmap annotation
# sample_ordered_50 <- data.frame(sample = str_extract(colnames(mat_50),pattern = c("V\\d+"))) %>% 
#   left_join(LUT_sample)
# 
# table(sample_ordered_50$type)
# 
# sample_ordered_50_conversion <- sample_ordered_50 %>% pull(type)
# column_ha_50 <- HeatmapAnnotation(sample = sample_ordered_50_conversion,
#                                   col = list(sample = c("NAWM"="cyan",
#                                                         "BASAL"="gray",
#                                                         "CORTEX"="green",
#                                                         "CA"="red",
#                                                         "A"="orange",
#                                                         "CI"="yellow",
#                                                         "RM"="blue")))
# 
# pdf("out/image/heatmap_avg_TSPO_50.pdf",width = 6,height = 4)
# Heatmap(mat_50,name = "avg exp TSPO",
#         top_annotation = column_ha_50,
#         col = viridis::viridis(20,option = "inferno"))  
# dev.off()
# 
# # assign only the spot that have more than 70 % of the specific cell type
# df_top_rank70 <- df_tot %>% 
#   group_by(barcodes,sample) %>% 
#   mutate(rank = rank(desc(prop),ties.method = "first")) %>% 
#   filter(rank == 1) %>% 
#   filter(prop > 0.7)
# 
# df_top_rank70 %>%
#   # average expression per sample, type and cell_type
#   group_by(sample,type,cell_type) %>% 
#   summarise(avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   # average expression per type fo the sample
#   group_by(type) %>% 
#   mutate(avg_type = mean(avg_TSPO)) %>%
#   ungroup() %>% 
#   mutate(type = fct_reorder(type,avg_type)) %>% 
#   # average expression per cell_type
#   group_by(cell_type) %>% 
#   mutate(avg_cell_type = mean(avg_TSPO)) %>%
#   ungroup() %>% 
#   mutate(cell_type = fct_reorder(cell_type,desc(avg_cell_type))) %>%
#   # ggplot(aes(x=type,y=avg_TSPO,col=cell_type)) +
#   ggplot(aes(x=type,y=avg_TSPO)) +
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
#   facet_wrap(~cell_type)+
#   theme(axis.text.x = element_text(hjust = 1,angle = 90),
#         strip.background = element_blank())+
#   stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)
# ggsave("out/image/scatter_avg_exp_TSPO_sampleCellType_Above70.pdf",width = 9,height = 9)
# 
# # render as an heatmap
# mat_70 <- df_top_rank70 %>%
#   # average expression per sample, type and cell_type
#   group_by(sample,type,cell_type) %>% 
#   summarise(avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   mutate(sample_type = paste0(sample,"_",type)) %>% 
#   dplyr::select(-c(sample,type)) %>% 
#   pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
#   column_to_rownames("cell_type")
# 
# # build the heatmap annotation
# sample_ordered_70 <- data.frame(sample = str_extract(colnames(mat_70),pattern = c("V\\d+"))) %>% 
#   left_join(LUT_sample)
# 
# table(sample_ordered_70$type)
# 
# sample_ordered_70_conversion <- sample_ordered_70 %>% pull(type)
# column_ha_70 <- HeatmapAnnotation(sample = sample_ordered_70_conversion,
#                                   col = list(sample = c("NAWM"="cyan",
#                                                         "BASAL"="gray",
#                                                         "CORTEX"="green",
#                                                         "CA"="red",
#                                                         "A"="orange",
#                                                         "CI"="yellow",
#                                                         "RM"="blue")))
# 
# pdf("out/image/heatmap_avg_TSPO_70.pdf",width = 6,height = 4)
# Heatmap(mat_70,name = "avg exp TSPO",
#         top_annotation = column_ha_70,
#         col = viridis::viridis(20,option = "inferno"))  
# dev.off()
