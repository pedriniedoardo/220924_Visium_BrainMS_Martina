# # libraries ---------------------------------------------------------------
# library(Seurat)
# library(SeuratData)
# library(ggplot2)
# library(patchwork)
# library(dplyr)
# library(tidyverse)
# library(hdf5r)
# library(limma)
# library(future)
# library(ComplexHeatmap)
# 
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
# 
# # read in the deconvolution result from each slice
# folder2 <- "out/table/"
# file_id2 <- dir(folder2) %>% 
#   str_subset("SPOTLight") %>% 
#   str_remove_all(pattern = "_SPOTLight.tsv")
# 
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
# # laod the segmentation from Sofia
# folder3 <- "data/segmentation/sofia/all/"
# file_id3 <- dir(folder3) %>% 
#   str_subset("sofia_segmentation") %>% 
#   str_remove_all(pattern = "_sofia_segmentation.csv")
# 
# # build a full matrix for the expression information on TSPO
# df_sofi <- lapply(file_id3,function(x){
#   
#   # keep track of the progress of the reading
#   print(x)
#   
#   df_meta <- read_csv(paste0("data/segmentation/sofia/all/",x,"_sofia_segmentation.csv")) %>% 
#     mutate(dataset = paste0("brain_",x)) %>% 
#     # some point do not have an annotation fiter them out
#     filter(!is.na(sofia_segmentation))
#   
#   return(df_meta)
#   
# }) %>% 
#   bind_rows()
# 
# dim(df_sofi)
# 
# # add the annotation of the sample type based on the slide
# LUT_sample <- read_csv("data/sample_classification.csv") %>% 
#   mutate(dataset = paste0("brain_",sample))
# 
# # join the two datasets
# # some slide are available in sofia's annotation but not in mine therefore use inner join
# df_tot <- left_join(df_dec,LUT_sample,by = "dataset") %>% 
#   left_join(df_exp,c("barcodes","dataset")) %>%
#   inner_join(df_sofi,c("barcodes","dataset")) %>% 
#   pivot_longer(cols = AST:VAS,names_to = "cell_type",values_to = "prop")
# 
# # fix the meta
# df_tot_fix <- df_tot %>% 
#   mutate(sofia_segmentation = case_when(sofia_segmentation == "vessels"~"vessel",
#                                         sofia_segmentation == "WM_1layer"~"WM_01_layer",
#                                         sofia_segmentation == "WM_2layer"~"WM_02_layer",
#                                         sofia_segmentation == "WM_3layer"~"WM_03_layer",
#                                         sofia_segmentation == "WM_4lyer"~"WM_04_layer",
#                                         sofia_segmentation == "WM_4layer"~"WM_04_layer",
#                                         sofia_segmentation == "WM_5layer"~"WM_05_layer",
#                                         sofia_segmentation == "WM_6layer"~"WM_06_layer",
#                                         sofia_segmentation == "WM_7layer"~"WM_07_layer",
#                                         sofia_segmentation == "WM_8layer"~"WM_08_layer",
#                                         sofia_segmentation == "WM_9layer"~"WM_09_layer",
#                                         sofia_segmentation == "WM_10layer"~"WM_10_layer",
#                                         T~sofia_segmentation)) %>% 
#   mutate(sofia_segmentation_short = case_when(str_detect(sofia_segmentation,pattern = "WM_")~"WM",
#                                               T~sofia_segmentation))
# 
# # check the metadata
# table(df_tot_fix$sofia_segmentation)
# table(df_tot_fix$sofia_segmentation_short)
# sum(table(df_tot_fix$sofia_segmentation_short))
# 
# # plotting ----------------------------------------------------------------
# # corr plot using the short annotation
# list_plot <- df_tot_fix %>% 
#   # split by tissue segmentaion
#   split(f = .$sofia_segmentation_short)
# 
# list_plot_corr <- pmap(list(names(list_plot),list_plot),function(name,x){
#     x %>%
#     ggplot(aes(x=prop,y=exp)) +
#     geom_point(alpha=0.2) +
#     facet_grid(dataset~cell_type)+geom_smooth(method = "lm")+
#     theme_bw()+
#     theme(strip.background = element_blank())+ggtitle(name)
#   })
# 
# # save the plots fron the list in different pages fo the same file
# pdf("out/image/scatterplot_cor_TSPO_sofiaSegmentationShort.pdf",width = 20,height = 20)
# list_plot_corr
# dev.off()
# 
# # corr plot using the long annotation  
# list_plot_long <- df_tot_fix %>% 
#   # split by tissue segmentaion
#   split(f = .$sofia_segmentation)
# 
# list_plot_corr_long <- pmap(list(names(list_plot_long),list_plot_long),function(name,x){
#   x %>%
#     ggplot(aes(x=prop,y=exp)) +
#     geom_point(alpha=0.2) +
#     facet_grid(dataset~cell_type)+geom_smooth(method = "lm")+
#     theme_bw()+
#     theme(strip.background = element_blank())+ggtitle(name)
# })
# 
# # save the plots fron the list in different pages fo the same file
# pdf("out/image/scatterplot_cor_TSPO_sofiaSegmentationLong.pdf",width = 20,height = 20)
# list_plot_corr_long
# dev.off()
# 
# # -------------------------------------------------------------------------
# mat_short <- df_tot_fix %>% 
#   group_by(dataset,cell_type,sofia_segmentation_short) %>%  
#   summarise(cor_TSPO = cor(prop,exp,method = "pearson")) %>%
#   mutate(dataset_segmentation = paste0(dataset,"_",sofia_segmentation_short)) %>%
#   ungroup() %>% 
#   dplyr::select(-c(dataset,sofia_segmentation_short)) %>%
#   pivot_wider(names_from = dataset_segmentation,values_from = cor_TSPO) %>% 
#   column_to_rownames("cell_type")
# 
# # mat_short %>%
# #   rownames_to_column() %>%
# #   write_tsv("out/table/test_tidyverse.tsv")
# 
# # filter out all the columns that are NA
# id <- apply(mat_short,MARGIN = 2,is.na) %>% 
#   colSums()==0
# 
# mat_short_filter <- mat_short[,id]
# 
# # build the annotation on the heatmap
# # read in the sample classification
# # LUT_sample <- read_csv("data/sample_classification.csv")
# sample_ordered <- data.frame(sample = str_extract(colnames(mat_short_filter),pattern = c("V\\d+"))) %>% 
#   left_join(LUT_sample)
# 
# table(sample_ordered$type)
# 
# sample_ordered_conversion <- sample_ordered %>% pull(type)
# area_ordered <- str_extract(colnames(mat_short_filter),pattern = c("WM|core|cortex|edge|vessel|active lesion"))
# 
# column_ha <- HeatmapAnnotation(sample = sample_ordered_conversion,
#                                area = area_ordered,
#                                col = list(area = c("WM"="beige",
#                                                     "core"="red",
#                                                     "cortex"="gray",
#                                                     "edge"="orange",
#                                                     "vessel"="green",
#                                                     "active lesion"="yellow"),
#                                           sample = c("NAWM"="cyan",
#                                                      "BASAL"="gray",
#                                                      "CORTEX"="green",
#                                                      "CA"="red",
#                                                      "A"="orange",
#                                                      "CI"="yellow",
#                                                      "RM"="blue")))
# 
# 
# pdf("out/image/heatmap_cor_TSPO_pearson_SofiaShort.pdf",width = 10,height = 4)
# Heatmap(mat_short_filter,name = "cor",
#         top_annotation = column_ha)  
# dev.off()
# 
# # -------------------------------------------------------------------------
# mat_short2 <- df_tot_fix %>% 
#   group_by(dataset,cell_type,sofia_segmentation_short) %>%  
#   summarise(cor_TSPO = cor(prop,exp,method = "spearman")) %>%
#   mutate(dataset_segmentation = paste0(dataset,"_",sofia_segmentation_short)) %>%
#   ungroup() %>% 
#   dplyr::select(-c(dataset,sofia_segmentation_short)) %>%
#   pivot_wider(names_from = dataset_segmentation,values_from = cor_TSPO) %>% 
#   column_to_rownames("cell_type")
# 
# # mat_short %>%
# #   rownames_to_column() %>%
# #   write_tsv("out/table/test_tidyverse.tsv")
# 
# # filter out all the columns that are NA
# id2 <- apply(mat_short2,MARGIN = 2,is.na) %>% 
#   colSums()==0
# 
# mat_short_filter2 <- mat_short2[,id2]
# 
# # build the annotation on the heatmap
# # read in the sample classification
# # LUT_sample <- read_csv("data/sample_classification.csv")
# sample_ordered2 <- data.frame(sample = str_extract(colnames(mat_short_filter2),pattern = c("V\\d+"))) %>% 
#   left_join(LUT_sample)
# 
# table(sample_ordered2$type)
# 
# sample_ordered_conversion2 <- sample_ordered2 %>% pull(type)
# area_ordered2 <- str_extract(colnames(mat_short_filter2),pattern = c("WM|core|cortex|edge|vessel|active lesion"))
# 
# column_ha2 <- HeatmapAnnotation(sample = sample_ordered_conversion2,
#                                area = area_ordered2,
#                                col = list(area = c("WM"="beige",
#                                                    "core"="red",
#                                                    "cortex"="gray",
#                                                    "edge"="orange",
#                                                    "vessel"="green",
#                                                    "active lesion"="yellow"),
#                                           sample = c("NAWM"="cyan",
#                                                      "BASAL"="gray",
#                                                      "CORTEX"="green",
#                                                      "CA"="red",
#                                                      "A"="orange",
#                                                      "CI"="yellow",
#                                                      "RM"="blue")))
# 
# 
# pdf("out/image/heatmap_cor_TSPO_spearman_SofiaShort.pdf",width = 10,height = 4)
# Heatmap(mat_short_filter2,name = "cor",
#         top_annotation = column_ha2)  
# dev.off()
# 
# # -------------------------------------------------------------------------
# # average expression per sofia's region
# df_plot <- df_tot_fix %>% 
#   group_by(sample,type,sofia_segmentation_short) %>% 
#   summarise(avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   # average sofia_segmentation_short
#   group_by(sofia_segmentation_short) %>% 
#   mutate(avg_sofia = mean(avg_TSPO)) %>%
#   ungroup() %>% 
#   mutate(sofia_segmentation_short = fct_reorder(sofia_segmentation_short,avg_sofia)) %>%
#   # average the type
#   group_by(type) %>% 
#   mutate(avg_type = mean(avg_TSPO)) %>%
#   ungroup() %>% 
#   mutate(type = fct_reorder(type,desc(avg_type)))
# 
# df_plot %>%
#   ggplot(aes(x=sofia_segmentation_short,y=avg_TSPO))+
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
#   stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+
#   # geom_text_repel(aes(label = sample))+
#   theme(axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("out/image/scatter_avg_TSPO_SofiaShort_tot.pdf",width = 5,height = 5)
# 
# # split by type
# df_plot %>%
#   ggplot(aes(x=sofia_segmentation_short,y=avg_TSPO))+
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
#   stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+
#   # geom_text_repel(aes(label = sample))+
#   theme(axis.text.x = element_text(hjust = 1,angle = 90),strip.background = element_blank())+facet_wrap(~type)
# ggsave("out/image/scatter_avg_TSPO_SofiaShort_tot_split.pdf",width = 9,height = 6)
# 
# df_plot %>%
#   filter(type == "A")
# 
# df_tot_fix %>% 
#   filter(type == "A") %>% 
#   group_by(sample,type,sofia_segmentation_short,cell_type) %>% 
#   summarise(n = n(),
#             median = median(exp),
#     avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   print(n = 30)
# 
# # is there an overall higher expression of TSPO in the slide overall?
# mat_avg <- df_tot_fix %>% 
#   group_by(sample,type,sofia_segmentation_short) %>% 
#   summarise(avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   mutate(sample_type = paste0(sample,"_",type)) %>%
#   ungroup() %>% 
#   dplyr::select(-c(sample,type)) %>%
#   pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
#   column_to_rownames("sofia_segmentation_short")
# 
# # mat_short %>%
# #   rownames_to_column() %>%
# #   write_tsv("out/table/test_tidyverse.tsv")
# 
# # filter out all the columns that are NA
# id_avg <- apply(mat_avg,MARGIN = 2,is.na) %>% 
#   colSums()>0
# 
# mat_avg_filter <- mat_avg[1:5,id_avg]
# 
# # build the annotation on the heatmap
# # read in the sample classification
# # LUT_sample <- read_csv("data/sample_classification.csv")
# sample_ordered_avg <- data.frame(sample = str_extract(colnames(mat_avg_filter),pattern = c("V\\d+"))) %>% 
#   left_join(LUT_sample)
# 
# table(sample_ordered_avg$type)
# 
# sample_ordered_conversion_avg <- sample_ordered_avg %>% pull(type)
# 
# column_ha_avg <- HeatmapAnnotation(sample = sample_ordered_conversion_avg,
#                                    col = list(sample = c("NAWM"="cyan",
#                                                          "BASAL"="gray",
#                                                          "CORTEX"="green",
#                                                          "CA"="red",
#                                                          "A"="orange","CI"="yellow",
#                                                          "RM"="blue")))
# 
# 
# pdf("out/image/heatmap_avg_TSPO_SofiaShort_tot.pdf",width = 6,height = 4)
# Heatmap(mat_avg_filter,name = "avg exp TSPO",
#         top_annotation = column_ha_avg,
#         col = viridis::viridis(20,option = "inferno"))  
# dev.off()
# 
# # -------------------------------------------------------------------------
# # rank the prop of the annotated cell type per spot
# df_top_rank <- df_tot_fix %>% 
#   group_by(barcodes,sample) %>% 
#   mutate(rank = rank(desc(prop),ties.method = "first")) %>% 
#   filter(rank == 1)
# 
# # render as an heatmap
# mat_tot_list <- df_top_rank %>%
#   # average expression per sample, type and cell_type
#   group_by(sample,type,cell_type,sofia_segmentation_short) %>% 
#   summarise(avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   split(f = .$sofia_segmentation_short) %>% 
#   lapply(function(x){
#     x %>%
#       mutate(sample_type = paste0(sample,"_",type)) %>% 
#       dplyr::select(-c(sample,type,sofia_segmentation_short)) %>% 
#       # pivot_wider(names_from = sample_type,values_from = avg_TSPO,values_fill = 0) %>% 
#       pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
#       column_to_rownames("cell_type")
#   })
# 
# list_ht <- pmap(list(mat_tot_list[-1],names(mat_tot_list[-1])),function(x,name){
#   # build the heatmap annotation
#   sample_ordered_tot <- data.frame(sample = str_extract(colnames(x),pattern = c("V\\d+"))) %>% 
#     left_join(LUT_sample)
#   
#   table(sample_ordered_tot$type)
#   
#   sample_ordered_tot_conversion <- sample_ordered_tot %>% pull(type)
#   column_ha_tot <- HeatmapAnnotation(sample = sample_ordered_tot_conversion,show_annotation_name = F,
#                                      col = list(sample = c("NAWM"="cyan",
#                                                            "BASAL"="gray",
#                                                            "CORTEX"="green",
#                                                            "CA"="red",
#                                                            "A"="orange",
#                                                            "CI"="yellow",
#                                                            "RM"="blue")))
#   
#   ht <- Heatmap(x,name = "avg exp TSPO",
#           top_annotation = column_ha_tot,
#           column_title = name,
#           col = viridis::viridis(20,option = "inferno"))  
#   
#   return(ht)
# })
# 
# pdf("out/image/heatmap_avg_TSPO_SofiaShort_CellType_tot.pdf",width = 15,height = 4)
# list_ht$core+list_ht$cortex+list_ht$edge+list_ht$vessel+list_ht$WM
# dev.off()
# 
# # -------------------------------------------------------------------------
# # rank the prop of the annotated cell type per spot
# df_top_rank2 <- df_tot_fix %>% 
#   group_by(barcodes,sample) %>% 
#   mutate(rank = rank(desc(prop),ties.method = "first")) %>% 
#   filter(rank == 1) %>% 
#   filter(prop > 0.5)
# 
# # render as an heatmap
# mat_tot_list2 <- df_top_rank2 %>%
#   # average expression per sample, type and cell_type
#   group_by(sample,type,cell_type,sofia_segmentation_short) %>% 
#   summarise(avg_TSPO = mean(exp)) %>% 
#   ungroup() %>% 
#   split(f = .$sofia_segmentation_short) %>% 
#   lapply(function(x){
#     x %>%
#       mutate(sample_type = paste0(sample,"_",type)) %>% 
#       dplyr::select(-c(sample,type,sofia_segmentation_short)) %>% 
#       # pivot_wider(names_from = sample_type,values_from = avg_TSPO,values_fill = 0) %>% 
#       pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
#       column_to_rownames("cell_type")
#   })
# 
# list_ht2 <- pmap(list(mat_tot_list2[-1],names(mat_tot_list2[-1])),function(x,name){
#   # build the heatmap annotation
#   sample_ordered_tot <- data.frame(sample = str_extract(colnames(x),pattern = c("V\\d+"))) %>% 
#     left_join(LUT_sample)
#   
#   table(sample_ordered_tot$type)
#   
#   sample_ordered_tot_conversion <- sample_ordered_tot %>% pull(type)
#   column_ha_tot <- HeatmapAnnotation(sample = sample_ordered_tot_conversion,show_annotation_name = F,
#                                      col = list(sample = c("NAWM"="cyan",
#                                                            "BASAL"="gray",
#                                                            "CORTEX"="green",
#                                                            "CA"="red",
#                                                            "A"="orange",
#                                                            "CI"="yellow",
#                                                            "RM"="blue")))
#   
#   ht <- Heatmap(x,name = "avg exp TSPO",
#                 top_annotation = column_ha_tot,
#                 column_title = name,
#                 col = viridis::viridis(20,option = "inferno"))  
#   
#   return(ht)
# })
# 
# pdf("out/image/heatmap_avg_TSPO_SofiaShort_CellType_50.pdf",width = 15,height = 4)
# list_ht2$core+list_ht2$cortex+list_ht2$edge+list_ht2$vessel+list_ht2$WM
# dev.off()
