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

# laod the segmentation from Sofia
folder3 <- "data/segmentation/sofia/all/draft_02/"
file_id3 <- dir(folder3) %>% 
  str_subset("sofia_segmentation") %>% 
  str_remove_all(pattern = "_sofia_segmentation_no_vessel.csv")

# build a full matrix for the expression information on TSPO
df_sofi <- lapply(file_id3,function(x){
  
  # keep track of the progress of the reading
  print(x)
  
  df_meta <- read_csv(paste0("data/segmentation/sofia/all/draft_02/",x,"_sofia_segmentation_no_vessel.csv")) %>% 
    mutate(dataset = paste0("brain_",x)) %>% 
    # some point do not have an annotation fiter them out
    filter(!is.na(sofia_segmentation))
  
  return(df_meta)
  
}) %>% 
  bind_rows()

dim(df_sofi)

# add the annotation of the sample type based on the slide
LUT_sample <- read_csv("data/sample_classification.csv") %>% 
  mutate(dataset = paste0("brain_",sample))

# join the two datasets
df_tot <- left_join(df_dec,LUT_sample,by = "dataset") %>% 
  left_join(df_exp,c("barcodes","dataset")) %>%
  inner_join(df_sofi,c("barcodes","dataset")) %>% 
  pivot_longer(cols = AST:VAS,names_to = "cell_type",values_to = "prop") %>% 
  # after chatting with sofia she decided to remove the BASAL sample
  filter(type != "BASAL")

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
sum(table(df_tot_fix$sofia_segmentation_short))

# plotting ----------------------------------------------------------------
# corr plot using the short annotation
list_plot <- df_tot_fix %>% 
  # split by tissue segmentaion
  split(f = .$sofia_segmentation_short)

list_plot_corr <- pmap(list(names(list_plot),list_plot),function(name,x){
  x %>%
    ggplot(aes(x=prop,y=exp)) +
    geom_point(alpha=0.2) +
    facet_grid(dataset~cell_type)+geom_smooth(method = "lm")+
    theme_bw()+
    theme(strip.background = element_blank())+ggtitle(name)
})

# save the plots fron the list in different pages fo the same file
pdf("out/image/sofia_TSPO_segmentation_draft02/scatterplot_cor_TSPO_sofiaSegmentationShort.pdf",width = 20,height = 20)
list_plot_corr
dev.off()

# corr plot using the long annotation  
list_plot_long <- df_tot_fix %>% 
  # split by tissue segmentaion
  split(f = .$sofia_segmentation)

list_plot_corr_long <- pmap(list(names(list_plot_long),list_plot_long),function(name,x){
  x %>%
    ggplot(aes(x=prop,y=exp)) +
    geom_point(alpha=0.2) +
    facet_grid(dataset~cell_type)+geom_smooth(method = "lm")+
    theme_bw()+
    theme(strip.background = element_blank())+ggtitle(name)
})

# save the plots fron the list in different pages fo the same file
pdf("out/image/sofia_TSPO_segmentation_draft02/scatterplot_cor_TSPO_sofiaSegmentationLong.pdf",width = 20,height = 20)
list_plot_corr_long
dev.off()

# -------------------------------------------------------------------------
mat_short <- df_tot_fix %>% 
  group_by(dataset,cell_type,sofia_segmentation_short) %>%  
  summarise(cor_TSPO = cor(prop,exp,method = "pearson")) %>%
  mutate(dataset_segmentation = paste0(dataset,"_",sofia_segmentation_short)) %>%
  ungroup() %>% 
  dplyr::select(-c(dataset,sofia_segmentation_short)) %>%
  pivot_wider(names_from = dataset_segmentation,values_from = cor_TSPO) %>% 
  column_to_rownames("cell_type")

# mat_short %>%
#   rownames_to_column() %>%
#   write_tsv("out/table/test_tidyverse.tsv")

# filter out all the columns that are NA
# id <- apply(mat_short,MARGIN = 2,is.na) %>% 
#   colSums()==0
id <- is.na(mat_short) %>%
  colSums() == 0

id2 <- is.na(mat_short[,id]) %>%
  rowSums() == 0

mat_short_filter <- mat_short[id2,id]

# build the annotation on the heatmap
# read in the sample classification
# LUT_sample <- read_csv("data/sample_classification.csv")
sample_ordered <- data.frame(sample = str_extract(colnames(mat_short_filter),pattern = c("V\\d+"))) %>% 
  left_join(LUT_sample)

table(sample_ordered$type)

sample_ordered_conversion <- sample_ordered %>% pull(type)
area_ordered <- str_extract(colnames(mat_short_filter),pattern = c("WM|core|cortex|edge|periplaque|active lesion"))

column_ha <- HeatmapAnnotation(sample = sample_ordered_conversion,
                               area = area_ordered,
                               col = list(area = c("WM"="beige",
                                                   "core"="red",
                                                   "cortex"="gray",
                                                   "edge"="orange",
                                                   "periplaque"="green",
                                                   "active lesion"="yellow"),
                                          sample = c("NAWM"="cyan",
                                                     "BASAL"="gray",
                                                     "CORTEX"="green",
                                                     "CA"="red",
                                                     "A"="orange",
                                                     "CI"="yellow",
                                                     "RM"="blue")))


pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_cor_TSPO_pearson_SofiaShort.pdf",width = 10,height = 4)
Heatmap(mat_short_filter,name = "cor",
        top_annotation = column_ha)  
dev.off()

# do the same but group by area
mat_short2 <- df_tot_fix %>% 
  group_by(cell_type,sofia_segmentation_short) %>%  
  summarise(cor_TSPO = cor(prop,exp,method = "pearson")) %>%
  # mutate(dataset_segmentation = paste0(dataset,"_",sofia_segmentation_short)) %>%
  # ungroup() %>% 
  # dplyr::select(-c(dataset,sofia_segmentation_short)) %>%
  # pivot_wider(names_from = dataset_segmentation,values_from = cor_TSPO) %>% 
  pivot_wider(names_from = sofia_segmentation_short,values_from = cor_TSPO) %>% 
  column_to_rownames("cell_type")

# mat_short %>%
#   rownames_to_column() %>%
#   write_tsv("out/table/test_tidyverse.tsv")

# filter out all the columns that are NA
# id <- apply(mat_short,MARGIN = 2,is.na) %>% 
#   colSums()==0
id2 <- is.na(mat_short2) %>%
  colSums() == 0

id22 <- is.na(mat_short2[,id2]) %>%
  rowSums() == 0

mat_short2_filter <- mat_short2[id22,id2]

pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_cor_TSPO_pearson_SofiaShort2.pdf",width = 4,height = 4)
Heatmap(mat_short2_filter,name = "cor")  
dev.off()

# -------------------------------------------------------------------------
mat_short_spearman <- df_tot_fix %>% 
  group_by(dataset,cell_type,sofia_segmentation_short) %>%  
  summarise(cor_TSPO = cor(prop,exp,method = "spearman")) %>%
  mutate(dataset_segmentation = paste0(dataset,"_",sofia_segmentation_short)) %>%
  ungroup() %>% 
  dplyr::select(-c(dataset,sofia_segmentation_short)) %>%
  pivot_wider(names_from = dataset_segmentation,values_from = cor_TSPO) %>% 
  column_to_rownames("cell_type")

# mat_short %>%
#   rownames_to_column() %>%
#   write_tsv("out/table/test_tidyverse.tsv")

# filter out all the columns that are NA
id_spearman <- is.na(mat_short_spearman) %>%
  colSums() == 0

id2_spearman <- is.na(mat_short_spearman[,id_spearman]) %>%
  rowSums() == 0

mat_short_spearman_filter <- mat_short_spearman[id2_spearman,id_spearman]

# build the annotation on the heatmap
# read in the sample classification
# LUT_sample <- read_csv("data/sample_classification.csv")
sample_ordered_spearman <- data.frame(sample = str_extract(colnames(mat_short_spearman_filter),pattern = c("V\\d+"))) %>% 
  left_join(LUT_sample)

table(sample_ordered_spearman$type)

sample_ordered_conversion_spearman <- sample_ordered_spearman %>% pull(type)
area_ordered_spearman <- str_extract(colnames(mat_short_spearman_filter),pattern = c("WM|core|cortex|edge|periplaque|active lesion"))

column_ha_spearman <- HeatmapAnnotation(sample = sample_ordered_conversion_spearman,
                                area = area_ordered_spearman,
                                col = list(area = c("WM"="beige",
                                                    "core"="red",
                                                    "cortex"="gray",
                                                    "edge"="orange",
                                                    "periplaque"="green",
                                                    "active lesion"="yellow"),
                                           sample = c("NAWM"="cyan",
                                                      "BASAL"="gray",
                                                      "CORTEX"="green",
                                                      "CA"="red",
                                                      "A"="orange",
                                                      "CI"="yellow",
                                                      "RM"="blue")))


pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_cor_TSPO_spearman_SofiaShort.pdf",width = 10,height = 4)
Heatmap(mat_short_spearman_filter,name = "cor",
        top_annotation = column_ha_spearman)  
dev.off()

# do the same by type
mat_short_spearman2 <- df_tot_fix %>% 
  group_by(cell_type,sofia_segmentation_short) %>%  
  summarise(cor_TSPO = cor(prop,exp,method = "spearman")) %>%
  # mutate(dataset_segmentation = paste0(dataset,"_",sofia_segmentation_short)) %>%
  # ungroup() %>% 
  # dplyr::select(-c(dataset,sofia_segmentation_short)) %>%
  # pivot_wider(names_from = dataset_segmentation,values_from = cor_TSPO) %>% 
  pivot_wider(names_from = sofia_segmentation_short,values_from = cor_TSPO) %>% 
  column_to_rownames("cell_type")

# mat_short %>%
#   rownames_to_column() %>%
#   write_tsv("out/table/test_tidyverse.tsv")

# filter out all the columns that are NA
id_spearman2 <- is.na(mat_short_spearman2) %>%
  colSums() == 0

id2_spearman2 <- is.na(mat_short_spearman2[,id_spearman2]) %>%
  rowSums() == 0

mat_short_spearman_filter2 <- mat_short_spearman2[id2_spearman2,id_spearman2]

pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_cor_TSPO_spearman_SofiaShort2.pdf",width = 4,height = 4)
Heatmap(mat_short_spearman_filter2,name = "cor")  
dev.off()


# average expression per sofia's region -----------------------------------
df_plot <- df_tot_fix %>% 
  group_by(sample,type,sofia_segmentation_short) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  # average sofia_segmentation_short
  group_by(sofia_segmentation_short) %>% 
  mutate(avg_sofia = mean(avg_TSPO)) %>%
  ungroup() %>% 
  mutate(sofia_segmentation_short = fct_reorder(sofia_segmentation_short,avg_sofia)) %>%
  # average the type
  group_by(type) %>% 
  mutate(avg_type = mean(avg_TSPO)) %>%
  ungroup() %>% 
  mutate(type = fct_reorder(type,dplyr::desc(avg_type)))

df_plot %>%
  ggplot(aes(x=sofia_segmentation_short,y=avg_TSPO))+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+
  # geom_text_repel(aes(label = sample))+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/sofia_TSPO_segmentation_draft02/scatter_avg_TSPO_SofiaShort_tot.pdf",width = 5,height = 5)

# split by type
df_plot %>%
  ggplot(aes(x=sofia_segmentation_short,y=avg_TSPO))+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+
  # geom_text_repel(aes(label = sample))+
  theme(axis.text.x = element_text(hjust = 1,angle = 90),strip.background = element_blank())+facet_wrap(~type)
ggsave("out/image/sofia_TSPO_segmentation_draft02/scatter_avg_TSPO_SofiaShort_tot_split.pdf",width = 9,height = 6)

# df_plot2 <- df_tot_fix %>% 
#   group_by(class_top,sample,type,sofia_segmentation_short) %>% 
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
#   mutate(type = fct_reorder(type,dplyr::desc(avg_type)))
# 
# df_plot2 %>%
#   ggplot(aes(x=sofia_segmentation_short,y=avg_TSPO))+
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
#   stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+
#   # geom_text_repel(aes(label = sample))+
#   theme(axis.text.x = element_text(hjust = 1,angle = 45))+
#   facet_grid(type~class_top)
# ggsave("out/image/sofia_TSPO_segmentation_draft02/scatter_avg_TSPO_SofiaShort_tot.pdf",width = 5,height = 5)

df_plot %>%
  filter(type == "A")

df_tot_fix %>% 
  filter(type == "A") %>% 
  group_by(sample,type,sofia_segmentation_short,cell_type) %>% 
  summarise(n = n(),
            median = median(exp),
            avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  print(n = 30)

# is there an overall higher expression of TSPO in the slide overall?
mat_avg <- df_tot_fix %>% 
  group_by(sample,type,sofia_segmentation_short) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  mutate(sample_type = paste0(sample,"_",type)) %>%
  ungroup() %>% 
  dplyr::select(-c(sample,type)) %>%
  pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  column_to_rownames("sofia_segmentation_short")

# mat_short %>%
#   rownames_to_column() %>%
#   write_tsv("out/table/test_tidyverse.tsv")

# # filter out all the columns that are NA
# id_avg <- apply(mat_avg,MARGIN = 2,is.na) %>% 
#   colSums()>0
# 
# mat_avg_filter <- mat_avg[1:5,id_avg]
id_avg <- is.na(mat_avg) %>% 
  colSums() < 7
id_avg2 <- is.na(mat_avg) %>% 
  rowSums() < 10

mat_avg_filter <- mat_avg[id_avg2,id_avg]

# build the annotation on the heatmap
# read in the sample classification
# LUT_sample <- read_csv("data/sample_classification.csv")
sample_ordered_avg <- data.frame(sample = str_extract(colnames(mat_avg_filter),pattern = c("V\\d+"))) %>% 
  left_join(LUT_sample)

table(sample_ordered_avg$type)

sample_ordered_conversion_avg <- sample_ordered_avg %>% pull(type)

column_ha_avg <- HeatmapAnnotation(sample = sample_ordered_conversion_avg,
                                   col = list(sample = c("NAWM"="cyan",
                                                         "BASAL"="gray",
                                                         "CORTEX"="green",
                                                         "CA"="red",
                                                         "A"="orange","CI"="yellow",
                                                         "RM"="blue")))


pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_avg_TSPO_SofiaShort_tot.pdf",width = 6,height = 4)
Heatmap(mat_avg_filter,name = "avg exp TSPO",
        top_annotation = column_ha_avg,
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# do the same as above but only for the type
mat_avg2 <- df_tot_fix %>% 
  group_by(type,sofia_segmentation_short) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  # mutate(sample_type = paste0(sample,"_",type)) %>%
  # ungroup() %>% 
  # dplyr::select(-c(sample,type)) %>%
  # pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  pivot_wider(names_from = type,values_from = avg_TSPO) %>% 
  column_to_rownames("sofia_segmentation_short")

# mat_short %>%
#   rownames_to_column() %>%
#   write_tsv("out/table/test_tidyverse.tsv")

# filter out all the columns that are NA
# id_avg2 <- apply(mat_avg2,MARGIN = 2,is.na) %>% 
#   colSums()>0

id_avg2 <- is.na(mat_avg2) %>% 
  colSums() < 7
id_avg22 <- is.na(mat_avg2) %>% 
  rowSums() < 5

mat_avg2_filter <- mat_avg2[id_avg22,id_avg2]

pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_avg_TSPO_SofiaShort_tot2.pdf",width = 4,height = 3)
Heatmap(mat_avg2_filter,name = "avg exp TSPO",
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# plot expression of TSPO top assignment ----------------------------------
df_tot_fix %>% 
  # average expression per sample, type and cell_type
  group_by(type,sofia_segmentation_short,sample,class_top) %>%
  summarise(avg_TSPO = mean(exp)) %>%
  ungroup() %>%
  # average expression per type fo the sample
  group_by(type) %>%
  mutate(avg_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(type = fct_reorder(type,avg_type)) %>%
  # average expression per classification fo the sample
  group_by(sofia_segmentation_short) %>%
  mutate(avg_segmentation = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(sofia_segmentation_short = fct_reorder(sofia_segmentation_short,avg_segmentation)) %>%
  # average expression per cell_type
  group_by(class_top) %>%
  mutate(avg_cell_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(class_top = fct_reorder(class_top,dplyr::desc(avg_cell_type))) %>%
  # ggplot(aes(x=type,y=avg_TSPO,col=cell_type)) +
  ggplot(aes(x=sofia_segmentation_short,y=avg_TSPO)) +
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  facet_grid(type~class_top)+
  theme(axis.text.x = element_text(hjust = 1,angle = 90),
        strip.background = element_blank())+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+scale_y_log10()
ggsave("out/image/sofia_TSPO_segmentation_draft02/scatter_avg_exp_TSPO_sampleCellType_TopAssigned_SofiaShort.pdf",width = 9,height = 9)

# render as an heatmap
mat_tot0 <- df_tot_fix %>%
  # average expression per sample, type and cell_type
  group_by(sofia_segmentation_short,type,class_top) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  mutate(segmentation_type = paste0(sofia_segmentation_short,"_",type)) %>% 
  dplyr::select(-c(sofia_segmentation_short,type)) %>% 
  pivot_wider(names_from = segmentation_type,values_from = avg_TSPO) %>% 
  column_to_rownames("class_top")

# build the heatmap annotation
area_ordered <- str_extract(colnames(mat_tot0),pattern = c("WM|core|cortex|edge|periplaque|active lesion"))

column_ha_tot0 <- HeatmapAnnotation(sample = str_extract(colnames(mat_tot0),pattern = c("CA|RM|A|CI|NAWM|CORTEX")),
                                    area = area_ordered,
                                    col = list(area = c("WM"="beige",
                                                        "core"="red",
                                                        "cortex"="gray",
                                                        "edge"="orange",
                                                        "periplaque"="green",
                                                        "active lesion"="yellow"),
                                               sample = c("NAWM"="cyan",
                                                          "BASAL"="gray",
                                                          "CORTEX"="green",
                                                          "CA"="red",
                                                          "A"="orange",
                                                          "CI"="yellow",
                                                          "RM"="blue")))

pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_avg_exp_TSPO_sampleCellType_TopAssigned_SofiaShort.pdf",width = 7,height = 4)
Heatmap(mat_tot0,name = "avg exp TSPO",
        top_annotation = column_ha_tot0,
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# do the same but compound the slice in sample types
mat_tot02 <- df_tot_fix %>%
  # average expression per sample, type and cell_type
  group_by(sofia_segmentation_short,class_top) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  # mutate(sample_type = paste0(sample,"_",type)) %>% 
  # dplyr::select(-c(sample,type)) %>% 
  # pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  pivot_wider(names_from = sofia_segmentation_short,values_from = avg_TSPO) %>% 
  column_to_rownames("class_top")

# build the heatmap annotation
pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_avg_exp_TSPO_sampleCellType_TopAssigned_SofiaShort2.pdf",width = 4,height = 4)
Heatmap(mat_tot02,name = "avg exp TSPO",
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# plot expression of TSPO top 50assignment --------------------------------
df_tot_fix %>% 
  # average expression per sample, type and cell_type
  group_by(type,sofia_segmentation_short,sample,class_50) %>%
  summarise(avg_TSPO = mean(exp)) %>%
  ungroup() %>%
  # average expression per type fo the sample
  group_by(type) %>%
  mutate(avg_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(type = fct_reorder(type,avg_type)) %>%
  # average expression per classification fo the sample
  group_by(sofia_segmentation_short) %>%
  mutate(avg_segmentation = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(sofia_segmentation_short = fct_reorder(sofia_segmentation_short,avg_segmentation)) %>%
  # average expression per cell_type
  group_by(class_50) %>%
  mutate(avg_cell_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(class_50 = fct_reorder(class_50,dplyr::desc(avg_cell_type))) %>%
  # ggplot(aes(x=type,y=avg_TSPO,col=cell_type)) +
  ggplot(aes(x=sofia_segmentation_short,y=avg_TSPO)) +
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  facet_grid(type~class_50)+
  theme(axis.text.x = element_text(hjust = 1,angle = 90),
        strip.background = element_blank())+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+scale_y_log10()
ggsave("out/image/sofia_TSPO_segmentation_draft02/scatter_avg_exp_TSPO_sampleCellType_Top50Assigned_SofiaShort.pdf",width = 9,height = 9)

# render as an heatmap
mat_tot500 <- df_tot_fix %>%
  # average expression per sample, type and cell_type
  group_by(sofia_segmentation_short,type,class_50) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  mutate(segmentation_type = paste0(sofia_segmentation_short,"_",type)) %>% 
  dplyr::select(-c(sofia_segmentation_short,type)) %>% 
  pivot_wider(names_from = segmentation_type,values_from = avg_TSPO) %>% 
  column_to_rownames("class_50")

# build the heatmap annotation
area_ordered50 <- str_extract(colnames(mat_tot500),pattern = c("WM|core|cortex|edge|periplaque|active lesion"))

column_ha_tot500 <- HeatmapAnnotation(sample = str_extract(colnames(mat_tot500),pattern = c("CA|RM|A|CI|NAWM|CORTEX")),
                                    area = area_ordered50,
                                    col = list(area = c("WM"="beige",
                                                        "core"="red",
                                                        "cortex"="gray",
                                                        "edge"="orange",
                                                        "periplaque"="green",
                                                        "active lesion"="yellow"),
                                               sample = c("NAWM"="cyan",
                                                          "BASAL"="gray",
                                                          "CORTEX"="green",
                                                          "CA"="red",
                                                          "A"="orange",
                                                          "CI"="yellow",
                                                          "RM"="blue")))

pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_avg_exp_TSPO_sampleCellType_Top50Assigned_SofiaShort.pdf",width = 7,height = 4)
Heatmap(mat_tot500,name = "avg exp TSPO",
        top_annotation = column_ha_tot500,
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# do the same but compound the slice in sample types
mat_tot5002 <- df_tot_fix %>%
  # average expression per sample, type and cell_type
  group_by(sofia_segmentation_short,class_50) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  # mutate(sample_type = paste0(sample,"_",type)) %>% 
  # dplyr::select(-c(sample,type)) %>% 
  # pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  pivot_wider(names_from = sofia_segmentation_short,values_from = avg_TSPO) %>% 
  column_to_rownames("class_50")

# build the heatmap annotation
pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_avg_exp_TSPO_sampleCellType_Top50Assigned_SofiaShort2.pdf",width = 4,height = 4)
Heatmap(mat_tot5002,name = "avg exp TSPO",
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# plot expression of TSPO top 50assignment --------------------------------
df_tot_fix %>% 
  # average expression per sample, type and cell_type
  group_by(type,sofia_segmentation_short,sample,class_70) %>%
  summarise(avg_TSPO = mean(exp)) %>%
  ungroup() %>%
  # average expression per type fo the sample
  group_by(type) %>%
  mutate(avg_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(type = fct_reorder(type,avg_type)) %>%
  # average expression per classification fo the sample
  group_by(sofia_segmentation_short) %>%
  mutate(avg_segmentation = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(sofia_segmentation_short = fct_reorder(sofia_segmentation_short,avg_segmentation)) %>%
  # average expression per cell_type
  group_by(class_70) %>%
  mutate(avg_cell_type = mean(avg_TSPO)) %>%
  ungroup() %>%
  mutate(class_70 = fct_reorder(class_70,dplyr::desc(avg_cell_type))) %>%
  # ggplot(aes(x=type,y=avg_TSPO,col=cell_type)) +
  ggplot(aes(x=sofia_segmentation_short,y=avg_TSPO)) +
  geom_point(position = position_jitter(width = 0.1),alpha = 0.5)+theme_bw()+
  facet_grid(type~class_70)+
  theme(axis.text.x = element_text(hjust = 1,angle = 90),
        strip.background = element_blank())+
  stat_summary(fun = "mean", colour = "red",geom = "crossbar",size = 0.1)+scale_y_log10()
ggsave("out/image/sofia_TSPO_segmentation_draft02/scatter_avg_exp_TSPO_sampleCellType_Top70Assigned_SofiaShort.pdf",width = 9,height = 9)

# render as an heatmap
mat_tot700 <- df_tot_fix %>%
  # average expression per sample, type and cell_type
  group_by(sofia_segmentation_short,type,class_70) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  mutate(segmentation_type = paste0(sofia_segmentation_short,"_",type)) %>% 
  dplyr::select(-c(sofia_segmentation_short,type)) %>% 
  pivot_wider(names_from = segmentation_type,values_from = avg_TSPO) %>% 
  column_to_rownames("class_70")

# build the heatmap annotation
area_ordered70 <- str_extract(colnames(mat_tot700),pattern = c("WM|core|cortex|edge|periplaque|active lesion"))

column_ha_tot700 <- HeatmapAnnotation(sample = str_extract(colnames(mat_tot700),pattern = c("CA|RM|A|CI|NAWM|CORTEX")),
                                      area = area_ordered70,
                                      col = list(area = c("WM"="beige",
                                                          "core"="red",
                                                          "cortex"="gray",
                                                          "edge"="orange",
                                                          "periplaque"="green",
                                                          "active lesion"="yellow"),
                                                 sample = c("NAWM"="cyan",
                                                            "BASAL"="gray",
                                                            "CORTEX"="green",
                                                            "CA"="red",
                                                            "A"="orange",
                                                            "CI"="yellow",
                                                            "RM"="blue")))

pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_avg_exp_TSPO_sampleCellType_Top70Assigned_SofiaShort.pdf",width = 7,height = 4)
Heatmap(mat_tot700,name = "avg exp TSPO",
        top_annotation = column_ha_tot700,
        col = viridis::viridis(20,option = "inferno"))  
dev.off()

# do the same but compound the slice in sample types
mat_tot7002 <- df_tot_fix %>%
  # average expression per sample, type and cell_type
  group_by(sofia_segmentation_short,class_70) %>% 
  summarise(avg_TSPO = mean(exp)) %>% 
  ungroup() %>% 
  # mutate(sample_type = paste0(sample,"_",type)) %>% 
  # dplyr::select(-c(sample,type)) %>% 
  # pivot_wider(names_from = sample_type,values_from = avg_TSPO) %>% 
  pivot_wider(names_from = sofia_segmentation_short,values_from = avg_TSPO) %>% 
  column_to_rownames("class_70")

# build the heatmap annotation
pdf("out/image/sofia_TSPO_segmentation_draft02/heatmap_avg_exp_TSPO_sampleCellType_Top70Assigned_SofiaShort2.pdf",width = 4,height = 4)
Heatmap(mat_tot7002,name = "avg exp TSPO",
        col = viridis::viridis(20,option = "inferno"))  
dev.off()
















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
#                 top_annotation = column_ha_tot,
#                 column_title = name,
#                 col = viridis::viridis(20,option = "inferno"))  
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
