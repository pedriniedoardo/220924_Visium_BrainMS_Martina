# AIM ---------------------------------------------------------------------
# try harmony processing for spatial datasets.
# the reference for the processing is suggested here https://www.10xgenomics.com/analysis-guides/correcting-batch-effects-in-visium-data
# The main limitation seems to be the fact that the processing is recommended on a technical replicate.
# some modidifcations are needed to adapt the code to the current object
# differently from the original analysis, I will try to regress out the nCount_Spatial and percent_mt variables

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
library(cowplot)

# read in the data --------------------------------------------------------
# read in the sample classification
LUT_sample <- read_csv("../../Visium_brainMS_single_sample/data/sample_classification.csv") %>%
  mutate(dataset = paste0("brain_",sample))

# read the data from the previous analysis
list_brain <- readRDS("../../Visium_brainMS_single_sample/out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")
list_brain_bayespace <- readRDS(file = "../../Visium_brainMS_single_sample/out/object/list_brain_all_BayesSpace1000_q10.rds")

# merge the individual objetcs to create a single count matrix
data.list.id <- names(list_brain)
data.combined.all <- merge(list_brain[[1]], y = list_brain[-1], add.cell.ids = data.list.id, project = "visium_brain")

# check the size of the dataset
data.combined.all

# confirm the total number of cells
lapply(list_brain,function(x){
  dim(x@assays$Spatial@counts)[2]
}) %>%
  unlist() %>%
  sum()

# QC
data.combined.all <- PercentageFeatureSet(data.combined.all,
                                          pattern = "^MT-|^Mt-|^mt-",
                                          col.name = "percent_mt")

# 4. Visualize
# First, visualize the data before running batch effect correction. Within this command, we will also be normalizing, scaling the data, and calculating gene and feature variance which will be used to run a PCA and UMAP.
# This need to be modified from the original line becouse the object has changed
data.combined.all <- SCTransform(data.combined.all, assay = "Spatial",vars.to.regress = c("nCount_Spatial","percent_mt"),verbose = TRUE) %>%
  RunPCA(verbose = T) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(verbose = T,resolution = 0.1) %>%
  RunUMAP(dims = 1:30)

DimPlot(data.combined.all, group.by = "dataset")

# 5. Run Harmony
# The next few commands are adapted from the Seurat vignette.
data.combined.all_harmony <- RunHarmony(data.combined.all, group.by.vars = "dataset") %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = seq(from=0.1,to=1,by=0.1),verbose = T)

# save the object
saveRDS(data.combined.all_harmony, file = "../out/object/100_data.combined.all_harmony_regressOut.rds")

# plot --------------------------------------------------------------------
# Martina asked to gather the periplaques categories in one.
data.combined.all_harmony@meta.data$manual_segmentation2 <- 
  data.combined.all_harmony@meta.data %>%
  mutate(manual_segmentation2 = case_when(manual_segmentation %in% c("periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6") ~ "WM",
                                          TRUE ~ manual_segmentation)) %>%
  pull(manual_segmentation2)

# martina suggested to split the edge also by its pathological class
data.combined.all_harmony@meta.data$manual_segmentation3 <- data.combined.all_harmony@meta.data %>%
  left_join(LUT_sample,by="dataset") %>%
  mutate(manual_segmentation3 = case_when(manual_segmentation2 == "edge" ~ paste0("edge_",type),
                                          TRUE ~ manual_segmentation2)) %>%
  pull(manual_segmentation3)

# plot the UMAP with all the resolutions runs
id_resolution <- str_subset(colnames(data.combined.all_harmony@meta.data),pattern = "SCT_snn_res") %>%
  sort()

# x <- id_resolution[1]
list_plot <- lapply(id_resolution,function(x){
  
  # define many colors
  # color_id <- alphabet(length(unique(data.combined.all_harmony[[x]] %>% unlist())))
  # show_col(color_id)
  # names(color_id) <- unique(data.combined.all_harmony[[x]] %>% unlist())
  
  # make the plot
  plot <- DimPlot(data.combined.all_harmony,
                  reduction = "umap",
                  group.by = x,
                  label = F,
                  raster = T)
  
  # scale_color_manual(values = color_id)
  
  # return plot
  return(plot)
})

wrap_plots(list_plot)
ggsave("../out/plot/100_UMAPCluster_resolutions_regressOut.pdf",width = 25,height = 15)

# decide the res 0.1
color_id <- tableau20(length(unique(data.combined.all_harmony$dataset)))
show_col(color_id)
names(color_id) <- unique(data.combined.all_harmony$dataset %>% unlist())

# nclust <- data.combined.all_harmony[["seurat_clusters"]] |> unique() |> nrow()
# cluster_cols <- DiscretePalette(nclust, palette = "polychrome")

nclust <- data.combined.all_harmony[["seurat_clusters"]] |> unique() |> nrow()
cluster_cols <- tableau20(nclust)
show_col(cluster_cols)

(DimPlot(data.combined.all_harmony, reduction = "umap",group.by = "dataset") + scale_color_manual(values = color_id)) +
  DimPlot(data.combined.all_harmony,
          group.by = "SCT_snn_res.0.1",
          shuffle = TRUE,
          cols = cluster_cols)

# save the plot
ggsave("../out/plot/100_UMAP_res0.1_regressOut.pdf",width = 10,height = 4)

# plot clusters per slide
DimPlot(data.combined.all_harmony,
        group.by = "SCT_snn_res.0.1",split.by = "dataset",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 3)

# show the clusters split by slide
# SpatialPlot(seu, pt.size.factor = 2) + 
#   plot_layout(guides='collect') &
#   theme(legend.position = "none") &
#   scale_fill_manual(values = cluster_cols)

# show the cluster by top anntation of the deconvolution approach
DimPlot(data.combined.all_harmony,
        group.by = "SCT_snn_res.0.1",split.by = "class_top",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 4)

DimPlot(data.combined.all_harmony,
        group.by = "SCT_snn_res.0.1",split.by = "manual_segmentation",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 4)

DimPlot(data.combined.all_harmony,
        group.by = "SCT_snn_res.0.1",split.by = "manual_segmentation2",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 3)

DimPlot(data.combined.all_harmony,
        group.by = "SCT_snn_res.0.1",split.by = "manual_segmentation3",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 3)

# plot the number of reads on the UMAP
(FeaturePlot(data.combined.all_harmony,features = "nCount_Spatial",order = T)+scale_color_viridis_c(option = "turbo",limits = c(0,5000),oob = scales::squish)) +(FeaturePlot(data.combined.all_harmony,features = "nFeature_Spatial",order = T)+scale_color_viridis_c(option = "turbo",limits = c(0,3000),oob = scales::squish))

# plot as a violin
VlnPlot(data.combined.all_harmony,features = c("nCount_Spatial"),group.by = "manual_segmentation2",pt.size = 0)+scale_y_continuous(trans = "log1p",breaks = c(0,1,10,100,1000,10000,100000))
# save this plot
ggsave("../out/plot/100_VlnPlot_nCountSpatial_res0.1_regressOut.pdf",width = 4,height = 3)

# hide the points
VlnPlot(data.combined.all_harmony,features = c("nCount_Spatial"),group.by = "manual_segmentation2")+scale_y_continuous(trans = "log1p",breaks = c(0,1,10,100,1000,10000,100000))

(FeaturePlot(data.combined.all_harmony,features = "nCount_SCT",order = T)+scale_color_viridis_c(option = "turbo")) +(FeaturePlot(data.combined.all_harmony,features = "nFeature_SCT",order = T)+scale_color_viridis_c(option = "turbo"))

# Dotplot custom ----------------------------------------------------------
shortlist_features_list_test <- list(
  WM = c("MBP", "PLP1", "MOBP"),
  WM_MT = c("MT-CO2","MT-ND4"),
  WM_IMM = c("AQP1","SERPINA3","APOE"),
  CORE = c("NEAT1", "SPARC", "HSP90AB1"),
  CORTEX = c("NRGN","GPM6A","CALM1")
)

# DotPlot(data.combined.all_harmony, features = shortlist_features_list_test, dot.scale = 8,cluster.idents = T,)

GOI <- unique(unlist(shortlist_features_list_test))
GOI

# specify the identity of the splitting variable
# define the grouping variable
Idents(data.combined.all_harmony) <- "SCT_snn_res.0.1"

# generate the dotplot using the unique gene name. this is used to extract the average expressions, the percentage of expression and classification
DefaultAssay(data.combined.all_harmony) <- "SCT"
# DefaultAssay(data.combined.all_harmony) <- "Spatial"

test <- DotPlot(data.combined.all_harmony, features = GOI,cluster.idents = T)

# pull the info for each gene for each cluster based on the signature grouping
df_test <- lapply(shortlist_features_list_test,function(x){
  test$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

head(df_test)

# plot custom dotplot. plot mapping radius
df_test %>%
  # force the order
  mutate(id = factor(id,levels = c(0,3,4,1,2))) %>% 
  mutate(cell_type = factor(cell_type,levels = c("WM","WM_MT","WM_IMM","CORE","CORTEX"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled)) +
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90)) +
  scale_color_gradient(low = "lightgrey",high = "blue")

ggsave("../out/plot/100_DotPlot_custom_res0.1_regressOut.pdf",width = 8,height = 4)

# top markers -------------------------------------------------------------
DefaultAssay(data.combined.all_harmony) <- "SCT"
# data.combined.all_harmony2 <- PrepSCTFindMarkers(data.combined.all_harmony)

# all_marks <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# top_markers <- all_marks |>
#   mutate(order_value = avg_log2FC * -log10(p_val_adj + 1e-300)) |>
#   group_by(cluster) |>
#   slice_max(n = 3, order_by = order_value)
# 
# DotPlot(seu, features = unique(top_markers$gene)) +
#   scale_x_discrete(guide = guide_axis(angle = 45))

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Idents(data.combined.all_harmony) <- "SCT_snn_res.0.1"
total_h.markers <- RunPrestoAll(data.combined.all_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Idents(data.combined.all_harmony2) <- "SCT_snn_res.0.1"
# total_h.markers2 <- RunPrestoAll(data.combined.all_harmony2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
total_h.markers %>%
  write_tsv("../out/table/100_FindAllMarkers_HarmonySample_res0.1_regressOut.tsv")

# # pick the top 100 markers per cluster
# total_h.markers %>%
#   group_by(cluster) %>%
#   dplyr::slice(1:100) %>%
#   write_tsv("../../out/table/129_FindAllMarkers_HarmonySample_res0.1_MG_subcluster_top100.tsv")
# 
# sobj_total_h.markers %>%
#   group_by(cluster) %>%
#   dplyr::slice(1:100) %>%
#   filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
#   filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
#   filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
#   write_tsv("../../out/table/129_FindAllMarkers_HarmonySample_res0.1_MG_subcluster_top100_noRIBOandMT.tsv")

# try plotting the top markers
top_specific_markers <- total_h.markers %>%
  # filter ribosomal and mt genes
  # filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(data.combined.all_harmony,
                       vars = unique(top_specific_markers$gene), 
                       group.by = "SCT_snn_res.0.1")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave("../out/plot/100_Ditto_HarmonySample_res0.1_MG_subcluster_regressOut.pdf",width = 15,height = 5)

# plot the proportions
df_summary <- data.combined.all_harmony@meta.data %>%
  group_by(dataset,manual_segmentation2,SCT_snn_res.0.1) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(dataset,manual_segmentation2) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

df_summary %>%
  ggplot(aes(x=manual_segmentation2,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2),shape=1)+
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1,angle = 45)) +
  facet_wrap(~SCT_snn_res.0.1,ncol = 1,scales = "free")+
  scale_y_sqrt()
ggsave("../out/plot/100_plot_clusterProp_res0.1_MG_subcluster_regressOut.pdf",height = 20,width = 4)

mat_zscore <- data.combined.all_harmony@meta.data %>%
  filter(!is.na(manual_segmentation2)) %>%
  group_by(manual_segmentation2,SCT_snn_res.0.1) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(manual_segmentation2) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  group_by(SCT_snn_res.0.1) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  select(manual_segmentation2,SCT_snn_res.0.1,zscore) %>%
  pivot_wider(names_from = manual_segmentation2,values_from = zscore) %>%
  column_to_rownames(var = "SCT_snn_res.0.1")

Heatmap(mat_zscore)

data.combined.all_harmony@meta.data %>%
  filter(!is.na(manual_segmentation2)) %>%
  group_by(manual_segmentation2,SCT_snn_res.0.1) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(manual_segmentation2) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  group_by(SCT_snn_res.0.1) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  write_tsv("../out/table/100_manualSegmentation2_cluster_proportion_regressOut.tsv")

data.combined.all_harmony@meta.data %>%
  filter(!is.na(manual_segmentation2)) %>%
  group_by(dataset,manual_segmentation2,SCT_snn_res.0.1) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(dataset) %>%
  mutate(tot_slide = sum(n)) %>%
  write_tsv("../out/table/100_slide_manualSegmentation2_cluster_proportion_regressOut.tsv")


# do the same also for the manual_segmentation3
# plot the proportions
df_summary2 <- data.combined.all_harmony@meta.data %>%
  # remove the NA classification
  # filter(!is.na(manual_segmentation3)) %>%
  group_by(dataset,manual_segmentation3,SCT_snn_res.0.1) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(dataset,manual_segmentation3) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

df_summary2 %>%
  ggplot(aes(x=manual_segmentation3,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2),shape=1)+
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1,angle = 45)) +
  facet_wrap(~SCT_snn_res.0.1,scales = "free")+
  scale_y_sqrt()
ggsave("../out/plot/100_plot_clusterProp_res0.1_MG_subcluster_regressOut_manulaSegmentation3.pdf",height = 8,width = 12)

# try also to split by  manual annotation rather then cluster id
df_summary2 %>%
  ggplot(aes(x=SCT_snn_res.0.1,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2),shape=1)+
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1,angle = 45)) +
  facet_wrap(~manual_segmentation3,scales = "free")+
  scale_y_sqrt()
ggsave("../out/plot/100_plot_clusterProp_res0.1_MG_subcluster_regressOut_manulaSegmentation32.pdf",height = 8,width = 12)


mat_zscore2 <- data.combined.all_harmony@meta.data %>%
  filter(!is.na(manual_segmentation3)) %>%
  group_by(manual_segmentation3,SCT_snn_res.0.1) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(manual_segmentation3) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  group_by(SCT_snn_res.0.1) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  select(manual_segmentation3,SCT_snn_res.0.1,zscore) %>%
  pivot_wider(names_from = manual_segmentation3,values_from = zscore) %>%
  column_to_rownames(var = "SCT_snn_res.0.1")

Heatmap(mat_zscore2)

# show the raw proportions of cells
mat_zscore22 <- data.combined.all_harmony@meta.data %>%
  filter(!is.na(manual_segmentation3)) %>%
  group_by(manual_segmentation3,SCT_snn_res.0.1) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(manual_segmentation3) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(manual_segmentation3,SCT_snn_res.0.1,prop) %>%
  pivot_wider(names_from = manual_segmentation3,values_from = prop) %>%
  column_to_rownames(var = "SCT_snn_res.0.1")

Heatmap(mat_zscore22,col = viridis::turbo(100))

# martina also wanted to see the values scaled by location
mat_zscore3 <- data.combined.all_harmony@meta.data %>%
  filter(!is.na(manual_segmentation3)) %>%
  group_by(manual_segmentation3,SCT_snn_res.0.1) %>%
  summarise(n = n(),.groups = "drop") %>%
  group_by(manual_segmentation3) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  # group_by(SCT_snn_res.0.1) %>%
  group_by(manual_segmentation3) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  select(manual_segmentation3,SCT_snn_res.0.1,zscore) %>%
  pivot_wider(names_from = SCT_snn_res.0.1,values_from = zscore) %>%
  column_to_rownames(var = "manual_segmentation3")

Heatmap(mat_zscore3)

# add the deconvolution score information per cluster ---------------------
# from each barcode I also have access to the deconvolution per cell type try to see if each cluster is particularly enriched for a specific cell type.

# split by cluster id
data.combined.all_harmony@meta.data %>%
  select(SCT_snn_res.0.1,AST,IMM,LYM,NEU,OLIGO,OPC,VAS,manual_segmentation3,dataset) %>%
  pivot_longer(cols = -c(SCT_snn_res.0.1,manual_segmentation3,dataset),names_to = "cell_type",values_to = "score") %>%
  group_by(SCT_snn_res.0.1,cell_type,dataset) %>%
  summarise(avg_score = mean(score),.groups = "drop") %>%
  ggplot(aes(x=cell_type,y=avg_score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),shape=1) +
  facet_wrap(~SCT_snn_res.0.1,scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),strip.background = element_blank())

# split by cell annotation
data.combined.all_harmony@meta.data %>%
  select(SCT_snn_res.0.1,AST,IMM,LYM,NEU,OLIGO,OPC,VAS,manual_segmentation3,dataset) %>%
  pivot_longer(cols = -c(SCT_snn_res.0.1,manual_segmentation3,dataset),names_to = "cell_type",values_to = "score") %>%
  group_by(SCT_snn_res.0.1,cell_type,dataset) %>%
  summarise(avg_score = mean(score),.groups = "drop") %>%
  ggplot(aes(x=SCT_snn_res.0.1,y=avg_score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),shape=1) +
  facet_wrap(~cell_type,scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),strip.background = element_blank())

mat_zscore4 <- data.combined.all_harmony@meta.data %>%
  filter(!is.na(manual_segmentation3)) %>%
  select(SCT_snn_res.0.1,AST,IMM,LYM,NEU,OLIGO,OPC,VAS,manual_segmentation3,dataset) %>%
  pivot_longer(cols = -c(SCT_snn_res.0.1,manual_segmentation3,dataset),names_to = "cell_type",values_to = "score") %>%
  group_by(SCT_snn_res.0.1,cell_type) %>%
  summarise(avg_score = mean(score),.groups = "drop") %>%
  group_by(SCT_snn_res.0.1) %>%
  mutate(zscore = (avg_score-mean(avg_score))/sd(avg_score)) %>%
  select(cell_type,SCT_snn_res.0.1,zscore) %>%
  pivot_wider(names_from = cell_type,values_from = zscore) %>%
  column_to_rownames(var = "SCT_snn_res.0.1")

Heatmap(mat_zscore4)


mat_zscore5 <- data.combined.all_harmony@meta.data %>%
  filter(!is.na(manual_segmentation3)) %>%
  select(SCT_snn_res.0.1,AST,IMM,LYM,NEU,OLIGO,OPC,VAS,manual_segmentation3,dataset) %>%
  pivot_longer(cols = -c(SCT_snn_res.0.1,manual_segmentation3,dataset),names_to = "cell_type",values_to = "score") %>%
  group_by(SCT_snn_res.0.1,cell_type) %>%
  summarise(avg_score = mean(score),.groups = "drop") %>%
  group_by(cell_type) %>%
  mutate(zscore = (avg_score-mean(avg_score))/sd(avg_score)) %>%
  select(cell_type,SCT_snn_res.0.1,zscore) %>%
  pivot_wider(names_from = SCT_snn_res.0.1,values_from = zscore) %>%
  column_to_rownames(var = "cell_type")

Heatmap(mat_zscore5)

# -------------------------------------------------------------------------
# try plotting the cluster id on the original object
# define the scale factor for the points. this is needed for older objects assay v3 running on Seurat 5.1.0
# ref: https://github.com/satijalab/seurat/issues/9049
spot_size <- list_brain$V01@images$slice1@scale.factors$fiducial/list_brain$V01@images$slice1@scale.factors$hires

SpatialDimPlot(list_brain$V01, group.by = "manual_segmentation",label = F,pt.size = spot_size)
SpatialDimPlot(list_brain_bayespace$V01, group.by = "BayesSpace",label = F,pt.size = spot_size)

meta_V01 <- list_brain$V01@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode)

meta_V01_bayespace <- list_brain_bayespace$V01@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode,BayesSpace)

meta_integrated <- data.combined.all_harmony@meta.data %>%
  rownames_to_column("sample_barcode") %>%
  separate(sample_barcode,into = c("sample_id","barcode"),sep = "_",remove = T) %>%
  filter(sample_id == "V01")

# join all the info
meta_V01_update <- meta_V01 %>%
  left_join(meta_V01_bayespace,by = "barcode") %>%
  left_join(meta_integrated,by = "barcode") %>%
  column_to_rownames("barcode")

test <- list_brain$V01
test@meta.data <- meta_V01_update

# plot the UMAP with all the annotations
p11 <- SpatialDimPlot(test, group.by = "BayesSpace",label = F,pt.size = spot_size)
p21 <- SpatialDimPlot(test, group.by = "manual_segmentation3",label = F,pt.size = spot_size)
p31 <- SpatialDimPlot(test, group.by = "SCT_snn_res.0.1",label = F,pt.size = spot_size)

Idents(test) <- "BayesSpace"
p12 <- SpatialDimPlot(test,
                      cells.highlight = CellsByIdentities(object = test),
                      facet.highlight = TRUE,pt.size = spot_size)

Idents(test) <- "manual_segmentation3"
p22 <- SpatialDimPlot(test,
                      cells.highlight = CellsByIdentities(object = test),
                      facet.highlight = TRUE,pt.size = spot_size)

Idents(test) <- "SCT_snn_res.0.1"
p32 <- SpatialDimPlot(test,
                      cells.highlight = CellsByIdentities(object = test),
                      facet.highlight = TRUE,pt.size = spot_size)

(p11 + p12) + plot_annotation("BayesSpace",theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("../out/plot/100_SpatialDimPlot_BayesSpace.pdf",width = 25,height = 10)
(p21 + p22) + plot_annotation("manual segmentation",theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("../out/plot/100_SpatialDimPlot_ManualSegmentation3.pdf",width = 25,height = 10)
(p31 + p32) + plot_annotation("integrated clustering",theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("../out/plot/100_SpatialDimPlot_Integratedres0.1.pdf",width = 25,height = 10)


# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# test
a <- c('potato', 'tomotto', 'chips', 'baloon')
b <- c('car', 'chips', 'bird', 'salt')

jaccard(a, b)

# build the dataset for the correlatino plot
df_crossing <- crossing(manual_segmentation = unique(meta_V01_update$manual_segmentation3),
                        integrated_clustering = unique(meta_V01_update$SCT_snn_res.0.1))
df_crossing

# build the scatter plot
df_jaccard_score <- pmap(list(id1 = df_crossing$manual_segmentation,
                              id2 = df_crossing$integrated_clustering), function(id1,id2){
                                
                                # calculate the jaccard score
                                a <- meta_V01_update %>%
                                  rownames_to_column("barcode") %>%
                                  filter(manual_segmentation3 == id1) %>% pull(barcode)
                                b <- meta_V01_update %>%
                                  rownames_to_column("barcode") %>%
                                  filter(SCT_snn_res.0.1 == id2) %>% pull(barcode)
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame(manual_segmentation = id1,
                                                 integrated_cluster = id2,
                                                 jaccard_score = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

head(df_jaccard_score)

# shape it as a matrix
mat_jaccard_score <- df_jaccard_score %>%
  pivot_wider(names_from = manual_segmentation,values_from = jaccard_score) %>%
  column_to_rownames("integrated_cluster")

mat_jaccard_score


# plot the matrix
ht_02 <- Heatmap(mat_jaccard_score,
                 name = "Jaccard score",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 100),
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize = 8),
                 column_names_side = "bottom",
                 column_names_gp = gpar(fontsize = 8),
                 row_dend_reorder = FALSE,
                 column_dend_reorder = FALSE,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 show_column_names = T,
                 show_row_names = T)

ht_02

df_crossing2 <- crossing(manual_segmentation = unique(meta_V01_update$manual_segmentation3),
                         bayespace = unique(meta_V01_update$BayesSpace))
df_crossing2

# build the scatter plot
df_jaccard_score2 <- pmap(list(id1 = df_crossing2$manual_segmentation,
                               id2 = df_crossing2$bayespace), function(id1,id2){
                                
                                # calculate the jaccard score
                                a <- meta_V01_update %>%
                                  rownames_to_column("barcode") %>%
                                  filter(manual_segmentation3 == id1) %>% pull(barcode)
                                b <- meta_V01_update %>%
                                  rownames_to_column("barcode") %>%
                                  filter(BayesSpace == id2) %>% pull(barcode)
                                
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame(manual_segmentation = id1,
                                                 BayeSpace = id2,
                                                 jaccard_score = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

head(df_jaccard_score2)

# shape it as a matrix
mat_jaccard_score2 <- df_jaccard_score2 %>%
  pivot_wider(names_from = manual_segmentation,values_from = jaccard_score) %>%
  column_to_rownames("BayeSpace")

# plot the matrix
ht_12 <- Heatmap(mat_jaccard_score2,
                 name = "Jaccard score",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 100),
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize = 8),
                 column_names_side = "bottom",
                 column_names_gp = gpar(fontsize = 8),
                 row_dend_reorder = FALSE,
                 column_dend_reorder = FALSE,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 show_column_names = T,
                 show_row_names = T)

ht_12



# save the annotation object ----------------------------------------------
saveRDS(data.combined.all_harmony, file = "../out/object/100_data.combined.all_harmony_regressOut_fullAnno.rds")

# save the metadata of the object
meta <- data.combined.all_harmony@meta.data %>%
  rownames_to_column("barcode")

meta %>%
  write_tsv("../out/table/100_meta_data.combined.all_harmony_regressOut_fullAnno.tsv")
