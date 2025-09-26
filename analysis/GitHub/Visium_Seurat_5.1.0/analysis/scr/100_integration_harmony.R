# AIM ---------------------------------------------------------------------
# try harmony processing for spatial datasets.
# the reference for the processing is suggested here https://www.10xgenomics.com/analysis-guides/correcting-batch-effects-in-visium-data
# The main limitation seems to be the fact that the processing is recommended on a technical replicate.
# some modidifcations are needed to adapt the code to the current object

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

# read in the data --------------------------------------------------------
# read the data from the previous analysis
list_brain <- readRDS("../../Visium_brainMS_single_sample/out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")

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

# 4. Visualize
# First, visualize the data before running batch effect correction. Within this command, we will also be normalizing, scaling the data, and calculating gene and feature variance which will be used to run a PCA and UMAP.
# This need to be modified from the original line becouse the object has changed
data.combined.all <- SCTransform(data.combined.all, assay = "Spatial", verbose = TRUE) %>%
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
saveRDS(data.combined.all_harmony, file = "../out/object/100_data.combined.all_harmony.rds")

# plot --------------------------------------------------------------------
# Martina asked to gather the periplaques categories in one.
data.combined.all_harmony@meta.data$manual_segmentation2 <- 
  data.combined.all_harmony@meta.data %>%
  mutate(manual_segmentation2 = case_when(manual_segmentation %in% c("periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6") ~ "WM",
                                          TRUE ~ manual_segmentation)) %>%
  pull(manual_segmentation2)

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
ggsave("../out/plot/100_UMAPCluster_resolutions.pdf",width = 25,height = 15)

# decide the res 0.1
color_id <- alphabet(length(unique(data.combined.all_harmony$dataset)))
show_col(color_id)
names(color_id) <- unique(data.combined.all_harmony$dataset %>% unlist())

nclust <- data.combined.all_harmony[["seurat_clusters"]] |> unique() |> nrow()
cluster_cols <- DiscretePalette(nclust, palette = "polychrome")

(DimPlot(data.combined.all_harmony, reduction = "umap",group.by = "dataset") + scale_color_manual(values = color_id)) +
  DimPlot(data.combined.all_harmony,
          group.by = "SCT_snn_res.0.1",
          shuffle = TRUE,
          cols = cluster_cols)

# plot clusters per slide
DimPlot(data.combined.all_harmony,
        group.by = "seurat_clusters",split.by = "dataset",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 3)

# show the clusters split by slide
# SpatialPlot(seu, pt.size.factor = 2) + 
#   plot_layout(guides='collect') &
#   theme(legend.position = "none") &
#   scale_fill_manual(values = cluster_cols)

# show the cluster by top anntation of the deconvolution approach
DimPlot(data.combined.all_harmony,
        group.by = "seurat_clusters",split.by = "class_top",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 4)

DimPlot(data.combined.all_harmony,
        group.by = "seurat_clusters",split.by = "manual_segmentation",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 4)

DimPlot(data.combined.all_harmony,
        group.by = "seurat_clusters",split.by = "manual_segmentation2",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 3)

# plot the number of reads on the UMAP
(FeaturePlot(data.combined.all_harmony,features = "nCount_Spatial",order = T)+scale_color_viridis_c(option = "turbo",limits = c(0,5000),oob = scales::squish)) +(FeaturePlot(data.combined.all_harmony,features = "nFeature_Spatial",order = T)+scale_color_viridis_c(option = "turbo",limits = c(0,3000),oob = scales::squish))

# plot as a violin
VlnPlot(data.combined.all_harmony,features = c("nCount_Spatial"),group.by = "manual_segmentation2")+scale_y_continuous(trans = "log1p",breaks = c(0,1,10,100,1000,10000,100000))

(FeaturePlot(data.combined.all_harmony,features = "nCount_SCT",order = T)+scale_color_viridis_c(option = "turbo")) +(FeaturePlot(data.combined.all_harmony,features = "nFeature_SCT",order = T)+scale_color_viridis_c(option = "turbo"))

# top markers -------------------------------------------------------------
DefaultAssay(data.combined.all_harmony) <- "SCT"
data.combined.all_harmony2 <- PrepSCTFindMarkers(data.combined.all_harmony)

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

# save the table of all markers
total_h.markers %>%
  write_tsv("../out/table/100_FindAllMarkers_HarmonySample_res0.1.tsv")

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
ggsave("../../out/image/129_Ditto_HarmonySample_res0.1_MG_subcluster.pdf",width = 15,height = 5)

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
ggsave("../../out/image/129_plot_clusterProp_res0.1_MG_subcluster.pdf",height = 20,width = 4)

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
