# read in the data --------------------------------------------------------
list_brain <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")

# In order to work with multiple slices in the same Seurat object, we provide the merge function.
data.list.id <- names(list_brain)
data.combined.all <- merge(list_brain[[1]], y = list_brain[-1], add.cell.ids = data.list.id, project = "visium_brain")

# This then enables joint dimensional reduction and clustering on the underlying RNA expression data.
DefaultAssay(data.combined.all) <- "SCT"

# VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
VariableFeatures(data.combined.all) <- purrr::map(list_brain, function(x){
  VariableFeatures(x)
}) %>%
  unlist() %>%
  unique()

data.combined.all <- RunPCA(data.combined.all, verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(verbose = FALSE,resolution = 0.1) %>%
  RunUMAP(dims = 1:30)

# Finally, the data can be jointly visualized in a single UMAP plot. SpatialDimPlot() and SpatialFeaturePlot() will by default plot all slices as columns and groupings/features as rows.

# define many colors
color_id <- alphabet(length(unique(data.combined.all$dataset)))
show_col(color_id)
names(color_id) <- unique(data.combined.all$dataset)

nclust <- data.combined.all[["seurat_clusters"]] |> unique() |> nrow()
cluster_cols <- DiscretePalette(nclust, palette = "polychrome")

(DimPlot(data.combined.all, reduction = "umap",group.by = "dataset") + scale_color_manual(values = color_id)) +
  DimPlot(data.combined.all,
          group.by = "seurat_clusters",
          shuffle = TRUE,
          cols = cluster_cols)

# plot clusters per slide
DimPlot(data.combined.all,
        group.by = "seurat_clusters",split.by = "dataset",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 3)

# show the clusters split by slide
# SpatialPlot(seu, pt.size.factor = 2) + 
#   plot_layout(guides='collect') &
#   theme(legend.position = "none") &
#   scale_fill_manual(values = cluster_cols)

# show the cluster by top anntation of the deconvolution approach
DimPlot(data.combined.all,
        group.by = "seurat_clusters",split.by = "class_top",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 4)

DimPlot(data.combined.all,
        group.by = "seurat_clusters",split.by = "manual_segmentation",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 4)
