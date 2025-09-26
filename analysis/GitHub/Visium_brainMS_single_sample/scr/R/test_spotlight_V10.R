# libraries ---------------------------------------------------------------
library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)

# old function spotlight --------------------------------------------------
SPOTlight_old_functions <- dir("scr/SPOTlight_old/")
lapply(paste0("scr/SPOTlight_old/",SPOTlight_old_functions), function(x){
  source(x)
})

# load the SP dataset -----------------------------------------------------
# this one was already preprocessed usin the seurat workflow
brain <- readRDS("../out_large/brain_V10.rds")

# load the sc reference ---------------------------------------------------
# notice that this one was aready preprocessed
reference <- readRDS("../out_large/CCA_SCTransformed.rds")
DimPlot(reference, group.by = "clusterCellType", label = TRUE)

# here is the recommended preprocessing
# set.seed(123)
# cortex_sc <- Seurat::SCTransform(cortex_sc, verbose = FALSE) %>% 
#   Seurat::RunPCA(., verbose = FALSE) %>%
#   Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)
# 
# Seurat::DimPlot(cortex_sc,
#                 group.by = "subclass",
#                 label = TRUE) + Seurat::NoLegend()



# 0.6 Compute marker genes ------------------------------------------------
# To determine the most important marker genes we can use the function Seurat::FindAllMarkers which will return the markers for each cluster.

# Seurat::Idents(reference) <- "clusterCellType"
# cluster_markers_all <- Seurat::FindAllMarkers(object = reference, 
#                                               assay = "SCT",
#                                               slot = "data",
#                                               verbose = TRUE, 
#                                               only.pos = TRUE)
# 
# saveRDS(object = cluster_markers_all,file = "out/object/reference_CCA_markers.rds")
cluster_markers_all <- readRDS(file = "out/object/reference_CCA_markers.rds")

# 0.6.1 SPOTlight Decomposition -------------------------------------------
set.seed(123)

spotlight_ls <- spotlight_deconvolution(
  se_sc = reference,
  counts_spatial = brain@assays$Spatial@counts,
  clust_vr = "clusterCellType", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 3000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

saveRDS(object = spotlight_ls, file = "out/object/spotlight_ls_V10.rds")


# Read RDS object ---------------------------------------------------------
# spotlight_ls <- readRDS(file = "out/object/spotlight_ls_V10.rds")

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

# # 0.6.2 Assess deconvolution ----------------------------------------------
# # Before even looking at the decomposed spots we can gain insight on how well the model performed by looking at the topic profiles for the cell types.
# 
# # The first thing we can do is look at how specific the topic profiles are for each cell type.
# 
# h <- NMF::coef(nmf_mod[[1]])
# 
# rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
# 
# topic_profile_plts <- dot_plot_profiles_fun(
#   h = h,
#   train_cell_clust = nmf_mod[[2]])
# 
# topic_profile_plts[[2]] + ggplot2::theme(
#   axis.text.x = ggplot2::element_text(angle = 90), 
#   axis.text = ggplot2::element_text(size = 12))

# Next we can take a look at the how the individual topic profiles of each cell within each cell-type behave.
# Here we expect that all the cells from the same cell type show a similar topic profile distribution, if not there might be a bit more substructure in that cluster and we may only be capturing one or the other.

# topic_profile_plts[[1]] + theme(axis.text.x = element_text(angle = 90), 
#                                 axis.text = element_text(size = 12))

# Lastly we can take a look at which genes are the most important for each topic and therefore get an insight into which genes are driving them.

# basis_spotlight <- data.frame(NMF::basis(nmf_mod[[1]]))
# 
# colnames(basis_spotlight) <- unique(stringr::str_wrap(nmf_mod[[2]], width = 30))
# 
# basis_spotlight %>%
#   dplyr::arrange(desc(IMM)) %>%
#   round(., 5) %>% 
#   DT::datatable(., filter = "top")

# 0.7 Visualization -------------------------------------------------------
# Join decomposition with metadata

# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(brain)

decon_df <- decon_mtrx %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

brain@meta.data <- brain@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

# save the meta after deconvolution
brain@meta.data %>% 
  rownames_to_column("barcodes") %>% 
  write_tsv("out/table/V10_SPOTLight.tsv")

# 0.7.1 Specific cell-types -----------------------------------------------
# we can use the standard Seurat::SpatialFeaturePlot to view predicted celltype proportions one at a time.

Seurat::SpatialFeaturePlot(
  object = brain,
  features = c("AST", "IMM", "NEU", "OLIGO","OPC","VAS"),
  alpha = c(0.1, 1))
ggsave("out/image/01_spotlight_deconvolution_panel_V10.pdf",width = 15,height = 10)

# 0.7.2 Spatial scatterpies -----------------------------------------------
# Plot spot composition of all the spots.

cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

spatial_scatterpie(se_obj = brain,
                   cell_types_all = c("AST","IMM","NEU","OLIGO","VAS"),
                   img_path = "../../../raw_data/spaceranger_out/10_SP1/outs/spatial/tissue_lowres_image.png",
                   pie_scale = 0.4)
ggsave("out/image/01_spotlight_deconvolution_piechart_V10.pdf",width = 5,height = 5)

# Plot spot composition of spots containing cell-types of interest

# spatial_scatterpie(se_obj = brain,
#                    cell_types_all = c("AST","IMM","NEU","OLIGO","VAS"),
#                    img_path = "../../../raw_data/spaceranger_out/01_SP1_novaseq/outs/spatial/tissue_lowres_image.png",
#                    cell_types_interest = "IMM",
#                    pie_scale = 0.4)

# # 0.7.3 Spatial interaction graph -----------------------------------------
# # Now that we know which cell types are found within each spot we can make a graph representing spatial interactions where cell types will have stronger edges between them the more often we find them within the same spot. To do this we will only need to run the function get_spatial_interaction_graph, this function prints the plot and returns the elements needed to plot it.
# 
# graph_ntw <- get_spatial_interaction_graph(decon_mtrx = decon_mtrx[,c("AST","IMM","NEU","OLIGO","VAS")])
# 
# # If you want to tune how the graph looks you can do the following or you can check out more options here:
# deg <- degree(graph_ntw, mode="all")
# 
# # Get color palette for difusion
# edge_importance <- E(graph_ntw)$importance
# 
# # Select a continuous palette
# qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]
# 
# # Create a color palette
# getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
# 
# # Get how many values we need
# grad_edge <- seq(0, max(edge_importance), 0.1)
# 
# # Generate extended gradient palette dataframe
# graph_col_df <- data.frame(value = as.character(grad_edge),
#                            color = getPalette(length(grad_edge)),
#                            stringsAsFactors = FALSE)
# # Assign color to each edge
# color_edge <- data.frame(value = as.character(round(edge_importance, 1)), stringsAsFactors = FALSE) %>%
#   dplyr::left_join(graph_col_df, by = "value") %>%
#   dplyr::pull(color)
# 
# # Open a pdf file
# plot(graph_ntw,
#      # Size of the edge
#      edge.width = edge_importance,
#      edge.color = color_edge,
#      # Size of the buble
#      vertex.size = deg,
#      vertex.color = "#cde394",
#      vertex.frame.color = "white",
#      vertex.label.color = "black",
#      vertex.label.family = "Ubuntu", # Font family of the label (e.g.“Times”, “Helvetica”)
#      layout = layout.circle)
# 
# # Lastly one can compute cell-cell correlations to see groups of cells that correlate positively or negatively.
# 
# # Remove cell types not predicted to be on the tissue
# decon_mtrx_sub <- decon_mtrx[, c("AST","IMM","NEU","OLIGO","VAS")]
# decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]
# 
# # Compute correlation
# decon_cor <- cor(decon_mtrx_sub,method = "spearman")
# 
# # Compute correlation P-value
# p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)
# 
# # Visualize
# ggcorrplot::ggcorrplot(
#   corr = decon_cor,
#   p.mat = p.mat[[1]],
#   hc.order = TRUE,
#   type = "full",
#   insig = "blank",
#   lab = TRUE,
#   outline.col = "lightgrey",
#   method = "square",
#   # colors = c("#4477AA", "white", "#BB4444"))
#   colors = c("#6D9EC1", "white", "#E46726"),
#   title = "Predicted cell-cell proportion correlation",
#   legend.title = "Correlation\n(Spearman)") +
#   ggplot2::theme(
#     plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
#     legend.text = ggplot2::element_text(size = 12),
#     legend.title = ggplot2::element_text(size = 15),
#     axis.text.x = ggplot2::element_text(angle = 90),
#     axis.text = ggplot2::element_text(size = 18, vjust = 0.5))


# # 0.8.1 Downsample data ---------------------------------------------------
# # If the dataset is very large we want to downsample it, both in terms of number of cells and number of genes, to train the model. To do this downsampling we want to keep a representative amount of cells per cluster and the most important genes. We show that this downsampling doesn’t affect the performance of the model and greatly speeds up the model training.
# 
# # Downsample scRNAseq to select gene set and number of cells to train the model
# se_sc_down <- downsample_se_obj(se_obj = cortex_sc,
#                                 clust_vr = "subclass",
#                                 cluster_markers = cluster_markers_all,
#                                 cl_n = 100,
#                                 hvg = 3000)
# 0.8.2 Train NMF model
# Once we have the data ready to pass to the model we can train it as shown below.
# 
# start_time <- Sys.time()
# nmf_mod_ls <- train_nmf(cluster_markers = cluster_markers_all, 
#                         se_sc = se_sc_down, 
#                         mtrx_spatial = anterior@assays$Spatial@counts,
#                         clust_vr = "subclass",
#                         ntop = NULL,
#                         hvg = 3000,
#                         transf = "uv",
#                         method = "nsNMF")
# 
# nmf_mod <- nmf_mod_ls[[1]]
# Extract matrices form the model:
#   
#   # get basis matrix W
#   w <- basis(nmf_mod)
# dim(w)
# 
# # get coefficient matrix H
# h <- coef(nmf_mod)
# dim(h)
# Look at cell-type specific topic profile
# 
# rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
# topic_profile_plts <- dot_plot_profiles_fun(
#   h = h,
#   train_cell_clust = nmf_mod_ls[[2]]
# )
# 
# topic_profile_plts[[2]] + theme(axis.text.x = element_text(angle = 90))
# 0.8.3 Spot Deconvolution
# # Extract count matrix
# spot_counts <- anterior@assays$Spatial@counts
# 
# # Subset to genes used to train the model
# spot_counts <- spot_counts[rownames(spot_counts) %in% rownames(w), ]
# Run spots through the basis to get the pertinent coefficients. To do this for every spot we are going to set up a system of linear equations where we need to find the coefficient, we will use non-negative least squares to determine the best coefficient fit.
# 
# ct_topic_profiles <- topic_profile_per_cluster_nmf(h = h,
#                                                    train_cell_clust = nmf_mod_ls[[2]])
# 
# decon_mtrx <- mixture_deconvolution_nmf(nmf_mod = nmf_mod,
#                                         mixture_transcriptome = spot_counts,
#                                         transf = "uv",
#                                         reference_profiles = ct_topic_profiles, 
#                                         min_cont = 0.09)