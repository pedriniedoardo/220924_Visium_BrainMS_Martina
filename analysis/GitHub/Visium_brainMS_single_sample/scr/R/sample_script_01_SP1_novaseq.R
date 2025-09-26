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

# setup parallel ----------------------------------------------------------
# library(future)
# plan("multiprocess",workers=12)
# plan()

# read in the data --------------------------------------------------------
# list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
# Load10X_Spatial(data.dir = data_dir)
data_dir <- "../../../raw_data/spaceranger_out/01_SP1_novaseq/outs/"
brain <- Load10X_Spatial(data.dir = data_dir)

# Data preprocessing ------------------------------------------------------
# The initial preprocessing steps that we perform on the spot by gene expression data are similar to a typical scRNA-seq experiment. We first need to normalize the data in order to account for variance in sequencing depth across data points. We note that the variance in molecular counts / spot can be substantial for spatial datasets, particularly if there are differences in cell density across the tissue. We see substantial heterogeneity here, which requires effective normalization.

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
ggsave("out/image/01_ncount.pdf",width = 10,height = 5)

# These plots demonstrate that the variance in molecular counts across spots is not just technical in nature, but also is dependent on the tissue anatomy. For example, regions of the tissue that are depleted for neurons (such as the cortical white matter), reproducibly exhibit lower molecular counts. As a result, standard approaches (such as the LogNormalize() function), which force each data point to have the same underlying ‘size’ after normalization, can be problematic.

# As an alternative, we recommend using sctransform (Hafemeister and Satija, Genome Biology 2019), which which builds regularized negative binomial models of gene expression in order to account for technical artifacts while preserving biological variance. For more details on sctransform, please see the paper here and the Seurat vignette here. sctransform normalizes the data, detects high-variance features, and stores the data in the SCT assay.

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain@assays

# How do results compare to log-normalization?
# To explore the differences in normalization methods, we examine how both the sctransform and log normalization results correlate with the number of UMIs. For this comparison, we first rerun sctransform to store values for all genes and run a log-normalization procedure via NormalizeData().
# rerun normalization to store sctransform residuals for all genes
brain <- SCTransform(brain, assay = "Spatial", return.only.var.genes = FALSE, verbose = T)
# also run standard log normalization for comparison
brain <- NormalizeData(brain, verbose = FALSE, assay = "Spatial")
# Computes the correlation of the log normalized data and sctransform residuals with the
# number of UMIs
brain <- GroupCorrelation(brain, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
brain <- GroupCorrelation(brain, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p1 <- GroupCorrelationPlot(brain, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(brain, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p1 + p2

# For the boxplots above, we calculate the correlation of each feature (gene) with the number of UMIs (the nCount_Spatial variable here). We then place genes into groups based on their mean expression, and generate boxplots of these correlations. You can see that log-normalization fails to adequately normalize genes in the first three groups, suggesting that technical factors continue to influence normalized expression estimates for highly expressed genes. In contrast, sctransform normalization substantially mitigates this effect.

# Gene expression visualization -------------------------------------------
# In Seurat, we have functionality to explore and interact with the inherently visual nature of spatial data. The SpatialFeaturePlot() function in Seurat extends FeaturePlot(), and can overlay molecular data on top of tissue histology. For example, in this data set of the mouse brain, the gene Hpca is a strong hippocampus marker and Ttr is a marker of the choroid plexus.

SpatialFeaturePlot(brain, features = c("SPP1", "MBP"))
ggsave("out/image/02_GOI.pdf",width = 10,height = 5)

# The default parameters in Seurat emphasize the visualization of molecular data. However, you can also adjust the size of the spots (and their transparency) to improve the visualization of the histology image, by changing the following parameters:
# pt.size.factor- This will scale the size of the spots. Default is 1.6
# alpha - minimum and maximum transparency. Default is c(1, 1).
# Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
p1 <- SpatialFeaturePlot(brain, features = c("SPP1"), pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = c("SPP1"), alpha = c(0.1, 1))
p1 + p2
ggsave("out/image/03_SPP1.pdf",width = 10,height = 5)

p1 <- SpatialFeaturePlot(brain, features = c("MBP"), pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = c("MBP"), alpha = c(0.1, 1))
p1 + p2
ggsave("out/image/03_MBP.pdf",width = 10,height = 5)

# Dimensionality reduction, clustering, and visualization -----------------
# We can then proceed to run dimensionality reduction and clustering on the RNA expression data, using the same workflow as we use for scRNA-seq analysis.

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# We can then visualize the results of the clustering either in UMAP space (with DimPlot()) or overlaid on the image with SpatialDimPlot().

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
ggsave("out/image/04_UMAP.pdf",width = 10,height = 5)

# As there are many colors, it can be challenging to visualize which voxel belongs to which cluster. We have a few strategies to help with this. Setting the label parameter places a colored box at the median of each cluster (see the plot above).

# You can also use the cells.highlight parameter to demarcate particular cells of interest on a SpatialDimPlot(). This can be very useful for distinguishing the spatial localization of individual clusters, as we show below:

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(0:8)), facet.highlight = TRUE, ncol = 3)
ggsave("out/image/05_PanelClusters.pdf",width = 15,height = 15)

# Interactive plotting ----------------------------------------------------
# We have also built in a number of interactive plotting capabilities. Both SpatialDimPlot() and SpatialFeaturePlot() now have an interactive parameter, that when set to TRUE, will open up the Rstudio viewer pane with an interactive Shiny plot. The example below demonstrates an interactive SpatialDimPlot() in which you can hover over spots and view the cell name and current identity class (analogous to the previous do.hover behavior).

SpatialDimPlot(brain, interactive = TRUE)

# For SpatialFeaturePlot(), setting interactive to TRUE brings up an interactive pane in which you can adjust the transparency of the spots, the point size, as well as the Assay and feature being plotted. After exploring the data, selecting the done button will return the last active plot as a ggplot object.

SpatialFeaturePlot(brain, features = "SPP1", interactive = TRUE)

# The LinkedDimPlot() function links the UMAP representation to the tissue image representation and allows for interactive selection. For example, you can select a region in the UMAP plot and the corresponding spots in the image representation will be highlighted.

LinkedDimPlot(brain)

# Identification of Spatially Variable Features ---------------------------
# Seurat offers two workflows to identify molecular features that correlate with spatial location within a tissue. The first is to perform differential expression based on pre-annotated anatomical regions within the tissue, which may be determined either from unsupervised clustering or prior knowledge. This strategy works will in this case, as the clusters above exhibit clear spatial restriction.

# de_markers <- FindMarkers(brain, ident.1 = 4, ident.2 = 6)
de_markers4 <- FindMarkers(brain, ident.1 = 4,logfc.threshold = 0,only.pos = T)
SpatialFeaturePlot(object = brain, features = rownames(de_markers4)[1:9], alpha = c(0.1, 1), ncol = 3)
ggsave("out/image/06_TopCluster4.pdf",width = 15,height = 15)

# An alternative approach, implemented in FindSpatiallyVariables(), is to search for features exhibiting spatial patterning in the absence of pre-annotation. The default method (method = 'markvariogram), is inspired by the Trendsceek, which models spatial transcriptomics data as a mark point process and computes a ‘variogram’, which identifies genes whose expression level is dependent on their spatial location. More specifically, this process calculates gamma(r) values measuring the dependence between two spots a certain “r” distance apart. By default, we use an r-value of ‘5’ in these analyses, and only compute these values for variable genes (where variation is calculated independently of spatial location) to save time.

# We note that there are multiple methods in the literature to accomplish this task, including SpatialDE, and Splotch. We encourage interested users to explore these methods, and hope to add support for them in the near future.
# see the suggestino on Using parellelization for findSpatiallyVariableFeatures #2856
t <- Sys.time()
brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "markvariogram",verbose = T)
Sys.time()-t
# Now we visualize the expression of the top 6 features identified by this measure.

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

# since the computation is quite intense save the object up to this poing
saveRDS(brain,"../out_large/brain_SP01.rds")

# Subset out anatomical regions -------------------------------------------
# As with single-cell objects, you can subset the object to focus on a subset of data. Here, we approximately subset the frontal brain. This process also facilitates the integration of these data with a cortical scRNA-seq dataset in the next section. First, we take a subset of clusters, and then further segment based on exact positions. After subsetting, we can visualize the cortical cells either on the full image, or a cropped image.

cortex <- subset(brain, idents = c(5))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400 | image_imagecol < 150))
# cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
# cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
# cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2

# Integration with single-cell data ---------------------------------------
# At ~50um, spots from the visium assay will encompass the expression profiles of multiple cells. For the growing list of systems where scRNA-seq data is available, users may be interested to ‘deconvolute’ each of the spatial voxels to predict the underlying composition of cell types. In preparing this vignette, we tested a wide variety of decovonlution and integration methods, using a reference scRNA-seq dataset of ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, generated with the SMART-Seq2 protocol. We consistently found superior performance using integration methods (as opposed to deconvolution methods), likely because of substantially different noise models that characterize spatial and single-cell datasets, and integration methods are specifiically designed to be robust to these differences. We therefore apply the ‘anchor’-based integration workflow introduced in Seurat v3, that enables the probabilistic transfer of annotations from a reference to a query set. We therefore follow the label transfer workflow introduced here, taking advantage of sctransform normalization, but anticipate new methods to be developed to accomplish this task.

# # in this case I am loading the dataset priduced by Martina
# reference <- readRDS("/media/edo/INTENSO/HSR/project_absinta/scRNAseq/211201_Nature_2021/analysis/GitHub/data/all20_integrated_clean_metadata.rds")
# Idents(reference) <- "pathology"
# # use only the cells from the same tiddue
# reference <- subset(reference, idents = "c_chronic_active")
# 
# DimPlot(reference,reduction = "curated")
# # add some missing general annotations
# meta <- reference@meta.data %>%
#   rownames_to_column(var = "barcodes") %>%
#   # add the macroclassification
#   mutate(clusterCellType = case_when(seurat_clusters %in% c(0,1,2,3,6) ~"OLIGO",
#                                      seurat_clusters %in% c(7,15) ~"NEU",
#                                      seurat_clusters %in% c(8) ~"OPC",
#                                      seurat_clusters %in% c(11,13) ~"VAS",
#                                      seurat_clusters %in% c(16) ~"LYM",
#                                      seurat_clusters %in% c(5,10,17) ~"IMM",
#                                      seurat_clusters %in% c(4,9,12,14) ~"AST"))
# 
# reference$clusterCellType <- meta$clusterCellType
# # confirm the addition fo the data
# head(reference@meta.data)
# 
# # note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k cells this speeds up SCTransform dramatically with no loss in performance
# # library(dplyr)
# reference <- SCTransform(reference, ncells = 3000, verbose = T) %>%
#   RunPCA(verbose = T) %>%
#   RunUMAP(dims = 1:30)
# # save the sctransformed reference
# saveRDS(reference,"../out_large/CCA_SCTransformed.rds")

# -------------------------------------------------------------------------
# After subsetting, we renormalize cortex
# brain <- SCTransform(brain, assay = "Spatial", verbose = T) %>%
#   RunPCA(verbose = T)

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(reference, group.by = "clusterCellType", label = TRUE)

anchors <- FindTransferAnchors(reference = reference, query = brain, normalization.method = "SCT")

predictions.assay <- TransferData(anchorset = anchors, refdata = reference$clusterCellType, prediction.assay = TRUE,
                                  weight.reduction = brain[["pca"]], dims = 1:30)
brain[["predictions"]] <- predictions.assay
# Now we get prediction scores for each spot for each class. Of particular interest in the frontal brain region are the laminar excitatory neurons. Here we can distinguish between distinct sequential layers of these neuronal subtypes, for example:

brain

DefaultAssay(brain) <- "predictions"
SpatialFeaturePlot(brain, features = c("OLIGO", "OPC", "AST", "IMM", "NEU", "VAS", "LYM"), pt.size.factor = 1.6, ncol = 3, crop = T)
ggsave("out/image/07_CCA_label_transfer_SP01.pdf",width = 15,height = 15)

# Based on these prediction scores, we can also predict cell types whose location is spatially restricted. We use the same methods based on marked point processes to define spatially variable features, but use the cell type prediction scores as the “marks” rather than gene expression.

brain <- FindSpatiallyVariableFeatures(brain, assay = "predictions", selection.method = "markvariogram",
                                       features = c("OLIGO", "OPC", "AST", "IMM", "NEU", "VAS", "LYM"), r.metric = 5, slot = "data")

top.clusters <- head(SpatiallyVariableFeatures(brain), 4)

SpatialPlot(object = brain, features = top.clusters, ncol = 2)

# Finally, we show that our integrative procedure is capable of recovering the known spatial localization patterns of both neuronal and non-neuronal subsets, including laminar excitatory, layer-1 astrocytes, and the cortical grey matter.

SpatialFeaturePlot(brain, features = c("OLIGO", "OPC", "AST", "IMM", "NEU", "VAS", "LYM"),
                   pt.size.factor = 1, ncol = 3, crop = FALSE, alpha = c(0.1, 1))
