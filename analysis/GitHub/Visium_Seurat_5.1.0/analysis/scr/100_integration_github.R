# REF ---------------------------------------------------------------------
# reference for the processing
# https://github.com/satijalab/seurat/issues/7094
list_brain <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")
data.list.id <- names(list_brain)
data.combined.all <- merge(list_brain[[1]], y = list_brain[-1], add.cell.ids = data.list.id, project = "visium_brain")

# Run PCA again on all samples in group
data.combined.all <- SCTransform(data.combined.all, assay = "Spatial", verbose = TRUE) %>%
  RunPCA(assay = "SCT", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(verbose = FALSE,resolution = 0.1) %>%
  RunUMAP(reduction = "pca", dims = 1:30)

DimPlot(data.combined.all,
        group.by = "seurat_clusters",split.by = "dataset",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 3)

# Integration -------------------------------------------------------------
p.features <- SelectIntegrationFeatures(list_brain, nfeatures = 3000, verbose = FALSE)

#Prep for integration
p.prepsct <- PrepSCTIntegration(object.list = list_brain,
                                anchor.features = p.features, verbose = FALSE)

# Find anchors
p.int.anchors <- FindIntegrationAnchors(object.list = p.prepsct, normalization.method = "SCT", verbose = FALSE,
                                        anchor.features = p.features)
# Integrate
parental.integrated <- IntegrateData(anchorset = p.int.anchors, normalization.method = "SCT", verbose = FALSE)

parental.integrated <- RunPCA(parental.integrated, verbose = FALSE)
parental.integrated <- FindNeighbors(parental.integrated, dims = 1:30)
parental.integrated <- FindClusters(parental.integrated, verbose = FALSE)
parental.integrated <- RunUMAP(parental.integrated, dims = 1:30)