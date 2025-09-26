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

# processing full dataset -------------------------------------------------
# in this case I am loading the dataset priduced by Martina
reference <- readRDS("../data/all20_integrated_clean_metadata.rds")
Idents(reference) <- "pathology"
# use only the cells from the same tiddue
reference <- subset(reference, idents = "c_chronic_active")

DimPlot(reference,reduction = "curated")
# add some missing general annotations
meta <- reference@meta.data %>%
  rownames_to_column(var = "barcodes") %>%
  # add the macroclassification
  mutate(clusterCellType = case_when(seurat_clusters %in% c(0,1,2,3,6) ~"OLIGO",
                                     seurat_clusters %in% c(7,15) ~"NEU",
                                     seurat_clusters %in% c(8) ~"OPC",
                                     seurat_clusters %in% c(11,13) ~"VAS",
                                     seurat_clusters %in% c(16) ~"LYM",
                                     seurat_clusters %in% c(5,10,17) ~"IMM",
                                     seurat_clusters %in% c(4,9,12,14) ~"AST"))

reference$clusterCellType <- meta$clusterCellType
# confirm the addition fo the data
head(reference@meta.data)

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k cells this speeds up SCTransform dramatically with no loss in performance
# library(dplyr)
reference <- SCTransform(reference, ncells = 3000, verbose = T) %>%
  RunPCA(verbose = T) %>%
  RunUMAP(dims = 1:30)
# save the sctransformed reference
saveRDS(reference,"../out_large/CCA_SCTransformed.rds")

# processing IMM subset ---------------------------------------------------
# in this case I am loading the dataset priduced by Martina
reference_IMM <- readRDS("../data/all20_immune.rds")
# Idents(reference) <- "pathology"
# # use only the cells from the same tiddue
# reference <- subset(reference, idents = "c_chronic_active")

DimPlot(reference_IMM,reduction = "umap",label = T)

# add some missing general annotations
meta_IMM <- reference_IMM@meta.data %>%
  rownames_to_column(var = "barcodes") %>%
  # add the macroclassification
  mutate(clusterCellType = case_when(seurat_clusters %in% c(0) ~"MG Homeo",
                                     seurat_clusters %in% c(1) ~"MIMS-foamy",
                                     seurat_clusters %in% c(2,5,6) ~"MG other",
                                     seurat_clusters %in% c(3) ~"MG Stress",
                                     seurat_clusters %in% c(4) ~"IMM high Mye",
                                     seurat_clusters %in% c(7) ~"T cells",
                                     seurat_clusters %in% c(8) ~"MIMS-iron",
                                     seurat_clusters %in% c(9) ~"MONO",
                                     seurat_clusters %in% c(10) ~"B Plasma"
                                     ))

reference_IMM$clusterCellType <- meta_IMM$clusterCellType
# confirm the addition fo the data
head(reference_IMM@meta.data)

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k cells this speeds up SCTransform dramatically with no loss in performance
# library(dplyr)
reference_IMM <- SCTransform(reference_IMM, ncells = 3000, verbose = T) %>%
  RunPCA(verbose = T) %>%
  RunUMAP(dims = 1:30)
# save the sctransformed reference
saveRDS(reference_IMM,"../out_large/IMM_ALL_SCTransformed.rds")

Idents(reference_IMM) <- "clusterCellType"
DimPlot(reference_IMM,reduction = "umap",label = T)

# processing ASTRO subset -------------------------------------------------
# in this case I am loading the dataset priduced by Martina
reference_ASTRO <- readRDS("../data/all20_astro.rds")
# Idents(reference) <- "pathology"
# # use only the cells from the same tiddue
# reference <- subset(reference, idents = "c_chronic_active")

DimPlot(reference_ASTRO,reduction = "umap",label = T)

# add some missing general annotations
meta_ASTRO <- reference_ASTRO@meta.data %>%
  rownames_to_column(var = "barcodes") %>%
  # add the macroclassification
  mutate(clusterCellType = case_when(seurat_clusters %in% c(0) ~"ASTRO nr",
                                     seurat_clusters %in% c(1,5) ~"ASTRO reactive",
                                     seurat_clusters %in% c(2,3,8) ~"ASTRO other",
                                     seurat_clusters %in% c(4) ~"ASTRO sen",
                                     seurat_clusters %in% c(6) ~"AIMS",
                                     seurat_clusters %in% c(7) ~"ASTRO peri"
  ))

reference_ASTRO$clusterCellType <- meta_ASTRO$clusterCellType
# confirm the addition fo the data
head(reference_ASTRO@meta.data)

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k cells this speeds up SCTransform dramatically with no loss in performance
# library(dplyr)
reference_ASTRO <- SCTransform(reference_ASTRO, ncells = 3000, verbose = T) %>%
  RunPCA(verbose = T) %>%
  RunUMAP(dims = 1:30)
# save the sctransformed reference
saveRDS(reference_ASTRO,"../out_large/ASTRO_ALL_SCTransformed.rds")

Idents(reference_ASTRO) <- "clusterCellType"
DimPlot(reference_ASTRO,reduction = "umap",label = T)
