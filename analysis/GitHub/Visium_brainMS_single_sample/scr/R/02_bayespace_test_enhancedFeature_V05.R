# libraries ---------------------------------------------------------------
#Load libraries
library(tidyverse)
library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(BayesSpace)

# read in the seurat object already processed -----------------------------
list_brain <- readRDS(file = "out/object/list_brain_all_BayesSpace1000_q05.rds")

# try enhanced resolution -------------------------------------------------
# select just on one slice
brain <- list_brain$V05
# increase the number of HVG
# brain <- SCTransform(brain, assay = "Spatial", verbose = T,variable.features.n = 5000)
# 
# Convert to SCE
diet.seurat <- Seurat::DietSeurat(brain, graphs = "pca") #slim down Seurat obj prior to conversion
sce <- as.SingleCellExperiment(diet.seurat) #convert seurat to SCE
colData(sce) <- cbind(colData(sce), brain@images$slice1@coordinates) #add spatial info to SCE

# BayesSpace Workflow
sce <- spatialPreprocess(sce, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
# sce <- qTune(sce, qs=seq(2, 10), platform="Visium")
# qPlot(sce)

sce <- spatialCluster(sce, nrep = 1000, burn.in = 100, q = 10) #quickly cluster via BayesSpace

# enhance spots
sce.enhanced <- spatialEnhance(sce, q=10, platform="Visium",
                               nrep=2000, gamma=3, verbose=TRUE,
                               jitter_scale=5.5, jitter_prior=0.3,
                               save.chain=TRUE,burn.in = 100)
# 
# # enhance the expresison of the features focus on the HVG
markers <- brain@assays$SCT@var.features
str_subset(markers,pattern = "TSPO")
str_subset(rownames(brain),pattern = "TSPO")
# sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
#                                 feature_names=markers,
#                                 nrounds=0)

# enhance only some markers of interest
sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                feature_names=c("TSPO"),
                                nrounds=0)

featurePlot(sce.enhanced, "TSPO")

# save the sample input and enhanced object
saveRDS(sce.enhanced,"out/object/test_sce.enhanced_V05.rds")
saveRDS(sce,"out/object/test_sce_V05.rds")

sce.enhanced <- readRDS("out/object/test_sce.enhanced_V05.rds")
sce <- readRDS("out/object/test_sce_V05.rds")
