# libraries ---------------------------------------------------------------
#Load libraries
library(tidyverse)
library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(BayesSpace)
library(patchwork)

# read in data ------------------------------------------------------------
# read in the processed data of the brain slices
list_brain <- readRDS("out/object/list_brain_all.rds")

# run the bayespace and export the seurat object --------------------------
# define the slices
file_id <- names(list_brain)

list_brain_bayespace <- lapply(file_id,function(x){
  # keep track of the progress of the reading
  print(x)
  # read in the data seurat object
  brain <- test <- list_brain[[x]]
  # SpatialDimPlot(brain, "seurat_clusters")
  
  #Convert to SCE
  diet.seurat <- Seurat::DietSeurat(brain, graphs = "pca") #slim down Seurat obj prior to conversion
  sce <- as.SingleCellExperiment(diet.seurat) #convert seurat to SCE
  colData(sce) <- cbind(colData(sce), brain@images$slice1@coordinates) #add spatial info to SCE
  
  #BayesSpace Workflow
  set.seed(21)
  
  sce <- spatialPreprocess(sce, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
  # sce <- qTune(sce, qs=seq(2, 10), platform="Visium")
  # qPlot(sce)
  
  sce <- spatialCluster(sce, nrep = 1000, burn.in = 100, q = 5) #quickly cluster via BayesSpace
  # sce <- spatialCluster(sce, nrep = 1000, burn.in = 100,q = 4)
  # clusterPlot(sce) #plot via BayesSpace
  saveRDS(sce,file = paste0("out/object/sce_BayesSpaceCluster1000_",x,"_q05.rds"))
  
  #Add BayesSpace results to Seurat
  brain@meta.data = cbind(brain@meta.data, BayesSpace = sce$spatial.cluster) #add BayesSpace clusters to Seurat obj
  # SpatialDimPlot(brain, "BayesSpace") #plot via Seurat
  
  return(brain)
}) %>% 
  setNames(file_id)

# save the list of all the brain slices with the bayespace annotation
saveRDS(list_brain_bayespace,file = "out/object/list_brain_all_BayesSpace1000_q05.rds")

list_brain_bayespace <- readRDS(file = "out/object/list_brain_all_BayesSpace1000_q05.rds")

# plotting ----------------------------------------------------------------
lapply(file_id,function(x){
  # track the progress
  print(x)
  # read in the seurat object -----------------------------------------------
  test <- list_brain_bayespace[[x]]
  
  # read in the metadata ----------------------------------------------------
  # wrangling
  pdf(paste0("out/image/10_BayesSpace1000_q05_",x,".pdf"),width = 15,height = 15)
  
  print(SpatialDimPlot(test, "BayesSpace") +
          plot_annotation(title = "BayesSpace"))
  
  Idents(test) <- "BayesSpace"
  print(SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE, ncol = 4,alpha = 0.5) +
          plot_annotation(title = "BayesSpace"))
  
  print(SpatialDimPlot(test, "seurat_clusters") +
          plot_annotation(title = "seurat clusters"))
  
  Idents(test) <- "seurat_clusters"
  print(SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE, ncol = 3,alpha = 0.5) +
          plot_annotation(title = "seurat clusters"))
  
  dev.off()
  
})
