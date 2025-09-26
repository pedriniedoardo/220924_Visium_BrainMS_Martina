# libraries ---------------------------------------------------------------
# Load libraries
library(tidyverse)
library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(BayesSpace)
library(scales)

# read in the seurat object already processed -----------------------------
list_brain <- readRDS(file = "out/object/list_brain_all_BayesSpace1000_q10.rds")

# set the genes of interest
rownames(list_brain[[1]])[rownames(list_brain[[1]]) %in% c("SERPINA3")]

# try enhanced resolution -------------------------------------------------
# run all the slices
# brain <- list_brain$V01
list_enhanced <- pmap(list(list_brain,names(list_brain)),function(brain,nm){
  # check the progression
  print(nm)
  
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
  str_subset(markers,pattern = "SERPINA3")
  str_subset(rownames(brain),pattern = "SERPINA3")
  # sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
  #                                 feature_names=markers,
  #                                 nrounds=0)
  
  # enhance only some markers of interest
  sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                  feature_names = c("SERPINA3"),
                                  nrounds=0)
  return(sce.enhanced)
})

# featurePlot(sce.enhanced, "TSPO")
# 
# # save the sample input and enhanced object
saveRDS(list_enhanced,"out/object/list_enhanced_ALL_SERPINA3.rds")

list_sl_subset <- list_enhanced[c("V01","V03","V05","V06","V12")]

# plot SERPINA3
pdf("out/image/list_SERPINA3_enhanced_subset.pdf")
pmap(list(list_sl_subset,names(list_sl_subset)),function(x,name){
  # enhance only some markers of interest
  # sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
  #                                 feature_names=c("TSPO"),
  #                                 nrounds=0)
  featurePlot(x, "SERPINA3")+ggtitle(name)
})
dev.off()

# plot the scale above ----------------------------------------------------
#
library(patchwork)
pdf("out/image/list_SERPINA3_panel_enhanced.pdf",width = 6,height = 5)
pmap(list(list_sl_subset,names(list_sl_subset)),function(x,name){
  
    # list_plot <- lapply(c("SPP1","C3","HMGB1","IGFBP5","CDKN1A","CDKN1B"),function(gene){
    list_plot <- lapply(c("SERPINA3"),function(gene){
      featurePlot(x,gene)+theme(legend.position = "top")
    })
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
  
})
dev.off()
