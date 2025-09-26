# libraries ---------------------------------------------------------------
#Load libraries
library(tidyverse)
library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(BayesSpace)
library(scales)

# read in the seurat object already processed -----------------------------
list_brain <- readRDS(file = "out/object/list_brain_all_BayesSpace1000_q05.rds")


# plot --------------------------------------------------------------------
pdf("out/image/00_list_slide.pdf")
pmap(list(list_brain,names(list_brain)),function(x,name){
  # enhance only some markers of interest
  # sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
  #                                 feature_names=c("TSPO"),
  #                                 nrounds=0)
  SpatialFeaturePlot(x, features = "nCount_Spatial",alpha=0,crop = F) + theme(legend.position = "none") + ggtitle(name)
})
dev.off()
