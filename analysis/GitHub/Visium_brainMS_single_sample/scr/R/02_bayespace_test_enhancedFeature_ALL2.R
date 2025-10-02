# libraries ---------------------------------------------------------------
#Load libraries
library(tidyverse)
library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(BayesSpace)
library(scales)
library(patchwork)

# read in the seurat object already processed -----------------------------
list_brain <- readRDS(file = "out/object/list_brain_all_BayesSpace1000_q05.rds")

# set the genes of interest
rownames(list_brain[[1]])[rownames(list_brain[[1]]) %in% c("SPP1","C3","HMGB1","IGFBP5","CHIT1","TSPO")]

# try enhanced resolution -------------------------------------------------
# run all the slices
# brain <- list_brain$V01
list_enhanced <- map(list_brain,function(brain){
  # print()
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
  # str_subset(markers,pattern = "TSPO")
  # str_subset(rownames(brain),pattern = "TSPO")
  # sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
  #                                 feature_names=markers,
  #                                 nrounds=0)
  
  # enhance only some markers of interest
  sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
                                  feature_names = c("SPP1","C3","HMGB1","IGFBP5","CHIT1","TSPO"),
                                  nrounds=0)
  return(sce.enhanced)
})

# featurePlot(sce.enhanced, "TSPO")
# 
# # save the sample input and enhanced object
saveRDS(list_enhanced,"out/object/list_enhanced_ALL2.rds")

# cap the values for all the slices ---------------------------------------

# plot CHIT1 from all the slices
# x <- list_enhanced$V01
# name <- "V01"

list_plot <- pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  # track the progresso of the plotting
  print(name)
  
  # if the gene is missiong from the matrix, add it as NA
  if("CHIT1" %in% rownames(x)){
    p <- featurePlot(x, "CHIT1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                     limits = c(0,2),
                                                     oob = scales::squish) +
      theme(legend.position = "top") +
      labs(title = name)
    
    return(p)
    
  } else {
    # remove the dim reduction from the object to allow the concatenation
    reducedDims(x) <- NULL
    # reducedDims(new_gene_sce)
    
    # 1. Create a vector of raw counts for the new gene across all cells
    # new_gene_counts <- rpois(n = ncol(sce), lambda = 5)
    new_gene_counts <- rep(NA,ncol(x))
    
    # 2. Create a new, single-gene SingleCellExperiment object
    new_gene_sce <- SingleCellExperiment(
      assays = list(logcounts = matrix(new_gene_counts, nrow = 1,
                                       dimnames = list("CHIT1", colnames(x))))
    )
    
    # 3. Combine the original and new objects
    combined_sce <- rbind(x, new_gene_sce)
    
    # produce the plot
    p <- featurePlot(combined_sce, "CHIT1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                     limits = c(0,2),
                                                     oob = scales::squish) +
      theme(legend.position = "top") +
      labs(title = name)
    
    return(p)
  }
  
})

pdf("out/image/list_panel_enhanced_CHIT1.pdf",width = 20,height = 20)
wrap_plots(list_plot)
dev.off()

# plot the data using 0 instead of NA
list_plot2 <- pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  # track the progresso of the plotting
  print(name)
  
  # if the gene is missiong from the matrix, add it as NA
  if("CHIT1" %in% rownames(x)){
    p <- featurePlot(x, "CHIT1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                     limits = c(0,1.5),
                                                     oob = scales::squish) +
      theme(legend.position = "top") +
      labs(title = name)
    
    return(p)
    
  } else {
    # remove the dim reduction from the object to allow the concatenation
    reducedDims(x) <- NULL
    # reducedDims(new_gene_sce)
    
    # 1. Create a vector of raw counts for the new gene across all cells
    # new_gene_counts <- rpois(n = ncol(sce), lambda = 5)
    new_gene_counts <- rep(0,ncol(x))
    
    # 2. Create a new, single-gene SingleCellExperiment object
    new_gene_sce <- SingleCellExperiment(
      assays = list(logcounts = matrix(new_gene_counts, nrow = 1,
                                       dimnames = list("CHIT1", colnames(x))))
    )
    
    # 3. Combine the original and new objects
    combined_sce <- rbind(x, new_gene_sce)
    
    # produce the plot
    p <- featurePlot(combined_sce, "CHIT1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                                limits = c(0,1.5),
                                                                oob = scales::squish) +
      theme(legend.position = "top") +
      labs(title = name)
    
    return(p)
  }
  
})

pdf("out/image/list_panel_enhanced_CHIT1_2.pdf",width = 20,height = 20)
wrap_plots(list_plot2)
dev.off()
