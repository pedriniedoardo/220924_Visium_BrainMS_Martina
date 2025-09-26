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
# saveRDS(sce,"out/object/test_sce_V05.rds")
# 
# sce.enhanced <- readRDS("out/object/test_sce.enhanced_V05.rds")
# sce <- readRDS("out/object/test_sce_V05.rds")

# list_enhanced <- readRDS("out/object/list_enhanced_ALL.rds")

pdf("out/image/list_TSPO_enhanced.pdf")
pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  # enhance only some markers of interest
  # sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
  #                                 feature_names=c("TSPO"),
  #                                 nrounds=0)
  featurePlot(x, "TSPO")+ggtitle(name)
})
dev.off()

list_plot <- lapply(c("SPP1","C3","HMGB1","IGFBP5","FGF1","CDKN1A"),function(gene){
  featurePlot(list_enhanced$V16,gene)+theme(legend.position = "top")
})
wrap_plots(list_plot,nrow = 1)

library(patchwork)
pdf("out/image/list_panel_enhanced.pdf",width = 30,height = 5)
pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  if(name=="V08"){
    list_plot <- lapply(c("SPP1","HMGB1","IGFBP5","FGF1","CDKN1B"),function(gene){
      featurePlot(x,gene)+theme(legend.position = "top")
    })
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
  } else if(name=="V14"){
    list_plot <- lapply(c("SPP1","C3","HMGB1","IGFBP5","FGF1","CDKN1B"),function(gene){
      featurePlot(x,gene)+theme(legend.position = "top")
    })
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
  } else{
    # list_plot <- lapply(c("SPP1","C3","HMGB1","IGFBP5","CDKN1A","CDKN1B"),function(gene){
    list_plot <- lapply(c("SPP1","C3","HMGB1","IGFBP5","FGF1","CDKN1A"),function(gene){
      featurePlot(x,gene)+theme(legend.position = "top")
    })
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
  }
  
})
dev.off()

# cap the values on the whole panel using CA as cut off -------------------
# use the following color
# "#F03B20"

p1 <- featurePlot(list_enhanced$V01, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                                 limits = c(0,2),
                                                                 oob = scales::squish)+
  theme(legend.position = "top")

p2 <- featurePlot(list_enhanced$V01, "C3")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                               limits = c(0,1.4),
                                                               oob = scales::squish)+
  theme(legend.position = "top")

p3 <- featurePlot(list_enhanced$V01, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                                  limits = c(0,0.6),
                                                                  oob = scales::squish)+
  theme(legend.position = "top")

p4 <- featurePlot(list_enhanced$V01, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                                   limits = c(0,0.6),
                                                                   oob = scales::squish)+
  theme(legend.position = "top")

p5 <- featurePlot(list_enhanced$V01, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                                 limits = c(0,0.4),
                                                                 oob = scales::squish)+
  theme(legend.position = "top")

p6 <- featurePlot(list_enhanced$V01, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                                   limits = c(0,0.6),
                                                                   oob = scales::squish)+
  theme(legend.position = "top")

list_plot<-list(p1,p2,p3,p4,p5,p6)
wrap_plots(list_plot,nrow = 1) + plot_annotation(name)

pdf("out/image/list_panel_enhanced1.pdf",width = 30,height = 5)
pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  if(name=="V08"){
    # c("SPP1","HMGB1","IGFBP5","FGF1","CDKN1B")
    p1 <- featurePlot(x, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                     limits = c(0,2),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    # p2 <- featurePlot(x, "C3")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                                limits = c(0,1.4),
    #                                                                oob = scales::squish)+
    #   theme(legend.position = "top")
    
    p3 <- featurePlot(x, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                      limits = c(0,0.6),
                                                      oob = scales::squish)+
      theme(legend.position = "top")
    
    p4 <- featurePlot(x, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    p5 <- featurePlot(x, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                     limits = c(0,0.4),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p6 <- featurePlot(x, "CDKN1B")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    # p6 <- featurePlot(x, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                    limits = c(0,0.6),
    #                                                    oob = scales::squish)+
    #   theme(legend.position = "top")
    
    list_plot<-list(p1,p3,p4,p5,p6)
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
    
  } else if(name=="V14"){
    # c("SPP1","C3","HMGB1","IGFBP5","FGF1","CDKN1B")
    p1 <- featurePlot(x, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                     limits = c(0,2),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p2 <- featurePlot(x, "C3")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                   limits = c(0,1.4),
                                                   oob = scales::squish)+
      theme(legend.position = "top")
    
    p3 <- featurePlot(x, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                      limits = c(0,0.6),
                                                      oob = scales::squish)+
      theme(legend.position = "top")
    
    p4 <- featurePlot(x, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    p5 <- featurePlot(x, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                     limits = c(0,0.4),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p6 <- featurePlot(x, "CDKN1B")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    # p6 <- featurePlot(x, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                    limits = c(0,0.6),
    #                                                    oob = scales::squish)+
    #   theme(legend.position = "top")
    
    list_plot<-list(p1,p2,p3,p4,p5,p6)
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
    
  } else{
    # list_plot <- lapply(c("SPP1","C3","HMGB1","IGFBP5","CDKN1A","CDKN1B"),function(gene){
    # c("SPP1","C3","HMGB1","IGFBP5","FGF1","CDKN1A")
    p1 <- featurePlot(x, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                     limits = c(0,2),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p2 <- featurePlot(x, "C3")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                   limits = c(0,1.4),
                                                   oob = scales::squish)+
      theme(legend.position = "top")
    
    p3 <- featurePlot(x, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                      limits = c(0,0.6),
                                                      oob = scales::squish)+
      theme(legend.position = "top")
    
    p4 <- featurePlot(x, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    p5 <- featurePlot(x, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                     limits = c(0,0.4),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    # p6 <- featurePlot(x, "CDKN1B")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                    limits = c(0,0.6),
    #                                                    oob = scales::squish)+
    #   theme(legend.position = "top")
    p6 <- featurePlot(x, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    list_plot<-list(p1,p2,p3,p4,p5,p6)
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
  }
  
})
dev.off()

pdf("out/image/list_panel_enhanced2.pdf",width = 30,height = 5)
pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  if(name=="V08"){
    # c("SPP1","HMGB1","IGFBP5","FGF1","CDKN1B")
    p1 <- featurePlot(x, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                     limits = c(0,2),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    # p2 <- featurePlot(x, "C3")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                                limits = c(0,1.4),
    #                                                                oob = scales::squish)+
    #   theme(legend.position = "top")
    
    p3 <- featurePlot(x, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                      limits = c(0,0.6),
                                                      oob = scales::squish)+
      theme(legend.position = "top")
    
    p4 <- featurePlot(x, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    p5 <- featurePlot(x, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                     limits = c(0,0.4),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p6 <- featurePlot(x, "CDKN1B")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    # p6 <- featurePlot(x, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                    limits = c(0,0.6),
    #                                                    oob = scales::squish)+
    #   theme(legend.position = "top")
    
    list_plot<-list(p1,p3,p4,p5,p6)
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
    
  } else if(name=="V14"){
    # c("SPP1","C3","HMGB1","IGFBP5","FGF1","CDKN1B")
    p1 <- featurePlot(x, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                     limits = c(0,2),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p2 <- featurePlot(x, "C3")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                   limits = c(0,1.4),
                                                   oob = scales::squish)+
      theme(legend.position = "top")
    
    p3 <- featurePlot(x, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                      limits = c(0,0.6),
                                                      oob = scales::squish)+
      theme(legend.position = "top")
    
    p4 <- featurePlot(x, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    p5 <- featurePlot(x, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                     limits = c(0,0.4),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p6 <- featurePlot(x, "CDKN1B")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    # p6 <- featurePlot(x, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                    limits = c(0,0.6),
    #                                                    oob = scales::squish)+
    #   theme(legend.position = "top")
    
    list_plot<-list(p1,p2,p3,p4,p5,p6)
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
    
  } else{
    # list_plot <- lapply(c("SPP1","C3","HMGB1","IGFBP5","CDKN1A","CDKN1B"),function(gene){
    # c("SPP1","C3","HMGB1","IGFBP5","FGF1","CDKN1A")
    p1 <- featurePlot(x, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                     limits = c(0,2),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p2 <- featurePlot(x, "C3")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                   limits = c(0,1.4),
                                                   oob = scales::squish)+
      theme(legend.position = "top")
    
    p3 <- featurePlot(x, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                      limits = c(0,0.6),
                                                      oob = scales::squish)+
      theme(legend.position = "top")
    
    p4 <- featurePlot(x, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    p5 <- featurePlot(x, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                     limits = c(0,0.4),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    # p6 <- featurePlot(x, "CDKN1B")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                    limits = c(0,0.6),
    #                                                    oob = scales::squish)+
    #   theme(legend.position = "top")
    p6 <- featurePlot(x, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = "black",
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    list_plot<-list(p1,p2,p3,p4,p5,p6)
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
  }
  
})
dev.off()

pdf("out/image/list_panel_enhanced3.pdf",width = 30,height = 5)
pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  if(name=="V08"){
    # c("SPP1","HMGB1","IGFBP5","FGF1","CDKN1B")
    p1 <- featurePlot(x, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                     limits = c(0,2),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    # p2 <- featurePlot(x, "C3")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                                limits = c(0,1.4),
    #                                                                oob = scales::squish)+
    #   theme(legend.position = "top")
    
    p3 <- featurePlot(x, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                      limits = c(0,0.6),
                                                      oob = scales::squish)+
      theme(legend.position = "top")
    
    p4 <- featurePlot(x, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    p5 <- featurePlot(x, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                     limits = c(0,0.4),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p6 <- featurePlot(x, "CDKN1B")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    # p6 <- featurePlot(x, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                    limits = c(0,0.6),
    #                                                    oob = scales::squish)+
    #   theme(legend.position = "top")
    
    list_plot<-list(p1,p3,p4,p5,p6)
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
    
  } else if(name=="V14"){
    # c("SPP1","C3","HMGB1","IGFBP5","FGF1","CDKN1B")
    p1 <- featurePlot(x, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                     limits = c(0,2),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p2 <- featurePlot(x, "C3")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                   limits = c(0,1.4),
                                                   oob = scales::squish)+
      theme(legend.position = "top")
    
    p3 <- featurePlot(x, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                      limits = c(0,0.6),
                                                      oob = scales::squish)+
      theme(legend.position = "top")
    
    p4 <- featurePlot(x, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    p5 <- featurePlot(x, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                     limits = c(0,0.4),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p6 <- featurePlot(x, "CDKN1B")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    # p6 <- featurePlot(x, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                    limits = c(0,0.6),
    #                                                    oob = scales::squish)+
    #   theme(legend.position = "top")
    
    list_plot<-list(p1,p2,p3,p4,p5,p6)
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
    
  } else{
    # list_plot <- lapply(c("SPP1","C3","HMGB1","IGFBP5","CDKN1A","CDKN1B"),function(gene){
    # c("SPP1","C3","HMGB1","IGFBP5","FGF1","CDKN1A")
    p1 <- featurePlot(x, "SPP1")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                     limits = c(0,2),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    p2 <- featurePlot(x, "C3")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                   limits = c(0,1.4),
                                                   oob = scales::squish)+
      theme(legend.position = "top")
    
    p3 <- featurePlot(x, "HMGB1")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                      limits = c(0,0.6),
                                                      oob = scales::squish)+
      theme(legend.position = "top")
    
    p4 <- featurePlot(x, "IGFBP5")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    p5 <- featurePlot(x, "FGF1")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                     limits = c(0,0.4),
                                                     oob = scales::squish)+
      theme(legend.position = "top")
    
    # p6 <- featurePlot(x, "CDKN1B")+scale_fill_gradient(low = "#F0F0F0",high = "#F03B20",
    #                                                    limits = c(0,0.6),
    #                                                    oob = scales::squish)+
    #   theme(legend.position = "top")
    p6 <- featurePlot(x, "CDKN1A")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                                       limits = c(0,0.6),
                                                       oob = scales::squish)+
      theme(legend.position = "top")
    
    list_plot<-list(p1,p2,p3,p4,p5,p6)
    wrap_plots(list_plot,nrow = 1) + plot_annotation(name)
  }
  
})
dev.off()

# library(patchwork)
# pdf("out/image/test.pdf",width = 15,height = 10)
# pmap(list(list_enhanced[-12],names(list_enhanced[-12])),function(x,name){
#     list_plot <- lapply(c("SPP1","C3","HMGB1","IGFBP5","CDKN1A","CDKN1B"),function(gene){
#       featurePlot(x,gene)
#     })
#     wrap_plots(list_plot) + plot_annotation(name)
# })
# dev.off()



# test --------------------------------------------------------------------
#
featurePlot(list_enhanced$V05, "TSPO")
featurePlot(list_enhanced$V05, "TSPO")+scale_fill_gradient(low = "#F0F0F0",high = muted("black"),
                                                           limits = c(0,0.3),
                                                           oob = scales::squish)+ggtitle("V05")

pdf("out/image/list_TSPO_enhanced_cap04.pdf")
pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  featurePlot(x, "TSPO")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                             limits = c(0,0.4),
                                             oob = scales::squish) +
    ggtitle(name)
})
dev.off()

pdf("out/image/list_TSPO_enhanced_cap03.pdf")
pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  featurePlot(x, "TSPO")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                             limits = c(0,0.3),
                                             oob = scales::squish) +
    ggtitle(name)
})
dev.off()

pdf("out/image/list_TSPO_enhanced_cap02.pdf")
pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  featurePlot(x, "TSPO")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                             limits = c(0,0.2),
                                             oob = scales::squish) +
    ggtitle(name)
})
dev.off()

pdf("out/image/list_TSPO_enhanced_cap01.pdf")
pmap(list(list_enhanced,names(list_enhanced)),function(x,name){
  featurePlot(x, "TSPO")+scale_fill_gradient(low = "#F0F0F0",high = muted("red"),
                                             limits = c(0,0.1),
                                             oob = scales::squish) +
    ggtitle(name)
})
dev.off()
