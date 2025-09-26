# libraries ---------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(NICHES)
library(viridis)
library(tidyverse)
library(ggrepel)
library(pals)
library(scales)

# read in the data --------------------------------------------------------
list_brain_NICHES <- readRDS("../out/object/list_brain_NICHES_all.rds")
LR_sig <- read_tsv("../data/135_net.up_cellID_SEN_AST.tsv")

# plot the panel defined by Martina
LR_sig_pairs <- c("C3—ITGB2","ANGPTL4—SDC4")

LUT_sample <- read_csv("../../Visium_brainMS_single_sample/data/sample_classification.csv") %>%
  mutate(dataset = paste0("brain_",sample))

# set the color parameters of the scale
library(RColorBrewer)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

library(scales)
show_col(SpatialColors(10))

# wrangling ---------------------------------------------------------------
# select just a subset of the slides
list_brain_NICHES_subset <- list_brain_NICHES[c("V05","V01","V07","V03","V08")]

# plot --------------------------------------------------------------------
# fetch all the plotting data for the current panel of genes and slide.
# the aim is to define the top and bottom value across the whole set of slides

# for scaled values
range_slide_scaled <- pmap(list(list_brain_NICHES_subset,names(list_brain_NICHES_subset)), function(sl,nm){
  df <- FetchData(object = sl$test,vars = LR_sig_pairs,slot = 'scale.data')
  # df <- FetchData(object = sl$test,vars = LR_sig_pairs,slot = 'data')
  
  df %>%
    pivot_longer(names_to = "gene",values_to = "value",everything()) %>%
    group_by(gene) %>%
    summarise(min_val = min(value),
              max_val = max(value),.groups = "drop") %>%
    mutate(slide = nm)
  
}) %>%
  bind_rows()

# for normalized values
range_slide_norm <- pmap(list(list_brain_NICHES_subset,names(list_brain_NICHES_subset)), function(sl,nm){
  # df <- FetchData(object = sl$test,vars = LR_sig_pairs,slot = 'scale.data')
  df <- FetchData(object = sl$test,vars = LR_sig_pairs,slot = 'data')
  
  df %>%
    pivot_longer(names_to = "gene",values_to = "value",everything()) %>%
    group_by(gene) %>%
    summarise(min_val = min(value),
              max_val = max(value),.groups = "drop") %>%
    mutate(slide = nm)
  
}) %>%
  bind_rows()

# define the range for the genes. scaled
list_range_scaled <- range_slide_scaled %>%
  group_by(gene) %>%
  summarise(min_val = min(min_val),
            max_val = max(max_val),.groups = "drop") %>%
  split(.$gene)

# define the range for the genes. scaled
list_range_norm <- range_slide_norm %>%
  group_by(gene) %>%
  summarise(min_val = min(min_val),
            max_val = max(max_val),.groups = "drop") %>%
  split(.$gene)

# list_brain_NICHES_subset$V01$test

# plot using the scaled values
lapply(list_range_scaled, function(p){
  # pick the pair
  pair <- p$gene
  
  # keep track of the processing
  print(pair)
  
  # loop the plot across the slides
  list_plot <- lapply(list_brain_NICHES_subset, function(sl){
    
    # pull the ranges and feature name
    limit_min <- list_range_scaled[[pair]]$min_val
    limit_max <- list_range_scaled[[pair]]$max_val -3
    pair_id <- list_range_scaled[[pair]]$gene
    
    test_test <- sl$test
    spot_size <- test_test@images$slice1@scale.factors$fiducial/test_test@images$slice1@scale.factors$hires
    
    # plot <- SpatialFeaturePlot(test_test,
    #                            features = LR_sig_pairs[1],
    #                            slot = 'data',pt.size.factor = spot_size*1.5)
    
    plot <- SpatialFeaturePlot(test_test,
                               features = pair_id,
                               slot = 'scale.data',pt.size.factor = spot_size*1.5)
    
    p1 <- plot +
      scale_fill_gradientn(colours = SpatialColors(n = 100),
                           breaks = c(limit_min, round(limit_max - (limit_max-limit_min)/2),limit_max),
                           labels = c("Low", "Medium", "High"),
                           limits = c(limit_min, limit_max),
                           # this one allows to cap the value above a certain threshold
                           oob = scales::squish)+
      theme(legend.position = "top")
    
    return(p1)
  })
  
  p2 <- wrap_plots(list_plot,ncol = 5)
  ggsave(plot = p2,paste0("../out/plot/112_NICHES_",pair,"_scaled_panelSigLR_AST.pdf"),width = 20,height = 4)
})

# plot using the normalized values values
lapply(list_range_norm, function(p){
  # pick the pair
  pair <- p$gene
  
  # keep track of the processing
  print(pair)
  
  # loop the plot across the slides
  list_plot <- lapply(list_brain_NICHES_subset, function(sl){
    
    # pull the ranges and feature name
    limit_min <- list_range_norm[[pair]]$min_val
    limit_max <- list_range_norm[[pair]]$max_val -10
    pair_id <- list_range_norm[[pair]]$gene
    
    test_test <- sl$test
    spot_size <- test_test@images$slice1@scale.factors$fiducial/test_test@images$slice1@scale.factors$hires
    
    plot <- SpatialFeaturePlot(test_test,
                               features = pair_id,
                               slot = 'data',pt.size.factor = spot_size*1.5)
    
    # plot <- SpatialFeaturePlot(test_test,
    #                            features = pair_id,
    #                            slot = 'scale.data',pt.size.factor = spot_size*1.5)
    
    p1 <- plot +
      scale_fill_gradientn(colours = SpatialColors(n = 100),
                           breaks = c(limit_min, round(limit_max - (limit_max-limit_min)/2),limit_max),
                           labels = c("Low", "Medium", "High"),
                           limits = c(limit_min, limit_max),
                           # this one allows to cap the value above a certain threshold
                           oob = scales::squish)+
      theme(legend.position = "top")
    
    return(p1)
  })
  
  p2 <- wrap_plots(list_plot,ncol = 5)
  ggsave(plot = p2,paste0("../out/plot/112_NICHES_",pair,"_norm_panelSigLR_AST.pdf"),width = 20,height = 4)
})
