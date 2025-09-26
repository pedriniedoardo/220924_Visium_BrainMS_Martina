# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)

# read in data ------------------------------------------------------------
# read in the deconvolution result from each slice
folder2 <- "out/table/"
file_id2 <- dir(folder2) %>% 
  str_subset("SPOTLight") %>% 
  str_subset("IMM",negate = T) %>% 
  str_subset("ASTRO",negate = T) %>% 
  str_remove_all(pattern = "_SPOTLight.tsv")

# build a full matrix for the expression information on TSPO
list_meta_dec <- lapply(file_id2,function(x){
  
  # keep track of the progress of the reading
  print(x)
  
  df_meta <- read_tsv(paste0("out/table/",x,"_SPOTLight.tsv")) %>% 
    mutate(dataset = paste0("brain_",x))
  
  return(df_meta)
  
}) %>% 
  setNames(file_id2)

# read in the processed data of the brain slices
list_brain <- readRDS("out/object/list_brain_all.rds")

# add the spotlight annotation --------------------------------------------
list_brain_spotlight <- lapply(file_id2,function(x){
  # check the progress
  print(x)
  
  # read in the seurat object -----------------------------------------------
  test <- list_brain[[x]]
  
  # read in the metadata ----------------------------------------------------
  # wrangling pick the top class per score
  df_class_top <- list_meta_dec[[x]] %>% 
    dplyr::select(barcodes,AST:VAS) %>% 
    pivot_longer(names_to = "class_top",values_to = "score",-barcodes) %>% 
    group_by(barcodes) %>% 
    arrange(desc(score)) %>% 
    dplyr::slice(1) %>% 
    # top_n(wt = score,n = 1) %>% 
    dplyr::select(-score)
  
  # wrangling pick the class for spots with top 50  or top 70 representation
  df_class_prop <- list_meta_dec[[x]] %>% 
    mutate(class_50 = case_when(AST > 0.5~"AST",
                                IMM > 0.5~"IMM",
                                LYM > 0.5~"LYM",
                                NEU > 0.5~"NEU",
                                OLIGO > 0.5~"OLIGO",
                                OPC > 0.5~"OPC",
                                VAS > 0.5~"VAS",
                                T~"undefined"),
           class_70 = case_when(AST > 0.7~"AST",
                                IMM > 0.7~"IMM",
                                LYM > 0.7~"LYM",
                                NEU > 0.7~"NEU",
                                OLIGO > 0.7~"OLIGO",
                                OPC > 0.7~"OPC",
                                VAS > 0.7~"VAS",
                                T~"undefined")
    )
  
  # join the two metadata
  meta_fix <- left_join(df_class_prop,df_class_top,"barcodes")
  # meta_fix
  
  # add the meta from the spotlight
  # class(list_brain$V01@meta.data)
  test@meta.data <- meta_fix %>%
    column_to_rownames("barcodes")
  
  return(test)
}) %>% 
  setNames(file_id2)

# save a list of all the preprocessed brain slices with spotligth annotation
saveRDS(list_brain_spotlight,"out/object/list_brain_all_spotlight.rds")
list_brain_spotlight <- readRDS("out/object/list_brain_all_spotlight.rds")

# loop all the plotting ---------------------------------------------------
# x <- file_id2[1]
lapply(file_id2,function(x){
  # check the progress
  print(x)
  
  # read in the seurat object -----------------------------------------------
  test <- list_brain_spotlight[[x]]
  
  # make the different classes as categorical
  # SpatialFeaturePlot(test, features = "IMM") + theme(legend.position = "right")
  pdf(paste0("out/image/10_class_deconvolution_",x,".pdf"),width = 15,height = 15)
  Idents(test) <- "class_50"
  print(SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE, ncol = 3,alpha = 0.5) +
          plot_annotation(title = "class_50"))

  Idents(test) <- "class_70"
  print(SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE, ncol = 3,alpha = 0.5) +
          plot_annotation(title = "class_70"))

  Idents(test) <- "class_top"
  print(SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE, ncol = 3,alpha = 0.5) +
          plot_annotation(title = "class_top"))
  dev.off()
})
