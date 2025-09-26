# AIM ---------------------------------------------------------------------
# try to run the NICHES pipeline on the all the slides as test

# libraries ---------------------------------------------------------------
# First, let's load dependencies.
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(NICHES)
library(viridis)
library(tidyverse)


# read in the data --------------------------------------------------------
# read the data from the previous analysis
list_brain <- readRDS("../../Visium_brainMS_single_sample/out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")
list_brain_bayespace <- readRDS(file = "../../Visium_brainMS_single_sample/out/object/list_brain_all_BayesSpace1000_q10.rds")

LUT_sample <- read_csv("../../Visium_brainMS_single_sample/data/sample_classification.csv") %>%
  mutate(dataset = paste0("brain_",sample))

# add the BayesSpace classification to all the datasets in the list
# sl <- list_brain$V01 
# sl_bayes <- list_brain_bayespace$V01 

list_brain_full <- pmap(list(list_brain,list_brain_bayespace[names(list_brain)],names(list_brain)), function(sl,sl_bayes,name){
  print(name)
  meta_full <- sl@meta.data %>%
    rownames_to_column("barcode") %>%
    left_join(sl_bayes@meta.data %>%
                rownames_to_column("barcode") %>%
                select(barcode,BayesSpace),by = "barcode") %>%
    column_to_rownames("barcode")
  
  # update the meta
  sl@meta.data <- meta_full
  
  return(sl)
})

# read in the LR database. this is the same one used for the single cell CCC analysis.
LR_database_list <- readRDS(file = "../data/CellChatDB.human.used.rds")
LR_database_df <- LR_database_list$interaction %>%
  select(ligand.symbol,receptor.symbol) %>%
  # modify the comma space with _
  mutate(ligand.symbol = str_replace_all(ligand.symbol,", ","_"),
         receptor.symbol = str_replace_all(receptor.symbol,", ","_"))

# read in the list of significant LR pairs detected in the single cell dataset
LR_sig <- read_tsv("../data/135_net.up_cellID_SEN_AST.tsv")

# load the default LR database
# db_human <- NICHES::LoadFantom5(species = "human")

# str(db_mouse)

# custom_LR_database
# data.frame. Each row is a ligand-receptor mechanism where the first column corresponds to the source genes that express the ligands subunits (separated by '_') and the second column corresponds to the receptor genes that express the receptor subunits (separated by '_').
# db_human_df <- data.frame(source.subunits = db_human$source.subunits[,1],
#                           target.subunits = db_human$target.subunits[,1])

# shape the cellchat database as the one in the NICHES
LR_database_df_fix <- LR_database_df %>%
  separate_rows(ligand.symbol,sep = "_") %>%
  separate_rows(receptor.symbol,sep = "_") %>%
  # make the pairs unique
  mutate(pair = paste(ligand.symbol,receptor.symbol,sep = "_")) %>%
  group_by(pair) %>%
  summarise() %>%
  separate(pair,into = c("ligand.symbol","receptor.symbol"),sep = "_")

# make sure the LR pairs of significance are present
LR_sig %>%
  select(ligand,receptor) %>%
  separate_rows(ligand,sep = "_") %>%
  separate_rows(receptor,sep = "_") %>%
  # make the pairs unique
  mutate(pair = paste(ligand,receptor,sep = "_")) %>%
  group_by(pair) %>%
  summarise() %>%
  separate(pair,into = c("ligand.symbol","receptor.symbol"),sep = "_") %>%
  left_join(LR_database_df_fix %>% mutate(test = 1),by = c("ligand.symbol","receptor.symbol")) %>%
  filter(is.na(test))

# one is missing
# ligand.symbol receptor.symbol  test
# PGD2-PTGDS    PTGDR              NA

# # wrangling ---------------------------------------------------------------
# # x <- list_brain_full$V01
# # nm <- "V01"
# list_brain_NICHES <- pmap(list(list_brain_full,names(list_brain_full)),function(x,nm){
#   # check the progress
#   print(nm)
#   
#   # subset the slide v01
#   # test <- list_brain_full$V01
#   test <- x
#   test <- UpdateSeuratObject(test) # JC: need to update seurat obj
#   
#   # Next, we will format the spatial coordinate metadata so that every cell has an explicitly labeled x and y coordinate.
#   
#   ## Format Spatial Coordinates and Normalize
#   test@meta.data$x <- test@images$slice1@coordinates$row
#   test@meta.data$y <- test@images$slice1@coordinates$col
#   
#   DefaultAssay(test) <- "Spatial"
#   
#   # NICHES can be run on imputed or non-imputed data. Here, we will use imputed data.
#   
#   ## Impute and Run NICHES
#   # the process below seems to be randomic, therefore to make the process reproducible, set a seed
#   set.seed(2144)
#   test <- SeuratWrappers::RunALRA(test)
#   
#   NICHES_output <- RunNICHES(object = test,
#                              LR.database = "custom",
#                              custom_LR_database = LR_database_df_fix,
#                              # species = "human",
#                              assay = "alra",
#                              position.x = 'x',
#                              position.y = 'y',
#                              k = 8,
#                              cell_types = "BayesSpace",
#                              min.cells.per.ident = 0,
#                              min.cells.per.gene = NULL,
#                              meta.data.to.map = c('dataset','BayesSpace'),
#                              CellToCell = F,CellToSystem = F,SystemToCell = F,
#                              CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = T)
#   
#   # NICHES outputs a list of objects. Each object contains a certain style of cell-system signaling atlas. Above, we have only calculated a single one of interest, namely, individual cellular microenvironment. We next isolate this output and embed using UMAP to visualize the microenvironemnt of each cell.
#   niche <- NICHES_output[['NeighborhoodToCell']]
#   Idents(niche) <- niche@meta.data$ReceivingType
#   
#   # Scale and visualize
#   niche <- ScaleData(niche)
#   niche <- FindVariableFeatures(niche,selection.method = "disp")
#   niche <- RunPCA(niche)
#   # ElbowPlot(niche,ndims = 50)
#   niche <- RunUMAP(niche,dims = 1:10)
#   
#   # Further, and perhaps more usefully, we can map over the output from NICHES onto the original spatial object as follows:
#   
#   # Add Niches output as an assay
#   niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
#   colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
#   test[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
#   DefaultAssay(test) <- "NeighborhoodToCell"
#   test <- ScaleData(test)
#   
#   return(list(test = test,niche = niche))
#   
# })
# 
# saveRDS(list_brain_NICHES,"../out/object/list_brain_NICHES_all_135_net.up_cellID_SEN_AST.rds")
# 
# test <- readRDS("../out/object/list_brain_NICHES_all.rds")
