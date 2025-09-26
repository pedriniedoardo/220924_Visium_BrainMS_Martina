# AIM ---------------------------------------------------------------------
# try to run the NICHES pipeline on the V01 as test

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
LR_sig <- read_tsv("../data/128_net.up_cellID_SEN.tsv")

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

# wrangling ---------------------------------------------------------------
# subset the slide v01
test <- list_brain_full$V01
test <- UpdateSeuratObject(test) # JC: need to update seurat obj

# for Seurat v5.1.0 there is the need to update the size of the spots
spot_size <- test@images$slice1@scale.factors$fiducial/test@images$slice1@scale.factors$hires

# make sure the data are already normalized
test@assays$SCT@data
SpatialFeaturePlot(test, features = c("SPP1", "CD44"),pt.size.factor = spot_size,slot = "data")

p1 <- DimPlot(test, reduction = "umap",group.by = 'BayesSpace', label = TRUE)
p2 <- SpatialDimPlot(test, label = TRUE,group.by = 'BayesSpace', label.size = 3,pt.size.factor = spot_size)
Idents(test) <- 'BayesSpace'
p3 <- (SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE,pt.size.factor = spot_size)+plot_annotation('BayesSpace Clusters'))
wrap_plots(list(p1,p2,p3))

# Next, we will format the spatial coordinate metadata so that every cell has an explicitly labeled x and y coordinate.

## Format Spatial Coordinates and Normalize
test@meta.data$x <- test@images$slice1@coordinates$row
test@meta.data$y <- test@images$slice1@coordinates$col

DefaultAssay(test) <- "Spatial"
# make sure the slot is already normalized
# test <- NormalizeData(test)
test@assays$Spatial@data

# NICHES can be run on imputed or non-imputed data. Here, we will use imputed data.

## Impute and Run NICHES
# the process below seems to be randomic, therefore to make the process reproducible, set a seed
set.seed(2144)
test <- SeuratWrappers::RunALRA(test)

# sample routine for a custom LR database ---------------------------------
# is it possible to change the database?
# LR.database	string. Default: "fantom5". Currently accepts "fantom5","omnipath", or "custom".
# ?LoadCustom
# LoadCustom()

# head(db_mouse_df)
# dim(db_mouse_df)

head(LR_database_df)
dim(LR_database_df)

NICHES_output <- RunNICHES(object = test,
                                LR.database = "custom",
                                custom_LR_database = LR_database_df_fix,
                                # species = "human",
                                assay = "alra",
                                position.x = 'x',
                                position.y = 'y',
                                k = 8,
                                cell_types = "BayesSpace",
                                min.cells.per.ident = 0,
                                min.cells.per.gene = NULL,
                                meta.data.to.map = c('dataset','BayesSpace'),
                                CellToCell = F,CellToSystem = F,SystemToCell = F,
                                CellToCellSpatial = F,CellToNeighborhood = F,NeighborhoodToCell = T)

# NICHES outputs a list of objects. Each object contains a certain style of cell-system signaling atlas. Above, we have only calculated a single one of interest, namely, individual cellular microenvironment. We next isolate this output and embed using UMAP to visualize the microenvironemnt of each cell.
niche <- NICHES_output[['NeighborhoodToCell']]
Idents(niche) <- niche@meta.data$ReceivingType

# Scale and visualize
niche <- ScaleData(niche)
niche <- FindVariableFeatures(niche,selection.method = "disp")
niche <- RunPCA(niche)
ElbowPlot(niche,ndims = 50)
niche <- RunUMAP(niche,dims = 1:10)

# as.character(niche$BayesSpace) == as.character(niche$ReceivingType)
p4 <- DimPlot(niche,reduction = 'umap',pt.size = 0.5,shuffle = T, label = T) +ggtitle('Cellular Microenvironment')

p5 <- DimPlot(niche,reduction = 'umap',pt.size = 0.5,shuffle = T,split.by = 'ReceivingType',ncol = 3)

p4+p5

# p2 <- DimPlot(niche,reduction = 'umap',pt.size = 0.5,shuffle = T,split.by = 'BayesSpace') +ggtitle('Cellular Microenvironment')
# We can already see, from this plot, some notable overlap between the microenvironments of celltypes 1 & 7 and celltypes 6 & 3. Let's explore this more deeply by finding signaling mechanisms specific to each celltype niche, plotting some of the results in heatmap form:

# Find markers
mark <- FindAllMarkers(niche,min.pct = 0.25,only.pos = T,test.use = "roc")
GOI_niche <- mark %>% group_by(cluster) %>% top_n(10,myAUC)

mark %>%
  filter(str_detect(gene,"SPP1"))

mark %>%
  filter(str_detect(gene,"HLA"))

DoHeatmap(niche,features = unique(GOI_niche$gene))+ 
  scale_fill_gradientn(colors = c("grey","white", "blue"))

# This confirms that celltypes 1 & 7 and 6 & 3 do indeed have some shared character.
# We can further confirm that identified celltype specific signaling mechanisms are indeed specific to tissue regions in which those cells are found, by plotting matched ligand and receptor pairs:

# Check that these make sense and print little plots
DefaultAssay(test) <- 'alra'
p1 <- SpatialFeaturePlot(test, crop = TRUE, features = "SPP1",slot = "data",min.cutoff =  'q1',
                         max.cutoff = 'q99',pt.size.factor = spot_size)+ggtitle("Ligand")+theme(legend.position = "right")
p2 <- SpatialFeaturePlot(test, crop = TRUE, features = "CD44",slot = "data",min.cutoff =  'q1',
                         max.cutoff = 'q99',pt.size.factor = spot_size)+ggtitle("Receptor")+theme(legend.position = "right")

ggpubr::ggarrange(p1,p2)

# Further, and perhaps more usefully, we can map over the output from NICHES onto the original spatial object as follows:

# Add Niches output as an assay
niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
test[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
DefaultAssay(test) <- "NeighborhoodToCell"
test <- ScaleData(test)

# Which allows direct visualization of niche interactions of interest in a spatial context:

# the spacer between ligand and receptro is the following 8212
# intToUtf8(8212)
# 3. Linux
# You can enter an em dash using Unicode input:
#   Press Ctrl + Shift + U.
# Type 2014 (hexadecimal code for 8212).
# Press Enter or Space.

# Plot celltype specific niche signaling
SpatialFeaturePlot(test,
                   features = c("SPP1—CD44"),
                   slot = 'scale.data',pt.size.factor = spot_size)

# plot the panel for all the significant LR pairs
LR_sig_pairs <- LR_sig %>%
  select(ligand,receptor) %>%
  separate_rows(ligand,sep = "_") %>%
  separate_rows(receptor,sep = "_") %>%
  # make the pairs unique
  mutate(pair = paste(ligand,receptor,sep = "_")) %>%
  group_by(pair) %>%
  summarise() %>%
  separate(pair,into = c("ligand.symbol","receptor.symbol"),sep = "_") %>%
  mutate(pair = paste0(ligand.symbol,intToUtf8(8212),receptor.symbol)) %>%
  pull(pair)

SpatialFeaturePlot(test,
                   features = LR_sig_pairs,
                   slot = 'scale.data',pt.size.factor = spot_size*2)
ggsave("../out/plot/110_NICHES_V01_panelSigLR.pdf",width = 20,height = 20)

# pull the scaled score to plot the values per manual annotation
df_exp <- test@assays$SCT@scale.data %>%
  .[c("SPP1","CD44"),] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene,names_to = "barcode",values_to = "scaled_norm_exp")

df_exp2 <- test@assays$NeighborhoodToCell@scale.data %>%
  .[c("SPP1—CD44","SPP1—ITGA5"),] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene,names_to = "barcode",values_to = "scaled_norm_exp")

# general meta
df_meta <- test@meta.data %>%
  rownames_to_column("barcode") %>%
  # add more info based on the manual segmentation
  # Martina asked to gather the periplaques categories in one.
  mutate(manual_segmentation2 = case_when(manual_segmentation %in% c("periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6") ~ "WM",
                                          TRUE ~ manual_segmentation)) %>%

  left_join(LUT_sample,by="dataset") %>%
  mutate(manual_segmentation3 = case_when(manual_segmentation2 == "edge" ~ paste0("edge_",type),
                                          TRUE ~ manual_segmentation2))

# add the meta info based on the barcodes
left_join(df_exp,df_meta,by = "barcode") %>%
  ggplot(aes(x = manual_segmentation3,y = scaled_norm_exp))+
  geom_violin()+
  facet_wrap(~gene,scales = "free")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

left_join(df_exp2,df_meta,by = "barcode") %>%
  ggplot(aes(x = manual_segmentation3,y = scaled_norm_exp))+
  geom_boxplot()+
  facet_wrap(~gene,scales = "free")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  


# add the spatial niche clustering to the main object
# dim(niche@meta.data)
# dim(test)
# 
# test@meta.data <- test@meta.data %>% rownames_to_column("barcodes") %>%
#   left_join(niche@meta.data %>%
#               select(barcodes = ReceivingCell,niches_cluster = seurat_clusters,ReceivingType),by = "barcodes") %>%
#   column_to_rownames("barcodes")
# 
# identical(test@meta.data$seurat_clusters,test@meta.data$ReceivingType)
# identical(test@meta.data$seurat_clusters,test@meta.data$niches_cluster)
# 
# p3 <- SpatialDimPlot(test, label = TRUE,group.by = 'niches_cluster', label.size = 3,pt.size.factor = spot_size)
# p4 <- SpatialDimPlot(test, label = TRUE,group.by = 'seurat_clusters', label.size = 3,pt.size.factor = spot_size)
# 
# p3 + p4
# 
# Idents(test) <- 'niches_cluster'
# p32 <- SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE,pt.size.factor = spot_size)+plot_annotation('Niche Clusters')
# 
# Idents(test) <- 'seurat_clusters'
# p42 <- SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE,pt.size.factor = spot_size)+plot_annotation('seurat clusters')
# 
# wrap_plots(list(p32,p42))
