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
# list_brain_NICHES <- readRDS("../out/object/list_brain_NICHES_all_135_net.up_cellID_SEN_AST.rds")
list_brain_NICHES <- readRDS("../out/object/list_brain_NICHES_all.rds")
LR_sig <- read_tsv("../data/138_net.up_cellID_SEN_VAS.tsv")

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

LUT_sample <- read_csv("../../Visium_brainMS_single_sample/data/sample_classification.csv") %>%
  mutate(dataset = paste0("brain_",sample))

# plot --------------------------------------------------------------------
# make the full plot
# sl <- list_brain_NICHES$V01
# nm <- "V01"

pmap(list(list_brain_NICHES,names(list_brain_NICHES)),function(sl,nm){
  
  # print name
  print(nm)
  
  # pull one slide to test
  # sl <- list_brain_NICHES$V01
  
  # start with one sample
  niche_test <- sl$niche
  # Idents(niche_test)
  test_test <- sl$test
  # Idents(niche) <- niche_test@meta.data$ReceivingType
  
  # spot size
  spot_size <- test_test@images$slice1@scale.factors$fiducial/test_test@images$slice1@scale.factors$hires
  
  # Further, and perhaps more usefully, we can map over the output from NICHES onto the original spatial object as follows:
  
  # # Add Niches output as an assay
  # niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], slot = 'data')
  # colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
  # test[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data )
  # DefaultAssay(test) <- "NeighborhoodToCell"
  # test <- ScaleData(test)
  
  # Which allows direct visualization of niche interactions of interest in a spatial context:
  
  # the spacer between ligand and receptro is the following 8212
  # intToUtf8(8212)
  # 3. Linux
  # You can enter an em dash using Unicode input:
  #   Press Ctrl + Shift + U.
  # Type 2014 (hexadecimal code for 8212).
  # Press Enter or Space.
  
  # make sure the specific slot is selected for plotting
  DefaultAssay(test_test) <- "NeighborhoodToCell"
  
  # # Plot celltype specific niche signaling
  # SpatialFeaturePlot(test_test,
  #                    features = c("SPP1—CD44"),
  #                    slot = 'scale.data',pt.size.factor = spot_size)
  
  
  plot <- SpatialFeaturePlot(test_test,
                             features = LR_sig_pairs,
                             slot = 'scale.data',pt.size.factor = spot_size*2)
  ggsave(plot = plot,paste0("../out/plot/110_NICHES_",nm,"_panelSigLR_VASSenmayo.pdf"),width = 20,height = 20)
})

# plot individual ---------------------------------------------------------
test_test <- list_brain_NICHES$V01$test
test_niche <- list_brain_NICHES$V01$niche

# plot dimensionality reduction
p4 <- DimPlot(test_niche,reduction = 'umap',pt.size = 0.5,shuffle = T, label = T) +ggtitle('Cellular Microenvironment')
p5 <- DimPlot(test_niche,reduction = 'umap',pt.size = 0.5,shuffle = T,split.by = 'ReceivingType',ncol = 3)

p4+p5

# print a specific ligand and receptor
# for Seurat v5.1.0 there is the need to update the size of the spots
spot_size <- test_test@images$slice1@scale.factors$fiducial/test_test@images$slice1@scale.factors$hires
DefaultAssay(test_test) <- 'alra'
p1 <- SpatialFeaturePlot(test_test, crop = TRUE, features = "SPP1",slot = "data",min.cutoff =  'q1',
                         max.cutoff = 'q99',pt.size.factor = spot_size)+ggtitle("Ligand")+theme(legend.position = "right")
p2 <- SpatialFeaturePlot(test_test, crop = TRUE, features = "CD44",slot = "data",min.cutoff =  'q1',
                         max.cutoff = 'q99',pt.size.factor = spot_size)+ggtitle("Receptor")+theme(legend.position = "right")

ggpubr::ggarrange(p1,p2)

# Find markers
df_mark <- FindAllMarkers(test_niche,min.pct = 0.25,only.pos = T,test.use = "roc")
df_GOI_niche <- df_mark %>% group_by(cluster) %>% top_n(10,myAUC)

df_mark %>%
  filter(str_detect(gene,"SPP1"))

df_mark %>%
  filter(str_detect(gene,"HLA"))

DoHeatmap(test_niche,features = unique(df_GOI_niche$gene))+ 
  scale_fill_gradientn(colors = c("grey","white", "blue"))

# -------------------------------------------------------------------------
# loop the extraction of the scaled expression for the ligand and receptor
# x <- list_brain_NICHES$V01

list_meta <- lapply(list_brain_NICHES, function(x){
  
  # extract the objects
  test_test <- x$test
  test_niche <- x$niche
  
  # define the GOI presents
  GOI <- LR_sig_pairs %>% str_split(pattern = intToUtf8(8212)) %>% unlist() %>% unique() %>%
    data.frame() %>%
    setNames("gene") %>%
    filter(gene %in% rownames(test_test@assays$SCT@scale.data)) %>%
    pull(gene)
  
  # pull the scaled score to plot the values per manual annotation
  df_exp <- test_test@assays$SCT@scale.data %>%
    .[GOI,] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene,names_to = "barcode",values_to = "scaled_norm_exp")
  
  # pull the scaled expression of the LR pairs
  df_exp2 <- test_test@assays$NeighborhoodToCell@scale.data %>%
    .[LR_sig_pairs,] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene,names_to = "barcode",values_to = "scaled_norm_exp")
  
  # general meta
  df_meta <- test_test@meta.data %>%
    rownames_to_column("barcode") %>%
    # add more info based on the manual segmentation
    # Martina asked to gather the periplaques categories in one.
    mutate(manual_segmentation2 = case_when(manual_segmentation %in% c("periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6") ~ "WM",
                                            TRUE ~ manual_segmentation)) %>%
    
    left_join(LUT_sample,by="dataset") %>%
    mutate(manual_segmentation3 = case_when(manual_segmentation2 == "edge" ~ paste0("edge_",type),
                                            TRUE ~ manual_segmentation2))
  
  # average the expression per gene and manual area
  df_avg_gene <- left_join(df_exp,df_meta,by = "barcode") %>%
    group_by(gene,dataset,type,manual_segmentation3) %>%
    summarise(avg_scaled_exp = mean(scaled_norm_exp))
  
  df_avg_pair <- left_join(df_exp2,df_meta,by = "barcode") %>%
    group_by(gene,dataset,type,manual_segmentation3) %>%
    summarise(avg_scaled_exp = mean(scaled_norm_exp))
  
  list_out <- list(df_avg_gene = df_avg_gene,df_avg_pair = df_avg_pair)
  return(list_out)
})

# make a unique meta
meta_full_gene <- lapply(list_meta,function(x){
  x$df_avg_gene
}) %>%
  bind_rows()

meta_full_pair <- lapply(list_meta,function(x){
  x$df_avg_pair
}) %>%
  bind_rows()

# plot the pairs
meta_full_pair %>%
  filter(!is.na(manual_segmentation3)) %>%
  ggplot(aes(x = manual_segmentation3,y = avg_scaled_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2,height = 0))+
  facet_wrap(~gene,scales = "free_y")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),strip.background = element_blank())

color_id <- tableau20(length(unique(meta_full_pair$dataset)))
# color_id <- rev(color_id)
# color_id <- c("#AEC7E8","#1F77B4","#FF7F0E","#FFBB78","#2CA02C","#98DF8A","#D62728","#FF9896" ,"#9467BD" ,"#C5B0D5" ,"#8C564B" ,"#C49C94" ,"#E377C2")
show_col(color_id)
names(color_id) <- unique(meta_full_pair$dataset %>% unlist())

meta_full_pair %>%
  filter(!is.na(manual_segmentation3)) %>%
  ggplot(aes(x = manual_segmentation3,y = avg_scaled_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2,height = 0),aes(col=dataset))+
  scale_color_manual(values = color_id)+
  facet_wrap(~gene,scales = "free_y")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),strip.background = element_blank())

# save also the full table
meta_full_pair %>%
  write_csv("../out/table/111_meta_full_pair_NICHES_scaled_VASSenmayo.csv")

# -------------------------------------------------------------------------
# loop the extraction of the normalized expression for the ligand and receptor
list_meta_norm <- lapply(list_brain_NICHES, function(x){
  
  # extract the objects
  test_test <- x$test
  test_niche <- x$niche
  
  # define the GOI presents
  GOI <- LR_sig_pairs %>% str_split(pattern = intToUtf8(8212)) %>% unlist() %>% unique() %>%
    data.frame() %>%
    setNames("gene") %>%
    filter(gene %in% rownames(test_test@assays$SCT@scale.data)) %>%
    pull(gene)
  
  # pull the scaled score to plot the values per manual annotation
  df_exp <- test_test@assays$SCT@data %>%
    .[GOI,] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene,names_to = "barcode",values_to = "scaled_norm_exp")
  
  # pull the scaled expression of the LR pairs
  df_exp2 <- test_test@assays$NeighborhoodToCell@data %>%
    .[LR_sig_pairs,] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene,names_to = "barcode",values_to = "scaled_norm_exp")
  
  # general meta
  df_meta <- test_test@meta.data %>%
    rownames_to_column("barcode") %>%
    # add more info based on the manual segmentation
    # Martina asked to gather the periplaques categories in one.
    mutate(manual_segmentation2 = case_when(manual_segmentation %in% c("periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6") ~ "WM",
                                            TRUE ~ manual_segmentation)) %>%
    
    left_join(LUT_sample,by="dataset") %>%
    mutate(manual_segmentation3 = case_when(manual_segmentation2 == "edge" ~ paste0("edge_",type),
                                            TRUE ~ manual_segmentation2))
  
  # average the expression per gene and manual area
  df_avg_gene <- left_join(df_exp,df_meta,by = "barcode") %>%
    group_by(gene,dataset,type,manual_segmentation3) %>%
    summarise(avg_norm_exp = mean(scaled_norm_exp))
  
  df_avg_pair <- left_join(df_exp2,df_meta,by = "barcode") %>%
    group_by(gene,dataset,type,manual_segmentation3) %>%
    summarise(avg_norm_exp = mean(scaled_norm_exp))
  
  list_out <- list(df_avg_gene = df_avg_gene,df_avg_pair = df_avg_pair)
  return(list_out)
})

# make a unique meta
meta_full_gene_norm <- lapply(list_meta_norm,function(x){
  x$df_avg_gene
}) %>%
  bind_rows()

meta_full_pair_norm <- lapply(list_meta_norm,function(x){
  x$df_avg_pair
}) %>%
  bind_rows()

# plot the pairs
meta_full_pair_norm %>%
  filter(!is.na(manual_segmentation3)) %>%
  ggplot(aes(x = manual_segmentation3,y = avg_norm_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2,height = 0))+
  facet_wrap(~gene,scales = "free_y")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),strip.background = element_blank())

# save also the full table
meta_full_pair_norm %>%
  write_csv("../out/table/111_meta_full_pair_NICHES_norm_VASSenmayo.csv")
