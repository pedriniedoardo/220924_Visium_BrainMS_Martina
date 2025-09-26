# AIM ---------------------------------------------------------------------
# finalize some custom plotting

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)
library(scales)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

show_col(SpatialColors(10))

# read in the data --------------------------------------------------------
# read in the processed data of the brain slices
list_brain <- readRDS(file = "out/object/list_brain_all_BayesSpace1000_q10.rds")
list_brain2 <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")


# plot panel --------------------------------------------------------------
limit <- max(list_brain$V01$nCount_Spatial)

SpatialFeaturePlot(list_brain$V01, features = "nCount_Spatial",crop = T) +
  scale_fill_gradientn(colours = SpatialColors(n = 100),
                       breaks = c(0, round(limit/2),limit),
                       labels = c("Low", "Medium", "High"),
                       limits = c(0, limit))+
  theme(legend.position = "top")

list_plot <- pmap(list(list_brain,names(list_brain)), function(sl,nm){
  
  limit <- max(sl$nCount_Spatial)
  
  p <- SpatialFeaturePlot(sl, features = "nCount_Spatial",crop = T) +
    scale_fill_gradientn(colours = SpatialColors(n = 100),
                         breaks = c(0, round(limit/2),limit),
                         labels = c("Low", "Medium", "High"),
                         limits = c(0, limit)) +
    theme(legend.position = "top") +
    ggtitle(nm)
  
  return(p)
})

wrap_plots(list_plot)
ggsave("out/image/99_plot_ncount_all_slides.pdf",width = 16,height = 16)


list_brain2$V01@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode,AST,IMM,LYM,NEU,OLIGO,OPC,VAS) %>%
  pivot_longer(cols = -barcode, names_to = "celltype", values_to = "value") %>%
  group_by(barcode) %>%
  summarise(tot = sum(value)) %>%
  summary()

# plot bayespace ----------------------------------------------------------
# martina asked to plot bayespace fllowing a custom coloring
# V01
test_v01 <- list_brain$V01
SpatialDimPlot(test_v01,group.by = "BayesSpace")
Idents(test_v01) <- "BayesSpace"
SpatialDimPlot(test_v01)
SpatialDimPlot(test_v01, cells.highlight = CellsByIdentities(object = test_v01), facet.highlight = TRUE)
# group the clusters
meta_v01 <- test_v01@meta.data
meta_v01_update <- meta_v01 %>%
  mutate(BayesSpace_group = case_when(
    BayesSpace %in% c("1","5","6","10") ~ "core",
    BayesSpace %in% c("9") ~ "cortex",
    BayesSpace %in% c("4") ~ "edge",
    BayesSpace %in% c("3","7","8","2") ~ "WM"
  ))

# update the plot
test_v01@meta.data <- meta_v01_update

#
Idents(test_v01) <- "BayesSpace_group"
SpatialDimPlot(test_v01)
ggsave("out/image/99_plot_bayespace_group_v01.pdf",width = 9,height = 9)
SpatialDimPlot(test_v01, cells.highlight = CellsByIdentities(object = test_v01), facet.highlight = TRUE)
ggsave("out/image/99_plot_bayespace_group_v01_split.pdf",width = 9,height = 9)


# V02
test_v02 <- list_brain$V02
SpatialDimPlot(test_v02,group.by = "BayesSpace")
Idents(test_v02) <- "BayesSpace"
SpatialDimPlot(test_v02)
SpatialDimPlot(test_v02, cells.highlight = CellsByIdentities(object = test_v02), facet.highlight = TRUE)
# group the clusters
meta_v02 <- test_v02@meta.data
meta_v02_update <- meta_v02 %>%
  mutate(BayesSpace_group = case_when(
    BayesSpace %in% c("2","5","6") ~ "core",
    BayesSpace %in% c("1","7") ~ "cortex",
    BayesSpace %in% c("3") ~ "edge",
    BayesSpace %in% c("9","10","8","4") ~ "WM"
  ))

# update the plot
test_v02@meta.data <- meta_v02_update

#
Idents(test_v02) <- "BayesSpace_group"
SpatialDimPlot(test_v02)
ggsave("out/image/99_plot_bayespace_group_v02.pdf",width = 9,height = 9)
SpatialDimPlot(test_v02, cells.highlight = CellsByIdentities(object = test_v02), facet.highlight = TRUE)
ggsave("out/image/99_plot_bayespace_group_v02_split.pdf",width = 9,height = 9)

# V07
test_v07 <- list_brain$V07
SpatialDimPlot(test_v07,group.by = "BayesSpace")
Idents(test_v07) <- "BayesSpace"
SpatialDimPlot(test_v07)
SpatialDimPlot(test_v07, cells.highlight = CellsByIdentities(object = test_v07), facet.highlight = TRUE)
# group the clusters
meta_v07 <- test_v07@meta.data
meta_v07_update <- meta_v07 %>%
  mutate(BayesSpace_group = case_when(
    BayesSpace %in% c("2","5","3","9") ~ "core",
    BayesSpace %in% c("1","7") ~ "cortex",
    BayesSpace %in% c("4") ~ "edge",
    BayesSpace %in% c("6","10","8") ~ "WM"
  ))

# update the plot
test_v07@meta.data <- meta_v07_update

#
Idents(test_v07) <- "BayesSpace_group"
SpatialDimPlot(test_v07)
ggsave("out/image/99_plot_bayespace_group_v07.pdf",width = 9,height = 9)
SpatialDimPlot(test_v07, cells.highlight = CellsByIdentities(object = test_v07), facet.highlight = TRUE)
ggsave("out/image/99_plot_bayespace_group_v07_split.pdf",width = 9,height = 9)


# plot and table the proportions per deconvolution ------------------------
# looking back at the script used for the SPOTlight deconvolution, I have assigned a value of 0 for the celltypes that had a proportion lower than 0.08. This was done following the vignette. This means that the summarization per spot might not add up to 1.
# Neverhless, this are all proportions.
list_brain2$V01@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode,AST,IMM,LYM,NEU,OLIGO,OPC,VAS) %>%
  pivot_longer(cols = -barcode, names_to = "celltype", values_to = "value") %>%
  group_by(barcode) %>%
  summarise(tot = sum(value)) %>%
  summary()

list_brain2$V01@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode,AST,IMM,LYM,NEU,OLIGO,OPC,VAS) %>%
  pivot_longer(cols = -barcode, names_to = "celltype", values_to = "value") %>%
  filter(value > 0) %>%
  summary()

# add the annotation suggested by Martina
# LUT_sample <- read_csv("data/sample_classification.csv") %>%
#   mutate(dataset = paste0("brain_",sample))
# 
# # Martina asked to gather the periplaques categories in one.
# # sl <- list_brain2$V01
# list_brain3 <- lapply(list_brain2, function(sl){
#   sl@meta.data$manual_segmentation2 <- 
#     sl@meta.data %>%
#     mutate(manual_segmentation2 = case_when(manual_segmentation %in% c("periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6") ~ "WM",
#                                             TRUE ~ manual_segmentation)) %>%
#     pull(manual_segmentation2)
#   
#   # martina suggested to split the edge also by its pathological class
#   sl@meta.data$manual_segmentation3 <-
#     sl@meta.data %>%
#     left_join(LUT_sample,by="dataset") %>%
#     mutate(manual_segmentation3 = case_when(manual_segmentation2 == "edge" ~ paste0("edge_",type),
#                                             TRUE ~ manual_segmentation2)) %>%
#     pull(manual_segmentation3)
#   
#   return(sl)
# })
# 
# saveRDS(list_brain3,"out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation3.rds")

list_brain3 <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation3.rds")

# pull the meta from all the slides
meta_all <- lapply(list_brain3, function(sl){
  sl@meta.data %>%
    rownames_to_column("barcode")
}) %>%
  bind_rows()

# calculate the averege prop of deconvolution per spot per menual segmentation3
# regardless of the slide
summary_01 <- meta_all %>%
  filter(!is.na(manual_segmentation3)) %>%
  select(barcode,manual_segmentation3,AST,IMM,LYM,NEU,OLIGO,OPC,VAS) %>%
  pivot_longer(cols = -c(barcode,manual_segmentation3), names_to = "celltype", values_to = "value") %>%
  group_by(manual_segmentation3,celltype) %>%
  summarise(mean_prop = mean(value),.groups = "drop")

write_tsv(summary_01,"out/table/99_deconvolution_summary_01.tsv")

hm01 <- summary_01 %>%
  group_by(celltype) %>%
  # mutate(z_score = scale(mean_prop)[,1])
  mutate(z_score = (mean_prop-mean(mean_prop))/sd(mean_prop)) %>%
  ungroup() %>%
  select(-mean_prop) %>%
  pivot_wider(names_from = "manual_segmentation3", values_from = "z_score") %>%
  column_to_rownames("celltype")

Heatmap(hm01,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

Heatmap(hm01,
        col = colorRamp2(
          breaks = seq(-2, 2, length.out = 100),
          colors = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100)
        ))

Heatmap(hm01,colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(100))

Heatmap(hm01,
        col = colorRamp2(
          breaks = seq(-2, 2, length.out = 100),
          colors = viridis::viridis(100,option = "turbo")
        ))


# add the slide info
summary_02 <- meta_all %>%
  filter(!is.na(manual_segmentation3)) %>%
  select(dataset,barcode,manual_segmentation3,AST,IMM,LYM,NEU,OLIGO,OPC,VAS) %>%
  pivot_longer(cols = -c(dataset,barcode,manual_segmentation3), names_to = "celltype", values_to = "value") %>%
  group_by(dataset,manual_segmentation3,celltype) %>%
  summarise(mean_prop = mean(value),.groups = "drop")

write_tsv(summary_02,"out/table/99_deconvolution_summary_02.tsv")


# raw prop
summary_02 %>%
  ggplot(aes(x=manual_segmentation3,y=mean_prop)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~celltype,scale="free") +
  geom_point(aes(y=mean_prop),position = position_jitter(width = 0.1),alpha = 0.5)+
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))

# scaled
summary_02 %>%
  group_by(celltype) %>%
  # mutate(z_score = scale(mean_prop)[,1])
  mutate(z_score = (mean_prop-mean(mean_prop))/sd(mean_prop)) %>%
  ungroup() %>%
  ggplot(aes(x=manual_segmentation3,y=z_score)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~celltype,scale="free") +
  geom_point(aes(y=z_score),position = position_jitter(width = 0.1),alpha = 0.5)+
  theme_bw() +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))


# -------------------------------------------------------------------------
# extract and save the plots with sofia's annotation
list_brain3 <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation3.rds")
# x <- list_brain3$V01
# name <- "V01"
pdf("out/image/99_panel_manual_segmentation_all_martina.pdf",width = 7,height = 6)
pmap(list(list_brain3,names(list_brain3)),function(x,name){
  # Idents(x) <- "manual_segmentation3"
  SpatialDimPlot(x, group.by = "manual_segmentation3",alpha=0.5,crop = F,pt.size.factor = 1) +
    # theme(legend.position = "none") +
    ggtitle(name)
})
dev.off()

# -------------------------------------------------------------------------
# extract and save the plots with sofia's annotation
list_brainBayes <- readRDS("out/object/list_brain_all_BayesSpace1000_q10.rds")
# x <- list_brainBayes$V01
# name <- "V01"
pdf("out/image/99_panel_bayespace_q10.pdf",width = 7,height = 6)
pmap(list(list_brainBayes,names(list_brainBayes)),function(x,name){
  # Idents(x) <- "manual_segmentation3"
  SpatialDimPlot(x, group.by = "BayesSpace",alpha=0.5,crop = F,pt.size.factor = 1) +
    # theme(legend.position = "none") +
    ggtitle(name)
})
dev.off()


