# AIM ---------------------------------------------------------------------
# attempt the integration following SIB processing

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(UpSetR)
library(pals)
library(scales)
library(clustree)
library(patchwork)
library(ComplexHeatmap)
library(matrixStats)
library(circlize)

# read in the data --------------------------------------------------------
list_brain <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")

# test common HVG ---------------------------------------------------------
# merge the tables
# pull the names of the individula objects
data.list.id <- names(list_brain)
data.combined.all <- merge(list_brain[[1]], y = list_brain[-1], add.cell.ids = data.list.id, project = "visium_brain")

# In order to perform dimensionality reduction, we first need to select variable features of all the slices. To get a good representation, we take the intersect (i.e. genes that are variable in all slices):
# intersect(VariableFeatures(list_brain$V01),
#           VariableFeatures(list_brain$V02))

list_var_features <- purrr::map(list_brain, function(x){
  VariableFeatures(x)
})

# check the instersection
upset(fromList(list_var_features), order.by = "freq",nsets = 1000)

# try to integrate with the common ones
# pull the common ones
top_common_hvf <- purrr::reduce(list_var_features, intersect)

test_common <- data.combined.all
VariableFeatures(test_common) <- top_common_hvf

# How many variable features do we have? Why did we select fewer genes than the default (check ?VariableFeatures)?

# Now that we have selected the most variable features, we can generate a PCA based on the normalized and scaled data of those:
test_common <- RunPCA(test_common, assay = "SCT", npcs = 50, verbose = FALSE)

# define many colors
color_id <- alphabet(length(unique(test_common$dataset)))
show_col(color_id)
names(color_id) <- unique(test_common$dataset)

# plot the PCA
DimPlot(test_common, reduction = "pca", group.by = "dataset")+scale_color_manual(values = color_id)

# Based on the PCA, we can create a UMAP to get a representation of all 50 dimensions in a two dimensional space:
test_common <- RunUMAP(test_common, reduction = "pca", dims = 1:10)

# plot the UMAP before the integration
DimPlot(test_common, reduction = "umap", group.by = "dataset")+scale_color_manual(values = color_id)

# full integration --------------------------------------------------------
# To integrate the slices, we first need to select integration features. These are genes that are variable in all slices. We then prepare the data for integration, find the integration anchors (i.e. spots that are within each others neigbourhoods), and integrate the data:

# You can safely ignore the warning: Warning: Different cells and/or features from existing assay SCT. See this issue.
features <- SelectIntegrationFeatures(list_brain)
list_brain <- PrepSCTIntegration(list_brain, anchor.features = features)

anchors <- FindIntegrationAnchors(
  list_brain,
  normalization.method = "SCT",
  anchor.features = features
)

seu <- IntegrateData(anchors, normalization.method = "SCT")

# Because we re-do the dimensionality reduction, we also again extract the variable features, run the PCA and the UMAP:
seu <- FindVariableFeatures(seu)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:50)
(DimPlot(seu, reduction = "umap",group.by = "dataset") + scale_color_manual(values = color_id))

# Identifying clusters
# Seurat implements a graph-based clustering approach. Distances between the spots are calculated based on previously identified PCs. Briefly, Seurat identifies clusters of spots by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First, it identifies k-nearest neighbors (KNN) and constructs the SNN graph. Then it optimizes the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B.

# The FindClusters function implements the procedure, and contains a resolution parameter that sets the granularity of the downstream clustering, with increased values leading to a greater number of clusters.

resolution_vector <- seq(0.1,1,0.1)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:50)
seu <- FindClusters(object = seu,
                    resolution = resolution_vector,
                    verbose=T)

# Some new columns appeared in the metadata data frame after the clustering, each representing the cluster ID per spot for a given resolution:
colnames(seu@meta.data)

# To get an overview of the clustering over the different resolutions, we can use clustree to get an idea:
clustree(seu, prefix = "integrated_snn_res.")
ggsave("out/image/UMAPCluster_tree_SIB.pdf",width = 15,height = 10)

# this is not (necessarily) the correct answer to the previous question!
res <- "integrated_snn_res.0.1"

# Now that we have selected a resolution, we can color both the UMAP and the slices accordingly. First we defnie some appropriate colors, then we plot the UMAP with DimPlot and the slices with SpatialPlot.
# 
# define a color palette based on the number of clusters
nclust <- seu[[res]] |> unique() |> nrow()
cluster_cols <- DiscretePalette(nclust, palette = "polychrome")

(DimPlot(seu, reduction = "umap",group.by = "dataset") + scale_color_manual(values = color_id)) +
DimPlot(seu,
        group.by = res,
        shuffle = TRUE,
        cols = cluster_cols)

# plot clusters per slide
DimPlot(seu,
        group.by = res,split.by = "dataset",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 3)

# show the clusters split by slide
# SpatialPlot(seu, pt.size.factor = 2) + 
#   plot_layout(guides='collect') &
#   theme(legend.position = "none") &
#   scale_fill_manual(values = cluster_cols)

# show the cluster by top anntation of the deconvolution approach
DimPlot(seu,
        group.by = res,split.by = "class_top",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 4)

DimPlot(seu,
        group.by = res,split.by = "manual_segmentation",
        shuffle = TRUE,
        cols = cluster_cols,ncol = 4)

# -------------------------------------------------------------------------
# calculate the relative proportion of cell in each cluster compared to different annotations.
prop_table_tot <- seu@meta.data %>%
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(integrated_snn_res.0.1,class_top) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(integrated_snn_res.0.1,class_top,prop) %>%
  pivot_wider(names_from = integrated_snn_res.0.1,values_from = prop,values_fill=0) %>%
  column_to_rownames("class_top") %>%
  as.matrix()

# try to plot the z score
zscore_table_tot <- prop_table_tot %>%
  as.data.frame() %>%
  rownames_to_column("class_top") %>%
  pivot_longer(names_to = "cluster",values_to = "prop",-class_top) %>%
  group_by(class_top) %>%
  mutate(z_score = (prop - mean(prop))/sd(prop)) %>%
  select(cluster,class_top,z_score) %>%
  pivot_wider(names_from = cluster,values_from = z_score) %>%
  column_to_rownames("class_top") %>%
  as.matrix()

rowSums(zscore_table_tot)
rowSds(zscore_table_tot)

# clip the z scores
# Define a color function with a specific range
color_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# pdf("../../out/image/129_heatmap_label_transfer_refWMImmune_queryIMMLYMsubset_raw_res0.9.pdf",height = 4,width = 4)
Heatmap(prop_table_tot,
        name = "prop", 
        column_title = "deconvolution top annotation",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(11))
# dev.off()
Heatmap(zscore_table_tot,
        name = "z score", 
        column_title = "deconvolution top annotation",
        cluster_rows = T,cluster_columns = T,col = color_fun)


# -------------------------------------------------------------------------
prop_table_tot2 <- seu@meta.data %>%
  # remove the spots where there is no annotation for the class manual segmentation
  filter(!is.na(manual_segmentation)) %>%
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(integrated_snn_res.0.1,manual_segmentation) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(integrated_snn_res.0.1,manual_segmentation,prop) %>%
  pivot_wider(names_from = integrated_snn_res.0.1,values_from = prop,values_fill=0) %>%
  column_to_rownames("manual_segmentation") %>%
  as.matrix()

# try to plot the z score
zscore_table_tot2 <- prop_table_tot2 %>%
  as.data.frame() %>%
  rownames_to_column("manual_segmentation") %>%
  pivot_longer(names_to = "cluster",values_to = "prop",-manual_segmentation) %>%
  group_by(manual_segmentation) %>%
  mutate(z_score = (prop - mean(prop))/sd(prop)) %>%
  select(cluster,manual_segmentation,z_score) %>%
  pivot_wider(names_from = cluster,values_from = z_score) %>%
  column_to_rownames("manual_segmentation") %>%
  as.matrix()

rowSums(zscore_table_tot2)
rowSds(zscore_table_tot2)

color_fun2 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# pdf("../../out/image/129_heatmap_label_transfer_refWMImmune_queryIMMLYMsubset_raw_res0.9.pdf",height = 4,width = 4)
Heatmap(prop_table_tot2,
        name = "prop", 
        column_title = "manual annotation",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(11))
# dev.off()
Heatmap(zscore_table_tot2,
        name = "z score", 
        column_title = "manual annotation",
        cluster_rows = T,cluster_columns = T,col = color_fun2)


# -------------------------------------------------------------------------
# try to check gene expression for different known markers by cluster
DefaultAssay(seu) <- "SCT"
shortlist_features_list_long <- list(
  IMMUNE = c("LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC","AIF1","HLA-DRA","TYROBP"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS =  c("SKAP1", "CD8A", "CD2"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"),
  SCHWANN = c("PMP22","MPZ","PRX"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
  STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
)

# make a shortlist of the markers
shortlist_features_list_short <- list(
  IMMUNE = c("CX3CR1","P2RY12","C3","CSF1R", "CD74","C1QB"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS =  c("SKAP1", "CD8A", "CD2"),
  OLIGO = c("MOG","MBP","MAG"),
  OPC = c("NLGN4X","OLIG1","OLIG2"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1"),
  NEURONS = c("CUX2","SYP", "NEFL","SYT1"),
  VAS = c("VWF","FLT1","CLDN5","PDGFRB"),
  SCHWANN = c("PMP22","MPZ","PRX"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
  STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long01 <- DotPlot(seu,
                       features = shortlist_features_list_long,
                       dot.scale = 8,
                       cluster.idents = T,
                       group.by = "integrated_snn_res.0.1",assay = "SCT") +
  RotatedAxis() +
  labs(title = "integrated_snn_res.0.1")+
  theme(strip.text = element_text(angle = 90))
# ggsave(plot=test_long01,"../../out/image/129_DotplotLong_res0.1_MG_subcluster.pdf",width = 30,height = 6)

# score the signatures, both long and short
seu <- Seurat::AddModuleScore(seu,
                              features = shortlist_features_list_long,
                              name = "_score_long")

df_rename_long <- data.frame(names = seu@meta.data %>%
                               colnames() %>%
                               str_subset("_score_long"),
                             rename = paste0("scoreLong_",names(shortlist_features_list_long)))

lookup_long <- df_rename_long$names
names(lookup_long) <- df_rename_long$rename

# rename the columns
seu@meta.data <- dplyr::rename(seu@meta.data,all_of(lookup_long))

# plot the scores from AddModuleScore
list_plot_02_long <- lapply(df_rename_long$rename,function(x){
  plot <- FeaturePlot(seu,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_long)
# ggsave("../../out/image/129_UMAPCluster_long_MG_subcluster.pdf",width = 22,height = 12)

# # same as above but as violin plot
# list_plot <- lapply(df_rename_long$rename, function(x){ 
#   test <- VlnPlot(object = seu,features = x, group.by = "integrated_snn_res.0.1",raster = T)
#   return(test)
# })
# 
# # make it a dataframe
# # x <- list_plot[[1]]
# df_violin <- lapply(list_plot,function(x){ 
#   df <- x[[1]]$data 
#   
#   # extract the name of the gene 
#   feature <- colnames(df)[1] 
#   
#   df %>% 
#     mutate(feature = feature) %>% 
#     setNames(c("value","ident","feature")) 
# }) %>% 
#   bind_rows()
# 
# head(df_violin) 
# 
# # set.seed(123)
# df_plot_violin <- df_violin %>%
#   group_by(ident,feature) %>%
#   # sample_n(size = 400,replace = F) %>%
#   ungroup()
# 
# df_plot_violin_summary <- df_plot_violin %>%
#   group_by(feature) %>%
#   summarise(med_score = median(value))
# 
# df_plot_violin %>%
#   ggplot(aes(y = ident, x = value)) +
#   geom_violin(scale = "width")+
#   #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) +
#   geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) +
#   facet_wrap(~feature,nrow = 1,scales = "free") +
#   theme_bw() +
#   geom_vline(data = df_plot_violin_summary,aes(xintercept = med_score),linetype="dashed",col="red") +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(hjust = 1,angle = 45))

# -------------------------------------------------------------------------
# save the object so fat
saveRDS(seu,"out/object/100_seu_SIB.rds")

