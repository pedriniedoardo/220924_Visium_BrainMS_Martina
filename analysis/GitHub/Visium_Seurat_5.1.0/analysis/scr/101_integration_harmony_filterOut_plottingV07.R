
# read in the data --------------------------------------------------------
data.combined.all_harmony <- readRDS(file = "../out/object/100_data.combined.all_harmony_regressOut_fullAnno.rds")

LUT_sample <- read_csv("../../Visium_brainMS_single_sample/data/sample_classification.csv") %>%
  mutate(dataset = paste0("brain_",sample))

# read the data from the previous analysis
list_brain <- readRDS("../../Visium_brainMS_single_sample/out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")
list_brain_bayespace <- readRDS(file = "../../Visium_brainMS_single_sample/out/object/list_brain_all_BayesSpace1000_q10.rds")


# -------------------------------------------------------------------------
# try plotting the cluster id on the original object
# define the scale factor for the points. this is needed for older objects assay v3 running on Seurat 5.1.0
# ref: https://github.com/satijalab/seurat/issues/9049
spot_size <- list_brain$V07@images$slice1@scale.factors$fiducial/list_brain$V01@images$slice1@scale.factors$hires

SpatialDimPlot(list_brain$V07, group.by = "manual_segmentation",label = F,pt.size = spot_size)
SpatialDimPlot(list_brain_bayespace$V07, group.by = "BayesSpace",label = F,pt.size = spot_size)

meta_V07 <- list_brain$V07@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode)

meta_V07_bayespace <- list_brain_bayespace$V07@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode,BayesSpace)

meta_integrated <- data.combined.all_harmony@meta.data %>%
  rownames_to_column("sample_barcode") %>%
  separate(sample_barcode,into = c("sample_id","barcode"),sep = "_",remove = T) %>%
  filter(sample_id == "V07")

# join all the info
meta_V07_update <- meta_V07 %>%
  left_join(meta_V07_bayespace,by = "barcode") %>%
  left_join(meta_integrated,by = "barcode") %>%
  column_to_rownames("barcode")

test <- list_brain$V07
test@meta.data <- meta_V07_update

# plot the UMAP with all the annotations
p17 <- SpatialDimPlot(test, group.by = "BayesSpace",label = F,pt.size = spot_size)
p27 <- SpatialDimPlot(test, group.by = "manual_segmentation3",label = F,pt.size = spot_size)
p37 <- SpatialDimPlot(test, group.by = "SCT_snn_res.0.1",label = F,pt.size = spot_size)

Idents(test) <- "BayesSpace"
p127 <- SpatialDimPlot(test,
                       cells.highlight = CellsByIdentities(object = test),
                       facet.highlight = TRUE,pt.size = spot_size)

Idents(test) <- "manual_segmentation3"
p227 <- SpatialDimPlot(test,
                       cells.highlight = CellsByIdentities(object = test),
                       facet.highlight = TRUE,pt.size = spot_size)

Idents(test) <- "SCT_snn_res.0.1"
p327 <- SpatialDimPlot(test,
                       cells.highlight = CellsByIdentities(object = test),
                       facet.highlight = TRUE,pt.size = spot_size)

(p17 + p17) + plot_annotation("BayesSpace",theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("../out/plot/101_SpatialDimPlot_BayesSpace_V07.pdf",width = 25,height = 10)
(p27 + p227) + plot_annotation("manual segmentation",theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("../out/plot/101_SpatialDimPlot_ManualSegmentation3_V07.pdf",width = 25,height = 10)
(p37 + p327) + plot_annotation("integrated clustering",theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("../out/plot/101_SpatialDimPlot_Integratedres0.1_V07.pdf",width = 25,height = 10)


# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# test
a <- c('potato', 'tomotto', 'chips', 'baloon')
b <- c('car', 'chips', 'bird', 'salt')

jaccard(a, b)

# build the dataset for the correlatino plot
df_crossing <- crossing(manual_segmentation = unique(meta_V07_update$manual_segmentation3),
                        integrated_clustering = unique(meta_V07_update$SCT_snn_res.0.1))
df_crossing

# build the scatter plot
df_jaccard_score <- pmap(list(id1 = df_crossing$manual_segmentation,
                              id2 = df_crossing$integrated_clustering), function(id1,id2){
                                
                                # calculate the jaccard score
                                a <- meta_V07_update %>%
                                  rownames_to_column("barcode") %>%
                                  filter(manual_segmentation3 == id1) %>% pull(barcode)
                                b <- meta_V07_update %>%
                                  rownames_to_column("barcode") %>%
                                  filter(SCT_snn_res.0.1 == id2) %>% pull(barcode)
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame(manual_segmentation = id1,
                                                 integrated_cluster = id2,
                                                 jaccard_score = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

head(df_jaccard_score)

# shape it as a matrix
mat_jaccard_score <- df_jaccard_score %>%
  pivot_wider(names_from = manual_segmentation,values_from = jaccard_score) %>%
  column_to_rownames("integrated_cluster")

mat_jaccard_score


# plot the matrix
ht_02 <- Heatmap(mat_jaccard_score,
                 name = "Jaccard score",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 100),
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize = 8),
                 column_names_side = "bottom",
                 column_names_gp = gpar(fontsize = 8),
                 row_dend_reorder = FALSE,
                 column_dend_reorder = FALSE,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 show_column_names = T,
                 show_row_names = T)

ht_02

df_crossing2 <- crossing(manual_segmentation = unique(meta_V07_update$manual_segmentation3),
                         bayespace = unique(meta_V07_update$BayesSpace))
df_crossing2

# build the scatter plot
df_jaccard_score2 <- pmap(list(id1 = df_crossing2$manual_segmentation,
                               id2 = df_crossing2$bayespace), function(id1,id2){
                                 
                                 # calculate the jaccard score
                                 a <- meta_V07_update %>%
                                   rownames_to_column("barcode") %>%
                                   filter(manual_segmentation3 == id1) %>% pull(barcode)
                                 b <- meta_V07_update %>%
                                   rownames_to_column("barcode") %>%
                                   filter(BayesSpace == id2) %>% pull(barcode)
                                 
                                 jaccard_score <- jaccard(a,b)
                                 
                                 # build a data.frame
                                 df <- data.frame(manual_segmentation = id1,
                                                  BayeSpace = id2,
                                                  jaccard_score = jaccard_score)
                                 return(df)
                               }) %>%
  bind_rows()

head(df_jaccard_score2)

# shape it as a matrix
mat_jaccard_score2 <- df_jaccard_score2 %>%
  pivot_wider(names_from = manual_segmentation,values_from = jaccard_score) %>%
  column_to_rownames("BayeSpace")

# plot the matrix
ht_12 <- Heatmap(mat_jaccard_score2,
                 name = "Jaccard score",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 100),
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize = 8),
                 column_names_side = "bottom",
                 column_names_gp = gpar(fontsize = 8),
                 row_dend_reorder = FALSE,
                 column_dend_reorder = FALSE,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 show_column_names = T,
                 show_row_names = T)

ht_12
