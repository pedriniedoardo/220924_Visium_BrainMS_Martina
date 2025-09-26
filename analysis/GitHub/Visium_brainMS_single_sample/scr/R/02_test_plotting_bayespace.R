# libraries ---------------------------------------------------------------
#Load libraries
library(tidyverse)
library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(BayesSpace)
library(patchwork)

# read in data ------------------------------------------------------------
# read in the processed data of the brain slices
list_brain <- readRDS(file = "out/object/list_brain_all_BayesSpace1000_q10.rds")
list_brain2 <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")

# generate an ad hoc plot for martina
test <- list_brain$V01

Idents(test) <- "BayesSpace"
SpatialDimPlot(test, cells.highlight = CellsByIdentities(object = test), facet.highlight = TRUE, ncol = 4,alpha = 0.5) +
  plot_annotation(title = "BayesSpace")

# plot the full UMAP
DimPlot(test,label = T)

# filter only the one shown in the shortlist panel
test_subset <- subset(test,subset = BayesSpace %in% c(1,4,9,7))

DimPlot(test_subset)
ggsave("out/image/test_subset.pdf",width = 5,height = 5)

test@meta.data

# plot markers ------------------------------------------------------------
# shortlist_features_list_long <- list(
#   IMMUNE = c("HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC","AIF1","HLA-DRA","TYROBP"),
#   B_CELLS = c("IGHG1", "CD38"),
#   T_CELLS =  c("SKAP1", "CD8A", "CD2"),
#   OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
#   ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
#   NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
#   ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR")
# )

shortlist_features_list_long <- list(
  IMMUNE = c("C3","CD74","C1QB","HLA-DRA"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS =  c("SKAP1", "CD8A", "CD2"),
  OLIGOLINEAGE = c("PLP1","OLIG1","OLIG2","MBP","MAG"),
  ASTRO = c("AQP4", "GFAP","VIM","APOE", "VCAN", "STAT3","SLC1A2","S100B"),
  NEURONS = c("THY1", "SLC17A7", "NRGN", "RORB", "SYP", "NEFL","SYT1"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","CLDN5")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
# source image
test_long_run03 <- DotPlot(test,
                           features = shortlist_features_list_long,
                           dot.scale = 10,
                           cluster.idents = T,
                           group.by = "BayesSpace",assay = "Spatial") +
  RotatedAxis() +
  labs(title = "BayesSpace")+
  theme(strip.text = element_text(angle = 90))

# generate the ggplot object
# force the order suggested by matrina
df_test_run03 <- lapply(shortlist_features_list_long,function(x){
  test_long_run03$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

# plot mapping radius
test_long_run03_alt <- df_test_run03 %>%
  # force the order
  mutate(id = factor(id,levels = c(4,5,9,1,6,3,7,8,2,10))) %>% 
  mutate(cell_type = factor(cell_type,levels = c("IMMUNE","B_CELLS","T_CELLS","NEURONS","ASTRO","OLIGOLINEAGE","ENDO"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90),strip.text = element_text(hjust = 0,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")


# jaccard score between clusters and manual annotation --------------------
meta_test <- list_brain2$V01@meta.data %>%
  rownames_to_column("barcode") %>%
  left_join(test@meta.data %>%
              rownames_to_column("barcode") %>%
              select(barcode,BayesSpace),by = "barcode")

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

#
df_crossing <- crossing(manual_annotation = unique(meta_test$manual_segmentation),
                        cluster_id = unique(meta_test$BayesSpace))
df_crossing

# build the scatter plot
# man_anno = "core"
# clu_id = "1"
df_jaccard_score <- pmap(list(anno = df_crossing$manual_annotation,
                              clu = df_crossing$cluster_id), function(anno,clu){
                                
                                # print(anno)
                                # calculate the jaccard score
                                a <- meta_test %>%
                                  filter(manual_segmentation == anno) %>% pull(barcode)
                                b <- meta_test %>%
                                  filter(BayesSpace == clu) %>% pull(barcode)
                                jaccard_score <- jaccard(a,b)

                                # build a data.frame
                                df <- data.frame(manual_annotation = anno,
                                                 cluster_id = clu,
                                                 jaccard_score = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

head(df_jaccard_score)

# shape it as a matrix
mat_jaccard_score <- df_jaccard_score %>%
  pivot_wider(names_from = manual_annotation,values_from = jaccard_score) %>%
  column_to_rownames("cluster_id")

mat_jaccard_score

# plot the matrix
ht_02 <- Heatmap(mat_jaccard_score,
                 name = "Jaccard score",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 20),
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
