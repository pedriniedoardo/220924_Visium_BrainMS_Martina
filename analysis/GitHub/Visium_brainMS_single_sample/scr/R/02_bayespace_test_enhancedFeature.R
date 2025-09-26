# libraries ---------------------------------------------------------------
#Load libraries
library(tidyverse)
library(Seurat)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(BayesSpace)

# read in the seurat object already processed -----------------------------
list_brain <- readRDS(file = "out/object/list_brain_all_BayesSpace1000_q05.rds")

# try enhanced resolution -------------------------------------------------
# select just on one slice
brain <- list_brain$V01
# 
# # Convert to SCE
# diet.seurat <- Seurat::DietSeurat(brain, graphs = "pca") #slim down Seurat obj prior to conversion
# sce <- as.SingleCellExperiment(diet.seurat) #convert seurat to SCE
# colData(sce) <- cbind(colData(sce), brain@images$slice1@coordinates) #add spatial info to SCE
# 
# # BayesSpace Workflow
# sce <- spatialPreprocess(sce, platform = "Visium", skip.PCA = T, log.normalize = F) #add BayesSpace metadata, without messing with PCA/logcounts
# # sce <- qTune(sce, qs=seq(2, 10), platform="Visium")
# # qPlot(sce)
# 
# sce <- spatialCluster(sce, nrep = 1000, burn.in = 100, q = 10) #quickly cluster via BayesSpace
# 
# # enhance spots
# sce.enhanced <- spatialEnhance(sce, q=10, platform="Visium",
#                                nrep=2000, gamma=3, verbose=TRUE,
#                                jitter_scale=5.5, jitter_prior=0.3,
#                                save.chain=TRUE,burn.in = 100)
# 
# # enhance the expresison of the features focus on the HVG
markers <- brain@assays$SCT@var.features
# sce.enhanced <- enhanceFeatures(sce.enhanced, sce,
#                                 feature_names=markers,
#                                 nrounds=0)
# save the sample input and enhanced object
# saveRDS(sce.enhanced,"out/object/test_sce.enhanced_V01.rds")
# saveRDS(sce,"out/object/test_sce_V01.rds")

sce.enhanced <- readRDS("out/object/test_sce.enhanced_V01.rds")
sce <- readRDS("out/object/test_sce_V01.rds")

# explore the enhanced object ---------------------------------------------
# By default, log-normalized expression (logcounts(sce)) is imputed, although other assays or arbitrary feature matrices can be specified.
logcounts(sce.enhanced)[markers, 1:5]

# notice the dimension of the expression matrix for the enhanced version
dim(sce.enhanced@assays@data$logcounts)
# compared with the non enhanced dataste
sce@assays@data$logcounts[1:10,1:10]
dim(sce@assays@data$logcounts)

# the main difference is the fact that the enhanced is not filling the values for the genes not asked for the enhancement. only the one specified
sce.enhanced@assays@data$logcounts[markers[1:10],1:10]

# Diagnostic measures from each predictive model, such as rmse when using xgboost, are added to the rowData of the enhanced dataset.
rowData(sce.enhanced)[markers, ]

# Spatial gene expression is visualized with featurePlot().
# featurePlot(sce.enhanced, "SPP1")|featurePlot(sce, "SPP1")
# featurePlot(sce.enhanced, "MBP")|featurePlot(sce, "MBP")

# Here, we compare the spatial expression of the imputed marker genes.
enhanced.plots <- purrr::map(markers[1:10], function(x) featurePlot(sce.enhanced, x))
# patchwork::wrap_plots(enhanced.plots, ncol=2)

# And we can compare to the spot-level expression.
spot.plots <- purrr::map(markers[1:10], function(x) featurePlot(sce, x))
# patchwork::wrap_plots(c(enhanced.plots, spot.plots), ncol=4)

# compare the identity of the spots after enhancement ---------------------
# fill a set of features
# head(colData(sce.enhanced))
dim(colData(sce.enhanced))
dim(colData(sce))

# I can use the columns spot.row and spot.col in the enhanced file to join the info from sce
colData(sce.enhanced) %>% 
  data.frame() %>% 
  group_by(spot.row,spot.col) %>% 
  summarise(n = n())

# join the two datasets
df_meta <- left_join(data.frame(colData(sce.enhanced)),data.frame(colData(sce)),by=c("spot.row"="row","spot.col"="col"),suffix = c(".enhanced",".ref"))

# is the spatial cluster from the enhanced spot the same of the original bayespace cluster from the reference object?
# df_meta %>%
#   mutate(spatial.cluster.enhanced = factor(spatial.cluster.enhanced),
#          BayesSpace = factor(BayesSpace)) %>%
#   group_by(BayesSpace,spatial.cluster.enhanced) %>% 
#   summarise(count = n()) %>% 
#   mutate(prop = count/sum(count)) %>% 
#   # ggplot(aes(x = BayesSpace,fill = spatial.cluster.enhanced)) + geom_bar()
#   ggplot(aes(x = BayesSpace,y=prop,fill = spatial.cluster.enhanced)) + geom_col()+theme_bw()

# df_meta %>%
#   mutate(spatial.cluster.enhanced = factor(spatial.cluster.enhanced),
#          BayesSpace = factor(BayesSpace)) %>%
#   ggplot(aes(x = BayesSpace,fill = spatial.cluster.enhanced)) + geom_bar() + theme_bw()

df_meta %>%
  mutate(spatial.cluster.enhanced = factor(spatial.cluster.enhanced),
         BayesSpace = factor(BayesSpace)) %>%
  ggplot(aes(x = BayesSpace,fill = spatial.cluster.enhanced)) + geom_bar(position = "fill") + theme_bw()
# the same spot gets more classifications

# modify the plotting -----------------------------------------------------
# compare the plotting fo the same feature
# regular option
(clusterPlot(sce,label = "BayesSpace")+ggtitle("V01_sce")) | (clusterPlot(sce.enhanced,label = "spatial.cluster")+ggtitle("V01_sce.enhanced"))
(featurePlot(sce, "SPP1")+ggtitle("V01_sce")) | (featurePlot(sce.enhanced, "SPP1")+ggtitle("V01_sce.enhanced"))
(featurePlot(sce.enhanced, "SPP1")+ggtitle("V01_sce.enhanced"))
ggsave("out/image/SPP1_enhanced.pdf",width = 6,height = 5)

(featurePlot(sce.enhanced, "SERPINA3")+ggtitle("V01_sce.enhanced"))

# try the same plots with different layers
(clusterPlot(sce, color="gray",label = "BayesSpace") +
    theme_bw() +
    labs(fill="BayesSpace\ncluster", title="V01_sce")) |
  (clusterPlot(sce.enhanced, color="gray",label = "spatial.cluster") +
     theme_bw() +
     labs(fill="enhanced\ncluster", title="V01_sce.enhanced"))

# try to change the shape of the unit -------------------------------------
test.ref_cluster <- clusterPlot(sce,label = "BayesSpace")
test.enhanced_cluster <- clusterPlot(sce.enhanced,label = "spatial.cluster")

# notice that even the ref object has more rows than spot because it has also the coordinated of the hexagon
test.ref_cluster$data %>% 
  group_by(x.pos,y.pos,spot,fill) %>% 
  summarise(n = n())

# in the enhanced one there are even more spots
test.enhanced_cluster$data %>% 
  # group_by(x.pos,y.pos,spot,fill,x.offset,y.offset) %>% 
  group_by(x.pos,y.pos,spot,fill) %>% 
  summarise(n = n())

# trimm the poin to the correct number
plot_ref_cluster <-test.ref_cluster$data %>% 
  arrange(x.pos,y.pos,spot) %>% 
  group_by(x.pos,y.pos,spot,fill) %>% 
  summarise(n = n()) %>% 
  # dim()
  mutate(fill = factor(fill,levels = 1:10))

p10 <- plot_ref_cluster %>%
  ggplot(aes(x=x.pos,y=y.pos)) +
  # geom_point(color='black',aes(fill=fill),shape=21) +
  geom_point(aes(col=fill)) +
  scale_y_continuous(trans = "reverse") +
  coord_fixed() +
  theme_void()
  

clusterPlot(sce,label = "BayesSpace") +
  theme_bw()

# trimm the points to the correct number
test.ref_cluster$data %>% 
  arrange(x.pos,y.pos,spot) %>% head(6)

# # according to the data, the radius of the spot is 1 in fact if x offset is 
# 0.5/sin(30*pi/180)
# # this is the radius of the section of the exagon
offset_x <- 0.3/cos(30*pi/180)
offset_y <- sqrt(offset_x^2-(offset_x/2)^2)

# what is the distance between two consicutive dots
test.ref_cluster$data %>% 
  group_by(x.pos,y.pos) %>% 
  summarise(n = n())

# # check to order of the subspot
# clusterPlot(sce.enhanced, color="gray",label = "subspot.idx") +
#   theme_bw() +
#   labs(fill="enhanced\ncluster", title="V01_sce.enhanced")

# fix the counting of the spots
plot_enhanced_cluster <- test.enhanced_cluster$data %>%
  separate(spot,into = c("spot","subspot"),sep = "\\.") %>% 
  # arrange(x.pos,y.pos,spot) %>% head(18)
  mutate(x.pos_2 = case_when(subspot == 6~(x.pos-offset_x),
                             subspot == 2~(x.pos-offset_x/2),
                             subspot == 1~(x.pos+offset_x/2),
                             subspot == 5~(x.pos+offset_x),
                             subspot == 3~(x.pos+offset_x/2),
                             subspot == 4~(x.pos-offset_x/2)
                             )) %>%
  mutate(y.pos_2 = case_when(subspot == 6~(y.pos),
                             subspot == 2~(y.pos+offset_y),
                             subspot == 1~(y.pos+offset_y),
                             subspot == 5~(y.pos),
                             subspot == 3~(y.pos-offset_y),
                             subspot == 4~(y.pos-offset_y)
                             )) %>%
  group_by(x.pos_2,y.pos_2,spot,fill,subspot) %>%
  # group_by(spot,fill) %>%
  summarise(n = n()) %>%
  arrange(spot) %>%
  mutate(fill = factor(fill,levels = 1:10))
# chack the dimensitons are correct
head(plot_enhanced_cluster)
dim(plot_enhanced_cluster)

p20 <- plot_enhanced_cluster %>%
  ggplot(aes(x=x.pos_2,y=y.pos_2)) +
  # geom_point(color='black',aes(fill=fill),shape=21,size = 1) +
  # geom_point(aes(col=subspot),size = 1) +
  geom_point(aes(col=fill),size = 1) +
  scale_y_continuous(trans = "reverse") +
  coord_fixed() +
  theme_void()

# clusterPlot(sce.enhanced,label = "spatial.cluster") +
#   theme_bw()

p10+ggtitle("V01_sce")|p20+ggtitle("V01_sce.enhanced")
(clusterPlot(sce,label = "BayesSpace")+ggtitle("V01_sce")) | (clusterPlot(sce.enhanced,label = "spatial.cluster")+ggtitle("V01_sce.enhanced"))

# do the same as above for the feature expression -------------------------
test.ref_feature <- featurePlot(sce, "SPP1")
test.enhanced_feature <- featurePlot(sce.enhanced, "SPP1")
list_plot <- map(c("SPP1","CHIT1","FTL","APOE","FTH1","C1QB","C5orf63"),.f = ~featurePlot(sce.enhanced,.x))
pdf("out/image/test_enhanced",width = 7,height = 6)
list_plot
dev.off()

featurePlot(sce.enhanced, "TSPO")

featurePlot(sce.enhanced,"C5orf63",is.enhanced = T)
featurePlot(sce.enhanced,"SPP1",is.enhanced = T)
featurePlot(sce,"C5orf63")

# trimm the poin to the correct number
plot_ref_feature <-test.ref_feature$data %>% 
  arrange(x.pos,y.pos,spot) %>% 
  group_by(x.pos,y.pos,spot,fill) %>% 
  summarise(n = n())
  # dim()
  # mutate(fill = factor(fill,levels = 1:10))

p11 <- plot_ref_feature %>%
  ggplot(aes(x=x.pos,y=y.pos)) +
  # geom_point(color='black',aes(fill=fill),shape=21) +
  geom_point(aes(col=fill)) +
  scale_y_continuous(trans = "reverse") +
  theme_void() +
  coord_fixed() +
  scale_color_gradient(low = "gray",high = "red")

# # according to the data, the radius of the spot is 1 in fact if x offset is 
# 0.5/sin(30*pi/180)
# # this is the radius of the section of the exagon
offset_x <- 0.3/cos(30*pi/180)
offset_y <- sqrt(offset_x^2-(offset_x/2)^2)

# # check to order of the subspot
# clusterPlot(sce.enhanced, color="gray",label = "subspot.idx") +
#   theme_bw() +
#   labs(fill="enhanced\ncluster", title="V01_sce.enhanced")

# fix the counting of the spots
plot_enhanced_feature <- test.enhanced_feature$data %>%
  separate(spot,into = c("spot","subspot"),sep = "\\.") %>% 
  # arrange(x.pos,y.pos,spot) %>% head(18)
  mutate(x.pos_2 = case_when(subspot == 6~(x.pos-offset_x),
                             subspot == 2~(x.pos-offset_x/2),
                             subspot == 1~(x.pos+offset_x/2),
                             subspot == 5~(x.pos+offset_x),
                             subspot == 3~(x.pos+offset_x/2),
                             subspot == 4~(x.pos-offset_x/2)
  )) %>%
  mutate(y.pos_2 = case_when(subspot == 6~(y.pos),
                             subspot == 2~(y.pos+offset_y),
                             subspot == 1~(y.pos+offset_y),
                             subspot == 5~(y.pos),
                             subspot == 3~(y.pos-offset_y),
                             subspot == 4~(y.pos-offset_y)
  )) %>%
  group_by(x.pos_2,y.pos_2,spot,fill,subspot) %>%
  # group_by(spot,fill) %>%
  summarise(n = n()) %>%
  arrange(spot)
  # mutate(fill = factor(fill,levels = 1:10))
# chack the dimensitons are correct
head(plot_enhanced_cluster)
dim(plot_enhanced_cluster)

p21 <- plot_enhanced_feature %>%
  ggplot(aes(x=x.pos_2,y=y.pos_2)) +
  # geom_point(color='black',aes(fill=fill),shape=21,size = 1) +
  # geom_point(aes(col=subspot),size = 1) +
  geom_point(aes(col=fill),size = 1) +
  scale_y_continuous(trans = "reverse")+coord_fixed()+
  theme_void() +
  scale_color_gradient(low = "gray",high = "red")

p11+ggtitle("V01_sce")|p21+ggtitle("V01_sce.enhanced")
(test.ref_feature+ggtitle("V01_sce")) | (test.enhanced_feature+ggtitle("V01_sce.enhanced"))
