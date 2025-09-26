# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(hdf5r)
library(limma)
library(future)

# setup parallel ----------------------------------------------------------
# library(future)
# plan("multiprocess",workers=12)
# plan()

# read in the data --------------------------------------------------------
# img <- Read10X_Image(image.dir = "../../../raw_data/spaceranger_out/01_SP1_novaseq/outs/spatial/",
#                      image.name = "tissue_lowres_image.png",filter.matrix = F)
# spobj <- Load10X_Spatial("../../../raw_data/spaceranger_out/01_SP1_novaseq/outs/",filename = "raw_feature_bc_matrix.h5",image = img)

brain_all <- Load10X_Spatial("../../../raw_data/spaceranger_out/01_SP1_novaseq/outs/",filename = "raw_feature_bc_matrix.h5",filter.matrix = F)

brain_short <- Load10X_Spatial("../../../raw_data/spaceranger_out/01_SP1_novaseq/outs/",filename = "filtered_feature_bc_matrix.h5",filter.matrix = T)

p1 <- SpatialFeaturePlot(brain_all, features = "nCount_Spatial",alpha=0,crop = F) + theme(legend.position = "none")
p2 <- SpatialFeaturePlot(brain_short, features = "nCount_Spatial",crop = F) + theme(legend.position = "right")
p3 <- SpatialFeaturePlot(brain_all, features = "nCount_Spatial",crop = F) + theme(legend.position = "right")

wrap_plots(p1, p2, p3)
ggsave("out/image/test_outside.pdf",width = 15,height = 5)

# wrangling ---------------------------------------------------------------
# which one are the spot in the background
# one is in tissue and 0 is background add it to metadata
meta_image <- brain_all@images$slice1@coordinates %>%
  rownames_to_column("barcode")

# meta dataset
meta_seq <- brain_all@meta.data %>%
  rownames_to_column("barcode")

# plotting counts per spot
left_join(meta_seq,meta_image,"barcode") %>% 
  mutate(tissue_cat = case_when(tissue ==0~"background",
                                T~"tissue")) %>% 
  ggplot(aes(x=tissue_cat,y = nCount_Spatial))+geom_violin()+scale_y_log10()+geom_boxplot(width=0.1,outlier.shape = NA)+theme_bw()
