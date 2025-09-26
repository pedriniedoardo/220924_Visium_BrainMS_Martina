brain <- readRDS(paste0("../out_large/brain_V05.rds"))
dim(df_exp)

DefaultAssay(brain) <- "SCT"
SpatialFeaturePlot(brain, features = c("TSPO"))
test_slot$SCT %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% "TSPO")

SpatialFeaturePlot(brain, features = c("TSPO"),slot = "scale.data")

brain@assays$Spatial
brain@assays$SCT@counts

DefaultAssay(brain) <- "Spatial"
test_slot$Spatial %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% "TSPO")
SpatialFeaturePlot(brain, features = c("TSPO"))
SpatialFeaturePlot(brain, features = c("TSPO"),slot = "scale.data")

?SpatialFeaturePlot

Idents(brain)
brain@meta.data


brain@assays$SCT@counts
brain@assays$SCT@data
brain@assays$SCT@data["TSPO",]

test <- left_join(brain@meta.data %>% rownames_to_column("barcodes"),df_sofi %>% 
                    filter(dataset == "brain_V05")) %>% 
  column_to_rownames("barcodes")

brain@meta.data <- test
Idents(brain) <- "sofia_segmentation"
test_slot <- AverageExpression(brain)
str(test_slot)
test_slot$Spatial %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% "TSPO")

test_slot$SCT %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% "TSPO")

df_sofi %>% 
  filter(dataset == "brain_V05")

test_slot2 <- AverageExpression(brain,slot ="scale.data")

test_slot2$SCT %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% "TSPO")

test_slot2 %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% "TSPO")