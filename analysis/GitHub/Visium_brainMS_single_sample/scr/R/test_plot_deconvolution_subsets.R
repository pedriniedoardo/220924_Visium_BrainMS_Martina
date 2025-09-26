df_tot_fix_IMM <- readRDS("out/object/df_tot_fix_IMM_dec.rds")
df_tot_fix_ASTRO <- readRDS("out/object/df_tot_fix_ASTRO_dec.rds")

# wrangling ---------------------------------------------------------------
df_full <- left_join(df_tot_fix_IMM %>% 
  dplyr::select(barcodes,dataset,class_top2,type,sofia_segmentation_short),
  df_tot_fix_ASTRO %>% 
  dplyr::select(barcodes,dataset,class_top3,type,sofia_segmentation_short),by = c("barcodes","dataset","type","sofia_segmentation_short")) %>% 
  mutate(class_final = case_when(class_top3=="IMM"~class_top2,
                                 class_top2=="AST"~class_top3,
                                 T~class_top2))

# heatmap 01 --------------------------------------------------------------
mat_120 <- df_full %>% 
  mutate(sofia_segmentation_short2 = case_when(sofia_segmentation_short%in% c("active lesion","core (RM)","core")~"core",
                                               T~sofia_segmentation_short)) %>% 
  dplyr::select(barcodes,dataset,class_final,type,sofia_segmentation_short2) %>% 
  group_by(sofia_segmentation_short2,class_final) %>% 
  summarise(count_spot = n()) %>% 
  ungroup() %>% 
  group_by(sofia_segmentation_short2) %>% 
  mutate(tot = sum(count_spot)) %>% 
  ungroup() %>% 
  mutate(prop = count_spot/tot) %>% 
  ungroup() %>% 
  group_by(sofia_segmentation_short2) %>% 
  mutate(norm = (prop-mean(prop))/sd(prop)) %>% 
  dplyr::select(sofia_segmentation_short2,class_final,norm) %>% 
  pivot_wider(names_from = class_final,values_from = norm) %>% 
  column_to_rownames("sofia_segmentation_short2")

Heatmap(mat_120,name = "Z score",col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")))

mat_122 <- df_full %>% 
  mutate(sofia_segmentation_short2 = case_when(sofia_segmentation_short%in% c("active lesion","core (RM)","core")~"core",
                                               T~sofia_segmentation_short)) %>% 
  dplyr::select(barcodes,dataset,class_final,type,sofia_segmentation_short2) %>% 
  group_by(sofia_segmentation_short2,class_final) %>% 
  summarise(count_spot = n()) %>% 
  ungroup() %>% 
  group_by(sofia_segmentation_short2) %>% 
  mutate(tot = sum(count_spot)) %>% 
  ungroup() %>% 
  mutate(prop = count_spot/tot) %>% 
  ungroup() %>% 
  group_by(class_final) %>% 
  mutate(norm = (prop-mean(prop))/sd(prop)) %>% 
  dplyr::select(sofia_segmentation_short2,class_final,norm) %>% 
  pivot_wider(names_from = class_final,values_from = norm) %>% 
  column_to_rownames("sofia_segmentation_short2")

library(circlize)
pdf("out/image/summary_SPOTLight_ASTRO_IMM_LocationCellType_heatmap.pdf",width = 5,height = 3)
Heatmap(mat_122,name = "Z score",col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")))
dev.off()

mat_123 <- df_full %>% 
  mutate(sofia_segmentation_short2 = case_when(sofia_segmentation_short%in% c("active lesion","core (RM)","core")~"core",
                                               T~sofia_segmentation_short)) %>% 
  dplyr::select(barcodes,dataset,class_final,type,sofia_segmentation_short2) %>% 
  group_by(sofia_segmentation_short2,class_final) %>% 
  summarise(count_spot = n()) %>% 
  ungroup() %>% 
  group_by(sofia_segmentation_short2) %>% 
  mutate(tot = sum(count_spot)) %>% 
  ungroup() %>% 
  mutate(prop = count_spot/tot) %>% 
  ungroup() %>% 
  group_by(class_final) %>% 
  mutate(norm = (prop-mean(prop))/sd(prop)) %>% 
  dplyr::select(sofia_segmentation_short2,class_final,norm) %>% 
  pivot_wider(names_from = sofia_segmentation_short2,values_from = norm) %>% 
  column_to_rownames("class_final")

library(circlize)
pdf("out/image/summary_SPOTLight_ASTRO_IMM_LocationCellType_heatmap2.pdf",width = 4,height = 5)
Heatmap(mat_123,name = "Z score",col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")))
dev.off()


# mat_123 <- df_full %>% 
#   dplyr::select(barcodes,dataset,class_final,type,sofia_segmentation_short) %>% 
#   group_by(sofia_segmentation_short,class_final) %>% 
#   summarise(count_spot = n()) %>% 
#   ungroup() %>% 
#   group_by(sofia_segmentation_short) %>% 
#   mutate(tot = sum(count_spot)) %>% 
#   ungroup() %>% 
#   mutate(prop = count_spot/tot) %>% 
#   ungroup() %>% 
#   dplyr::select(sofia_segmentation_short,class_final,prop) %>% 
#   pivot_wider(names_from = class_final,values_from = prop) %>% 
#   column_to_rownames("sofia_segmentation_short")
# Heatmap(mat_123,name = "prop")

# -------------------------------------------------------------------------
mat_02 <- df_full %>% 
  dplyr::select(barcodes,dataset,class_final,type,sofia_segmentation_short) %>% 
  group_by(type,class_final) %>% 
  summarise(count_spot = n()) %>% 
  ungroup() %>% 
  group_by(type) %>% 
  mutate(tot = sum(count_spot)) %>% 
  ungroup() %>% 
  mutate(prop = count_spot/tot) %>% 
  ungroup() %>% 
  group_by(class_final) %>% 
  mutate(norm = (prop-mean(prop))/sd(prop)) %>% 
  dplyr::select(type,class_final,norm) %>% 
  pivot_wider(names_from = class_final,values_from = norm) %>% 
  column_to_rownames("type")


# pdf("out/image/summary_SPOTLight_ASTRO_TypeCellType_heatmap.pdf",width = 5,height = 4)
# Heatmap(mat_01_filter,name = "prop",
#         top_annotation = column_ha)
Heatmap(mat_02,name = "Z score")
# dev.off()


