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
library(ComplexHeatmap)
library(Matrix)
library(data.table)
library(gt)
library(SPOTlight)
library(RColorBrewer)
library(ggridges)
library(cowplot)
library(ggside)
library(circlize)
library(scales)

# read in the data --------------------------------------------------------
# read in the data already scored by senescence
list_brain <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE.rds")

# read in a list of new metadata
folder <- "data/segmentation/martina_new/"
file <- dir(folder)

tot_metadata <- lapply(file, function(x){
  read_csv(paste0(folder,x))
}) %>%
  setNames(file %>% str_remove_all(pattern = ",csv")) %>% 
  bind_rows(.id ="file") %>% 
  separate(file,into = c("id1","id2","slide"),remove = F) %>% 
  mutate(slide_fix = paste0("V",slide))

# is the annotation uniform
table(tot_metadata$manual_segmentation)

# resplit the table into slides
list_meta <- split(tot_metadata,f = tot_metadata$slide_fix)

# define the slides
names(list_brain)
names(list_meta)

# intersect the two informations
int_slides <- intersect(names(list_brain),names(list_meta))

# x <- "V01"
list_out <- lapply(int_slides,function(x){
  # set the brain object
  test <- list_brain[[x]]
  
  # read in the novel segmantation run by sofia
  test_meta <- list_meta[[x]] %>% 
    dplyr::select(-c("file","slide_fix","id1","id2","slide"))
  
  # check that the categories are coherent
  table(test_meta$manual_segmentation)
  dim(test_meta)
  
  # add the annotation provided by Sofia
  current_meta <- test@meta.data %>% 
    rownames_to_column("Barcode")
  
  dim(current_meta)
  
  # add the new annotation provide by sofia
  full_meta <- left_join(current_meta,test_meta,"Barcode") 
  dim(full_meta)
  
  # add the metadata to the original object. it is just one column
  test$manual_segmentation <- full_meta$manual_segmentation
  
  # plot the slide with the new segmentation
  plot_segmentation <- Seurat::SpatialDimPlot(
    object = test,group.by = "manual_segmentation",
    alpha = c(1))
  # ggsave("out/image/panel_V01_manual_segmentation.pdf",width = 6,height = 5)
  
  # try to check the gradints of all the different deconvolution scores. order the categories
  df_plot <- full_meta %>% 
    dplyr::select(Barcode,manual_segmentation,AST:VAS) %>% 
    dplyr::select(-contains("SCT_snn_res")) %>% 
    pivot_longer(names_to = "cell_type",values_to = "score",-c(Barcode,manual_segmentation))
  
  return(list(plot = plot_segmentation,
              df_plot = df_plot))
}) %>% 
  setNames(int_slides)

# check the colnames
lapply(list_out,function(x){
  table(x$df_plot$cell_type)
})

# extract all the meta and save them
tot_scores <- lapply(list_out,function(x){
  x$df_plot
}) %>% 
  bind_rows(.id = "slide")

tot_scores %>% 
  write_tsv("out/table/gradient_cellType_martina.tsv")

# make the summary of the gradients
tot_scores_summary <- tot_scores %>% 
  group_by(slide,manual_segmentation,cell_type) %>% 
  summarise(avg = mean(score),
            sd = sd(score),
            med = median(score))

tot_scores_summary %>% 
  write_tsv("out/table/gradient_cellType_summary_martina.tsv")

# make the plot
pdf("out/image/panel_all_manual_segmentation_gradient_cellType_unscaled_martina.pdf",width = 20,height = 20)
tot_scores %>% 
  split(tot_scores$slide) %>% 
  lapply(function(df_plot){
    name <- unique(df_plot$slide)
    df_plot %>% 
      ggplot(aes(x=manual_segmentation,y=score))+geom_violin()+
      geom_boxplot(width=0.1, outlier.shape = NA, position=position_dodge(0.75))+
      theme_bw()+
      theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
      facet_wrap(~cell_type)+
      geom_hline(yintercept = 0,linetype="dashed",alpha=0.5)+ggtitle(name)
  })
dev.off()

pdf("out/image/panel_all_manual_segmentation_gradient_cellType_scaled_martina.pdf",width = 20,height = 20)
tot_scores %>% 
  split(tot_scores$slide) %>% 
  lapply(function(df_plot){
    name <- unique(df_plot$slide)
    df_plot %>% 
      ggplot(aes(x=manual_segmentation,y=score))+geom_violin()+
      geom_boxplot(width=0.1, outlier.shape = NA, position=position_dodge(0.75))+
      theme_bw()+
      theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
      facet_wrap(~cell_type,scales = "free")+
      geom_hline(yintercept = 0,linetype="dashed",alpha=0.5)+ggtitle(name)
  })
dev.off()

# plotting test -----------------------------------------------------------
# load the summary and try to plot different situations
test <- read_tsv("out/table/gradient_cellType_summary_martina.tsv")
table(test$manual_segmentation)
table(test$manual_segmentation,test$slide)

# add some metaannotation
LUT_sample <- data.frame(slide = c("V01","V02","V03","V04","V05","V06","V08","V09","V10","V11","V12","V14"),
                         condition = c("CAL","CAL","remyel","non-lesional","active","CI","non-lesional","CI","CI","CI","non-lesional","non-lesional")) %>%
  mutate(annotation = paste(condition,slide))

# build a matrix for each signture
id_cellType <- unique(test$cell_type)
# signa <- "XIMERAKIS_ASC_fixed1"
list_mat <- lapply(id_cellType,function(ctype){
  mat <- left_join(test,LUT_sample,by = "slide") %>%
    mutate(area = manual_segmentation) %>% 
    # filer out the spots that do not have an area
    filter(!is.na(area)) %>% 
    dplyr::filter(cell_type == ctype) %>%
    dplyr::select(annotation,area,avg) %>%
    pivot_wider(names_from = area,values_from = avg) %>%
    # filter one sample
    # filter(annotation != "CI V07") %>%
    column_to_rownames("annotation") %>%
    # reorder the sample as recommended
    .[c("active V05","CAL V01","CAL V02","CI V06","CI V09","CI V10","CI V11","remyel V03","non-lesional V04","non-lesional V08","non-lesional V12","non-lesional V14"),] %>% 
    t() %>% 
    .[c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex"),]
  return(mat)
}) %>% 
  setNames(id_cellType)

# save the object
saveRDS(list_mat,"out/object/matrix_cellType_gradient_martina.rds")

show_col(c("blue",muted("blue")))

# plot the heatmaps
pdf("out/image/heatmap_gradient_cellType.pdf",width = 5,height = 4)
pmap(list(list_mat,names(list_mat)), function(x,name){
  Heatmap(as.matrix(x),cluster_columns = F,cluster_rows = F,col = viridis::plasma(10),name = "score",na_col = "white",column_title = name)
})
dev.off()

# try to map size to the average
p1 <- list_mat$IMM %>% 
  data.frame() %>% 
  rownames_to_column("area") %>% 
  pivot_longer(names_to = "tissue",values_to = "score",-area) %>% 
  filter(!is.na(score)) %>% 
  # order the y values
  mutate(area = factor(area,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex"))) %>% 
  # filter(ID %in% c("MICROGLIA")) %>%
  ggplot(aes(x=tissue,y=fct_rev(area),size = score,fill=score)) +
  geom_point(shape = 21) +
  scale_radius()+
  theme_cowplot()+
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('')+
  # scale_fill_gradientn(colours = viridis::viridis(20,option = "plasma"), limits = c(0,0.1), oob = scales::squish, name = 'score')+
  scale_fill_gradientn(colours = viridis::viridis(20,option = "plasma"), name = 'score')+
  theme(plot.margin = margin(2, 1, 1, 1, "cm"))+
  theme(legend.position = "right")

# pdf("out/image/test_gradient_astro_hm_dot_martina.pdf",width = 5,height = 4.5)
p1
# dev.off()

# plot the gratient with a line
list1 <- list_mat

lapply(list1,function(x){
  x %>% 
    data.frame() %>% 
    rownames_to_column("area") %>% 
    pivot_longer(names_to = "tissue",values_to = "score",-area) %>% 
    filter(!is.na(score)) %>% 
    # order the y values
    mutate(area = factor(area,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex"))) %>% 
    mutate(tissue = toupper(tissue)) %>% 
    mutate(tissue_class = str_extract(tissue,pattern = "ACTIVE|CAL|CI|REMYEL|NON.LESIONAL"))
}) %>% 
  bind_rows(.id = "cell_type") %>% 
  ggplot(aes(x=area,y=score,group=tissue))+geom_line()+facet_grid(cell_type~tissue_class,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/test_gradient_cellType_linePlot_martina.pdf",width = 10,height = 12)

# test plotting on the slide ----------------------------------------------
# list_lut <- readRDS("out/object/list_SIT.rds")
# slide <- "V01"

list_out <- lapply(int_slides,function(x){
  # set the brain object
  test <- list_brain[[x]]
  
  # read in the novel segmantation run by sofia
  test_meta <- list_meta[[x]] %>% 
    dplyr::select(-c("file","slide_fix","id1","id2","slide"))
  
  # check that the categories are coherent
  table(test_meta$manual_segmentation)
  dim(test_meta)
  
  # add the annotation provided by Sofia
  current_meta <- test@meta.data %>% 
    rownames_to_column("Barcode")
  
  dim(current_meta)
  
  # add the new annotation provide by sofia
  full_meta <- left_join(current_meta,test_meta,"Barcode") 
  dim(full_meta)
  
  # add the metadata to the original object. it is just one column
  test$manual_segmentation <- full_meta$manual_segmentation
  
  # plot the slide with the new segmentation
  plot_segmentation <- Seurat::SpatialDimPlot(
    object = test,group.by = "manual_segmentation",
    alpha = c(1))
  # ggsave("out/image/panel_V01_manual_segmentation.pdf",width = 6,height = 5)
  
  # try to check the gradints of all the different deconvolution scores. order the categories
  df_plot <- full_meta %>% 
    dplyr::select(Barcode,manual_segmentation,AST:VAS) %>% 
    dplyr::select(-contains("SCT_snn_res")) %>% 
    pivot_longer(names_to = "cell_type",values_to = "score",-c(Barcode,manual_segmentation))
  
  return(list(plot = plot_segmentation,
              df_plot = df_plot))
}) %>% 
  setNames(int_slides)


pdf("out/image/00_list_slide_cellType_gradient.pdf",width = 15,height = 15)
# run the analysis only on the common samples
lapply(int_slides,function(slide){
  # set the brain object
  test <- list_brain[[slide]]
  
  # read in the novel segmantation run by sofia
  test_meta <- list_meta[[slide]] %>% 
    dplyr::select(-c("file","slide_fix","id1","id2","slide"))
  
  # add the annotation provided by Sofia
  current_meta <- test@meta.data %>% 
    rownames_to_column("Barcode")
  
  # add the new annotation provide by sofia
  full_meta <- left_join(current_meta,test_meta,"Barcode") 
  dim(full_meta)
  
  # update the metadata
  test@meta.data <- full_meta %>% 
    column_to_rownames("Barcode")
  
  p1 <- Seurat::SpatialDimPlot(object = test,group.by = "manual_segmentation",alpha = 0)+ theme(legend.position = "none") + ggtitle(slide)
  p2 <- Seurat::SpatialDimPlot(object = test,group.by = "manual_segmentation",alpha = 1)
  # p1 <- SpatialFeaturePlot(brain, features = "nCount_Spatial",alpha=0,crop = F) + theme(legend.position = "none") + ggtitle(slide)
  p3 <- Seurat::SpatialFeaturePlot(object = test,features = id_cellType,alpha = 0.8)
  
  layout <- 
  "ACC
  BCC
  "
  p1 + p2 + p3 + plot_layout(design = layout)
})
dev.off()


# correlate the celltype score --------------------------------------------
# redefine the LUT of the slides
LUT_sample <- data.frame(slide = c("V01","V02","V03","V04","V05","V06","V08","V09","V10","V11","V12","V14"),
                         condition = c("CAL","CAL","remyel","non-lesional","active","CI","non-lesional","CI","CI","CI","non-lesional","non-lesional")) %>%
  mutate(annotation = paste(condition,slide))

# read in the score per gradient
df_cellType_gradient <- read_tsv("out/table/gradient_cellType_summary_martina.tsv")
df_senescence_gradient <- read_tsv("out/table/gradient_senescence_summary_martina.tsv")

# make it as a matrix
df_cellType_gradient2 <- df_cellType_gradient %>% 
  filter(!is.na(manual_segmentation)) %>% 
  mutate(row = paste0(slide,"_",manual_segmentation)) %>% 
  dplyr::select(row,cell_type,avg) %>% 
  pivot_wider(names_from = cell_type,values_from = avg) %>% 
  arrange(row)

df_senescence_gradient2 <- df_senescence_gradient %>% 
  filter(!is.na(manual_segmentation)) %>%
  mutate(row = paste0(slide,"_",manual_segmentation)) %>% 
  dplyr::select(row,signature,avg) %>% 
  pivot_wider(names_from = signature,values_from = avg) %>% 
  arrange(row)

# run the correlation matrix
mat_cor <- cor(df_cellType_gradient2 %>% 
                 column_to_rownames("row"),
               df_senescence_gradient2 %>%
                 column_to_rownames("row"),method = "spearman")

# plot the heatmap
pdf("out/image/00_corr_heatmap_gradient_ALL.pdf",width = 6,height = 6)
draw(Heatmap(mat_cor,column_title = "ALL"),heatmap_legend_side = "left",padding = unit(c(30, 2, 2, 2), "mm"))
dev.off()

# check some specific correlations
df_all <- left_join(df_cellType_gradient2,df_senescence_gradient2,by="row") %>% 
  separate(row,into = c("slide","area"),sep = "_",remove = F) %>% 
  left_join(LUT_sample,by = "slide")

mat_cor["AST","Inhibits1"]
cor(df_all$AST,df_all$Inhibits1)

df_all %>% 
  ggplot(aes(x=Inhibits1,y=AST))+geom_point()+geom_smooth(method = "lm")
df_all %>% 
  ggplot(aes(x=senmayo1,y=IMM))+geom_point()+geom_smooth(method = "lm")

# -------------------------------------------------------------------------
# check specific local correlations
id_sample <- df_all %>%
  filter(condition=="CAL") %>% 
  pull(row)


df_cellType_gradient3 <- df_cellType_gradient %>% 
  filter(!is.na(manual_segmentation)) %>% 
  mutate(row = paste0(slide,"_",manual_segmentation)) %>% 
  dplyr::select(row,cell_type,avg) %>% 
  pivot_wider(names_from = cell_type,values_from = avg) %>%
  dplyr::filter(row %in% id_sample) %>% 
  arrange(row)

df_senescence_gradient3 <- df_senescence_gradient %>% 
  filter(!is.na(manual_segmentation)) %>%
  mutate(row = paste0(slide,"_",manual_segmentation)) %>% 
  dplyr::select(row,signature,avg) %>% 
  pivot_wider(names_from = signature,values_from = avg) %>% 
  dplyr::filter(row %in% id_sample) %>% 
  arrange(row)

# run the correlation matrix
mat_cor3 <- cor(df_cellType_gradient3 %>% 
                 column_to_rownames("row"),
               df_senescence_gradient3 %>%
                 column_to_rownames("row"),method = "spearman")

# plot the heatmap
pdf("out/image/00_corr_heatmap_gradient_CAL.pdf",width = 6,height = 6)
draw(Heatmap(mat_cor3,column_title = "CAL"),heatmap_legend_side = "left",padding = unit(c(30, 2, 2, 2), "mm"))
dev.off()

# check some specific correlations
df_all3 <- left_join(df_cellType_gradient3,df_senescence_gradient3,by="row") %>% 
  separate(row,into = c("slide","area"),sep = "_",remove = F) %>% 
  left_join(LUT_sample,by = "slide")

mat_cor3["LYM","XIMERAKIS_PC_fixed1"]
cor(df_all3$LYM,df_all3$XIMERAKIS_PC_fixed1,method = "spearman")

df_all3 %>% 
  ggplot(aes(x=XIMERAKIS_PC_fixed1,y=LYM))+geom_point()+geom_smooth(method = "lm")

# correlate the prop of SIT -----------------------------------------------
# redefine the LUT of the slides
# LUT_sample <- data.frame(slide = c("V01","V02","V03","V04","V05","V06","V08","V09","V10","V11","V12","V14"),
#                          condition = c("CAL","CAL","remyel","non-lesional","active","CI","non-lesional","CI","CI","CI","non-lesional","non-lesional")) %>%
#   mutate(annotation = paste(condition,slide))

# read in the score per gradient
# df_cellType_gradient <- read_tsv("out/table/gradient_cellType_summary_martina.tsv")
df_senescence_gradient <- read_tsv("out/table/gradient_senescence_summary_martina.tsv")
df_SIT_gradient <- read_tsv("out/table/gradient_SIT_summary_martina.tsv")

# make it as a matrix
df_SIT_gradient2 <- df_SIT_gradient %>% 
  filter(!is.na(manual_segmentation)) %>%
  mutate(row = paste0(slide,"_",manual_segmentation)) %>% 
  dplyr::select(row,SENEQUANTILE,prop_senescence) %>% 
  pivot_wider(names_from = SENEQUANTILE,values_from = prop_senescence) %>% 
  arrange(row)

df_cellType_gradient2 <- df_cellType_gradient %>%
  filter(!is.na(manual_segmentation)) %>%
  mutate(row = paste0(slide,"_",manual_segmentation)) %>%
  dplyr::select(row,cell_type,avg) %>%
  pivot_wider(names_from = cell_type,values_from = avg) %>%
  filter(row %in% df_SIT_gradient2$row) %>% 
  arrange(row)

# not all the samples are present

# run the correlation matrix
mat_cor4 <- cor(df_cellType_gradient2 %>%
                  column_to_rownames("row"),
                df_SIT_gradient2 %>%
                  column_to_rownames("row"),method = "spearman")

# plot the heatmap
pdf("out/image/00_corr_heatmapSITcellType_gradient_ALL.pdf",width = 6,height = 6)
draw(Heatmap(mat_cor4,column_title = "ALL"),heatmap_legend_side = "left",padding = unit(c(30, 2, 2, 2), "mm"))
dev.off()

# read in the score per gradient
# df_cellType_gradient <- read_tsv("out/table/gradient_cellType_summary_martina.tsv")
# df_senescence_gradient <- read_tsv("out/table/gradient_senescence_summary_martina.tsv")
# df_SIT_gradient <- read_tsv("out/table/gradient_SIT_summary_martina.tsv")

# focus on the CAL sample only
id_sample <- df_all %>%
  filter(condition=="CAL") %>% 
  pull(row)

# make it as a matrix
df_SIT_gradient3 <- df_SIT_gradient %>% 
  filter(!is.na(manual_segmentation)) %>%
  mutate(row = paste0(slide,"_",manual_segmentation)) %>% 
  dplyr::select(row,SENEQUANTILE,prop_senescence) %>% 
  pivot_wider(names_from = SENEQUANTILE,values_from = prop_senescence) %>% 
  dplyr::filter(row %in% id_sample) %>% 
  arrange(row)

df_cellType_gradient3 <- df_cellType_gradient %>%
  filter(!is.na(manual_segmentation)) %>%
  mutate(row = paste0(slide,"_",manual_segmentation)) %>%
  dplyr::select(row,cell_type,avg) %>%
  pivot_wider(names_from = cell_type,values_from = avg) %>%
  filter(row %in% df_SIT_gradient2$row) %>% 
  dplyr::filter(row %in% id_sample) %>%
  arrange(row)

# not all the samples are present

# run the correlation matrix
mat_cor5 <- cor(df_cellType_gradient3 %>%
                  column_to_rownames("row"),
                df_SIT_gradient3 %>%
                  column_to_rownames("row"),method = "spearman")

# plot the heatmap
pdf("out/image/00_corr_heatmapSITcellType_gradient_CAL.pdf",width = 6,height = 6)
draw(Heatmap(mat_cor5,column_title = "CAL"),heatmap_legend_side = "left",padding = unit(c(30, 2, 2, 2), "mm"))
dev.off()