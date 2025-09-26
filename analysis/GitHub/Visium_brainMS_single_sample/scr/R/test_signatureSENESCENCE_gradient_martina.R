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
table(tot_metadata$manual_segmentation,tot_metadata$slide_fix)
table(tot_metadata$id1)
table(tot_metadata$slide_fix)

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
  
  # try to check the gradints of all the different signatures scored. order the categories
  df_plot <- full_meta %>% 
    dplyr::select(Barcode,CASELLA_UP_fixed1:manual_segmentation) %>% 
    dplyr::select(-contains("SCT_snn_res")) %>% 
    pivot_longer(names_to = "signature",values_to = "score",-c(Barcode,manual_segmentation))
  
  return(list(plot = plot_segmentation,
              df_plot = df_plot))
}) %>% 
  setNames(int_slides)

# check the colnames
lapply(list_out,function(x){
  table(x$df_plot$signature)
})

# extract and save the plots with sofia's annotation
pdf("out/image/panel_manual_segmentation_all_martina.pdf",width = 7,height = 6)
pmap(list(list_out,names(list_out)),function(x,name){
  x$plot+ggtitle(name)
})
dev.off()

# extract all the meta and save them
tot_scores <- lapply(list_out,function(x){
  x$df_plot
}) %>% 
  bind_rows(.id = "slide")

tot_scores %>% 
  write_tsv("out/table/gradient_senescence_martina.tsv")

# make the summary of the gradients
tot_scores_summary <- tot_scores %>% 
  group_by(slide,manual_segmentation,signature) %>% 
  summarise(avg = mean(score),
            sd = sd(score),
            med = median(score))

tot_scores_summary %>% 
  write_tsv("out/table/gradient_senescence_summary_martina.tsv")

# make the plot
pdf("out/image/panel_all_manual_segmentation_gradient_senescence_unscaled_martina.pdf",width = 20,height = 20)
tot_scores %>% 
  split(tot_scores$slide) %>% 
  lapply(function(df_plot){
    name <- unique(df_plot$slide)
    df_plot %>% 
      ggplot(aes(x=manual_segmentation,y=score))+geom_violin()+
      geom_boxplot(width=0.1, outlier.shape = NA, position=position_dodge(0.75))+
      theme_bw()+
      theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
      facet_wrap(~signature)+
      geom_hline(yintercept = 0,linetype="dashed",alpha=0.5)+ggtitle(name)
  })
dev.off()

pdf("out/image/panel_all_manual_segmentation_gradient_senescence_scaled_martina.pdf",width = 20,height = 20)
tot_scores %>% 
  split(tot_scores$slide) %>% 
  lapply(function(df_plot){
    name <- unique(df_plot$slide)
    df_plot %>% 
      ggplot(aes(x=manual_segmentation,y=score))+geom_violin()+
      geom_boxplot(width=0.1, outlier.shape = NA, position=position_dodge(0.75))+
      theme_bw()+
      theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
      facet_wrap(~signature,scales = "free")+
      geom_hline(yintercept = 0,linetype="dashed",alpha=0.5)+ggtitle(name)
  })
dev.off()

# plotting test -----------------------------------------------------------
# load the summary and try to plot different situations
test <- read_tsv("out/table/gradient_senescence_summary_martina.tsv")
table(test$manual_segmentation)
table(test$manual_segmentation,test$slide)

# add some metaannotation
LUT_sample <- data.frame(slide = c("V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V14"),
                  condition = c("CAL","CAL","remyel","non-lesional","active","CI","CI","non-lesional","CI","CI","CI","non-lesional","non-lesional")) %>%
  mutate(annotation = paste(condition,slide))

# build a matrix for each signture
id_signatures <- unique(test$signature)
# signa <- "XIMERAKIS_ASC_fixed1"
list_mat <- lapply(id_signatures,function(signa){
  mat <- left_join(test,LUT_sample,by = "slide") %>%
    mutate(area = manual_segmentation) %>% 
    # filer out the spots that do not have an area
    dplyr::filter(!is.na(area)) %>% 
    dplyr::filter(signature == signa) %>%
    dplyr::select(annotation,area,avg) %>%
    pivot_wider(names_from = area,values_from = avg) %>%
    # filter one sample
    # filter(annotation != "CI V07") %>%
    column_to_rownames("annotation") %>%
    # reorder the sample as recommended
    .[c("active V05","CAL V01","CAL V02","CI V06","CI V07","CI V09","CI V10","CI V11","remyel V03","non-lesional V04","non-lesional V08","non-lesional V12","non-lesional V14"),] %>% 
    t() %>% 
    .[c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex"),]
  return(mat)
}) %>% 
  setNames(id_signatures)

# save the object
saveRDS(list_mat,"out/object/matrix_senescence_gradient_martina.rds")

show_col(c("blue",muted("blue")))

# Heatmap(as.matrix(mat),cluster_columns = F,cluster_rows = F,col = viridis::turbo(10))
# Heatmap(as.matrix(mat_astro),cluster_columns = F,cluster_rows = F,col = viridis::plasma(10))
col_fun3 <- colorRamp2(breaks = seq(min(list_mat$XIMERAKIS_ASC_fixed1,na.rm = T), 0.1, length = 10), colors = viridis::plasma(10))
# col_fun4 <- colorRamp2(breaks = seq(min(mat_astro,na.rm = T), 0.1, length = 10), colors = viridis::plasma(10))

pdf("out/image/test_gradient_astro_hm_color_test_cap_martina.pdf",width = 4,height = 3)
Heatmap(as.matrix(list_mat$XIMERAKIS_ASC_fixed1),cluster_columns = F,cluster_rows = F,col = col_fun3,name = "score",na_col = "white")
dev.off()

pdf("out/image/test_gradient_astro_hm_color_martina.pdf",width = 4,height = 3)
Heatmap(as.matrix(list_mat$XIMERAKIS_ASC_fixed1),cluster_columns = F,cluster_rows = F,col = viridis::plasma(10),name = "score",na_col = "white")
dev.off()

pdf("out/image/test_gradient_endo_hm_color_martina.pdf",width = 4,height = 3)
Heatmap(as.matrix(list_mat$XIMERAKIS_EC_fixed1),cluster_columns = F,cluster_rows = F,col = viridis::plasma(10),name = "score",na_col = "white")
dev.off()

pdf("out/image/test_gradient_endo_hm_color2_martina.pdf",width = 4,height = 3)
Heatmap(as.matrix(list_mat$XIMERAKIS_EC_fixed1),cluster_columns = F,cluster_rows = F,col = viridis::plasma(10),name = "score",na_col = "white")
dev.off()

pdf("out/image/test_gradient_micro_hm_color_martina.pdf",width = 4,height = 3)
Heatmap(as.matrix(list_mat$XIMERAKIS_MG_fixed1),cluster_columns = F,cluster_rows = F,col = viridis::plasma(10),name = "score",na_col = "white")
dev.off()

pdf("out/image/test_gradient_micro_hm_color2_martina.pdf",width = 4,height = 3)
Heatmap(as.matrix(list_mat$XIMERAKIS_MG_fixed1),cluster_columns = F,cluster_rows = F,col = viridis::plasma(10),name = "score",na_col = "white")
dev.off()

pdf("out/image/test_gradient_opc_hm_color_martina.pdf",width = 4,height = 3)
Heatmap(as.matrix(list_mat$XIMERAKIS_OPC_fixed1),cluster_columns = F,cluster_rows = F,col = viridis::plasma(10),name = "score",na_col = "white")
dev.off()

pdf("out/image/test_gradient_opc_hm_color2_martina.pdf",width = 4,height = 3)
Heatmap(as.matrix(list_mat$XIMERAKIS_OPC_fixed1),cluster_columns = F,cluster_rows = F,col = viridis::plasma(10),name = "score",na_col = "white")
dev.off()

# Heatmap(as.matrix(mat),cluster_columns = F,cluster_rows = F,col = colorRamp2(c(0, 0.1), c("lightgray", muted("red"))),na_col = "white")

# try to map size to the average
p1 <- list_mat$XIMERAKIS_ASC_fixed1 %>% 
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
  scale_fill_gradientn(colours = viridis::viridis(20,option = "plasma"), limits = c(0,0.1), oob = scales::squish, name = 'score')+
  theme(plot.margin = margin(2, 1, 1, 1, "cm"))+
  theme(legend.position = "right")

pdf("out/image/test_gradient_astro_hm_dot_martina.pdf",width = 5,height = 4.5)
p1
dev.off()

p2 <- list_mat$XIMERAKIS_EC_fixed1 %>% 
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
  scale_fill_gradientn(colours = viridis::viridis(20,option = "plasma"))+
  theme(plot.margin = margin(2, 1, 1, 1, "cm"))+
  theme(legend.position = "right")

pdf("out/image/test_gradient_endo2_hm_dot_martina.pdf",width = 5,height = 4.5)
p2
dev.off()

p3 <- list_mat$XIMERAKIS_MG_fixed1 %>% 
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
  scale_fill_gradientn(colours = viridis::viridis(20,option = "plasma"))+
  theme(plot.margin = margin(2, 1, 1, 1, "cm"))+
  theme(legend.position = "right")

pdf("out/image/test_gradient_micro2_hm_dot_martina.pdf",width = 5,height = 4.5)
p3
dev.off()

p4 <- list_mat$XIMERAKIS_OPC_fixed1 %>% 
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
  scale_fill_gradientn(colours = viridis::viridis(20,option = "plasma"))+
  theme(plot.margin = margin(2, 1, 1, 1, "cm"))+
  theme(legend.position = "right")

pdf("out/image/test_gradient_opc2_hm_dot_martina.pdf",width = 5,height = 4.5)
p4
dev.off()

# plot lines
list_mat$XIMERAKIS_ASC_fixed1 %>% 
  data.frame() %>%
  rownames_to_column("area") %>% 
  pivot_longer(names_to = "tissue",values_to = "score",-area) %>% 
  filter(!is.na(score)) %>% 
  # order the y values
  mutate(area = factor(area,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex"))) %>% 
  mutate(tissue_class = str_extract(tissue,pattern = "active|CAL|CI|remyel|non.lesional")) %>% 
  ggplot(aes(x=area,y=score,group=tissue))+geom_line()+facet_grid(~tissue_class)+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/test_gradient_astro_linePlot_martina.pdf",width = 10,height = 3)

list_mat$XIMERAKIS_ASC_fixed1 %>% 
  data.frame() %>%
  rownames_to_column("area") %>% 
  pivot_longer(names_to = "tissue",values_to = "score",-area) %>% 
  filter(!is.na(score)) %>% 
  # order the y values
  mutate(area = factor(area,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex"))) %>% 
  mutate(tissue_class = str_extract(tissue,pattern = "active|CAL|CI|remyel|non.lesional")) %>% 
  ggplot(aes(x=area,y=score,group=tissue))+geom_line()+
  # facet_wrap(~tissue_class,ncol = 1,scales = "free_x")+
  facet_wrap(~tissue_class,ncol = 1)+
  theme_cowplot()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/test_gradient_astro_linePlot_martina2.pdf",width = 4,height = 10)

# save the table
list_mat$XIMERAKIS_ASC_fixed1 %>% 
  data.frame() %>%
  rownames_to_column("area") %>% 
  pivot_longer(names_to = "tissue",values_to = "score",-area) %>% 
  filter(!is.na(score)) %>% 
  # order the y values
  mutate(area = factor(area,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex"))) %>% 
  mutate(tissue_class = str_extract(tissue,pattern = "active|CAL|CI|remyel|non.lesional")) %>%
  write_tsv("out/table/gradient_XIMERAKIS_ASC_fixed1_martina.tsv")

# 
# list1 <- list(mat_astro = list_mat$XIMERAKIS_ASC_fixed1,
#               mat_endo = list_mat$XIMERAKIS_EC_fixed1,
#               mat_micro = list_mat$XIMERAKIS_MG_fixed1,
#               mat_opc = list_mat$XIMERAKIS_OPC_fixed1)
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
  bind_rows(.id = "signature") %>% 
  ggplot(aes(x=area,y=score,group=tissue))+geom_line()+facet_grid(signature~tissue_class,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/test_gradient_linePlot_martina.pdf",width = 10,height = 35)


list2 <- list(mat_astro = list_mat$XIMERAKIS_ASC_fixed1,
              mat_endo = list_mat$XIMERAKIS_EC_fixed1,
              mat_micro = list_mat$XIMERAKIS_MG_fixed1,
              mat_opc = list_mat$XIMERAKIS_OPC_fixed1)

# lapply(list2,function(x){
#   x %>% 
#     data.frame() %>% 
#     rownames_to_column("area") %>% 
#     pivot_longer(names_to = "tissue",values_to = "score",-area) %>% 
#     filter(!is.na(score)) %>% 
#     # order the y values
#     mutate(area = factor(area,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex"))) %>% 
#     mutate(tissue = toupper(tissue)) %>% 
#     mutate(tissue_class = str_extract(tissue,pattern = "ACTIVE|CAL|CI|REMYEL|NON.LESIONAL"))
# }) %>% 
#   bind_rows(.id = "signature") %>% 
#   ggplot(aes(x=area,y=score,group=tissue))+geom_line()+facet_grid(signature~tissue_class,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("out/image/test_gradient_linePlot2.pdf",width = 10,height = 12)

# test plotting on the slide ----------------------------------------------
list_brain <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE.rds")

# pull the medatada from the median estimate of the V01
test <- list_brain$V01

# read in the segmantation provided by sofia
# read in a list of new metadata
folder <- "data/segmentation/martina_new/"
file <- dir(folder)

tot_metadata <- lapply(file, function(x){
  read_csv(paste0(folder,x))
}) %>%
  setNames(file %>% str_remove_all(pattern = ".csv")) %>% 
  bind_rows(.id ="file") %>% 
  separate(file,into = c("id1","id2","slide"),remove = F) %>% 
  mutate(slide_fix = paste0("V",slide))

test_meta <- tot_metadata %>% 
  filter(slide_fix == "V01")

# read in the estimates summarized of the senesncece values
LUT_sigantures <- read_tsv("out/table/gradient_senescence_summary_martina.tsv")
test_meta_senescence <- LUT_sigantures %>% 
  mutate(signature = paste0(signature,".","avg")) %>% 
  filter(slide == "V01") %>% 
  dplyr::select(manual_segmentation,signature,avg) %>% 
  pivot_wider(names_from = signature,values_from = avg)

# join the tables
test_meta_senescence_full <- left_join(test_meta,test_meta_senescence,by="manual_segmentation")

meta_test <- test@meta.data %>% 
  rownames_to_column("Barcode")

# join the full table with all the annotations
meta_test_full <- left_join(meta_test,test_meta_senescence_full,"Barcode") %>% 
  column_to_rownames("Barcode")

# swap the metadata
test@meta.data <- meta_test_full

SpatialDimPlot(test,group.by = "manual_segmentation")+
  # SpatialDimPlot(test,group.by = "XIMERAKIS_ASC_fixed1.avg")+scale_fill_gradient()
  (SpatialFeaturePlot(test,features = "XIMERAKIS_ASC_fixed1.avg")+scale_fill_gradientn(colours = viridis::viridis(20,option = "plasma")))
ggsave("out/image/test_panel_gradient_astro_V01_martina.pdf",width = 10,height = 4)

SpatialDimPlot(test,group.by = "manual_segmentation")+
  (SpatialFeaturePlot(test,features = "XIMERAKIS_ASC_fixed1.avg")+scale_fill_gradientn(colours = viridis::viridis(20,option = "turbo")))
ggsave("out/image/test_panel_gradient_astro_V01_2_martina.pdf",width = 10,height = 4)

SpatialDimPlot(test,group.by = "manual_segmentation")+
  (SpatialFeaturePlot(test,features = "XIMERAKIS_ASC_fixed1.avg"))
ggsave("out/image/test_panel_gradient_astro_V01_3_martina.pdf",width = 10,height = 4)

# plot the gradients of cell type -----------------------------------------


