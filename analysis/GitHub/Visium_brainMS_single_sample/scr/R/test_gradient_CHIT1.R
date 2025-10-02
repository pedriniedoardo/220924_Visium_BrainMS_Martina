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
library(broom)

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

# check if CHIT1 is expressed or not
df_gene_exp <- lapply(list_brain, function(x){
  data.frame(CHIT1_present = "CHIT1" %in% rownames(x))
}) %>%
  bind_rows(.id = "slide")

slide_CHIT1 <- df_gene_exp %>%
  filter(CHIT1_present == T) %>%
  pull(slide)

# keep only the slides that have CHIT1
int_slides2 <- intersect(int_slides,slide_CHIT1)

# x <- "V14"
list_out <- lapply(int_slides2,function(x){
  # track the progress
  print(x)
  
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
  
  # add the expression of CHIT1. use the normaized espression
  dfCHIT11 <- FetchData(test,"CHIT1",slot="data") %>%
    dplyr::rename("CHIT1_norm" = "CHIT1") %>%
    rownames_to_column("Barcode")
  
  dfCHIT12 <- FetchData(test,"CHIT1",slot="count") %>%
    dplyr::rename("CHIT1_raw" = "CHIT1") %>%
    rownames_to_column("Barcode")
  
  # add the new annotation provide by sofia
  full_meta <- left_join(current_meta,test_meta,"Barcode") %>%
    left_join(dfCHIT11,"Barcode") %>%
    left_join(dfCHIT12,"Barcode")
  dim(full_meta)
  
  # add the metadata to the original object. it is just one column
  test$manual_segmentation <- full_meta$manual_segmentation
  test$CHIT1_raw <- full_meta$CHIT1_raw
  test$CHIT1_norm <- full_meta$CHIT1_norm
  
  # plot the slide with the new segmentation
  plot_segmentation <- Seurat::SpatialDimPlot(
    object = test,group.by = "manual_segmentation",
    alpha = c(1))
  # ggsave("out/image/panel_V01_manual_segmentation.pdf",width = 6,height = 5)
  
  # try to check the gradints of all the different signatures scored. order the categories
  df_plot <- full_meta %>% 
    dplyr::select(Barcode,manual_segmentation,CHIT1_norm,CHIT1_raw) %>% 
    pivot_longer(names_to = "CHIT1",values_to = "exp",-c(Barcode,manual_segmentation))
  
  return(list(plot = plot_segmentation,
              df_plot = df_plot))
}) %>% 
  setNames(int_slides2)

# extract all the meta and save them
tot_scores <- lapply(list_out,function(x){
  x$df_plot
}) %>% 
  bind_rows(.id = "slide")

tot_scores %>% 
  write_tsv("out/table/gradient_CHIT1_martina.tsv")

# make the summary of the gradients
tot_scores_summary <- tot_scores %>% 
  group_by(slide,manual_segmentation,CHIT1) %>% 
  summarise(avg = mean(exp),
            sd = sd(exp),
            med = median(exp))

tot_scores_summary %>% 
  write_tsv("out/table/gradient_CHIT1_summary_martina.tsv")

# plotting test -----------------------------------------------------------
# load the summary and try to plot different situations
test <- read_tsv("out/table/gradient_CHIT1_summary_martina.tsv")
table(test$manual_segmentation)
table(test$manual_segmentation,test$slide)

# add some metaannotation
LUT_sample <- data.frame(slide = c("V01","V02","V03","V04","V05","V06","V08","V09","V10","V11","V12","V14"),
                         condition = c("CAL","CAL","remyel","non-lesional","active","CI","non-lesional","CI","CI","CI","non-lesional","non-lesional")) %>%
  mutate(annotation = paste(condition,slide))

df_test <- test %>%
  left_join(LUT_sample,by = "slide") %>%
  filter(!is.na(manual_segmentation)) %>% 
  # order the y values
  mutate(area = factor(manual_segmentation,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex")))%>%
  mutate(tissue_class = str_extract(condition,pattern = "active|CAL|CI|remyel|non-lesional")) %>%
  mutate(tissue = toupper(tissue_class))

df_test %>%
  ggplot(aes(x=area,y=avg,group=slide))+geom_line()+facet_grid(CHIT1~tissue,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/test_gradient_CHIT1_linePlot_martina.pdf",width = 12,height = 6)

####
# is there any correlation between the CHIT1 avg expressino in the different areas and the cell type deconvulution?
# load the gradient file from the same sampels
test_cellType <- read_tsv("out/table/gradient_cellType_summary_martina.tsv")

df_test_cellType <- test_cellType %>%
  left_join(LUT_sample,by = "slide") %>%
  filter(!is.na(manual_segmentation)) %>% 
  # order the y values
  mutate(area = factor(manual_segmentation,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex")))%>%
  mutate(tissue_class = str_extract(condition,pattern = "active|CAL|CI|remyel|non-lesional")) %>%
  mutate(tissue = toupper(tissue_class))

df_test_cellType %>%
  ggplot(aes(x=area,y=avg,group=slide))+geom_line()+facet_grid(cell_type~tissue,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("out/image/test_gradient_celltype_linePlot_martina.pdf",width = 12,height = 12)

# test a correlation approach
df_full <- left_join(df_test,df_test_cellType,by=c("slide","manual_segmentation","condition","annotation","tissue_class","tissue"),suffix = c(".CHIT1",".cellType"))

# plot the correlations
df_full %>%
  filter(CHIT1 == "CHIT1_norm") %>%
  ggplot(aes(x=avg.CHIT1,y=avg.cellType))+geom_point()+facet_wrap(tissue~cell_type,scales = "free",ncol=7)+geom_smooth(method = "lm")+theme_bw()+theme(strip.background = element_blank())
ggsave("out/image/test_correaltion_celltype_CHIT1_scatterPlot_martina.pdf",width = 18,height = 12)

# try a correlation matrix
mat_corr <- df_full %>%
  filter(CHIT1 == "CHIT1_norm") %>%
  group_by(tissue,cell_type) %>%
  summarise(cor = cor(avg.CHIT1,avg.cellType,method = "pearson")) %>%
  pivot_wider(names_from = cell_type,values_from = cor) %>%
  column_to_rownames("tissue")

pdf("out/image/test_correaltion_celltype_CHIT1_hm_martina.pdf",width = 5,height = 3)
Heatmap(mat_corr)
dev.off()

# pull the p values for the correlation test
mat_corr_pvalue <- df_full %>%
  filter(CHIT1 == "CHIT1_norm") %>%
  group_by(tissue,cell_type) %>%
  nest() %>% 
  mutate(test = map(data,function(x){
    cor.test(x$avg.CHIT1,x$avg.cellType) %>% 
      tidy()
  })) %>% 
  unnest(test) %>% 
  ungroup() %>% 
  mutate(padj = p.adjust(p.value,method = "fdr")) %>% 
  arrange(padj)

# save the full table
write_tsv(mat_corr_pvalue,"out/table/mat_corr_pvalue_gradient_CHIT1.tsv")

# try another option for the plotting
# try to show the trend side by side of the deconvolution and the CHIT1 expression in CAL samples
p1 <- df_test_cellType %>% 
  filter(tissue %in% c("CAL"),
         cell_type %in% c("IMM")) %>% 
  ggplot(aes(x=area,y=avg,group=slide))+geom_line()+facet_grid(cell_type~tissue,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))

p2 <- df_test %>%
  filter(tissue %in% c("CAL")) %>%
  filter(CHIT1 %in% c("CHIT1_norm")) %>%
  ggplot(aes(x=area,y=avg,group=slide))+geom_line()+facet_grid(~tissue,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))

p3 <- df_full %>%
  filter(CHIT1 == "CHIT1_norm") %>%
  filter(tissue %in% c("CAL"),
         cell_type %in% c("IMM")) %>% 
  ggplot(aes(x=avg.CHIT1,y=avg.cellType))+geom_point()+
  # facet_wrap(tissue~cell_type,scales = "free",ncol=7)+
  geom_smooth(method = "lm")+theme_bw()+theme(strip.background = element_blank())

layout <-"AACCC
BBCCC"
p1+p2+p3+ plot_layout(design = layout)
ggsave("out/image/plot_correlation_CHIT1_deconvolution_CAL_IMM.pdf", width = 9,height = 6)
