# AIM ---------------------------------------------------------------------
# compile a metadata table to share fro the spatial dataset

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(cowplot)

# read in the data --------------------------------------------------------
# read in the slide classification
list_brain3 <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation3.rds")

# original 
LUT <- read_csv("data/sample_classification.csv")

# add martina's updated classification
LUT_update <- read_csv("data/sample_classification_martina.csv")

LUT_sample <- data.frame(slide = c("V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V14"),
                         condition = c("CAL","CAL","remyel","non-lesional","active","CI","CI","non-lesional","CI","CI","CI","non-lesional","non-lesional")) %>%
  mutate(annotation = paste(condition,slide))

# # metadata to replicate the picture 5 b -----------------------------------
# # pull the meta from all the slides
# meta_all <- lapply(list_brain3, function(sl){
#   sl@meta.data %>%
#     rownames_to_column("barcode")
# }) %>%
#   bind_rows()
# 
# test <- meta_all %>%
#   # filter(dataset != "brain_V07") %>%
#   separate("dataset",into = c("tissue","slide"),sep = "_",remove = F) %>%
#   # add metadata original
#   left_join(LUT,by = c("slide" = "sample")) %>%
#   # add metadata updated
#   left_join(LUT_update,by = c("slide" = "sample")) %>%
#   group_by(slide,type,type_update_martina,manual_segmentation3) %>%
#   summarise(n_spot = n()) %>%
#   ungroup() %>%
#   group_by(slide) %>%
#   mutate(n_spot_tot = sum(n_spot))
# 
# # confirm the numbers of the plot
# 
# 
# # metadata replicate the picture 5 d --------------------------------------
# 
# summary_02 <- meta_all %>%
#   filter(!is.na(manual_segmentation3)) %>%
#   select(dataset,barcode,manual_segmentation3,AST,IMM,LYM,NEU,OLIGO,OPC,VAS) %>%
#   pivot_longer(cols = -c(dataset,barcode,manual_segmentation3), names_to = "celltype", values_to = "value") %>%
#   group_by(dataset,manual_segmentation3,celltype) %>%
#   summarise(mean_prop = mean(value),.groups = "drop")
# 
# # write_tsv(summary_02,"out/table/99_deconvolution_summary_02.tsv")
# 
# 
# # metadata replicate picture 6 b ------------------------------------------
# #
# LUT_sample <- data.frame(slide = c("V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12","V14"),
#                          condition = c("CAL","CAL","remyel","non-lesional","active","CI","CI","non-lesional","CI","CI","CI","non-lesional","non-lesional")) %>%
#   mutate(annotation = paste(condition,slide))
# 
# LUT_sample %>%
#   left_join(LUT_update,by = c("slide" = "sample"))
# 
# tot_scores_summary2 <- meta_all %>%
#   select(barcode,dataset,CASELLA_UP_fixed1:Inhibits1,manual_segmentation) %>%
#   pivot_longer(names_to = "signature",values_to = "score",-c(barcode,manual_segmentation,dataset)) %>%
#   group_by(dataset,manual_segmentation,signature) %>% 
#   summarise(avg = mean(score),
#             sd = sd(score),
#             med = median(score)) %>%
#   separate("dataset",into = c("tissue","slide"),sep = "_",remove = F) %>%
#   # add metadata updated
#   left_join(LUT_update,by = c("slide" = "sample")) %>%
#   # add metadata hardcoded
#   left_join(LUT_sample,by = c("slide"))
# 
# test1 <- tot_scores_summary2 %>%
#   filter(!is.na(manual_segmentation)) %>%
#   # filter(dataset != "brain_V07") %>%
#   filter(signature == "XIMERAKIS_ASC_fixed1") %>%
#   filter(!is.na(avg)) %>%
#   mutate(manual_segmentation = factor(manual_segmentation,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex")))
# 
# test1 %>%
#   group_by(condition,slide) %>%
#   summarise(n = n()) %>%
#   pivot_wider(names_from = condition,values_from = n)
# 
# test1 %>%  
#   # mutate(tissue_class = str_extract(tissue,pattern = "active|CAL|CI|remyel|non.lesional")) %>% 
#   ggplot(aes(x=manual_segmentation,y=avg,group=annotation))+geom_line()+
#   # facet_wrap(~tissue_class,ncol = 1,scales = "free_x")+
#   facet_wrap(~condition,ncol = 1)+
#   theme_cowplot()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# 
# test1 %>%
#   select(annotation,avg,manual_segmentation) %>%
#   pivot_wider(names_from = manual_segmentation,values_from = avg)

# share metadata ----------------------------------------------------------
# pull the full metadata
meta_all <- lapply(list_brain3, function(sl){
  sl@meta.data %>%
    rownames_to_column("barcode")
}) %>%
  bind_rows()

# build the medatada to share
test_share_meta <- meta_all %>%
  filter(dataset != "brain_V07") %>%
  select(barcode,dataset,CASELLA_UP_fixed1:Inhibits1,AST,IMM,LYM,NEU,OLIGO,OPC,VAS,manual_segmentation,manual_segmentation3) %>%
  separate("dataset",into = c("tissue","slide"),sep = "_",remove = F) %>%
  # add metadata updated
  left_join(LUT_update,by = c("slide" = "sample")) %>%
  # add metadata hardcoded fro figure 6 b
  left_join(LUT_sample,by = c("slide"))

# save the table with the info that martina suggested to keep
test_share_meta %>%
  select(-c(type_update_martina,condition,annotation)) %>%
  saveRDS("out/object/99_metadata_spatial_share.rds")

# make sure the annotation for manual segmentaiton makes sense
# notice that desipite the condition paremanter is not matching with the type_update_martina, there is no core or edge annotation of spots for the two misslabelled slides for type_update_martina.
test_share_meta %>%
  group_by(slide,type_update_martina,condition,manual_segmentation,manual_segmentation3) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  select(-manual_segmentation3) %>%
  pivot_wider(names_from = manual_segmentation,values_from = n) %>%
  arrange(type_update_martina)

test_share_meta %>% dim()

# confirm I can replicate the plots
# 5 b
test5b <- test_share_meta %>%
  separate("dataset",into = c("tissue","slide"),sep = "_",remove = F) %>%
  # add metadata updated
  group_by(slide,type_update_martina,manual_segmentation3) %>%
  summarise(n_spot = n()) %>%
  ungroup() %>%
  group_by(slide) %>%
  mutate(n_spot_tot = sum(n_spot))

# test5b %>%
#   write_csv("test5b.csv")

# confirm the counts
test5b2 <- test5b %>%
  ungroup() %>%
  group_by(type_update_martina,manual_segmentation3) %>%
  summarise(tot = sum(n_spot))

# test5b2 %>%
#   write_csv("test5b2.csv")

# 5 d
test_share_meta %>%
  filter(!is.na(manual_segmentation3)) %>%
  select(dataset,barcode,manual_segmentation3,AST,IMM,LYM,NEU,OLIGO,OPC,VAS) %>%
  pivot_longer(cols = -c(dataset,barcode,manual_segmentation3), names_to = "celltype", values_to = "value") %>%
  group_by(dataset,manual_segmentation3,celltype) %>%
  summarise(mean_prop = mean(value),.groups = "drop")

# 6 b
test4 <- test_share_meta %>%
  select(barcode,dataset,CASELLA_UP_fixed1:Inhibits1,manual_segmentation,type_update_martina,condition,annotation) %>%
  # select(barcode,dataset,CASELLA_UP_fixed1:Inhibits1,manual_segmentation) %>%
  # pivot_longer(names_to = "signature",values_to = "score",-c(barcode,manual_segmentation,dataset)) %>%
  pivot_longer(names_to = "signature",values_to = "score",-c(barcode,manual_segmentation,dataset,type_update_martina,condition,annotation)) %>%
  group_by(dataset,manual_segmentation,signature,type_update_martina,condition,annotation) %>% 
  summarise(avg = mean(score),
            sd = sd(score),
            med = median(score)) %>%
  separate("dataset",into = c("tissue","slide"),sep = "_",remove = F) %>%
  filter(!is.na(manual_segmentation)) %>%
  # filter(dataset != "brain_V07") %>%
  filter(signature == "XIMERAKIS_ASC_fixed1") %>%
  filter(!is.na(avg)) %>%
  mutate(manual_segmentation = factor(manual_segmentation,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex")))

test4 %>%
  group_by(condition,slide) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = condition,values_from = n)

test4 %>%
  group_by(slide,type_update_martina,condition) %>%
  summarise()

test4 %>%  
  # mutate(tissue_class = str_extract(tissue,pattern = "active|CAL|CI|remyel|non.lesional")) %>% 
  ggplot(aes(x=manual_segmentation,y=avg,group=annotation))+geom_line()+
  # facet_wrap(~tissue_class,ncol = 1,scales = "free_x")+
  facet_wrap(~condition,ncol = 1)+
  theme_cowplot()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))

test4 %>%  
  # mutate(tissue_class = str_extract(tissue,pattern = "active|CAL|CI|remyel|non.lesional")) %>% 
  ggplot(aes(x=manual_segmentation,y=avg,group=annotation))+geom_line()+
  # facet_wrap(~tissue_class,ncol = 1,scales = "free_x")+
  facet_wrap(~type_update_martina,ncol = 1)+
  theme_cowplot()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
