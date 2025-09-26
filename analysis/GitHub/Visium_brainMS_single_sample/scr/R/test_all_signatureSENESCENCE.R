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

# read in the data --------------------------------------------------------
list_spotlight <- readRDS("out/object/list_brain_all_spotlight.rds")
# skip V05
# read in the sigantures
list_sig <- readRDS("data/signature_new/senescence_pathways.rds")

# select a specific signature
# list_sig$senmayo
# 
# dataset<-list_spotlight$V08
# name_dataset<-"V08"
# signature<-"senmayo"

# score the signatures and save the meta
list_brain <- pmap(list(list_spotlight,names(list_spotlight)), function(dataset,name_dataset){
  list_df_meta <- lapply(names(list_sig),function(signature){
    # extract the dataframe of genes in the signature
    signature.genes.df <- list_sig[[signature]]
    
    # pull the genes
    signature.genes <- signature.genes.df %>%
      pull(Genes) %>%
      unique()
    
    if(name_dataset %in% c("V10","V11")){
      # score the module
      data.combined <- AddModuleScore(dataset,
                                      features = list(signature.genes),
                                      name=signature,nbin = 15)
    }else if(name_dataset %in% c("V12")){
      # score the module
      data.combined <- AddModuleScore(dataset,
                                      features = list(signature.genes),
                                      name=signature,nbin = 19)
    }else if(name_dataset %in% c("V14")){
      # score the module
      data.combined <- AddModuleScore(dataset,
                                      features = list(signature.genes),
                                      name=signature,nbin = 10)
    }else if(name_dataset %in% c("V15")){
      # score the module
      data.combined <- AddModuleScore(dataset,
                                      features = list(signature.genes),
                                      name=signature,nbin = 9)
    }else if(name_dataset %in% c("V05")){
      # score the module
      data.combined <- AddModuleScore(dataset,
                                      features = list(signature.genes),
                                      name=signature,nbin = 22)
    }else if(name_dataset %in% c("V08")){
      # score the module
      data.combined <- AddModuleScore(dataset,
                                      features = list(signature.genes),
                                      name=signature,nbin = 7)
    }else{
      # score the module
      data.combined <- AddModuleScore(dataset,
                                      features = list(signature.genes),
                                      name=signature)
    }
    # confirm the addition of the score for the module
    # data.combined@meta.data
    
    df_meta <- data.combined@meta.data %>%
      rownames_to_column("barcode")
    
    return(df_meta)
  })
  # put all the signature together
  df_meta_full <- purrr::reduce(list_df_meta,left_join,by=c("barcode","orig.ident","nCount_Spatial","nFeature_Spatial","nCount_SCT","nFeature_SCT","SCT_snn_res.0.8","seurat_clusters","AST","IMM","LYM","NEU","OLIGO","OPC","VAS","res_ss","dataset","class_50","class_70","class_top"))
  
  # save the table with the scores
  df_meta_full %>% 
    write_tsv(paste0("out/table/modules_SENESCENCE/Module_score_SENESCENCE_full_",name_dataset,".tsv"))
  
  # add the meta to the original object
  dataset@meta.data <- df_meta_full %>% column_to_rownames("barcode")
  return(dataset)
})

saveRDS(list_brain,"out/object/list_brain_all_spotlight_SENESCENCE.rds")
list_brain <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE.rds")

# plotting ----------------------------------------------------------------
pdf("out/image/panel_SENESCENCE_all.pdf",width = 40,height = 40)
pmap(list(list_brain,names(list_brain)),function(x,name){
  # save the plot
  plot <- Seurat::SpatialFeaturePlot(
    object = x,
    features = c("CASELLA_UP_fixed1","CLASSICAL_SASP1","EPC_TEAM_SENESCENCE_DOWN_NO_IF1","EPC_TEAM_SENESCENCE_UP_NO_IF1","FRIDMAN_SENESCENCE_UP1","HERNANDEZ_UP_fixed1","PURCELL_fixed1","REACTOME_SENESCENCE_SASP1","SAEPHIA_CURATED_SASP1","SENCAN_fixed1","senmayo1","XIMERAKIS_ASC_fixed1","XIMERAKIS_EC_fixed1","XIMERAKIS_MAC_fixed1","XIMERAKIS_MG_fixed1","XIMERAKIS_mNEUR_fixed1","XIMERAKIS_OLG_fixed1","XIMERAKIS_OPC_fixed1","XIMERAKIS_PC_fixed1","Induces1","Inhibits1"),
    alpha = c(0.1, 1)) + plot_annotation(name)
  
  # print(plot)
})
dev.off()

pdf("out/image/panel_SENESCENCE_all.png",width = 40,height = 40)
pmap(list(list_brain,names(list_brain)),function(x,name){
  # save the plot
  plot <- Seurat::SpatialFeaturePlot(
    object = x,
    features = c("CASELLA_UP_fixed1","CLASSICAL_SASP1","EPC_TEAM_SENESCENCE_DOWN_NO_IF1","EPC_TEAM_SENESCENCE_UP_NO_IF1","FRIDMAN_SENESCENCE_UP1","HERNANDEZ_UP_fixed1","PURCELL_fixed1","REACTOME_SENESCENCE_SASP1","SAEPHIA_CURATED_SASP1","SENCAN_fixed1","senmayo1","XIMERAKIS_ASC_fixed1","XIMERAKIS_EC_fixed1","XIMERAKIS_MAC_fixed1","XIMERAKIS_MG_fixed1","XIMERAKIS_mNEUR_fixed1","XIMERAKIS_OLG_fixed1","XIMERAKIS_OPC_fixed1","XIMERAKIS_PC_fixed1","Induces1","Inhibits1"),
    alpha = c(0.1, 1)) + plot_annotation(name)
  
  # print(plot)
})
dev.off()

Seurat::SpatialFeaturePlot(
  object = list_brain$V01,
  features = c("CASELLA_UP_fixed1","CLASSICAL_SASP1","EPC_TEAM_SENESCENCE_DOWN_NO_IF1","EPC_TEAM_SENESCENCE_UP_NO_IF1","FRIDMAN_SENESCENCE_UP1","HERNANDEZ_UP_fixed1","PURCELL_fixed1","REACTOME_SENESCENCE_SASP1","SAEPHIA_CURATED_SASP1","SENCAN_fixed1","senmayo1","XIMERAKIS_ASC_fixed1","XIMERAKIS_EC_fixed1","XIMERAKIS_MAC_fixed1","XIMERAKIS_MG_fixed1","XIMERAKIS_mNEUR_fixed1","XIMERAKIS_OLG_fixed1","XIMERAKIS_OPC_fixed1","XIMERAKIS_PC_fixed1","Induces1","Inhibits1"),
  alpha = c(0.1, 1))+plot_annotation("test")
ggsave("out/image/panel_SENESCENCE_V01.pdf",width = 40,height = 40)

Seurat::SpatialFeaturePlot(
  object = list_brain$V02,
  features = c("CASELLA_UP_fixed1","CLASSICAL_SASP1","EPC_TEAM_SENESCENCE_DOWN_NO_IF1","EPC_TEAM_SENESCENCE_UP_NO_IF1","FRIDMAN_SENESCENCE_UP1","HERNANDEZ_UP_fixed1","PURCELL_fixed1","REACTOME_SENESCENCE_SASP1","SAEPHIA_CURATED_SASP1","SENCAN_fixed1","senmayo1","XIMERAKIS_ASC_fixed1","XIMERAKIS_EC_fixed1","XIMERAKIS_MAC_fixed1","XIMERAKIS_MG_fixed1","XIMERAKIS_mNEUR_fixed1","XIMERAKIS_OLG_fixed1","XIMERAKIS_OPC_fixed1","XIMERAKIS_PC_fixed1","Induces1","Inhibits1"),
  alpha = c(0.1, 1))
ggsave("out/image/panel_SENESCENCE_V02.pdf",width = 40,height = 40)


Seurat::SpatialFeaturePlot(
  object = list_brain$V05,
  features = c("CASELLA_UP_fixed1","CLASSICAL_SASP1","EPC_TEAM_SENESCENCE_DOWN_NO_IF1","EPC_TEAM_SENESCENCE_UP_NO_IF1","FRIDMAN_SENESCENCE_UP1","HERNANDEZ_UP_fixed1","PURCELL_fixed1","REACTOME_SENESCENCE_SASP1","SAEPHIA_CURATED_SASP1","SENCAN_fixed1","senmayo1","XIMERAKIS_ASC_fixed1","XIMERAKIS_EC_fixed1","XIMERAKIS_MAC_fixed1","XIMERAKIS_MG_fixed1","XIMERAKIS_mNEUR_fixed1","XIMERAKIS_OLG_fixed1","XIMERAKIS_OPC_fixed1","XIMERAKIS_PC_fixed1","Induces1","Inhibits1"),
  alpha = c(0.1, 1))
ggsave("out/image/panel_SENESCENCE_V05.pdf",width = 40,height = 40)

# what is the distibutino of the score based on the proportion of different cell types
list_brain$V01@meta.data %>%
  # ordet the cluster by the median value of the score
  group_by(class_top) %>% 
  mutate(score_median_senmayo = median(senmayo1)) %>%
  ungroup() %>% 
  mutate(class_top = fct_reorder(class_top,score_median_senmayo)) %>%
  # mutate(BraakStage=as.factor(BraakStage)) %>% 
  # ggplot(aes(x=signature_score1,y=fct_reorder(seurat_clusters_fix,score_median)))+
  ggplot(aes(x=senmayo1,y=class_top,fill=after_stat(x)))+
  # geom_boxplot()
  # geom_density_ridges(alpha=0.5) +
  geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2,alpha=0.5,vline_linetype="dashed")+
  scale_fill_viridis_c(option = "turbo")+
  # facet_wrap(~pathology)+
  theme_ridges()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))

# check the correlation between senmayo and propo of cell type
test <- list_brain$V01@meta.data %>%
  rownames_to_column("barcode") %>% 
  dplyr::select(barcode,AST:VAS,senmayo1) %>% 
  pivot_longer(names_to = "class",values_to = "prop_class",AST:VAS)

# classical scatter plot
test %>% 
  ggplot(aes(x=prop_class,y = senmayo1))+geom_point(alpha=0.1)+facet_wrap(~class)+theme_cowplot()+theme(strip.background = element_blank())+geom_smooth(method = "lm")

test %>%
  # filter(class=="VAS") %>%
  filter(prop_class>0) %>% 
  ggplot(aes(x=prop_class,y=senmayo1)) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  )+scale_fill_viridis_c(option = "turbo")+facet_wrap(~class)+theme_cowplot()+theme(strip.background = element_blank())

# test fix the axes
list_plot_density <- lapply(c("AST","IMM","LYM","NEU","OLIGO","OPC","VAS"),function(x){
  test %>%
    filter(class==x) %>%
    filter(prop_class>0) %>% 
    ggplot(aes(x=prop_class,y=senmayo1)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_x_continuous(limits = c(0,1),expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.06,0.06),expand = c(0,0)) +
    theme(
      legend.position='none'
    )+scale_fill_viridis_c(option = "turbo")+facet_wrap(~class)+theme_cowplot()+theme(strip.background = element_blank())
})

patchwork::wrap_plots(list_plot_density)

# try with the side sitribution
test %>%
  filter(class=="VAS") %>%
  filter(prop_class>0) %>% 
  ggplot(aes(x=prop_class,y=senmayo1)) +
  geom_point(alpha=0.1)+facet_wrap(~class)+theme_cowplot()+theme(strip.background = element_blank())+geom_smooth(method = "lm")+
  geom_xsidedensity(aes(y = after_stat(density)), position = "stack",xfill="gray") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack",yfill="gray")

test %>%
  # filter(class=="VAS") %>%
  filter(prop_class>0) %>% 
  # ggplot(aes(x=prop_class,y=senmayo1,col=class)) +
  ggplot(aes(x=prop_class,y=senmayo1)) +
  geom_point(alpha=0.1)+theme_cowplot()+theme(strip.background = element_blank())+geom_smooth(method = "lm")+
  geom_xsidedensity(aes(y = after_stat(density),xfill=class),alpha=0.5) +
  geom_ysidedensity(aes(x = after_stat(density),yfill=class),alpha=0.5)

test %>%
  # filter(class=="VAS") %>%
  filter(prop_class>0) %>% 
  # ordet the cluster by the median value of the score
  group_by(class) %>% 
  mutate(score_median_prop = median(prop_class)) %>%
  ungroup() %>% 
  mutate(class = fct_reorder(class,score_median_prop)) %>% # pull(class) %>% levels()
  # ggplot(aes(x=prop_class,y=senmayo1,col=class)) +
  ggplot(aes(x=prop_class,y=senmayo1)) +
  geom_point(alpha=0.1)+theme_cowplot()+theme(strip.background = element_blank())+
  geom_xsideboxplot(aes(y = class,fill=class), orientation = "y") +
  scale_xsidey_discrete()

test %>%
  # filter(class=="VAS") %>%
  filter(prop_class>0) %>% 
  # ordet the cluster by the median value of the score
  group_by(class) %>% 
  mutate(score_median_senmayo = median(senmayo1)) %>%
  ungroup() %>% 
  mutate(class = fct_reorder(class,score_median_senmayo)) %>% # pull(class) %>% levels()
  # ggplot(aes(x=prop_class,y=senmayo1,col=class)) +
  ggplot(aes(y=prop_class,x=senmayo1)) +
  geom_point(alpha=0.1,aes(col=class))+theme_cowplot()+theme(strip.background = element_blank())+
  geom_xsideboxplot(aes(y = class,fill=class), orientation = "y",outlier.alpha = 0.1) +
  scale_xsidey_discrete()+
  geom_ysidedensity(aes(col=class))+
  theme(ggside.panel.scale = .2)

# SpatialFeaturePlot(dataset, features = "nCount_Spatial") + theme(legend.position = "right")                
