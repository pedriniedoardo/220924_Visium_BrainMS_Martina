##########################################################
############## SENESCENCE INDEX TOOL (SIT) ###############
##########################################################

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(homologene)
library(RColorBrewer)
library(scales)
library(msigdbr)
library(clusterProfiler)
library(GSVA)
library(limma)
library(cowplot)

# library(org.Mm.eg.db)

# read in the seignatures -------------------------------------------------
# get the specific pathways
# get a shortlist of the terms
m_df <- msigdbr(species = "Homo sapiens",category = "C2")

# filter only the terms fo interest and pull the genes
GSred <- m_df %>% 
  filter(gs_name %in% c("KEGG_CELL_CYCLE","REACTOME_CELL_CYCLE" ,"WP_CELL_CYCLE")) %>% 
  split(f = .$gs_name) %>% 
  map(function(x){
    x %>% 
      pull(gene_symbol)
  })

# define the specific signature score -------------------------------------
# reference markers from the authors
senescence_marker_mouse <- c("Cdkn2a","Cdkn1a",  "Serpine1" ,"Cdkn1b", "Cdkn2d", "Cdkn2b")

# convert the human gene neames
df2 <- homologene(senescence_marker_mouse, inTax = 10090, outTax = 9606) %>% 
  setNames(c("mouse_gene","human_gene","mouse_ID","human_ID"))
df2 
senescence_marker <- df2$human_gene

# read in the data --------------------------------------------------------
# read in the data already scored by senescence
list_brain <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE.rds")
# Seurat.object_all <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15.rds")
# Seurat.object <- list_brain$V01
# make it slimmer to speed up the computation
# Seurat.object <- subset(Seurat.object_all,subset = ID == "RR16_BASELINE_0h")

# DimPlot(Seurat.object)
# Seurat::SpatialDimPlot(object = Seurat.object,group.by = "seurat_clusters")
# VlnPlot(Seurat.object, senescence_marker)

# generate a lut for both data and segmentations
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

LUT_sample <- tot_metadata %>% 
  group_by(slide_fix,file) %>% 
  summarise() %>% 
  # use only segmentation for which there is a slide file
  filter(slide_fix %in% names(list_brain)) %>% 
  # skip slide 8 remove also slide 10
  filter(!slide_fix %in% c("V08","V10","V11"))

# run the scoring ---------------------------------------------------------
# setdiff(names(list_brain),LUT_sample$slide_fix)
# setdiff(LUT_sample$slide_fix,names(list_brain))

# slide <- "V01"
# segmentation <- "sofia_segmentation_01.csv"

list_lut <- pmap(list(LUT_sample$slide_fix,LUT_sample$file),function(slide,segmentation){
  # keep track fo the pregress
  print(slide)
  
  # read in the object
  Seurat.object <- list_brain[[slide]]
  
  # add the score fo the subset of markers
  if(slide %in% c("V10","V11")){
    # score the module
    # Seurat.object <- AddModuleScore(Seurat.object, features = list(senescence_marker))
    Seurat.object <- AddModuleScore(Seurat.object,
                                    features = list(senescence_marker),nbin = 15)
  }else if(slide %in% c("V12")){
    # score the module
    Seurat.object <- AddModuleScore(Seurat.object,
                                    features = list(senescence_marker),nbin = 19)
  }else if(slide %in% c("V14")){
    # score the module
    Seurat.object <- AddModuleScore(Seurat.object,
                                    features = list(senescence_marker),nbin = 10)
  }else if(slide %in% c("V15")){
    # score the module
    Seurat.object <- AddModuleScore(Seurat.object,
                                    features = list(senescence_marker),nbin = 9)
  }else if(slide %in% c("V05")){
    # score the module
    Seurat.object <- AddModuleScore(Seurat.object,
                                    features = list(senescence_marker),nbin = 22)
  }else if(slide %in% c("V08")){
    # score the module
    Seurat.object <- AddModuleScore(Seurat.object,
                                    features = list(senescence_marker),nbin = 7)
  }else{
    # score the module
    Seurat.object <- AddModuleScore(Seurat.object,
                                    features = list(senescence_marker))
  }
  # Seurat.object@meta.data
  # change the name of the AddModuleScore to Senescence_signature
  Seurat.object$Cluster1-> Seurat.object$Senescence_signature
  # Seurat.object@meta.data
  
  # -------------------------------------------------------------------------
  # this is too mem intensive
  # exprdf <- as.matrix(Seurat.object@assays$RNA@data)
  # gsvascore1 <- GSVA::gsva(expr = exprdf, GSred[1],  mx.diff=T, method="ssgsea", parallel.sz=24)
  
  # test use less cores to reduce RAM memory
  gsvascore <- GSVA::gsva(expr = Seurat.object@assays$Spatial@data, GSred,  mx.diff=T, method="ssgsea", parallel.sz=8)
  # free the memory
  # gc()
  gsvascore <- normalizeQuantiles(gsvascore)
  # save the object
  # saveRDS(gsvascore,"out/object/gsvascore_V01.rds")
  # gsvascore <- readRDS("../../out/object/gsvascore.rds")
  
  # Kegg signature, add the scores to the object
  Seurat.object$KEGG_CELL_CYCLE <- gsvascore["KEGG_CELL_CYCLE", colnames(Seurat.object)]
  Seurat.object$Cell_cycle_arrest_signatureK <- (-(Seurat.object$KEGG_CELL_CYCLE))
  
  # Reactome signature add the scores to the object
  Seurat.object$REACTOME_CELL_CYCLE <- gsvascore["REACTOME_CELL_CYCLE", colnames(Seurat.object)]
  Seurat.object$Cell_cycle_arrest_signatureR <- -(Seurat.object$REACTOME_CELL_CYCLE)
  
  # WP signature add the scores to the object
  Seurat.object$WP_CELL_CYCLE <- gsvascore["WP_CELL_CYCLE", colnames(Seurat.object)]
  Seurat.object$Cell_cycle_arrest_signatureWP <- -(Seurat.object$WP_CELL_CYCLE)
  
  ########  SENESCENCE SCORE ########
  senenorm <- scales::rescale(x=Seurat.object$Senescence_signature, c(0,1))
  keggnorm <- scales::rescale(x=Seurat.object$Cell_cycle_arrest_signatureK, c(0,1))
  reactomenorm <- scales::rescale(x=Seurat.object$Cell_cycle_arrest_signatureR, c(0,1))
  wpnorm <- scales::rescale(x=Seurat.object$Cell_cycle_arrest_signatureWP, c(0,1))
  
  Seurat.object$Senescence_scoreK <- senenorm+keggnorm
  Seurat.object$Senescence_scoreWP <- senenorm+wpnorm
  Seurat.object$Senescence_scoreR <- senenorm+reactomenorm
  
  ##########  Identification of senescent cells
  quantile <- cbind(as.character(cut(as.numeric(as.character(Seurat.object$Senescence_scoreK)),as.numeric(quantile(Seurat.object$Senescence_scoreK)))),
                    as.character(cut(as.numeric(as.character(Seurat.object$Senescence_scoreR)),as.numeric(quantile(Seurat.object$Senescence_scoreR)))),
                    as.character(cut(as.numeric(as.character(Seurat.object$Senescence_scoreWP)),as.numeric(quantile(Seurat.object$Senescence_scoreWP)))))
  
  # rename the levels of the quantiles using the level nomenclature
  df_quantiles <- lapply(1:ncol(quantile),function(x){
    data.frame(x = factor(quantile[,x],labels = c("quantile1","quantile2","quantile3","quantile4")))
  }) %>% 
    bind_cols() %>% 
    # rename the columns
    setNames(c("Senescence_scoreKQUANTILE", "Senescence_scoreRQUANTILE", "Senescence_scoreWPQUANTILE")) %>% 
    # add the quantiles
    mutate(barcode = colnames(Seurat.object))
  
  # ref metadata
  df_meta_ref <- Seurat.object@meta.data
  
  # merge the metadata
  df_meta_full <- df_meta_ref %>% 
    rownames_to_column("barcode") %>% 
    # add the quantile information
    left_join(df_quantiles,by = "barcode") %>% 
    # add the criterion of selection for the senescent cells. all the cells should be in the top quantiles for all the signatures
    mutate(SENEQUANTILE = case_when(Senescence_scoreKQUANTILE == "quantile4" &
                                      Senescence_scoreRQUANTILE == "quantile4" &
                                      Senescence_scoreWPQUANTILE == "quantile4" ~ "YES",
                                    T ~ "NO"))
  
  # # check
  # df_meta_full %>% 
  #   ggplot(aes(x=1,y=Senescence_scoreK))+
  #   geom_boxplot()+
  #   geom_point(aes(col=Senescence_scoreKQUANTILE),position=position_jitter(width = 0.2),alpha=0.5)
  
  # update the metadata in the seurat object
  Seurat.object@meta.data <- df_meta_full %>% 
    column_to_rownames("barcode")
  
  # extract the metadata to be more flexible with the plotting
  # save the current meta add also the coordinates of the UMAP
  df_umap <- Seurat.object@reductions$umap@cell.embeddings %>% 
    data.frame() %>% 
    rownames_to_column()
  
  df_meta <- Seurat.object@meta.data %>% 
    data.frame() %>% 
    rownames_to_column()
  
  df_meta_full <- left_join(df_umap,df_meta,"rowname")
  
  folder <- "data/segmentation/martina_new/"
  # count the senescent cells per sample
  df_meta_SIT <- df_meta_full %>% 
    left_join(read_csv(paste0(folder,segmentation,".csv")),by = c("rowname"="Barcode"))
  
  write_csv(df_meta_SIT,file = paste0("out/table/df_meta_SIT_",slide,".csv"))
  
  return(df_meta_SIT)
}) %>% 
  setNames(LUT_sample$slide_fix)
# save the list object
saveRDS(list_lut,"out/object/list_SIT.rds")

# -------------------------------------------------------------------------
LUT_sample2 <- data.frame(slide = c("V01","V02","V03","V04","V05","V06","V08","V09","V10","V11","V12","V14"),
                         condition = c("CAL","CAL","remyel","non-lesional","active","CI","non-lesional","CI","CI","CI","non-lesional","non-lesional")) %>% 
  mutate(dataset = paste0("brain_",slide))

df_summary <- list_lut %>%
  lapply(function(x){
    x %>%
      dplyr::select(dataset,manual_segmentation,SENEQUANTILE) 
  }) %>% 
  bind_rows() %>% 
  group_by(dataset,manual_segmentation,SENEQUANTILE) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(dataset,manual_segmentation) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop_senescence = n/tot) %>% 
  left_join(LUT_sample2,by = "dataset")

df_summary %>% 
  write_tsv("out/table/gradient_SIT_summary_martina.tsv")


df_summary %>%
  filter(SENEQUANTILE=="YES") %>% 
  mutate(condition = factor(condition,levels = c("non-lesional","CI","active","CAL","remyel"))) %>% 
  ggplot(aes(x=condition,y=prop_senescence))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/SIT_tissue.pdf",width = 4,height = 3)

df_summary %>%
  filter(SENEQUANTILE=="YES",!is.na(manual_segmentation)) %>% 
  mutate(condition = factor(condition,levels = c("non-lesional","CI","active","CAL","remyel"))) %>%
  mutate(manual_segmentation = factor(manual_segmentation,levels = c("core","edge","periplaque 1","periplaque 2","periplaque 3","periplaque 4","periplaque 5","periplaque 6","cortex"))) %>% 
  ggplot(aes(x=manual_segmentation,y=prop_senescence))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  theme_bw() +
  facet_wrap(~condition,scales = "free")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("out/image/SIT_tissue_split.pdf",width = 12,height = 6)

# full map
df_meta_full %>%
  arrange(SENEQUANTILE) %>% 
  ggplot(aes(x=UMAP_1,y=UMAP_2,col=SENEQUANTILE))+geom_point(size = 0.1,alpha=0.5)+theme_cowplot()+
  scale_color_manual(values = c("NO"="#F15A2B" , "YES"="#1D75BC"))+
  guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))
# ggsave("../../out/image/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_SIT.pdf",height = 10,width = 11)

# try to plot the spatial features
(Seurat::SpatialDimPlot(object = Seurat.object,group.by = "SENEQUANTILE",alpha = 0) + theme(legend.position = "none"))+
  (Seurat::SpatialDimPlot(object = Seurat.object,group.by = "SENEQUANTILE",alpha = 0) + theme(legend.position = "none"))+
(Seurat::SpatialDimPlot(object = Seurat.object,group.by = "SENEQUANTILE",alpha = 0.7)+scale_fill_manual(values = c("YES"="#F15A2B" , "NO"="#1D75BC")))+
  (Seurat::SpatialDimPlot(object = Seurat.object,group.by = "SENEQUANTILE",alpha = 0.7)+scale_fill_manual(values = c("NO"="#F15A2B" , "YES"="#1D75BC")))
# read the new annotation
list_lut <- readRDS("out/object/list_SIT.rds")

pdf("out/image/00_list_slide_SIT_gradient.pdf",width = 15,height = 5)
# run the analysis only on the common samples
lapply(LUT_sample$slide_fix,function(slide){
  # read in the object
  brain <- list_brain[[slide]]
  # read in the new annotation
  meta_brain <- list_lut[[slide]] %>% 
    column_to_rownames("rowname")
  
  # update the metadata
  brain@meta.data <- meta_brain
  p1 <- Seurat::SpatialDimPlot(object = brain,group.by = "SENEQUANTILE",alpha = 0)+ theme(legend.position = "none") + ggtitle(slide)
  p2 <- Seurat::SpatialDimPlot(object = brain,group.by = "manual_segmentation",alpha = 1)
  # p1 <- SpatialFeaturePlot(brain, features = "nCount_Spatial",alpha=0,crop = F) + theme(legend.position = "none") + ggtitle(slide)
  p3 <- Seurat::SpatialDimPlot(object = brain,group.by = "SENEQUANTILE",alpha = 0.7)+scale_fill_manual(values = c("NO"="#F15A2B" , "YES"="#1D75BC"))
  p1+p2+p3
})
dev.off()


