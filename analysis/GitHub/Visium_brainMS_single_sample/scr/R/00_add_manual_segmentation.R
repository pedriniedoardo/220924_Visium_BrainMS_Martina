# AIM ---------------------------------------------------------------------
# add the manual segmentation information from Martina to the spatial data


# libraries ---------------------------------------------------------------



# read in the data --------------------------------------------------------
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
  
  # read in the novel segmentation run by Martina
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
  full_meta <- left_join(current_meta,test_meta,"Barcode") %>%
    select(c("Barcode","orig.ident","nCount_Spatial","nFeature_Spatial","nCount_SCT","nFeature_SCT","seurat_clusters","AST","IMM","LYM","NEU","OLIGO","OPC","VAS","res_ss","dataset","class_50","class_70","class_top","CASELLA_UP_fixed1","CLASSICAL_SASP1","EPC_TEAM_SENESCENCE_DOWN_NO_IF1","EPC_TEAM_SENESCENCE_UP_NO_IF1","FRIDMAN_SENESCENCE_UP1","HERNANDEZ_UP_fixed1","PURCELL_fixed1","REACTOME_SENESCENCE_SASP1","SAEPHIA_CURATED_SASP1","SENCAN_fixed1","senmayo1","XIMERAKIS_ASC_fixed1","XIMERAKIS_EC_fixed1","XIMERAKIS_MAC_fixed1","XIMERAKIS_MG_fixed1","XIMERAKIS_mNEUR_fixed1","XIMERAKIS_OLG_fixed1","XIMERAKIS_OPC_fixed1","XIMERAKIS_PC_fixed1","Induces1","Inhibits1","manual_segmentation")) %>%
    column_to_rownames("Barcode")
  
  dim(full_meta)
  
  # swap the metadata
  test@meta.data <- full_meta

  # return the object with full annotation
  return(test)
}) %>% 
  setNames(int_slides)

# check the colnames
lapply(list_out,function(x){
  dim(x@meta.data)
})

# save the annotated objects
saveRDS(list_out,"out/object/list_brain_all_spotlight_SENESCENCE_manualSegmentation.rds")
