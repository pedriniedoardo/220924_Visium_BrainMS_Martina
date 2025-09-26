# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(patchwork)

# read in data ------------------------------------------------------------
# define all the files available
folder <- "../out_large/"
file_id <- dir(folder) %>% 
  str_subset("^brain_V") %>% 
  str_remove_all(pattern = ".rds")
#
file_id_index <- lapply(str_split(file_id,"_"),function(x){x[2]}) %>% unlist()

# build a full matrix for the expression information on TSPO
list_brain <- lapply(file_id,function(x){
  # keep track of the progress of the reading
  print(x)
  # read in the data
  brain <- readRDS(paste0("../out_large/",x,".rds"))
  # # extract the expression of the feature
  # df_exp <- brain@assays$SCT@data["TSPO",] %>% 
  #   data.frame(exp = .) %>% 
  #   rownames_to_column("barcodes") %>% 
  #   mutate(gene = "TSPO") %>% 
  #   mutate(dataset = x)
  # return(df_exp)
  return(brain)
}) %>% 
  setNames(file_id_index)

# save a list of all the preprocessed brain slices
saveRDS(list_brain,"out/object/list_brain_all.rds")
