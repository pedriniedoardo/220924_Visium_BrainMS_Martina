# libraries ---------------------------------------------------------------
library(tidyverse)

# senescence --------------------------------------------------------------
signatures <- dir("data/siganture_new/senescence/") %>%
  str_subset(pattern = ".txt")
name_sign <- signatures %>%
  str_remove_all(pattern = "_review.txt|_signature_short_fixed")

# the classical sigantures
pathway_1 <- pmap(list(signatures,name_sign),function(x,y){
  read_tsv(paste0("data/signature_new/senescence/",x)) %>%
    mutate(Pathway = y) %>%
    mutate(Genes = human_gene) %>%
    dplyr::select(Pathway,Genes)
}) %>%
  setNames(name_sign)

# extract the dataframe of genes in the signature
pathway_2 <- read_csv("data/signature_new/senescence/cellAge/cellAge_fix.csv") %>%
  mutate(Pathway = paste0("CellAge_",senescence_effect),
         Genes = gene_name) %>%
  filter(senescence_effect != "Unclear") %>%
  split(f = .$senescence_effect) %>%
  lapply(.,function(x){
    x %>%
      dplyr::select(Pathway,Genes)
  })

# pull the genes
pathways_senescence <- c(pathway_1,pathway_2)

# save the object
saveRDS(object = pathways_senescence,file = "data/signature_new/senescence_pathways.rds")

# # inflammation ------------------------------------------------------------
# signatures2 <- dir("data/signature_new/inflammation/") %>%
#   str_subset(pattern = ".txt")
# name_sign2 <- signatures2 %>%
#   str_remove_all(pattern = "_review.txt|_signature_short_fixed")
# 
# # the classical sigantures
# pathway_12 <- pmap(list(signatures2,name_sign2),function(x,y){
#   read_tsv(paste0("data/signature_new/inflammation/",x)) %>%
#     mutate(Pathway = y) %>%
#     mutate(Genes = human_gene) %>%
#     dplyr::select(Pathway,Genes)
# }) %>%
#   setNames(name_sign2)
# 
# # pull the genes
# pathways_inflammation <- c(pathway_12)
# 
# # save the object
# saveRDS(object = pathways_inflammation,file = "data/signature_new/inflammation_pathways.rds")
# 
