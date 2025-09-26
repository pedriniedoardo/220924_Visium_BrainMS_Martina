# libraries ---------------------------------------------------------------
library(tidyverse)

# read in the data --------------------------------------------------------
# read in the data already scored by senescence
list_brain <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE.rds")

df <- lapply(list_brain,function(x){
  # total unique genes pre slide
  unique_gene <- dim(x)[1]
  
  x@meta.data %>%
    summarise(n_spot = n(),
              tot_nCount = sum(nCount_Spatial)) |> 
    mutate(avg_count_spot = tot_nCount/n_spot) |> 
    mutate(unique_gene = unique_gene)
    
}) |> 
  bind_rows(.id = "slide")

write_tsv(df,file = "out/table/summary_total_spots.tsv")
