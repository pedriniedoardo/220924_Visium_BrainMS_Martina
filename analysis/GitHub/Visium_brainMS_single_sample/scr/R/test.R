file <- dir("scr/R/") %>% 
  str_subset(pattern = "test_spotlight") %>% 
  str_subset(pattern = "14|15|16")
  
lapply(file, function(x){
  print(x)
  
  file_dir <- paste0("scr/R/",x)
  print(file_dir)
  
  source(file_dir)
  
  remove(list = ls())
})

