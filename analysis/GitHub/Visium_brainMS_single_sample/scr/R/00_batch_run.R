# run all the scritps
visium_scripts <- dir("scr/R/") %>% 
  str_subset(pattern = "sample_script_V") %>% 
  str_subset(pattern = "01|02|03|04|05|06|07|08|09|10|11|12|13",negate = T)
  

lapply(paste0("scr/R/",visium_scripts), function(x){
  source(x)
  print(x)
  gc()
  gc()
})


visium_scripts <- dir("scr/R/") %>% 
  str_subset(pattern = "test_spotlight") %>% 
  str_subset(pattern = "01|02|08|13",negate = T)


lapply(paste0("scr/R/",visium_scripts), function(x){
  source(x)
  print(x)
  gc()
  gc()
})
