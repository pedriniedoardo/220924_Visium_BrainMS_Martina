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
library(factoextra)
library(magick)
library(enrichR)

# read in the data --------------------------------------------------------
# read in the data already scored by senescence
list_brain <- readRDS("out/object/list_brain_all_spotlight_SENESCENCE.rds")

# aggragate the expression per sample
list_exp <- pmap(list(list_brain,names(list_brain)),function(x,nam){
    test <- AggregateExpression(x,group.by = "orig.ident")
    test$SCT %>%
    data.frame() %>%
    mutate(slide=nam) %>%
    rownames_to_column("gene")
})

# make it a table of expression
df_exp <- list_exp %>%
bind_rows() %>%
# fill with 0 if NA
pivot_wider(names_from = slide,values_from = all,values_fill = 0)

df_summary_gene <- df_exp %>%
pivot_longer(names_to = "sample",values_to = "exp",-gene) %>%
group_by(gene) %>%
summarise(tot = sum(exp),n_zero = sum(exp==0))

df_summary_gene %>%
filter(n_zero==0)

df_summary_gene %>%
filter(n_zero!=0)

df_summary_gene %>%
filter(gene=="TSPO")

# keep only the genes that are expressed in all the sample, n_zero = 0
keep_genes <- df_summary_gene %>%
filter(n_zero == 0) %>%
pull(gene)

df_filter <- df_exp %>%
pivot_longer(names_to = "sample",values_to = "exp",-gene) %>%
filter(gene %in% keep_genes) %>%
# add the summary of the total number of reads per sample
group_by(sample) %>%
mutate(tot = sum(exp)) %>%
mutate(norm_count = exp/tot*10^6)

# scale the data by row
df_filter_scale <- df_filter %>%
dplyr::select(-c(exp,tot)) %>%
pivot_wider(names_from = sample,,values_from = norm_count) %>%
column_to_rownames("gene") %>%
t(.) %>%
scale(center=T, scale=T) %>%
t(.)

# save the total normalized
saveRDS(df_filter_scale,"./out/object/df_filter_scale.rds")

# confirm the center and sd of the features
df_filter_scale %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(names_to = "sample",values_to = "zscore",-gene) %>% 
  group_by(gene) %>% 
  summarise(avg = mean(zscore)) %>% 
  ggplot(aes(x=avg))+geom_histogram()+theme_bw()
df_filter_scale %>% 
  data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(names_to = "sample",values_to = "zscore",-gene) %>% 
  group_by(gene) %>% 
  summarise(sd = sd(zscore)) %>% 
  ggplot(aes(x=sd))+geom_histogram()+theme_bw()

# used the filtered scaled dataset
test <- df_filter_scale %>%
  t() %>%
  cor()

test[1:5,1:5]

# save the total correlation matrix
saveRDS(test,"./out/object/test_corr_matrix.rds")

# Perform hierarchical clustering on the correlation matrix
dist_matrix <- as.dist(1 - test)  # Convert correlation matrix to distance matrix
glimpse(dist_matrix)

hc_result <- hclust(dist_matrix)  # Hierarchical clustering
hc_result

# -------------------------------------------------------------------------
# Determine the number of modules (clusters)
# reference: http://www.sthda.com/english/articles/29-cluster-validation-essentials/96-determiningthe-optimal-number-of-clusters-3-must-know-methods/
# Elbow method 
fviz_nbclust(test, kmeans, method = "wss") +
  # geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# set.seed(123)
# fviz_nbclust(test, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+   labs(subtitle = "Gap statistic method")

num_modules <- 7

# Cut the dendrogram to obtain modules
module_labels <- cutree(hc_result, k = num_modules)
str(module_labels)

module_labels_fix <- paste0("m_",module_labels)
str(module_labels_fix)

# Convert module labels to colors
# module_colors <- colorRamp2(seq(min(module_labels), max(module_labels)),
#                             colors = c("red", "green", "blue"))(module_labels)

# build the annotation object 
# column_ha <- HeatmapAnnotation(treat = sample_ordered, 
#                                col = list(treat = c("VEH" = "green", "BMP9" = "gray")))  
column_ha <- HeatmapAnnotation(module = module_labels_fix)
column_ha

pdf("out/image/heatmap_correaltion_matrix.pdf",width = 20,height = 20)
# Create a heatmap with module assignments
Heatmap(t(test),
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 8),
        row_title = "gene",
        column_title = "gene",
        clustering_distance_rows = dist_matrix,
        clustering_distance_columns = dist_matrix,
        # clustering_method_rows = "none",
        # clustering_method_columns = "none",
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE,
        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        # row_annotation = anno_text(labels = module_labels,
        #                            col = module_colors,
        #                            fontface = "bold",
        #                            size = 8),
        # row_annotation_side = "left",
        # row_annotation_width = unit(1, "cm"),
        show_column_names = FALSE,
        show_row_names = FALSE,top_annotation = column_ha)
dev.off()

# determin the mudules identity
# pull the features belonging to the same module
df_modules <- data.frame(feature = names(module_labels),module = module_labels_fix) %>% 
  arrange(module)
head(df_modules)

# save the table
df_modules %>%
write_tsv("./out/table/df_correlation_modules.tsv")

module_id <- df_modules %>%
dplyr::filter(feature=="TSPO") %>%
pull(module)

####
# run enrichr with the list of genes in the module
# DB selection ------------------------------------------------------------
dbs <- listEnrichrDbs()
#filter fo the db of interest
dbs %>%
  filter(str_detect(libraryName,pattern = "Atlas"))

dbs %>%
  filter(str_detect(libraryName,pattern = "Cell"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2016","HDSigDB_Human_2021","Azimuth_Cell_Types_2021","GO_Biological_Process_2023","Descartes_Cell_Types_and_Tissue_2021","CellMarker_Augmented_2021")

# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs
list_genes <- df_modules %>%
  filter(module == module_id) %>% 
  pull(feature) %>% 
  list("m_2" = .)

# define the background
# background <- df_modules$feature

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr <- lapply(list_genes,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>% 
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

list_enrichr %>%
  write_tsv("./out/table/enrichR_m2.tsv")

# see the top term from azimuth cell type
list_enrichr %>%
  as.tibble() %>% 
  # group_by(annotation) %>% summarise()
  filter(annotation == "Azimuth_Cell_Types_2021") %>% 
  # select(Term,Overlap,Genes,Adjusted.P.value,Odds.Ratio) %>% 
  data.frame()  %>% 
  head()

# list  <- read_tsv("out/table/enrichR_out_all_filtered_clusters_CTRL_vs_CSF_harmonyMartina.tsv")

plot_list_UP <- list_enrichr %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")

  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_UP)
ggsave("./out/image/enrichR_m2.pdf",width = 7,height = 25,limitsize = FALSE)

