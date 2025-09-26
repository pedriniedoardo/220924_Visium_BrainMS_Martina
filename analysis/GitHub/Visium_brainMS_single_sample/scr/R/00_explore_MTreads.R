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
library(ggrepel)
library(cowplot)

# setup parallel ----------------------------------------------------------
# library(future)
# plan("multiprocess",workers=12)
# plan()

# read in the data --------------------------------------------------------
# list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
# Load10X_Spatial(data.dir = data_dir)
data_dir <- "../../../raw_data/spaceranger_out/"
sample_id <- dir(data_dir)

# extrapolate the percent MT per sample
list_meta <- lapply(sample_id,function(x){
  brain <- Load10X_Spatial(data.dir = paste0(data_dir,x,"/outs/"))
  brain$percent.mt <- PercentageFeatureSet(brain, pattern = "^MT-")
  brain@meta.data %>%
    mutate(sample_id = x)
})

# average the esprssion per gene in each slide
# x <- "01_SP1_novaseq"
list_cpm <- lapply(sample_id,function(x){
  brain <- Load10X_Spatial(data.dir = paste0(data_dir,x,"/outs/"))
  
  # aggregare all the raw counts
  test <- AggregateExpression(brain,slot = "count")
  
  # make the counts as CPM
  tot <- colSums(test$Spatial)
  
  # make it a dataframe
  df <- test$Spatial %>%
    data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(cpm = all/tot*10^6) %>%
    mutate(sample = x)
  
  return(df)
})

# -------------------------------------------------------------------------
# plot the prop of MT
list_meta %>%
  bind_rows() %>%
  ggplot(aes(x=percent.mt)) + geom_histogram()+facet_wrap(~sample_id)+theme_bw()+theme(strip.background = element_blank())

# plot the gene expression
df_cpm <- list_cpm %>%
  bind_rows() %>%
  group_by(gene) %>%
  summarise(avg_cpm = mean(cpm))

# filter onlyt he MT genes
df_cpm_mt <- 
  df_cpm %>%
  filter(str_detect(gene,pattern = "^MT-"))

df_cpm %>%
  # filter(avg_cpm > 0) %>%
  ggplot(aes(y=avg_cpm,x="brain spatial",label = gene)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),shape = 1,alpha = 0.01) +
  geom_point(data = df_cpm_mt,position = position_jitter(width = 0.2),shape = 1,alpha = 0.01) +
  ggrepel::geom_text_repel(data = df_cpm_mt,aes(y=avg_cpm),force = 30,min.segment.length = 0) +
  theme_bw() +
  scale_y_continuous(trans = "log1p")

df_cpm %>%
  # filter(avg_cpm > 0) %>%
  ggplot(aes(y=avg_cpm,x="brain spatial",label = gene)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),shape = 1,alpha = 0.01) +
  geom_point(data = df_cpm_mt,position = position_jitter(width = 0.2),shape = 1,alpha = 0.01) +
  ggrepel::geom_text_repel(
    data = df_cpm_mt,
    aes(y=avg_cpm),
    force = 30,
    min.segment.length = 0,
    # --- Add these lines for curved connections ---
    segment.curvature = -0.1, # Adjust this value (e.g., -0.1 to 0.5) for desired curve
    segment.angle = 20,       # Adjust this value (e.g., 0 to 180) for angle
    segment.ncp = 3           # Number of control points for bezier curve (default is 5, 3 is good for simpler curves)
    # -----------------------------------------------
  ) +
  theme_bw() +
  scale_y_continuous(trans = "log1p")

df_cpm %>%
  # filter(avg_cpm > 0) %>%
  ggplot(aes(y=avg_cpm,x="brain spatial",label = gene)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_point(position = position_jitter(width = 0.2),shape = 1,alpha = 0.01) +
  geom_point(data = df_cpm_mt,shape = 1,col="red") +
  ggrepel::geom_text_repel(
    data = df_cpm_mt,
    aes(y=avg_cpm),
    force = 20,
    min.segment.length = 0,
    segment.curvature = -0.1,
    segment.angle = 20,
    segment.ncp = 3,segment.alpha = 0.5,
    # # --- Add or adjust these values to increase distance ---
    # box.padding = 0.5,  # Increase this value (default is 0.25)
    point.padding = 0.2 # Increase this value (default is 0.25)
    # # ----------------------------------------------------
  ) +
  theme_bw() +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,50000))  

ggsave("out/image/00_explore_MTReads_01.pdf",width = 4,height = 8)

df_cpm %>%
  # filter(avg_cpm > 0) %>%
  ggplot(aes(y=avg_cpm,x="brain spatial",label = gene)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_point(position = position_jitter(width = 0.2),shape = 1,alpha = 0.01) +
  geom_point(data = df_cpm_mt,shape = 1,col="red") +
  ggrepel::geom_text_repel(
    data = df_cpm_mt,
    aes(y=avg_cpm),
    force = 20,
    min.segment.length = 0,
    segment.alpha = 0.5,
    # # --- Add or adjust these values to increase distance ---
    # box.padding = 0.5,  # Increase this value (default is 0.25)
    point.padding = 0.2 # Increase this value (default is 0.25)
    # # ----------------------------------------------------
  ) +
  theme_bw() +
  scale_y_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,50000))  


df_cpm %>%
  # filter(avg_cpm > 0) %>%
  ggplot(aes(x=avg_cpm,label = gene)) +
  geom_histogram() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_x_continuous(trans = "log1p",breaks = c(0,10,100,1000,10000,50000))  

ggsave("out/image/00_explore_MTReads_01.pdf",width = 4,height = 8)

ggplot(df_cpm, aes(x = avg_cpm)) +
  geom_histogram(aes(y = ..density..), fill = "gray", color = "black") + # Using density for y-axis
  # geom_histogram(aes(y = ..density..)) +
  # You can also use y = ..count.. for raw counts on y-axis
  # geom_histogram() + # Original way using counts
  
  # Add geom_text_repel for the labels
  # Use the filtered 'df_cpm_labeled' for the data argument
  ggrepel::geom_text_repel(
    data = df_cpm_mt, # Use the subset of data for labels
    aes(label = gene, y = 0), # Labels on x-axis, adjusting y to appear "below" points or near axis
    # You might need to adjust y for positioning, e.g., y = 0 for labels near x-axis.
    # If you want them floating, you might need to make 'y' a calculated aesthetic
    # that represents a suitable vertical position for the label.
    # For histograms, there isn't a direct 'y' value for the 'gene' like in a scatter.
    # A common trick is to set 'y' to a constant or use `position = position_nudge_y` if needed.
    # For now, setting y=0 places them at the bottom, and ggrepel adjusts.
    direction = "y",         # Constrain movement to y-direction
    nudge_y = -0.05,         # Nudge labels slightly downwards (adjust as needed)
    force = 20,               # Amount of force to push labels away from each other
    box.padding = 0.5,       # Padding around the label box
    point.padding = 0.5,     # Padding around the point (here, the 'point' is conceptually the x-value)
    min.segment.length = 0.1, # Minimum length of the segment line
    segment.color = "grey50", # Color of the connecting line
    segment.size = 0.2,       # Thickness of the connecting line
    segment.alpha = 0.6,      # Transparency of the connecting line
    size=3
  ) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000, 50000)) +
  labs(
    title = "Distribution of Gene CPM with Top Gene Labels",
    x = "Average CPM (log1p scale)",
    y = "Density" # Or "Count" if you used y = ..count..
  )
ggsave("out/image/00_explore_MTReads_02.pdf",width = 12,height = 6)
