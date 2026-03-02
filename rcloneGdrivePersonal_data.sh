#!/bin/bash
. /home/edo/micromamba/bin/activate;
conda activate env_rclone;

rclone copy -PL /mnt/SSD02/work/HSR/Absinta/spatial_visium/220924_Visium_BrainMS_Martina/raw_data \
gdrive_edo:work/analysis/HSR/absinta/spatial/220924_Visium_BrainMS_Martina/data/raw_data

rclone copy -PL /mnt/SSD02/work/HSR/Absinta/spatial_visium/220924_Visium_BrainMS_Martina/analysis/GitHub/data \
gdrive_edo:work/analysis/HSR/absinta/spatial/220924_Visium_BrainMS_Martina/data/data_common

rclone copy -PL /mnt/SSD02/work/HSR/Absinta/spatial_visium/220924_Visium_BrainMS_Martina/analysis/GitHub/Visium_Seurat_5.1.0/data \
gdrive_edo:work/analysis/HSR/absinta/spatial/220924_Visium_BrainMS_Martina/data/Visium_Seurat_5.1.0/data

rclone copy -PL /mnt/SSD02/work/HSR/Absinta/spatial_visium/220924_Visium_BrainMS_Martina/analysis/GitHub/Visium_brainMS_single_sample/data \
gdrive_edo:work/analysis/HSR/absinta/spatial/220924_Visium_BrainMS_Martina/data/Visium_brainMS_single_sample/data
