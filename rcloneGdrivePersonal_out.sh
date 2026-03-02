#!/bin/bash
. /home/edo/micromamba/bin/activate;
conda activate env_rclone;

rclone copy -PL /mnt/SSD02/work/HSR/Absinta/spatial_visium/220924_Visium_BrainMS_Martina/analysis/GitHub/out_large \
gdrive_edo:work/analysis/HSR/absinta/spatial/220924_Visium_BrainMS_Martina/out/out_large

rclone copy -PL /mnt/SSD02/work/HSR/Absinta/spatial_visium/220924_Visium_BrainMS_Martina/analysis/GitHub/Visium_Seurat_5.1.0/out \
gdrive_edo:work/analysis/HSR/absinta/spatial/220924_Visium_BrainMS_Martina/out/Visium_Seurat_5.1.0/out

rclone copy -PL /mnt/SSD02/work/HSR/Absinta/spatial_visium/220924_Visium_BrainMS_Martina/analysis/GitHub/Visium_brainMS_single_sample/out \
gdrive_edo:work/analysis/HSR/absinta/spatial/220924_Visium_BrainMS_Martina/out/Visium_brainMS_single_sample/out
