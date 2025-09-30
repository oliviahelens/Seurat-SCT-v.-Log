============================
Seurat PBMC SCT v. log 
============================

Purpose:
Compare Seurat log-normalization with SCTransform (SCT) on the standard PBMC3k dataset.  
This repo follows the Seurat tutorial flow with added steps:
- Unique object names so both pipelines can run in one R session
- Manual labels (on log-normalized data only) and Azimuth auto-labels for side-by-side comparison

Environment info:
 - R: 4.4.3  
 - Seurat: 5.3.0  
 - SeuratData: 0.2.2.9002  
 - glmGamPoi: 1.18.0  
 - dplyr: 1.1.4  
 - ggplot2: 3.5.2  
 - sctransform: 0.4.2  
 - patchwork: 1.3.2

Minimal install:
install.packages(c("Seurat","SeuratObject","glmGamPoi","dplyr","ggplot2","patchwork")

Repository layout: 
R/ 
  seurat_pbmc_main.R 
  seurat_pbmc_log_workflow.R 
  seurat_pbmc_SCT_workflow.R 
graphs/              		(selected output figures) 
data/filtered_gene_bc_matrices/hg19/
<project>.Rproj 
README.txt

Start here: 
Open the entry script: seurat_pbmc_main.R 
At the bottom of that file set one of: 
  mode <- "log" 
  mode <- "SCT" 
Run the script (source it or run line-by-line). The script loads libraries and data, 
does QC + preprocessing, then calls the matching workflow script.