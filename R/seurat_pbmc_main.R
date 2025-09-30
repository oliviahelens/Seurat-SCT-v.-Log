# ============================
# Seurat PBMC main 
# ============================

# Purpose:
# This script is the single entry point for running either the
# log-normalized or SCTransform (SCT) Seurat workflows with the standard 
# PBMC3k dataset. It runs through the tasks specified in the Seurat 
# tutorial with several additional steps.

# Notes:
#   - This script loads shared libraries and PBMC data, then runs standard
#      preprocessing to filter the data. 
#   - Based on `mode` entered at the end of the file, it calls either
#        source("seurat_pbmc_log_workflow.R")
#     or source("seurat_pbmc_SCT_workflow.R")
#     to execute the chosen pipeline.
#   - You can do this manually instead, if you want.
#   - Objects from each workflow are named uniquely so you can run both 
#     pipelines in the same R session if you like.
#   - This library also includes both the manual data labels for the pbmc_log
#     UMAP and the Azimuth automatic labeling for both pbmc_log and pbmc_sct for 
#     comparison.

# Environment info:
#   R: 4.4.3
#   Seurat: 5.3.0
#   SeuratData: 0.2.2.9002
#   glmGamPoi: 1.18.0
#   dplyr: 1.1.4
#   ggplot2: 3.5.2
#   sctransform: 0.4.2
#   patchwork: 1.3.2

# Load libraries
library(Seurat)
library(SeuratData)
library(glmGamPoi)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Examine a three genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# Dense versus Sparse sizing

# Calculate dense size
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

# Calculate sparse size
sparse.size <- object.size(pbmc.data)
sparse.size

# Calculate ratio of dense vs. sparse size
dense.size/sparse.size

#_____________Quality control_____________
# Calculate the percentage of reads mapping to mitochondrial genes for each cell
# and store it as a new metadata column ("percent.mt") in the Seurat object.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = .1, ncol = 3)

# QC scatterplots:
# Titles are Pearson correlation coefficients
# Left: percent.mt (percent mitochondrial RNA) vs. nCount_RNA (total RNA molecules per cell).
#       Used to spot cells with high mitochondrial content (potentially stressed or dying).
# Right: nFeature_RNA (number of detected genes) vs. nCount_RNA.
#        Shows the relationship between total molecules and detected genes;
#        points far from the diagonal can indicate doublets or technical artifacts.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)

# Filter out low-quality cells:
#   - Keep cells with more than 200 detected genes (nFeature_RNA > 200)
#   - Exclude potential doublets with more than 2500 genes (nFeature_RNA < 2500)
#   - Remove cells with >5% mitochondrial transcripts (percent.mt < 5)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#_____________Choose workflow_____________
mode <- "SCT"   # "SCT" or "log"
if (mode == "SCT") {
  source("seurat_pbmc_SCT_workflow.R")
} else {
  source("seurat_pbmc_log_workflow.R")
}
