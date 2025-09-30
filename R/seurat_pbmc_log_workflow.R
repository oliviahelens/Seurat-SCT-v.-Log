# ============================
# Seurat PBMC log workflow 
# ============================

#_____________**Normalize data using Log transform**_____________
pbmc_log <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#_____________Feature selection on the log-normalized pipeline object_____________
# Identify genes with the most variable expression across cells
#   - Uses variance-stabilizing transformation ("vst") method
#   - Keeps the top 2,000 most variable genes for downstream analysis (PCA, clustering)
pbmc_log <- FindVariableFeatures(pbmc_log, selection.method = "vst", nfeatures = 2000)

# Extract the names of the 10 most variable genes
top10_log <- head(VariableFeatures(pbmc_log), 10)

# Create a scatter plot of all genes showing mean expression vs. variance,
plot1_log_variablefeature <- VariableFeaturePlot(pbmc_log)

# A labels to highlight the top 10 variable genes
plot2_log_variablefeature <- LabelPoints(plot = plot1_log_variablefeature, points = top10_log, repel = TRUE)

# Display the labeled version
print(plot2_log_variablefeature)


#_____________Scale data for the log-normalized pipeline object_____________
# Select all genes for scaling
all.genes <- rownames(pbmc_log)

# Scale and center expression values:
#   - Mean of each gene across cells set to 0
#   - Variance set to 1
#   - Regress out mitochondrial percentage (percent.mt) to remove technical effects
pbmc_log <- ScaleData(
  pbmc_log,
  features = all.genes,
  vars.to.regress = "percent.mt"
)

# Verify that ScaleData() centered and scaled each variable gene
#   - Mean expression per gene should be ~0
#   - Standard deviation (variance^0.5) per gene should be ~1
# Get the list of variable features
var.genes <- VariableFeatures(pbmc_log)

# Extract their scaled expression matrix
mat <- GetAssayData(pbmc_log, slot = "scale.data")[var.genes, ]

# Inspect the first 5 variable genes
apply(mat[1:5, ], 1, mean)  # should be close to 0
apply(mat[1:5, ], 1, sd)    # should be close to 1

#_____________PCA visualization on the log-normalized pipeline object_____________

# Run PCA using the current variable features (from FindVariableFeatures)
pbmc_log <- RunPCA(pbmc_log, features = VariableFeatures(pbmc_log), verbose = FALSE)

# Inspect top loadings (genes) for the first 5 PCs
print(pbmc_log[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize which genes drive PCs 1 and 2
pbmc_log_visualize_PC1_PC2_driving_genes <- VizDimLoadings(pbmc_log, dims = 1:2, reduction = "pca")
print(pbmc_log_visualize_PC1_PC2_driving_genes)

# Scatter plot of cells in PC space (PC1 vs PC2)
pbmc_log_scatterplot_PC1_PC2 <- DimPlot(pbmc_log, reduction = "pca")
print(pbmc_log_scatterplot_PC1_PC2)

# Heatmap of the top genes for selected PCs
pbmc_log_heatmap_selected_PCs <- DimHeatmap(pbmc_log, dims = 1:9, cells = 500, balanced = TRUE, fast = TRUE)
print(pbmc_log_heatmap_selected_PCs)

# Choosing how many PCs to keep downstream via Elbow plot to choose a knee point 
pbmc_log_elbow_plot <- ElbowPlot(pbmc_log, ndims = 50)
print(pbmc_log_elbow_plot)

#_____________Cluster the cells from the log-normalized pipeline object_____________

# Build a Shared Nearest Neighbor (SNN) graph using the first 10 PCs
pbmc_log <- FindNeighbors(pbmc_log, dims = 1:10)

# Perform community detection on the SNN graph to assign cluster IDs
#   - The resolution parameter controls cluster granularity (higher = more clusters)
pbmc_log <- FindClusters(pbmc_log, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc_log), 5)

#_____________Run non-linear dimensional reduction (UMAP) on the log-normalized pipeline object (sorry Lior)_____________

# This computes a 2D embedding using the first 10 PCs and stores it in the 
# "umap" slot of pbmc_log
pbmc_log <- RunUMAP(pbmc_log, dims = 1:10)

# Plot the UMAP embedding to visualize cells in two dimensions
#   - Each point is a cell and colors represent clusters or identities
pbmc_log_unlabeled_UMAP <- DimPlot(pbmc_log, reduction = "umap")
print(pbmc_log_unlabeled_UMAP)

#_____________Finding differentially expressed features (cluster biomarkers) via the log-normalized pipeline object_____________

# Find all markers of cluster 2
cluster2_log.markers <- FindMarkers(pbmc_log, ident.1 = 2)
head(cluster2_log.markers, n = 5)

# Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5_log.markers <- FindMarkers(pbmc_log, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5_log.markers, n = 5)

# Find marker genes for each cluster vs. all other cells.
# only.pos = TRUE means we only keep genes that are up-regulated in the cluster
pbmc_log.markers <- FindAllMarkers(pbmc_log, only.pos = TRUE)

# Group the results by cluster and keep only genes with average log2 fold change > 1
pbmc_log.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Identify marker genes that best distinguish cluster 0 from all other cells
#   - Uses ROC test to measure each gene’s ability to classify cells as cluster 0
#   - logfc.threshold = 0.25 filters out genes with small fold changes (<0.25 log2FC)
#   - only.pos = TRUE keeps only genes up-regulated in cluster 0
cluster0_log.markers <- FindMarkers(
  pbmc_log,           
  ident.1 = 0,         # target cluster to compare against all others
  logfc.threshold = 0.25,
  test.use = "roc",
  only.pos = TRUE
)

# View the top results to see the highest AUC (classification power) markers
head(cluster0.markers)

# Visualizing marker expression via scaled/normalized expression values
#   - MS4A1 / CD79A – B-cell markers (B-cell receptor components)
pbmc_log_vlnplot_ms4a1_cd79a <- VlnPlot(pbmc_log, features = c("MS4A1", "CD79A"))
print(pbmc_log_vlnplot_ms4a1_cd79a)

# Visualizing marker expression via raw counts
#   - NKG7 – strongly expressed in NK cells and cytotoxic T cells
#   - PF4 – platelet/megakaryocyte marker
pbmc_log_vlnplot_nkg7_pf4 <- VlnPlot(pbmc_log, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
print(pbmc_log_vlnplot_nkg7_pf4)

# Visualize expression of key marker genes across the UMAP embedding
#   - FeaturePlot overlays each gene’s normalized expression on the UMAP layout
#   - Cells are colored by expression level (low = gray, high = purple/blue scale)
#   - This plots common immune cell markers:
#       MS4A1  (B cells) 
#       GNLY   (NK cells / cytotoxic T)
#       CD3E   (T cells)
#       CD14   (monocytes)
#       FCER1A (dendritic cells)
#       FCGR3A (NK/monocytes)
#       LYZ    (monocytes)
#       PPBP   (platelets)
#       CD8A   (CD8+ T cells)
pbmc_log_immunecells_featureplot <- FeaturePlot(pbmc_log, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                                   "CD8A"))
print(pbmc_log_immunecells_featureplot)

# DoHeatmap() generates an expression heatmap for given cells and features 
# This plots the top 20 markers for each cluster

pbmc_log.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_log
pbmc_log_heatmap_top10features <- DoHeatmap(pbmc_log, features = top10_log$gene) + NoLegend()
print(pbmc_log_heatmap_top10features)

#_____________Assigning cell type identity to clusters using the log-normalized pipeline object_____________
# Use canonical marker genes to rename Seurat cluster IDs with known PBMC cell types

new.cluster.ids <- c("Naive CD4 T",    # cluster 0: IL7R, CCR7 markers
                     "CD14+ Mono",     # cluster 1: CD14, LYZ
                     "Memory CD4 T",   # cluster 2: IL7R, S100A4
                     "B",              # cluster 3: MS4A1
                     "CD8 T",          # cluster 4: CD8A
                     "FCGR3A+ Mono",   # cluster 5: FCGR3A, MS4A7
                     "NK",             # cluster 6: GNLY, NKG7
                     "DC",             # cluster 7: FCER1A, CST3 (dendritic cells)
                     "Platelet")       # cluster 8: PPBP

# Map these new names to the existing cluster levels of the Seurat object
names(new.cluster.ids) <- levels(pbmc_log)

# Rename the cluster identities in the Seurat object to the descriptive cell types
pbmc_log <- RenameIdents(pbmc_log, new.cluster.ids)

# Plot the UMAP with the updated cell type labels instead of numeric cluster IDs
pbmc_log_manual_labeled_UMAP <- DimPlot(pbmc_log, reduction = "umap", label = TRUE, pt.size = 0.5)
print(pbmc_log_manual_labeled_UMAP)


#_____________Automated cell-type labeling with Azimuth using the log-normalized pipeline object_____________

# Load an Azimuth PBMC reference (SCT-normalized single-cell atlas)
InstallData("pbmcsca")            
ref_log <- LoadData("pbmcsca")    

# Prepare log-normalized view for the reference
DefaultAssay(ref_log) <- "RNA"
ref_log <- NormalizeData(ref_log, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
ref_log <- FindVariableFeatures(ref_log, verbose = FALSE)
ref_log <- ScaleData(ref_log, verbose = FALSE)
ref_log <- RunPCA(ref_log, verbose = FALSE)         

# Ensure query uses log
DefaultAssay(pbmc_log) <- "RNA"

# Find transfer anchors
anchors_log <- FindTransferAnchors(
  reference = ref_log,
  query = pbmc_log,
  normalization.method = "LogNormalize",
  reference.reduction = "pca"  
)

# Transfer reference cell type labels to the query object

# Check column names to properly choose the column in the reference that holds the cell type labels
colnames(ref_log@meta.data)

# Use the anchor set to map cell type labels from the reference to our PBMC query.
#   - anchorset: the anchor relationships found by FindTransferAnchors()
#   - refdata:   the metadata column in the reference that holds the cell type labels
#                (here "CellType"; change if your reference uses a different column)
#   - query:     the query Seurat object we want to annotate (pbmc_sct)

labels_log <- setNames(ref_log@meta.data$CellType, colnames(ref_log))

pred_log <- TransferData(
  anchorset = anchors_log,
  refdata   = labels_log,
  query     = pbmc_log
)

# Check column names 
colnames(pred_log@meta.data)

# UMAP plot
pbmc_log_azimuth_labeled_UMAP <- DimPlot(pred_log, reduction = "umap", group.by = "predicted.id", label = TRUE)
print(pbmc_log_azimuth_labeled_UMAP)
