# ============================
# Seurat PBMC SCT workflow 
# ============================

#_____________Normalize data using SCTransform _____________
pbmc_sct <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

#_____________PCA visualization on the SCTransform pipeline object_____________

# Run PCA on the variable features selected during SCTransform (no need to specify features)
pbmc_sct <- RunPCA(pbmc_sct, verbose = FALSE)

# Inspect top loadings (genes) for the first 5 principal components
print(pbmc_sct[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize which genes drive PCs 1 and 2
pbmc_sct_visualize_PC1_PC2_driving_genes <- VizDimLoadings(pbmc_sct, dims = 1:2, reduction = "pca")
print(pbmc_sct_visualize_PC1_PC2_driving_genes)

# Scatter plot of cells in PC space (PC1 vs PC2)
pbmc_sct_scatterplot_PC1_PC2 <- DimPlot(pbmc_sct, reduction = "pca")
print(pbmc_sct_scatterplot_PC1_PC2)

# Heatmap of the top genes for selected PCs
pbmc_sct_heatmap_selected_PCs <- DimHeatmap(pbmc_sct, dims = 1:9, cells = 500, balanced = TRUE, fast = TRUE)
print(pbmc_sct_heatmap_selected_PCs)

# Determine how many PCs to retain downstream using an elbow plot
pbmc_sct_elbow_plot <-ElbowPlot(pbmc_sct, ndims = 50)
print(pbmc_sct_elbow_plot)

#_____________Cluster the cells from the SCTransform pipeline object_____________

# Build a Shared Nearest Neighbor (SNN) graph using the first 10 PCs
pbmc_sct <- FindNeighbors(pbmc_sct, dims = 1:10)

# Perform community detection to assign cluster IDs
pbmc_sct <- FindClusters(pbmc_sct, resolution = 0.5)

# Show cluster IDs of the first 5 cells
head(Idents(pbmc_sct), 5)

#_____________Run non-linear dimensional reduction (UMAP)using the SCT-normalized pipeline object(sorry Lior)_____________

# Compute a 2D UMAP embedding using the first 10 PCs
pbmc_sct <- RunUMAP(pbmc_sct, dims = 1:10)

# Plot the UMAP embedding to visualize cells in two dimensions
pbmc_sct_unlabeled_UMAP <- DimPlot(pbmc_sct, reduction = "umap")
print(pbmc_sct_unlabeled_UMAP)

#_____________Finding differentially expressed features (cluster biomarkers) using the SCT-normalized pipeline object_____________

# Set default to the SCT assay
DefaultAssay(pbmc_sct) <- "SCT"

# Find all markers of cluster 2
cluster2_sct.markers <- FindMarkers(pbmc_sct, ident.1 = 2)
head(cluster2_sct.markers, n = 5)

# Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5_sct.markers <- FindMarkers(pbmc_sct, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5_sct.markers, n = 5)

# Find marker genes for each cluster vs. all other cells.
pbmc_sct.markers <- FindAllMarkers(pbmc_sct, only.pos = TRUE)

# Keep strong markers (avg_log2FC > 1) within each cluster
pbmc_sct.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Group the results by cluster and keep only genes with average log2 fold change > 1
cluster0_sct.markers <- FindMarkers(
  pbmc_sct,
  ident.1 = 0,
  logfc.threshold = 0.25,  # ignore tiny changes
  test.use = "roc",
  only.pos = TRUE
)

# View the top results to see the highest AUC (classification power) markers
head(cluster0_sct.markers)

# Visualizing marker expression via SCT normalized expression values 
#   - MS4A1 / CD79A – B-cell markers (B-cell receptor components)
pbmc_sct_vlnplot_ms4a1_cd79a <- VlnPlot(pbmc_sct, features = c("MS4A1", "CD79A"))
print(pbmc_sct_vlnplot_ms4a1_cd79a)

# Visualizing marker expression via raw counts
#   - NKG7 – strongly expressed in NK cells and cytotoxic T cells
#   - PF4 – platelet/megakaryocyte marker
pbmc_sct_vlnplot_nkg7_pf4 <- VlnPlot(pbmc_sct, features = c("NKG7", "PF4"), assay = "RNA", slot = "counts", log = TRUE)
print(pbmc_sct_vlnplot_nkg7_pf4)

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
pbmc_sct_immunecells_featureplot <- FeaturePlot(
  pbmc_sct,
  features = c("MS4A1","GNLY","CD3E","CD14","FCER1A","FCGR3A","LYZ","PPBP","CD8A")
)
print(pbmc_sct_immunecells_featureplot)

# DoHeatmap() generates an expression heatmap for given cells and features
# This plots the top 20 markers for each cluster
top10_sct <- pbmc_sct.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

pbmc_sct_heatmap_top10features <- DoHeatmap(pbmc_sct, features = top10_sct$gene) + NoLegend()
print(pbmc_sct_heatmap_top10features)

#_____________Automated cell-type labeling with Azimuth using the SCT-normalized pipeline object_____________

# Load an Azimuth PBMC reference (SCT-normalized single-cell atlas)
InstallData("pbmcsca")            
ref_sct <- LoadData("pbmcsca")    

# Create an SCT assay on the reference
ref_sct <- SCTransform(ref_sct, vst.flavor = "v2", verbose = FALSE)
DefaultAssay(ref_sct) <- "SCT"

# Give the reference an SCT-based PCA
ref_sct <- RunPCA(ref_sct, verbose = FALSE)

# Ensure query uses SCT
DefaultAssay(pbmc_sct) <- "SCT"

# Find transfer anchors
anchors_sct <- FindTransferAnchors(
  reference = ref_sct,  
  query = pbmc_sct,
  normalization.method = "SCT",
  reference.reduction = "pca"  
)

# Transfer reference cell type labels to the query object

# Check column names to properly choose the column in the reference that holds the cell type labels
colnames(ref_sct@meta.data)

# Use the anchor set to map cell type labels from the reference to our PBMC query.
#   - anchorset: the anchor relationships found by FindTransferAnchors()
#   - refdata:   the metadata column in the reference that holds the cell type labels
#                (here "CellType"; change if your reference uses a different column)
#   - query:     the query Seurat object we want to annotate (pbmc_sct)

labels_sct <- setNames(ref_sct@meta.data$CellType, colnames(ref_sct))

pred_sct <- TransferData(
  anchorset = anchors_sct,
  refdata   = labels_sct,
  query     = pbmc_sct
)

# Check column names 
colnames(pred_sct@meta.data)

# UMAP plot
pbmc_sct_azimuth_labeled_UMAP <- DimPlot(pred_sct, reduction = "umap", group.by = "predicted.id", label = TRUE)
print(pbmc_sct_azimuth_labeled_UMAP)
