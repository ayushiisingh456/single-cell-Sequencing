# ==========================================
# Single-cell RNA-seq Analysis
# Refined and Organized Script
# ==========================================

# -------------------------------
# Load Libraries
# -------------------------------
if(!require(Seurat)) install.packages("Seurat")
if(!require(anndata)) install.packages("anndata")
if(!require(org.Hs.eg.db)) BiocManager::install("org.Hs.eg.db")
if(!require(AnnotationDbi)) BiocManager::install("AnnotationDbi")
if(!require(dplyr)) install.packages("dplyr")
if(!require(EnhancedVolcano)) install.packages("EnhancedVolcano")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggrepel)) install.packages("ggrepel")

library(anndata)
library(Seurat)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)

# -------------------------------
# Q1: Quality Control (QC)
# -------------------------------

# Load data
ad <- read_h5ad("C:/Users/Lenovo/OneDrive/Desktop/Cancer_Project/data/57c2a91b-97ee-4b88-8b8b-5806e5f9e3be.h5ad")
seurat_obj <- CreateSeuratObject(
  counts = t(ad$X),
  meta.data = ad$obs
)

# Convert gene IDs to symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(seurat_obj),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
valid <- !is.na(gene_symbols)
seurat_obj <- seurat_obj[valid, ]
rownames(seurat_obj) <- gene_symbols[valid]

# QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Violin plot before filtering
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/Fig1_QC_violin_before.png")

# Scatter plots
FeatureScatter(seurat_obj, "nCount_RNA", "percent.mt")
ggsave("figures/Fig2_QC_scatter1.png")
FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA")
ggsave("figures/Fig2_QC_scatter2.png")

# Filtering
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                              nFeature_RNA < 6000 &
                              percent.mt < 15)

# Violin plot after filtering
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("figures/Fig1_QC_violin_after.png")

# -------------------------------
# Q2: Cluster Stability (PCA & Resolution)
# -------------------------------

# Normalize and scale
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# PCA
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)
ggsave("figures/Fig3_PCA_diagnostics.png")

# Clustering at different resolutions
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
resolutions <- c(0.4, 0.6, 0.8)
for (res in resolutions) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
  DimPlot(seurat_obj, label = TRUE) +
    ggtitle(paste("Resolution:", res)) 
}

# UMAP visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
ggsave("figures/Fig4_UMAP_clusters.png")

# -------------------------------
# Q3: Marker Identification & Cell Type Assignment
# -------------------------------

# Find all markers
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(markers, "tables/cluster_markers.csv", row.names = FALSE)

# Top 10 markers per cluster heatmap
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)
DoHeatmap(seurat_obj, features = top10$gene, group.by = "seurat_clusters") +
  NoLegend()
ggsave("figures/Fig5_marker_heatmap.png")

# DotPlot for key markers
DotPlot(seurat_obj, features = c("EPCAM", "CD3D", "MS4A1"))
ggsave("figures/Fig5_marker_dotplot.png")

# Rename clusters with cell type identities
new_ids <- c(
  "Neoplastic_Epithelial_Cells",
  "Progenitor_Like_Cells",
  "Proliferative_Epithelial_Progenitor",
  "Quiescent_Epithelial_Cells",
  "Differentiated_Epithelial_Cells",
  "Cycling_Epithelial_Cells",
  "Stem_Like_Epithelial_Cells",
  "Stress_Resilient_Epithelial_Cells",
  "EMT_Like_Epithelial_Cells"
)
names(new_ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new_ids)
DimPlot(seurat_obj, label = TRUE)
ggsave("figures/Fig4_UMAP_celltypes.png")

# -------------------------------
# Q4: Differential Expression (DE)
# -------------------------------

Idents(seurat_obj) <- "seurat_clusters"
de_results <- FindMarkers(seurat_obj,
  ident.1 = "Proliferative_Epithelial_Progenitor",
  ident.2 = "Differentiated_Epithelial_Cells"
)

# Save DE table
dir.create("tables", showWarnings = FALSE)
write.csv(de_results, "tables/differential_expression.csv")

# Volcano plot using EnhancedVolcano
EnhancedVolcano(
  de_results,
  lab = rownames(de_results),
  x = "avg_log2FC",
  y = "p_val_adj",
  pCutoff = 0.05,
  FCcutoff = 0.25,
  pointSize = 2,
  labSize = 3
)
ggsave("figures/Fig6_DE_volcano.png")

# Custom ggplot volcano with top genes labeled
de_results$negLogP <- -log10(de_results$p_val_adj)
logFC_cutoff <- 0.25
pval_cutoff <- 0.05

top_genes <- de_results %>%
  filter(p_val_adj < 0.01 & abs(avg_log2FC) > 1)

ggplot(de_results, aes(x = avg_log2FC, y = negLogP)) +
  geom_point(color = "grey", alpha = 0.5) +
  geom_point(data = top_genes, color = "red") +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), size = 3) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
  labs(title = "Differential Expression Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()
ggsave("figures/Fig6_DE_volcano_custom.png")

# -------------------------------
# Q5: Study-linked Insight
# Example: Highlight rare cluster
# -------------------------------

rare_cluster <- subset(seurat_obj, idents = "Stem_Like_Epithelial_Cells")
DimPlot(rare_cluster, reduction = "umap", label = TRUE)
ggsave("figures/Fig5_rare_cluster.png")

# -------------------------------
# Save session info
# -------------------------------
dir.create("session_info", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "session_info/sessionInfo.txt")
