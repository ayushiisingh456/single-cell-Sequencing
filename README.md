# Single-cell RNA-seq Analysis of Epithelial Tumor Cells

This repository contains a complete single-cell RNA sequencing (scRNA-seq) analysis pipeline implemented in R using the Seurat framework. The analysis focuses on epithelial tumor cell populations and follows standard best practices for quality control, clustering, cell type annotation, and differential expression analysis.

---

## Overview of Analysis

The workflow addresses the following research objectives:

1. **Quality Control (QC):**  
   Assessment and filtering of low-quality cells based on gene counts, UMI counts, and mitochondrial gene expression.

2. **Dimensionality Reduction & Clustering Stability:**  
   Principal Component Analysis (PCA) and evaluation of clustering stability across multiple resolution parameters.

3. **Marker Gene Identification & Cell Type Annotation:**  
   Identification of cluster-specific marker genes and assignment of biologically meaningful epithelial cell identities.

4. **Differential Expression Analysis:**  
   Comparative gene expression analysis between biologically relevant epithelial subpopulations.

5. **Study-linked Insight:**  
   Identification and visualization of a rare stem-like epithelial cell cluster.

---

## Repository Structure
scRNAseq-project/
│
├── script/
│ └── analysis.R # Complete scRNA-seq analysis pipeline
│
├── figures/
│ ├── Fig1_QC_violin_before.png
│ ├── Fig1_QC_violin_after.png
│ ├── Fig2_QC_scatter1.png
│ ├── Fig2_QC_scatter2.png
│ ├── Fig3_PCA_diagnostics.png
│ ├── Fig4_UMAP_clusters.png
│ ├── Fig4_UMAP_celltypes.png
│ ├── Fig5_marker_heatmap.png
│ ├── Fig5_marker_dotplot.png
│ ├── Fig5_rare_cluster.png
│ └── Fig6_DE_volcano.png
│
├── tables/
│ ├── cluster_markers.csv # Marker genes for all clusters
│ └── differential_expression.csv # DE results between selected cell groups
│
├── report/
│ └── scRNAseq_analysis_report.pdf # 2–4 page analysis report
│
├── session_info/
│ └── sessionInfo.txt # R session and package versions
│
└── README.md

---

## Methods Summary

- Input data was imported from a `.h5ad` (AnnData) file and converted into a Seurat object.
- Gene identifiers were mapped from Ensembl IDs to gene symbols.
- Cells were filtered based on:
  - Number of detected genes (`nFeature_RNA`)
  - Total UMI counts (`nCount_RNA`)
  - Percentage of mitochondrial gene expression (`percent.mt`)
- PCA was used for dimensionality reduction, with the elbow plot guiding PC selection.
- Clustering was performed using graph-based methods at multiple resolutions, with resolution 0.6 selected for downstream analyses.
- Cell types were assigned based on canonical marker gene expression.
- Differential expression analysis was performed between proliferative progenitor-like and differentiated epithelial cell populations.
- Results were visualized using UMAP, heatmaps, dot plots, and volcano plots.

---

## Software & Reproducibility

All analyses were conducted in R using Seurat and associated visualization and annotation packages.  
Exact package versions and R session details are provided in `session_info/sessionInfo.txt` to ensure reproducibility.

---

## Notes

This repository is intended for educational and research purposes and demonstrates a complete, reproducible scRNA-seq analysis workflow suitable for coursework, projects, or portfolio presentation.
