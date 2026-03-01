scRNAseq-analysis/
│
├── README.md
│
├── data/
│   ├── raw/            # (optional, if allowed)
│   └── processed/      # filtered Seurat object
│
├── scripts/
│   ├── 01_QC.R
│   ├── 02_Normalization_PCA.R
│   ├── 03_Clustering_UMAP.R
│   ├── 04_Markers_CellType.R
│   └── 05_DifferentialExpression.R
│
├── figures/
│   ├── Fig1_QC_violin.png
│   ├── Fig2_QC_scatter.png
│   ├── Fig3_PCA.png
│   ├── Fig4_UMAP.png
│   ├── Fig5_Markers.png
│   └── Fig6_Volcano.png
│
├── results/
│   ├── markers_all_clusters.csv
│   └── DE_epithelial_vs_fibroblast.csv
│
└── session_info.txt
