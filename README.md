# Single-cell RNA-seq Analysis

## Overview
This repository contains a single-cell RNA sequencing analysis performed using Seurat in R.
The analysis includes quality control, normalization, dimensionality reduction, clustering,
cell type annotation, and differential expression analysis.

## Dataset
- Technology: 10x Genomics 3' v3
- Tissue: Urinary bladder (organoid)
- Organism: Human
- Data format: h5ad (converted to Seurat)

## Analysis Workflow
1. Quality control and filtering
2. Normalization and variable feature selection
3. PCA and clustering stability analysis
4. UMAP visualization and cell type annotation
5. Differential expression analysis
6. Biological insight extraction

## Repository Structure
- `scripts/` – R scripts for each analysis step
- `figures/` – All figures generated in the analysis
- `results/` – Marker genes and differential expression tables

## Key Results
- Distinct epithelial and stromal populations were identified
- Stable clustering observed across PCA and resolution parameters
- Differential expression highlights lineage-specific transcriptional programs

## Software
- R
- Seurat
- ggplot2

## Author
Ayushi Singh
