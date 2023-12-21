# Single_Cell_RNA_Analysis
### Overview 
A single cell RNA analysis using data obtained from Hong et al, 2022 paper - **Cure of syngeneic carcinomas with targeted IL-12 through obligate reprogramming of lymphoid and myeloid immunity** (https://pubmed.ncbi.nlm.nih.gov/35260537/) The focus is on T/NK cells. This analyses include re-normalising data, performing clustering and dimensional reduction as well as marker identification.

# Libraries
| Data Manipulation | Single Cell Analysis | Visualization |
|-------------------|----------------------|---------------|
| dplyr             | Seurat               | ggplot2       |
| tidyr             | edgeR                | viridis       |
| stringr           | scCustomize          | Patchwork     |
|                   | limma                | RColorBrewer  |
|                   | statmod              | UpSetR        |
|                   |                      | Pheatmap      |
|                   |                      | Glimma        |


# Single Cell RNA Analysis Workflow

## 1. Data Importation and Initial QC
- Import single-cell datasets.
- QC Metrics to assess the quality of the data in the initial stages.

### QC Metrics and Filtration

1. **nFeature_RNA**: Number of genes per cell. Cells with low counts may be of poor quality or empty droplets; high counts may indicate doublets.
2. **nCount_RNA**: Total number of reads per cell.
3. **percent.mt**: Percentage of mitochondrial genes. High mitochondrial content can indicate cell stress or death.
4. **percent.ribo**: Percentage of transcripts that are ribosomal.
5. **Cells cycle**

## 2. Initial Normalization
- NormalizeData
- FindVariableFeatures
- SCTransform Normalization
- Cell Type Scoring
- FeatureScatter

## 3. Dimensionality Reduction
- RunPCA: Principal component analysis to reduce dimensionality.
- ElbowPlot: Determine the number of significant principal components.
- FindNeighbors & FindClusters: Cluster cells based on their PCA scores.
- RunUMAP: Visualize clusters using Uniform Manifold Approximation and Projection.

## 4. Gene Expression and Marker Detection
- FeaturePlot
- FindAllMarkers
- FindMarkers

## 5. Visualization
- Dimplot
- VlnPlot
- RidgePlot
- DotPlot

# Pseudobulk Analysis

## Overview

This workflow outlines the steps for conducting a Pseudobulk analysis using `Seurat` and `edgeR` libraries.

## Steps

### 1. PB Generation
- Details for PB generation go here.

### 2. Create a DGEList
- `DGEList`: Function to create an object to store the read counts and experimental data.

### 3. Filter Lowly Expressed Genes
- `filterByExpr`: Function to filter out lowly expressed genes.

### 4. QC Metrics
- Perform Quality Control using:
  - `logCPM`
  - `RLE`
  - `MDS`
  - `PCA`

### 5. Design Matrix for Expression of Selected Genes
- `model.matrix`: Function to design the matrix for differential expression analysis.

### 6. Estimate Variance
- `voom`: Function to estimate the mean-variance relationship.

### 7. Linear Modelling to Assess Differential Expression (DE)
- `lmFit`: Fits linear models to expression data.
- `makeContrasts`: Creates contrasts for differential expression.
- `contrasts.fit`: Fits the contrasts.
- `eBayes`: Empirical Bayes statistics to assess DE.
- `topTable`: Table of the top-ranked DE genes.

### 8. Visualisation with Glimma
- `glMDSPlot`: Multidimensional scaling plot.
- `glVolcano`: Volcano plot for DE genes.
- `glMAPlot`: MA-plot (M (log ratio) vs A (mean average)).

## Installation

To install the required libraries, use the following R commands:


