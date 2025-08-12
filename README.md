# Microglia_subtypes_mic_scRNA 
scRNA_seq on microglia from wild-derived HIV's mouse model

Title:
Microglia_subtypes_mic_scRNA: Identification and Characterization of Microglial Subtypes from scRNA-seq Data

Description:
The Microglia_subtypes_mic_scRNA R package provides a reproducible workflow for processing, analyzing, and visualizing single-cell RNA sequencing (scRNA-seq) data of microglia, with a specific application to studies on antiretroviral drug (AZT)-induced neuropathic pain.

The package implements end-to-end scRNA-seq analysis, including:

Raw Data Processing: Import of Cell Ranger outputs (.h5, matrix.mtx) with automated QC, doublet removal, and filtering of low-quality cells.

Normalization & Integration: Single-cell–specific normalization methods (SCTransform, log-normalization) and batch integration using Seurat, Harmony, or Liger.

Dimensionality Reduction & Clustering: PCA, UMAP/t-SNE projections, and community detection for identifying microglia subpopulations.

Cell Type Annotation: Semi-automated annotation using marker gene databases (SingleR, scCATCH) and curated microglia-specific markers.

Differential Expression & Functional Enrichment: DEG identification between conditions (WT vs. AZT), pain-related gene prioritization, and pathway analysis (IPA, GSEA).

Publication-Ready Visualizations: Violin plots, dot plots, cluster heatmaps, and subtype abundance bar charts.

Target Audience:
Researchers studying microglia biology, neuroinflammation, neuropathic pain, HIV-associated neurological disorders, or anyone interested in scRNA-seq workflows for CNS immune cell populations.


