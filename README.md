# Microglia_subtypes_mic_scRNA 
scRNA_seq on microglia from wild-derived HIV's mouse model

Identification and Characterization of Microglial Subtypes from scRNA-seq Data

The Microglia_subtypes_mic_scRNA R package provides a reproducible workflow for processing, analyzing, and visualizing single-cell RNA sequencing (scRNA-seq) data of microglia, with a specific focus on understanding antiretroviral drug (AZT)-induced neuropathic pain.

This package was developed for the analyses presented in:

Elqumsan M, Li X, Tang S-J. Identification of Dysregulated Microglia Subtypes Induced by Antiretroviral Drug AZT Using Single-Cell RNA Sequencing. [Manuscript in preparation].


## Features
- Raw Data Processing
  - Import Cell Ranger output (.h5, .mtx)
  - Quality control: filter low-quality cells and remove doublets

- Normalization & Integration
   - SCTransform-based normalization
   - Batch correction via Harmony or Liger

- Dimensional Reduction & Clustering
   - PCA, UMAP
   - Graph-based clustering to identify subpopulations

- Cell Type Annotation
   -  Marker-based annotation
   - Automated annotation via SingleR and scCATCH

- Differential Expression & Functional Enrichment
   - Subtype-specific DEG analysis
   - Pathway enrichment via IPA and GSEA

- Publication-Ready Visualizations
   - UMAP plots, violin plots, dot plots, heatmaps, and subtype abundance barplots
