# Microglia_subtypes_mic_scRNA :  Identification of Dysregulated Microglia Subtypes Induced by Antiretroviral Drug AZT Using Single-Cell RNA Sequencing. [Manuscript in preparation]

## Developers 

 - Mustafa AbuElqumsan (ORCID [0000-0002-1018-1410](https://orcid.org/0000-0002-1018-1410))  
 - Xin Liu, 
 - Shao-Jun Tang ( ORCID [0000-0002-6076-5481] (https://orcid.org/0000-0002-6076-5481))

## Description    
scRNA_seq on microglia from wild-derived HIV's mouse model , Identification and Characterization of Microglial Subtypes from scRNA-seq Data :

The Microglia_subtypes_mic_scRNA R package provides a reproducible workflow for processing, analyzing, and visualizing single-cell RNA sequencing (scRNA-seq) data of microglia, with a specific focus on understanding antiretroviral drug (AZT)-induced neuropathic pain.


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

## Installation      

```
#### Install from GitHub

if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("elqumsan/Microglia_subtypes_mic_scRNA")
```
```
#### Usage
library(Microglia_subtypes_mic_scRNA)
```
```
##### Load example data
data_dir <- system.file("extdata", package = "Microglia_subtypes_mic_scRNA")
```
```
## Run complete workflow
results <- run_microglia_analysis(
    data_path = data_dir,
    integration_method = "Harmony",
    annotation_method = "SingleR"
)
```
```
## Visualize results
plot_umap(results)
plot_subtype_abundance(results)
```
## Requirements
- R (>= 4.5.1)
- Key dependencies: Seurat, SingleR, scCATCH, harmony, liger, DoubletFinder, SoupX
  For a full list of dependencies, see DESCRIPTION file.

## DataSet Information 
THe package is designed to process scRNA-seq datasets derived from microglia.
Example datasets are provided for demonstration and testing purposes

## Running All Analyses
To reproduce the entire analysis pipeline from the manuscript 

```
source(system.file("scripts", "full_analysis.R", package = "Microglia_subtypes_mic_scRNA"))
```
## Materials & Methods
Detailed analysis steps, QC thresholds, and parameter setting are described in the supplementary materials section of the manuscript.

## Citation 
If you use this package, please cite:
AbuElqumsan M, Li X, Tang S-J. Identification of Dysregulated Microglia Subtypes Induced by Antiretroviral Drug AZT Using Single-Cell RNA Sequencing. [Manuscript in preparation]

## License
????
## contact 
For questions or bug reports, please contact:
- Mustafa AbuElqumssan --- mustafa.abuelqumsan@stonybrookmedicine.edu


