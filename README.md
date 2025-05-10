# Single-nucleus Profiling Reveals Cell Type-specific Transcriptomic Signatures of Inner Ear Sensory Organs across Aging

## Preamable
- This directory contains illustration of the overall workflow of single-nucleus profiling of inner ear across aging.
- Samples of the sequencing data were sensory organs including cochleae and utricles from the C57BL/6 strains of 3m, 12m and 24m.
- All the scripts were runned on Linux and R.4.0.3/R.4.2.0/R.4.2.2, unless specially mentioned.
- The snRNA-seq data can be assessed on the GEO database through the accession number, GSE274279.
- Code contributors: Feng Zhang and Yunjie Li

## Code
### Cellranger 
````unix
step01_cellranger_snRNA.sh
````

### QC, analysis and combination
````unix
step02_QC_analysis/
````
- Scripts in this folder hold the function of following steps:
  - quality control
  - normalisation through SCTransform
  - clustering through runUmap
  - combination of the data of different ages
 
### Only utricle analysis
````unix
step03_Only_Utricle/
````
- Scripts in this folder perform an overall analysis of only utricle data:
  - clustring umap
  - clusters identification
  - markers illustration

### Only cochlea analysis
````unix
step04_only_cochlea/
````
- Scripts in this folder perform an overall analysis of only cochlea data:
  - clustring umap
  - clusters identification
  - markers illustration

### Combine analysis
````unix
step05_combine_analysis.R
````
- This script includes the workflow of analyzing the combination of cochlea and utricle data:
  - integrating two data using 'BEER' ('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
 
### HC aging analysis
````unix
step06_HC_aging/
````
- Scripts in this folder depict the overall workflow how we presented the transcriptional features of cochlear and utricular epithelia across ages:
  - DEG analysis across three ages of both cochlea and utricle
  - DEG analysis between two types HC of cochlea
  - DEG analysis across four types HC of utricle

### Identity of Cell Types
````unix
step07_ReCluster.R
step08_rename_clusters.R
````
- These two scripts re-adjust the clusters and update cluster IDs.

### Macrophages analysis
````unix
step09_Macrophages_aging/
````
- Scripts in this folder illustrate the workflow focusing on the analysis of macrophages of both cochlea and utricle.

### Figure adjustment
````unix
step10_figure_adjust.R
````
- This script presents the detailed adjustment of most of the figures including: UMAP, heatmaps, dotplots, violin plots, volcano plots and so on.
