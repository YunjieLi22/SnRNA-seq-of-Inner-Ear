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
  - runUmap
  - combination of the data of different ages



