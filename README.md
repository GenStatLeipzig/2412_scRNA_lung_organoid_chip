# Mechanical strain exacerbates Pseudomonas infection in an organoid-based pneumonia-on-a-chip model

This repository contains code for single-cell RNA sequencing data from a study investigating how mechanical strain affects *Pseudomonas aeruginosa* infection in an organoid-based lung-on-chip model.


## Study Overview

We developed a pneumonia-on-a-chip (POC) model using human primary pulmonary microvascular endothelial cells (HPMVECs) and organoid-derived alveolar epithelial cells (ODAECs) to study ventilator-associated pneumonia. The model incorporates physiological (5%) and hyperphysiological (10%) mechanical strain to mimic breathing conditions.

## Repository Structure

```
├── scripts/
│   ├── 01_setup.R                    # Project setup and directories
│   ├── 03_nf_scrnaseq.R             # nf-core/scrnaseq pipeline execution
│   ├── 04_scDblFinder.qmd            # initial doublet estimation, extended later with demuxify
│   ├── 06_clusterbasedQC.qmd        # Quality control and integration
│   ├── 08_markerdotplot_rpca.R        # celltype identification support
│   ├── 10_integrate_with_karen.qmd   # Integration with published organoid data
│   ├── 11_markerplot_integrated_with_karen.R  # Marker gene visualization
│   ├── 13_1_integrate_with_karen_higher_resolution.qmd   # Integration with unstretched data
│   ├── 18_check_doublet.R           # Doublet detection and consensus
│   ├── 21_experiments_vs_celltype_barplot.R # Cell type distribution
│   └── Figure_Plots.R # reproduction of manuscript figures
├── demuxify/                  # doublet detection
└── README.md

```

## Methods Overview

### Single-cell RNA Sequencing Pipeline
- **Platform:** 10x Genomics Chromium Next GEM Single Cell 3' v3.1
- **Processing:** nf-core/scrnaseq pipeline with CellRanger
- **Analysis:** Seurat v5 in R with multiple integration methods (RPCA, Harmony, scVI)
- **Quality Control:** Comprehensive doublet detection using multiple algorithms implemented in demuxify
- **Cell Type Annotation:** GPT-4.1 assisted annotation with manual curation


## Requirements

- R ≥4.0 with packages: Seurat, data.table, ggplot2, patchwork
- Nextflow and nf-core/scrnaseq pipeline
- Python environment with scVI for integration (optional)

## Citation

Hoffmann K, Behrendt U, Pennitz P, et al. Mechanical strain exacerbates *Pseudomonas* infection in an organoid-based pneumonia-on-a-chip model.  (in review).



---

**Note:** This repository contains bioinformatic analysis code. Wet lab protocols and chip fabrication methods are detailed in the manuscript.
