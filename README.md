# multiomic-spatial-integration
Integrated single-nucleus RNA-seq and GeoMx WTA spatial transcriptomics workflow for cell-type-resolved spatial deconvolution.

🚧 Work in progress

![R](https://img.shields.io/badge/R-4.x-blue)
![Python](https://img.shields.io/badge/Python-3.x-yellow)
![Seurat](https://img.shields.io/badge/Seurat-v5-green)
![Pyro](https://img.shields.io/badge/Pyro-Bayesian%20modeling-orange)
![Status](https://img.shields.io/badge/status-in%20progress-orange)

---

## Overview

This repository contains a modular multi-omic workflow for integrating annotated single-nucleus RNA-seq (snRNA-seq) data with NanoString GeoMx Whole Transcriptome Atlas (WTA) Digital Spatial Profiling data.

The goal is to estimate cell-type-resolved spatial abundance in human brain tissue and test how inferred cellular composition varies by disease status, amyloid pathology, and vascular/parenchymal microenvironment.

The example biological use case focuses on AD/CAA versus control human brain tissue, with emphasis on vascular and parenchymal amyloid-associated remodeling.

---

## Workflow Summary

The workflow consists of five major stages:

1. snRNA-seq preprocessing, QC, clustering, subclustering, and cell-type annotation  
2. GeoMx WTA processing, QC, and AnnData object construction  
3. Regression-based cell-type signature inference and Bayesian spatial deconvolution
   Note: training signatures separately for AD vs Control, however, If there is only a unified reference. this can also be made with this repo.
5. Extraction, normalization, and annotation of spatial cell-type abundances  
6. Statistical analysis of inferred cell-type proportions and spatial visualization  

This structure follows the integrated workflow described in the internal project documentation. :contentReference[oaicite:0]{index=0}

---

## Repository Structure

```text
multiomic-spatial-integration/
├── R/
│   ├── seurat_reference_utils.R          # Prepare annotated snRNA-seq reference
│   ├── signature_export_utils.R          # Export cell-type signatures
│   ├── proportion_stats_utils.R          # Beta mixed models for inferred proportions
│   ├── contrast_utils.R                  # Disease, amyloid, and weighted contrasts
│   ├── enrichment_utils.R                # Pathway sanity-check inputs
│   └── plotting_utils.R                  # Heatmaps, dotplots, boxplots
│
├── python/
│   ├── geomx_anndata_utils.py            # Build AnnData from GeoMx WTA matrices
│   ├── regression_utils.py               # Cell2Location regression helpers
│   ├── spacejam_pyro_utils.py            # Pyro/SpaceJam model helpers
│   ├── abundance_extraction_utils.py     # Extract and normalize spot factors
│   └── plotting_utils.py                 # Exploratory spatial plots
│
├── models/
│   └── LocationModelWTAMultiExperimentHierarchicalGeneLevel_Modified.py
│
├── notebooks/
│   ├── 01_wta_to_anndata.ipynb
│   ├── 02_regression_signatures.ipynb
│   ├── 03_spacejam_cell2location.ipynb
│   ├── 04_extract_cell_proportions.ipynb
│   └── 05_plotting_and_stats.ipynb
│
├── scripts/
│   ├── run_sn_reference_export.R
│   ├── run_spatial_stats.R
│   └── run_pathway_sanity_checks.R
│
├── docs/
│   ├── workflow_overview.md
│   ├── geomx_anndata_structure.md
│   ├── bayesian_deconvolution_notes.md
│   ├── pyro_model_modifications.md
│   └── statistical_modeling.md
│
├── README.md
└── .gitignore
```

## Stage 1: snRNA-seq Reference Preparation
Annotated snRNA-seq data are used as the cellular reference for spatial deconvolution.
Main steps include:

- Ambient RNA correction when appropriate
- QC filtering
- doublet removal
- SCTransform normalization
- Harmony integration
- clustering and manual annotation
- targeted subclustering of astrocytes, microglia, and vascular cells
- Pearson correlation-based merging of highly similar subclusters export of annotated reference signatures


The reference includes major brain cell populations such as astrocytes, microglia, oligodendrocytes, OPCs, excitatory/inhibitory neurons, and vascular-associated populations. 

## Stage 2: GeoMx WTA AnnData Construction
GeoMx WTA matrices are converted into a unified AnnData object.
The AnnData structure preserves:
- raw WTA gene counts
- negative probe counts
- ROI metadata
- disease status
- pathology status
- anatomical region
- scan/slide identifiers
- nuclei counts


Raw counts are stored in:

```Python
adata.layers["counts"]
```
Negative probe counts are stored in:
```Python
adata.obsm["negProbes"]
```
This structure ensures compatibility with downstream Bayesian spatial modeling. 


## Stage 3: Bayesian Spatial Deconvolution
The workflow uses a two-stage modeling strategy:

1.- Regression-based inference of cell-type signatures from snRNA-seq

2.- Bayesian spatial deconvolution of GeoMx WTA data using a SpaceJam / Cell2Location-style model implemented in Pyro


The spatial model incorporates:

- inferred snRNA-seq signatures
- GeoMx WTA counts
- negative probe background
- ROI nuclei counts
- scan/slide-level batch structure


The output is a cell-type abundance estimate per ROI. 

Model Adaptation
The spatial model was adapted from a PyMC3/Theano implementation into a Pyro/PyTorch framework to support:
- GPU acceleration
- stochastic variational inference
- WTA-scale gene panels
- explicit negative probe modeling
- nuclei-count-informed priors
- stable Negative Binomial likelihoods
- memory-safe saving of posterior outputs

These modifications preserve the conceptual SpaceJam structure while enabling scalable GeoMx WTA analysis. 


## Stage 4: Abundance Extraction and Annotation
After model training, posterior spatial factors are extracted and converted into:

- absolute abundance estimates
- relative abundance estimates
- ROI × cell-type long-format tables


These outputs are merged with GeoMx metadata and factor-to-cell-type annotation mappings for downstream visualization and statistical testing. 

## Stage 5: Statistical Modeling
Spatial cell-type proportions are modeled using beta mixed-effects models in R.
The core model structure is:
``` R
rel_abundance ~ contrast_variable + (1 | Scan_ID)
```
using:
``` R
glmmTMB::beta_family(link = "logit")
```
This is appropriate because inferred cell-type abundances are proportions bounded between 0 and 1.
Core contrasts include:

1.- Disease effect: AD/CAA versus control in amyloid-free ROIs


2.- Amyloid effect: amyloid-positive versus amyloid-free ROIs within AD/CAA


3.- Max Amyloid Pathology:  amyloid-positive in AD/CAA versus amyloid-free ROIs within Controls


4.- Weighted AD overall effect: pathology-weighted AD/CAA versus control


Multiple testing correction is performed across cell types using Benjamini-Hochberg FDR. 

---

## Outputs
Common outputs include:
```
*_inferred_signatures.csv
*_spot_factors_abs.pt
*_spot_factors_rel.pt
*_roi_celltype_abundance_long.csv
*_factor_celltype_annotation.csv
*_beta_mixed_model_results.csv
*_contrast_heatmap.pdf
```
Example Use Case
This repository was developed for AD/CAA spatial multi-omic integration in human brain tissue.
The workflow links:
- transcriptionally defined snRNA-seq cell states
- GeoMx WTA spatial ROIs
- vascular and parenchymal pathology annotations
- Bayesian-inferred spatial cell-type composition
- beta mixed-effects statistical inference


This enables compartment-resolved analysis of gliovascular and neuroimmune remodeling in AD/CAA.

### Notes

- This repository is intended as a modular research framework, not a single-click pipeline.
- Sensitive raw data and patient metadata should not be committed.
- Large model outputs should be stored outside Git or tracked with Git LFS.
- The modified Pyro model is intended for research use and should be validated for each dataset.



### Author
Enrique Chimal
PhD Candidate – Medical Neuroscience

This README positions the repo as a full **computational biology framework**, not just a collection of scripts.
