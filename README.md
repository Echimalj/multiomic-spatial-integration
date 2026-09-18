# Gene Sextant: Multi-omic Spatial Integration
A modular computational framework for transforming single-cell and spatial transcriptomic data into cellular, statistical, mechanistic, and translational evidence.

🚧 Work in progress (Current version 0.5.0)

![R](https://img.shields.io/badge/R-4.x-blue)
![Python](https://img.shields.io/badge/Python-3.x-yellow)
![Seurat](https://img.shields.io/badge/Seurat-v5-green)
![Pyro](https://img.shields.io/badge/Pyro-Bayesian%20modeling-orange)
![Status](https://img.shields.io/badge/status-in%20progress-orange)

---

## Overview

**Gene Sextant** is a modular computational framework for integrating single-nucleus RNA sequencing (snRNA-seq) with spatial transcriptomics to resolve how molecular and cellular changes vary across anatomical and pathological microenvironments.

The framework was developed using NanoString GeoMx Whole Transcriptome Atlas (WTA) data and human snRNA-seq references, with an initial application to Alzheimer's disease (AD) and cerebral amyloid angiopathy (CAA).

Rather than treating spatial transcriptomic analysis as a single endpoint, Gene Sextant progressively transforms spatial molecular measurements into multiple layers of biological evidence:

**Data → Cellular Composition → Statistical Evidence → Mechanistic Hypotheses → Candidate Prioritization**

The central goal is to move beyond asking:

> **What changes in disease?**

toward asking:

> **Where does it change, which cells are associated with the change, under
> which pathological conditions does it occur, and how can that evidence
> inform biological interpretation and candidate prioritization?**

---

## Gene Sextant Architecture

Gene Sextant is organized into five conceptual modules.

```text
Human snRNA-seq                              Spatial transcriptomics
 Nuclei × Genes                                  ROIs × Genes
       │                                              │
       └──────────────────────┬───────────────────────┘
                              ▼
                 ┌─────────────────────────┐
                 │ MODULE 1                │
                 │ QC & PREPROCESSING      │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ MODULE 2                │
                 │ CELLULAR COMPOSITION    │
                 │ ROIs × Cell types       │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ MODULE 3                │
                 │ STATISTICAL EVIDENCE    │
                 │ Contrast × Feature      │
                 │ × Region                │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ MODULE 4                │
                 │ MECHANISTIC             │
                 │ INTERPRETATION          │
                 │ Cells • Genes           │
                 │ Pathways • Networks     │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ MODULE 5                │
                 │ BIOMARKER & TARGET      │
                 │ PRIORITIZATION          │
                 └────────────┬────────────┘
                              │
                              ▼
                 Spatially contextualized
                    biological evidence
                              │
                              ▼
                    Precision neuroscience
```

The modules are conceptually sequential but computationally modular. 
Individual analyses can be run independently once their required upstream data products are available.

---
## Module 1 - QC & Preprocessing

Question: Are the single-cell and spatial datasets suitable for integration?

Inputs
- annotated snRNA-seq reference data
- GeoMx WTA count matrices
- negative-probe measurements
- ROI metadata
- nuclei counts
- disease, pathology, anatomical-region, and scan identifiers


Approach

The snRNA-seq reference undergoes QC, clustering, annotation, and targeted subclustering. GeoMx WTA measurements are converted into a standardized AnnData representation containing expression counts, negative probes, and ROI-level metadata.

Outputs
- snRNA-seq reference
   Nuclei × Genes

- Spatial expression
   ROIs × Genes

- ROI metadata
   ROIs × Covariates

These standardized inputs feed the cellular-composition module.

## Module 2 - Cellular Composition

Question: Which cell populations contribute to each spatial microenvironment?

Gene Sextant uses a two-stage reference-based deconvolution strategy:

1.- regression-derived cell-type expression signatures from annotated snRNA-seq;
2.- Bayesian spatial deconvolution of GeoMx WTA counts using a SpaceJam/Cell2Location-style hierarchical model implemented in Pyro.

The spatial model incorporates:

- snRNA-seq-derived expression signatures
- WTA counts
- negative-probe background
- nuclei-count-informed priors
- scan/slide structure

Posterior abundance estimates are extracted as both absolute and relative cell-type abundance.

Validation

Deconvolution is evaluated using complementary strategies including:

* held-out synthetic pseudobulk recovery
* marker-abundance concordance
* regression-signature identity
* leave-one-scan-out robustness
* compositional sensitivity analyses
* cross-method benchmarking

Outputs
-ROI × Cell type abundance
   P(cell type | ROI)

Relative abundance estimates provide the primary input for downstream spatial statistical modeling.

## Module 3 - Statistical Evidence

Question: Which cellular or molecular features change with disease and local pathology across spatial microenvironments?

Gene Sextant uses scan-aware mixed-effects models to evaluate features across arterial, capillary, and parenchymal contexts.

For relative cell-type abundance:
```text
Y_ROI ~ Beta(μ, φ)

logit(μ) = Biological effects + (1 | Scan_ID)
```

Models are fitted using glmmTMB, planned contrasts are estimated with emmeans, and multiple testing is controlled using Benjamini-Hochberg FDR.

*Biological contrasts*

Four complementary contrasts separate disease and local pathology effects:

1.-Disease effect
   AD/CAA vs Control among amyloid-negative ROIs.
2.-Amyloid effect
   Amyloid-positive vs amyloid-negative ROIs within AD/CAA.
3.-Maximum-pathology effect
   AD/CAA amyloid-positive vs Control amyloid-negative.
4.-Overall disease effect
   Pathology-weighted AD/CAA vs Control.

The same statistical architecture can be applied to different feature spaces, including:
* ROI × Cell abundances
* ROI × Genes
* ROI × Pathway activities
* ROI × Candidate biomarkers

Outputs
These analyses produce standardized statistical effect matrices:
Contrast × Feature × Region

Effect size • Direction • FDR

## Module 4 - Mechanistic Interpretation

Question: How can statistically supported spatial changes be organized into testable biological hypotheses?

Module 4 integrates evidence across multiple biological scales rather than interpreting individual significant features in isolation.

Current analyses include:

* cell-type abundance changes
* spatial gene-expression effects
* pathway activity associated with cellular composition
* lineage-associated molecular perturbations
* cell-type co-occurrence structure
* anatomical comparisons across arteries, capillaries, and parenchyma

This module connects statistical associations to candidate mechanisms while preserving anatomical and pathological context.

```text
Statistical effect matrices
          │
          ├── Cell changes
          ├── Gene changes
          ├── Pathway changes
          └── Co-occurrence structure
          │
          ▼
Spatially resolved,
testable mechanistic hypotheses
```

These analyses identify associations and generate mechanistic hypotheses; they do not by themselves establish causal relationships.

## Module 5 - Biomarker & Target Prioritization

Question: Can human spatial multi-omic evidence refine the interpretation of candidate biomarkers and therapeutic targets?

Module 5 applies the Gene Sextant evidence framework to externally defined candidate gene panels.

```text
Candidate panel
      │
      ▼
Reference attribution
      │
      ▼
Spatial perturbation
Arteries | Capillaries | Parenchyma
      │
      ▼
Lineage attribution
      │
      ▼
Fine-cell-type associations
      │
      ▼
Integrated evidence
      │
      ├───────────────┐
      ▼               ▼
Biomarker         Therapeutic
interpretation    target prioritization
```

Current panel applications:

The repository currently supports analyses of:

* candidate fluid biomarker panels
* AGORA therapeutic targets
* single- and repeatedly nominated AGORA targets
* matrisome/ECM biomarker panels

Panel-agnostic architecture

Module 5 uses a common analytical interface while preserving source-specific
metadata throughout the workflow.

A panel therefore contains a standardized core:
- gene
- panel
- category
- headline
- source

while retaining arbitrary panel-specific annotations such as:

```text
AGORA
├── total_nominations
├── nominating_teams
├── programs
└── pharos_class

Matrisome
├── MatrisomeDivision
└── MatrisomeCategory
```

These annotations remain attached to the candidate throughout reference attribution, spatial modeling, lineage analysis, fine-cell-type analysis, and integrated evidence reporting.

This allows new biomarker or therapeutic candidate collections to use the same statistical framework without discarding their original provenance.


## Repository Structure

```text
multiomic-spatial-integration/
├── R/
│   ├── signature_export_utils.R          # Export cell-type signatures from snRNA-seq references
│   ├── proportion_stats_utils.R          # Beta mixed models for inferred proportions
│   ├── contrast_utils.R                  # Disease, amyloid, and weighted contrasts
│   ├── enrichment_utils.R                # Marker gene extraction and enrichment helpers
│   ├── pathway_proportion_utils.R        # Correlate cell-type abundance with spatial expression and enrichment
│   ├── marker_concordance_utils.R        # Validate inferred abundances using independent marker-gene expression
│   ├── robustness_utils.R                # Leave-one-scan-out robustness analyses for spatial contrasts
│   ├── cooccurrence_utils.R              # CLR-based cell-type co-occurrence and spatial correlation analysis
│   └── plotting_utils.R                  # Heatmaps, dotplots, boxplots, and summary visualizations
│   #Note: Several viz_xxx.R files are still being validated and mapped to the workflow
│
├── python/
│   ├── geomx_anndata_utils.py            # Build AnnData objects from GeoMx WTA data
│   ├── regression_utils.py               # Cell2Location regression signature inference
│   ├── spacejam_pyro_utils.py            # Bayesian SpaceJam / Pyro model helpers
│   ├── abundance_extraction_utils.py     # Extract and normalize inferred abundances
│   ├── pseudobulk_utils.py               # Generate synthetic pseudobulk ROIs for model validation
│   ├── pseudobulk_validation_utils.py    # Evaluate deconvolution recovery on synthetic datasets
│   └── plotting_utils.py                 # Exploratory spatial plotting utilities
│
├── models/
│   └── LocationModelWTAMultiExperimentHierarchicalGeneLevel_Modified.py
│
├── notebooks/
│   ├── 01_wta_to_anndata.py              # Convert GeoMx WTA outputs into AnnData format
│   ├── 02_regression_signatures.py       # Learn cell-type gene signatures via regression (preferred)
│   ├── 03_spacejam_cell2location.py      # Run Bayesian spatial deconvolution (Cell2Location / SpaceJam)
│   ├── 04_extract_cell_proportions.py    # Extract and normalize inferred cell-type abundances
│   ├── 05_plotting_and_stats.py          # Exploratory visualization and generate analysis-ready tables
│   ├── 06_pseudobulk_validation.py       # Validate deconvolution accuracy using synthetic pseudobulk mixtures
│   ├── run_pseudobulk_validation_gpu.sh  # Submit GPU job for pseudobulk validation workflow
│   ├── run_regression_validation_gpu.sh  # Submit GPU job for Cell2Location regression model training
│   ├── run_spacejam_gpu.sh               # Submit GPU job for Bayesian SpaceJam spatial deconvolution
│   └──Examples/
│       ├── WTA to AnnData.ipynb
│       ├── Regression.ipynb
│       ├── Cell2Location + SpaceJam.ipynb
│       ├── Cell Proportions after SpaceJam.ipynb
│       ├── Cell2location_Plotting_and_Stats.ipynb
│       └── Stats_Integration.Rmd
│
├── resources/
│   ├── README.md                        #Documentation for external reference resources
│   └── download_gmt_resources.sh        #Download MSigDB and pathway GMT files
│
├── scripts/
│   ├── check_marker_factor_mapping.R          # Compare regression signatures against independent marker sets
│   ├── run_sn_reference_export.R              # Export annotated snRNA-seq reference
│   ├── run_spatial_stats.R                    # Run spatial statistical analyses and contrasts
│   ├── export_pathway_gene_list.R            # Export marker-based pathway inputs
│   ├── run_pathway_proportion_link.R          # Associate cell-type abundance with pathway activity
│   ├── run_marker_concordance_check.R         # Validate inferred abundances using marker concordance
│   ├── run_contrast_robustness_check.R        # Perform leave-one-scan-out robustness analyses
│   ├── run_all_signature_marker_validation.py # Comprehensive validation of all inferred cell-type signatures
│   └── run_all_figures.R                      # Generate publication-ready figures from analysis outputs
│
├── docs/
│   ├── workflow_overview.md
│   ├── geomx_anndata_structure.md
│   ├── bayesian_deconvolution_notes.md
│   ├── pyro_model_modifications.md
│   ├── statistical_modeling.md
│   └── Stats_updated.md
│
├── data/
│   └── nuc_count.csv    # ROI-level nuclei counts for abundance normalization and compositional analyses
│   #Note, Rest of Raw data will be made available upon publication
│
├── results/
│   ├── AD_CAA_cluster_markers.txt        # Seurat-derived marker table for the reference cell types
│   ├── all_marker_validation/            # Integrated marker-concordance and signature-identity QC
│   ├── cell_proportions/                 # ROI-level absolute and relative inferred abundances
│   ├── geomx_exports/                    # GeoMx expression and negative-probe exports for downstream R analyses
│   ├── marker_concordance/               # Real-data marker-abundance concordance results
│   ├── pathway_proportion_link/          # Gene rankings and pathway enrichment linked to cell-type abundance
│   ├── pseudobulk_validation/            # Held-out synthetic-mixture deconvolution validation
│   ├── regression_model/                 # Condition-specific regression signatures and factor-order checks
│   ├── spacejam/                         # SpaceJam model outputs, manifests, metadata, and training diagnostics
│   └── spatial_stats/                    # Mixed-effects contrasts, heterogeneity, co-occurrence, and robustness 
│
├── README.md
├── CHANGELOG
└── .gitignore
```
# Note: Files under 'notebooks'
These files represent the modular, production-ready pipeline for spatial transcriptomics integration.
Although organized under notebooks/, they are implemented as .py scripts for:
   * reproducibility
   * version control stability
   * compatibility with HPC and batch execution

Each script:
   * uses utils
   * mirrors the logical steps of the analysis pipeline
   * can be run independently or as part of the full workflow

### Example notebooks:
This folder contains the original interactive notebooks and R Markdown files used during method development.

**Since V0.02 the Notebooks have had some bugs identified and fixed directly to the Notebooks and the example files are illustrative of the intital versions.**

These files:
* provide step-by-step exploratory workflows
* include intermediate checks and visualizations
* reflect the iterative development process
  
Recommended usage
* Use the .py scripts in notebooks → reproducible pipeline execution
* Use the Examples/ folder for:
  
      - understanding the workflow
      - exploratory analysis
      - adapting to new datasets


## Module 1 specifics:
###  snRNA-seq Reference Preparation
Annotated snRNA-seq data are used as the cellular reference for spatial deconvolution.
Main steps include:

- Ambient RNA correction when appropriate
- QC filtering
- doublet removal
- SCTransform normalization
- Harmony integration
- clustering and manual annotation
- targeted subclustering of astrocytes, microglia, and vascular cells
- Pearson correlation-based merging of highly similar subclusters


The reference includes major brain cell populations such as astrocytes, microglia, oligodendrocytes, OPCs, excitatory/inhibitory neurons, and vascular-associated populations. 

This repository uses the preprocessing logic developed in `neuro-snRNAseq-tools` to generate the annotated snRNA-seq reference. The original project-specific script included SoupX correction, QC, DoubletFinder, SCTransform, Harmony integration, manual annotation, targeted subclustering, Pearson-correlation-guided merging, and expression aggregation for spatial integration. 

#### Dependency: neuro-snRNAseq-tools

This script depends on the `neuro-snRNAseq-tools` repository.

Clone it locally:

```bash
git clone https://github.com/echimalj/neuro-snRNAseq-tools.git
```

### GeoMx WTA AnnData Construction
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

## Module 2 specifics:
### Bayesian Spatial Deconvolution
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


### Abundance Extraction and Annotation
After model training, posterior spatial factors are extracted and converted into:

- absolute abundance estimates
- relative abundance estimates
- ROI × cell-type long-format tables


These outputs are merged with GeoMx metadata and factor-to-cell-type annotation mappings for downstream visualization and statistical testing. 


### Validation framework
Multiple orthogonal validation strategies are implemented to evaluate inferred cell-type abundances.

Current validation includes:

• pseudobulk recovery

• independent marker-gene concordance

• regression-signature identity

• leave-one-scan robustness

• compositional analysis

• ROI nuclei normalization
   Because inferred abundances are compositional, the framework includes dedicated analyses to distinguish
      • true biological depletion
      from
      • apparent increases caused by compositional redistribution.

      ROI nuclei counts can also be incorporated to normalize inferred abundance on a per-cell basis.

These analyses distinguish technical failures from biologically expected subtype overlap and provide confidence in downstream interpretation.

### Deconvolution benchmarking (WORK IN PROGRESS)
A modular benchmarking framework compares the Bayesian model against multiple deconvolution algorithms, including
   - SpatialDecon
   - MuSiC
   - Bisque
   - DWLS
   - BayesPrism   
   - SPOTlight
   - RCTD   
   - STdeconvolve

with standardized outputs for method agreement and reproducibility.


## Module 3 specifics:
### Statistical Modeling
Spatially inferred cell-type proportions are modeled using **beta mixed-effects models** in R.
The core model structure is:
``` R
rel_abundance ~ contrast_variable + (1 | Scan_ID)
```
using:
``` R
glmmTMB::beta_family(link = "logit")
```
This framework is appropriate because cell-type abundances are continuous proportions bounded between 0 and 1 and exhibit heteroscedasticity. A random intercept for Scan_ID accounts for repeated measurements at the ROI/sample level.

--
Multi-contrast framework

To disentangle disease and pathology effects, we implement four complementary contrasts:

### 1. Disease effect
 AD/CAA vs Control using amyloid-free ROIs only
→ isolates disease-driven changes independent of amyloid pathology

### 2. Amyloid effect (within AD)
 Amyloid vs AmyloidFree within AD/CAA
→ captures pathology-specific effects within diseased tissue

### 3. Max pathology effect
 AD/CAA Amyloid vs Control AmyloidFree
→ represents the strongest pathological contrast

### 4. Weighted overall AD effect 
 Pathology-weighted AD/CAA vs Control

This model accounts for heterogeneous amyloid burden across regions by defining:
```text
AD_overall = (1 − w) × AD_AmyloidFree + w × AD_Amyloid
```
where:
```text
w = proportion of amyloid-positive ROIs within AD/CAA (region-specific)
```
→ provides a biologically realistic estimate of the overall disease effect

### Effect size and inference
- Effect sizes are reported as log2 odds ratios
- Statistical contrasts are computed using emmeans
- Multiple testing correction is performed using Benjamini–Hochberg FDR, applied within region and contrast

### Design strengths

This modeling framework:
- Separates disease and pathology effects
- Accounts for spatial heterogeneity in amyloid burden
- Provides both controlled and extreme contrasts
- Enables biologically interpretable comparisons across microenvironments

For full details, see:
``` text
docs/statistical_modeling.md
```

## Module 4 specifics:
## Pathway association analysis
Rather than performing enrichment on marker genes, the framework correlates inferred cell-type abundance with spatial gene expression, ranks genes by abundance association, and performs GSEA/fGSEA on those ranked lists.

This identifies biological pathways associated with changes in spatial abundance.

Individual cellular effects were subsequently integrated with spatial gene expression, pathway activity, lineage associations, and cellular co-occurrence.

```text
Cells + Genes + Pathways + Networks
                ↓
Spatially resolved mechanistic hypotheses
```

## Module 5 specifics:
The current architecture generalizes candidate analysis across different panel types.

Instead of building separate workflows for each candidate collection, Gene Sextant preserves panel-specific metadata while applying a shared analytical pipeline.

```text
                    Gene Sextant
                         │
       ┌─────────────────┼─────────────────┐
       ▼                 ▼                 ▼
  Biomarkers           AGORA          Matrisome
       │                 │                 │
       └─────────────────┼─────────────────┘
                         ▼
                Shared evidence model
                         ▼
          Spatially contextualized results
```
This design separates the analytical framework from the candidate source, making the architecture extensible to future biomarker panels, therapeutic target collections, and other biologically defined gene sets.


---

## Outputs
The pipeline generates outputs for each major analysis stage, including:

```text
Regression model
├── *_inferred_signatures.csv
├── *_mean_cluster_expr.csv

Bayesian deconvolution
├── *_spot_factors_abs.pt
├── *_spot_factors_rel.pt
├── *_celltype_abundance_long.csv

Validation
├── all_marker_validation/
├── marker_concordance/
├── pseudobulk_validation/

Statistical analysis
├── spatial_stats/
├── pathway_proportion_link/
```
Most downstream analyses use the standardized abundance tables in
`results/cell_proportions/`, allowing independent statistical workflows
to be reproduced without rerunning Bayesian deconvolution.

## Example Use Case

This framework was developed for spatial multi-omic integration of
Alzheimer's disease (AD) and cerebral amyloid angiopathy (CAA) using
human postmortem brain tissue.

The workflow integrates:

- annotated single-nucleus RNA-seq reference atlases
- GeoMx Whole Transcriptome Atlas spatial transcriptomics
- Bayesian cell-type deconvolution
- regression-derived cell-type signatures
- ROI-level statistical modeling
- marker-based biological validation
- pathway enrichment linked to spatial cell-type abundance

In the current application, this enables compartment-resolved analysis of vascular, glial, and neuronal remodeling across parenchymal and vascular amyloid microenvironments.

Although developed for AD/CAA, the framework is readily adaptable to any study combining annotated sc/snRNA-seq references with GeoMx WTA or related spatial transcriptomic platforms.

### Notes

- The repository is organized as a modular computational framework rather than a single-click pipeline.
- Each analysis stage can be executed independently using the corresponding scripts.
- Validation modules are designed to distinguish technical artifacts from biologically meaningful subtype overlap.
- Sensitive raw data and patient metadata should never be committed.
- Large trained model objects are best managed outside Git or through Git LFS.
- The modified Bayesian SpaceJam implementation is intended for research use and should be validated for new datasets.

### Current Status

| Component                                    | Status        |
| -------------------------------------------- | ------------- |
| Module 1 — QC & preprocessing                | ✅ Implemented |
| Module 2 — Cellular composition              | ✅ Implemented |
| Module 3 — Statistical evidence              | ✅ Implemented |
| Module 4 — Mechanistic interpretation        | ✅ Implemented |
| Module 5 — Biomarker & target prioritization | ✅ Implemented |
| Pseudobulk validation                        | ✅ Implemented |
| Marker/signature validation                  | ✅ Implemented |
| Cross-method benchmarking                    | 🟡 Ongoing    |
| Publication/reporting figures                | 🟡 Ongoing    |


**Deconvolution methods**

| Method | Status |
|--------|--------|
| Cell2location (baseline) | ✅ Validated |
| SpatialDecon | ✅ Implemented & benchmarked |
| DWLS | ✅ Implemented & benchmarked |
| MuSiC | ✅ Implemented & benchmarked|
| Bisque | ✅ Implemented & benchmarked |
| BayesPrism | 🟡 Pending |
| SPOTlight | 🟡 Pending |
| RCTD | 🟡 Pending |
| STdeconvolve | 🟡 Pending |
| CIBERSORTx | 🟡 Input preparation implemented |


## Highlights

- Integration of human snRNA-seq references with GeoMx WTA spatial transcriptomics
- Regression-derived cell-type expression signatures
- Bayesian spatial deconvolution using a modified SpaceJam/Pyro framework
- GPU-accelerated stochastic variational inference
- Synthetic pseudobulk validation with known cellular composition
- Orthogonal marker, signature, robustness, and cross-method validation
- Scan-aware mixed-effects modeling across disease, pathology, and anatomical context
- Standardized `Contrast × Feature × Region` statistical effect matrices
- Spatial pathway and cell-type co-occurrence analysis
- Multi-scale integration of cellular, molecular, pathway, and spatial evidence
- Panel-agnostic biomarker and therapeutic target validation
- Preservation of source-specific candidate metadata throughout evidence integration
- Human spatial evidence for contextualizing established and emerging therapeutic candidates

### Authors:
Enrique Chimal
PhD Candidate – Medical Neuroscience - Indiana University School of Medicine

Subah Hussain
PhD Student - Genetics and Genomics - Baylor College of Medicine

Gene Sextant transforms spatial heterogeneity into interpretable cellular, statistical, mechanistic, and translational evidence.
