# Bayesian Deconvolution Notes

This pipeline uses Cell2Location / SpaceJam to infer cell-type composition in spatial transcriptomics data.

---

## Concept

Spatial transcriptomics measures:
- Mixed gene expression from multiple cell types

Goal:
- Infer contribution of each cell type per ROI

---

## Model Overview

Observed counts are modeled as:

Gene expression = Σ(cell type signatures × abundance) + noise

---

## Key Advantages

- Accounts for:
  - Technical noise
  - Sampling variability
  - Gene-specific dispersion

- Produces:
  - Posterior distributions (uncertainty-aware)

---

## Inputs

- AnnData (GeoMx WTA)
- Reference signatures (snRNA-seq)

---

## Outputs

- Cell-type abundance per ROI
- Relative proportions
- Uncertainty estimates

---

## Preferred Approach

Regression-based signature learning is used because it:

- Reduces bias from compositional effects
- Improves robustness across datasets
