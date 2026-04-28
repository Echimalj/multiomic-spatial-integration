# Pyro Model Modifications

This repository uses a modified version of the Cell2Location model:

LocationModelWTAMultiExperimentHierarchicalGeneLevel_Modified.py

---

## Motivation

The original model was adapted to better handle:

- GeoMx WTA data structure
- Multi-experiment datasets
- Gene-level hierarchical variation

---

## Key Modifications

### 1. Multi-experiment handling
- Improved handling of batch-specific variation

### 2. Gene-level hierarchy
- Allows gene-specific variance modeling

### 3. Improved convergence stability
- Adjusted priors and optimization

---

## Impact

- Improved robustness across ROIs
- Better biological interpretability
- Reduced overfitting in low-count regions

---

## Notes

These modifications are experimental and tailored for GeoMx WTA data.
