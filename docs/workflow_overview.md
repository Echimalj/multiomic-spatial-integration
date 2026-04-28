# Workflow Overview

This repository implements a full multi-omic integration pipeline combining:

- Single-nucleus RNA-seq (snRNA-seq)
- Spatial transcriptomics (GeoMx WTA)
- Bayesian deconvolution (Cell2Location / SpaceJam)
- Statistical modeling (beta mixed-effects)

---

## Pipeline Summary

```text
snRNA-seq (Seurat)
        ↓
Reference annotation + subclustering
        ↓
Signature extraction * See notes in Notebook *
        ↓
GeoMx WTA → AnnData
        ↓
Regression-based signature learning
        ↓
Bayesian spatial deconvolution
        ↓
Cell-type abundance estimation
        ↓
Exploratory visualization (Python)
        ↓
Formal statistical modeling (R)
```

---

## Key Components
1. snRNA-seq Reference
  - Built using neuro-snRNAseq-tools
  - Includes:
      * SCTransform normalization
      * Harmony integration
      * Manual annotation
      * Subclustering + merging

Output:
- Annotated Seurat object
- Aggregated signatures (optional)

  
2. Spatial Data Processing
- GeoMx WTA → AnnData

- Includes:
    * ROI metadata
    * Negative probe matrices
    * Gene-level counts
      
3. Signature Learning (Preferred)
- Regression-based (Notebook 02)
- Accounts for:
    * Technical noise
    * Compositional bias

- Alternative:
   *Aggregation (Seurat)

4. Bayesian Deconvolution
- Cell2Location / SpaceJam
  
Outputs:
- Cell-type abundance per ROI
- Posterior distributions
  
5. Statistical Modeling
- Beta mixed-effects models (glmmTMB)
- Controls:
    * ROI-level correlation
    * Disease status
    * Pathology
      
Key Outputs:
- Cell-type abundance tables
- Spatial cell-type maps
- Statistical results
- Pathway validation inputs
