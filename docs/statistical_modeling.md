# Statistical Modeling

Formal statistical testing is performed using beta mixed-effects models.

---

## Why Beta Models?

Cell-type proportions:
- Continuous
- Bounded between 0 and 1

→ Beta distribution is appropriate

---

## Model

```r
rel_abundance ~ condition + (1 | ROI)
```

Implemented with 
```r
glmmTMB::beta_family(link = "logit")
```

---
Fixed Effects
- Disease status
- Pathology
- Region (optional)

  
Random Effects
- ROI / sample ID

Accounts for:
- Repeated measurements
- Biological variability


Post-hoc Testing
- Performed using ```emmeans```
- Pairwise contrasts:
     * WT vs AD
     * AD vs AD+CAA
     * etc.
  
Multiple Testing
- Benjamini-Hochberg correction

---  
Notes
Exploratory statistics in Python (Notebook 05) are not used for final inference.
