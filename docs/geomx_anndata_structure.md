# GeoMx AnnData Structure

This document describes the structure of the AnnData object generated from GeoMx WTA data.

---

## Main Components

### X matrix
- Gene expression matrix
- Shape: (ROIs × genes)

---

### obs (ROI metadata)

Includes:
- `ROI`
- `disease_status`
- `pathology`
- `region`
- `Scan_ID`

---

### var (genes)

Includes:
- Gene IDs
- Gene symbols
- Annotation fields

---

### obsm

#### `negProbes`
- Negative control probe counts
- Used for background estimation

---

## Notes

- ROIs are the fundamental spatial unit
- Data is not single-cell — requires deconvolution
- Metadata alignment is critical for downstream modeling
