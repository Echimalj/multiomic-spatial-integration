# Review of Subah Branch
| File | Decision | Status |
|---|---|---|
| `notebooks/04b_extract_roi_total_counts.py` | Merge after edits | Approved and tested |
| `scripts/run_calibration_gate_check.R` | Merge | Approved and validated |
| `scripts/run_absolute_abundance_contrasts.R` | Merge | Approved and validated |
| `scripts/check_absolute_abundance_outputs.R` | Merge | Approved and validated |
| `R/viz_model_agreement.R` | Merge | Refactored, generalized, and validated |
| `scripts/run_model_agreement_figures.R` | Merge | Refactored, generalized, and validated |

## `notebooks/04b_extract_roi_total_counts.py`

**Purpose** 
Extracts per-ROI total counts, detected-gene counts, and DESeq2-style median-of-ratios size factors from the GeoMx AnnData object for the absolute-abundance and calibration workflows.

**Files reviewed together**
- `notebooks/04b_extract_roi_total_counts.py`
- New ROI-count utilities added to `python/abundance_extraction_utils.py`

**Issues identified**
- The original script was not portable because `python.abundance_extraction_utils` could not be imported when running the script directly.
- Input and output paths depended on the current working directory.
- The output directory was assumed to exist.
- The notebook depended on `compute_roi_total_counts()`, which was not yet present in the integration branch. (This was my bad, I was reviewing them appart haha ECJ )
- The original implementation did not explicitly protect against using transformed rather than raw counts.
- Size factors were not centered to a geometric mean of 1.

**Changes made**
- Added repository-root detection using `Path(__file__).resolve().parents[1]`.
- Added the repository root to `sys.path` before importing local Python utilities.
- Replaced relative paths with repository-root-relative input and output paths.
- Added input existence and AnnData dimension checks.
- Added output-table validation and automatic output-directory creation.
- Added QC summary statistics to the script log.
- Manually integrated the following utilities without overwriting existing branch changes:
  - `_build_roi_id_frame()`
  - `compute_median_of_ratios_size_factors()`
  - `compute_roi_total_counts()`
- Updated `compute_roi_total_counts()` to prefer `adata.layers["counts"]` when available.
- Added validation that the selected matrix contains finite, nonnegative, integer-like values.
- Normalized size factors so their geometric mean equals 1.
- Preserved the existing output column name `size_factor_mor` for downstream compatibility.

**Validation performed**
- Script executed successfully from the repository root.
- Log written to `results/Script_results/`.
- Output created at:
  - `results/cell_proportions/roi_total_counts.csv`
- AnnData shape: `190 × 18,676`.
- Count source: `adata.layers["counts"]`.
- Count matrix was integer-valued, finite, nonnegative, and had a minimum value of 1.
- Output contained 190 rows and 5 columns.
- No duplicated `ROI_ID` values.
- No missing values.
- No nonpositive total counts or size factors.
- All 18,676 genes contributed to the size-factor reference.
- Detected-gene count was 18,676 for every ROI, consistent with the GeoMx matrix having a minimum value of 1.

**Output summary before geometric-mean rescaling**
- Total counts: 18,905–169,461
- Median total counts: 46,656.5
- Size factors: 0.616–3.492
- Initial geometric mean of size factors: 0.9606

**Decision** 
**Merge after edits — approved.**

The script is now portable, reproducible, and appropriately validated. The scientific implementation is acceptable for this GeoMx count matrix, and the added safeguards reduce the risk of silently using transformed data.

## `scripts/run_calibration_gate_check.R`

**Purpose** 
Evaluates whether disease-associated differences in ROI cellularity require calibration before interpreting absolute-abundance models. Fits a calibration model, estimates the disease-associated scaling factor, and determines whether downstream analyses require correction before confirmatory interpretation.

**Issues identified**
- The original script depended on the working directory.
- Input files were not validated before model fitting.
- Calibration metadata were not propagated to downstream analyses.
- The decision summary could be overwritten later in the script.
- Logging and reproducibility information were limited.

**Changes made**
- Removed dependence on `setwd()` by using repository-relative paths.
- Added validation of required input files.
- Improved console logging and summary output.
- Saved the fitted calibration model and session information.
- Preserved calibration metadata including:
  - disease coefficient
  - standard error
  - multiplicative scale factor
  - inverse correction factor
  - calibration decision

**Validation performed**
- Calibration model executed successfully.
- Estimated disease coefficient: **0.2455** (log2 scale).
- Estimated multiplicative shift: **1.185×**.
- Estimated inverse correction factor: **0.8435**.
- Calibration decision: **proceed**.
- No correction required under the predefined calibration threshold.

**Decision** 
**Merge — approved.**

The calibration workflow is reproducible, well documented, and appropriately determines whether absolute-abundance analyses require calibration before confirmatory interpretation.

---

## `scripts/run_absolute_abundance_contrasts.R`

**Purpose** 
Runs the complete absolute-abundance differential analysis using both normalization strategies, propagates calibration metadata to all outputs, and generates combined summaries for downstream interpretation.

**Issues identified**
- The script depended on the working directory.
- Calibration metadata were not propagated to exported results.
- Returned contrast objects were not unpacked correctly.
- Helper functions were defined after first use.
- `glmmTMB::gaussian()` was used instead of `stats::gaussian()`.
- Silent model failures were difficult to diagnose.
- Vector/scalar logical evaluation inside `mutate()` produced incorrect confirmatory flags.

**Changes made**
- Removed dependence on `setwd()`.
- Added calibration-gate integration throughout the workflow.
- Added validation of joins and duplicate ROI identifiers.
- Correctly unpacked returned contrast objects.
- Moved helper functions before first use.
- Replaced `glmmTMB::gaussian()` with `stats::gaussian()`.
- Replaced silent `try()` calls with `tryCatch()` for improved diagnostics.
- Corrected vectorized logical evaluation for confirmatory eligibility.
- Added calibration metadata to all exported summaries.

**Validation performed**
- Successfully completed analyses using both:
  - `size_factor_mor`
  - `total_counts`
- Generated:
  - individual contrast tables
  - combined summaries
  - weighted summaries
  - model diagnostics
- Calibration metadata were correctly propagated to every exported result.
- No model-fitting failures occurred.

**Decision** 
**Merge — approved.**

The pipeline is reproducible, calibration-aware, and produces consistent absolute-abundance results across both normalization strategies.

---

## `scripts/check_absolute_abundance_outputs.R`

**Purpose** 
Performs quality control of the complete absolute-abundance workflow by validating calibration-gate outputs, model diagnostics, normalization robustness, and agreement with the proportion-based analyses.

**Issues identified**
- Effect-size column selection preferred the generic `estimate` column over explicitly named effect-size columns.
- The proportion-analysis summary contained an `estimate` column that was entirely missing for the Amyloid and Disease contrasts, causing incorrect `NA` agreement statistics.
- Effect-size selection depended only on the existence of a column rather than whether it contained usable values.

**Changes made**
- Implemented automatic selection of the first non-empty effect-size column.
- Prioritized explicitly named effect-size columns (`log2_OR`, `log2FC`, etc.) over generic `estimate` columns.
- Added diagnostic reporting of the selected effect-size column.
- Expanded QC reporting to include:
  - calibration-gate validation
  - fit diagnostics
  - normalization sensitivity
  - absolute-versus-proportion agreement
  - closure-artifact assessment
  - reproducibility metadata

**Validation performed**
- Overall QC status: **PASSED**.
- All calibration-gate decisions were **proceed**.
- All results were eligible for confirmatory interpretation.
- No flagged model-fitting problems.
- Excellent agreement between normalization strategies:

| Contrast | Spearman | Direction agreement |
|---|---:|---:|
| Amyloid | 0.990 | 96.4% |
| Disease | 0.966 | 99.3% |
| MaxPathology | 0.995 | 93.5% |
| Overall | 0.991 | 100% |

Absolute-versus-proportion agreement:

| Contrast | Spearman |
|---|---:|
| Amyloid | 0.294 |
| Disease | 0.534 |
| MaxPathology | 0.130 |
| Overall | 0.641 |

Closure-artifact assessment:
- Possible closure artifacts: **59**
- Strong closure-artifact candidates: **13**

All QC checklist items passed successfully.

The QC framework now validates the complete absolute-abundance workflow, correctly handles heterogeneous summary formats, and provides reproducible assessment of calibration, normalization robustness, and agreement with proportion-based analyses. The generalized model-agreement framework further confirms that the remaining differences between absolute-abundance and proportion-based analyses reflect methodological and biological differences rather than implementation artifacts, while preserving a reusable comparison infrastructure for future analytical models.

**Decision** 
**Merge — approved.**

## `R/viz_model_agreement.R` and `scripts/run_model_agreement_figures.R`

The QC framework now validates the complete absolute-abundance workflow, correctly handles heterogeneous summary formats, and provides reproducible assessment of calibration, normalization robustness, and agreement with proportion-based analyses. The remaining differences between absolute-abundance and proportion analyses appear to reflect methodological and biological differences rather than implementation errors.

The objective of this review was **not** to merge every file from the feature branch unchanged, but rather to integrate the validated scientific and computational improvements into the main analysis framework while preserving a coherent repository architecture.

Several files in the feature branch were developed specifically to compare analyses performed **with and without fibroblast exclusion**. During review, we concluded that fibroblast exclusion represents a distinct biological sensitivity analysis rather than the primary analytical workflow. Because excluding fibroblasts changes the biological question being asked, these workflows were intentionally **not** merged into the main pipeline.

Instead, we adopted a "daughter function" strategy.

Rather than modifying or extending the retired fibroblast-comparison functions directly, we created new generalized implementations that reuse the validated computational concepts while removing assumptions specific to fibroblast exclusion.

Specifically:

- `R/viz_model_comparison.R`
  → replaced by
  `R/viz_model_agreement.R`

- `scripts/run_model_comparison_figures.R`
  → replaced by
  `scripts/run_model_agreement_figures.R`

The new implementation generalizes model comparison to support arbitrary spatial-analysis pipelines, including:

- relative abundance
- absolute abundance (offset-adjusted)
- absolute abundance (raw library-size offset)
- future normalization or statistical models

without embedding assumptions about fibroblast exclusion.

The original fibroblast-exclusion scripts remain valuable as exploratory sensitivity analyses but were intentionally kept outside the primary workflow because they address a different biological question.

This approach allowed the repository to retain the validated statistical methodology while avoiding propagation of workflow-specific assumptions into the core analysis framework.

**Purpose**

Provides a generalized visualization framework for comparing matched spatial-transcriptomic effect estimates across multiple analytical pipelines.

The implementation replaces the retired fibroblast-comparison framework with a model-agnostic comparison system capable of evaluating agreement between relative-abundance, absolute-abundance, and future statistical models while preserving backward compatibility with existing visualization workflows.

**Files reviewed together**

- `R/viz_model_agreement.R`
- `scripts/run_model_agreement_figures.R`

**Retired parent functions**

- `R/viz_model_comparison.R`
- `scripts/run_model_comparison_figures.R`

**Related workflows intentionally not merged**

- `R/exclude_fibroblast_utils.R`
- `scripts/run_spatial_stats_no_fibroblast.R`

**Rationale**

The original visualization framework was designed specifically to compare analyses before and after fibroblast exclusion.

During review, we concluded that fibroblast exclusion represents a biologically distinct sensitivity analysis rather than a preprocessing step that should be embedded within the primary workflow.

Because excluding fibroblasts fundamentally changes the biological hypothesis being tested, directly extending the original implementation would have unnecessarily coupled future analyses to a fibroblast-signature comparison paradigm.

Instead, we retained the validated computational methodology while replacing the fibroblast-specific abstractions with a generalized model-agreement framework applicable to any pair of spatial-analysis models.

This strategy minimizes duplicated visualization code while allowing future analytical methods to reuse the same comparison infrastructure.

**Issues identified**

- Visualization logic assumed fibroblast exclusion as the comparison target.
- Comparison utilities only supported `log2_OR` effect sizes.
- Absolute-versus-relative comparison duplicated substantial plotting logic.
- Multiple supported effect-size columns produced ambiguous comparisons.
- QC summaries unnecessarily required effect-size detection despite using only adjusted p-values.
- Several tidyselect expressions used deprecated `.data` syntax.
- The comparison runner depended on fibroblast-specific configuration.

**Changes made**

- Replaced the fibroblast-specific comparison framework with a generalized model-agreement implementation.
- Added validation of required summary columns.
- Added duplicate join-key validation before model comparison.
- Generalized effect-size handling to support:
  - `log2_OR`
  - `log2_fold_change`
- Added explicit effect-column selection for summaries containing multiple supported effect-size statistics.
- Standardized comparison plots using:
  - symmetric axes
  - identity reference line
  - zero-reference lines
  - standardized significance categories
- Added facet-level agreement summaries reporting:
  - Pearson correlation
  - Spearman correlation
  - same-direction agreement
  - matched comparison counts
- Added descriptive QC summaries reporting significant-result frequency across analytical models.
- Removed deprecated tidyselect syntax.
- Implemented metadata-driven model configuration supporting arbitrary comparison pairs.
- Refactored `plot_absolute_vs_relative_fc_scatter()` into a lightweight backward-compatible wrapper that delegates all plotting to the generalized implementation.

**Validation performed**

Successfully generated:

- Descriptive significance-rate QC figure.
- Relative abundance vs absolute abundance (offset-adjusted).
- Relative abundance vs absolute abundance (raw library-size offset).
- Absolute abundance (offset-adjusted) vs absolute abundance (raw library-size offset).

For every comparison:

- 552 matched contrasts joined successfully.
- Zero unmatched rows.
- Correlation summaries generated for every analysis type.
- Figures written successfully without warnings or errors.

Backward compatibility was additionally validated by successfully regenerating:

- `absolute_vs_relative_fc_scatter`

for both absolute-abundance normalization strategies through the compatibility wrapper.

**Decision**

**Merge — approved.**

The visualization framework is now generalized, reusable, and substantially easier to maintain. The refactoring removes duplicated plotting logic, supports heterogeneous effect-size definitions, preserves backward compatibility with existing workflows, and provides a common infrastructure for future model-comparison analyses without coupling the primary pipeline to fibroblast-specific sensitivity analyses.


---

# Implemented Changes

The following improvements were incorporated into `main` following review of
Subah's branch.

## Pathway recurrence analysis

- Renamed the recurrence workflow for clarity:
  - `run_pathway_specificity_summary.R` →
    `run_pathway_recurrence_summary.R`
  - `run_pathway_sanity_check.R` →
    `export_pathway_gene_list.R`

- Replaced `compute_pathway_specificity()` with
  `summarize_pathway_recurrence()`.

- Recurrence is now computed from unique
  **pathway × cell-type × stratum** associations, preventing duplicate rows
  from inflating recurrence counts.

- Added comprehensive input validation and empty-result handling.

- Introduced `load_pathway_enrichment_results()` as a reusable loader for
  pathway enrichment outputs.

- Replaced the fixed recurrence threshold with a proportion-based threshold
  (`MAX_FREQUENCY_PROP`), allowing the cutoff to scale automatically with the
  number of tested cell-type–stratum combinations.

- Added pathway annotations:
  - `n_significant`
  - `recurrence_fraction`
  - `n_lineages`
  - `cell_lineages`
  - `broad_category_match`
  - `is_housekeeping`

- Added lineage-aware summaries to distinguish recurrence across related
  biological compartments without altering the recurrence statistic itself.

- Collapsed vascular-associated cell types
  (Endothelial, Pericyte, SMC, VLMC, Fibroblast) into a common
  `"Vascular"` lineage for improved interpretability.

- Added `plot_pathway_recurrence_dotplot()` with improved labeling,
  NES coloring, and recurrence summaries.

- Updated console logging to report recurrence and represented lineage groups.

---

## Co-occurrence network visualization

The context-specific co-occurrence network implemented in `main` supersedes
the exploratory network visualization in Subah's branch.

Major improvements include:

- Fixed node coordinates across biological contexts by computing a
  region-specific union graph and reusing its layout in every panel.

- Preserved isolated cell types by explicitly including all vertices in the
  network.

- Added lineage-based node coloring.

- Scaled node size according to union degree.

- Encoded edge direction and strength using both color and width.

- Normalized network coordinates to allow direct visual comparison between
  panels.

- Added descriptive figure annotations documenting the exploratory nature of
  the visualization and the limited number of independent scans per context.

The original exploratory network remains useful for quick visualization of
significant pairwise correlations, but the context-specific implementation
provides substantially improved biological interpretability.
