# ============================================================
# Absolute cell-type abundance contrasts (total-count offset)
#
# Mirrors run_spatial_stats_no_fibroblast.R's shape: same 4-contrast
# framework as run_spatial_stats.R, reused here against log(abs_abundance)
# with a per-ROI total-count offset instead of rel_abundance.
#
# Requires scripts/run_calibration_gate_check.R to have already run --
# every output row is stamped with the gate's gate_status, and
# Disease/Overall/MaxPathology should not be reported as confirmatory
# unless gate_status == "proceed".

# The gate is currently diagnostic only. No calibration correction is
# automatically applied when gate_status == "proceed_with_correction".
# Such outputs remain non-confirmatory unless a separate, explicit
# correction sensitivity analysis is run.
#
# Runs twice:
#   - offset_col = "size_factor_mor"  (primary, DESeq2-style)  -> results/spatial_stats_absolute/
#   - offset_col = "total_counts"     (raw-sum sensitivity)    -> results/spatial_stats_absolute_rawlibsize/
# Divergence between the two is itself diagnostic, not discarded.
# ============================================================

library(glmmTMB)
library(emmeans)
library(dplyr)
library(readr)

source("R/contrast_utils.R")
source("R/absolute_abundance_utils.R")

gate_summary_file <- "results/spatial_stats_absolute/calibration_gate_summary.csv"

if (!file.exists(gate_summary_file)) {
  stop(
    "Calibration gate summary not found at ",
    gate_summary_file,
    ". Run scripts/run_calibration_gate_check.R first.",
    call. = FALSE
  )
}

# ============================================================
# Parse the gate decision file into a per-contrast-type lookup
# ============================================================

gate_summary <- readr::read_csv(
  gate_summary_file,
  show_col_types = FALSE
)

if (nrow(gate_summary) != 1) {
  stop(
    "Calibration gate summary must contain exactly one row; found ",
    nrow(gate_summary),
    ".",
    call. = FALSE
  )
}

required_gate_columns <- c(
  "decision",
  "disease_coefficient_log2",
  "disease_se_log2",
  "biological_threshold_log2",
  "region_pattern_correlation",
  "correction_applied"
)

missing_gate_columns <- setdiff(
  required_gate_columns,
  colnames(gate_summary)
)

if (length(missing_gate_columns) > 0) {
  stop(
    "Calibration gate summary is missing required column(s): ",
    paste(missing_gate_columns, collapse = ", "),
    call. = FALSE
  )
}

allowed_gate_decisions <- c(
  "proceed",
  "proceed_with_correction",
  "do_not_proceed"
)

cross_condition_gate_status <- gate_summary$decision[[1]]

if (
  length(cross_condition_gate_status) != 1 ||
    is.na(cross_condition_gate_status) ||
    !cross_condition_gate_status %in% allowed_gate_decisions
) {
  stop(
    "Invalid calibration-gate decision: ",
    cross_condition_gate_status,
    call. = FALSE
  )
}

if (isTRUE(gate_summary$correction_applied[[1]])) {
  stop(
    "The gate summary claims that a correction was applied, but this ",
    "contrast script does not currently implement or verify corrected ",
    "abundance values.",
    call. = FALSE
  )
}

gate_status_lookup <- c(
  Amyloid = "proceed",
  Disease = cross_condition_gate_status,
  Overall = cross_condition_gate_status,
  MaxPathology = cross_condition_gate_status
)

message("Calibration gate decisions loaded:")
print(gate_status_lookup)

message(
  "Calibration correction applied: FALSE"
)

if (cross_condition_gate_status == "proceed_with_correction") {
  warning(
    "Calibration gate returned 'proceed_with_correction'. ",
    "The current script will fit and save uncorrected cross-condition ",
    "contrasts, but they will be marked confirmatory_eligible = FALSE.",
    call. = FALSE
  )
}

if (cross_condition_gate_status == "do_not_proceed") {
  warning(
    "Calibration gate returned 'do_not_proceed'. ",
    "Cross-condition contrasts will be generated for diagnostic purposes ",
    "only and marked confirmatory_eligible = FALSE.",
    call. = FALSE
  )
}

# ============================================================
# Load abundance + total-count data (shared across both offset runs)
# ============================================================

abs_long_file <- "results/cell_proportions/roi_celltype_abundance_long_abs.csv"
total_counts_file <- "results/cell_proportions/roi_total_counts.csv"

required_input_files <- c(
  abs_long_file,
  total_counts_file
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    "Absolute-abundance input file(s) not found:\n",
    paste0("  - ", missing_input_files, collapse = "\n"),
    "\nRun this script from the repository root.",
    call. = FALSE
  )
}

abs_df_raw <- readr::read_csv(
  abs_long_file,
  show_col_types = FALSE
)

total_counts_df <- readr::read_csv(
  total_counts_file,
  show_col_types = FALSE
)

abs_only <- setdiff(abs_df_raw$ROI_ID, total_counts_df$ROI_ID)
counts_only <- setdiff(total_counts_df$ROI_ID, abs_df_raw$ROI_ID)

if (length(abs_only) > 0 || length(counts_only) > 0) {
  stop(
    "Incomplete join between abs_abundance ROIs and total_counts ROIs.\n",
    "In abs_abundance but not total_counts (", length(abs_only), "): ",
    paste(utils::head(abs_only, 10), collapse = ", "), "\n",
    "In total_counts but not abs_abundance (", length(counts_only), "): ",
    paste(utils::head(counts_only, 10), collapse = ", "),
    call. = FALSE
  )
}

duplicated_count_rois <- unique(
  total_counts_df$ROI_ID[
    duplicated(total_counts_df$ROI_ID)
  ]
)

if (length(duplicated_count_rois) > 0) {
  stop(
    "Duplicate ROI rows detected in total_counts_df: ",
    paste(
      utils::head(duplicated_count_rois, 10),
      collapse = ", "
    ),
    if (length(duplicated_count_rois) > 10) ", ..." else "",
    call. = FALSE
  )
}

merged_df <- dplyr::inner_join(
  abs_df_raw,
  total_counts_df,
  by = "ROI_ID",
  relationship = "many-to-one"
)

# ============================================================
# Helper: run all four contrasts for one offset choice, stamp gate_status,
# write the same file set as spatial_stats_no_fibroblast (per-contrast CSVs
# + region weights + fit status + combined summary).
# ============================================================

# ============================================================
# Calibration-gate stamping helpers
# ============================================================

stamp_gate <- function(contrast_df, type_label) {
  gate_status <- unname(gate_status_lookup[[type_label]])

  if (
    length(gate_status) != 1 ||
      is.na(gate_status)
  ) {
    stop(
      "No valid calibration-gate status found for contrast type '",
      type_label,
      "'.",
      call. = FALSE
    )
  }

  is_amyloid <- identical(type_label, "Amyloid")

  contrast_df |>
    dplyr::mutate(
      gate_status = gate_status,

      confirmatory_eligible =
        if (is_amyloid) {
          TRUE
        } else {
          gate_status == "proceed"
        },

      calibration_correction_applied = FALSE,

      calibration_disease_coefficient_log2 =
        if (is_amyloid) {
          NA_real_
        } else {
          gate_summary$disease_coefficient_log2[[1]]
        },

      calibration_disease_se_log2 =
        if (is_amyloid) {
          NA_real_
        } else {
          gate_summary$disease_se_log2[[1]]
        },

      calibration_region_pattern_correlation =
        if (is_amyloid) {
          NA_real_
        } else {
          gate_summary$region_pattern_correlation[[1]]
        }
    )
}


stamp_formatted_gate <- function(summary_df, source_df) {
  gate_columns <- c(
    "gate_status",
    "confirmatory_eligible",
    "calibration_correction_applied",
    "calibration_disease_coefficient_log2",
    "calibration_disease_se_log2",
    "calibration_region_pattern_correlation"
  )

  missing_gate_columns <- setdiff(
    gate_columns,
    colnames(source_df)
  )

  if (length(missing_gate_columns) > 0) {
    stop(
      "stamp_formatted_gate(): source dataframe is missing gate columns: ",
      paste(missing_gate_columns, collapse = ", "),
      call. = FALSE
    )
  }

  for (column_name in gate_columns) {
    summary_df[[column_name]] <- rep(
      source_df[[column_name]][[1]],
      nrow(summary_df)
    )
  }

  summary_df
}


# ============================================================
# Run and write absolute-abundance pipeline
# ============================================================

run_and_write_absolute_pipeline <- function(offset_col, output_dir) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  df <- prepare_absolute_abundance_data(
    df = merged_df,
    offset_col = offset_col
  )

  message(
    "Running absolute-abundance models with offset column: ",
    offset_col
  )

  # ----------------------------------------------------------
  # Run models
  # ----------------------------------------------------------

  amyloid_res <- run_amyloid_effect_absolute(df)
  disease_res <- run_disease_effect_absolute(df)
  overall_res <- run_weighted_overall_effect_absolute(df)
  maxpath_res <- run_max_pathology_effect_absolute(df)

  # ----------------------------------------------------------
  # Stamp calibration-gate metadata onto raw contrasts
  # ----------------------------------------------------------

  amyloid_contrasts <- stamp_gate(
    amyloid_res$contrasts,
    "Amyloid"
  )

  disease_contrasts <- stamp_gate(
    disease_res$contrasts,
    "Disease"
  )

  overall_contrasts <- stamp_gate(
    overall_res$contrasts,
    "Overall"
  )

  maxpath_contrasts <- stamp_gate(
    maxpath_res$contrasts,
    "MaxPathology"
  )

  # ----------------------------------------------------------
  # Write raw contrast outputs
  # ----------------------------------------------------------

  readr::write_csv(
    amyloid_contrasts,
    file.path(
      output_dir,
      "amyloid_effect_contrasts.csv"
    )
  )

  readr::write_csv(
    disease_contrasts,
    file.path(
      output_dir,
      "disease_effect_contrasts.csv"
    )
  )

  readr::write_csv(
    overall_contrasts,
    file.path(
      output_dir,
      "overall_effect_contrasts.csv"
    )
  )

  readr::write_csv(
    maxpath_contrasts,
    file.path(
      output_dir,
      "max_pathology_effect_contrasts.csv"
    )
  )

  # ----------------------------------------------------------
  # Write estimated means
  # ----------------------------------------------------------

  readr::write_csv(
    amyloid_res$means,
    file.path(
      output_dir,
      "amyloid_effect_means.csv"
    )
  )

  readr::write_csv(
    disease_res$means,
    file.path(
      output_dir,
      "disease_effect_means.csv"
    )
  )

  readr::write_csv(
    overall_res$means,
    file.path(
      output_dir,
      "overall_effect_means.csv"
    )
  )

  readr::write_csv(
    maxpath_res$means,
    file.path(
      output_dir,
      "max_pathology_effect_means.csv"
    )
  )

  # ----------------------------------------------------------
  # Write diagnostic and auxiliary outputs
  # ----------------------------------------------------------

  readr::write_csv(
    overall_res$weights,
    file.path(
      output_dir,
      "overall_effect_region_weights.csv"
    )
  )

  readr::write_csv(
    amyloid_res$fit_status,
    file.path(
      output_dir,
      "amyloid_effect_fit_status.csv"
    )
  )

  readr::write_csv(
    maxpath_res$fit_status,
    file.path(
      output_dir,
      "max_pathology_fit_status.csv"
    )
  )

  # ----------------------------------------------------------
  # Format standardized summaries
  # ----------------------------------------------------------

  amyloid_df <- format_absolute_contrast_summary(
    amyloid_contrasts,
    contrast_type = "Amyloid"
  )

  disease_df <- format_absolute_contrast_summary(
    disease_contrasts,
    contrast_type = "Disease"
  )

  overall_df <- format_absolute_contrast_summary(
    overall_contrasts |>
      dplyr::filter(
        .data$contrast == "AD_overall_vs_Control"
      ),
    contrast_type = "Overall"
  )

  maxpath_df <- format_absolute_contrast_summary(
    maxpath_contrasts,
    contrast_type = "MaxPathology"
  )

  # ----------------------------------------------------------
  # Restore calibration-gate columns after summary formatting
  # ----------------------------------------------------------

  amyloid_df <- stamp_formatted_gate(
    amyloid_df,
    amyloid_contrasts
  )

  disease_df <- stamp_formatted_gate(
    disease_df,
    disease_contrasts
  )

  overall_df <- stamp_formatted_gate(
    overall_df,
    overall_contrasts
  )

  maxpath_df <- stamp_formatted_gate(
    maxpath_df,
    maxpath_contrasts
  )

  combined_contrasts <- dplyr::bind_rows(
    amyloid_df,
    disease_df,
    overall_df,
    maxpath_df
  )

  readr::write_csv(
    combined_contrasts,
    file.path(
      output_dir,
      "combined_spatial_contrast_summary.csv"
    )
  )

  message(
    "Absolute-abundance pipeline completed for offset_col = '",
    offset_col,
    "'. Outputs written to: ",
    output_dir
  )

  invisible(
    list(
      amyloid = amyloid_res,
      disease = disease_res,
      overall = overall_res,
      max_pathology = maxpath_res,
      combined_summary = combined_contrasts
    )
  )
}
# ============================================================
# Primary: median-of-ratios size factor
# ============================================================

run_and_write_absolute_pipeline(
  offset_col = "size_factor_mor",
  output_dir = "results/spatial_stats_absolute"

)

# ============================================================
# Sensitivity variant: raw total-count offset (mirrors the Pyro model's
# internal l_r exactly, but is dominated by highly-expressed genes --
# divergence from the median-of-ratios run is itself diagnostic)
# ============================================================

run_and_write_absolute_pipeline(
  offset_col = "total_counts",
  output_dir = "results/spatial_stats_absolute_rawlibsize"
)

message("Absolute-abundance contrasts completed successfully.")
