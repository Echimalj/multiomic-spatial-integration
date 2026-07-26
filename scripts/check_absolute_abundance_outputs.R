#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

options(dplyr.summarise.inform = FALSE)

# ============================================================
# Paths
# ============================================================

project_root <- "/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/multiomic-spatial-integration"
results_root <- file.path(project_root, "results")

primary_dir <- file.path(results_root, "spatial_stats_absolute")
sensitivity_dir <- file.path(results_root, "spatial_stats_absolute_rawlibsize")
qc_dir <- file.path(results_root, "absolute_abundance_qc")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# Add another candidate here if your proportion results live elsewhere.
proportion_candidates <- c(
  file.path(results_root, "spatial_stats", "combined_spatial_contrast_summary.csv"),
  file.path(results_root, "spatial_stats_proportion", "combined_spatial_contrast_summary.csv"),
  file.path(results_root, "combined_spatial_contrast_summary.csv")
)

# ============================================================
# Helpers
# ============================================================

write_qc <- function(x, filename) {
  readr::write_csv(x, file.path(qc_dir, filename), na = "")
}

read_checked <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path, call. = FALSE)
  }

  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

first_existing <- function(paths) {
  found <- paths[file.exists(paths)]
  if (length(found) == 0) NA_character_ else found[[1]]
}

find_col <- function(df, candidates, required = TRUE, label = NULL) {
  hit <- candidates[candidates %in% names(df)]

  if (length(hit) > 0) {
    return(hit[[1]])
  }

  if (required) {
    stop(
      "Could not find ", ifelse(is.null(label), "required column", label),
      ". Tried: ", paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }

  NA_character_
}

normalize_logical <- function(x) {
  if (is.logical(x)) return(x)

  x <- tolower(trimws(as.character(x)))

  dplyr::case_when(
    x %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)

  if (
    sum(keep) < 3 ||
      dplyr::n_distinct(x[keep]) < 2 ||
      dplyr::n_distinct(y[keep]) < 2
  ) {
    return(NA_real_)
  }

  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

standardize_summary <- function(df, label) {
  type_col <- find_col(
    df,
    c("type", "contrast_type", "analysis_type"),
    label = paste(label, "type")
  )

  region_col <- find_col(
    df,
    c("region", "Region"),
    label = paste(label, "region")
  )

  celltype_col <- find_col(
    df,
    c("celltype", "cell_type", "CellType", "cellType"),
    label = paste(label, "cell type")
  )

  contrast_col <- find_col(
    df,
    c("contrast", "comparison", "comparison_name"),
    required = FALSE
  )

  effect_candidates <- c(
    "log2_fold_change",
    "log2FC",
    "log2_fc",
    "log2_OR",
    "log2OR",
    "estimate"
  )

  available_effect_cols <- effect_candidates[
    effect_candidates %in% names(df)
  ]

  effect_col <- NULL

  for (candidate in available_effect_cols) {
    values <- suppressWarnings(as.numeric(df[[candidate]]))

    if (any(is.finite(values))) {
      effect_col <- candidate
      break
    }
  }

  if (is.null(effect_col)) {
    stop(
      sprintf(
        "No usable effect-size column found for %s. Checked: %s",
        label,
        paste(effect_candidates, collapse = ", ")
      )
    )
  }

  message(
    sprintf(
      "[%s] Using effect-size column: %s",
      label,
     effect_col
    )
  )

  p_col <- find_col(
    df,
    c("p.value", "p_value", "pvalue", "p"),
    required = FALSE
  )

  padj_col <- find_col(
    df,
    c("p_adj", "padj", "adjusted_p_value", "adj_p_value", "FDR", "q_value"),
    label = paste(label, "adjusted p-value")
  )

  df |>
    mutate(
      type_std = as.character(.data[[type_col]]),
      region_std = as.character(.data[[region_col]]),
      celltype_std = as.character(.data[[celltype_col]]),
      contrast_std = if (!is.na(contrast_col)) {
        as.character(.data[[contrast_col]])
      } else {
        as.character(.data[[type_col]])
      },
      effect_std = suppressWarnings(as.numeric(.data[[effect_col]])),
      p_value_std = if (!is.na(p_col)) {
        suppressWarnings(as.numeric(.data[[p_col]]))
      } else {
        NA_real_
      },
      p_adj_std = suppressWarnings(as.numeric(.data[[padj_col]]))
    )
}

# ============================================================
# Required files
# ============================================================

required_files <- c(
  "amyloid_effect_contrasts.csv",
  "amyloid_effect_means.csv",
  "disease_effect_contrasts.csv",
  "disease_effect_means.csv",
  "overall_effect_contrasts.csv",
  "overall_effect_means.csv",
  "overall_effect_region_weights.csv",
  "max_pathology_effect_contrasts.csv",
  "max_pathology_effect_means.csv",
  "max_pathology_fit_status.csv",
  "combined_spatial_contrast_summary.csv"
)

optional_files <- c(
  "amyloid_effect_fit_status.csv",
  "disease_effect_fit_status.csv",
  "overall_effect_fit_status.csv"
)

check_directory <- function(directory, normalization) {
  tibble(
    normalization = normalization,
    filename = c(required_files, optional_files),
    required = c(
      rep(TRUE, length(required_files)),
      rep(FALSE, length(optional_files))
    ),
    path = file.path(directory, c(required_files, optional_files))
  ) |>
    mutate(exists = file.exists(path))
}

file_inventory <- bind_rows(
  check_directory(primary_dir, "size_factor_mor"),
  check_directory(sensitivity_dir, "total_counts")
)

write_qc(file_inventory, "file_inventory.csv")

missing_required <- file_inventory |>
  filter(required, !exists)

if (nrow(missing_required) > 0) {
  stop(
    "Required files are missing:\n",
    paste(missing_required$path, collapse = "\n"),
    call. = FALSE
  )
}

# ============================================================
# Read outputs
# ============================================================

read_outputs <- function(directory) {
  list(
    amyloid = read_checked(file.path(directory, "amyloid_effect_contrasts.csv")),
    disease = read_checked(file.path(directory, "disease_effect_contrasts.csv")),
    overall = read_checked(file.path(directory, "overall_effect_contrasts.csv")),
    maxpath = read_checked(file.path(directory, "max_pathology_effect_contrasts.csv")),
    combined = read_checked(file.path(directory, "combined_spatial_contrast_summary.csv")),
    amyloid_means = read_checked(file.path(directory, "amyloid_effect_means.csv")),
    disease_means = read_checked(file.path(directory, "disease_effect_means.csv")),
    overall_means = read_checked(file.path(directory, "overall_effect_means.csv")),
    maxpath_means = read_checked(file.path(directory, "max_pathology_effect_means.csv")),
    overall_weights = read_checked(file.path(directory, "overall_effect_region_weights.csv")),
    maxpath_fit = read_checked(file.path(directory, "max_pathology_fit_status.csv")),
    amyloid_fit = {
      p <- file.path(directory, "amyloid_effect_fit_status.csv")
      if (file.exists(p)) read_checked(p) else NULL
    },
    disease_fit = {
      p <- file.path(directory, "disease_effect_fit_status.csv")
      if (file.exists(p)) read_checked(p) else NULL
    },
    overall_fit = {
      p <- file.path(directory, "overall_effect_fit_status.csv")
      if (file.exists(p)) read_checked(p) else NULL
    }
  )
}

primary <- read_outputs(primary_dir)
sensitivity <- read_outputs(sensitivity_dir)

# ============================================================
# Row counts
# ============================================================

row_counts <- bind_rows(
  tibble(
    normalization = "size_factor_mor",
    table = names(primary),
    n_rows = map_int(primary, ~ if (is.null(.x)) 0L else nrow(.x)),
    n_columns = map_int(primary, ~ if (is.null(.x)) 0L else ncol(.x))
  ),
  tibble(
    normalization = "total_counts",
    table = names(sensitivity),
    n_rows = map_int(sensitivity, ~ if (is.null(.x)) 0L else nrow(.x)),
    n_columns = map_int(sensitivity, ~ if (is.null(.x)) 0L else ncol(.x))
  )
)

write_qc(row_counts, "row_counts.csv")

# ============================================================
# Calibration gate
# ============================================================

gate_columns <- c(
  "gate_status",
  "confirmatory_eligible",
  "calibration_correction_applied",
  "calibration_disease_coefficient_log2",
  "calibration_disease_se_log2",
  "calibration_region_pattern_correlation"
)

gate_column_inventory <- bind_rows(
  tibble(
    normalization = "size_factor_mor",
    column = gate_columns,
    present = gate_columns %in% names(primary$combined)
  ),
  tibble(
    normalization = "total_counts",
    column = gate_columns,
    present = gate_columns %in% names(sensitivity$combined)
  )
)

write_qc(gate_column_inventory, "gate_column_inventory.csv")

summarize_gate <- function(df, normalization) {
  if (!all(gate_columns %in% names(df))) {
    return(tibble(
      normalization = normalization,
      type = NA_character_,
      n = nrow(df),
      gate_statuses = NA_character_,
      all_confirmatory_eligible = NA,
      correction_values = NA_character_
    ))
  }

  type_col <- find_col(df, c("type", "contrast_type", "analysis_type"))

  df |>
    mutate(
      type = as.character(.data[[type_col]]),
      confirmatory_eligible = normalize_logical(confirmatory_eligible),
      calibration_correction_applied = normalize_logical(
        calibration_correction_applied
      )
    ) |>
    group_by(type) |>
    summarise(
      normalization = normalization,
      n = n(),
      gate_statuses = paste(sort(unique(gate_status)), collapse = "; "),
      all_confirmatory_eligible = all(confirmatory_eligible, na.rm = TRUE),
      correction_values = paste(
        sort(unique(as.character(calibration_correction_applied))),
        collapse = "; "
      ),
      disease_coefficient_values = paste(
        unique(round(calibration_disease_coefficient_log2, 6)),
        collapse = "; "
      ),
      disease_se_values = paste(
        unique(round(calibration_disease_se_log2, 6)),
        collapse = "; "
      ),
      region_correlation_values = paste(
        unique(round(calibration_region_pattern_correlation, 6)),
        collapse = "; "
      ),
      .groups = "drop"
    ) |>
    select(normalization, everything())
}

gate_summary <- bind_rows(
  summarize_gate(primary$combined, "size_factor_mor"),
  summarize_gate(sensitivity$combined, "total_counts")
)

write_qc(gate_summary, "gate_summary.csv")

# ============================================================
# Fit diagnostics
# ============================================================

summarize_fit <- function(df, analysis, normalization) {
  if (is.null(df)) {
    return(tibble(
      normalization = normalization,
      analysis = analysis,
      status = "file_not_available",
      n_combinations = 0L
    ))
  }

  status_col <- find_col(
    df,
    c("status", "fit_status", "model_status", "result_status"),
    required = FALSE
  )

  if (!is.na(status_col)) {
    return(
      df |>
        mutate(status = as.character(.data[[status_col]])) |>
        count(status, name = "n_combinations") |>
        mutate(
          normalization = normalization,
          analysis = analysis
        ) |>
        select(normalization, analysis, status, n_combinations)
    )
  }

  fit_ok_col <- find_col(
    df,
    c("fit_ok", "success", "model_ok"),
    required = FALSE
  )

  if (!is.na(fit_ok_col)) {
    return(
      df |>
        mutate(
          fit_ok_std = normalize_logical(.data[[fit_ok_col]]),
          status = case_when(
            fit_ok_std ~ "success",
            is.na(fit_ok_std) ~ "unknown",
            TRUE ~ "failed"
          )
        ) |>
        count(status, name = "n_combinations") |>
        mutate(
          normalization = normalization,
          analysis = analysis
        ) |>
        select(normalization, analysis, status, n_combinations)
    )
  }

  tibble(
    normalization = normalization,
    analysis = analysis,
    status = "status_column_not_found",
    n_combinations = nrow(df)
  )
}

fit_summary <- bind_rows(
  summarize_fit(primary$amyloid_fit, "Amyloid", "size_factor_mor"),
  summarize_fit(sensitivity$amyloid_fit, "Amyloid", "total_counts"),
  summarize_fit(primary$disease_fit, "Disease", "size_factor_mor"),
  summarize_fit(sensitivity$disease_fit, "Disease", "total_counts"),
  summarize_fit(primary$overall_fit, "Overall", "size_factor_mor"),
  summarize_fit(sensitivity$overall_fit, "Overall", "total_counts"),
  summarize_fit(primary$maxpath_fit, "MaxPathology", "size_factor_mor"),
  summarize_fit(sensitivity$maxpath_fit, "MaxPathology", "total_counts")
)

write_qc(fit_summary, "fit_summary.csv")

extract_fit_problems <- function(df, analysis, normalization) {
  if (is.null(df)) return(tibble())

  status_col <- find_col(
    df,
    c("status", "fit_status", "model_status"),
    required = FALSE
  )

  fit_ok_col <- find_col(
    df,
    c("fit_ok", "success", "model_ok"),
    required = FALSE
  )

  hessian_col <- find_col(
    df,
    c("pd_hessian", "positive_definite_hessian"),
    required = FALSE
  )

  out <- df |>
    mutate(
      analysis = analysis,
      normalization = normalization
    )

  if (!is.na(status_col)) {
    return(
      out |>
        filter(
          !tolower(as.character(.data[[status_col]])) %in%
            c("success", "ok", "fit_ok")
        )
    )
  }

  if (!is.na(fit_ok_col)) {
    fit_ok <- normalize_logical(out[[fit_ok_col]])
    hessian_bad <- if (!is.na(hessian_col)) {
      normalize_logical(out[[hessian_col]]) %in% FALSE
    } else {
      rep(FALSE, nrow(out))
    }

    return(out[fit_ok %in% FALSE | hessian_bad, , drop = FALSE])
  }

  out[0, , drop = FALSE]
}

fit_problems <- bind_rows(
  extract_fit_problems(primary$amyloid_fit, "Amyloid", "size_factor_mor"),
  extract_fit_problems(sensitivity$amyloid_fit, "Amyloid", "total_counts"),
  extract_fit_problems(primary$maxpath_fit, "MaxPathology", "size_factor_mor"),
  extract_fit_problems(sensitivity$maxpath_fit, "MaxPathology", "total_counts")
)

write_qc(fit_problems, "fit_problems.csv")

# ============================================================
# Compare size_factor_mor versus total_counts
# ============================================================

primary_std <- standardize_summary(
  primary$combined,
  "size_factor_mor combined summary"
)

sensitivity_std <- standardize_summary(
  sensitivity$combined,
  "total_counts combined summary"
)

join_keys <- c("type_std", "region_std", "celltype_std", "contrast_std")

primary_compare <- primary_std |>
  select(all_of(join_keys), effect_std, p_value_std, p_adj_std) |>
  rename(
    effect_mor = effect_std,
    p_value_mor = p_value_std,
    p_adj_mor = p_adj_std
  )

sensitivity_compare <- sensitivity_std |>
  select(all_of(join_keys), effect_std, p_value_std, p_adj_std) |>
  rename(
    effect_total_counts = effect_std,
    p_value_total_counts = p_value_std,
    p_adj_total_counts = p_adj_std
  )

duplicate_primary_keys <- primary_compare |>
  count(across(all_of(join_keys)), name = "n") |>
  filter(n > 1)

duplicate_sensitivity_keys <- sensitivity_compare |>
  count(across(all_of(join_keys)), name = "n") |>
  filter(n > 1)

write_qc(duplicate_primary_keys, "duplicate_primary_comparison_keys.csv")
write_qc(duplicate_sensitivity_keys, "duplicate_sensitivity_comparison_keys.csv")

normalization_comparison <- inner_join(
  primary_compare,
  sensitivity_compare,
  by = join_keys
) |>
  mutate(
    effect_difference = effect_mor - effect_total_counts,
    absolute_effect_difference = abs(effect_difference),
    same_direction = case_when(
      is.na(effect_mor) | is.na(effect_total_counts) ~ NA,
      effect_mor == 0 & effect_total_counts == 0 ~ TRUE,
      TRUE ~ sign(effect_mor) == sign(effect_total_counts)
    ),
    significant_mor = !is.na(p_adj_mor) & p_adj_mor < 0.05,
    significant_total_counts =
      !is.na(p_adj_total_counts) & p_adj_total_counts < 0.05,
    significant_both = significant_mor & significant_total_counts,
    significance_agreement =
      significant_mor == significant_total_counts
  )

write_qc(normalization_comparison, "normalization_comparison.csv")

normalization_summary <- normalization_comparison |>
  group_by(type_std) |>
  summarise(
    n_compared = n(),
    spearman_effect_correlation = safe_cor(
      effect_mor,
      effect_total_counts,
      "spearman"
    ),
    pearson_effect_correlation = safe_cor(
      effect_mor,
      effect_total_counts,
      "pearson"
    ),
    median_absolute_effect_difference =
      safe_median(absolute_effect_difference),
    maximum_absolute_effect_difference =
      safe_max(absolute_effect_difference),
    direction_agreement = safe_mean(same_direction),
    significance_agreement = safe_mean(significance_agreement),
    n_significant_mor = sum(significant_mor, na.rm = TRUE),
    n_significant_total_counts =
      sum(significant_total_counts, na.rm = TRUE),
    n_significant_both = sum(significant_both, na.rm = TRUE),
    .groups = "drop"
  ) |>
  rename(type = type_std)

write_qc(normalization_summary, "normalization_summary.csv")

offset_sensitive_results <- normalization_comparison |>
  filter(
    same_direction %in% FALSE |
      significance_agreement %in% FALSE |
      absolute_effect_difference > 0.5
  ) |>
  arrange(desc(absolute_effect_difference))

write_qc(offset_sensitive_results, "offset_sensitive_results.csv")

robust_significant_results <- normalization_comparison |>
  filter(significant_both, same_direction %in% TRUE) |>
  arrange(type_std, region_std, p_adj_mor)

write_qc(robust_significant_results, "robust_significant_results.csv")

primary_only_significant <- normalization_comparison |>
  filter(significant_mor, !significant_total_counts) |>
  arrange(p_adj_mor)

write_qc(primary_only_significant, "significant_primary_only.csv")

total_counts_only_significant <- normalization_comparison |>
  filter(!significant_mor, significant_total_counts) |>
  arrange(p_adj_total_counts)

write_qc(
  total_counts_only_significant,
  "significant_total_counts_only.csv"
)

direction_disagreements <- normalization_comparison |>
  filter(same_direction %in% FALSE) |>
  arrange(desc(absolute_effect_difference))

write_qc(
  direction_disagreements,
  "direction_disagreements_between_offsets.csv"
)

largest_offset_differences <- normalization_comparison |>
  arrange(desc(absolute_effect_difference)) |>
  slice_head(n = 50)

write_qc(
  largest_offset_differences,
  "largest_offset_sensitivity_differences.csv"
)

# ============================================================
# Compare absolute abundance versus proportion
# ============================================================

proportion_file <- first_existing(proportion_candidates)
proportion_comparison_completed <- FALSE

absolute_vs_proportion <- tibble()
absolute_vs_proportion_summary <- tibble()
possible_closure_artifacts <- tibble()
strong_closure_artifact_candidates <- tibble()

if (is.na(proportion_file)) {
  writeLines(
    c(
      "Absolute-versus-proportion comparison was not run.",
      "No proportion-based combined summary was found at:",
      paste0("  - ", proportion_candidates),
      "",
      "Update 'proportion_candidates' near the top of this script and rerun."
    ),
    file.path(qc_dir, "proportion_comparison_not_run.txt")
  )
} else {
  message("Using proportion summary: ", proportion_file)

  proportion_raw <- read_checked(proportion_file)
  proportion_std <- standardize_summary(
    proportion_raw,
    "proportion combined summary"
  )

  absolute_for_prop <- primary_std |>
    select(all_of(join_keys), effect_std, p_value_std, p_adj_std) |>
    rename(
      absolute_effect = effect_std,
      absolute_p_value = p_value_std,
      absolute_p_adj = p_adj_std
    )

  proportion_for_compare <- proportion_std |>
    select(all_of(join_keys), effect_std, p_value_std, p_adj_std) |>
    rename(
      proportion_effect = effect_std,
      proportion_p_value = p_value_std,
      proportion_p_adj = p_adj_std
    )

  absolute_vs_proportion <- inner_join(
    absolute_for_prop,
    proportion_for_compare,
    by = join_keys
  ) |>
    mutate(
      same_direction = case_when(
        is.na(absolute_effect) | is.na(proportion_effect) ~ NA,
        absolute_effect == 0 & proportion_effect == 0 ~ TRUE,
        TRUE ~ sign(absolute_effect) == sign(proportion_effect)
      ),
      absolute_significant =
        !is.na(absolute_p_adj) & absolute_p_adj < 0.05,
      proportion_significant =
        !is.na(proportion_p_adj) & proportion_p_adj < 0.05,
      significance_agreement =
        absolute_significant == proportion_significant
    )

  write_qc(absolute_vs_proportion, "absolute_vs_proportion.csv")

  absolute_vs_proportion_summary <- absolute_vs_proportion |>
    group_by(type_std) |>
    summarise(
      n_compared = n(),
      spearman_effect_correlation = safe_cor(
        absolute_effect,
        proportion_effect,
        "spearman"
      ),
      pearson_effect_correlation = safe_cor(
        absolute_effect,
        proportion_effect,
        "pearson"
      ),
      direction_agreement = safe_mean(same_direction),
      significance_agreement = safe_mean(significance_agreement),
      n_absolute_significant =
        sum(absolute_significant, na.rm = TRUE),
      n_proportion_significant =
        sum(proportion_significant, na.rm = TRUE),
      n_significant_both =
        sum(absolute_significant & proportion_significant, na.rm = TRUE),
      .groups = "drop"
    ) |>
    rename(type = type_std)

  write_qc(
    absolute_vs_proportion_summary,
    "absolute_vs_proportion_summary.csv"
  )

  possible_closure_artifacts <- absolute_vs_proportion |>
    filter(proportion_significant, !absolute_significant) |>
    arrange(proportion_p_adj)

  write_qc(
    possible_closure_artifacts,
    "possible_closure_artifacts.csv"
  )

  strong_closure_artifact_candidates <- absolute_vs_proportion |>
    filter(
      proportion_significant,
      !absolute_significant,
      same_direction %in% FALSE
    ) |>
    arrange(proportion_p_adj)

  write_qc(
    strong_closure_artifact_candidates,
    "strong_closure_artifact_candidates.csv"
  )

  proportion_comparison_completed <- TRUE
}

# ============================================================
# QC status and text report
# ============================================================

required_files_ok <- nrow(missing_required) == 0
gate_columns_ok <- all(gate_column_inventory$present)

all_gate_proceed <- (
  nrow(gate_summary) > 0 &&
    all(tolower(gate_summary$gate_statuses) == "proceed", na.rm = TRUE)
)

all_confirmatory_eligible <- (
  nrow(gate_summary) > 0 &&
    all(gate_summary$all_confirmatory_eligible %in% TRUE)
)

normalization_metrics_ok <- (
  nrow(normalization_summary) > 0 &&
    all(
      is.na(normalization_summary$spearman_effect_correlation) |
        normalization_summary$spearman_effect_correlation >= 0.80
    ) &&
    all(
      is.na(normalization_summary$direction_agreement) |
        normalization_summary$direction_agreement >= 0.90
    )
)

qc_status <- case_when(
  !required_files_ok ~ "FAILED",
  !gate_columns_ok ~ "REVIEW",
  !all_gate_proceed ~ "REVIEW",
  !all_confirmatory_eligible ~ "REVIEW",
  nrow(fit_problems) > 0 ~ "REVIEW",
  !normalization_metrics_ok ~ "REVIEW",
  TRUE ~ "PASSED"
)

qc_status_table <- tibble(
  check = c(
    "Required files present",
    "Calibration-gate columns present",
    "All gate decisions proceed",
    "All results confirmatory eligible",
    "No flagged fit problems",
    "Normalization sensitivity acceptable",
    "Proportion comparison completed"
  ),
  passed = c(
    required_files_ok,
    gate_columns_ok,
    all_gate_proceed,
    all_confirmatory_eligible,
    nrow(fit_problems) == 0,
    normalization_metrics_ok,
    proportion_comparison_completed
  )
)

write_qc(qc_status_table, "qc_status.csv")

pct <- function(x) {
  ifelse(is.na(x), "NA", paste0(round(100 * x, 1), "%"))
}

normalization_lines <- normalization_summary |>
  mutate(
    line = paste0(
      type,
      ": Spearman=", round(spearman_effect_correlation, 3),
      "; direction agreement=", pct(direction_agreement),
      "; significance agreement=", pct(significance_agreement),
      "; significant with both=", n_significant_both
    )
  ) |>
  pull(line)

report_lines <- c(
  "============================================================",
  "Absolute-Abundance QC Report",
  "============================================================",
  paste0("Generated UTC: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("Project root: ", project_root),
  paste0("QC directory: ", qc_dir),
  paste0("Overall QC status: ", qc_status),
  "",
  "ROW COUNTS",
  "----------",
  capture.output(
    print(
      row_counts |>
        filter(table %in% c("amyloid", "disease", "overall", "maxpath", "combined")),
      n = Inf
    )
  ),
  "",
  "CALIBRATION GATE",
  "----------------",
  capture.output(print(gate_summary, n = Inf)),
  "",
  "FIT DIAGNOSTICS",
  "---------------",
  paste0("Flagged fit problems: ", nrow(fit_problems)),
  "",
  "NORMALIZATION SENSITIVITY",
  "-------------------------",
  normalization_lines,
  paste0("Offset-sensitive results: ", nrow(offset_sensitive_results)),
  paste0("Direction disagreements: ", nrow(direction_disagreements)),
  paste0("Robust significant results: ", nrow(robust_significant_results)),
  "",
  "ABSOLUTE VERSUS PROPORTION",
  "--------------------------",
  if (proportion_comparison_completed) {
    c(
      paste0("Proportion summary: ", proportion_file),
      capture.output(print(absolute_vs_proportion_summary, n = Inf)),
      paste0("Possible closure artifacts: ", nrow(possible_closure_artifacts)),
      paste0(
        "Strong closure-artifact candidates: ",
        nrow(strong_closure_artifact_candidates)
      )
    )
  } else {
    "Comparison skipped because no proportion summary was found."
  },
  "",
  "QC CHECKLIST",
  "------------",
  capture.output(print(qc_status_table, n = Inf)),
  "",
  "INTERPRETATION",
  "--------------",
  "Highest confidence: significant with both offsets, same direction, successful fit, gate=proceed.",
  "Offset-sensitive: significant with one offset only, direction reversal, or absolute effect difference > 0.5.",
  "Possible closure artifact: significant in the proportion analysis but not in the absolute-abundance analysis.",
  "Strong closure candidate: proportion-only significance plus a direction reversal.",
  "",
  "============================================================"
)

writeLines(report_lines, file.path(qc_dir, "qc_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))

message("")
message("============================================================")
message("Absolute-abundance QC completed")
message("Overall QC status: ", qc_status)
message("Outputs: ", qc_dir)
message("Main report: ", file.path(qc_dir, "qc_summary.txt"))
message("Matched normalization rows: ", nrow(normalization_comparison))
message("Offset-sensitive results: ", nrow(offset_sensitive_results))
message("Robust significant results: ", nrow(robust_significant_results))
message("Flagged fit problems: ", nrow(fit_problems))

if (proportion_comparison_completed) {
  message("Possible closure artifacts: ", nrow(possible_closure_artifacts))
  message(
    "Strong closure-artifact candidates: ",
    nrow(strong_closure_artifact_candidates)
  )
} else {
  message("Absolute-versus-proportion comparison was skipped.")
}

message("============================================================")





