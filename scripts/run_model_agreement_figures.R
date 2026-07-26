# Daugther function of retired function: scripts/run_model_comparison_figures.R
# ============================================================
# Cross-model agreement figures
#
# Compares matched spatial cell-type contrast estimates across:
#   1. Relative abundance
#   2. Absolute abundance with adjusted offset
#   3. Absolute abundance with raw-library-size offset
#
# Fibroblast-excluded analyses are intentionally not included.
# ============================================================

library(ggplot2)
library(dplyr)

source("R/viz_theme.R")
source("R/viz_model_agreement.R")


# ============================================================
# Configuration
# ============================================================

fig_root <- file.path(
  "results",
  "figures",
  "model_agreement"
)

dir.create(
  fig_root,
  recursive = TRUE,
  showWarnings = FALSE
)

model_info <- tibble::tribble(
  ~model_id, ~model_label, ~stats_dir, ~effect_col,

  "relative",
  "Relative abundance",
  file.path("results", "spatial_stats"),
  "log2_OR",

  "absolute_adjusted",
  "Absolute abundance (offset-adjusted)",
  file.path("results", "spatial_stats_absolute"),
  "log2_fold_change",

  "absolute_rawlibsize",
  "Absolute abundance (raw libsize offset)",
  file.path("results", "spatial_stats_absolute_rawlibsize"),
  "log2_fold_change"
) |>
  dplyr::mutate(
    combined_csv = file.path(
      stats_dir,
      "combined_spatial_contrast_summary.csv"
    )
  )

# ============================================================
# Helpers
# ============================================================

if_exists <- function(paths, expr) {
  missing_paths <- paths[!file.exists(paths)]

  if (length(missing_paths) > 0) {
    message(
      "Skipping (missing input): ",
      paste(missing_paths, collapse = ", ")
    )
    return(invisible(NULL))
  }

  force(expr)
}


get_model_value <- function(id, column) {

  value <- model_info |>
    dplyr::filter(.data$model_id == id) |>
    dplyr::pull(dplyr::all_of(column))

  if (length(value) != 1) {
    stop(
      "Expected exactly one model configuration for `",
      id,
      "`.",
      call. = FALSE
    )
  }

  value
}


# ============================================================
# Descriptive significance-rate QC
# ============================================================

named_model_csvs <- stats::setNames(
  model_info$combined_csv,
  model_info$model_label
)

present_model_csvs <- named_model_csvs[
  file.exists(named_model_csvs)
]

if (length(present_model_csvs) > 0) {
  plot_model_significance_rate_qc(
    model_csvs = present_model_csvs,
    output_dir = fig_root
  )
} else {
  message(
    "Skipping significance-rate QC: ",
    "no model summaries were found."
  )
}


# ============================================================
# Pairwise model comparisons
# ============================================================

comparison_info <- tibble::tribble(
  ~model_a, ~model_b, ~filename,
  "relative",
  "absolute_adjusted",
  "relative_vs_absolute_adjusted",

  "relative",
  "absolute_rawlibsize",
  "relative_vs_absolute_rawlibsize",

  "absolute_adjusted",
  "absolute_rawlibsize",
  "absolute_adjusted_vs_absolute_rawlibsize"
)


for (i in seq_len(nrow(comparison_info))) {
  model_a_id <- comparison_info$model_a[[i]]
  model_b_id <- comparison_info$model_b[[i]]

  csv_a <- get_model_value(
    model_a_id,
    "combined_csv"
  )

  csv_b <- get_model_value(
    model_b_id,
    "combined_csv"
  )

  label_a <- get_model_value(
    model_a_id,
    "model_label"
  )

  label_b <- get_model_value(
    model_b_id,
    "model_label"
  )

  effect_col_a <- get_model_value(
    model_a_id,
    "effect_col"
  )

  effect_col_b <- get_model_value(
    model_b_id,
    "effect_col"
  )

  message(
    "\nGenerating model-agreement figure for ",
    label_a,
    " versus ",
    label_b,
    "..."
  )

  if_exists(
    c(csv_a, csv_b),
    plot_model_effect_agreement(
        model_a_csv = csv_a,
        model_b_csv = csv_b,
        label_a = label_a,
        label_b = label_b,
        output_dir = fig_root,
        effect_col_a = effect_col_a,
        effect_col_b = effect_col_b,
        filename = comparison_info$filename[[i]]
    )
  )
}


message(
  "\nAll available model-agreement figures written under:\n",
  normalizePath(
    fig_root,
    winslash = "/",
    mustWork = FALSE
  )
)
