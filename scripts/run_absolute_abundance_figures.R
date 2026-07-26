# ============================================================
# Figures for the absolute cell-type abundance (total-count offset) pipeline
#
# Covers, for both offset variants (results/spatial_stats_absolute and
# results/spatial_stats_absolute_rawlibsize):
#   - per-pipeline effect-size heatmaps + forest plots, reusing
#     plot_contrast_effsize_heatmap() / plot_contrast_forest()
#     (R/viz_spatial_stats.R) unchanged
#   - the calibration gate diagnostic (results/spatial_stats_absolute only --
#     the gate is run once, upstream of both offset variants)
#   - the flagship absolute-vs-relative log2FC scatter, against
#     results/spatial_stats (the least-confounded proportion-based pipeline)
# ============================================================

# ============================================================
# Libraries
# ============================================================

library(ggplot2)
library(dplyr)

# ============================================================
# Visualization utilities
# ============================================================

source("R/viz_theme.R")
source("R/viz_model_agreement.R")
source("R/viz_spatial_stats.R")
source("R/viz_absolute_abundance.R")


# ============================================================
# Configuration
# ============================================================

fig_root <- file.path(
    "results",
    "figures",
    "abs-abundance-figures"
)

dir.create(
    fig_root,
    recursive = TRUE,
    showWarnings = FALSE
)

# helper: run a figure call only if its input exists
if_exists <- function(path, expr) {
  if (all(file.exists(path))) {
    force(expr)
  } else {
    message("Skipping (missing input): ", paste(path, collapse = ", "))
  }
}

contrast_csvs <- c(
  Amyloid = "amyloid_effect_contrasts.csv",
  Disease = "disease_effect_contrasts.csv",
  Overall = "overall_effect_contrasts.csv",
  MaxPathology = "max_pathology_effect_contrasts.csv"
)

abundance_csv <- file.path(
  "results",
  "cell_proportions",
  "spatial_celltype_proportions_for_R.csv"
)

pipeline_info <- tibble::tribble(
  ~model_label,                    ~stats_dir,                                    ~figure_subdir,
  "Absolute (offset-adjusted)",    "results/spatial_stats_absolute",              "spatial_stats_absolute",
  "Absolute (raw libsize offset)", "results/spatial_stats_absolute_rawlibsize",   "spatial_stats_absolute_rawlibsize"
)

calibration_ratio_csv <- file.path(
  "results",
  "spatial_stats_absolute",
  "calibration_ratio_table.csv"
)

calibration_summary_csv <- file.path(
  "results",
  "spatial_stats_absolute",
  "calibration_gate_summary.csv"
)

relative_summary_csv <- file.path(
  "results",
  "spatial_stats",
  "combined_spatial_contrast_summary.csv"
)

# ============================================================
# Spatial statistics

     # Per-pipeline heatmaps + forest plots
# ============================================================

for (i in seq_len(nrow(pipeline_info))) {
  model_label <- pipeline_info$model_label[[i]]
  stats_dir <- pipeline_info$stats_dir[[i]]
  out_dir <- file.path(
    fig_root,
    pipeline_info$figure_subdir[[i]]
  )

  message(
    "\nGenerating spatial-statistics figures for ",
    model_label,
    "..."
  )

  combined_csv <- file.path(
    stats_dir,
    "combined_spatial_contrast_summary.csv"
  )

  if_exists(
    c(
      combined_csv,
      abundance_csv
    ),
    plot_contrast_effsize_heatmap(
      combined_csv = combined_csv,
      output_dir = out_dir,
      abundance_csv = abundance_csv,
      effect_col = "log2_fold_change",
      effect_label = "log2 fold-change"
    )
  )

# Forest plots are temporarily disabled because the validated main-branch
# R/viz_spatial_stats.R does not currently define plot_contrast_forest().
#
# for (contrast_label in names(contrast_csvs)) {
#   contrast_csv <- file.path(
#     stats_dir,
#     contrast_csvs[[contrast_label]]
#   )
#
#   if_exists(
#     contrast_csv,
#     plot_contrast_forest(
#       contrast_csv = contrast_csv,
#       contrast_label = contrast_label,
#       output_dir = out_dir
#     )
#   )
# }
  

  message(
    "Finished spatial-statistics figures for ",
    model_label,
    "."
  )
}

# ============================================================
# Calibration gate diagnostic
# ============================================================

if_exists(
  c(
    calibration_ratio_csv,
    calibration_summary_csv
  ),
  plot_calibration_gate_diagnostic(
    calib_csv = calibration_ratio_csv,
    gate_summary_csv = calibration_summary_csv,
    output_dir = file.path(
      fig_root,
      "spatial_stats_absolute"
    )
  )
)

# ============================================================
# Absolute-versus-relative comparison
# Flagship: absolute vs. relative log2FC scatter, against the
# least-confounded proportion-based pipeline (raw beta-regression)
# ============================================================

for (i in seq_len(nrow(pipeline_info))) {
  model_label <- pipeline_info$model_label[[i]]
  stats_dir <- pipeline_info$stats_dir[[i]]
  out_dir <- file.path(
    fig_root,
    pipeline_info$figure_subdir[[i]]
  )

  absolute_summary_csv <- file.path(
    stats_dir,
    "combined_spatial_contrast_summary.csv"
  )

  message(
    "\nGenerating absolute-versus-relative comparison for ",
    model_label,
    "..."
  )

  if_exists(
    c(
      absolute_summary_csv,
      relative_summary_csv
    ),
    plot_absolute_vs_relative_fc_scatter(
      absolute_csv = absolute_summary_csv,
      relative_csv = relative_summary_csv,
      output_dir = out_dir
    )
  )
}

print(warnings())

message(
  "\nAll available absolute-abundance figures written under:\n",
  normalizePath(
    fig_root,
    winslash = "/",
    mustWork = FALSE
  )
)
