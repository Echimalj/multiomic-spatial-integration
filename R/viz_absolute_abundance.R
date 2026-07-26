# ============================================================
# Visualization utilities for the absolute-abundance pipeline
#
# Requires:
#   - R/viz_theme.R
#
# Provides:
#   - plot_calibration_gate_diagnostic()
#   - plot_absolute_vs_relative_fc_scatter()
#
# The standard spatial-statistics heatmaps and forest plots remain in:
#   - R/viz_spatial_stats.R
# ============================================================

#' Figures for the absolute cell-type abundance (total-count offset) pipeline
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2 and dplyr.
#'
#' @keywords internal
NULL






# ============================================================
# Internal helpers
# ============================================================

#' Read and validate a CSV used by an absolute-abundance figure
#'
#' @param path Path to the CSV file.
#' @param required_cols Character vector of required column names.
#' @param context Short label used in error messages.
#'
#' @return A tibble.
#' @keywords internal
.read_absolute_figure_csv <- function(path,
                                      required_cols,
                                      context = basename(path)) {
  if (!file.exists(path)) {
    stop(
      context,
      ": input file does not exist: ",
      path,
      call. = FALSE
    )
  }

  dat <- readr::read_csv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  )

  missing_cols <- setdiff(required_cols, names(dat))

  if (length(missing_cols) > 0L) {
    stop(
      context,
      ": missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(dat) == 0L) {
    stop(
      context,
      ": input file contains no rows: ",
      path,
      call. = FALSE
    )
  }

  dat
}


#' Check that join keys uniquely identify rows
#'
#' @param dat Data frame to inspect.
#' @param keys Character vector of key columns.
#' @param context Label used in error messages.
#'
#' @return The input data frame, invisibly.
#' @keywords internal
.validate_unique_figure_keys <- function(dat, keys, context) {
  duplicated_keys <- dat |>
    dplyr::count(
      dplyr::across(dplyr::all_of(keys)),
      name = ".n"
    ) |>
    dplyr::filter(.n > 1L)

  if (nrow(duplicated_keys) > 0L) {
    stop(
      context,
      ": join keys are not unique. Found ",
      nrow(duplicated_keys),
      " duplicated key combination(s) using: ",
      paste(keys, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(dat)
}


# ============================================================
# Calibration gate
# ============================================================
#'
#' The figure a reviewer needs before trusting anything else from the
#' absolute-abundance pipeline: per-ROI `calibration_ratio` (total
#' `abs_abundance` / `total_counts`) by `disease_status`, faceted by region,
#' annotated with the fitted `disease_status` coefficient and the
#' plausibility threshold used to make the gate decision.
#'
#' @param calib_csv Path to `calibration_ratio_table.csv`
#'   (`compute_calibration_ratio_table()` output).
#' @param gate_summary_csv Path to `calibration_gate_summary.csv`
#'   (`run_calibration_gate_check()` output, written by
#'   `scripts/run_calibration_gate_check.R`).
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly), or NULL if inputs are missing.
#' @export
plot_calibration_gate_diagnostic <- function(calib_csv,
                                             gate_summary_csv,
                                             output_dir) {
  input_paths <- c(calib_csv, gate_summary_csv)

  if (!all(file.exists(input_paths))) {
    missing_paths <- input_paths[!file.exists(input_paths)]

    message(
      "plot_calibration_gate_diagnostic: skipping; missing input(s): ",
      paste(missing_paths, collapse = ", ")
    )

    return(invisible(NULL))
  }

  calib_required <- c(
    "disease_status",
    "region",
    "calibration_ratio"
  )

  summary_required <- c(
    "decision",
    "disease_coefficient_log2",
    "biological_threshold_log2"
  )

  calib_df <- .read_absolute_figure_csv(
    path = calib_csv,
    required_cols = calib_required,
    context = "Calibration-ratio table"
  )

  gate_summary <- .read_absolute_figure_csv(
    path = gate_summary_csv,
    required_cols = summary_required,
    context = "Calibration-gate summary"
  )

  if (nrow(gate_summary) > 1L) {
    message(
      "plot_calibration_gate_diagnostic: gate summary contains ",
      nrow(gate_summary),
      " rows; using the first row."
    )
  }

  decision <- as.character(gate_summary$decision[[1]])
  coefficient_log2 <- as.numeric(
    gate_summary$disease_coefficient_log2[[1]]
  )
  threshold_log2 <- as.numeric(
    gate_summary$biological_threshold_log2[[1]]
  )

  if (
    !is.finite(coefficient_log2) ||
      !is.finite(threshold_log2) ||
      is.na(decision) ||
      !nzchar(decision)
  ) {
    stop(
      "Calibration-gate summary contains an invalid decision, coefficient, ",
      "or biological threshold.",
      call. = FALSE
    )
  }

  calib_df <- calib_df |>
    dplyr::mutate(
      calibration_ratio = as.numeric(calibration_ratio),
      disease_status = factor(
        disease_status,
        levels = c("Control", "AD-CAA")
      )
    ) |>
    dplyr::filter(
      is.finite(calibration_ratio),
      !is.na(disease_status),
      !is.na(region)
    )

  if (nrow(calib_df) == 0L) {
    message(
      "plot_calibration_gate_diagnostic: skipping; ",
      "no finite calibration ratios remain after filtering."
    )

    return(invisible(NULL))
  }

  message(
    "plot_calibration_gate_diagnostic: plotting ",
    nrow(calib_df),
    " ROI(s) across ",
    dplyr::n_distinct(calib_df$region),
    " region(s)."
  )

  subtitle_text <- paste0(
    "Disease coefficient = ",
    formatC(coefficient_log2, digits = 3, format = "f"),
    " log2; biological threshold = ",
    formatC(threshold_log2, digits = 3, format = "f"),
    "; decision = ",
    decision
  )

  p <- ggplot2::ggplot(
    calib_df,
    ggplot2::aes(
      x = disease_status,
      y = calibration_ratio,
      color = disease_status
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "grey65",
      linewidth = 0.5
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      height = 0,
      alpha = 0.6,
      size = 1.4
    ) +
    ggplot2::geom_boxplot(
      width = 0.4,
      alpha = 0,
      outlier.shape = NA,
      color = "grey30"
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(region)
    ) +
    ggplot2::scale_color_manual(
      values = okabe_ito_palette(2),
      name = NULL,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Calibration ratio by disease status",
      subtitle = subtitle_text,
      x = NULL,
      y = "Calibration ratio\n(total absolute abundance / total counts)",
      caption = "Dashed horizontal line indicates a calibration ratio of 1."
    ) +
    theme_analysis() +
    ggplot2::theme(
      legend.position = "none"
    )

  save_figure(
    plot = p,
    filename = "calibration_gate_diagnostic",
    output_dir = output_dir,
    width = 10,
    height = 5
  )

  invisible(p)
}


# ============================================================
# Absolute-versus-relative comparison (backward-compatible wrapper)
# ============================================================
#'
#' Compare absolute- and relative-abundance effect estimates.
#'
#' Backward-compatible wrapper around `plot_model_effect_agreement()`
#' configured to compare the standard absolute-abundance pipeline
#' (`log2_fold_change`) against the standard relative-abundance pipeline
#' (`log2_OR`). This preserves the historical API while delegating all
#' plotting logic to the generic model-agreement framework.
#'
#' The resulting scatterplot joins matched cell type × region × contrast ×
#' analysis-type results and compares absolute abundance effects (x-axis)
#' against relative abundance effects (y-axis), highlighting agreement in
#' both effect size and statistical significance.
#'
#' @param absolute_csv Path to the absolute-abundance
#'   `combined_spatial_contrast_summary.csv`.
#' @param relative_csv Path to the relative-abundance
#'   `combined_spatial_contrast_summary.csv`.
#'   Defaults to
#'   `results/spatial_stats/combined_spatial_contrast_summary.csv`.
#' @param output_dir Directory where figures will be written.
#' @param padj_cutoff Adjusted p-value threshold used to define significant
#'   results.
#'
#' @return Invisibly returns the ggplot object produced by
#'   `plot_model_effect_agreement()`, or `NULL` if required inputs are
#'   unavailable.
#'
#' @seealso [plot_model_effect_agreement()]
#' @export

plot_absolute_vs_relative_fc_scatter <- function(
    absolute_csv,
    relative_csv =
      "results/spatial_stats/combined_spatial_contrast_summary.csv",
    output_dir,
    padj_cutoff = 0.05
) {
  plot_model_effect_agreement(
    model_a_csv = absolute_csv,
    model_b_csv = relative_csv,
    label_a = "Absolute abundance",
    label_b = "Relative abundance",
    output_dir = output_dir,
    effect_col_a = "log2_fold_change",
    effect_col_b = "log2_OR",
    padj_cutoff = padj_cutoff,
    filename = "absolute_vs_relative_fc_scatter"
  )
}
