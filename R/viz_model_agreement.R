#' This function is the daughter of the retired function: R/viz_model_comparison.R
#' Model-agreement visualizations for spatial cell-type contrasts
#'
#' Compares matched contrast estimates across spatial-analysis pipelines.
#' Supported effect columns include `log2_OR` and `log2_fold_change`.
#'
#' Requires:
#'   - R/viz_theme.R
#'   - ggplot2
#'   - dplyr
#'
#' @keywords internal
NULL


# ============================================================
# Internal helpers
# ============================================================

#' Detect the effect-size column in a spatial contrast summary
#'
#' @param df Data frame containing spatial contrast results.
#' @return A list with `col` and `label`.
#' @keywords internal
.effect_size_column_info <- function(
    df,
    effect_col = NULL
) {
  candidates <- c(
    log2_OR = "log2 OR",
    log2_fold_change = "log2 fold-change"
  )

  if (!is.null(effect_col)) {
    if (
      length(effect_col) != 1L ||
        is.na(effect_col) ||
        !nzchar(effect_col)
    ) {
      stop(
        "`effect_col` must be NULL or one non-empty column name.",
        call. = FALSE
      )
    }

    if (!effect_col %in% names(candidates)) {
      stop(
        paste0(
          "Unsupported effect column `",
          effect_col,
          "`. Expected one of: ",
          paste(names(candidates), collapse = ", "),
          "."
        ),
        call. = FALSE
      )
    }

    if (!effect_col %in% colnames(df)) {
      stop(
        paste0(
          "Requested effect column `",
          effect_col,
          "` was not found in the input."
        ),
        call. = FALSE
      )
    }

    return(list(
      col = effect_col,
      label = unname(candidates[[effect_col]])
    ))
  }

  present <- names(candidates)[
    names(candidates) %in% colnames(df)
  ]

  if (length(present) == 0L) {
    stop(
      paste0(
        "No supported effect-size column found. Expected one of: ",
        paste(names(candidates), collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  if (length(present) > 1L) {
    stop(
      paste0(
        "Multiple supported effect-size columns found: ",
        paste(present, collapse = ", "),
        ". Supply `effect_col` explicitly."
      ),
      call. = FALSE
    )
  }

  list(
    col = present,
    label = unname(candidates[[present]])
  )
}


#' Validate a combined spatial-contrast summary
#'
#' @param df Input data frame.
#' @param source_label Label used in error messages.
#' @return Effect-column metadata.
#' @keywords internal
.validate_model_summary <- function(df, source_label = "input") {
  required <- c(
    "celltype",
    "region",
    "contrast",
    "type",
    "p_adj"
  )

  missing_cols <- setdiff(required, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      paste0(
        source_label,
        " is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Prepare one model summary for pairwise joining
#'
#' @param csv Path to a combined spatial-contrast summary.
#' @param suffix Either `"a"` or `"b"`.
#' @param source_label Model label used in messages.
#' @return A standardized data frame.
#' @keywords internal
.prepare_model_summary <- function(
    csv,
    suffix,
    source_label,
    effect_col = NULL
) {
  df <- utils::read.csv(
    csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  .validate_model_summary(
    df,
    source_label = source_label
  )

  effect_info <- .effect_size_column_info(
    df,
    effect_col = effect_col
  )

  out <- df |>
    dplyr::select(
      "celltype",
      "region",
      "contrast",
      "type",
      effect = dplyr::all_of(effect_info$col),
      "p_adj"
    ) |>
    dplyr::mutate(
      effect = suppressWarnings(as.numeric(.data$effect)),
      p_adj = suppressWarnings(as.numeric(.data$p_adj))
    )

  names(out)[names(out) == "effect"] <- paste0("effect_", suffix)
  names(out)[names(out) == "p_adj"] <- paste0("p_adj_", suffix)

  list(
    data = out,
    effect_col = effect_info$col,
    effect_label = effect_info$label
  )
}


#' Create facet-level model-agreement statistics
#'
#' @param joined Matched model-comparison data.
#' @return Data frame with correlation and agreement summaries.
#' @keywords internal
.summarize_model_agreement <- function(joined) {
  joined |>
    dplyr::group_by(.data$type) |>
    dplyr::summarise(
      n = sum(
        is.finite(.data$effect_a) &
          is.finite(.data$effect_b)
      ),
      pearson_r = if (n >= 3) {
        stats::cor(
          .data$effect_a,
          .data$effect_b,
          method = "pearson",
          use = "complete.obs"
        )
      } else {
        NA_real_
      },
      spearman_rho = if (n >= 3) {
        stats::cor(
          .data$effect_a,
          .data$effect_b,
          method = "spearman",
          use = "complete.obs"
        )
      } else {
        NA_real_
      },
      same_direction_pct = {
        valid <- is.finite(.data$effect_a) &
          is.finite(.data$effect_b)

        if (sum(valid) == 0) {
          NA_real_
        } else {
          100 * mean(
            sign(.data$effect_a[valid]) ==
              sign(.data$effect_b[valid])
          )
        }
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      annotation = paste0(
        "n = ", .data$n,
        "\nPearson r = ", formatC(.data$pearson_r, digits = 2, format = "f"),
        "\nSpearman rho = ", formatC(.data$spearman_rho, digits = 2, format = "f"),
        "\nSame direction = ",
        formatC(.data$same_direction_pct, digits = 1, format = "f"),
        "%"
      )
    )
}

#' Validate uniqueness of model-comparison join keys
#'
#' @param df Data frame to validate.
#' @param keys Character vector of key columns.
#' @param context Label used in error messages.
#' @return The input data frame invisibly.
#' @keywords internal
.validate_model_join_keys <- function(
    df,
    keys,
    context = "Model summary"
) {
  duplicate_keys <- df |>
    dplyr::count(
      dplyr::across(dplyr::all_of(keys)),
      name = "n"
    ) |>
    dplyr::filter(.data$n > 1L)

  if (nrow(duplicate_keys) > 0L) {
    stop(
      context,
      " contains ",
      nrow(duplicate_keys),
      " duplicated join-key combination(s) across: ",
      paste(keys, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(df)
}

# ============================================================
# Pairwise effect-size comparison
# ============================================================

#' Compare matched effect estimates between two models
#'
#' @param model_a_csv Path to the first combined contrast summary.
#' @param model_b_csv Path to the second combined contrast summary.
#' @param label_a,label_b Human-readable model labels.
#' @param output_dir Figure output directory.
#' @param padj_cutoff Adjusted-p-value significance threshold.
#' @param filename Optional output filename without extension.
#' @return The ggplot object invisibly, or NULL if inputs are absent.
#' @export
plot_model_effect_agreement <- function(
    model_a_csv,
    model_b_csv,
    label_a,
    label_b,
    output_dir,
    effect_col_a = NULL,
    effect_col_b = NULL,
    padj_cutoff = 0.05,
    filename = NULL
) {
  if (
    length(padj_cutoff) != 1L ||
      !is.finite(padj_cutoff) ||
      padj_cutoff <= 0 ||
      padj_cutoff >= 1
  ) {
    stop(
      "`padj_cutoff` must be one finite number between 0 and 1.",
      call. = FALSE
    )
  }
  if (!file.exists(model_a_csv) || !file.exists(model_b_csv)) {
    message(
      "plot_model_effect_agreement: missing input for ",
      label_a,
      " vs ",
      label_b,
      "."
    )
    return(invisible(NULL))
  }

  model_a <- .prepare_model_summary(
    csv = model_a_csv,
    suffix = "a",
    source_label = label_a,
    effect_col = effect_col_a
  )

  model_b <- .prepare_model_summary(
    csv = model_b_csv,
    suffix = "b",
    source_label = label_b,
    effect_col = effect_col_b
  )

  join_keys <- c(
    "celltype",
    "region",
    "contrast",
    "type"
  )

  .validate_model_join_keys(
    model_a$data,
    keys = join_keys,
    context = label_a
  )

  .validate_model_join_keys(
    model_b$data,
    keys = join_keys,
    context = label_b
  )

  joined <- dplyr::inner_join(
    model_a$data,
    model_b$data,
    by = join_keys
  )

  a_only <- dplyr::anti_join(
    model_a$data,
    model_b$data,
    by = join_keys
  )

  b_only <- dplyr::anti_join(
    model_b$data,
    model_a$data,
    by = join_keys
  )

  message(
    "plot_model_effect_agreement: joined ",
    nrow(joined),
    " row(s); ",
    nrow(a_only),
    " ",
    label_a,
    "-only row(s); ",
    nrow(b_only),
    " ",
    label_b,
    "-only row(s)."
  )

  if (nrow(joined) == 0) {
    message(
      "plot_model_effect_agreement: no overlapping rows for ",
      label_a,
      " vs ",
      label_b,
      "."
    )
    return(invisible(NULL))
  }

  joined <- joined |>
    dplyr::filter(
      is.finite(.data$effect_a),
      is.finite(.data$effect_b)
    ) |>
    dplyr::mutate(
      sig_a = !is.na(.data$p_adj_a) &
        .data$p_adj_a < padj_cutoff,
      sig_b = !is.na(.data$p_adj_b) &
        .data$p_adj_b < padj_cutoff
    )

  if (nrow(joined) == 0) {
    message(
      "plot_model_effect_agreement: no paired finite estimates for ",
      label_a,
      " vs ",
      label_b,
      "."
    )
    return(invisible(NULL))
  }

  effect_limit <- max(
    abs(c(
      joined$effect_a,
      joined$effect_b
    )),
    na.rm = TRUE
  )

  if (!is.finite(effect_limit)) {
    message(
      "plot_model_effect_agreement: ",
      "unable to determine finite plotting limits for ",
      label_a,
      " vs ",
      label_b,
      "."
    )

    return(invisible(NULL))
  }

  if (effect_limit == 0) {
    effect_limit <- 1
  } else {
    effect_limit <- effect_limit * 1.05
  }

  category_levels <- c(
    "Neither",
    paste0(label_a, " only"),
    paste0(label_b, " only"),
    "Both"
  )

  joined <- joined |>
    dplyr::mutate(
      significance = dplyr::case_when(
        .data$sig_a & .data$sig_b ~ "Both",
        .data$sig_a ~ paste0(label_a, " only"),
        .data$sig_b ~ paste0(label_b, " only"),
        TRUE ~ "Neither"
      ),
      significance = factor(
        .data$significance,
        levels = category_levels
      )
    )

  agreement_stats <- .summarize_model_agreement(joined)

  annotation_df <- agreement_stats |>
    dplyr::mutate(
      x = -effect_limit * 0.95,
      y = effect_limit * 0.95
    )

  category_colors <- stats::setNames(
    c("grey70", okabe_ito_palette(3)),
    category_levels
  )

  p <- ggplot2::ggplot(
    joined,
    ggplot2::aes(
      x = .data$effect_a,
      y = .data$effect_b,
      color = .data$significance
    )
  ) +
     ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.5
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dotted",
      color = "grey70",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dotted",
      color = "grey70",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(
      size = 1.7,
      alpha = 0.75
    ) +
    ggplot2::geom_text(
      data = annotation_df,
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        label = .data$annotation
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.2
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(type)
    ) +
    ggplot2::coord_equal(
      xlim = c(-effect_limit, effect_limit),
      ylim = c(-effect_limit, effect_limit)
    ) +
    ggplot2::scale_color_manual(
      values = category_colors,
      name = paste0(
        "Significance\nBH p < ",
        padj_cutoff
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = paste0(
        label_a,
        " versus ",
        label_b
      ),
      subtitle = paste0(
        "Matched cell type-region-contrast estimates; ",
        "solid diagonal denotes identical effect estimates"
      ),
      x = paste0(
        label_a,
        " (",
        model_a$effect_label,
        ")"
      ),
      y = paste0(
        label_b,
        " (",
        model_b$effect_label,
        ")"
      )
    ) +
    theme_analysis()

  if (is.null(filename)) {
    safe_a <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      label_a
    )

    safe_b <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      label_b
    )

    filename <- paste0(
      "model_effect_agreement_",
      safe_a,
      "_vs_",
      safe_b
    )
  }

  save_figure(
    p,
    filename,
    output_dir,
    width = 12,
    height = 8
  )

  invisible(p)
}


#' Summarize significant-result frequency across models
#'
#' This is intended as a descriptive QC comparison. A lower or higher
#' significant-hit rate should not be interpreted by itself as evidence
#' that one model is preferable.
#'
#' @param model_csvs Named character vector of combined-summary paths.
#' @param output_dir Figure output directory.
#' @param padj_cutoff Adjusted-p-value threshold.
#' @return The ggplot object invisibly, or NULL.
#' @export
plot_model_significance_rate_qc <- function(
    model_csvs,
    output_dir,
    padj_cutoff = 0.05
) {
  if (is.null(names(model_csvs)) ||
      any(names(model_csvs) == "")) {
    stop(
      "`model_csvs` must be a named character vector.",
      call. = FALSE
    )
  }

  present <- model_csvs[file.exists(model_csvs)]

  if (length(present) == 0) {
    message(
      "plot_model_significance_rate_qc: ",
      "none of the input files exist."
    )
    return(invisible(NULL))
  }

  model_data <- lapply(
    names(present),
    function(model_label) {
      df <- utils::read.csv(
        present[[model_label]],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      .validate_model_summary(
        df,
        source_label = model_label
      )

      df |>
        dplyr::select(
          "celltype",
          "region",
          "contrast",
          "type",
          "p_adj"
        ) |>
        dplyr::mutate(
          model = model_label,
          p_adj = suppressWarnings(
            as.numeric(.data$p_adj)
          )
        )
    }
  )

  df <- dplyr::bind_rows(model_data)

  summary_df <- df |>
    dplyr::group_by(
      .data$model,
      .data$region,
      .data$type
    ) |>
    dplyr::summarise(
      n_tested = dplyr::n(),
      n_sig = sum(
        !is.na(.data$p_adj) &
          .data$p_adj < padj_cutoff
      ),
      pct_sig = 100 * .data$n_sig / .data$n_tested,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      model = factor(
        .data$model,
        levels = names(model_csvs)
      )
    )

  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = .data$region,
      y = .data$pct_sig,
      fill = .data$model
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.7
    ) +
    ggplot2::facet_wrap(~ type) +
    ggplot2::scale_fill_manual(
      values = okabe_ito_palette(length(present)),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Significant-result frequency across models",
      subtitle = paste0(
        "Descriptive QC summary; BH-adjusted p < ",
        padj_cutoff
      ),
      x = "Region",
      y = "Results significant (%)"
    ) +
    theme_analysis() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 30,
        hjust = 1
      )
    )

  save_figure(
    p,
    "model_significance_rate_qc",
    output_dir,
    width = 11,
    height = 7
  )

  invisible(p)
}

