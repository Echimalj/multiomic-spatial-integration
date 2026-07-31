#' Figures for the deconvolution method comparison
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2. There is no ground
#' truth for the real ROIs, so these visualize agreement with the
#' Cell2location baseline and between methods, not accuracy.
#'
#' @keywords internal
NULL


#' Method x method correlation heatmap
#'
#' @param cross_cor_csv Path to `cross_method_correlation.csv` (a square matrix
#'   written with row names).
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly).
#' @export
plot_cross_method_heatmap <- function(cross_cor_csv, output_dir) {
  mat <- as.matrix(utils::read.csv(cross_cor_csv, row.names = 1, check.names = FALSE))

  long <- data.frame(
    method_a = rep(rownames(mat), times = ncol(mat)),
    method_b = rep(colnames(mat), each = nrow(mat)),
    corr = as.vector(mat),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(long, ggplot2::aes(x = method_a, y = method_b, fill = corr)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", corr)), size = 3, color = "grey15") +
    scale_fill_diverging(name = "corr.", limits = c(-1, 1)) +
    ggplot2::labs(
      title = "Agreement between deconvolution methods",
      subtitle = "Correlation of full ROI x cell-type proportion vectors",
      x = NULL, y = NULL
    ) +
    theme_analysis() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  save_figure(p, "cross_method_correlation", output_dir, width = 8, height = 7)
  invisible(p)
}


#' Per-method agreement with the baseline (distribution over cell types)
#'
#' Reads `concordance_vs_cell2location.csv` and draws, per method, the
#' distribution of per-cell-type correlations against the Cell2location
#' baseline. A method whose box sits high and tight agrees closely across
#' cell types.
#'
#' @param concordance_csv Path to `concordance_vs_cell2location.csv`.
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly).
#' @export
plot_method_agreement_box <- function(concordance_csv, output_dir) {
  df <- utils::read.csv(concordance_csv, stringsAsFactors = FALSE)

  med <- stats::aggregate(pearson_r ~ method, df, median, na.rm = TRUE)
  method_order <- med$method[order(med$pearson_r)]
  df$method <- factor(df$method, levels = method_order)

  methods_n <- length(method_order)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = pearson_r, y = method, fill = method)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = okabe_ito_palette(methods_n), guide = "none") +
    ggplot2::labs(
      title = "Per-method agreement with Cell2location",
      subtitle = "Distribution of per-cell-type correlations across ROIs",
      x = "Correlation with baseline", y = NULL
    ) +
    theme_analysis()

  save_figure(p, "method_agreement_distribution", output_dir, width = 8, height = 6)
  invisible(p)
}


#' Load every method's long proportions plus the Cell2location baseline
#'
#' @keywords internal
.load_combined_proportions <- function(method_output_dir, baseline_csv,
                                       roi_id_col = "ROI_ID",
                                       celltype_col = "celltype",
                                       proportion_col = "rel_abundance") {
  files <- list.files(method_output_dir, pattern = "_proportions\\.csv$", full.names = TRUE)
  method_dfs <- lapply(files, utils::read.csv, stringsAsFactors = FALSE)
  combined <- do.call(rbind, method_dfs)

  baseline <- utils::read.csv(baseline_csv, stringsAsFactors = FALSE)
  baseline <- data.frame(
    method = "Cell2location",
    ROI_ID = baseline[[roi_id_col]],
    celltype = baseline[[celltype_col]],
    proportion = baseline[[proportion_col]],
    stringsAsFactors = FALSE
  )

  rbind(combined, baseline)
}


#' Per-cell-type agreement heatmap (method x cell-type correlation vs baseline)
#'
#' More diagnostic than the agreement boxplot: shows *which* cell types drive
#' each method's agreement or disagreement with Cell2location.
#'
#' @param concordance_csv Path to `concordance_vs_cell2location.csv`.
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly).
#' @export
plot_concordance_heatmap <- function(concordance_csv, output_dir) {
  df <- utils::read.csv(concordance_csv, stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = celltype, fill = pearson_r)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_diverging(name = "corr. vs\nbaseline", limits = c(-1, 1)) +
    ggplot2::labs(
      title = "Per-cell-type agreement with Cell2location",
      x = NULL, y = NULL
    ) +
    theme_analysis(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  save_figure(p, "concordance_heatmap", output_dir, width = 8, height = 10)
  invisible(p)
}


#' Cross-method spread of mean cell-type proportion
#'
#' For each cell type, the distribution of per-method mean proportions (one
#' point per method). Cell types with a tight box are method-robust (all tools
#' agree); wide boxes flag method-sensitive cell types to interpret cautiously.
#' Methods are not individually color-coded (there can be > 8, beyond the safe
#' categorical limit) -- the spread itself is the message.
#'
#' @param method_output_dir Directory holding `<method>_proportions.csv`.
#' @param baseline_csv Cell2location long proportions CSV.
#' @param output_dir Figures output directory.
#' @param top_n If set, show only the `top_n` most method-sensitive cell types.
#' @return The ggplot object (invisibly).
#' @export
plot_celltype_method_spread <- function(method_output_dir, baseline_csv, output_dir,
                                        top_n = NULL) {
  combined <- .load_combined_proportions(method_output_dir, baseline_csv)
  mean_by <- stats::aggregate(proportion ~ method + celltype, combined, mean)

  spread <- stats::aggregate(proportion ~ celltype, mean_by, function(x) diff(range(x)))
  colnames(spread)[2] <- "range"
  order_ct <- spread$celltype[order(spread$range)]
  if (!is.null(top_n)) order_ct <- utils::tail(order_ct, top_n)

  mean_by <- mean_by[mean_by$celltype %in% order_ct, ]
  mean_by$celltype <- factor(mean_by$celltype, levels = order_ct)

  p <- ggplot2::ggplot(mean_by, ggplot2::aes(x = proportion, y = celltype)) +
    ggplot2::geom_boxplot(outlier.shape = NA, color = "grey55", fill = "grey95") +
    ggplot2::geom_jitter(height = 0.15, size = 1.1, alpha = 0.7, color = "#0072B2") +
    ggplot2::labs(
      title = "Cross-method spread of mean cell-type proportion",
      subtitle = "Each point = one method's mean; wide spread = method-sensitive cell type",
      x = "Mean proportion (per method)", y = NULL
    ) +
    theme_analysis()

  save_figure(p, "celltype_method_spread", output_dir, width = 8, height = 10)
  invisible(p)
}


#' Per-method distance-from-baseline bar
#'
#' Ranked mean absolute proportion difference vs. the Cell2location baseline
#' (lower = closer).
#'
#' @param mad_csv Path to `mean_abs_diff_vs_cell2location.csv`.
#' @param output_dir Figures output directory.
#' @return The ggplot object (invisibly).
#' @export
plot_mean_abs_diff_bar <- function(mad_csv, output_dir) {
  df <- utils::read.csv(mad_csv, stringsAsFactors = FALSE)
  df <- df[order(df$mean_abs_diff), ]
  df$method <- factor(df$method, levels = rev(df$method))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_abs_diff, y = method, fill = mean_abs_diff)) +
    ggplot2::geom_col(width = 0.7) +
    scale_fill_magnitude(name = NULL) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
      title = "Distance from Cell2location baseline",
      subtitle = "Mean absolute proportion difference (lower = closer)",
      x = "Mean |proportion difference|", y = NULL
    ) +
    theme_analysis()

  save_figure(p, "mean_abs_diff_bar", output_dir, width = 7, height = 5)
  invisible(p)
}


#' Hierarchical clustering dendrogram of methods
#'
#' Clusters methods on 1 - (cross-method correlation), showing which tools
#' group together. Base-graphics (no extra package dependency); written to PDF
#' and PNG directly.
#'
#' @param cross_cor_csv Path to `cross_method_correlation.csv`.
#' @param output_dir Figures output directory.
#' @return Invisibly, the hclust object.
#' @export
plot_method_dendrogram <- function(cross_cor_csv, output_dir) {
  mat <- as.matrix(utils::read.csv(cross_cor_csv, row.names = 1, check.names = FALSE))
  d <- stats::as.dist(1 - mat)
  hc <- stats::hclust(d, method = "average")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(output_dir, "method_dendrogram.pdf")
  grDevices::pdf(pdf_path, width = 8, height = 6)
  graphics::plot(hc, main = "Method clustering (1 - correlation)", xlab = "", sub = "")
  grDevices::dev.off()

  png_path <- file.path(output_dir, "method_dendrogram.png")
  grDevices::png(png_path, width = 8 * 200, height = 6 * 200, res = 200)
  graphics::plot(hc, main = "Method clustering (1 - correlation)", xlab = "", sub = "")
  grDevices::dev.off()

  message("Wrote ", pdf_path, " and .png")
  invisible(hc)
}


#' Per-cell-type scatter of each method against the baseline
#'
#' Direct agreement view: for a handful of cell types, each method's per-ROI
#' proportion (y) vs. the Cell2location proportion (x), faceted by cell type x
#' method, with a y=x reference line.
#'
#' @param method_output_dir Directory holding `<method>_proportions.csv`.
#' @param baseline_csv Cell2location long proportions CSV.
#' @param output_dir Figures output directory.
#' @param celltypes Cell types to show; if NULL, the `n_celltypes` most
#'   abundant (by baseline mean) are used.
#' @param n_celltypes Number of cell types when `celltypes` is NULL.
#' @return The ggplot object (invisibly).
#' @export

plot_celltype_scatter_vs_baseline <- function(
    method_output_dir,
    baseline_csv,
    output_dir,
    celltypes = NULL,
    n_celltypes = 6
) {
  combined <- .load_combined_proportions(
    method_output_dir,
    baseline_csv
  )

  baseline <- combined[
    combined$method == "Cell2location",
    ,
    drop = FALSE
  ]

  others <- combined[
    combined$method != "Cell2location",
    ,
    drop = FALSE
  ]

  if (is.null(celltypes)) {
    baseline_for_ranking <- baseline[
      is.finite(baseline$proportion),
      ,
      drop = FALSE
    ]

    ct_mean <- stats::aggregate(
      proportion ~ celltype,
      baseline_for_ranking,
      mean
    )

    ct_mean <- ct_mean[
      order(-ct_mean$proportion),
      ,
      drop = FALSE
    ]

    celltypes <- utils::head(
      ct_mean$celltype,
      n_celltypes
    )
  }

  b <- baseline[
    baseline$celltype %in% celltypes,
    c(
      "ROI_ID",
      "celltype",
      "proportion"
    ),
    drop = FALSE
  ]

  colnames(b)[3] <- "baseline"

  o <- others[
    others$celltype %in% celltypes,
    ,
    drop = FALSE
  ]

  merged <- merge(
    o,
    b,
    by = c(
      "ROI_ID",
      "celltype"
    )
  )

  n_before <- nrow(merged)

  merged <- merged[
    is.finite(merged$baseline) &
      is.finite(merged$proportion),
    ,
    drop = FALSE
  ]

  n_removed <- n_before - nrow(merged)

  if (n_removed > 0) {
    message(
      "plot_celltype_scatter_vs_baseline: excluded ",
      n_removed,
      " rows without paired finite estimates."
    )
  }

  if (nrow(merged) == 0) {
    warning(
      "plot_celltype_scatter_vs_baseline: no paired finite estimates remain."
    )

    return(
      invisible(NULL)
    )
  }

  merged$celltype <- factor(
    merged$celltype,
    levels = celltypes
  )

panel_cor <- merged |>
  dplyr::group_by(
    celltype,
    method
  ) |>
  dplyr::summarise(
    n_pairs = sum(
      is.finite(baseline) &
        is.finite(proportion)
    ),
    rho = if (
      n_pairs >= 3 &&
        stats::sd(baseline, na.rm = TRUE) > 0 &&
        stats::sd(proportion, na.rm = TRUE) > 0
    ) {
      suppressWarnings(
        stats::cor(
          baseline,
          proportion,
          method = "spearman",
          use = "complete.obs"
        )
      )
    } else {
      NA_real_
    },
    .groups = "drop"
  )

panel_cor$label <- ifelse(
  is.finite(panel_cor$rho),
  paste0(
    "rho = ",
    sprintf("%.2f", panel_cor$rho)
  ),
  "rho = NA"
)

p <- ggplot2::ggplot(
  merged,
  ggplot2::aes(
    x = baseline,
    y = proportion
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "22",
    color = "grey40",
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    alpha = 0.4,
    size = 0.6,
    color = "#0072B2"
  ) +
  ggplot2::geom_text(
    data = panel_cor,
    mapping = ggplot2::aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.2,
    size = 2.5,
    fontface = "italic",
    colour = "grey25"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      NA
    ),
    expand = ggplot2::expansion(
      mult = c(0.05, 0.05)
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      NA
    ),
    expand = ggplot2::expansion(
      mult = c(0.05, 0.08)
    )
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(
      celltype,
      method
    ),
    scales = "free_y",
    ncol = length(
      unique(merged$method)
    )
  ) +
    ggplot2::labs(
      title = paste0(
        "Agreement between Cell2location and alternative ",
        "deconvolution methods"
      ),
      subtitle = paste0(
        "Top ",
        length(celltypes),
        " cell types by Cell2location abundance; ",
        "dashed line indicates perfect agreement (y = x)"
      ),
      x = "Cell2location proportion",
      y = "Alternative-method proportion"
    ) +
    theme_analysis(
      base_size = 8
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 5
      ),
      axis.text.y = ggplot2::element_text(
        size = 5
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 7
      )
    )

  save_figure(
    p,
    "celltype_scatter_vs_baseline",
    output_dir,
    width = 12,
    height = 11
  )

  invisible(p)
}

#' Method x cell-type mean-composition heatmap
#'
#' Loads every `<method>_proportions.csv` plus the Cell2location baseline and
#' draws the mean inferred proportion per method x cell type -- surfacing
#' systematic composition differences between methods at a glance.
#'
#' @param method_output_dir Directory holding `<method>_proportions.csv`.
#' @param baseline_csv Cell2location long proportions CSV.
#' @param output_dir Figures output directory.
#' @param roi_id_col,celltype_col,proportion_col Baseline column names.
#' @return The ggplot object (invisibly).
#' @export
plot_method_composition_heatmap <- function(
    method_output_dir,
    baseline_csv,
    output_dir,
    roi_id_col = "ROI_ID",
    celltype_col = "celltype",
    proportion_col = "rel_abundance"
) {

  files <- list.files(
    method_output_dir,
    pattern = "_proportions\\.csv$",
    full.names = TRUE
  )

  method_dfs <- lapply(
    files,
    utils::read.csv,
    stringsAsFactors = FALSE
  )

  combined <- do.call(
    rbind,
    method_dfs
  )

  baseline <- utils::read.csv(
    baseline_csv,
    stringsAsFactors = FALSE
  )

  baseline <- data.frame(
    method = "Cell2location",
    ROI_ID = baseline[[roi_id_col]],
    celltype = baseline[[celltype_col]],
    proportion = baseline[[proportion_col]],
    stringsAsFactors = FALSE
  )

  combined <- rbind(
    combined,
    baseline
  )

  # Remove non-finite proportions before summarizing
  n_before <- nrow(combined)

  combined <- combined[
    is.finite(combined$proportion),
    ,
    drop = FALSE
  ]

  n_removed <- n_before - nrow(combined)

  if (n_removed > 0) {
    message(
      "plot_method_composition_heatmap: excluded ",
      n_removed,
      " rows with non-finite proportions."
    )
  }

  mean_comp <- stats::aggregate(
    proportion ~ method + celltype,
    combined,
    function(x) {
      mean(
        x,
        na.rm = TRUE
      )
    }
  )
  # Order cell types by cross-method disagreement
  celltype_spread <- stats::aggregate(
    proportion ~ celltype,
    mean_comp,
    function(x) {
      max(x, na.rm = TRUE) -
        min(x, na.rm = TRUE)
    }
  )
  
  celltype_order <- celltype_spread$celltype[
    order(
      celltype_spread$proportion,
      decreasing = TRUE
    )
  ]
  
  mean_comp$celltype <- factor(
    mean_comp$celltype,
    levels = rev(celltype_order)
  )

  p <- ggplot2::ggplot(
    mean_comp,
    ggplot2::aes(
      x = method,
      y = celltype,
      fill = proportion
    )
  ) +
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::scale_fill_viridis_c(
      name = "Mean proportion",
      trans = "sqrt",
      limits = c(0, NA)
    ) +
    ggplot2::labs(
      title = "Mean inferred cell-type composition by method",
      subtitle = "Cell types ordered by cross-method variability",
      x = NULL,
      y = NULL
    ) +
    theme_analysis(
      base_size = 8
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
		vjust = 1,
		face = "bold"
      )
    )

  save_figure(
    p,
    "method_composition_heatmap",
    output_dir,
    width = 8,
    height = 10
  )

  invisible(p)
}

# ============================================================
# Mean cell-type abundance and method-bias heatmaps
#
# Panel A:
#   Mean estimated proportion for each cell type and method.
#
# Panel B:
#   Row-scaled z-score across methods. This highlights methods
#   that estimate relatively more or less of a given cell type.
#
# The same clustered row and column ordering is used in both
# panels so they can be compared directly.
# ============================================================

plot_method_celltype_mean_heatmaps <- function(
    deconv_dir,
    baseline_csv,
    output_dir,
    cluster_rows = TRUE,
    cluster_columns = TRUE
) {

  required_packages <- c(
    "ggplot2",
    "dplyr",
    "readr",
    "tidyr",
    "grid"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # ----------------------------------------------------------
  # Find method proportion files
  # ----------------------------------------------------------

  method_files <- list.files(
    deconv_dir,
    pattern = "_proportions\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  method_files <- method_files[
    !grepl(
      "Cell2location",
      basename(method_files),
      ignore.case = TRUE
    )
  ]

  if (length(method_files) == 0) {
    stop(
      "No method proportion files were found in ",
      deconv_dir,
      ".",
      call. = FALSE
    )
  }

  if (!file.exists(baseline_csv)) {
    stop(
      "Cell2location baseline file not found: ",
      baseline_csv,
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Load Cell2location baseline
  # ----------------------------------------------------------

  baseline_raw <- readr::read_csv(
    baseline_csv,
    show_col_types = FALSE
  )

  required_baseline_columns <- c(
    "ROI_ID",
    "celltype",
    "rel_abundance"
  )

  missing_baseline_columns <- setdiff(
    required_baseline_columns,
    names(baseline_raw)
  )

  if (length(missing_baseline_columns) > 0) {
    stop(
      "Cell2location baseline is missing required columns: ",
      paste(missing_baseline_columns, collapse = ", "),
      call. = FALSE
    )
  }

  baseline <- baseline_raw |>
    dplyr::transmute(
      method = "Cell2location",
      ROI_ID = as.character(.data$ROI_ID),
      celltype = as.character(.data$celltype),
      proportion = as.numeric(.data$rel_abundance)
    )

  # ----------------------------------------------------------
  # Load comparison methods
  # ----------------------------------------------------------

  method_list <- lapply(
    method_files,
    function(method_file) {

      method_df <- readr::read_csv(
        method_file,
        show_col_types = FALSE
      )

      required_method_columns <- c(
        "ROI_ID",
        "celltype",
        "proportion"
      )

      missing_method_columns <- setdiff(
        required_method_columns,
        names(method_df)
      )

      if (length(missing_method_columns) > 0) {
        warning(
          "Skipping ",
          basename(method_file),
          ": missing columns ",
          paste(missing_method_columns, collapse = ", "),
          call. = FALSE
        )

        return(NULL)
      }

      filename_method <- sub(
        "_proportions\\.csv$",
        "",
        basename(method_file)
      )

      method_label <- filename_method

      if ("method" %in% names(method_df)) {
        observed_labels <- unique(
          method_df$method[
            !is.na(method_df$method) &
              nzchar(method_df$method)
          ]
        )

        if (length(observed_labels) == 1) {
          method_label <- observed_labels
        }
      }

      method_df |>
        dplyr::transmute(
          method = method_label,
          ROI_ID = as.character(.data$ROI_ID),
          celltype = as.character(.data$celltype),
          proportion = as.numeric(.data$proportion)
        )
    }
  )

  method_list <- method_list[
    !vapply(
      method_list,
      is.null,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(method_list) == 0) {
    stop(
      "No valid comparison-method files could be loaded.",
      call. = FALSE
    )
  }

  all_proportions <- dplyr::bind_rows(
    c(
      list(baseline),
      method_list
    )
  )

  # ----------------------------------------------------------
  # Calculate mean abundance by method and cell type
  #
  # For methods with failed ROIs, the mean is calculated using
  # the finite estimates available for that method/cell type.
  # ----------------------------------------------------------

  mean_abundance <- all_proportions |>
    dplyr::group_by(
      .data$celltype,
      .data$method
    ) |>
    dplyr::summarise(
      n_evaluable = sum(
        is.finite(.data$proportion)
      ),
      mean_proportion = if (
        all(!is.finite(.data$proportion))
      ) {
        NA_real_
      } else {
        mean(
          .data$proportion[
            is.finite(.data$proportion)
          ],
          na.rm = TRUE
        )
      },
      .groups = "drop"
    )

  # Save the values underlying the figure.
  readr::write_csv(
    mean_abundance,
    file.path(
      deconv_dir,
      "all_methods_mean_proportion_by_celltype.csv"
    )
  )

  # ----------------------------------------------------------
  # Convert to cell type x method matrix
  # ----------------------------------------------------------

  mean_wide <- mean_abundance |>
    dplyr::select(
      celltype,
      method,
      mean_proportion
    ) |>
    tidyr::pivot_wider(
      names_from = method,
      values_from = mean_proportion
    ) |>
    dplyr::arrange(
      .data$celltype
    )

  mean_matrix <- as.matrix(
    mean_wide[
      ,
      setdiff(
        names(mean_wide),
        "celltype"
      ),
      drop = FALSE
    ]
  )

  storage.mode(mean_matrix) <- "numeric"

  rownames(mean_matrix) <- mean_wide$celltype

  # ----------------------------------------------------------
  # Row-scaled z-scores
  # ----------------------------------------------------------

  z_matrix <- matrix(
    NA_real_,
    nrow = nrow(mean_matrix),
    ncol = ncol(mean_matrix),
    dimnames = dimnames(mean_matrix)
  )

  for (i in seq_len(nrow(mean_matrix))) {

    values <- mean_matrix[i, ]
    finite_values <- values[is.finite(values)]

    if (length(finite_values) == 0) {
      next
    }

    row_mean <- mean(finite_values)
    row_sd <- stats::sd(finite_values)

    if (
      !is.finite(row_sd) ||
      row_sd == 0
    ) {
      z_matrix[i, is.finite(values)] <- 0
    } else {
      z_matrix[i, ] <- (
        values - row_mean
      ) / row_sd
    }
  }

  # ----------------------------------------------------------
  # Determine clustered ordering
  #
  # The z-score matrix is used for clustering so abundant cell
  # types do not dominate solely because of their magnitude.
  # ----------------------------------------------------------

  row_order <- rownames(mean_matrix)
  column_order <- colnames(mean_matrix)

  clustering_matrix <- z_matrix
  clustering_matrix[
    !is.finite(clustering_matrix)
  ] <- 0

  if (
    cluster_rows &&
    nrow(clustering_matrix) > 1
  ) {
    row_hclust <- stats::hclust(
      stats::dist(clustering_matrix),
      method = "complete"
    )

    row_order <- rownames(clustering_matrix)[
      row_hclust$order
    ]
  }

  if (
    cluster_columns &&
    ncol(clustering_matrix) > 1
  ) {
    column_hclust <- stats::hclust(
      stats::dist(
        t(clustering_matrix)
      ),
      method = "complete"
    )

    column_order <- colnames(clustering_matrix)[
      column_hclust$order
    ]
  }

  # ----------------------------------------------------------
  # Prepare plotting data
  # ----------------------------------------------------------

  raw_plot_df <- mean_abundance |>
    dplyr::mutate(
      celltype = factor(
        .data$celltype,
        levels = rev(row_order)
      ),
      method = factor(
        .data$method,
        levels = column_order
      )
    )

  z_plot_df <- as.data.frame(
    as.table(z_matrix),
    stringsAsFactors = FALSE
  )

  names(z_plot_df) <- c(
    "celltype",
    "method",
    "z_score"
  )

  z_plot_df <- z_plot_df |>
    dplyr::mutate(
      celltype = factor(
        .data$celltype,
        levels = rev(row_order)
      ),
      method = factor(
        .data$method,
        levels = column_order
      )
    )

  # ----------------------------------------------------------
  # Panel A: raw mean proportions
  #
  # The square-root transformation affects only the display
  # scale; the underlying values remain raw mean proportions.
  # ----------------------------------------------------------

  p_raw <- ggplot2::ggplot(
    raw_plot_df,
    ggplot2::aes(
      x = .data$method,
      y = .data$celltype,
      fill = .data$mean_proportion
    )
  ) +
    ggplot2::geom_tile(
      linewidth = 0.25,
      color = "white"
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      trans = "sqrt",
      na.value = "grey90",
      name = "Mean\nproportion"
    ) +
    ggplot2::labs(
      title = "A. Mean estimated proportion",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(
      base_size = 9
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = ggplot2::element_text(
        size = 6.5
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 11
      ),
      legend.position = "right"
    )

  # ----------------------------------------------------------
  # Panel B: row-scaled method bias
  # ----------------------------------------------------------

  finite_z <- z_plot_df$z_score[
    is.finite(z_plot_df$z_score)
  ]

  z_limit <- if (length(finite_z) == 0) {
    2
  } else {
    max(
      abs(finite_z),
      na.rm = TRUE
    )
  }

  if (
    !is.finite(z_limit) ||
    z_limit == 0
  ) {
    z_limit <- 2
  }

  p_z <- ggplot2::ggplot(
    z_plot_df,
    ggplot2::aes(
      x = .data$method,
      y = .data$celltype,
      fill = .data$z_score
    )
  ) +
    ggplot2::geom_tile(
      linewidth = 0.25,
      color = "white"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(
        -z_limit,
        z_limit
      ),
      oob = scales::squish,
      na.value = "grey90",
      name = "Row\nz-score"
    ) +
    ggplot2::labs(
      title = "B. Relative method bias",
      subtitle = "Scaled within each cell type",
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(
      base_size = 9
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = ggplot2::element_text(
        size = 6.5
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 11
      ),
      plot.subtitle = ggplot2::element_text(
        size = 8
      ),
      legend.position = "right"
    )

  # ----------------------------------------------------------
  # Draw the two ggplots side by side using base grid
  #
  # This avoids adding patchwork or cowplot as dependencies.
  # ----------------------------------------------------------

  draw_two_panels <- function() {

    grid::grid.newpage()

    figure_layout <- grid::grid.layout(
      nrow = 1,
      ncol = 2,
      widths = grid::unit(
        c(1, 1),
        "null"
      )
    )

    grid::pushViewport(
      grid::viewport(
        layout = figure_layout
      )
    )

    print(
      p_raw,
      vp = grid::viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
      )
    )

    print(
      p_z,
      vp = grid::viewport(
        layout.pos.row = 1,
        layout.pos.col = 2
      )
    )

    grid::popViewport()
  }

  pdf_file <- file.path(
    output_dir,
    "method_celltype_mean_and_bias_heatmaps.pdf"
  )

  png_file <- file.path(
    output_dir,
    "method_celltype_mean_and_bias_heatmaps.png"
  )

  grDevices::pdf(
    pdf_file,
    width = 13,
    height = 10,
    useDingbats = FALSE
  )

  draw_two_panels()

  grDevices::dev.off()

  grDevices::png(
    png_file,
    width = 13,
    height = 10,
    units = "in",
    res = 300
  )

  draw_two_panels()

  grDevices::dev.off()

  message(
    "Wrote ",
    pdf_file
  )

  message(
    "Wrote ",
    png_file
  )

  invisible(
    list(
      raw_plot = p_raw,
      zscore_plot = p_z,
      mean_abundance = mean_abundance,
      mean_matrix = mean_matrix,
      z_matrix = z_matrix,
      row_order = row_order,
      method_order = column_order
    )
  )
}

