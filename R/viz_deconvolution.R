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
