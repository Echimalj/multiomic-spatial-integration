# ============================================================
# Gene-panel visualization
# ============================================================

#' Plot panel-gene reference expression across broad lineages
#'
#' Produces a condition-faceted heatmap of reference signature expression.
#' Values are z-scaled across lineages separately for each gene and condition,
#' emphasizing the relative lineage preference of each gene rather than
#' absolute differences in expression scale among genes.
#'
#' @param reference_csv Path to
#'   `reference_expression_by_lineage_all_conditions.csv`.
#' @param output_dir Figure output directory.
#' @param headline_first Place genes marked as headline genes first.
#'
#' @return The ggplot object invisibly.
#' @export
plot_gene_panel_reference_heatmap <- function(
    reference_csv,
    output_dir,
    headline_first = TRUE
) {
  if (!file.exists(reference_csv)) {
    stop(
      "plot_gene_panel_reference_heatmap: file does not exist: ",
      reference_csv,
      call. = FALSE
    )
  }

  df <- utils::read.csv(
    reference_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "gene",
    "lineage",
    "mean_expression",
    "condition",
    "category",
    "headline"
  )

  missing_cols <- setdiff(
    required_cols,
    names(df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "plot_gene_panel_reference_heatmap: missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df$mean_expression <- suppressWarnings(
    as.numeric(df$mean_expression)
  )

  df <- df[
    is.finite(df$mean_expression) &
      !is.na(df$gene) &
      nzchar(df$gene) &
      !is.na(df$lineage) &
      nzchar(df$lineage),
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0L) {
    message(
      "plot_gene_panel_reference_heatmap: no finite expression rows."
    )

    return(invisible(NULL))
  }

  # Z-scale across lineages within each gene and condition.
  plot_df <- df |>
    dplyr::group_by(
      condition,
      gene
    ) |>
    dplyr::mutate(
      reference_z = {
        x <- mean_expression
        finite_x <- is.finite(x)

        if (
          sum(finite_x) < 2L ||
          stats::sd(
            x[finite_x],
            na.rm = TRUE
          ) == 0
        ) {
          rep(
            0,
            dplyr::n()
          )
        } else {
          as.numeric(
            scale(x)
          )
        }
      }
    ) |>
    dplyr::ungroup()

  gene_metadata <- plot_df |>
    dplyr::distinct(
      gene,
      category,
      headline
    )

  if (headline_first) {
    gene_order <- gene_metadata |>
      dplyr::arrange(
        dplyr::desc(headline),
        category,
        gene
      ) |>
      dplyr::pull(gene)
  } else {
    gene_order <- gene_metadata |>
      dplyr::arrange(
        category,
        gene
      ) |>
      dplyr::pull(gene)
  }

  lineage_order <- c(
    "Vascular",
    "Astrocyte",
    "Microglia",
    "OPC",
    "Oligodendrocyte",
    "ExcitatoryNeuron",
    "InhibitoryNeuron",
    "Other"
  )

  lineage_order <- lineage_order[
    lineage_order %in%
      unique(plot_df$lineage)
  ]

  remaining_lineages <- setdiff(
    sort(unique(plot_df$lineage)),
    lineage_order
  )

  lineage_order <- c(
    lineage_order,
    remaining_lineages
  )

  condition_order <- c(
    "Control",
    "AD+CAA"
  )

  condition_order <- condition_order[
    condition_order %in%
      unique(plot_df$condition)
  ]

  plot_df <- plot_df |>
    dplyr::mutate(
      gene = factor(
        gene,
        levels = rev(gene_order)
      ),
      lineage = factor(
        lineage,
        levels = lineage_order
      ),
      condition = factor(
        condition,
        levels = condition_order
      )
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = lineage,
      y = gene,
      fill = reference_z
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.15
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(category),
      cols = ggplot2::vars(condition),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Relative\nexpression\nz-score"
    ) +
    ggplot2::labs(
      title = "Reference lineage expression of gene-panel targets",
      subtitle = paste0(
        "Expression is z-scaled across lineages within each gene and ",
        "condition; positive values indicate relative lineage enrichment"
      ),
      x = NULL,
      y = NULL,
      caption = paste0(
        "Lineage values represent the mean inferred reference expression ",
        "across fine cell types assigned to each broad lineage."
      )
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
        size = 7
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 8
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        colour = "grey80"
      ),
      legend.position = "right",
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  save_figure(
    p,
    "gene_panel_reference_lineage_heatmap",
    output_dir,
    width = 12,
    height = 15
  )

  invisible(p)
}

#' Plot spatial perturbation effects for a gene panel
#'
#' @param spatial_csv Path to `spatial_perturbation_by_region.csv`.
#' @param output_dir Figure output directory.
#' @param padj_cutoff Adjusted significance threshold.
#'
#' @return The ggplot object invisibly.
#' @export
plot_gene_panel_spatial_perturbation <- function(
    spatial_csv,
    output_dir,
    padj_cutoff = 0.05
) {
  if (!file.exists(spatial_csv)) {
    stop(
      "plot_gene_panel_spatial_perturbation: file does not exist: ",
      spatial_csv,
      call. = FALSE
    )
  }

  df <- utils::read.csv(
    spatial_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

panel_label <- if (
  "panel" %in% names(df) &&
    length(unique(stats::na.omit(df$panel))) == 1L
) {
  unique(stats::na.omit(df$panel))[[1]]
} else {
  "Gene panel"
}

panel_title <- gsub(
  "_",
  " ",
  panel_label,
  fixed = TRUE
)

  required_cols <- c(
    "gene",
    "category",
    "headline",
    "region",
    "contrast",
    "estimate",
    "p_adj",
    "model_status"
  )

  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0L) {
    stop(
      "plot_gene_panel_spatial_perturbation: missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  plot_df <- df |>
    dplyr::filter(
      model_status == "ok",
      is.finite(estimate)
    ) |>
    dplyr::mutate(
      significant = is.finite(p_adj) & p_adj < padj_cutoff,
      significance_symbol = dplyr::case_when(
        is.finite(p_adj) & p_adj < 0.001 ~ "***",
        is.finite(p_adj) & p_adj < 0.01  ~ "**",
        is.finite(p_adj) & p_adj < 0.05  ~ "*",
        TRUE ~ ""
      )
    )

  gene_order <- plot_df |>
    dplyr::group_by(gene, category, headline) |>
    dplyr::summarise(
      max_abs_effect = max(abs(estimate), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(headline),
      category,
      dplyr::desc(max_abs_effect),
      gene
    ) |>
    dplyr::pull(gene)

  plot_df <- plot_df |>
    dplyr::mutate(
      gene = factor(gene, levels = rev(gene_order)),
      region = factor(
        region,
        levels = c(
          "Arteries",
          "Capillaries",
          "Parenchyma"
        )
      ),
      contrast = factor(
        contrast,
        levels = c(
          "Disease_effect",
          "Amyloid_effect",
          "MaxPathology_effect"
        ),
        labels = c(
          "Disease",
          "Amyloid",
          "Maximum pathology"
        )
      )
    )

  fill_limit <- max(
    abs(plot_df$estimate),
    na.rm = TRUE
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = contrast,
      y = gene,
      fill = estimate
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = significance_symbol
      ),
      size = 2.5
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(category),
      cols = ggplot2::vars(region),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      name = "Expression\neffect"
    ) +
    ggplot2::labs(
      title = paste0( "Spatial perturbation of ", panel_title, " genes"),
      subtitle = paste0(
        "Gaussian mixed models: expression ~ group × region + ",
        "(1 | Scan_ID)"
      ),
      x = NULL,
      y = NULL,
      caption = paste0(
        "* adjusted p < 0.05; ** < 0.01; *** < 0.001. ",
        "BH correction was applied across genes within each ",
        "region × contrast."
      )
    ) +
    ggplot2::theme_minimal(
      base_size = 9
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(
        size = 7
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 8
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        colour = "grey80"
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )

  save_figure(
    p,
    "gene_panel_spatial_perturbation_heatmap",
    output_dir,
    width = 13,
    height = 15
  )

  invisible(p)
}

#' Plot lineage-associated perturbation for a gene panel
#'
#' Uses globally adjusted p-values when at least one interaction passes the
#' requested FDR threshold. If no interaction passes global FDR, nominal
#' p-values are used for exploratory display and the figure caption explicitly
#' states that the displayed significance symbols are unadjusted.
#'
#' @param lineage_csv Path to `lineage_attributed_perturbation.csv`.
#' @param output_dir Figure output directory.
#' @param fdr_cutoff Global adjusted-p threshold.
#' @param nominal_cutoff Nominal p-value threshold used only when no global
#'   FDR-significant interactions are present.
#'
#' @return The ggplot object invisibly.
#' @export
plot_gene_panel_lineage_attribution <- function(
    lineage_csv,
    output_dir,
    fdr_cutoff = 0.05,
    nominal_cutoff = 0.05
) {
  if (!file.exists(lineage_csv)) {
    stop(
      "plot_gene_panel_lineage_attribution: file does not exist: ",
      lineage_csv,
      call. = FALSE
    )
  }

  df <- utils::read.csv(
    lineage_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "gene",
    "category",
    "headline",
    "lineage",
    "amyloid_interaction",
    "p.value",
    "p_adj_global",
    "model_status"
  )

  missing_cols <- setdiff(
    required_cols,
    names(df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "plot_gene_panel_lineage_attribution: missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  plot_df <- df |>
    dplyr::filter(
      model_status == "ok",
      is.finite(amyloid_interaction),
      is.finite(p.value)
    )

  if (nrow(plot_df) == 0L) {
    message(
      "plot_gene_panel_lineage_attribution: no evaluable models."
    )
    return(invisible(NULL))
  }

  has_fdr_significant <- any(
    is.finite(plot_df$p_adj_global) &
      plot_df$p_adj_global < fdr_cutoff
  )

  if (has_fdr_significant) {
    significance_mode <- "global_fdr"

    plot_df <- plot_df |>
      dplyr::mutate(
        significance_value = p_adj_global,
        significant = is.finite(p_adj_global) &
          p_adj_global < fdr_cutoff,
        significance_symbol = dplyr::case_when(
          is.finite(p_adj_global) & p_adj_global < 0.001 ~ "***",
          is.finite(p_adj_global) & p_adj_global < 0.01 ~ "**",
          is.finite(p_adj_global) & p_adj_global < 0.05 ~ "*",
          TRUE ~ ""
        )
      )

    subtitle_text <- paste0(
      "Symbols indicate globally FDR-adjusted significance across all ",
      "gene × lineage interaction tests"
    )

    caption_text <- paste0(
      "* global BH-adjusted p < 0.05; ** < 0.01; *** < 0.001. ",
      "The interaction estimates the pathology-associated change in the ",
      "expression–lineage relationship."
    )
  } else {
    significance_mode <- "nominal"

    plot_df <- plot_df |>
      dplyr::mutate(
        significance_value = p.value,
        significant = p.value < nominal_cutoff,
        significance_symbol = dplyr::case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01 ~ "**",
          p.value < 0.05 ~ "*",
          TRUE ~ ""
        )
      )

    subtitle_text <- paste0(
      "No interaction passed global FDR; symbols show nominal p-values ",
      "for exploratory interpretation"
    )

    caption_text <- paste0(
      "Exploratory display only: * nominal p < 0.05; ** < 0.01; ",
      "*** < 0.001. These symbols are not multiple-testing adjusted. ",
      "No gene × lineage interaction passed global BH FDR < ",
      fdr_cutoff,
      "."
    )
  }

  # Keep all genes that have at least one displayed association.
  displayed_genes <- plot_df |>
    dplyr::group_by(
      gene
    ) |>
    dplyr::summarise(
      any_displayed = any(significant),
      strongest_p = min(significance_value, na.rm = TRUE),
      maximum_effect = max(
        abs(amyloid_interaction),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::filter(
      any_displayed
    )

  if (nrow(displayed_genes) == 0L) {
    message(
      "plot_gene_panel_lineage_attribution: no associations passed the ",
      "selected display threshold."
    )
    return(invisible(NULL))
  }

  plot_df <- plot_df |>
    dplyr::filter(
      gene %in% displayed_genes$gene
    )

  gene_order <- plot_df |>
    dplyr::group_by(
      gene,
      category,
      headline
    ) |>
    dplyr::summarise(
      strongest_p = min(
        significance_value,
        na.rm = TRUE
      ),
      max_abs_effect = max(
        abs(amyloid_interaction),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(headline),
      category,
      strongest_p,
      dplyr::desc(max_abs_effect),
      gene
    ) |>
    dplyr::pull(
      gene
    )

  lineage_order <- c(
    "Vascular",
    "Astrocyte",
    "Microglia",
    "OPC",
    "Oligodendrocyte",
    "ExcitatoryNeuron",
    "InhibitoryNeuron"
  )

  lineage_order <- lineage_order[
    lineage_order %in%
      unique(plot_df$lineage)
  ]

  plot_df <- plot_df |>
    dplyr::mutate(
      gene = factor(
        gene,
        levels = rev(gene_order)
      ),
      lineage = factor(
        lineage,
        levels = lineage_order
      )
    )

  fill_limit <- max(
    abs(plot_df$amyloid_interaction),
    na.rm = TRUE
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = lineage,
      y = gene,
      fill = amyloid_interaction
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = significance_symbol
      ),
      size = 3
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(category),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(
        -fill_limit,
        fill_limit
      ),
      name = "Amyloid\ninteraction"
    ) +
    ggplot2::labs(
      title = "Lineage-associated perturbation of gene-panel expression",
      subtitle = subtitle_text,
      x = NULL,
      y = NULL,
      caption = caption_text
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
        size = 7
      ),
      strip.text.y = ggplot2::element_text(
        face = "bold",
        angle = 0
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        colour = "grey80"
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 7
      )
    )

  filename <- if (
    significance_mode == "global_fdr"
  ) {
    "gene_panel_lineage_attribution_heatmap"
  } else {
    "gene_panel_lineage_attribution_heatmap_exploratory_nominal"
  }

  figure_height <- max(
    6,
    0.28 * length(unique(plot_df$gene)) + 3
  )

  save_figure(
    p,
    filename,
    output_dir,
    width = 10,
    height = figure_height
  )

  message(
    "Lineage-attribution significance mode: ",
    if (significance_mode == "global_fdr") {
      "global FDR"
    } else {
      "nominal p-values because no test passed global FDR"
    },
    "."
  )

  invisible(p)
}

#' Plot fine-cell-type discovery correlations
#'
#' Displays scan-aware gene-expression versus cell-type-proportion
#' correlations separately in amyloid-free and amyloid-positive AD/CAA ROIs.
#'
#' Tile fill represents Spearman correlation. A black dot marks associations
#' passing global BH FDR across all gene x cell-type x stratum tests.
#'
#' @param correlation_csv Path to
#'   `fine_celltype_discovery_correlations.csv`.
#' @param output_dir Figure output directory.
#' @param fdr_cutoff Global FDR threshold.
#' @param show_only_evaluable Restrict to models with status `"ok"`.
#'
#' @return The ggplot object invisibly.
#' @export
plot_gene_panel_fine_celltype_correlations <- function(
    correlation_csv,
    output_dir,
    fdr_cutoff = 0.05,
    show_only_evaluable = TRUE
) {
  if (!file.exists(correlation_csv)) {
    stop(
      "plot_gene_panel_fine_celltype_correlations: file does not exist: ",
      correlation_csv,
      call. = FALSE
    )
  }

  df <- utils::read.csv(
    correlation_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "gene",
    "category",
    "headline",
    "celltype",
    "lineage",
    "stratum",
    "correlation",
    "p_adj_global",
    "model_status"
  )

  missing_cols <- setdiff(
    required_cols,
    names(df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "plot_gene_panel_fine_celltype_correlations: missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  plot_df <- df

  if (show_only_evaluable) {
    plot_df <- plot_df |>
      dplyr::filter(
        model_status == "ok"
      )
  }

  plot_df <- plot_df |>
    dplyr::filter(
      is.finite(correlation)
    ) |>
    dplyr::mutate(
      significant_global =
        is.finite(p_adj_global) &
        p_adj_global < fdr_cutoff
    )

  if (nrow(plot_df) == 0L) {
    message(
      "plot_gene_panel_fine_celltype_correlations: ",
      "no evaluable correlation rows."
    )

    return(invisible(NULL))
  }

  n_significant <- sum(
    plot_df$significant_global,
    na.rm = TRUE
  )

  # ----------------------------------------------------------
  # Order genes by panel category and strongest correlation
  # ----------------------------------------------------------

  gene_order <- plot_df |>
    dplyr::group_by(
      gene,
      category,
      headline
    ) |>
    dplyr::summarise(
      max_abs_correlation = max(
        abs(correlation),
        na.rm = TRUE
      ),
      min_global_fdr = if (
        all(!is.finite(p_adj_global))
      ) {
        NA_real_
      } else {
        min(
          p_adj_global,
          na.rm = TRUE
        )
      },
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(headline),
      category,
      min_global_fdr,
      dplyr::desc(max_abs_correlation),
      gene
    ) |>
    dplyr::pull(
      gene
    )

  # ----------------------------------------------------------
  # Order cell types within broad lineage
  # ----------------------------------------------------------

  lineage_order <- c(
    "Vascular",
    "Astrocyte",
    "Microglia",
    "OPC",
    "Oligodendrocyte",
    "ExcitatoryNeuron",
    "InhibitoryNeuron",
    "Other"
  )

  celltype_order <- plot_df |>
    dplyr::distinct(
      celltype,
      lineage
    ) |>
    dplyr::mutate(
      lineage = factor(
        lineage,
        levels = lineage_order
      )
    ) |>
    dplyr::arrange(
      lineage,
      celltype
    ) |>
    dplyr::pull(
      celltype
    )

  stratum_order <- c(
    "AmyloidFree",
    "Amyloid"
  )

  plot_df <- plot_df |>
    dplyr::mutate(
      gene = factor(
        gene,
        levels = gene_order
      ),
      celltype = factor(
        celltype,
        levels = rev(celltype_order)
      ),
      lineage = factor(
        lineage,
        levels = lineage_order
      ),
      stratum = factor(
        stratum,
        levels = stratum_order,
        labels = c(
          "Amyloid-free AD/CAA",
          "Amyloid-positive AD/CAA"
        )
      )
    )

  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = gene,
      y = celltype,
      fill = correlation
    )
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.1
    ) +
    ggplot2::geom_point(
      data = plot_df |>
        dplyr::filter(
          significant_global
        ),
      shape = 16,
      size = 0.75,
      colour = "black"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(lineage),
      cols = ggplot2::vars(stratum),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman\ncorrelation"
    ) +
    ggplot2::labs(
      title = "Fine-cell-type associations of gene-panel expression",
      subtitle = paste0(
        "AD/CAA ROIs stratified by local amyloid pathology; ",
        n_significant,
        " associations pass global BH FDR < ",
        fdr_cutoff
      ),
      x = NULL,
      y = NULL,
      caption = paste0(
        "Tile fill represents Spearman correlation between normalized ",
        "gene expression and estimated cell-type proportion. ",
        "Black dots indicate global BH-adjusted p < ",
        fdr_cutoff,
        ". P-values were obtained from scan-aware mixed models."
      )
    ) +
    ggplot2::theme_minimal(
      base_size = 8
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),

      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6
      ),

      axis.text.y = ggplot2::element_text(
        size = 6
      ),

      strip.text.x = ggplot2::element_text(
        face = "bold",
        size = 9
      ),

      strip.text.y = ggplot2::element_text(
        face = "bold",
        angle = 0,
        size = 7
      ),

      strip.background = ggplot2::element_rect(
        fill = "grey95",
        colour = "grey80"
      ),

      legend.position = "right",

      plot.title = ggplot2::element_text(
        face = "bold"
      ),

      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = 7
      )
    )

  save_figure(
    p,
    "gene_panel_fine_celltype_correlation_heatmap",
    output_dir,
    width = 18,
    height = 15
  )

  invisible(p)
}
