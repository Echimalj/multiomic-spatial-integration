#' Figures for the pathway-proportion linkage analysis
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2.
#'
#' @keywords internal
NULL


#' GSEA dot plot of top pathways for one cell type
#'
#' Reads an fgsea result table (`pathway`, `NES`, `padj`, `size`) and plots
#' the top pathways by adjusted p-value: NES on the x-axis, dot size by set
#' size, and dot color by NES.
#'
#' @param enrichment_csv Path to a `pathway_enrichment.csv`.
#' @param celltype_label Label used in the title and output filename.
#' @param output_dir Figure output directory.
#' @param top_n Number of pathways, ranked by adjusted p-value, to show.
#'
#' @return The ggplot object invisibly, or `NULL` if there is nothing to plot.
#'
#' @export
plot_gsea_dotplot <- function(
    enrichment_csv,
    celltype_label,
    output_dir,
    top_n = 15L
) {
  if (
    length(enrichment_csv) != 1L ||
    is.na(enrichment_csv) ||
    !file.exists(enrichment_csv)
  ) {
    stop(
      "plot_gsea_dotplot: enrichment file not found: ",
      enrichment_csv,
      call. = FALSE
    )
  }

  if (
    length(top_n) != 1L ||
    !is.numeric(top_n) ||
    !is.finite(top_n) ||
    top_n < 1 ||
    top_n != floor(top_n)
  ) {
    stop(
      "plot_gsea_dotplot: `top_n` must be one positive integer.",
      call. = FALSE
    )
  }

  df <- utils::read.csv(
    enrichment_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "pathway",
    "NES",
    "padj",
    "size"
  )

  missing_cols <- setdiff(
    required_cols,
    names(df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "plot_gsea_dotplot: enrichment table is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df <- df[
    !is.na(df$pathway) &
      nzchar(as.character(df$pathway)) &
      !is.na(df$NES) &
      is.finite(df$NES) &
      !is.na(df$padj) &
      is.finite(df$padj),
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0L) {
    message(
      "plot_gsea_dotplot: nothing to plot for ",
      celltype_label,
      "."
    )

    return(
      invisible(NULL)
    )
  }

  df <- df[
    order(
      df$padj,
      -abs(df$NES),
      df$pathway,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  df <- utils::head(
    df,
    top_n
  )

  df <- df[
    order(df$NES),
    ,
    drop = FALSE
  ]

  df$pathway <- factor(
    df$pathway,
    levels = unique(df$pathway)
  )

  nes_limit <- max(
    abs(df$NES),
    na.rm = TRUE
  )

  if (!is.finite(nes_limit) || nes_limit == 0) {
    nes_limit <- 1
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = NES,
      y = pathway
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = size,
        color = NES
      )
    ) +
    scale_color_diverging(
      name = "NES",
      limits = c(
        -nes_limit,
        nes_limit
      )
    ) +
    ggplot2::scale_size_continuous(
      name = "Set size"
    ) +
    ggplot2::labs(
      title = paste0(
        "Pathways covarying with ",
        celltype_label,
        " abundance"
      ),
      subtitle = paste0(
        "Top ",
        min(top_n, nrow(df)),
        " by adjusted p-value"
      ),
      x = "Normalized enrichment score",
      y = NULL
    ) +
    theme_analysis()

  save_figure(
    p,
    paste0(
      "gsea_",
      celltype_label
    ),
    output_dir,
    width = 9,
    height = 6
  )

  invisible(p)
}


#' Volcano plot of genes ranked by association with a cell type's proportion
#'
#' @param ranked_csv Path to a `ranked_genes.csv` containing `gene`,
#'   `statistic`, `p.value`, and `p_adj`.
#' @param celltype_label Label used in the title and output filename.
#' @param output_dir Figure output directory.
#' @param padj_cutoff Adjusted p-value threshold used for point coloring.
#' @param label_top Number of significant genes to label at each end of the
#'   association axis.
#'
#' @return The ggplot object invisibly.
#'
#' @export
plot_ranked_gene_volcano <- function(
    ranked_csv,
    celltype_label,
    output_dir,
    padj_cutoff = 0.05,
    label_top = 10L
) {
  if (
    length(ranked_csv) != 1L ||
    is.na(ranked_csv) ||
    !file.exists(ranked_csv)
  ) {
    stop(
      "plot_ranked_gene_volcano: ranked-gene file not found: ",
      ranked_csv,
      call. = FALSE
    )
  }

  if (
    length(padj_cutoff) != 1L ||
    !is.numeric(padj_cutoff) ||
    !is.finite(padj_cutoff) ||
    padj_cutoff <= 0 ||
    padj_cutoff >= 1
  ) {
    stop(
      "plot_ranked_gene_volcano: `padj_cutoff` must be one finite ",
      "numeric value between 0 and 1.",
      call. = FALSE
    )
  }

  if (
    length(label_top) != 1L ||
    !is.numeric(label_top) ||
    !is.finite(label_top) ||
    label_top < 0 ||
    label_top != floor(label_top)
  ) {
    stop(
      "plot_ranked_gene_volcano: `label_top` must be one ",
      "non-negative integer.",
      call. = FALSE
    )
  }

  df <- utils::read.csv(
    ranked_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "gene",
    "statistic",
    "p.value",
    "p_adj"
  )

  missing_cols <- setdiff(
    required_cols,
    names(df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "plot_ranked_gene_volcano: ranked-gene table is missing ",
      "required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df <- df[
    !is.na(df$statistic) &
      is.finite(df$statistic) &
      !is.na(df$p.value) &
      is.finite(df$p.value) &
      df$p.value >= 0 &
      !is.na(df$p_adj) &
      is.finite(df$p_adj),
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0L) {
    message(
      "plot_ranked_gene_volcano: nothing to plot for ",
      celltype_label,
      "."
    )

    return(
      invisible(NULL)
    )
  }

  positive_p <- df$p.value[df$p.value > 0]

  minimum_positive_p <- if (length(positive_p) > 0L) {
    min(positive_p)
  } else {
    .Machine$double.xmin
  }

  plotting_p <- pmax(
    df$p.value,
    minimum_positive_p
  )

  df$neglog_p <- -log10(plotting_p)

  df$direction <- ifelse(
    df$p_adj >= padj_cutoff,
    "n.s.",
    ifelse(
      df$statistic > 0,
      "pos. assoc.",
      "neg. assoc."
    )
  )

  significant <- df[
    df$p_adj < padj_cutoff,
    ,
    drop = FALSE
  ]

  significant <- significant[
    order(
      -significant$neglog_p,
      -abs(significant$statistic),
      significant$gene,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  top_positive <- utils::head(
    significant[
      significant$statistic > 0,
      ,
      drop = FALSE
    ],
    label_top
  )

  top_negative <- utils::head(
    significant[
      significant$statistic < 0,
      ,
      drop = FALSE
    ],
    label_top
  )

  to_label <- rbind(
    top_positive,
    top_negative
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = statistic,
      y = neglog_p,
      color = direction
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.5,
      size = 1
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "neg. assoc." = "#2166AC",
        "n.s." = "grey75",
        "pos. assoc." = "#B2182B"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = paste0(
        "Genes associated with ",
        celltype_label,
        " abundance"
      ),
      x = "Spearman correlation with proportion",
      y = expression(-log[10]~"(p)")
    ) +
    theme_analysis()

  if (
    requireNamespace(
      "ggrepel",
      quietly = TRUE
    ) &&
      nrow(to_label) > 0L
  ) {
    p <- p +
      ggrepel::geom_text_repel(
        data = to_label,
        ggplot2::aes(label = gene),
        size = 2.5,
        max.overlaps = 20,
        show.legend = FALSE
      )
  }

  save_figure(
    p,
    paste0(
      "volcano_",
      celltype_label
    ),
    output_dir,
    width = 8,
    height = 6
  )

  invisible(p)
}


#' Plot low-recurrence pathway associations
#'
#' Plots pathway–cell-type associations selected by
#' `top_low_recurrence_associations_by_stratum()`.
#'
#' Each point represents one pathway–cell-type association. The x-axis shows
#' normalized enrichment score, point color shows NES, and point size shows
#' adjusted significance as `-log10(padj)`.
#'
#' Pathways shown here have already been filtered based on their recurrence
#' across all tested cell-type–stratum combinations and the curated broad
#' pathway-category flag. Low recurrence should be interpreted as evidence
#' of context restriction, not necessarily as evidence of greater biological
#' credibility.
#'
#' @param pathway_df Output from
#'   `top_low_recurrence_associations_by_stratum()`. Expected columns are
#'   `stratum`, `celltype`, `pathway`, `padj`, `NES`, and `n_significant`.
#' @param output_dir Figure output directory.
#' @param label_width Maximum number of pathway-name characters retained
#'   before truncation.
#' @param combined If `TRUE`, also write a combined figure faceted by
#'   stratum.
#'
#' @return A named list of ggplot objects, one per stratum and optionally one
#'   named `"combined"`, invisibly. Returns `NULL` when there is nothing to
#'   plot.
#'
#' @export
plot_pathway_recurrence_dotplot <- function(
    pathway_df,
    output_dir,
    label_width = 55L,
    combined = FALSE
) {
  if (!is.data.frame(pathway_df)) {
    stop(
      "plot_pathway_recurrence_dotplot: `pathway_df` must be a ",
      "data frame.",
      call. = FALSE
    )
  }

  if (
    length(label_width) != 1L ||
    !is.numeric(label_width) ||
    !is.finite(label_width) ||
    label_width < 1 ||
    label_width != floor(label_width)
  ) {
    stop(
      "plot_pathway_recurrence_dotplot: `label_width` must be one ",
      "positive integer.",
      call. = FALSE
    )
  }

  if (
    length(combined) != 1L ||
    !is.logical(combined) ||
    is.na(combined)
  ) {
    stop(
      "plot_pathway_recurrence_dotplot: `combined` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  required_cols <- c(
    "stratum",
    "celltype",
    "pathway",
    "padj",
    "NES",
    "n_significant"
  )

  missing_cols <- setdiff(
    required_cols,
    names(pathway_df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "plot_pathway_recurrence_dotplot: input is missing required ",
      "column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df <- pathway_df[
    !is.na(pathway_df$stratum) &
      nzchar(as.character(pathway_df$stratum)) &
      !is.na(pathway_df$celltype) &
      nzchar(as.character(pathway_df$celltype)) &
      !is.na(pathway_df$pathway) &
      nzchar(as.character(pathway_df$pathway)) &
      !is.na(pathway_df$NES) &
      is.finite(pathway_df$NES) &
      !is.na(pathway_df$padj) &
      is.finite(pathway_df$padj) &
      pathway_df$padj >= 0,
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0L) {
    message(
      "plot_pathway_recurrence_dotplot: nothing to plot."
    )

    return(
      invisible(NULL)
    )
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  pathway_text <- as.character(df$pathway)

  df$pathway_short <- ifelse(
    nchar(pathway_text) > label_width,
    paste0(
      substr(
        pathway_text,
        1L,
        label_width
      ),
      "..."
    ),
    pathway_text
  )

  df$label <- paste0(
    df$pathway_short,
    "  [",
    df$celltype,
    "]"
  )

  positive_padj <- df$padj[df$padj > 0]

  minimum_positive_padj <- if (length(positive_padj) > 0L) {
    min(positive_padj)
  } else {
    .Machine$double.xmin
  }

  df$neglog_padj <- -log10(
    pmax(
      df$padj,
      minimum_positive_padj
    )
  )

  nes_limit <- max(
    abs(df$NES),
    na.rm = TRUE
  )

  if (!is.finite(nes_limit) || nes_limit == 0) {
    nes_limit <- 1
  }

  make_dotplot <- function(
      plot_df,
      title,
      subtitle
  ) {
    plot_df <- plot_df[
      order(
        plot_df$NES,
        plot_df$padj,
        plot_df$label,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    plot_df$label <- factor(
      plot_df$label,
      levels = unique(plot_df$label)
    )

    ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = NES,
        y = label
      )
    ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "grey60"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          size = neglog_padj,
          color = NES
        )
      ) +
      scale_color_diverging(
        name = "NES",
        limits = c(
          -nes_limit,
          nes_limit
        )
      ) +
      ggplot2::scale_size_continuous(
        name = expression(-log[10]~"(padj)")
      ) +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        x = "Normalized enrichment score",
        y = NULL
      ) +
      theme_analysis(
        base_size = 9
      )
  }

  results <- list()

  stratum_names <- unique(
    as.character(df$stratum)
  )

  for (stratum_name in stratum_names) {
    stratum_df <- df[
      df$stratum == stratum_name,
      ,
      drop = FALSE
    ]

    p <- make_dotplot(
      stratum_df,
      title = paste0(
        "Low-recurrence pathway associations (",
        stratum_name,
        ")"
      ),
      subtitle = "pathway [associated cell type]"
    )

    safe_stratum_name <- gsub(
      "[^[:alnum:]_-]+",
      "_",
      stratum_name
    )

    save_figure(
      p,
      paste0(
        "pathway_recurrence_dotplot_",
        safe_stratum_name
      ),
      output_dir,
      width = 9,
      height = max(
        4,
        0.35 * nrow(stratum_df)
      )
    )

    results[[stratum_name]] <- p
  }

  if (combined) {
    p_all <- make_dotplot(
      df,
      title = "Low-recurrence pathway associations",
      subtitle = "pathway [associated cell type]"
    ) +
      ggplot2::facet_wrap(
        ~stratum,
        scales = "free_y",
        ncol = 1
      )

    n_strata <- length(
      stratum_names
    )

    save_figure(
      p_all,
      "pathway_recurrence_dotplot_combined",
      output_dir,
      width = 10,
      height = max(
        8,
        1.8 * n_strata,
        0.25 * nrow(df)
      )
    )

    results[["combined"]] <- p_all
  }

  invisible(results)
}


#' Plot cell-type-restricted pathway hits
#'
#' Backward-compatible wrapper around
#' `plot_pathway_recurrence_dotplot()`.
#'
#' This function preserves the previous API used by the pathway-specificity
#' workflow. New code should use `plot_pathway_recurrence_dotplot()`.
#'
#' @param top_specific_df Pathway-association data frame.
#' @param output_dir Figure output directory.
#' @param label_width Maximum number of pathway-name characters retained
#'   before truncation.
#' @param combined If `TRUE`, also write a combined figure faceted by
#'   stratum.
#'
#' @return The result of `plot_pathway_recurrence_dotplot()`, invisibly.
#'
#' @seealso plot_pathway_recurrence_dotplot
#'
#' @export
plot_specificity_dotplot <- function(
    top_specific_df,
    output_dir,
    label_width = 55L,
    combined = FALSE
) {
  plot_pathway_recurrence_dotplot(
    pathway_df = top_specific_df,
    output_dir = output_dir,
    label_width = label_width,
    combined = combined
  )
}
