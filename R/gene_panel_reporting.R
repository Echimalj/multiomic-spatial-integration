# ============================================================
# Gene-panel reporting utilities
# ============================================================

#' Format a gene-panel label for reports
#'
#' @param panel_name Internal panel identifier.
#'
#' @return Human-readable panel label.
#' @export
format_gene_panel_report_label <- function(panel_name) {
  known_labels <- c(
    Hajjar_2026 = "Hajjar 2026 biomarker panel",
    AGORA_High_Nomination = "AGORA high-nomination targets"
  )

  if (
    length(panel_name) == 1L &&
      panel_name %in% names(known_labels)
  ) {
    return(
      unname(known_labels[[panel_name]])
    )
  }

  gsub(
    "_",
    " ",
    panel_name,
    fixed = TRUE
  )
}


#' Rank genes using cross-module evidence
#'
#' The ranking score is intentionally transparent:
#'
#' - Spatial FDR-significant evidence: 3 points
#' - Spatial nominal-only evidence: 1 point
#' - Fine-cell-type global FDR evidence: 3 points
#' - Lineage global FDR evidence: 3 points
#' - Exploratory nominal lineage evidence: 1 point
#' - Panel-designated headline target: 1 point
#'
#' Reference attribution is reported but is not scored because expression
#' in a reference lineage is annotation evidence rather than perturbation
#' evidence.
#'
#' @param integrated_summary Data frame from
#'   `integrated_gene_panel_summary.csv`.
#'
#' @return Ranked gene table.
#' @export
rank_gene_panel_priorities <- function(
    integrated_summary
) {
  required_cols <- c(
    "gene",
    "panel",
    "category",
    "headline",
    "spatial_evidence",
    "lineage_evidence",
    "fine_celltype_evidence"
  )

  missing_cols <- setdiff(
    required_cols,
    names(integrated_summary)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "rank_gene_panel_priorities: missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  ranked <- integrated_summary |>
    dplyr::mutate(
      headline = dplyr::coalesce(
        as.logical(headline),
        FALSE
      ),

      spatial_score = dplyr::case_when(
        spatial_evidence == "FDR-significant" ~ 3L,
        spatial_evidence == "Nominal only" ~ 1L,
        TRUE ~ 0L
      ),

      lineage_score = dplyr::case_when(
        lineage_evidence == "Global FDR-significant" ~ 3L,
        lineage_evidence == "Exploratory nominal" ~ 1L,
        TRUE ~ 0L
      ),

      fine_celltype_score = dplyr::case_when(
        fine_celltype_evidence == "Global FDR-significant" ~ 3L,
        TRUE ~ 0L
      ),

      headline_score = ifelse(
        headline,
        1L,
        0L
      ),

      evidence_score =
        spatial_score +
        lineage_score +
        fine_celltype_score +
        headline_score,

      n_supported_modules =
        as.integer(spatial_score > 0L) +
        as.integer(lineage_score > 0L) +
        as.integer(fine_celltype_score > 0L),

      spatial_strength = dplyr::if_else(
        is.finite(strongest_spatial_fdr) &
          strongest_spatial_fdr > 0,
        -log10(strongest_spatial_fdr),
        NA_real_
      ),

      fine_strength = pmax(
        dplyr::if_else(
          is.finite(p_adj_global_AmyloidFree) &
            p_adj_global_AmyloidFree > 0,
          -log10(p_adj_global_AmyloidFree),
          NA_real_
        ),
        dplyr::if_else(
          is.finite(p_adj_global_Amyloid) &
            p_adj_global_Amyloid > 0,
          -log10(p_adj_global_Amyloid),
          NA_real_
        ),
        na.rm = TRUE
      )
    )

  ranked$fine_strength[
    !is.finite(ranked$fine_strength)
  ] <- NA_real_

  ranked <- ranked |>
    dplyr::arrange(
      dplyr::desc(evidence_score),
      dplyr::desc(n_supported_modules),
      dplyr::desc(headline),
      dplyr::desc(spatial_strength),
      dplyr::desc(fine_strength),
      gene
    ) |>
    dplyr::mutate(
      evidence_rank = dplyr::row_number(),
      evidence_tier = dplyr::case_when(
        evidence_score >= 7L ~ "Tier A: convergent evidence",
        evidence_score >= 5L ~ "Tier B: strong evidence",
        evidence_score >= 3L ~ "Tier C: single-module FDR evidence",
        evidence_score >= 1L ~ "Tier D: exploratory evidence",
        TRUE ~ "Tier E: no detected evidence"
      )
    ) |>
    dplyr::relocate(
      evidence_rank,
      evidence_tier,
      evidence_score,
      n_supported_modules,
      spatial_score,
      lineage_score,
      fine_celltype_score,
      headline_score,
      .after = source
    )

  ranked
}


#' Summarize one analyzed gene panel
#'
#' @param integrated_summary Integrated per-gene table.
#' @param reference_missing Optional missing-reference table.
#'
#' @return One-row panel summary.
#' @export
summarize_gene_panel <- function(
    integrated_summary,
    reference_missing = NULL
) {
  panel_names <- unique(
    integrated_summary$panel
  )

  panel_names <- panel_names[
    !is.na(panel_names)
  ]

  if (length(panel_names) != 1L) {
    stop(
      "summarize_gene_panel: expected exactly one panel label.",
      call. = FALSE
    )
  }

  panel_name <- panel_names[[1]]

  n_missing_reference <- if (
    is.null(reference_missing) ||
      nrow(reference_missing) == 0L
  ) {
    0L
  } else {
    dplyr::n_distinct(
      reference_missing$gene
    )
  }

  data.frame(
    panel = panel_name,
    panel_label = format_gene_panel_report_label(
      panel_name
    ),
    n_panel_genes = nrow(
      integrated_summary
    ),
    n_missing_reference_genes = n_missing_reference,
    n_spatial_fdr_genes = sum(
      integrated_summary$spatial_evidence ==
        "FDR-significant",
      na.rm = TRUE
    ),
    n_spatial_nominal_only_genes = sum(
      integrated_summary$spatial_evidence ==
        "Nominal only",
      na.rm = TRUE
    ),
    n_lineage_global_fdr_genes = sum(
      integrated_summary$lineage_evidence ==
        "Global FDR-significant",
      na.rm = TRUE
    ),
    n_lineage_exploratory_genes = sum(
      integrated_summary$lineage_evidence ==
        "Exploratory nominal",
      na.rm = TRUE
    ),
    n_fine_celltype_fdr_genes = sum(
      integrated_summary$fine_celltype_evidence ==
        "Global FDR-significant",
      na.rm = TRUE
    ),
    n_headline_genes = sum(
      integrated_summary$headline,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}


#' Write a Markdown report for a gene panel
#'
#' @param integrated_summary Integrated summary table.
#' @param ranked_genes Ranked gene table.
#' @param panel_summary One-row panel summary.
#' @param panel_result_dir Panel-specific result directory.
#' @param panel_figure_dir Panel-specific figure directory.
#' @param output_file Markdown output path.
#' @param top_n Number of ranked genes displayed.
#'
#' @return Output path invisibly.
#' @export
write_gene_panel_markdown_report <- function(
    integrated_summary,
    ranked_genes,
    panel_summary,
    panel_result_dir,
    panel_figure_dir,
    output_file,
    top_n = 15L
) {
  panel_name <- panel_summary$panel[[1]]
  panel_label <- panel_summary$panel_label[[1]]

  top_genes <- utils::head(
    ranked_genes,
    top_n
  )

  markdown_escape <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    gsub("\\|", "\\\\|", x)
  }

  format_number <- function(x, digits = 3L) {
    ifelse(
      is.finite(x),
      formatC(
        x,
        digits = digits,
        format = "fg",
        flag = "#"
      ),
      ""
    )
  }

  top_table <- data.frame(
    Rank = top_genes$evidence_rank,
    Gene = top_genes$gene,
    Tier = top_genes$evidence_tier,
    Score = top_genes$evidence_score,
    Reference_lineage =
      top_genes$reference_dominant_lineage,
    Spatial_region =
      top_genes$strongest_spatial_region,
    Spatial_contrast =
      top_genes$strongest_spatial_contrast,
    Spatial_effect = format_number(
      top_genes$strongest_spatial_estimate
    ),
    Spatial_FDR = format_number(
      top_genes$strongest_spatial_fdr
    ),
    Fine_celltype_AmyloidFree =
      top_genes$celltype_AmyloidFree,
    Fine_celltype_Amyloid =
      top_genes$celltype_Amyloid,
    stringsAsFactors = FALSE
  )

  table_lines <- c(
    paste(
      names(top_table),
      collapse = " | "
    ),
    paste(
      rep("---", ncol(top_table)),
      collapse = " | "
    )
  )

  for (i in seq_len(nrow(top_table))) {
    table_lines <- c(
      table_lines,
      paste(
        markdown_escape(
          unlist(
            top_table[i, ],
            use.names = FALSE
          )
        ),
        collapse = " | "
      )
    )
  }

  figure_files <- c(
    "gene_panel_reference_lineage_heatmap.png",
    "gene_panel_spatial_perturbation_heatmap.png",
    "gene_panel_lineage_attribution_heatmap.png",
    "gene_panel_lineage_attribution_heatmap_exploratory_nominal.png",
    "gene_panel_fine_celltype_correlation_heatmap.png"
  )

  existing_figures <- figure_files[
    file.exists(
      file.path(
        panel_figure_dir,
        figure_files
      )
    )
  ]

  relative_figure_paths <- file.path(
    "..",
    "..",
    "figures",
    "gene_panel_validation",
    panel_name,
    existing_figures
  )

  figure_labels <- c(
    gene_panel_reference_lineage_heatmap.png =
      "Reference lineage attribution",
    gene_panel_spatial_perturbation_heatmap.png =
      "Spatial perturbation",
    gene_panel_lineage_attribution_heatmap.png =
      "Lineage-associated perturbation",
    gene_panel_lineage_attribution_heatmap_exploratory_nominal.png =
      "Lineage-associated perturbation — exploratory nominal display",
    gene_panel_fine_celltype_correlation_heatmap.png =
      "Fine-cell-type discovery correlations"
  )

  figure_section <- character(0)

  for (i in seq_along(existing_figures)) {
    filename <- existing_figures[[i]]

    figure_section <- c(
      figure_section,
      paste0(
        "### ",
        figure_labels[[filename]]
      ),
      "",
      paste0(
        "![",
        figure_labels[[filename]],
        "](",
        relative_figure_paths[[i]],
        ")"
      ),
      ""
    )
  }

  spatial_sentence <- paste0(
    panel_summary$n_spatial_fdr_genes,
    " of ",
    panel_summary$n_panel_genes,
    " genes had at least one FDR-significant spatial perturbation."
  )

  lineage_sentence <- if (
    panel_summary$n_lineage_global_fdr_genes > 0L
  ) {
    paste0(
      panel_summary$n_lineage_global_fdr_genes,
      " genes had a globally FDR-significant lineage interaction."
    )
  } else {
    paste0(
      "No lineage interaction passed global FDR; ",
      panel_summary$n_lineage_exploratory_genes,
      " genes had nominal exploratory lineage evidence."
    )
  }

  fine_sentence <- paste0(
    panel_summary$n_fine_celltype_fdr_genes,
    " of ",
    panel_summary$n_panel_genes,
    " genes had at least one globally FDR-significant ",
    "fine-cell-type association."
  )

  report_lines <- c(
    paste0("# ", panel_label),
    "",
    paste0(
      "_Generated: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M %Z"
      ),
      "_"
    ),
    "",
    "## Panel overview",
    "",
    paste0(
      "- Panel identifier: `",
      panel_name,
      "`"
    ),
    paste0(
      "- Genes analyzed: ",
      panel_summary$n_panel_genes
    ),
    paste0(
      "- Genes missing from the reference signatures: ",
      panel_summary$n_missing_reference_genes
    ),
    paste0(
      "- Headline genes: ",
      panel_summary$n_headline_genes
    ),
    "",
    "## Evidence summary",
    "",
    paste0("- ", spatial_sentence),
    paste0("- ", lineage_sentence),
    paste0("- ", fine_sentence),
    "",
    "## Ranking framework",
    "",
    "The evidence score is a transparent prioritization tool:",
    "",
    "- Spatial FDR evidence: 3 points",
    "- Spatial nominal-only evidence: 1 point",
    "- Fine-cell-type global FDR evidence: 3 points",
    "- Lineage global FDR evidence: 3 points",
    "- Exploratory nominal lineage evidence: 1 point",
    "- Panel-designated headline target: 1 point",
    "",
    "The score is intended for prioritization and does not represent an independent statistical test.",
    "",
    paste0(
      "## Top ",
      min(top_n, nrow(top_genes)),
      " ranked genes"
    ),
    "",
    table_lines,
    "",
    "## Figures",
    "",
    figure_section,
    "## Output files",
    "",
    paste0(
      "- Integrated summary: `",
      file.path(
        panel_result_dir,
        "integrated_gene_panel_summary.csv"
      ),
      "`"
    ),
    paste0(
      "- Ranked priorities: `",
      file.path(
        panel_result_dir,
        "ranked_gene_priorities.csv"
      ),
      "`"
    ),
    "",
    "## Interpretation note",
    "",
    paste0(
      "Reference attribution indicates the lineage in which a gene is ",
      "preferentially expressed. Spatial perturbation tests region-specific ",
      "expression changes. Lineage interactions assess whether the ",
      "expression–lineage association changes with local amyloid pathology. ",
      "Fine-cell-type correlations are scan-aware associations and should ",
      "not be interpreted as direct cellular causation."
    )
  )

  dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    report_lines,
    output_file
  )

  message(
    "Wrote ",
    output_file
  )

  invisible(output_file)
}
