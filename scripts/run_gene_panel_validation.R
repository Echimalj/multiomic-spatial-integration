#!/usr/bin/env Rscript

# ============================================================
# Gene-panel validation
#
# Initial implementation:
#   Reference expression and cell-type attribution using
#   Control and AD+CAA inferred signature matrices.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source("R/gene_panel_utils.R")
source("R/pathway_proportion_utils.R")
source("R/association_utils.R")

# ============================================================
# Configuration
# ============================================================

panel_files <- c(
  "resources/gene_panels/hajjar_2026_endothelial_biomarkers.csv",
  "resources/gene_panels/agora_high_nomination_targets.csv"
)

signature_files <- c(
  Control = "results/regression_model/Control_inferred_signatures.csv",
  `AD+CAA` = "results/regression_model/AD+CAA_inferred_signatures.csv"
)
expression_csv <-
  "results/geomx_exports/CAA-AD_expression_wide.csv"

proportions_csv <-
  "results/cell_proportions/spatial_celltype_proportions_for_R.csv"

output_root <- "results/gene_panel_validation"

dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Helpers
# ============================================================

safe_panel_name <- function(x) {
  x <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    x
  )

  gsub(
    "^_+|_+$",
    "",
    x
  )
}

load_signature_matrix <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Signature file does not exist: ",
      path,
      call. = FALSE
    )
  }

  signature <- as.matrix(
    utils::read.csv(
      path,
      row.names = 1,
      check.names = FALSE
    )
  )

  storage.mode(signature) <- "numeric"

  signature
}

# ============================================================
# Load shared spatial inputs
# ============================================================

message("Loading ROI expression...")

expression_raw <- load_roi_expression(
  expression_csv
)

expression_normalized <- normalize_expression_cpm(
  expression_raw
)

message(
  "Expression: ",
  nrow(expression_normalized),
  " ROIs x ",
  ncol(expression_normalized),
  " genes."
)

message("Loading cell-type proportions and ROI metadata...")

proportions_long <- utils::read.csv(
  proportions_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

roi_meta <- build_roi_meta(
  proportions_long
)

message(
  "ROI metadata: ",
  nrow(roi_meta),
  " ROIs; ",
  length(unique(roi_meta$Scan_ID)),
  " scans."
)
# ============================================================
# Run each gene panel
# ============================================================

for (panel_file in panel_files) {

  message(
    "\n============================================================"
  )
  message("Loading gene panel: ", panel_file)
  message(
    "============================================================"
  )

  panel <- load_gene_panel(
    panel_file
  )

  panel_names <- unique(
    panel$panel
  )

  if (length(panel_names) != 1L) {
    stop(
      "Expected one panel label in ",
      panel_file,
      ", but found: ",
      paste(panel_names, collapse = ", "),
      call. = FALSE
    )
  }

  panel_name <- panel_names[[1]]

  panel_output_dir <- file.path(
    output_root,
    safe_panel_name(panel_name)
  )

  dir.create(
    panel_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  all_per_lineage <- list()
  all_summary <- list()
  all_missing <- list()

  for (condition_name in names(signature_files)) {

    signature_file <- signature_files[[condition_name]]

    message(
      "\nReference attribution: ",
      panel_name,
      " / ",
      condition_name
    )

    signature_mat <- load_signature_matrix(
      signature_file
    )

    attribution <- attribute_genes_to_celltypes(
      signature_mat = signature_mat,
      genes = panel$gene
    )

    per_lineage <- attribution$per_lineage |>
      dplyr::mutate(
        condition = condition_name,
        panel = panel_name,
        .before = 1
      ) |>
      dplyr::left_join(
        panel |>
          dplyr::select(
            gene,
            category,
            headline,
            source
          ),
        by = "gene"
      )

    summary_df <- attribution$summary |>
      dplyr::mutate(
        condition = condition_name,
        panel = panel_name,
        .before = 1
      ) |>
      dplyr::left_join(
        panel |>
          dplyr::select(
            gene,
            category,
            headline,
            source
          ),
        by = "gene"
      )

    missing_df <- if (
      length(attribution$missing_genes) > 0L
    ) {
      data.frame(
        panel = rep(
          panel_name,
          length(attribution$missing_genes)
        ),
        condition = rep(
          condition_name,
          length(attribution$missing_genes)
        ),
        gene = attribution$missing_genes,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        panel = character(0),
        condition = character(0),
        gene = character(0),
        stringsAsFactors = FALSE
      )
    }

    all_per_lineage[[condition_name]] <- per_lineage
    all_summary[[condition_name]] <- summary_df
    all_missing[[condition_name]] <- missing_df

    readr::write_csv(
      per_lineage,
      file.path(
        panel_output_dir,
        paste0(
          "reference_expression_by_lineage_",
          safe_panel_name(condition_name),
          ".csv"
        )
      )
    )

    readr::write_csv(
      summary_df,
      file.path(
        panel_output_dir,
        paste0(
          "reference_attribution_summary_",
          safe_panel_name(condition_name),
          ".csv"
        )
      )
    )

    message(
      "Detected ",
      nrow(summary_df),
      " of ",
      nrow(panel),
      " panel genes."
    )

    if (length(attribution$missing_genes) > 0L) {
      message(
        "Missing genes: ",
        paste(
          attribution$missing_genes,
          collapse = ", "
        )
      )
    }
  }

  combined_per_lineage <- dplyr::bind_rows(
    all_per_lineage
  )

  combined_summary <- dplyr::bind_rows(
    all_summary
  )

  combined_missing <- dplyr::bind_rows(
    all_missing
  )

  # Compare dominant lineage assignments across conditions.
  condition_comparison <- combined_summary |>
    dplyr::select(
      gene,
      condition,
      dominant_lineage,
      dominant_celltype,
      tau_specificity
    ) |>
    tidyr::pivot_wider(
      names_from = condition,
      values_from = c(
        dominant_lineage,
        dominant_celltype,
        tau_specificity
      ),
      names_sep = "_"
    )

  if (
    all(
      c(
        "dominant_lineage_Control",
        "dominant_lineage_AD+CAA"
      ) %in% names(condition_comparison)
    )
  ) {
    condition_comparison <- condition_comparison |>
      dplyr::mutate(
        lineage_shift =
          !is.na(dominant_lineage_Control) &
          !is.na(`dominant_lineage_AD+CAA`) &
          dominant_lineage_Control !=
          `dominant_lineage_AD+CAA`
      )
  }

  readr::write_csv(
    combined_per_lineage,
    file.path(
      panel_output_dir,
      "reference_expression_by_lineage_all_conditions.csv"
    )
  )

  readr::write_csv(
    combined_summary,
    file.path(
      panel_output_dir,
      "reference_attribution_summary_all_conditions.csv"
    )
  )

  readr::write_csv(
    condition_comparison,
    file.path(
      panel_output_dir,
      "reference_attribution_condition_comparison.csv"
    )
  )

  readr::write_csv(
    combined_missing,
    file.path(
      panel_output_dir,
      "missing_reference_genes.csv"
    )
  )

  message(
    "\nReference attribution completed for panel: ",
    panel_name
  )

  # ==========================================================
  # Spatial perturbation
  # ==========================================================

  message(
    "\nRunning spatial perturbation models for panel: ",
    panel_name
  )

  spatial_perturbation <-
    test_gene_panel_spatial_perturbation(
      expr_mat = expression_normalized,
      roi_meta = roi_meta,
      genes = panel$gene,
      roi_id_col = "ROI_ID",
      disease_col = "disease_status",
      pathology_col = "pathology",
      region_col = "region",
      scan_col = "Scan_ID",
      adjust_scope = "contrast_region"
    )

  spatial_perturbation <- spatial_perturbation |>
    dplyr::mutate(
      panel = panel_name,
      .before = 1
    ) |>
    dplyr::left_join(
      panel |>
        dplyr::select(
          gene,
          category,
          headline,
          source
        ),
      by = "gene"
    )

  readr::write_csv(
    spatial_perturbation,
    file.path(
      panel_output_dir,
      "spatial_perturbation_by_region.csv"
    )
  )

  message(
    "Spatial perturbation models completed: ",
    sum(spatial_perturbation$model_status == "ok"),
    " contrast rows; ",
    sum(spatial_perturbation$model_status != "ok"),
    " failed or unevaluable rows."
  )

  # ==========================================================
  # Composition-attributed perturbation
  # ==========================================================

  message(
    "\nRunning lineage-associated perturbation models for panel: ",
    panel_name
  )

  lineage_attribution <- test_gene_panel_lineage_attribution(
    expr_mat = expression_normalized,
    proportions_long = proportions_long,
    roi_meta = roi_meta,
    genes = panel$gene,
    roi_id_col = "ROI_ID",
    celltype_col = "celltype",
    proportion_col = "rel_abundance",
    disease_col = "disease_status",
    pathology_col = "pathology",
    region_col = "region",
    scan_col = "Scan_ID",
    min_n = 20L,
    min_scans = 3L
  )

  lineage_attribution <- lineage_attribution |>
    dplyr::mutate(
      panel = panel_name,
      .before = 1
    ) |>
    dplyr::left_join(
      panel |>
        dplyr::select(
          gene,
          category,
          headline,
          source
        ),
      by = "gene"
    )

  readr::write_csv(
    lineage_attribution,
    file.path(
      panel_output_dir,
      "lineage_attributed_perturbation.csv"
    )
  )

  message(
    "Lineage attribution completed: ",
    sum(lineage_attribution$model_status == "ok"),
    " evaluable models; ",
    sum(lineage_attribution$model_status != "ok"),
    " failed or unevaluable models."
  )

  # ==========================================================
  # Fine-cell-type discovery correlations
  # ==========================================================

  message(
    "\nRunning fine-cell-type discovery correlations for panel: ",
    panel_name
  )

  fine_celltype_correlations <-
    correlate_gene_panel_with_celltypes(
      expr_mat = expression_normalized,
      proportions_long = proportions_long,
      roi_meta = roi_meta,
      genes = panel$gene,
      roi_id_col = "ROI_ID",
      celltype_col = "celltype",
      proportion_col = "rel_abundance",
      disease_col = "disease_status",
      pathology_col = "pathology",
      scan_col = "Scan_ID",
      disease_subset = "AD-CAA",
      method = "spearman",
      min_n = 8L,
      min_scans = 3L
    )

  fine_celltype_correlations <- fine_celltype_correlations |>
    dplyr::mutate(
      panel = panel_name,
      .before = 1
    ) |>
    dplyr::left_join(
      panel |>
        dplyr::select(
          gene,
          category,
          headline,
          source
        ),
      by = "gene"
    )

  readr::write_csv(
    fine_celltype_correlations,
    file.path(
      panel_output_dir,
      "fine_celltype_discovery_correlations.csv"
    )
  )

  message(
    "Fine-cell-type correlations completed: ",
    sum(fine_celltype_correlations$model_status == "ok"),
    " evaluable tests; ",
    sum(fine_celltype_correlations$model_status != "ok"),
    " failed or unevaluable tests."
  )

  # ==========================================================
  # Integrated per-gene summary
  # ==========================================================

  message(
    "\nBuilding integrated gene-panel summary for: ",
    panel_name
  )

  gene_panel_summary <- build_gene_panel_summary(
    panel = panel,
    reference_summary = combined_summary,
    spatial_results = spatial_perturbation,
    lineage_results = lineage_attribution,
    fine_results = fine_celltype_correlations,
    reference_condition = "Control",
    spatial_fdr_cutoff = 0.05,
    lineage_fdr_cutoff = 0.05,
    fine_fdr_cutoff = 0.05
  )

  readr::write_csv(
    gene_panel_summary,
    file.path(
      panel_output_dir,
      "integrated_gene_panel_summary.csv"
    )
  )

  message(
    "Integrated summary written: ",
    nrow(gene_panel_summary),
    " genes; ",
    sum(
      gene_panel_summary$spatial_evidence ==
        "FDR-significant",
      na.rm = TRUE
    ),
    " with spatial FDR evidence; ",
    sum(
      gene_panel_summary$fine_celltype_evidence ==
        "Global FDR-significant",
      na.rm = TRUE
    ),
    " with fine-cell-type FDR evidence."
  )
}
message(
  "\nGene-panel reference validation completed successfully."
)
