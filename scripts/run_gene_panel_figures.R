#!/usr/bin/env Rscript

# ============================================================
# Generate gene-panel validation figures
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("R/viz_theme.R")
source("R/viz_gene_panel.R")

result_root <- "results/gene_panel_validation"
figure_root <- "results/figures/gene_panel_validation"

if (!dir.exists(result_root)) {
  stop(
    "Gene-panel result directory does not exist: ",
    result_root,
    call. = FALSE
  )
}

panel_dirs <- list.dirs(
  result_root,
  recursive = FALSE,
  full.names = TRUE
)

if (length(panel_dirs) == 0L) {
  stop(
    "No panel-specific result directories found under ",
    result_root,
    ".",
    call. = FALSE
  )
}

for (panel_dir in panel_dirs) {

  panel_name <- basename(
    panel_dir
  )

  message(
    "\nGenerating figures for panel: ",
    panel_name
  )

  panel_figure_dir <- file.path(
    figure_root,
    panel_name
  )

  # ----------------------------------------------------------
  # Reference-expression heatmap
  # ----------------------------------------------------------

  reference_csv <- file.path(
    panel_dir,
    "reference_expression_by_lineage_all_conditions.csv"
  )

  if (file.exists(reference_csv)) {
    plot_gene_panel_reference_heatmap(
      reference_csv = reference_csv,
      output_dir = panel_figure_dir
    )
  } else {
    message(
      "Skipping reference heatmap; missing input: ",
      reference_csv
    )
  }

  # ----------------------------------------------------------
  # Spatial-perturbation heatmap
  # ----------------------------------------------------------

  spatial_csv <- file.path(
    panel_dir,
    "spatial_perturbation_by_region.csv"
  )

  if (file.exists(spatial_csv)) {
    plot_gene_panel_spatial_perturbation(
      spatial_csv = spatial_csv,
      output_dir = panel_figure_dir,
      padj_cutoff = 0.05
    )
  } else {
    message(
      "Skipping spatial perturbation heatmap; missing input: ",
      spatial_csv
    )
  }

  # ----------------------------------------------------------
  # Lineage-associated perturbation heatmap
  # ----------------------------------------------------------

  lineage_csv <- file.path(
    panel_dir,
    "lineage_attributed_perturbation.csv"
  )

  if (file.exists(lineage_csv)) {
    plot_gene_panel_lineage_attribution(
      lineage_csv = lineage_csv,
      output_dir = panel_figure_dir,
      fdr_cutoff = 0.05,
      nominal_cutoff = 0.05
    )
  } else {
    message(
      "Skipping lineage-attribution heatmap; missing input: ",
      lineage_csv
    )
  }

  # ----------------------------------------------------------
  # Fine-cell-type discovery correlation heatmap
  # ----------------------------------------------------------

  fine_celltype_csv <- file.path(
    panel_dir,
    "fine_celltype_discovery_correlations.csv"
  )

  if (file.exists(fine_celltype_csv)) {
    plot_gene_panel_fine_celltype_correlations(
      correlation_csv = fine_celltype_csv,
      output_dir = panel_figure_dir,
      fdr_cutoff = 0.05
    )
  } else {
    message(
      "Skipping fine-cell-type correlation heatmap; missing input: ",
      fine_celltype_csv
    )
  }


}

message(
  "\nAll available gene-panel figures written under ",
  figure_root
)
