#!/usr/bin/env Rscript

# ============================================================
# Build gene-panel summaries, rankings, and Markdown reports
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source("R/gene_panel_reporting.R")

result_root <- "results/gene_panel_validation"

figure_root <- file.path(
  "results",
  "figures",
  "gene_panel_validation"
)

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
    "No panel result directories found under ",
    result_root,
    ".",
    call. = FALSE
  )
}

panel_summaries <- list()

for (panel_dir in panel_dirs) {
  panel_name <- basename(
    panel_dir
  )

  message(
    "\nBuilding reporting outputs for panel: ",
    panel_name
  )

  integrated_file <- file.path(
    panel_dir,
    "integrated_gene_panel_summary.csv"
  )

  missing_file <- file.path(
    panel_dir,
    "missing_reference_genes.csv"
  )

  if (!file.exists(integrated_file)) {
    message(
      "Skipping panel; missing input: ",
      integrated_file
    )
    next
  }

  integrated_summary <- readr::read_csv(
    integrated_file,
    show_col_types = FALSE
  )

  reference_missing <- if (
    file.exists(missing_file)
  ) {
    readr::read_csv(
      missing_file,
      show_col_types = FALSE
    )
  } else {
    data.frame(
      panel = character(0),
      condition = character(0),
      gene = character(0)
    )
  }

  ranked_genes <- rank_gene_panel_priorities(
    integrated_summary
  )

  ranked_file <- file.path(
    panel_dir,
    "ranked_gene_priorities.csv"
  )

  readr::write_csv(
    ranked_genes,
    ranked_file
  )

  message(
    "Wrote ",
    ranked_file
  )

  panel_summary <- summarize_gene_panel(
    integrated_summary = integrated_summary,
    reference_missing = reference_missing
  )

  panel_summaries[[panel_name]] <- panel_summary

  report_file <- file.path(
    panel_dir,
    "report.md"
  )

  write_gene_panel_markdown_report(
    integrated_summary = integrated_summary,
    ranked_genes = ranked_genes,
    panel_summary = panel_summary,
    panel_result_dir = panel_dir,
    panel_figure_dir = file.path(
      figure_root,
      panel_name
    ),
    output_file = report_file,
    top_n = 15L
  )
}

combined_panel_summary <- dplyr::bind_rows(
  panel_summaries
)

summary_file <- file.path(
  result_root,
  "panel_summary.csv"
)

readr::write_csv(
  combined_panel_summary,
  summary_file
)

message(
  "\nWrote ",
  summary_file,
  " with ",
  nrow(combined_panel_summary),
  " panel(s)."
)

message(
  "\nGene-panel reporting completed successfully."
)
