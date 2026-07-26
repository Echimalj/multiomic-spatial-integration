# ============================================================
# Pathway recurrence summary
#
# Summarizes how broadly significant pathways recur across cell-type and
# stratum combinations. Frequently recurring pathways may represent shared
# biology, pathway redundancy, broad transcriptional programs,
# compositional structure, or technical effects.
#
# Low-recurrence, non-housekeeping pathway associations are highlighted as
# candidate cell-type-restricted signals for follow-up interpretation.
# ============================================================

library(dplyr)

source("R/viz_theme.R")
source("R/pathway_proportion_utils.R")
source("R/viz_pathway.R")

pathway_dir <- file.path(
  "results",
  "pathway_proportion_link"
)

fig_dir <- file.path(
  "results",
  "figures",
  "pathway_recurrence"
)

PADJ_CUTOFF <- 0.05

# Maximum fraction of tested celltype–stratum combinations in which a
# pathway may be significant and still be considered "low recurrence".
MAX_FREQUENCY_PROP <- 0.10

TOP_N_PER_STRATUM <- 15L

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

combined <- load_pathway_enrichment_results(pathway_dir)

message(sprintf(
  "Loaded %d pathway-association rows.",
  nrow(combined)
))

recurrence <- summarize_pathway_recurrence(
  combined,
  padj_cutoff = PADJ_CUTOFF
)

n_combinations <- recurrence$n_tested_combinations

max_frequency <- max(
  1L,
  floor(MAX_FREQUENCY_PROP * n_combinations)
)

message(sprintf(
  "Detected %d unique celltype-stratum combinations.",
  n_combinations
))

n_sig <- nrow(recurrence$significant)

n_broad <- sum(
  recurrence$significant$is_housekeeping,
  na.rm = TRUE
)

if (n_sig > 0L) {
  message(sprintf(
    "%d / %d (%.1f%%) significant pathway associations match curated broad pathway families.",
    n_broad,
    n_sig,
    100 * n_broad / n_sig
  ))
} else {
  message("No significant pathway associations detected.")
}

top_low_recurrence <-
  top_low_recurrence_associations_by_stratum(
    recurrence,
    max_frequency = max_frequency,
    top_n = TOP_N_PER_STRATUM
  )

if (nrow(top_low_recurrence) > 0L) {

  n_strata <- dplyr::n_distinct(top_low_recurrence$stratum)

  message(sprintf(
    paste0(
      "%d low-recurrence associations retained ",
      "across %d %s ",
      "(max recurrence = %d/%d combinations; %.2f%%)."
    ),
    nrow(top_low_recurrence),
    n_strata,
    ifelse(n_strata == 1L, "stratum", "strata"),
    max_frequency,
    n_combinations,
    100 * max_frequency / n_combinations
  ))

} else {

  message(
    "No low-recurrence pathway associations passed the recurrence filter."
  )

}

frequency_file <- file.path(
  pathway_dir,
  "pathway_recurrence_summary.csv"
)

top_file <- file.path(
  pathway_dir,
  "top_low_recurrence_associations_by_stratum.csv"
)

utils::write.csv(
  recurrence$frequency,
  frequency_file,
  row.names = FALSE
)

utils::write.csv(
  top_low_recurrence,
  top_file,
  row.names = FALSE
)

if (nrow(top_low_recurrence) > 0L) {

  message("\nTop low-recurrence association per stratum:")

  for (stratum_name in unique(top_low_recurrence$stratum)) {

    row <- top_low_recurrence[
      top_low_recurrence$stratum == stratum_name,
      ,
      drop = FALSE
    ][1, ]

    message(sprintf(
      paste0(
        "  %-28s %-20s %-50s ",
        "NES=%6.2f padj=%8.2e ",
        "(recurrence=%d; lineages=%s)"
      ),
      stratum_name,
      row$celltype,
      substr(row$pathway, 1, 48),
      row$NES,
      row$padj,
      row$n_significant,
      row$cell_lineages
    ))
  }

   plot_pathway_recurrence_dotplot(
     top_low_recurrence,
     fig_dir
   )

} else {

  message(
    "\nNo pathway associations passed the recurrence filter; ",
    "no summary plot was generated."
  )

}

message("\nWrote: ", frequency_file)
message("Wrote: ", top_file)
message("Pathway recurrence summary completed.")
