# ============================================================
# Contrast model robustness check (leave-one-Scan-out)
#
# With only a handful of distinct Scan_IDs backing the random intercept
# in each beta mixed-effects contrast model, a single scan can have an
# outsized influence on an effect estimate. This refits each of the four
# contrast functions once per excluded Scan_ID and flags region x
# celltype x contrast results whose significance call flips, or whose
# estimate swings sharply, when any single scan is dropped.
# ============================================================

library(glmmTMB)
library(emmeans)
library(dplyr)

source("R/contrast_utils.R")
source("R/robustness_utils.R")

# ============================================================
# Config -- adjust paths for your environment
# ============================================================

input_file <- "results/cell_proportions/spatial_celltype_proportions_for_R.csv"
output_dir <- "results/spatial_stats/robustness"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load and clean data
# ============================================================

df <- read.csv(input_file, stringsAsFactors = FALSE)

df <- prepare_spatial_proportion_data(
  df = df,
  abundance_col = "rel_abundance",
  disease_col = "disease_status",
  pathology_col = "pathology",
  region_col = "region",
  scan_col = "Scan_ID",
  celltype_col = "celltype"
)

message(sprintf("%d distinct Scan_ID(s) in the full dataset.", dplyr::n_distinct(df$Scan_ID)))

# ============================================================
# Run leave-one-Scan-out for each contrast function
# ============================================================

contrast_fns <- list(
  amyloid_effect = run_amyloid_effect,
  disease_effect = run_disease_effect,
  overall_effect = run_weighted_overall_effect,
  max_pathology_effect = run_max_pathology_effect
)

expected_contrast <- c(
  amyloid_effect = "Amyloid_effect",
  disease_effect = "Disease_effect",
  overall_effect = "AD_overall_vs_Control",
  max_pathology_effect = "AD_Amyloid_vs_Control"
)


for (name in names(contrast_fns)) {

  message(sprintf("Running leave-one-Scan-out for %s...", name))

  loo <- run_leave_one_scan_out(
    df = df,
    model_fn = contrast_fns[[name]],
    abundance_col = "rel_abundance"
  )

  summary_df <- summarize_robustness(loo$full, loo$leave_one_out)

  if (
    "contrast" %in% colnames(summary_df) &&
      name %in% names(expected_contrast)
  ) {
  
    n_before <- nrow(summary_df)

    summary_df <- summary_df[
      summary_df$contrast == expected_contrast[[name]],
      ,
      drop = FALSE
    ]

    if (nrow(summary_df) != n_before) {
      message(
        sprintf(
          "  Filtered robustness results from %d to %d rows (%s).",
          n_before,
          nrow(summary_df),
          expected_contrast[[name]]
        )
      )
    }
  }

  write.csv(
    summary_df,
    file = file.path(output_dir, paste0(name, "_robustness.csv")),
    row.names = FALSE
  )

  if (nrow(summary_df) == 0) {
    message("  No leave-one-out refits succeeded; nothing to summarize.")
    next
  }

  flagged <- summary_df[summary_df$any_sig_flip, ]

  if (nrow(flagged) > 0) {
    message(sprintf(
      "  %d of %d region x celltype result(s) flip significance when a single scan is dropped.",
      nrow(flagged), nrow(summary_df)
    ))
  } else {
    message("  No significance flips under leave-one-Scan-out.")
  }
}

message("Contrast robustness check completed.")
