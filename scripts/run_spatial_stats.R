# ============================================================
# Spatial statistical modeling (beta mixed-effects)
# Multi-contrast framework:
# - Amyloid effect
# - Disease effect
# - Overall effect (weighted)
# - Max pathology comparison
# ============================================================

library(glmmTMB)
library(emmeans)
library(dplyr)
library(tidyr)
library(readr)

source("R/contrast_utils.R")
source("R/plotting_utils.R")

# Inputs / outputs

input_file <- "results/cell_proportions/spatial_celltype_proportions_for_R.csv"
output_dir <- "results/spatial_stats"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load and clean data
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

# ============================================================
# Run contrasts
# ============================================================

amyloid_res <- run_amyloid_effect(
  df = df,
  abundance_col = "rel_abundance"
)

disease_res <- run_disease_effect(
  df = df,
  abundance_col = "rel_abundance"
)

overall_res <- run_weighted_overall_effect(
  df = df,
  abundance_col = "rel_abundance"
)

maxpath_res <- run_max_pathology_effect(
  df = df,
  abundance_col = "rel_abundance"
)

# ============================================================
# Save individual contrast outputs
# ============================================================

write.csv(
  amyloid_res$contrasts,
  file = file.path(output_dir, "amyloid_effect_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  disease_res$contrasts,
  file = file.path(output_dir, "disease_effect_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  overall_res$contrasts,
  file = file.path(output_dir, "overall_effect_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  maxpath_res$contrasts,
  file = file.path(output_dir, "max_pathology_effect_contrasts.csv"),
  row.names = FALSE
)

write.csv(
  overall_res$weights,
  file = file.path(output_dir, "overall_effect_region_weights.csv"),
  row.names = FALSE
)

# ============================================================
# Combined summary table
# ============================================================

amyloid_df <- format_contrast_summary(
  amyloid_res$contrasts,
  contrast_type = "Amyloid"
)

disease_df <- format_contrast_summary(
  disease_res$contrasts,
  contrast_type = "Disease"
)

overall_df <- format_contrast_summary(
  overall_res$contrasts %>% filter(contrast == "AD_overall_vs_Control"),
  contrast_type = "Overall"
)

maxpath_df <- format_contrast_summary(
  maxpath_res$contrasts,
  contrast_type = "MaxPathology"
)

combined_contrasts <- bind_rows(
  amyloid_df,
  disease_df,
  overall_df,
  maxpath_df
)

write.csv(
  combined_contrasts,
  file = file.path(output_dir, "combined_spatial_contrast_summary.csv"),
  row.names = FALSE
)

message("Spatial statistical modeling completed successfully.")
