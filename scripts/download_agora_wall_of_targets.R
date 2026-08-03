#!/usr/bin/env Rscript

# ============================================================
# Download the AGORA Wall of Targets
#
# Outputs:
#   1. Full AGORA Wall of Targets
#   2. Tier 1 high-nomination panel
#   3. Download metadata
#
# Tier 1 is defined as targets with at least 3 independent
# AGORA nominations. This prioritizes repeatedly nominated AD
# targets while keeping the downstream spatial, lineage, and
# fine-cell-type analyses comparable in scale to the Hajjar panel.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
})

source("R/agora_utils.R")

# ============================================================
# Configuration
# ============================================================

output_dir <- file.path(
  "resources",
  "gene_panels"
)

full_output_file <- file.path(
  output_dir,
  "agora_wall_of_targets.csv"
)

tier1_output_file <- file.path(
  output_dir,
  "agora_high_nomination_targets.csv"
)

metadata_output_file <- file.path(
  output_dir,
  "agora_wall_of_targets_download_metadata.csv"
)

minimum_tier1_nominations <- 3L
headline_nomination_threshold <- 4L

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Download the complete AGORA target resource
# ============================================================

message("Downloading AGORA nominated targets...")

targets <- download_agora_targets(
  page_size = 100L,
  max_retries = 3L,
  retry_wait = 2,
  verbose = TRUE
)

downloaded_at <- attr(
  targets,
  "downloaded_at"
)

api_total_elements <- attr(
  targets,
  "api_total_elements"
)

if (is.null(downloaded_at)) {
  downloaded_at <- format(
    Sys.time(),
    tz = "UTC",
    usetz = TRUE
  )
}

if (is.null(api_total_elements)) {
  api_total_elements <- nrow(
    targets
  )
}

# ============================================================
# Validate the complete download
# ============================================================

required_columns <- c(
  "gene",
  "panel",
  "category",
  "headline",
  "source",
  "total_nominations"
)

missing_columns <- setdiff(
  required_columns,
  names(targets)
)

if (length(missing_columns) > 0L) {
  stop(
    "Downloaded AGORA table is missing required column(s): ",
    paste(
      missing_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

invalid_gene_rows <- is.na(targets$gene) |
  !nzchar(
    trimws(
      as.character(targets$gene)
    )
  )

if (any(invalid_gene_rows)) {
  stop(
    "Downloaded AGORA table contains ",
    sum(invalid_gene_rows),
    " missing or empty gene symbol(s).",
    call. = FALSE
  )
}

if (anyDuplicated(targets$gene)) {
  duplicated_genes <- unique(
    targets$gene[
      duplicated(targets$gene)
    ]
  )

  stop(
    "Downloaded AGORA table contains duplicated gene symbol(s): ",
    paste(
      utils::head(
        duplicated_genes,
        10
      ),
      collapse = ", "
    ),
    call. = FALSE
  )
}

# ============================================================
# Write the complete Wall of Targets
# ============================================================

readr::write_csv(
  targets,
  full_output_file
)

message(
  "Wrote ",
  full_output_file,
  " with ",
  nrow(targets),
  " unique genes."
)

# ============================================================
# Create the Tier 1 high-nomination panel
# ============================================================

tier1_targets <- create_agora_tier1_panel(
  agora_targets = targets,
  minimum_nominations = minimum_tier1_nominations,
  panel_name = "AGORA_High_Nomination",
  headline_threshold = headline_nomination_threshold
)

readr::write_csv(
  tier1_targets,
  tier1_output_file
)

message(
  "Wrote ",
  tier1_output_file,
  " with ",
  nrow(tier1_targets),
  " Tier 1 targets."
)

# ============================================================
# Write reproducibility metadata
# ============================================================

nomination_distribution <- table(
  targets$total_nominations,
  useNA = "ifany"
)

tier1_nomination_distribution <- table(
  tier1_targets$total_nominations,
  useNA = "ifany"
)

metadata <- data.frame(
  resource = "AGORA Wall of Targets",
  source = "AGORA AD Knowledge Portal",
  api_endpoint = paste0(
    "https://agora.adknowledgeportal.org/",
    "api/v1/comparison-tools/targets"
  ),
  downloaded_at_utc = downloaded_at,
  api_total_elements = as.integer(
    api_total_elements
  ),
  retained_unique_genes = nrow(
    targets
  ),
  tier1_definition = paste0(
    "total_nominations >= ",
    minimum_tier1_nominations
  ),
  tier1_headline_definition = paste0(
    "total_nominations >= ",
    headline_nomination_threshold
  ),
  tier1_retained_genes = nrow(
    tier1_targets
  ),
  full_nomination_distribution = paste(
    names(nomination_distribution),
    as.integer(nomination_distribution),
    sep = "=",
    collapse = "; "
  ),
  tier1_nomination_distribution = paste(
    names(tier1_nomination_distribution),
    as.integer(tier1_nomination_distribution),
    sep = "=",
    collapse = "; "
  ),
  stringsAsFactors = FALSE
)

readr::write_csv(
  metadata,
  metadata_output_file
)

message(
  "Wrote download metadata to ",
  metadata_output_file,
  "."
)

# ============================================================
# Console summary
# ============================================================

message(
  "\nAGORA download completed successfully."
)

message(
  "Full Wall of Targets: ",
  nrow(targets),
  " genes."
)

message(
  "Tier 1 panel: ",
  nrow(tier1_targets),
  " genes with at least ",
  minimum_tier1_nominations,
  " nominations."
)

message(
  "Tier 1 headline targets: ",
  sum(
    tier1_targets$headline,
    na.rm = TRUE
  ),
  " genes with at least ",
  headline_nomination_threshold,
  " nominations."
)

message(
  "Tier 1 nomination distribution: ",
  paste(
    names(tier1_nomination_distribution),
    as.integer(tier1_nomination_distribution),
    sep = "=",
    collapse = ", "
  ),
  "."
)
