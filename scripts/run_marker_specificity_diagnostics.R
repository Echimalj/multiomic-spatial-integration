source("R/pathway_proportion_utils.R")
source("R/marker_concordance_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

expression_csv <- "results/geomx_exports/CAA-AD_expression_wide.csv"
proportions_csv <- "results/cell_proportions/roi_celltype_abundance_long.csv"
markers_file <- "results/AD_CAA_cluster_markers.txt"

output_dir <- "results/marker_specificity_diagnostics"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

partial_spearman <- function(x, y, z) {
  keep <- is.finite(x) & is.finite(y) & is.finite(z)

  x <- x[keep]
  y <- y[keep]
  z <- z[keep]

  if (length(x) < 5) {
    return(NA_real_)
  }

  if (
    stats::sd(x) == 0 ||
    stats::sd(y) == 0 ||
    stats::sd(z) == 0
  ) {
    return(NA_real_)
  }

  # Partial Spearman correlation:
  # rank variables, regress x and y on ranked covariate,
  # then correlate residuals.
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  rz <- rank(z, ties.method = "average")

  x_resid <- stats::residuals(stats::lm(rx ~ rz))
  y_resid <- stats::residuals(stats::lm(ry ~ rz))

  stats::cor(
    x_resid,
    y_resid,
    method = "pearson",
    use = "complete.obs"
  )
}


get_best_matches <- function(cor_mat) {
  inferred_types <- rownames(cor_mat)

  out <- lapply(inferred_types, function(ct) {
    vals <- cor_mat[ct, ]

    if (all(is.na(vals))) {
      return(data.frame(
        inferred_celltype = ct,
        best_marker_celltype = NA_character_,
        best_rho = NA_real_,
        self_rho = NA_real_,
        self_is_best = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    best_idx <- which.max(replace(vals, is.na(vals), -Inf))
    best_marker <- colnames(cor_mat)[best_idx]

    self_rho <- if (ct %in% colnames(cor_mat)) {
      unname(cor_mat[ct, ct])
    } else {
      NA_real_
    }

    data.frame(
      inferred_celltype = ct,
      best_marker_celltype = best_marker,
      best_rho = unname(vals[best_idx]),
      self_rho = self_rho,
      self_is_best = identical(ct, best_marker),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out) |>
    arrange(desc(best_rho))
}


write_matrix <- function(x, filename) {
  write.csv(
    x,
    file.path(output_dir, filename),
    row.names = TRUE,
    quote = FALSE
  )
}


plot_heatmap <- function(mat, filename, main) {
  pdf(
    file.path(output_dir, filename),
    width = 12,
    height = 12
  )

  heatmap(
    mat,
    Rowv = NA,
    Colv = NA,
    scale = "none",
    margins = c(10, 10),
    main = main
  )

  dev.off()
}


# -------------------------------------------------------------------------
# Load and normalize GeoMx expression
# -------------------------------------------------------------------------

expr_raw <- load_roi_expression(expression_csv)
expr_norm <- normalize_expression_cpm(expr_raw)

cat(
  "Expression matrix:",
  nrow(expr_norm), "ROIs x",
  ncol(expr_norm), "genes\n"
)


# -------------------------------------------------------------------------
# Build expression-aware marker lists
# -------------------------------------------------------------------------

markers_df <- read.delim(
  markers_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_marker_cols <- c(
  "p_val_adj",
  "avg_log2FC",
  "cluster",
  "gene"
)

missing_marker_cols <- setdiff(
  required_marker_cols,
  colnames(markers_df)
)

if (length(missing_marker_cols) > 0) {
  stop(
    "Marker file is missing required columns: ",
    paste(missing_marker_cols, collapse = ", ")
  )
}

marker_table <- markers_df |>
  filter(
    is.finite(p_val_adj),
    p_val_adj <= 0.05,
    gene %in% colnames(expr_norm)
  ) |>
  distinct(cluster, gene, .keep_all = TRUE) |>
  group_by(cluster) |>
  arrange(desc(avg_log2FC), .by_group = TRUE) |>
  slice_head(n = 20) |>
  ungroup()

marker_lists_df <- marker_table |>
  group_by(cluster) |>
  summarise(
    genes = list(gene),
    n_markers = n(),
    .groups = "drop"
  )

marker_lists <- stats::setNames(
  marker_lists_df$genes,
  marker_lists_df$cluster
)

write.csv(
  marker_table,
  file.path(output_dir, "expression_aware_marker_table.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  marker_lists_df |>
    select(cluster, n_markers),
  file.path(output_dir, "marker_counts_by_celltype.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat(
  "Expression-aware marker programs:",
  length(marker_lists), "\n"
)


# -------------------------------------------------------------------------
# Compute marker-score matrix
# -------------------------------------------------------------------------

score_list <- lapply(names(marker_lists), function(ct) {
  genes <- marker_lists[[ct]]

  score <- try(
    compute_marker_score(
      expr_norm,
      genes,
      min_markers = 3
    ),
    silent = TRUE
  )

  if (inherits(score, "try-error")) {
    warning(
      "Marker score failed for ", ct, ": ",
      as.character(score)
    )
    return(NULL)
  }

  score
})

names(score_list) <- names(marker_lists)

score_list <- score_list[
  !vapply(score_list, is.null, logical(1))
]

if (length(score_list) == 0) {
  stop("No marker-score programs could be calculated.")
}

score_mat <- do.call(cbind, score_list)
colnames(score_mat) <- names(score_list)


# -------------------------------------------------------------------------
# Build abundance matrix
# -------------------------------------------------------------------------

prop <- read.csv(
  proportions_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_prop_cols <- c(
  "ROI_ID",
  "celltype",
  "rel_abundance"
)

missing_prop_cols <- setdiff(
  required_prop_cols,
  colnames(prop)
)

if (length(missing_prop_cols) > 0) {
  stop(
    "Proportion file is missing required columns: ",
    paste(missing_prop_cols, collapse = ", ")
  )
}

abund_wide <- prop |>
  select(ROI_ID, celltype, rel_abundance) |>
  pivot_wider(
    id_cols = ROI_ID,
    names_from = celltype,
    values_from = rel_abundance,
    values_fill = 0,
    values_fn = sum
  )

abund_wide <- as.data.frame(abund_wide)

rownames(abund_wide) <- as.character(abund_wide$ROI_ID)
abund_wide$ROI_ID <- NULL

abund_mat <- as.matrix(abund_wide)
storage.mode(abund_mat) <- "numeric"


# -------------------------------------------------------------------------
# Align ROIs
# -------------------------------------------------------------------------

shared_rois <- intersect(
  rownames(score_mat),
  rownames(abund_mat)
)

if (length(shared_rois) == 0) {
  stop(
    "No shared ROI identifiers between marker-score ",
    "and abundance matrices."
  )
}

score_mat <- score_mat[shared_rois, , drop = FALSE]
abund_mat <- abund_mat[shared_rois, , drop = FALSE]

cat(
  "Score matrix:",
  nrow(score_mat), "x", ncol(score_mat), "\n"
)

cat(
  "Abundance matrix:",
  nrow(abund_mat), "x", ncol(abund_mat), "\n"
)

cat("Shared ROIs:", length(shared_rois), "\n")


# -------------------------------------------------------------------------
# Correlation structure within each data type
# -------------------------------------------------------------------------

abundance_cor <- stats::cor(
  abund_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

marker_score_cor <- stats::cor(
  score_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

write_matrix(
  abundance_cor,
  "abundance_abundance_spearman.csv"
)

write_matrix(
  marker_score_cor,
  "marker_score_marker_score_spearman.csv"
)


# -------------------------------------------------------------------------
# Raw abundance-marker correlation
# -------------------------------------------------------------------------

raw_cor <- stats::cor(
  abund_mat,
  score_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

write_matrix(
  raw_cor,
  "abundance_marker_raw_spearman.csv"
)

raw_best <- get_best_matches(raw_cor)

write.csv(
  raw_best,
  file.path(output_dir, "raw_best_marker_matches.csv"),
  row.names = FALSE,
  quote = FALSE
)


# -------------------------------------------------------------------------
# Partial correlations controlling overall excitatory abundance
# -------------------------------------------------------------------------

excitatory_types <- grep(
  "^ExNeuron",
  colnames(abund_mat),
  value = TRUE
)

if (length(excitatory_types) < 2) {
  stop(
    "Fewer than two ExNeuron abundance columns were found."
  )
}

total_excitatory <- rowSums(
  abund_mat[, excitatory_types, drop = FALSE],
  na.rm = TRUE
)

partial_cor <- matrix(
  NA_real_,
  nrow = ncol(abund_mat),
  ncol = ncol(score_mat),
  dimnames = list(
    colnames(abund_mat),
    colnames(score_mat)
  )
)

for (inferred_ct in rownames(partial_cor)) {

  # For an excitatory subtype, avoid controlling for a total that
  # contains the target itself. Use all other excitatory subtypes.
  if (inferred_ct %in% excitatory_types) {
    other_excitatory <- setdiff(
      excitatory_types,
      inferred_ct
    )

    covariate <- rowSums(
      abund_mat[, other_excitatory, drop = FALSE],
      na.rm = TRUE
    )
  } else {
    covariate <- total_excitatory
  }

  for (marker_ct in colnames(partial_cor)) {
    partial_cor[inferred_ct, marker_ct] <- partial_spearman(
      x = abund_mat[, inferred_ct],
      y = score_mat[, marker_ct],
      z = covariate
    )
  }
}

write_matrix(
  partial_cor,
  "abundance_marker_partial_spearman.csv"
)

partial_best <- get_best_matches(partial_cor)

write.csv(
  partial_best,
  file.path(output_dir, "partial_best_marker_matches.csv"),
  row.names = FALSE,
  quote = FALSE
)


# -------------------------------------------------------------------------
# Excitatory-only matrices and summaries
# -------------------------------------------------------------------------

shared_excitatory <- intersect(
  excitatory_types,
  colnames(score_mat)
)

raw_excitatory <- raw_cor[
  shared_excitatory,
  shared_excitatory,
  drop = FALSE
]

partial_excitatory <- partial_cor[
  shared_excitatory,
  shared_excitatory,
  drop = FALSE
]

abundance_excitatory_cor <- abundance_cor[
  shared_excitatory,
  shared_excitatory,
  drop = FALSE
]

marker_excitatory_cor <- marker_score_cor[
  shared_excitatory,
  shared_excitatory,
  drop = FALSE
]

write_matrix(
  raw_excitatory,
  "excitatory_raw_abundance_marker.csv"
)

write_matrix(
  partial_excitatory,
  "excitatory_partial_abundance_marker.csv"
)

write_matrix(
  abundance_excitatory_cor,
  "excitatory_abundance_correlation.csv"
)

write_matrix(
  marker_excitatory_cor,
  "excitatory_marker_score_correlation.csv"
)


# -------------------------------------------------------------------------
# Compare raw and partial self-correlations
# -------------------------------------------------------------------------

comparison <- full_join(
  raw_best |>
    select(
      inferred_celltype,
      raw_best_marker = best_marker_celltype,
      raw_best_rho = best_rho,
      raw_self_rho = self_rho,
      raw_self_is_best = self_is_best
    ),
  partial_best |>
    select(
      inferred_celltype,
      partial_best_marker = best_marker_celltype,
      partial_best_rho = best_rho,
      partial_self_rho = self_rho,
      partial_self_is_best = self_is_best
    ),
  by = "inferred_celltype"
) |>
  mutate(
    self_rho_change = partial_self_rho - raw_self_rho
  ) |>
  arrange(desc(self_rho_change))

write.csv(
  comparison,
  file.path(output_dir, "raw_vs_partial_comparison.csv"),
  row.names = FALSE,
  quote = FALSE
)

excitatory_comparison <- comparison |>
  filter(grepl("^ExNeuron", inferred_celltype))

write.csv(
  excitatory_comparison,
  file.path(output_dir, "excitatory_raw_vs_partial_comparison.csv"),
  row.names = FALSE,
  quote = FALSE
)


# -------------------------------------------------------------------------
# Generate heatmaps
# -------------------------------------------------------------------------

plot_heatmap(
  abundance_excitatory_cor,
  "excitatory_abundance_correlation_heatmap.pdf",
  "Excitatory abundance correlations"
)

plot_heatmap(
  marker_excitatory_cor,
  "excitatory_marker_score_correlation_heatmap.pdf",
  "Excitatory marker-score correlations"
)

plot_heatmap(
  raw_excitatory,
  "excitatory_raw_abundance_marker_heatmap.pdf",
  "Raw abundance-marker correlations"
)

plot_heatmap(
  partial_excitatory,
  "excitatory_partial_abundance_marker_heatmap.pdf",
  "Partial abundance-marker correlations"
)


# -------------------------------------------------------------------------
# Write text summary
# -------------------------------------------------------------------------

summary_file <- file.path(
  output_dir,
  "marker_specificity_summary.txt"
)

sink(summary_file)

cat("MARKER SPECIFICITY DIAGNOSTIC SUMMARY\n")
cat("=====================================\n\n")

cat("ROIs:", length(shared_rois), "\n")
cat("Abundance factors:", ncol(abund_mat), "\n")
cat("Marker programs:", ncol(score_mat), "\n")
cat("Excitatory subtypes:", length(shared_excitatory), "\n\n")

cat("RAW CORRELATIONS\n")
cat("----------------\n")
cat(
  "Self is best:",
  sum(raw_best$self_is_best, na.rm = TRUE),
  "of", nrow(raw_best), "\n"
)
cat(
  "Positive self rho:",
  sum(raw_best$self_rho > 0, na.rm = TRUE),
  "of", nrow(raw_best), "\n\n"
)

cat("PARTIAL CORRELATIONS\n")
cat("--------------------\n")
cat(
  "Self is best:",
  sum(partial_best$self_is_best, na.rm = TRUE),
  "of", nrow(partial_best), "\n"
)
cat(
  "Positive self rho:",
  sum(partial_best$self_rho > 0, na.rm = TRUE),
  "of", nrow(partial_best), "\n\n"
)

cat("EXCITATORY RAW VS PARTIAL\n")
cat("-------------------------\n")
print(
  excitatory_comparison,
  row.names = FALSE
)

cat("\nTARGET SUBTYPES\n")
cat("---------------\n")

print(
  comparison |>
    filter(
      inferred_celltype %in% c(
        "ExNeuron6",
        "ExNeuron10",
        "ExNeuron12"
      )
    ),
  row.names = FALSE
)

sink()

cat("\nDiagnostics completed successfully.\n")
cat("Outputs written to:", output_dir, "\n")
cat("Summary:", summary_file, "\n")
