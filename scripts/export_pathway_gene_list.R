# ============================================================
# Export pathway gene list
# ============================================================

source("R/enrichment_utils.R")

repo_root <- normalizePath(".", mustWork = TRUE)

input_file <- file.path(
  repo_root,
  "results",
  "AD_CAA_cluster_markers.txt"
)

output_file <- file.path(
  repo_root,
  "results",
  "pathway_gene_list.txt"
)

if (!file.exists(input_file)) {
  stop(
    "Marker file not found: ",
    input_file,
    call. = FALSE
  )
}

markers <- utils::read.delim(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(markers) == 0L) {
  stop(
    "Marker table contains no rows: ",
    input_file,
    call. = FALSE
  )
}

genes <- extract_top_genes(markers)
genes <- format_gene_list(genes)
genes <- unique(genes[!is.na(genes) & nzchar(genes)])

if (length(genes) == 0L) {
  stop(
    "No valid genes were extracted from the marker table.",
    call. = FALSE
  )
}

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

utils::write.table(
  genes,
  file = output_file,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

message(sprintf(
  "Exported %d unique genes to %s",
  length(genes),
  output_file
))
