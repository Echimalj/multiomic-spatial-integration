#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Biobase)
  library(BisqueRNA)
})

source("deconvolution_comparison/R/common_utils.R")
source("deconvolution_comparison/R/adapter_bisque.R")

project_dir <- getwd()

output_dir <- Sys.getenv(
  "BISQUE_OUTPUT_DIR",
  unset = file.path(
    project_dir,
    "results",
    "deconvolution_comparison"
  )
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Load single-nucleus reference
# ------------------------------------------------------------

message("Loading reference...")

reference_dir <- file.path(
  project_dir,
  "data"
)

ref <- load_reference_mtx(
  counts_mtx = file.path(
    reference_dir,
    "AD_CAA_counts.mtx"
  ),
  metadata_csv = file.path(
    reference_dir,
    "AD_CAA_meta.csv"
  ),
  genes_csv = file.path(
    reference_dir,
    "AD_CAA_genes.csv"
  ),
  celltype_col = "New_Idents",
  sample_col = "orig.ident"
)

if (is.null(ref$counts) ||
    is.null(ref$cell_type) ||
    is.null(ref$sample)) {
  stop(
    paste(
      "Reference loader must return counts,",
      "cell_type, and sample."
    )
  )
}

if (length(ref$cell_type) != ncol(ref$counts)) {
  stop(
    "Length of ref$cell_type does not match reference cells."
  )
}

if (length(ref$sample) != ncol(ref$counts)) {
  stop(
    "Length of ref$sample does not match reference cells."
  )
}

biological_samples <- unique(
  as.character(ref$sample)
)

message(
  "Unique biological samples: ",
  paste(
    biological_samples,
    collapse = ", "
  )
)

message(
  "Reference: ",
  nrow(ref$counts),
  " genes x ",
  ncol(ref$counts),
  " cells."
)

message(
  "Cell types: ",
  length(unique(ref$cell_type))
)

message(
  "Biological samples: ",
  length(biological_samples)
)

message("Sample distribution:")

print(
  table(ref$sample)
)

message("Cell-type distribution:")

print(
  table(ref$cell_type)
)

# ------------------------------------------------------------
# Load GeoMx mixture
# ------------------------------------------------------------

message("Loading GeoMx mixture...")

mixture_file <- file.path(
  project_dir,
  "results",
  "geomx_exports",
  "CAA-AD_expression_wide.csv"
)

if (!file.exists(mixture_file)) {
  stop(
    "Mixture file does not exist: ",
    mixture_file
  )
}

mixture <- load_mixture(
  mixture_file
)

if (!is.matrix(mixture)) {
  mixture <- as.matrix(mixture)
}

storage.mode(mixture) <- "numeric"

message(
  "Mixture: ",
  nrow(mixture),
  " genes x ",
  ncol(mixture),
  " ROIs."
)

# ------------------------------------------------------------
# Gene overlap diagnostics
# ------------------------------------------------------------

shared_genes <- intersect(
  rownames(mixture),
  rownames(ref$counts)
)

message(
  "Shared genes: ",
  length(shared_genes)
)

if (length(shared_genes) == 0) {
  stop(
    "No shared genes were found between the mixture and reference."
  )
}

if (length(shared_genes) < 1000) {
  warning(
    "Only ",
    length(shared_genes),
    " shared genes were found."
  )
}

# ------------------------------------------------------------
# Memory diagnostic
#
# The current adapter converts the aligned reference to a dense
# matrix before constructing the ExpressionSet. This calculation
# reports the approximate size of that dense matrix.
# ------------------------------------------------------------

estimated_dense_gb <- (
  length(shared_genes) *
    ncol(ref$counts) *
    8
) / 1024^3

message(
  "Approximate dense reference size after gene alignment: ",
  round(estimated_dense_gb, 2),
  " GB, excluding temporary copies and object overhead."
)

if (estimated_dense_gb > 20) {
  warning(
    paste0(
      "Bisque may require substantially more than ",
      round(estimated_dense_gb, 1),
      " GB because adapter_bisque.R currently converts ",
      "the reference matrix to a dense matrix."
    )
  )
}

# ------------------------------------------------------------
# Run BisqueRNA
# ------------------------------------------------------------

message("Running BisqueRNA...")

start_time <- Sys.time()

bisque_prop <- run_bisque(
  ref = ref,
  mixture = mixture
)

elapsed <- difftime(
  Sys.time(),
  start_time,
  units = "mins"
)

if (is.null(bisque_prop)) {
  stop(
    paste(
      "BisqueRNA returned NULL.",
      "Confirm that BisqueRNA and Biobase are installed."
    )
  )
}

bisque_prop <- as.matrix(
  bisque_prop
)

storage.mode(bisque_prop) <- "numeric"

if (!is.numeric(bisque_prop)) {
  stop(
    "BisqueRNA output is not numeric."
  )
}

if (is.null(rownames(bisque_prop))) {
  stop(
    "BisqueRNA output has no ROI row names."
  )
}

if (is.null(colnames(bisque_prop))) {
  stop(
    "BisqueRNA output has no cell-type column names."
  )
}

message(
  "BisqueRNA elapsed time: ",
  round(
    as.numeric(elapsed),
    2
  ),
  " minutes."
)

# ------------------------------------------------------------
# Validate output dimensions and values
# ------------------------------------------------------------

raw_sums <- rowSums(
  bisque_prop,
  na.rm = TRUE
)

message(
  "BisqueRNA result: ",
  nrow(bisque_prop),
  " ROIs x ",
  ncol(bisque_prop),
  " cell types."
)

finite_values <- bisque_prop[
  is.finite(bisque_prop)
]

if (length(finite_values) > 0) {
  message(
    "Raw value range: ",
    paste(
      range(finite_values),
      collapse = " to "
    )
  )
} else {
  warning(
    "BisqueRNA output contains no finite values."
  )
}

message(
  "Zero or failed ROI rows: ",
  sum(
    !is.finite(raw_sums) |
      raw_sums <= 0
  )
)

message(
  "Negative coefficients: ",
  sum(
    bisque_prop < 0,
    na.rm = TRUE
  )
)

expected_rois <- ncol(mixture)
expected_celltypes <- length(
  unique(ref$cell_type)
)

if (nrow(bisque_prop) != expected_rois) {
  stop(
    "Expected ",
    expected_rois,
    " ROI rows, but BisqueRNA returned ",
    nrow(bisque_prop),
    "."
  )
}

if (ncol(bisque_prop) != expected_celltypes) {
  warning(
    "Expected ",
    expected_celltypes,
    " cell types, but BisqueRNA returned ",
    ncol(bisque_prop),
    "."
  )
}

if (anyDuplicated(rownames(bisque_prop))) {
  stop(
    "BisqueRNA output contains duplicated ROI names."
  )
}

if (anyDuplicated(colnames(bisque_prop))) {
  stop(
    "BisqueRNA output contains duplicated cell-type names."
  )
}

# Reorder ROIs to match the input mixture when all names are present.
if (setequal(
  rownames(bisque_prop),
  colnames(mixture)
)) {
  bisque_prop <- bisque_prop[
    colnames(mixture),
    ,
    drop = FALSE
  ]
} else {
  missing_rois <- setdiff(
    colnames(mixture),
    rownames(bisque_prop)
  )

  extra_rois <- setdiff(
    rownames(bisque_prop),
    colnames(mixture)
  )

  if (length(missing_rois) > 0) {
    warning(
      "ROIs missing from BisqueRNA output: ",
      paste(
        head(missing_rois, 10),
        collapse = ", "
      )
    )
  }

  if (length(extra_rois) > 0) {
    warning(
      "Unexpected ROIs in BisqueRNA output: ",
      paste(
        head(extra_rois, 10),
        collapse = ", "
      )
    )
  }
}

# ------------------------------------------------------------
# Standardize and save
# ------------------------------------------------------------

long <- standardize_proportions(
  prop_mat = bisque_prop,
  method = "BisqueRNA",
  normalize = TRUE
)

save_method_output(
  long_df = long,
  method = "BisqueRNA",
  output_dir = output_dir
)

message(
  "Standardized output: ",
  nrow(long),
  " rows; ",
  length(unique(long$ROI_ID)),
  " ROIs; ",
  length(unique(long$celltype)),
  " cell types."
)

output_file <- file.path(
  output_dir,
  "BisqueRNA_proportions.csv"
)

message(
  "Output written to: ",
  output_file
)

message(
  "BisqueRNA run completed successfully."
)
