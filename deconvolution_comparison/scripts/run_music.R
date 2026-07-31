suppressPackageStartupMessages({
  library(Matrix)
  library(MuSiC)
  library(SingleCellExperiment)
})

source("deconvolution_comparison/R/common_utils.R")
source("deconvolution_comparison/R/adapter_music.R")

project_dir <- getwd()

output_dir <- Sys.getenv(
  "MUSIC_OUTPUT_DIR",
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

message("Loading reference...")

reference_dir <- "data"

ref <- load_reference_mtx(
  counts_mtx = file.path(reference_dir, "AD_CAA_counts.mtx"),
  metadata_csv = file.path(reference_dir, "AD_CAA_meta.csv"),
  genes_csv = file.path(reference_dir, "AD_CAA_genes.csv"),
  celltype_col = "New_Idents",
  sample_col = "orig.ident"
)

message(
  "Unique biological samples: ",
  paste(levels(ref$sample), collapse = ", ")
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
  length(unique(ref$sample))
)

message(
  "Sample distribution:"
)

print(table(ref$sample))

message("Loading GeoMx mixture...")

mixture <- load_mixture(
  "results/geomx_exports/CAA-AD_expression_wide.csv"
)

message(
  "Mixture: ",
  nrow(mixture),
  " genes x ",
  ncol(mixture),
  " ROIs."
)

message("Running MuSiC...")

start_time <- Sys.time()

message(
  "Shared genes: ",
  length(intersect(
    rownames(mixture),
    rownames(ref$counts)
  ))
)

music_prop <- run_music(
  ref = ref,
  mixture = mixture
)

elapsed <- difftime(
  Sys.time(),
  start_time,
  units = "mins"
)

stopifnot(
  is.matrix(music_prop),
  is.numeric(music_prop)
)

message(
  "MuSiC elapsed time: ",
  round(as.numeric(elapsed), 2),
  " minutes."
)

if (is.null(music_prop))
  stop("MuSiC returned NULL.")

raw_sums <- rowSums(
  music_prop,
  na.rm = TRUE
)

message(
  "MuSiC result: ",
  nrow(music_prop),
  " ROIs x ",
  ncol(music_prop),
  " cell types."
)

message(
  "Raw value range: ",
  paste(
    range(
      music_prop,
      na.rm = TRUE
    ),
    collapse = " to "
  )
)

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
    music_prop < 0,
    na.rm = TRUE
  )
)

stopifnot(
  nrow(music_prop) == 190,
  ncol(music_prop) == 46,
  !anyDuplicated(rownames(music_prop)),
  !anyDuplicated(colnames(music_prop))
)

long <- standardize_proportions(
  prop_mat = music_prop,
  method = "MuSiC",
  normalize = TRUE
)

save_method_output(
  long_df = long,
  method = "MuSiC",
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

message("MuSiC run completed successfully.")
