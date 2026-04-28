# ============================================================
# snRNA-seq Reference Preparation for Multi-omic Spatial Integration
# ============================================================
#
# Repository:
# multiomic-spatial-integration
#
# Upstream Framework:
# This script is adapted from the modular preprocessing framework developed in:
# neuro-snRNAseq-tools https://github.com/Echimalj/neuro-snRNAseq-tools
#
# The goal is to reuse standardized snRNA-seq preprocessing utilities from
# neuro-snRNAseq-tools and apply them here to generate an annotated reference
# suitable for GeoMx WTA spatial deconvolution.
#
# ------------------------------------------------------------
# Description:
# This script prepares an annotated single-nucleus RNA-seq reference from
# human AD/CAA and control samples for downstream integration with GeoMx
# Whole Transcriptome Atlas (WTA) spatial transcriptomics.
#
# The original project-specific workflow included:
# - SoupX ambient RNA correction
# - sample metadata annotation
# - QC filtering
# - DoubletFinder-based doublet removal
# - SCTransform normalization
# - Harmony integration
# - cluster marker detection
# - manual cell-type annotation
# - targeted subclustering of astrocytes, microglia, and vascular cells
# - Pearson-correlation-guided merging of highly similar subclusters
# - aggregation of expression by cell type and disease condition
#
# This cleaned version keeps the biological workflow intact while replacing
# hard-coded, sample-specific code with reusable functions and sample-sheet
# driven processing from neuro-snRNAseq-tools.
#
# ------------------------------------------------------------
# Inputs:
# - sample sheet CSV with:
#     sample_id
#     sample_path
#     condition / FDX
#     batch
#     species
#     use_soupx
#     orig.ident
#
# - Cell Ranger output folders for each sample
#
# ------------------------------------------------------------
# Outputs:
# - annotated Seurat reference object
# - QC summaries
# - doublet summaries
# - cluster marker tables
# - subcluster marker tables
# - Pearson correlation matrix between subclusters
# - final cell-type labels
# - aggregated expression matrix for spatial integration
#
# ------------------------------------------------------------
# Example downstream use:
# The final aggregated expression matrix is used as input for
# Cell2Location / SpaceJam-style Bayesian spatial deconvolution of GeoMx
# WTA data.
#
# ------------------------------------------------------------
# Author:
# Enrique Chimal
# ============================================================

# ============================================================
# Link to neuro-snRNAseq-tools repository
# ============================================================

# Option 1: set manually (recommended for HPC / reproducibility)
neuro_tools_path <- "/path/to/neuro-snRNAseq-tools"

# Option 2: relative path (if repos are cloned side-by-side)
# neuro_tools_path <- "../neuro-snRNAseq-tools"

# Option 3: environment variable (advanced / portable)
# neuro_tools_path <- Sys.getenv("NEURO_SN_RNASEQ_TOOLS")

if (!dir.exists(neuro_tools_path)) {
  stop(
    "neuro-snRNAseq-tools path not found. Please set 'neuro_tools_path' correctly.",
    call. = FALSE
  )
}

neuro_R <- file.path(neuro_tools_path, "R")

# ============================================================
# Load reusable snRNA-seq utilities
# ============================================================

source(file.path(neuro_R, "check_dependencies.R"))
source(file.path(neuro_R, "load_samples.R"))
source(file.path(neuro_R, "metadata_utils.R"))
source(file.path(neuro_R, "qc_utils.R"))
source(file.path(neuro_R, "doublet_utils.R"))
source(file.path(neuro_R, "preprocess_utils.R"))
source(file.path(neuro_R, "integration_utils.R"))
source(file.path(neuro_R, "marker_utils.R"))
source(file.path(neuro_R, "annotation_utils.R"))
source(file.path(neuro_R, "subcluster_utils.R"))
source(file.path(neuro_R, "checkpoint_utils.R"))

check_required_packages(c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tibble",
  "ggplot2",
  "patchwork",
  "sctransform",
  "harmony",
  "clustree"
))

check_required_packages(c(
  "SoupX",
  "DoubletFinder"
))

# ============================================================
# Project configuration
# ============================================================

project_dir <- "/path/to/multiomic-spatial-integration"

sample_sheet_file <- file.path(project_dir, "data", "human_sample_sheet.csv")

output_dir <- file.path(project_dir, "results", "sn_reference")
checkpoint_dir <- file.path(project_dir, "checkpoints")
figure_dir <- file.path(project_dir, "figures", "sn_reference")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load sample sheet
# ============================================================

sample_sheet <- read.csv(sample_sheet_file)

required_cols <- c(
  "sample_id",
  "sample_path",
  "FDX",
  "batch",
  "species",
  "use_soupx",
  "orig.ident"
)

missing_cols <- setdiff(required_cols, colnames(sample_sheet))

if (length(missing_cols) > 0) {
  stop(
    "Sample sheet is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

sample_sheet

# ============================================================
# Load samples
# ============================================================

seurat_list <- load_samples_from_sheet(
  sample_sheet = sample_sheet,
  min_cells = 3,
  min_features = 200,
  use_soupx = TRUE
)

seurat_list <- annotate_samples_from_sheet(
  seurat_list = seurat_list,
  sample_sheet = sample_sheet
)

AD_CAA <- merge_seurat_samples(
  seurat_list = seurat_list,
  project_name = "AD_CAA_2025"
)

summarize_cells_by_group(AD_CAA, group_col = "orig.ident")
summarize_cells_by_group(AD_CAA, group_col = "FDX")

## CheckPoint
AD_CAA <- checkpoint(
  object = AD_CAA,
  file = file.path(checkpoint_dir, "AD_CAA_01_loaded_merged.rds"),
  reload = TRUE
)

# ============================================================
# QC metrics and filtering
# ============================================================

AD_CAA <- add_qc_metrics(
  seu = AD_CAA,
  species = "human",
  add_ribo = FALSE,
  add_hb = FALSE
)

qc_plot <- plot_basic_qc(AD_CAA)

ggplot2::ggsave(
  filename = file.path(figure_dir, "AD_CAA_basic_qc.pdf"),
  plot = qc_plot,
  width = 10,
  height = 5
)

AD_CAA_before_qc <- AD_CAA

AD_CAA <- filter_cells_basic(
  seu = AD_CAA,
  min_features = 200,
  max_mt = 1
)

qc_summary <- summarize_filtering(
  before = AD_CAA_before_qc,
  after = AD_CAA,
  group_col = "orig.ident"
)

write.csv(
  qc_summary,
  file = file.path(output_dir, "AD_CAA_QC_filtering_summary.csv"),
  row.names = FALSE
)

## CheckPoint
AD_CAA <- checkpoint(
  object = AD_CAA,
  file = file.path(checkpoint_dir, "AD_CAA_02_after_qc.rds"),
  reload = TRUE
)
# ============================================================
# Doublet detection and singlet filtering
# ============================================================

AD_CAA <- annotate_doublets_by_sample(
  seu = AD_CAA,
  split_by = "orig.ident",
  assay = "RNA",
  sct = FALSE,
  resolution = 0.1
)

doublet_summary <- summarize_doublets(
  AD_CAA,
  group_col = "orig.ident"
)

write.csv(
  doublet_summary,
  file = file.path(output_dir, "AD_CAA_doublet_summary.csv"),
  row.names = FALSE
)

AD_CAA <- keep_singlets(AD_CAA)

## CheckPoint
AD_CAA <- checkpoint(
  object = AD_CAA,
  file = file.path(checkpoint_dir, "AD_CAA_03_singlets_only.rds"),
  reload = TRUE
)

