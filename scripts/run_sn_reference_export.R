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
