#' Signature export utilities
#'
#' Helper functions for exporting snRNA-seq reference signatures for
#' downstream spatial deconvolution.
#'
#' @keywords internal
NULL

#' Add combined cell class and condition label
#'
#' @param seu Seurat object.
#' @param cellclass_col Metadata column containing cell class labels.
#' @param condition_col Metadata column containing condition labels.
#' @param output_col Name of combined metadata column.
#' @param sep Separator.
#'
#' @return Seurat object.
#' @export
add_cellclass_condition_label <- function(seu,
                                          cellclass_col = "final_cellclass",
                                          condition_col = "FDX",
                                          output_col = "cellclass_condition",
                                          sep = "_") {
  required_cols <- c(cellclass_col, condition_col)
  missing_cols <- setdiff(required_cols, colnames(seu@meta.data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  seu[[output_col]] <- paste0(
    seu[[cellclass_col]][, 1],
    sep,
    seu[[condition_col]][, 1]
  )

  return(seu)
}

#' Export aggregate expression signatures
#'
#' @param seu Seurat object.
#' @param group_col Metadata column to aggregate by.
#' @param assay Assay to aggregate.
#' @param layer Expression layer to export.
#' @param file Output file.
#'
#' @return Aggregated expression matrix.
#' @export
export_aggregate_signatures <- function(seu,
                                        group_col = "cellclass_condition",
                                        assay = "RNA",
                                        layer = "data",
                                        file = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.", call. = FALSE)
  }

  if (!group_col %in% colnames(seu@meta.data)) {
    stop("group_col not found in metadata: ", group_col, call. = FALSE)
  }

  agg <- Seurat::AggregateExpression(
    seu,
    group.by = group_col,
    assays = assay,
    return.seurat = TRUE
  )

  mat <- Seurat::GetAssayData(
    agg,
    assay = assay,
    layer = layer
  )

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    write.table(
      mat,
      file = file,
      sep = "\t",
      quote = FALSE
    )
  }

  return(mat)
}

#' Export collapsed aggregate signatures
#'
#' @param seu Seurat object.
#' @param collapsed_col Collapsed cell class column.
#' @param condition_col Condition column.
#' @param output_col Combined collapsed-condition column.
#' @param assay Assay to aggregate.
#' @param layer Layer to export.
#' @param file Output file.
#'
#' @return Aggregated matrix.
#' @export
export_collapsed_aggregate_signatures <- function(seu,
                                                  collapsed_col = "final_collapsed_cellclass",
                                                  condition_col = "FDX",
                                                  output_col = "collapsed_condition",
                                                  assay = "RNA",
                                                  layer = "data",
                                                  file = NULL) {
  seu <- add_cellclass_condition_label(
    seu = seu,
    cellclass_col = collapsed_col,
    condition_col = condition_col,
    output_col = output_col
  )

  export_aggregate_signatures(
    seu = seu,
    group_col = output_col,
    assay = assay,
    layer = layer,
    file = file
  )
}

#' Export Seurat counts, metadata, and genes for AnnData conversion
#'
#' @param seu Seurat object.
#' @param output_dir Output directory.
#' @param assay Assay to export.
#' @param counts_layer Counts layer.
#' @param prefix File prefix.
#'
#' @return Invisibly returns output file paths.
#' @export
export_seurat_to_anndata_inputs <- function(seu,
                                            output_dir,
                                            assay = "RNA",
                                            counts_layer = "counts",
                                            prefix = "sn_reference") {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  counts <- Seurat::GetAssayData(
    seu,
    assay = assay,
    layer = counts_layer
  )

  counts_file <- file.path(output_dir, paste0(prefix, "_counts.mtx"))
  metadata_file <- file.path(output_dir, paste0(prefix, "_metadata.csv"))
  genes_file <- file.path(output_dir, paste0(prefix, "_genes.csv"))

  Matrix::writeMM(counts, counts_file)

  write.csv(
    seu@meta.data,
    file = metadata_file,
    quote = TRUE
  )

  write.csv(
    data.frame(gene_id = rownames(counts)),
    file = genes_file,
    row.names = FALSE,
    quote = TRUE
  )

  invisible(list(
    counts = counts_file,
    metadata = metadata_file,
    genes = genes_file
  ))
}

#' Save reference metadata summary
#'
#' @param seu Seurat object.
#' @param group_cols Metadata columns to summarize.
#' @param file Output CSV file.
#'
#' @return Summary data frame.
#' @export
save_reference_metadata_summary <- function(seu,
                                            group_cols = c("final_cellclass", "FDX"),
                                            file = NULL) {
  missing_cols <- setdiff(group_cols, colnames(seu@meta.data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df <- as.data.frame(
    table(seu@meta.data[, group_cols])
  )

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    write.csv(df, file = file, row.names = FALSE)
  }

  return(df)
}
