#' Enrichment preparation utilities
#' @keywords internal
NULL

#' Extract top genes per celltype
#'
#' @param markers DEG table
#' @param padj_cutoff threshold
#'
#' @export
extract_top_genes <- function(markers,
                              padj_cutoff = 0.05) {

  markers$gene[markers$p_val_adj < padj_cutoff]
}

#' Format gene list for enrichment tools
#'
#' @param genes vector
#'
#' @export
format_gene_list <- function(genes) {
  unique(na.omit(genes))
}
