#' Pathway-linkage utilities for cell-type proportions
#'
#' Cell-type proportions are compositional (they sum to ~1 within an ROI),
#' so pathway/GSEA tools can't be run on them directly. These utilities
#' instead correlate a cell type's per-ROI proportion with per-ROI gene
#' expression, rank genes by that association, and run rank-based gene set
#' enrichment (fgsea) on the ranked list -- i.e. "what pathways covary with
#' this cell type's abundance shifts."
#'
#' @keywords internal
NULL


#' Load a wide ROI x gene expression matrix
#'
#' @param expression_csv Path to a CSV with ROI ids as row names and genes
#'   as columns (see `python.geomx_anndata_utils.export_expression_csv()`).
#'
#' @return Numeric matrix, rows = ROI, columns = gene.
#' @export
load_roi_expression <- function(expression_csv) {
  expr_df <- read.csv(expression_csv, row.names = 1, check.names = FALSE)
  as.matrix(expr_df)
}


#' Per-ROI log2-CPM normalization
#'
#' Normalizes each ROI (row) by its own total count before log-transforming,
#' so correlation with cell-type proportion isn't confounded by sequencing
#' depth.
#'
#' @param expr_mat Numeric matrix, rows = ROI, columns = gene (raw counts).
#' @param pseudocount Added before taking log2.
#'
#' @return Numeric matrix, same shape as `expr_mat`.
#' @export
normalize_expression_cpm <- function(expr_mat, pseudocount = 1) {
  lib_sizes <- rowSums(expr_mat)

  if (any(lib_sizes == 0)) {
    stop(
      "normalize_expression_cpm: ", sum(lib_sizes == 0),
      " ROI(s) have zero total counts; remove them before normalizing.",
      call. = FALSE
    )
  }

  cpm <- sweep(expr_mat, 1, lib_sizes, "/") * 1e6

  log2(cpm + pseudocount)
}


#' Rank genes by association with a cell type's proportion
#'
#' Joins `proportions_df` (filtered to one cell type) to `expr_mat` by ROI
#' id, then computes a single vectorized correlation of every gene column
#' against the proportion vector (not a per-gene loop), converting
#' coefficients to p-values via the standard t-approximation
#' (`t = r * sqrt(n-2) / sqrt(1-r^2)`), BH-adjusted across genes.
#'
#' @param expr_mat Numeric matrix, rows = ROI, columns = gene (typically
#'   already normalized via `normalize_expression_cpm()`).
#' @param proportions_df Long-format cell-type proportion table (e.g.
#'   `results/cell_proportions/spatial_celltype_proportions_for_R.csv`).
#' @param celltype Cell type to test (must appear in `celltype_col`).
#' @param roi_id_col Column in `proportions_df` identifying the ROI,
#'   matching `rownames(expr_mat)`.
#' @param proportion_col Column with the relative abundance value.
#' @param celltype_col Column with cell-type labels.
#' @param method Correlation method passed to `stats::cor()`.
#' @param min_shared_rois Minimum overlapping ROIs required between
#'   `expr_mat` and `proportions_df` before proceeding.
#'
#' @return Data frame `gene, statistic, p.value, p_adj`, sorted by
#'   `statistic` descending.
#' @export
rank_genes_by_celltype_association <- function(expr_mat,
                                               proportions_df,
                                               celltype,
                                               roi_id_col = "ROI_ID",
                                               proportion_col = "rel_abundance",
                                               celltype_col = "celltype",
                                               method = "spearman",
                                               min_shared_rois = 5) {
  if (!roi_id_col %in% colnames(proportions_df)) {
    stop(
      "rank_genes_by_celltype_association: roi_id_col '", roi_id_col,
      "' not found in proportions_df. Available columns: ",
      paste(colnames(proportions_df), collapse = ", "),
      call. = FALSE
    )
  }

  dat_ct <- proportions_df[proportions_df[[celltype_col]] == celltype, ]

  if (nrow(dat_ct) == 0) {
    stop(
      "rank_genes_by_celltype_association: no rows found for celltype '",
      celltype, "' in column '", celltype_col, "'.",
      call. = FALSE
    )
  }

  shared_rois <- intersect(rownames(expr_mat), dat_ct[[roi_id_col]])

  if (length(shared_rois) < min_shared_rois) {
    stop(
      "rank_genes_by_celltype_association: only ", length(shared_rois),
      " shared ROI(s) between expr_mat (", nrow(expr_mat),
      " rows) and proportions_df (", nrow(dat_ct), " rows for '", celltype,
      "') via roi_id_col = '", roi_id_col, "' -- check that expr_mat's ",
      "row names and proportions_df[[roi_id_col]] use the same ROI ",
      "identifier convention.",
      call. = FALSE
    )
  }

  rownames(dat_ct) <- dat_ct[[roi_id_col]]
  dat_ct <- dat_ct[shared_rois, ]
  expr_sub <- expr_mat[shared_rois, , drop = FALSE]

  proportion_vec <- dat_ct[[proportion_col]]

  r <- as.numeric(stats::cor(expr_sub, proportion_vec, method = method))
  n <- length(proportion_vec)

  t_stat <- r * sqrt(n - 2) / sqrt(1 - r^2)
  p_value <- 2 * stats::pt(-abs(t_stat), df = n - 2)

  result <- data.frame(
    gene = colnames(expr_sub),
    statistic = r,
    p.value = p_value,
    stringsAsFactors = FALSE
  )

  result$p_adj <- stats::p.adjust(result$p.value, method = "BH")

  result[order(-result$statistic), ]
}


#' Run rank-based gene set enrichment on a ranked gene list
#'
#' @param ranked_genes Data frame from `rank_genes_by_celltype_association()`
#'   (`gene`, `statistic` columns used).
#' @param gmt_file Path to a GMT gene-set file (e.g. MSigDB Hallmark/KEGG/
#'   Reactome/GO). If `NULL`, enrichment is skipped (the ranked gene list is
#'   still useful on its own, e.g. via manual Enrichr upload).
#' @param min_size Minimum pathway size passed to `fgsea::fgsea()`.
#' @param max_size Maximum pathway size passed to `fgsea::fgsea()`.
#'
#' @return Data frame of fgsea results sorted by adjusted p-value, or `NULL`
#'   if `gmt_file` is `NULL`.
#' @export
run_pathway_enrichment <- function(ranked_genes,
                                   gmt_file,
                                   min_size = 15,
                                   max_size = 500) {
  if (is.null(gmt_file)) {
    message("run_pathway_enrichment: gmt_file is NULL, skipping enrichment.")
    return(NULL)
  }

  pathways <- fgsea::gmtPathways(gmt_file)

  ranks <- stats::setNames(ranked_genes$statistic, ranked_genes$gene)
  ranks <- sort(ranks, decreasing = TRUE)

  res <- fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    minSize = min_size,
    maxSize = max_size
  )

  res <- as.data.frame(res)
  res[order(res$padj), ]
}

# ============================================================
# Pathway recurrence summaries
# ============================================================


#' Broad pathway families commonly recurring across enrichment analyses
#'
#' Manually curated pathway-name terms used to flag broad pathway families
#' that frequently recur across many enrichment analyses. These categories
#' may represent shared biology, pathway redundancy, broad transcriptional
#' programs, compositional structure, or technical effects.
#'
#' The resulting flag is an interpretive aid only. It should not be treated
#' as evidence that a pathway is artifactual or biologically uninformative.
#'
#' @keywords internal
.default_broad_pathway_terms <- c(
  "Respiratory",
  "Electron Transport",
  "Ribosom",
  "rRNA",
  "Translation",
  "Spliceosome",
  "Processing Of Capped",
  "mRNA Splicing",
  "Nonsense.Mediated",
  "SRP.dependent",
  "Citric Acid",
  "Olfactory",
  "Selenocysteine",
  "Selenoamino",
  "Influenza",
  "Complex I Biogenesis",
  "Cristae Formation"
)


#' Match pathways to broad pathway families
#'
#' @param pathways Character vector of pathway names.
#' @param broad_terms Character vector of case-insensitive regular-expression
#'   terms defining broad pathway families.
#'
#' @return Character vector containing the first matching term for each
#'   pathway, or `NA_character_` when no term matches.
#'
#' @keywords internal
.match_broad_pathway_term <- function(
    pathways,
    broad_terms = .default_broad_pathway_terms
) {
  pathways <- as.character(pathways)

  if (length(broad_terms) == 0L) {
    return(rep(NA_character_, length(pathways)))
  }

  matched_term <- rep(NA_character_, length(pathways))

  for (term in broad_terms) {
    is_match <- is.na(matched_term) &
      !is.na(pathways) &
      grepl(
        term,
        pathways,
        ignore.case = TRUE,
        perl = TRUE
      )

    matched_term[is_match] <- term
  }

  matched_term
}


#' Load pathway-enrichment results from a linkage directory
#'
#' Recursively loads every `pathway_enrichment.csv` file under a directory
#' structured as:
#'
#' ```
#' <pathway_dir>/<stratum>/<celltype>/pathway_enrichment.csv
#' ```
#'
#' Each result is tagged with the stratum and cell type inferred from its
#' relative path.
#'
#' @param pathway_dir Root pathway-linkage results directory.
#'
#' @return Data frame containing `stratum`, `celltype`, `pathway`, `pval`,
#'   `padj`, `NES`, and `size`.
#'
#' @export
load_pathway_enrichment_results <- function(pathway_dir) {
  if (
    length(pathway_dir) != 1L ||
    is.na(pathway_dir) ||
    !nzchar(pathway_dir)
  ) {
    stop(
      "load_pathway_enrichment_results: `pathway_dir` must be a ",
      "single non-empty path.",
      call. = FALSE
    )
  }

  if (!dir.exists(pathway_dir)) {
    stop(
      "load_pathway_enrichment_results: directory not found: ",
      pathway_dir,
      call. = FALSE
    )
  }

  pathway_root <- normalizePath(
    pathway_dir,
    winslash = "/",
    mustWork = TRUE
  )

  files <- list.files(
    pathway_root,
    pattern = "^pathway_enrichment\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0L) {
    stop(
      "load_pathway_enrichment_results: no pathway_enrichment.csv ",
      "files found under: ",
      pathway_root,
      call. = FALSE
    )
  }

  required_cols <- c(
    "pathway",
    "pval",
    "padj",
    "NES",
    "size"
  )

  read_one <- function(file_path) {
    normalized_file <- normalizePath(
      file_path,
      winslash = "/",
      mustWork = TRUE
    )

    relative_path <- substring(
      normalized_file,
      nchar(pathway_root) + 2L
    )

    path_parts <- strsplit(
      relative_path,
      "/",
      fixed = TRUE
    )[[1L]]

    if (length(path_parts) < 3L) {
      stop(
        "load_pathway_enrichment_results: expected file layout ",
        "'<stratum>/<celltype>/pathway_enrichment.csv', but found: ",
        relative_path,
        call. = FALSE
      )
    }

    result <- utils::read.csv(
      normalized_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    missing_cols <- setdiff(
      required_cols,
      names(result)
    )

    if (length(missing_cols) > 0L) {
      stop(
        "load_pathway_enrichment_results: file is missing required ",
        "column(s): ",
        paste(missing_cols, collapse = ", "),
        "\nFile: ",
        normalized_file,
        call. = FALSE
      )
    }

    if (nrow(result) == 0L) {
      return(NULL)
    }

    result$stratum <- path_parts[[1L]]
    result$celltype <- path_parts[[2L]]

    result[
      ,
      c(
        "stratum",
        "celltype",
        "pathway",
        "pval",
        "padj",
        "NES",
        "size"
      ),
      drop = FALSE
    ]
  }

  result_list <- lapply(files, read_one)

  result_list <- Filter(
    Negate(is.null),
    result_list
  )

  if (length(result_list) == 0L) {
    stop(
      "load_pathway_enrichment_results: all pathway-enrichment ",
      "files were empty under: ",
      pathway_root,
      call. = FALSE
    )
  }

  combined <- do.call(
    rbind,
    result_list
  )

  rownames(combined) <- NULL

  combined
}


#' Infer broad cell lineages from cell-type labels
#'
#' Exact mappings supplied through `lineage_map` take precedence over the
#' default pattern-based assignments. Unmatched cell-type labels are retained
#' as their own lineage.
#'
#' @param celltypes Character vector of cell-type labels.
#' @param lineage_map Optional named character vector mapping exact cell-type
#'   labels to broader lineage labels.
#'
#' @return Character vector of lineage labels.
#'
#' @keywords internal
.infer_cell_lineage <- function(
    celltypes,
    lineage_map = NULL
) {
  celltypes <- as.character(celltypes)

  lineages <- rep(
    NA_character_,
    length(celltypes)
  )

  if (!is.null(lineage_map)) {
    if (
      !is.character(lineage_map) ||
      is.null(names(lineage_map)) ||
      any(is.na(names(lineage_map))) ||
      any(!nzchar(names(lineage_map))) ||
      anyDuplicated(names(lineage_map))
    ) {
      stop(
        ".infer_cell_lineage: `lineage_map` must be a named ",
        "character vector with unique, non-empty names.",
        call. = FALSE
      )
    }

    map_index <- match(
      celltypes,
      names(lineage_map)
    )

    has_exact_match <- !is.na(map_index)

    lineages[has_exact_match] <-
      unname(
        lineage_map[map_index[has_exact_match]]
      )
  }

  assign_pattern <- function(
      pattern,
      lineage
  ) {
    unresolved <- is.na(lineages) |
      !nzchar(lineages)

    matches <- unresolved &
      !is.na(celltypes) &
      nzchar(celltypes) &
      grepl(
        pattern,
        celltypes,
        ignore.case = TRUE,
        perl = TRUE
      )

    lineages[matches] <<- lineage
  }

  assign_pattern("^Astro", "Astrocyte")
  assign_pattern("^Microglia", "Microglia")
  assign_pattern("^(OPC|Oligodendrocyte precursor)", "OPC")
  assign_pattern("^Oligodendro", "Oligodendrocyte")
  assign_pattern("^(ExNeuron|Excitatory)", "Excitatory neuron")
  assign_pattern("^(InNeuron|Inhibitory)", "Inhibitory neuron")
  assign_pattern("^Endothelial", "Vascular")
  assign_pattern("^Pericyte", "Vascular")
  assign_pattern("^Fibroblast", "Vascular")
  assign_pattern("^VLMC", "Vascular")
  assign_pattern("^(SMC|Smooth[ _-]?muscle)", "Vascular")
  assign_pattern("^(Macrophage|Monocyte)", "Myeloid")
  assign_pattern("^(T[ _-]?cell|Tcell)", "T cell")
  assign_pattern("^(B[ _-]?cell|Bcell)", "B cell")

  unresolved <- is.na(lineages) |
    !nzchar(lineages)

  lineages[unresolved] <- celltypes[unresolved]

  lineages[
    is.na(lineages) |
      !nzchar(lineages)
  ] <- "Unknown"

  lineages
}

#' Summarize pathway recurrence across cell-type-stratum combinations
#'
#' Counts the number of distinct cell-type-stratum combinations in which each
#' pathway is significant. The function also annotates whether recurrence is
#' restricted to one biological lineage or distributed across multiple
#' lineages.
#'
#' @param combined_df Data frame containing pathway-enrichment results.
#' @param padj_cutoff Adjusted p-value threshold used to define significance.
#' @param broad_terms Named character vector or other object accepted by
#'   `.match_broad_pathway_term()`.
#' @param lineage_map Optional named character vector mapping exact cell-type
#'   labels to broader biological lineages. Exact mappings override the
#'   default pattern-based assignments.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `frequency`: one row per significant pathway containing:
#'       \itemize{
#'         \item `n_significant`: number of significant cell-type–stratum
#'           contexts.
#'         \item `recurrence_fraction`: proportion of tested contexts in which
#'           the pathway is significant.
#'         \item `n_lineages`: number of distinct biological lineage groups
#'           represented by those significant contexts.
#'         \item `cell_lineages`: semicolon-separated lineage groups in which
#'           the pathway is significant.
#'         \item `broad_category_match`: matched broad pathway category, if any.
#'         \item `is_housekeeping`: logical flag indicating whether the pathway
#'           belongs to one of the curated broad pathway families.
#'       }
#'     \item `significant`: significant pathway-context associations annotated
#'       with pathway recurrence and lineage information.
#'     \item `n_tested_combinations`: number of unique cell-type–stratum
#'       combinations evaluated.
#'   }
#'
#' @export
summarize_pathway_recurrence <- function(
    combined_df,
    padj_cutoff = 0.05,
    broad_terms = .default_broad_pathway_terms,
    lineage_map = NULL
) {
  if (!is.data.frame(combined_df)) {
    stop(
      "summarize_pathway_recurrence: `combined_df` must be a data frame.",
      call. = FALSE
    )
  }

  required_cols <- c(
    "pathway",
    "celltype",
    "stratum",
    "padj"
  )

  missing_cols <- setdiff(
    required_cols,
    names(combined_df)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "summarize_pathway_recurrence: missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(combined_df$padj)) {
    stop(
      "summarize_pathway_recurrence: `padj` must be numeric.",
      call. = FALSE
    )
  }

  if (
    length(padj_cutoff) != 1L ||
      !is.numeric(padj_cutoff) ||
      !is.finite(padj_cutoff) ||
      padj_cutoff <= 0 ||
      padj_cutoff >= 1
  ) {
    stop(
      "summarize_pathway_recurrence: `padj_cutoff` must be one ",
      "finite numeric value between 0 and 1.",
      call. = FALSE
    )
  }

  valid_context <- !is.na(combined_df$celltype) &
    nzchar(as.character(combined_df$celltype)) &
    !is.na(combined_df$stratum) &
    nzchar(as.character(combined_df$stratum))

  context_ids <- paste(
    combined_df$celltype[valid_context],
    combined_df$stratum[valid_context],
    sep = "::"
  )

  n_tested_combinations <- length(
    unique(context_ids)
  )

  significant <- combined_df[
    !is.na(combined_df$padj) &
      is.finite(combined_df$padj) &
      combined_df$padj < padj_cutoff &
      !is.na(combined_df$pathway) &
      nzchar(as.character(combined_df$pathway)) &
      valid_context,
    ,
    drop = FALSE
  ]

  frequency_template <- data.frame(
    pathway = character(),
    n_significant = integer(),
    recurrence_fraction = numeric(),
    n_lineages = integer(),
    cell_lineages = character(),
    broad_category_match = character(),
    is_housekeeping = logical(),
    stringsAsFactors = FALSE
  )

  if (nrow(significant) == 0L) {
    significant$cell_lineage <- character(0)
    significant$n_significant <- integer(0)
    significant$recurrence_fraction <- numeric(0)
    significant$n_lineages <- integer(0)
    significant$cell_lineages <- character(0)
    significant$broad_category_match <- character(0)
    significant$is_housekeeping <- logical(0)

    return(
      list(
        frequency = frequency_template,
        significant = significant,
        n_tested_combinations = n_tested_combinations
      )
    )
  }

  # Keep the most significant result when duplicate pathway-context rows exist.
  significant <- significant[
    order(
      significant$pathway,
      significant$celltype,
      significant$stratum,
      significant$padj,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  association_key <- paste(
    significant$pathway,
    significant$celltype,
    significant$stratum,
    sep = "::"
  )

  significant <- significant[
    !duplicated(association_key),
    ,
    drop = FALSE
  ]

  significant$cell_lineage <- .infer_cell_lineage(
    significant$celltype,
    lineage_map = lineage_map
  )

  # Count significant cell-type-stratum contexts for each pathway.
  recurrence_counts <- aggregate(
    x = list(
      n_significant = rep.int(
        1L,
        nrow(significant)
      )
    ),
    by = list(
      pathway = significant$pathway
    ),
    FUN = sum
  )

  # Count the number of unique lineages represented by each pathway.
  lineage_pairs <- unique(
    significant[
      ,
      c(
        "pathway",
        "cell_lineage"
      ),
      drop = FALSE
    ]
  )

  lineage_counts <- aggregate(
    x = list(
      n_lineages = rep.int(
        1L,
        nrow(lineage_pairs)
      )
    ),
    by = list(
      pathway = lineage_pairs$pathway
    ),
    FUN = sum
  )

  # Store the names of all lineages represented by each pathway.
  lineage_names <- aggregate(
    cell_lineage ~ pathway,
    data = lineage_pairs,
    FUN = function(x) {
      paste(
        sort(unique(x)),
        collapse = "; "
      )
    }
  )

  names(lineage_names)[
    names(lineage_names) == "cell_lineage"
  ] <- "cell_lineages"

  recurrence_counts <- merge(
      recurrence_counts,
      lineage_counts,
      by = "pathway",
      all.x = TRUE,
      sort = FALSE
    )

recurrence_counts <- merge(
  recurrence_counts,
  lineage_names,
  by = "pathway",
  all.x = TRUE,
  sort = FALSE
)

  if (n_tested_combinations > 0L) {
    recurrence_counts$recurrence_fraction <-
      recurrence_counts$n_significant /
      n_tested_combinations
  } else {
    recurrence_counts$recurrence_fraction <- NA_real_
  }

  recurrence_counts$broad_category_match <-
    .match_broad_pathway_term(
      recurrence_counts$pathway,
      broad_terms = broad_terms
    )

  recurrence_counts$is_housekeeping <-
    !is.na(recurrence_counts$broad_category_match)

  recurrence_counts <- recurrence_counts[
    order(
      -recurrence_counts$n_significant,
      recurrence_counts$pathway
    ),
    ,
    drop = FALSE
  ]

  # Add pathway-level recurrence and lineage annotations to every
  # significant pathway-context association.
  significant <- merge(
    significant,
    recurrence_counts,
    by = "pathway",
    all.x = TRUE,
    sort = FALSE
  )

  significant <- significant[
    order(
      significant$stratum,
      significant$celltype,
      significant$padj,
      significant$pathway,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  rownames(recurrence_counts) <- NULL
  rownames(significant) <- NULL

  list(
    frequency = recurrence_counts,
    significant = significant,
    n_tested_combinations = n_tested_combinations
  )
}


#' Select low-recurrence pathway associations by stratum
#'
#' Filters significant pathway associations to pathways that:
#'
#' - do not match the manually curated broad pathway categories; and
#' - are significant in no more than `max_frequency` distinct
#'   cell-type–stratum combinations.
#'
#' The top associations within each stratum are selected by adjusted p-value,
#' followed by decreasing absolute normalized enrichment score when `NES`
#' is available.
#'
#' The returned rows represent pathway–cell-type associations, not
#' necessarily unique pathway names.
#'
#' @param recurrence_result Output from
#'   `summarize_pathway_recurrence()`.
#' @param max_frequency Maximum number of distinct cell-type–stratum
#'   combinations in which a pathway may be significant.
#' @param top_n Maximum number of associations retained per stratum.
#'
#' @return Data frame sorted by stratum and ranking criteria. An empty data
#'   frame with the expected columns is returned when no associations pass
#'   the filters.
#'
#' @export
top_low_recurrence_associations_by_stratum <- function(
    recurrence_result,
    max_frequency = 10L,
    top_n = 15L
) {
  if (
    !is.list(recurrence_result) ||
    !"significant" %in% names(recurrence_result) ||
    !is.data.frame(recurrence_result$significant)
  ) {
    stop(
      "top_low_recurrence_associations_by_stratum: ",
      "`recurrence_result` must be the output of ",
      "`summarize_pathway_recurrence()`.",
      call. = FALSE
    )
  }

  if (
    length(max_frequency) != 1L ||
    !is.numeric(max_frequency) ||
    !is.finite(max_frequency) ||
    max_frequency < 1 ||
    max_frequency != floor(max_frequency)
  ) {
    stop(
      "top_low_recurrence_associations_by_stratum: ",
      "`max_frequency` must be one positive integer.",
      call. = FALSE
    )
  }

  if (
    length(top_n) != 1L ||
    !is.numeric(top_n) ||
    !is.finite(top_n) ||
    top_n < 1 ||
    top_n != floor(top_n)
  ) {
    stop(
      "top_low_recurrence_associations_by_stratum: ",
      "`top_n` must be one positive integer.",
      call. = FALSE
    )
  }

  significant <- recurrence_result$significant

required_cols <- c(
  "stratum",
  "pathway",
  "padj",
  "cell_lineage",
  "n_significant",
  "n_lineages",
  "cell_lineages",
  "is_housekeeping"
)

  missing_cols <- setdiff(
    required_cols,
    names(significant)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "top_low_recurrence_associations_by_stratum: ",
      "significant-results table is missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  candidates <- significant[
    !is.na(significant$is_housekeeping) &
      !significant$is_housekeeping &
      !is.na(significant$n_significant) &
      significant$n_significant <= max_frequency,
    ,
    drop = FALSE
  ]

  if (nrow(candidates) == 0L) {
    return(
      candidates
    )
  }

  split_candidates <- split(
    candidates,
    candidates$stratum,
    drop = TRUE
  )

  rank_one_stratum <- function(stratum_df) {
    if ("NES" %in% names(stratum_df)) {
      stratum_df <- stratum_df[
        order(
          stratum_df$padj,
          -abs(stratum_df$NES),
          stratum_df$pathway,
          na.last = TRUE
        ),
        ,
        drop = FALSE
      ]
    } else {
      stratum_df <- stratum_df[
        order(
          stratum_df$padj,
          stratum_df$pathway,
          na.last = TRUE
        ),
        ,
        drop = FALSE
      ]
    }

    utils::head(
      stratum_df,
      top_n
    )
  }

  output <- do.call(
    rbind,
    lapply(
      split_candidates,
      rank_one_stratum
    )
  )

  rownames(output) <- NULL

  output[
    order(
      output$stratum,
      output$padj,
      output$pathway,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
}
