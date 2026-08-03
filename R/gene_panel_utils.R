# ============================================================
# Gene-panel validation utilities
# ============================================================

#' Load and validate a gene-panel definition
#'
#' @param panel_csv Path to a CSV containing at least a `gene` column.
#'
#' @return Data frame containing one row per unique gene.
#' @export
load_gene_panel <- function(panel_csv) {
  if (
    length(panel_csv) != 1L ||
      !is.character(panel_csv) ||
      !nzchar(panel_csv)
  ) {
    stop(
      "load_gene_panel: `panel_csv` must be one non-empty file path.",
      call. = FALSE
    )
  }

  if (!file.exists(panel_csv)) {
    stop(
      "load_gene_panel: file does not exist: ",
      panel_csv,
      call. = FALSE
    )
  }

  panel <- utils::read.csv(
    panel_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"gene" %in% names(panel)) {
    stop(
      "load_gene_panel: panel must contain a `gene` column.",
      call. = FALSE
    )
  }

  panel$gene <- trimws(
    as.character(panel$gene)
  )

  panel <- panel[
    !is.na(panel$gene) &
      nzchar(panel$gene),
    ,
    drop = FALSE
  ]

  duplicated_genes <- unique(
    panel$gene[
      duplicated(panel$gene)
    ]
  )

  if (length(duplicated_genes) > 0L) {
    warning(
      "load_gene_panel: removing ",
      length(duplicated_genes),
      " duplicated gene(s): ",
      paste(
        utils::head(duplicated_genes, 10),
        collapse = ", "
      ),
      call. = FALSE
    )

    panel <- panel[
      !duplicated(panel$gene),
      ,
      drop = FALSE
    ]
  }

  if (!"panel" %in% names(panel)) {
    panel$panel <- tools::file_path_sans_ext(
      basename(panel_csv)
    )
  }

  panel$panel <- as.character(
    panel$panel
  )

  rownames(panel) <- NULL

  panel
}


#' Build one validated metadata row per ROI
#'
#' @param proportions_long Long-format cell-type proportion table.
#' @param roi_id_col ROI identifier column.
#' @param disease_col Disease-status column.
#' @param pathology_col Pathology-status column.
#' @param region_col Spatial-region column.
#' @param scan_col Biological scan/sample column.
#'
#' @return Data frame with one row per ROI and ROI IDs as row names.
#' @export
build_roi_meta <- function(
    proportions_long,
    roi_id_col = "ROI_ID",
    disease_col = "disease_status",
    pathology_col = "pathology",
    region_col = "region",
    scan_col = "Scan_ID"
) {
  if (!is.data.frame(proportions_long)) {
    stop(
      "build_roi_meta: `proportions_long` must be a data frame.",
      call. = FALSE
    )
  }

  required_cols <- c(
    roi_id_col,
    disease_col,
    pathology_col,
    region_col,
    scan_col
  )

  missing_cols <- setdiff(
    required_cols,
    names(proportions_long)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "build_roi_meta: missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  meta <- unique(
    proportions_long[
      ,
      required_cols,
      drop = FALSE
    ]
  )

  roi_ids <- as.character(
    meta[[roi_id_col]]
  )

  invalid_roi <- is.na(roi_ids) |
    !nzchar(roi_ids)

  if (any(invalid_roi)) {
    stop(
      "build_roi_meta: ROI identifiers contain missing or empty values.",
      call. = FALSE
    )
  }

  duplicated_roi_ids <- unique(
    roi_ids[
      duplicated(roi_ids) |
        duplicated(roi_ids, fromLast = TRUE)
    ]
  )

  if (length(duplicated_roi_ids) > 0L) {
    stop(
      "build_roi_meta: conflicting metadata rows were found for ",
      length(duplicated_roi_ids),
      " ROI(s), including: ",
      paste(
        utils::head(duplicated_roi_ids, 5),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  rownames(meta) <- roi_ids

  meta
}


#' Map fine cell types to broad cell lineages
#'
#' @param celltypes Character vector of fine cell-type names.
#' @param patterns Optional named list mapping regex patterns to lineage labels.
#'
#' @return Named character vector, one lineage per cell type.
#' @export
map_celltype_to_lineage <- function(
    celltypes,
    patterns = NULL
) {
  celltypes <- as.character(
    celltypes
  )

  if (is.null(patterns)) {
    patterns <- list(
      "^Microglia" = "Microglia",
      "^Astrocyte" = "Astrocyte",
      "^ExNeuron" = "ExcitatoryNeuron",
      "^InhNeuron" = "InhibitoryNeuron",
      "^Oligodendrocyte" = "Oligodendrocyte",
      "^OPC" = "OPC",
      "^(Endothelial|SMC|Pericyte|Fibroblast|VLMC)" = "Vascular"
    )
  }

  if (
    !is.list(patterns) ||
      is.null(names(patterns)) ||
      any(!nzchar(names(patterns)))
  ) {
    stop(
      "map_celltype_to_lineage: `patterns` must be a named list.",
      call. = FALSE
    )
  }

  lineage <- rep(
    "Other",
    length(celltypes)
  )

  names(lineage) <- celltypes

  for (pattern in names(patterns)) {
    hit <- grepl(
      pattern,
      celltypes
    )

    lineage[
      hit &
        lineage == "Other"
    ] <- as.character(
      patterns[[pattern]]
    )
  }

  lineage
}

#' Attribute panel genes to reference cell types and lineages
#'
#' For each panel gene present in the reference signature matrix, identifies
#' the dominant fine cell type and dominant broad lineage. Lineage expression
#' is calculated as the mean signature value across fine cell types assigned
#' to that lineage.
#'
#' The tau specificity score ranges from 0 to 1:
#'
#' - values near 0 indicate broad expression across lineages;
#' - values near 1 indicate lineage-restricted expression.
#'
#' @param signature_mat Numeric matrix with genes in rows and cell types in
#'   columns.
#' @param genes Character vector of panel genes.
#' @param lineage_map Optional named character vector mapping cell types to
#'   lineages. When NULL, `map_celltype_to_lineage()` is used.
#'
#' @return A list with:
#'   - `per_lineage`: gene-by-lineage mean expression;
#'   - `summary`: dominant lineage, specificity, and dominant cell type;
#'   - `missing_genes`: requested genes absent from the signature matrix.
#'
#' @export
attribute_genes_to_celltypes <- function(
    signature_mat,
    genes,
    lineage_map = NULL
) {
  if (
    !is.matrix(signature_mat) &&
      !is.data.frame(signature_mat)
  ) {
    stop(
      "attribute_genes_to_celltypes: `signature_mat` must be a matrix ",
      "or data frame.",
      call. = FALSE
    )
  }

  signature_mat <- as.matrix(
    signature_mat
  )

  if (!is.numeric(signature_mat)) {
    stop(
      "attribute_genes_to_celltypes: `signature_mat` must be numeric.",
      call. = FALSE
    )
  }

  if (
    is.null(rownames(signature_mat)) ||
      is.null(colnames(signature_mat))
  ) {
    stop(
      "attribute_genes_to_celltypes: `signature_mat` must have gene row ",
      "names and cell-type column names.",
      call. = FALSE
    )
  }

  if (anyDuplicated(rownames(signature_mat))) {
    stop(
      "attribute_genes_to_celltypes: duplicated gene row names were found.",
      call. = FALSE
    )
  }

  if (anyDuplicated(colnames(signature_mat))) {
    stop(
      "attribute_genes_to_celltypes: duplicated cell-type column names ",
      "were found.",
      call. = FALSE
    )
  }

  genes <- unique(
    trimws(
      as.character(genes)
    )
  )

  genes <- genes[
    !is.na(genes) &
      nzchar(genes)
  ]

  if (length(genes) == 0L) {
    stop(
      "attribute_genes_to_celltypes: `genes` contains no valid gene names.",
      call. = FALSE
    )
  }

  celltypes <- colnames(
    signature_mat
  )

  if (is.null(lineage_map)) {
    lineage_map <- map_celltype_to_lineage(
      celltypes
    )
  }

  if (
    is.null(names(lineage_map)) ||
      anyDuplicated(names(lineage_map))
  ) {
    stop(
      "attribute_genes_to_celltypes: `lineage_map` must be a uniquely ",
      "named character vector.",
      call. = FALSE
    )
  }

  missing_celltypes <- setdiff(
    celltypes,
    names(lineage_map)
  )

  if (length(missing_celltypes) > 0L) {
    stop(
      "attribute_genes_to_celltypes: lineage mapping is missing ",
      length(missing_celltypes),
      " cell type(s), including: ",
      paste(
        utils::head(missing_celltypes, 10),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  lineage_map <- as.character(
    lineage_map[celltypes]
  )

  names(lineage_map) <- celltypes

  present_genes <- intersect(
    genes,
    rownames(signature_mat)
  )

  missing_genes <- setdiff(
    genes,
    present_genes
  )

  if (length(present_genes) == 0L) {
    stop(
      "attribute_genes_to_celltypes: none of the requested panel genes ",
      "are present in the signature matrix.",
      call. = FALSE
    )
  }

  lineages <- sort(
    unique(
      lineage_map
    )
  )

  per_lineage_rows <- vector(
    "list",
    length(present_genes)
  )

  summary_rows <- vector(
    "list",
    length(present_genes)
  )

  names(per_lineage_rows) <- present_genes
  names(summary_rows) <- present_genes

  for (gene_name in present_genes) {
    celltype_expression <- as.numeric(
      signature_mat[
        gene_name,
        celltypes,
        drop = TRUE
      ]
    )

    names(celltype_expression) <- celltypes

    lineage_expression <- vapply(
      lineages,
      function(lineage_name) {
        values <- celltype_expression[
          lineage_map == lineage_name
        ]

        if (all(!is.finite(values))) {
          return(NA_real_)
        }

        mean(
          values,
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    per_lineage_rows[[gene_name]] <- data.frame(
      gene = gene_name,
      lineage = names(lineage_expression),
      mean_expression = as.numeric(lineage_expression),
      stringsAsFactors = FALSE
    )

    finite_lineages <- lineage_expression[
      is.finite(lineage_expression)
    ]

    finite_celltypes <- celltype_expression[
      is.finite(celltype_expression)
    ]

    dominant_lineage <- if (
      length(finite_lineages) > 0L
    ) {
      names(finite_lineages)[
        which.max(finite_lineages)
      ]
    } else {
      NA_character_
    }

    dominant_celltype <- if (
      length(finite_celltypes) > 0L
    ) {
      names(finite_celltypes)[
        which.max(finite_celltypes)
      ]
    } else {
      NA_character_
    }

    max_lineage_expression <- if (
      length(finite_lineages) > 0L
    ) {
      max(finite_lineages)
    } else {
      NA_real_
    }

    tau_specificity <- if (
      length(finite_lineages) > 1L &&
        is.finite(max_lineage_expression) &&
        max_lineage_expression > 0
    ) {
      sum(
        1 - finite_lineages / max_lineage_expression
      ) /
        (length(finite_lineages) - 1)
    } else {
      NA_real_
    }

    summary_rows[[gene_name]] <- data.frame(
      gene = gene_name,
      dominant_lineage = dominant_lineage,
      tau_specificity = tau_specificity,
      dominant_celltype = dominant_celltype,
      stringsAsFactors = FALSE
    )
  }

  per_lineage <- do.call(
    rbind,
    per_lineage_rows
  )

  summary <- do.call(
    rbind,
    summary_rows
  )

  rownames(per_lineage) <- NULL
  rownames(summary) <- NULL

  list(
    per_lineage = per_lineage,
    summary = summary,
    missing_genes = missing_genes
  )
}

#' Add the three-level disease/pathology group
#'
#' Creates:
#' - Control_AmyloidFree
#' - AD_AmyloidFree
#' - AD_Amyloid
#'
#' @keywords internal
.add_gene_panel_group <- function(
    roi_meta,
    disease_col = "disease_status",
    pathology_col = "pathology"
) {
  required_cols <- c(
    disease_col,
    pathology_col
  )

  missing_cols <- setdiff(
    required_cols,
    names(roi_meta)
  )

  if (length(missing_cols) > 0L) {
    stop(
      ".add_gene_panel_group: missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  disease <- as.character(
    roi_meta[[disease_col]]
  )

  pathology <- as.character(
    roi_meta[[pathology_col]]
  )

  valid_disease <- disease %in% c(
    "Control",
    "AD-CAA"
  )

  valid_pathology <- pathology %in% c(
    "AmyloidFree",
    "Amyloid"
  )

  if (any(!valid_disease | !valid_pathology)) {
    stop(
      ".add_gene_panel_group: unrecognized disease/pathology values were found.",
      call. = FALSE
    )
  }

  group <- ifelse(
    disease == "Control",
    "Control_AmyloidFree",
    ifelse(
      pathology == "AmyloidFree",
      "AD_AmyloidFree",
      "AD_Amyloid"
    )
  )

  roi_meta$group <- factor(
    group,
    levels = c(
      "Control_AmyloidFree",
      "AD_AmyloidFree",
      "AD_Amyloid"
    )
  )

  roi_meta
}


#' Test panel-gene expression across spatial microenvironments
#'
#' For each gene, fits a Gaussian mixed model:
#'
#' `expression ~ group * region + (1 | Scan_ID)`
#'
#' and estimates three contrasts separately within Arteries, Capillaries,
#' and Parenchyma:
#'
#' - Disease_effect:
#'   AD amyloid-free minus Control amyloid-free
#' - Amyloid_effect:
#'   AD amyloid-positive minus AD amyloid-free
#' - MaxPathology_effect:
#'   AD amyloid-positive minus Control amyloid-free
#'
#' @param expr_mat ROI x gene normalized-expression matrix.
#' @param roi_meta One-row-per-ROI metadata from `build_roi_meta()`.
#' @param genes Character vector of panel genes.
#' @param roi_id_col ROI identifier column in `roi_meta`.
#' @param disease_col Disease-status column.
#' @param pathology_col Pathology-status column.
#' @param region_col Spatial-region column.
#' @param scan_col Scan/sample column used as the random intercept.
#' @param adjust_scope Multiple-testing correction scope. `"contrast_region"`
#'   adjusts across genes separately within each contrast x region combination.
#'
#' @return Data frame containing gene, region, contrast, estimate, standard
#'   error, confidence interval, p-value, adjusted p-value, ROI count, scan
#'   count, and model status.
#' @export
test_gene_panel_spatial_perturbation <- function(
    expr_mat,
    roi_meta,
    genes,
    roi_id_col = "ROI_ID",
    disease_col = "disease_status",
    pathology_col = "pathology",
    region_col = "region",
    scan_col = "Scan_ID",
    adjust_scope = "contrast_region"
) {
  required_packages <- c(
    "glmmTMB",
    "emmeans"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_packages) > 0L) {
    stop(
      "test_gene_panel_spatial_perturbation: missing package(s): ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  expr_mat <- as.matrix(
    expr_mat
  )

  if (!is.numeric(expr_mat)) {
    stop(
      "test_gene_panel_spatial_perturbation: `expr_mat` must be numeric.",
      call. = FALSE
    )
  }

  if (
    is.null(rownames(expr_mat)) ||
      is.null(colnames(expr_mat))
  ) {
    stop(
      "test_gene_panel_spatial_perturbation: `expr_mat` must have ROI row ",
      "names and gene column names.",
      call. = FALSE
    )
  }

  required_meta_cols <- c(
    roi_id_col,
    disease_col,
    pathology_col,
    region_col,
    scan_col
  )

  missing_meta_cols <- setdiff(
    required_meta_cols,
    names(roi_meta)
  )

  if (length(missing_meta_cols) > 0L) {
    stop(
      "test_gene_panel_spatial_perturbation: missing metadata column(s): ",
      paste(missing_meta_cols, collapse = ", "),
      call. = FALSE
    )
  }

  genes <- unique(
    trimws(
      as.character(genes)
    )
  )

  genes <- genes[
    !is.na(genes) &
      nzchar(genes)
  ]

  present_genes <- intersect(
    genes,
    colnames(expr_mat)
  )

  if (length(present_genes) == 0L) {
    stop(
      "test_gene_panel_spatial_perturbation: none of the requested genes ",
      "are present in the expression matrix.",
      call. = FALSE
    )
  }

  roi_ids <- as.character(
    roi_meta[[roi_id_col]]
  )

  shared_rois <- intersect(
    rownames(expr_mat),
    roi_ids
  )

  if (length(shared_rois) == 0L) {
    stop(
      "test_gene_panel_spatial_perturbation: no shared ROI IDs between ",
      "expression and metadata.",
      call. = FALSE
    )
  }

  meta_index <- match(
    shared_rois,
    roi_ids
  )

  meta <- roi_meta[
    meta_index,
    ,
    drop = FALSE
  ]

  rownames(meta) <- shared_rois

  meta <- .add_gene_panel_group(
    meta,
    disease_col = disease_col,
    pathology_col = pathology_col
  )

  meta[[region_col]] <- factor(
    as.character(meta[[region_col]]),
    levels = c(
      "Arteries",
      "Capillaries",
      "Parenchyma"
    )
  )

  meta[[scan_col]] <- factor(
    meta[[scan_col]]
  )

  formula_text <- paste0(
    "expression_value ~ group * ",
    region_col,
    " + (1 | ",
    scan_col,
    ")"
  )

  model_formula <- stats::as.formula(
    formula_text
  )

  contrast_list <- list(
    Disease_effect = c(
      -1,
      1,
      0
    ),
    Amyloid_effect = c(
      0,
      -1,
      1
    ),
    MaxPathology_effect = c(
      -1,
      0,
      1
    )
  )

  result_rows <- list()
  result_index <- 1L

  for (gene_name in present_genes) {
    dat <- meta

    dat$expression_value <- as.numeric(
      expr_mat[
        shared_rois,
        gene_name
      ]
    )

    finite_rows <- is.finite(
      dat$expression_value
    )

    dat <- dat[
      finite_rows,
      ,
      drop = FALSE
    ]

    n_rois <- nrow(dat)
    n_scans <- length(
      unique(dat[[scan_col]])
    )

    if (
      n_rois < 8L ||
        n_scans < 2L ||
        stats::sd(dat$expression_value) == 0
    ) {
      result_rows[[result_index]] <- data.frame(
        gene = gene_name,
        region = NA_character_,
        contrast = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        statistic = NA_real_,
        df = NA_real_,
        p.value = NA_real_,
        n_rois = n_rois,
        n_scans = n_scans,
        model_status = if (
          stats::sd(dat$expression_value) == 0
        ) {
          "constant_expression"
        } else {
          "insufficient_data"
        },
        stringsAsFactors = FALSE
      )

      result_index <- result_index + 1L
      next
    }

    fit <- try(
      glmmTMB::glmmTMB(
        formula = model_formula,
        data = dat,
        family = gaussian()
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      result_rows[[result_index]] <- data.frame(
        gene = gene_name,
        region = NA_character_,
        contrast = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        statistic = NA_real_,
        df = NA_real_,
        p.value = NA_real_,
        n_rois = n_rois,
        n_scans = n_scans,
        model_status = "model_failed",
        stringsAsFactors = FALSE
      )

      result_index <- result_index + 1L
      next
    }

    emm <- try(
      emmeans::emmeans(
        fit,
        specs = stats::as.formula(
          paste0(
            "~ group | ",
            region_col
          )
        )
      ),
      silent = TRUE
    )

    if (inherits(emm, "try-error")) {
      result_rows[[result_index]] <- data.frame(
        gene = gene_name,
        region = NA_character_,
        contrast = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        statistic = NA_real_,
        df = NA_real_,
        p.value = NA_real_,
        n_rois = n_rois,
        n_scans = n_scans,
        model_status = "emmeans_failed",
        stringsAsFactors = FALSE
      )

      result_index <- result_index + 1L
      next
    }

    contrast_result <- try(
      emmeans::contrast(
        emm,
        method = contrast_list
      ),
      silent = TRUE
    )

    if (inherits(contrast_result, "try-error")) {
      result_rows[[result_index]] <- data.frame(
        gene = gene_name,
        region = NA_character_,
        contrast = NA_character_,
        estimate = NA_real_,
        std_error = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        statistic = NA_real_,
        df = NA_real_,
        p.value = NA_real_,
        n_rois = n_rois,
        n_scans = n_scans,
        model_status = "contrast_failed",
        stringsAsFactors = FALSE
      )

      result_index <- result_index + 1L
      next
    }

    contrast_df <- as.data.frame(
      summary(
        contrast_result,
        infer = c(
          TRUE,
          TRUE
        )
      )
    )

        n_contrasts <- nrow(
      contrast_df
    )

    # Identify the region column returned by emmeans.
    region_output_col <- intersect(
      c(
        region_col,
        "region",
        make.names(region_col)
      ),
      names(contrast_df)
    )

    if (length(region_output_col) == 0L) {
      stop(
        "test_gene_panel_spatial_perturbation: emmeans output for gene `",
        gene_name,
        "` does not contain the expected region column. Columns returned: ",
        paste(
          names(contrast_df),
          collapse = ", "
        ),
        call. = FALSE
      )
    }

    region_output_col <- region_output_col[[1]]

    statistic_values <- if (
      "z.ratio" %in% names(contrast_df)
    ) {
      as.numeric(
        contrast_df$z.ratio
      )
    } else if (
      "t.ratio" %in% names(contrast_df)
    ) {
      as.numeric(
        contrast_df$t.ratio
      )
    } else {
      rep(
        NA_real_,
        n_contrasts
      )
    }

    df_values <- if (
      "df" %in% names(contrast_df)
    ) {
      as.numeric(
        contrast_df$df
      )
    } else {
      rep(
        NA_real_,
        n_contrasts
      )
    }

    conf_low_values <- if (
      "lower.CL" %in% names(contrast_df)
    ) {
      as.numeric(
        contrast_df$lower.CL
      )
    } else if (
      "asymp.LCL" %in% names(contrast_df)
    ) {
      as.numeric(
        contrast_df$asymp.LCL
      )
    } else {
      rep(
        NA_real_,
        n_contrasts
      )
    }

    conf_high_values <- if (
      "upper.CL" %in% names(contrast_df)
    ) {
      as.numeric(
        contrast_df$upper.CL
      )
    } else if (
      "asymp.UCL" %in% names(contrast_df)
    ) {
      as.numeric(
        contrast_df$asymp.UCL
      )
    } else {
      rep(
        NA_real_,
        n_contrasts
      )
    }

    output <- data.frame(
      gene = rep(
        gene_name,
        n_contrasts
      ),
      region = as.character(
        contrast_df[[region_output_col]]
      ),
      contrast = as.character(
        contrast_df$contrast
      ),
      estimate = as.numeric(
        contrast_df$estimate
      ),
      std_error = as.numeric(
        contrast_df$SE
      ),
      conf_low = conf_low_values,
      conf_high = conf_high_values,
      statistic = statistic_values,
      df = df_values,
      p.value = as.numeric(
        contrast_df$p.value
      ),
      n_rois = rep(
        n_rois,
        n_contrasts
      ),
      n_scans = rep(
        n_scans,
        n_contrasts
      ),
      model_status = rep(
        "ok",
        n_contrasts
      ),
      stringsAsFactors = FALSE
    )

    result_rows[[result_index]] <- output
    result_index <- result_index + 1L
  }

  out <- do.call(
    rbind,
    result_rows
  )

  rownames(out) <- NULL

  out$p_adj <- NA_real_

  valid_rows <- out$model_status == "ok" &
    is.finite(out$p.value)

  if (adjust_scope == "contrast_region") {
    split_index <- interaction(
      out$contrast,
      out$region,
      drop = TRUE
    )

    adjusted_groups <- split(
      which(valid_rows),
      split_index[valid_rows]
    )

    for (group_rows in adjusted_groups) {
      out$p_adj[group_rows] <- stats::p.adjust(
        out$p.value[group_rows],
        method = "BH"
      )
    }
  } else if (adjust_scope == "global") {
    out$p_adj[valid_rows] <- stats::p.adjust(
      out$p.value[valid_rows],
      method = "BH"
    )
  } else {
    stop(
      "test_gene_panel_spatial_perturbation: unsupported `adjust_scope`: ",
      adjust_scope,
      call. = FALSE
    )
  }

  out
}

#' Aggregate fine cell-type proportions into broad lineages
#'
#' @param proportions_long Long cell-type proportion table.
#' @param lineage_map Named vector mapping fine cell types to lineages.
#' @param roi_id_col ROI identifier column.
#' @param celltype_col Fine cell-type column.
#' @param proportion_col Proportion column.
#'
#' @return Numeric matrix with ROIs in rows and lineages in columns.
#' @export
build_lineage_proportion_matrix <- function(
    proportions_long,
    lineage_map = NULL,
    roi_id_col = "ROI_ID",
    celltype_col = "celltype",
    proportion_col = "rel_abundance"
) {
  required_cols <- c(
    roi_id_col,
    celltype_col,
    proportion_col
  )

  missing_cols <- setdiff(
    required_cols,
    names(proportions_long)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "build_lineage_proportion_matrix: missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  duplicated_pairs <- duplicated(
    proportions_long[
      c(roi_id_col, celltype_col)
    ]
  )

  if (any(duplicated_pairs)) {
    stop(
      "build_lineage_proportion_matrix: duplicated ROI x celltype rows found.",
      call. = FALSE
    )
  }

  celltypes <- unique(
    as.character(proportions_long[[celltype_col]])
  )

  if (is.null(lineage_map)) {
    lineage_map <- map_celltype_to_lineage(
      celltypes
    )
  }

  missing_mapping <- setdiff(
    celltypes,
    names(lineage_map)
  )

  if (length(missing_mapping) > 0L) {
    stop(
      "build_lineage_proportion_matrix: missing lineage mappings for: ",
      paste(
        utils::head(missing_mapping, 10),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  df <- proportions_long
  df$lineage <- unname(
    lineage_map[
      as.character(df[[celltype_col]])
    ]
  )

  df$proportion_value <- suppressWarnings(
    as.numeric(df[[proportion_col]])
  )

  df <- df[
    !is.na(df$lineage) &
      is.finite(df$proportion_value),
    ,
    drop = FALSE
  ]

  aggregated <- stats::aggregate(
    proportion_value ~ ROI_ID + lineage,
    data = data.frame(
      ROI_ID = as.character(df[[roi_id_col]]),
      lineage = df$lineage,
      proportion_value = df$proportion_value,
      stringsAsFactors = FALSE
    ),
    FUN = sum
  )

  wide <- stats::reshape(
    aggregated,
    idvar = "ROI_ID",
    timevar = "lineage",
    direction = "wide"
  )

  rownames(wide) <- wide$ROI_ID
  wide$ROI_ID <- NULL

  colnames(wide) <- sub(
    "^proportion_value\\.",
    "",
    colnames(wide)
  )

  lineage_mat <- as.matrix(wide)
  storage.mode(lineage_mat) <- "numeric"

  # A missing lineage in one ROI means zero estimated abundance.
  lineage_mat[is.na(lineage_mat)] <- 0

  row_sums <- rowSums(lineage_mat)

  if (any(!is.finite(row_sums) | row_sums <= 0)) {
    stop(
      "build_lineage_proportion_matrix: one or more ROIs have invalid ",
      "total lineage abundance.",
      call. = FALSE
    )
  }

  message(
    "Lineage matrix: ",
    nrow(lineage_mat),
    " ROIs x ",
    ncol(lineage_mat),
    " lineages; row-sum range ",
    paste(
      round(range(row_sums), 4),
      collapse = " to "
    ),
    "."
  )

  lineage_mat
}

#' Test lineage-associated gene perturbation in amyloid-positive ROIs
#'
#' Within AD/CAA ROIs, fits one model per gene and lineage:
#'
#' `expression ~ pathology * lineage_proportion + region + (1 | Scan_ID)`
#'
#' The interaction coefficient tests whether the relationship between gene
#' expression and lineage abundance differs between amyloid-free and
#' amyloid-positive microenvironments.
#'
#' @param expr_mat ROI x gene normalized-expression matrix.
#' @param proportions_long Long-format fine cell-type proportions.
#' @param roi_meta One-row-per-ROI metadata.
#' @param genes Character vector of panel genes.
#' @param lineage_map Optional fine-cell-type to lineage mapping.
#' @param roi_id_col ROI identifier column.
#' @param celltype_col Fine cell-type column.
#' @param proportion_col Proportion column.
#' @param disease_col Disease-status column.
#' @param pathology_col Pathology column.
#' @param region_col Region column.
#' @param scan_col Scan/sample column.
#' @param min_n Minimum complete ROI count.
#' @param min_scans Minimum number of scans.
#'
#' @return One row per gene x lineage model.
#' @export
test_gene_panel_lineage_attribution <- function(
    expr_mat,
    proportions_long,
    roi_meta,
    genes,
    lineage_map = NULL,
    roi_id_col = "ROI_ID",
    celltype_col = "celltype",
    proportion_col = "rel_abundance",
    disease_col = "disease_status",
    pathology_col = "pathology",
    region_col = "region",
    scan_col = "Scan_ID",
    min_n = 20L,
    min_scans = 3L
) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "test_gene_panel_lineage_attribution: package `glmmTMB` is required.",
      call. = FALSE
    )
  }

  expr_mat <- as.matrix(expr_mat)

  if (
    !is.numeric(expr_mat) ||
      is.null(rownames(expr_mat)) ||
      is.null(colnames(expr_mat))
  ) {
    stop(
      "test_gene_panel_lineage_attribution: `expr_mat` must be a numeric ",
      "ROI x gene matrix with row and column names.",
      call. = FALSE
    )
  }

  required_meta <- c(
    roi_id_col,
    disease_col,
    pathology_col,
    region_col,
    scan_col
  )

  missing_meta <- setdiff(
    required_meta,
    names(roi_meta)
  )

  if (length(missing_meta) > 0L) {
    stop(
      "test_gene_panel_lineage_attribution: missing metadata column(s): ",
      paste(missing_meta, collapse = ", "),
      call. = FALSE
    )
  }

  lineage_mat <- build_lineage_proportion_matrix(
    proportions_long = proportions_long,
    lineage_map = lineage_map,
    roi_id_col = roi_id_col,
    celltype_col = celltype_col,
    proportion_col = proportion_col
  )

  genes <- unique(
    trimws(as.character(genes))
  )

  genes <- genes[
    !is.na(genes) &
      nzchar(genes)
  ]

  present_genes <- intersect(
    genes,
    colnames(expr_mat)
  )

  if (length(present_genes) == 0L) {
    stop(
      "test_gene_panel_lineage_attribution: no requested genes are present.",
      call. = FALSE
    )
  }

  ad_roi_ids <- as.character(
    roi_meta[[roi_id_col]][
      roi_meta[[disease_col]] == "AD-CAA"
    ]
  )

  shared_rois <- Reduce(
    intersect,
    list(
      rownames(expr_mat),
      rownames(lineage_mat),
      ad_roi_ids
    )
  )

  if (length(shared_rois) < min_n) {
    stop(
      "test_gene_panel_lineage_attribution: only ",
      length(shared_rois),
      " shared AD/CAA ROIs.",
      call. = FALSE
    )
  }

  meta_index <- match(
    shared_rois,
    roi_meta[[roi_id_col]]
  )

  meta <- roi_meta[
    meta_index,
    ,
    drop = FALSE
  ]

  rownames(meta) <- shared_rois

  pathology <- factor(
    as.character(meta[[pathology_col]]),
    levels = c(
      "AmyloidFree",
      "Amyloid"
    )
  )

  region <- factor(
    as.character(meta[[region_col]]),
    levels = c(
      "Arteries",
      "Capillaries",
      "Parenchyma"
    )
  )

  scan <- factor(
    meta[[scan_col]]
  )

  lineages <- colnames(lineage_mat)
  results <- list()
  result_index <- 1L

  for (gene_name in present_genes) {
    expression_value <- as.numeric(
      expr_mat[
        shared_rois,
        gene_name
      ]
    )

    for (lineage_name in lineages) {
      lineage_proportion <- as.numeric(
        lineage_mat[
          shared_rois,
          lineage_name
        ]
      )

      dat <- data.frame(
        expression_value = expression_value,
        pathology = pathology,
        lineage_proportion = lineage_proportion,
        region = region,
        scan = scan,
        stringsAsFactors = FALSE
      )

      complete <- is.finite(dat$expression_value) &
        is.finite(dat$lineage_proportion) &
        !is.na(dat$pathology) &
        !is.na(dat$region) &
        !is.na(dat$scan)

      dat <- dat[
        complete,
        ,
        drop = FALSE
      ]

      mean_lineage_abundance <- mean(
        dat$lineage_proportion,
        na.rm = TRUE
      )
      
      sd_lineage_abundance <- stats::sd(
        dat$lineage_proportion,
        na.rm = TRUE
      )

      n_nonzero_rois <- sum(
        dat$lineage_proportion > 0,
        na.rm = TRUE
      )

      dat$lineage_z <- if (
        is.finite(sd_lineage_abundance) &&
          sd_lineage_abundance > 0
      ) {
        as.numeric(
          scale(dat$lineage_proportion)
        )
      } else {
        rep(
          NA_real_,
          nrow(dat)
        )
      }

      n_rois <- nrow(dat)
      n_scans <- nlevels(
        droplevels(dat$scan)
      )

      status <- "ok"

      if (n_rois < min_n) {
        status <- "insufficient_rois"
      } else if (n_scans < min_scans) {
        status <- "insufficient_scans"
      } else if (stats::sd(dat$expression_value) == 0) {
        status <- "constant_expression"
      } else if (
        !is.finite(sd_lineage_abundance) ||
          sd_lineage_abundance <= 0
      ) {
        status <- "constant_lineage_proportion"
      }

      if (status != "ok") {
        results[[result_index]] <- data.frame(
          gene = gene_name,
          lineage = lineage_name,
          mean_lineage_abundance = mean_lineage_abundance,
          sd_lineage_abundance = sd_lineage_abundance,
          n_nonzero_rois = n_nonzero_rois,
          baseline_slope = NA_real_,
          baseline_slope_se = NA_real_,
          amyloid_interaction = NA_real_,
          interaction_se = NA_real_,
          interaction_z = NA_real_,
          p.value = NA_real_,
          amyloid_slope = NA_real_,
          n_rois = n_rois,
          n_scans = n_scans,
          model_status = status,
          stringsAsFactors = FALSE
        )

        result_index <- result_index + 1L
        next
      }

      fit <- try(
        glmmTMB::glmmTMB(
          expression_value ~
            pathology * lineage_z +
            region +
            (1 | scan),
          data = dat,
          family = gaussian()
        ),
        silent = TRUE
      )

      if (inherits(fit, "try-error")) {
        results[[result_index]] <- data.frame(
          gene = gene_name,
          lineage = lineage_name,
          mean_lineage_abundance = mean_lineage_abundance,
          sd_lineage_abundance = sd_lineage_abundance,
          n_nonzero_rois = n_nonzero_rois,
          baseline_slope = NA_real_,
          baseline_slope_se = NA_real_,
          amyloid_interaction = NA_real_,
          interaction_se = NA_real_,
          interaction_z = NA_real_,
          p.value = NA_real_,
          amyloid_slope = NA_real_,
          n_rois = n_rois,
          n_scans = n_scans,
          model_status = "model_failed",
          stringsAsFactors = FALSE
        )

        result_index <- result_index + 1L
        next
      }
      
	  fit_convergence <- tryCatch(
	    fit$fit$convergence,
	    error = function(e) NA_integer_
	  )

	  positive_definite_hessian <- tryCatch(
	    isTRUE(fit$sdr$pdHess),
	    error = function(e) FALSE
	  )

	  if (
	    !isTRUE(fit_convergence == 0L) ||
	      !positive_definite_hessian
	  ) {
	    results[[result_index]] <- data.frame(
	      gene = gene_name,
	      lineage = lineage_name,
          mean_lineage_abundance = mean_lineage_abundance,
          sd_lineage_abundance = sd_lineage_abundance,
          n_nonzero_rois = n_nonzero_rois,
	      baseline_slope = NA_real_,
	      baseline_slope_se = NA_real_,
	      amyloid_interaction = NA_real_,
	      interaction_se = NA_real_,
	      interaction_z = NA_real_,
	      p.value = NA_real_,
	      amyloid_slope = NA_real_,
	      n_rois = n_rois,
	      n_scans = n_scans,
	      model_status = "convergence_problem",
	      stringsAsFactors = FALSE
	    )

	    result_index <- result_index + 1L
	    next
	  }

      coefficient_table <- summary(fit)$coefficients$cond

      baseline_term <- "lineage_z"
      interaction_term <- paste0(
        "pathologyAmyloid:lineage_z"
      )

      # glmmTMB may reverse the interaction-term order.
      alternate_interaction_term <- paste0(
        "lineage_z:pathologyAmyloid"
      )

      if (
        !interaction_term %in% rownames(coefficient_table) &&
          alternate_interaction_term %in% rownames(coefficient_table)
      ) {
        interaction_term <- alternate_interaction_term
      }

      baseline_slope <- if (
        baseline_term %in% rownames(coefficient_table)
      ) {
        coefficient_table[
          baseline_term,
          "Estimate"
        ]
      } else {
        NA_real_
      }

      baseline_se <- if (
        baseline_term %in% rownames(coefficient_table)
      ) {
        coefficient_table[
          baseline_term,
          "Std. Error"
        ]
      } else {
        NA_real_
      }

      interaction_estimate <- if (
        interaction_term %in% rownames(coefficient_table)
      ) {
        coefficient_table[
          interaction_term,
          "Estimate"
        ]
      } else {
        NA_real_
      }

      interaction_se <- if (
        interaction_term %in% rownames(coefficient_table)
      ) {
        coefficient_table[
          interaction_term,
          "Std. Error"
        ]
      } else {
        NA_real_
      }

      interaction_z <- if (
        interaction_term %in% rownames(coefficient_table)
      ) {
        coefficient_table[
          interaction_term,
          "z value"
        ]
      } else {
        NA_real_
      }

      interaction_p <- if (
        interaction_term %in% rownames(coefficient_table)
      ) {
        coefficient_table[
          interaction_term,
          "Pr(>|z|)"
        ]
      } else {
        NA_real_
      }

      results[[result_index]] <- data.frame(
        gene = gene_name,
        lineage = lineage_name,
        mean_lineage_abundance = mean_lineage_abundance,
        sd_lineage_abundance = sd_lineage_abundance,
        n_nonzero_rois = n_nonzero_rois,
        baseline_slope = as.numeric(baseline_slope),
        baseline_slope_se = as.numeric(baseline_se),
        amyloid_interaction = as.numeric(interaction_estimate),
        interaction_se = as.numeric(interaction_se),
        interaction_z = as.numeric(interaction_z),
        p.value = as.numeric(interaction_p),
        amyloid_slope = as.numeric(
          baseline_slope + interaction_estimate
        ),
              n_rois = n_rois,
        n_scans = n_scans,
        model_status = if (
          is.finite(interaction_p)
        ) {
          "ok"
        } else {
          "interaction_missing"
        },
        stringsAsFactors = FALSE
      )

      result_index <- result_index + 1L
    }
  }

  out <- do.call(
    rbind,
    results
  )

  rownames(out) <- NULL

  valid <- out$model_status == "ok" &
    is.finite(out$p.value)

  out$p_adj_global <- NA_real_
  out$p_adj_within_lineage <- NA_real_
  out$p_adj_within_gene <- NA_real_

  out$p_adj_global[valid] <- stats::p.adjust(
    out$p.value[valid],
    method = "BH"
  )

  lineage_groups <- split(
    which(valid),
    out$lineage[valid]
  )

  for (rows in lineage_groups) {
    out$p_adj_within_lineage[rows] <- stats::p.adjust(
      out$p.value[rows],
      method = "BH"
    )
  }

  gene_groups <- split(
    which(valid),
    out$gene[valid]
  )

  for (rows in gene_groups) {
    out$p_adj_within_gene[rows] <- stats::p.adjust(
      out$p.value[rows],
      method = "BH"
    )
  }

  out
}

#' Correlate panel-gene expression with fine cell-type abundance
#'
#' Within AD/CAA ROIs, tests each panel gene against every fine cell-type
#' proportion separately within amyloid-free and amyloid-positive strata.
#'
#' Spearman correlation is reported as the descriptive effect size. The
#' p-value is obtained from a scan-aware mixed model through
#' `scan_aware_association()`.
#'
#' @param expr_mat ROI x gene normalized-expression matrix.
#' @param proportions_long Long-format fine cell-type proportion table.
#' @param roi_meta One-row-per-ROI metadata.
#' @param genes Character vector of panel genes.
#' @param lineage_map Optional named mapping of fine cell types to lineages.
#' @param roi_id_col ROI identifier column.
#' @param celltype_col Fine cell-type column.
#' @param proportion_col Proportion column.
#' @param disease_col Disease-status column.
#' @param pathology_col Pathology-status column used for stratification.
#' @param scan_col Scan/sample column.
#' @param disease_subset Disease group to retain. Defaults to `"AD-CAA"`.
#' @param method Correlation method.
#' @param min_n Minimum complete ROI count per test.
#' @param min_scans Minimum independent scan count.
#'
#' @return One row per gene x cell type x pathology stratum.
#' @export
correlate_gene_panel_with_celltypes <- function(
    expr_mat,
    proportions_long,
    roi_meta,
    genes,
    lineage_map = NULL,
    roi_id_col = "ROI_ID",
    celltype_col = "celltype",
    proportion_col = "rel_abundance",
    disease_col = "disease_status",
    pathology_col = "pathology",
    scan_col = "Scan_ID",
    disease_subset = "AD-CAA",
    method = "spearman",
    min_n = 8L,
    min_scans = 3L
) {
  if (!exists("scan_aware_association", mode = "function")) {
    stop(
      "correlate_gene_panel_with_celltypes: ",
      "`scan_aware_association()` is unavailable. ",
      "Source R/association_utils.R first.",
      call. = FALSE
    )
  }

  expr_mat <- as.matrix(
    expr_mat
  )

  if (
    !is.numeric(expr_mat) ||
      is.null(rownames(expr_mat)) ||
      is.null(colnames(expr_mat))
  ) {
    stop(
      "correlate_gene_panel_with_celltypes: `expr_mat` must be a numeric ",
      "ROI x gene matrix with row and column names.",
      call. = FALSE
    )
  }

  required_prop_cols <- c(
    roi_id_col,
    celltype_col,
    proportion_col
  )

  missing_prop_cols <- setdiff(
    required_prop_cols,
    names(proportions_long)
  )

  if (length(missing_prop_cols) > 0L) {
    stop(
      "correlate_gene_panel_with_celltypes: missing proportion column(s): ",
      paste(missing_prop_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_meta_cols <- c(
    roi_id_col,
    disease_col,
    pathology_col,
    scan_col
  )

  missing_meta_cols <- setdiff(
    required_meta_cols,
    names(roi_meta)
  )

  if (length(missing_meta_cols) > 0L) {
    stop(
      "correlate_gene_panel_with_celltypes: missing metadata column(s): ",
      paste(missing_meta_cols, collapse = ", "),
      call. = FALSE
    )
  }

  duplicated_pairs <- duplicated(
    proportions_long[
      c(
        roi_id_col,
        celltype_col
      )
    ]
  )

  if (any(duplicated_pairs)) {
    stop(
      "correlate_gene_panel_with_celltypes: duplicated ROI x celltype rows ",
      "were found.",
      call. = FALSE
    )
  }

  genes <- unique(
    trimws(
      as.character(genes)
    )
  )

  genes <- genes[
    !is.na(genes) &
      nzchar(genes)
  ]

  present_genes <- intersect(
    genes,
    colnames(expr_mat)
  )

  if (length(present_genes) == 0L) {
    stop(
      "correlate_gene_panel_with_celltypes: no requested genes are present ",
      "in the expression matrix.",
      call. = FALSE
    )
  }

  celltypes <- sort(
    unique(
      as.character(
        proportions_long[[celltype_col]]
      )
    )
  )

  if (is.null(lineage_map)) {
    lineage_map <- map_celltype_to_lineage(
      celltypes
    )
  }

  missing_lineages <- setdiff(
    celltypes,
    names(lineage_map)
  )

  if (length(missing_lineages) > 0L) {
    stop(
      "correlate_gene_panel_with_celltypes: missing lineage mapping for: ",
      paste(
        utils::head(missing_lineages, 10),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  meta <- roi_meta[
    roi_meta[[disease_col]] == disease_subset,
    ,
    drop = FALSE
  ]

  if (nrow(meta) == 0L) {
    stop(
      "correlate_gene_panel_with_celltypes: no metadata rows matched ",
      "disease_subset = `",
      disease_subset,
      "`.",
      call. = FALSE
    )
  }

  meta_roi_ids <- as.character(
    meta[[roi_id_col]]
  )

  rownames(meta) <- meta_roi_ids

  strata <- c(
    "AmyloidFree",
    "Amyloid"
  )

  strata <- strata[
    strata %in%
      unique(as.character(meta[[pathology_col]]))
  ]

  result_rows <- list()
  result_index <- 1L

  for (stratum_name in strata) {
    stratum_roi_ids <- as.character(
      meta[[roi_id_col]][
        meta[[pathology_col]] == stratum_name
      ]
    )

    for (celltype_name in celltypes) {
      celltype_df <- proportions_long[
        as.character(proportions_long[[celltype_col]]) ==
          celltype_name,
        c(
          roi_id_col,
          proportion_col
        ),
        drop = FALSE
      ]

      celltype_roi_ids <- as.character(
        celltype_df[[roi_id_col]]
      )

      rownames(celltype_df) <- celltype_roi_ids

      shared_rois <- Reduce(
        intersect,
        list(
          rownames(expr_mat),
          stratum_roi_ids,
          celltype_roi_ids
        )
      )

      proportion_vector <- as.numeric(
        celltype_df[
          shared_rois,
          proportion_col
        ]
      )

      scan_vector <- as.character(
        meta[
          shared_rois,
          scan_col
        ]
      )

      mean_proportion <- mean(
        proportion_vector,
        na.rm = TRUE
      )

      sd_proportion <- stats::sd(
        proportion_vector,
        na.rm = TRUE
      )

      n_nonzero_rois <- sum(
        proportion_vector > 0,
        na.rm = TRUE
      )

      for (gene_name in present_genes) {
        expression_vector <- as.numeric(
          expr_mat[
            shared_rois,
            gene_name
          ]
        )

        association <- scan_aware_association(
          x = expression_vector,
          y = proportion_vector,
          scan = scan_vector,
          method = method,
          min_n = min_n,
          min_scans = min_scans
        )

        result_rows[[result_index]] <- data.frame(
          gene = gene_name,
          celltype = celltype_name,
          lineage = unname(
            lineage_map[[celltype_name]]
          ),
          stratum = stratum_name,
          correlation = association$estimate,
          model_slope = association$slope,
          slope_se = association$slope_se,
          statistic = association$statistic,
          p.value = association$p.value,
          n_rois = association$n,
          n_scans = association$n_scans,
          mean_celltype_proportion = mean_proportion,
          sd_celltype_proportion = sd_proportion,
          n_nonzero_rois = n_nonzero_rois,
          scan_aware = association$scan_aware,
          model_status = association$status,
          stringsAsFactors = FALSE
        )

        result_index <- result_index + 1L
      }
    }
  }

  out <- do.call(
    rbind,
    result_rows
  )

  rownames(out) <- NULL

  valid <- out$model_status == "ok" &
    is.finite(out$p.value)

  out$p_adj_global <- NA_real_
  out$p_adj_within_stratum <- NA_real_
  out$p_adj_within_gene_stratum <- NA_real_

  out$p_adj_global[valid] <- stats::p.adjust(
    out$p.value[valid],
    method = "BH"
  )

  stratum_groups <- split(
    which(valid),
    out$stratum[valid]
  )

  for (rows in stratum_groups) {
    out$p_adj_within_stratum[rows] <- stats::p.adjust(
      out$p.value[rows],
      method = "BH"
    )
  }

  gene_stratum_groups <- split(
    which(valid),
    interaction(
      out$gene[valid],
      out$stratum[valid],
      drop = TRUE
    )
  )

  for (rows in gene_stratum_groups) {
    out$p_adj_within_gene_stratum[rows] <- stats::p.adjust(
      out$p.value[rows],
      method = "BH"
    )
  }

  out
}

#' Build an integrated summary table for a gene-panel analysis
#'
#' Combines reference attribution, spatial perturbation, lineage-interaction,
#' and fine-cell-type correlation results into one row per panel gene.
#'
#' @param panel Gene-panel data frame from `load_gene_panel()`.
#' @param reference_summary Combined reference-attribution summary.
#' @param spatial_results Spatial perturbation results.
#' @param lineage_results Lineage-attributed perturbation results.
#' @param fine_results Fine-cell-type discovery correlation results.
#' @param reference_condition Condition used as the primary reference.
#' @param spatial_fdr_cutoff Adjusted-p threshold for spatial results.
#' @param lineage_fdr_cutoff Global adjusted-p threshold for lineage results.
#' @param fine_fdr_cutoff Global adjusted-p threshold for fine-cell results.
#'
#' @return One-row-per-gene integrated summary data frame.
#' @export
build_gene_panel_summary <- function(
    panel,
    reference_summary,
    spatial_results,
    lineage_results,
    fine_results,
    reference_condition = "Control",
    spatial_fdr_cutoff = 0.05,
    lineage_fdr_cutoff = 0.05,
    fine_fdr_cutoff = 0.05
) {
  required_panel_cols <- c(
    "gene",
    "panel",
    "category",
    "headline",
    "source"
  )

  missing_panel_cols <- setdiff(
    required_panel_cols,
    names(panel)
  )

  if (length(missing_panel_cols) > 0L) {
    stop(
      "build_gene_panel_summary: panel is missing column(s): ",
      paste(missing_panel_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # Reference attribution
  # ----------------------------------------------------------

  reference_primary <- reference_summary |>
    dplyr::filter(
      condition == reference_condition
    ) |>
    dplyr::select(
      gene,
      reference_dominant_lineage = dominant_lineage,
      reference_dominant_celltype = dominant_celltype,
      reference_tau_specificity = tau_specificity
    )

  reference_condition_comparison <- reference_summary |>
    dplyr::select(
      gene,
      condition,
      dominant_lineage,
      dominant_celltype
    ) |>
    tidyr::pivot_wider(
      names_from = condition,
      values_from = c(
        dominant_lineage,
        dominant_celltype
      ),
      names_sep = "_"
    )

  if (
    all(
      c(
        "dominant_lineage_Control",
        "dominant_lineage_AD+CAA"
      ) %in% names(reference_condition_comparison)
    )
  ) {
    reference_condition_comparison <-
      reference_condition_comparison |>
      dplyr::mutate(
        reference_lineage_shift =
          !is.na(dominant_lineage_Control) &
          !is.na(`dominant_lineage_AD+CAA`) &
          dominant_lineage_Control !=
          `dominant_lineage_AD+CAA`
      ) |>
      dplyr::select(
        gene,
        reference_lineage_shift
      )
  } else {
    reference_condition_comparison <- data.frame(
      gene = unique(reference_summary$gene),
      reference_lineage_shift = NA,
      stringsAsFactors = FALSE
    )
  }

  # ----------------------------------------------------------
  # Strongest spatial perturbation per gene
  # ----------------------------------------------------------

  spatial_valid <- spatial_results |>
    dplyr::filter(
      model_status == "ok",
      is.finite(p.value),
      is.finite(estimate)
    )

  spatial_best <- spatial_valid |>
    dplyr::group_by(gene) |>
    dplyr::slice_min(
      order_by = p_adj,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      gene,
      strongest_spatial_region = region,
      strongest_spatial_contrast = contrast,
      strongest_spatial_estimate = estimate,
      strongest_spatial_p = p.value,
      strongest_spatial_fdr = p_adj,
      strongest_spatial_fdr_significant =
        is.finite(p_adj) &
        p_adj < spatial_fdr_cutoff
    )

  spatial_counts <- spatial_valid |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      n_spatial_fdr_significant = sum(
        is.finite(p_adj) &
          p_adj < spatial_fdr_cutoff
      ),
      .groups = "drop"
    )

  # ----------------------------------------------------------
  # Strongest lineage interaction per gene
  # ----------------------------------------------------------

  lineage_valid <- lineage_results |>
    dplyr::filter(
      model_status == "ok",
      is.finite(p.value),
      is.finite(amyloid_interaction)
    )

  lineage_best <- lineage_valid |>
    dplyr::group_by(gene) |>
    dplyr::slice_min(
      order_by = p.value,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      gene,
      strongest_lineage = lineage,
      strongest_lineage_interaction = amyloid_interaction,
      strongest_lineage_nominal_p = p.value,
      strongest_lineage_within_gene_fdr = p_adj_within_gene,
      strongest_lineage_global_fdr = p_adj_global,
      strongest_lineage_global_significant =
        is.finite(p_adj_global) &
        p_adj_global < lineage_fdr_cutoff,
      strongest_lineage_exploratory_nominal =
        p.value < 0.05 &
        (
          !is.finite(p_adj_global) |
            p_adj_global >= lineage_fdr_cutoff
        ),
      strongest_lineage_mean_abundance =
        mean_lineage_abundance,
      strongest_lineage_sd_abundance =
        sd_lineage_abundance
    )

  # ----------------------------------------------------------
  # Strongest fine-cell-type association per stratum
  # ----------------------------------------------------------

  fine_valid <- fine_results |>
    dplyr::filter(
      model_status == "ok",
      is.finite(p.value),
      is.finite(correlation)
    )

  fine_best <- fine_valid |>
    dplyr::group_by(
      gene,
      stratum
    ) |>
    dplyr::slice_min(
      order_by = p_adj_global,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      gene,
      stratum,
      celltype,
      lineage,
      correlation,
      p.value,
      p_adj_global
    ) |>
    tidyr::pivot_wider(
      names_from = stratum,
      values_from = c(
        celltype,
        lineage,
        correlation,
        p.value,
        p_adj_global
      ),
      names_sep = "_"
    )

  fine_counts <- fine_valid |>
    dplyr::group_by(
      gene,
      stratum
    ) |>
    dplyr::summarise(
      n_fine_global_fdr_significant = sum(
        is.finite(p_adj_global) &
          p_adj_global < fine_fdr_cutoff
      ),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = stratum,
      values_from = n_fine_global_fdr_significant,
      names_prefix = "n_fine_fdr_significant_",
      values_fill = 0
    )

  # ----------------------------------------------------------
  # Join all modules
  # ----------------------------------------------------------

  summary_table <- panel |>
    dplyr::select(
      gene,
      panel,
      category,
      headline,
      source
    ) |>
    dplyr::left_join(
      reference_primary,
      by = "gene"
    ) |>
    dplyr::left_join(
      reference_condition_comparison,
      by = "gene"
    ) |>
    dplyr::left_join(
      spatial_best,
      by = "gene"
    ) |>
    dplyr::left_join(
      spatial_counts,
      by = "gene"
    ) |>
    dplyr::left_join(
      lineage_best,
      by = "gene"
    ) |>
    dplyr::left_join(
      fine_best,
      by = "gene"
    ) |>
    dplyr::left_join(
      fine_counts,
      by = "gene"
    )

  # Add explicit evidence-tier labels.
  summary_table <- summary_table |>
    dplyr::mutate(
      spatial_evidence = dplyr::case_when(
        strongest_spatial_fdr_significant ~
          "FDR-significant",
        is.finite(strongest_spatial_p) &
          strongest_spatial_p < 0.05 ~
          "Nominal only",
        is.finite(strongest_spatial_p) ~
          "Not significant",
        TRUE ~
          "Not evaluated"
      ),
      lineage_evidence = dplyr::case_when(
        strongest_lineage_global_significant ~
          "Global FDR-significant",
        strongest_lineage_exploratory_nominal ~
          "Exploratory nominal",
        is.finite(strongest_lineage_nominal_p) ~
          "Not significant",
        TRUE ~
          "Not evaluated"
      ),
      fine_celltype_evidence = dplyr::case_when(
        (
          dplyr::coalesce(
            n_fine_fdr_significant_AmyloidFree,
            0L
          ) +
            dplyr::coalesce(
              n_fine_fdr_significant_Amyloid,
              0L
            )
        ) > 0 ~
          "Global FDR-significant",
        TRUE ~
          "No global FDR association"
      )
    )

  summary_table
}
