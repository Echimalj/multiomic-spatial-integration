#' Cell-type co-occurrence utilities
#'
#' Cell-type proportions are compositional (they sum to ~1 within an ROI),
#' so naive Pearson correlation between cell-type columns is biased toward
#' negative values purely from that sum-to-one constraint ("closure"),
#' independent of any real biological relationship. These utilities apply a
#' centered-log-ratio (CLR) transform before correlating, which is the
#' standard fix for compositional data.
#'
#' Composition also differs structurally by vascular compartment (e.g. SMCs
#' are near-exclusive to Arteries), so correlations are computed separately
#' within each `region` stratum rather than pooled — pooling would mostly
#' recover "which region is this ROI" rather than genuine co-occurrence.
#' Stratifying further by disease/pathology `group` is a natural extension
#' but is left out here given the smaller per-cell sample sizes that would
#' result.
#'
#' @keywords internal
NULL


#' Centered log-ratio transform
#'
#' Row-wise CLR: for each row (ROI), subtracts the mean of the log-proportions
#' across columns (cell types) from each log-proportion. Removes the
#' sum-to-one constraint that induces spurious negative correlations in raw
#' compositional data.
#'
#' @param prop_mat Matrix or data frame of proportions (rows = ROIs,
#'   columns = cell types).
#' @param eps Small value proportions are clipped to before taking logs.
#'
#' @return Data frame of CLR-transformed values, same shape as `prop_mat`.
#' @export
compute_clr <- function(prop_mat, eps = 1e-6) {
  prop_mat <- as.matrix(prop_mat)
  prop_mat <- pmin(pmax(prop_mat, eps), 1 - eps)

  log_mat <- log(prop_mat)
  row_geom_mean_log <- rowMeans(log_mat)

  clr_mat <- log_mat - row_geom_mean_log

  as.data.frame(clr_mat)
}

#' Add disease/pathology context labels
#'
#' Defines the three observed ROI contexts:
#' Control amyloid-free, AD amyloid-free, and AD amyloid-positive.
#'
#' @keywords internal
.add_cooccurrence_context <- function(df) {
  df |>
    dplyr::mutate(
      context = dplyr::case_when(
        .data$disease_status == "Control" &
          .data$pathology == "AmyloidFree" ~ "Control_AF",

        .data$disease_status == "AD-CAA" &
          .data$pathology == "AmyloidFree" ~ "AD_AF",

        .data$disease_status == "AD-CAA" &
          .data$pathology == "Amyloid" ~ "AD_Abeta",

        TRUE ~ NA_character_
      ),
      context = factor(
        .data$context,
        levels = c(
          "Control_AF",
          "AD_AF",
          "AD_Abeta"
        )
      )
    ) |>
    dplyr::filter(!is.na(.data$context))
}

#' Run CLR-transformed cell-type co-occurrence analysis
#'
#' Pivots the long proportion table to one row per ROI (`Scan_ID`) x
#' `stratify_by` with one column per cell type, applies `compute_clr()`,
#' and computes all pairwise Pearson correlations between cell types within
#' each `stratify_by` stratum.
#'
#' @param df Cleaned dataframe (as returned by `prepare_spatial_proportion_data`).
#' @param abundance_col Relative abundance column.
#' @param stratify_by Column to stratify correlations by. Defaults to
#'   `"region"`.
#' @param min_n Minimum number of complete ROIs required within a stratum
#'   before computing correlations for it.
#'
#' @return Long-format data frame: stratify_by column, celltype_1,
#'   celltype_2, r, p.value, p_adj, n.
#' @export
run_celltype_cooccurrence <- function(df,
                                      abundance_col = "rel_abundance",
                                      stratify_by = "region",
                                      min_n = 3) {
  wide_dat <- df |>
    dplyr::select(dplyr::all_of(c("Scan_ID", stratify_by, "celltype", abundance_col))) |>
    dplyr::filter(
      !is.na(.data[[abundance_col]]),
      !is.na(.data[[stratify_by]]),
      !is.na(.data$Scan_ID),
      !is.na(.data$celltype)
    ) |>
    dplyr::group_by(.data$Scan_ID, .data[[stratify_by]], .data$celltype) |>
    dplyr::summarise(
      abundance = mean(.data[[abundance_col]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(c("Scan_ID", stratify_by)),
      names_from = "celltype",
      values_from = abundance
    )

  celltype_cols <- setdiff(colnames(wide_dat), c("Scan_ID", stratify_by))

  if (length(celltype_cols) < 2) {
    return(data.frame())
  }

  strata <- unique(wide_dat[[stratify_by]])
  result_list <- list()

  for (grp in strata) {

    sub <- wide_dat[wide_dat[[stratify_by]] == grp, celltype_cols, drop = FALSE]

    sub <- sub[, colSums(!is.na(sub)) > 0, drop = FALSE]

    if (nrow(sub) < min_n) next
    if (ncol(sub) < 2) next

    sub[is.na(sub)] <- 0

    clr_mat <- compute_clr(sub)

    pairs <- utils::combn(colnames(clr_mat), 2, simplify = FALSE)

    pair_results <- lapply(pairs, function(pair) {
      ct1 <- pair[1]
      ct2 <- pair[2]

      test <- try(
        stats::cor.test(clr_mat[[ct1]], clr_mat[[ct2]], method = "pearson"),
        silent = TRUE
      )

      if (inherits(test, "try-error")) return(NULL)

      data.frame(
        celltype_1 = ct1,
        celltype_2 = ct2,
        r = unname(test$estimate),
        p.value = test$p.value,
        n = nrow(clr_mat),
        stringsAsFactors = FALSE
      )
    })

    pair_df <- dplyr::bind_rows(pair_results)

    if (nrow(pair_df) == 0) next

    pair_df[[stratify_by]] <- grp
    result_list[[as.character(grp)]] <- pair_df
  }

  cooc_df <- dplyr::bind_rows(result_list)

  if (nrow(cooc_df) == 0) {
    return(cooc_df)
  }

  cooc_df <- cooc_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(stratify_by))) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(c(stratify_by, "celltype_1", "celltype_2", "r", "p.value", "p_adj", "n")))

  cooc_df
}

#' Run CLR co-occurrence analysis by region and disease/pathology context
#'
#' Cell-type proportions are averaged within each independent
#' Scan_ID x region x context before CLR transformation and correlation.
#'
#' @param df Cleaned spatial proportion dataframe.
#' @param abundance_col Relative-abundance column.
#' @param min_n Minimum number of independent Scan_IDs required.
#'
#' @return Long-format pairwise correlation table.
#' @export
run_celltype_cooccurrence_by_context <- function(
    df,
    abundance_col = "rel_abundance",
    min_n = 3
) {
  required_cols <- c(
    "Scan_ID",
    "region",
    "celltype",
    "disease_status",
    "pathology",
    abundance_col
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(df)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  context_df <- .add_cooccurrence_context(df)

  wide_dat <- context_df |>
    dplyr::select(
      .data$Scan_ID,
      .data$region,
      .data$context,
      .data$celltype,
      abundance = dplyr::all_of(abundance_col)
    ) |>
    dplyr::filter(
      !is.na(.data$abundance),
      !is.na(.data$Scan_ID),
      !is.na(.data$region),
      !is.na(.data$context),
      !is.na(.data$celltype)
    ) |>
    dplyr::group_by(
      .data$Scan_ID,
      .data$region,
      .data$context,
      .data$celltype
    ) |>
    dplyr::summarise(
      abundance = mean(
        .data$abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols = c(
        "Scan_ID",
        "region",
        "context"
      ),
      names_from = "celltype",
      values_from = "abundance"
    )

  metadata_cols <- c(
    "Scan_ID",
    "region",
    "context"
  )

  celltype_cols <- setdiff(
    colnames(wide_dat),
    metadata_cols
  )

  if (length(celltype_cols) < 2) {
    return(data.frame())
  }

  strata <- wide_dat |>
    dplyr::distinct(
      .data$region,
      .data$context
    )

  result_list <- list()

  for (k in seq_len(nrow(strata))) {
    region_name <- strata$region[k]
    context_name <- strata$context[k]

    sub <- wide_dat[
      wide_dat$region == region_name &
        wide_dat$context == context_name,
      celltype_cols,
      drop = FALSE
    ]

    sub <- sub[
      ,
      colSums(!is.na(sub)) > 0,
      drop = FALSE
    ]

    n_scans <- nrow(sub)

    if (n_scans < min_n) {
      message(
        "Skipping ",
        region_name,
        " / ",
        context_name,
        ": only ",
        n_scans,
        " independent Scan_ID(s)."
      )
      next
    }

    if (ncol(sub) < 2) {
      next
    }

    sub[is.na(sub)] <- 0

    clr_mat <- compute_clr(sub)

    variable_cols <- vapply(
      clr_mat,
      function(x) {
        stats::sd(x, na.rm = TRUE) > 0
      },
      logical(1)
    )

    clr_mat <- clr_mat[
      ,
      variable_cols,
      drop = FALSE
    ]

    if (ncol(clr_mat) < 2) {
      next
    }

    pairs <- utils::combn(
      colnames(clr_mat),
      2,
      simplify = FALSE
    )

    pair_results <- lapply(
      pairs,
      function(pair) {
        ct1 <- pair[1]
        ct2 <- pair[2]

        test <- try(
          stats::cor.test(
            clr_mat[[ct1]],
            clr_mat[[ct2]],
            method = "pearson"
          ),
          silent = TRUE
        )

        if (inherits(test, "try-error")) {
          return(NULL)
        }

        data.frame(
          region = as.character(region_name),
          context = as.character(context_name),
          celltype_1 = ct1,
          celltype_2 = ct2,
          r = unname(test$estimate),
          p.value = test$p.value,
          n_scans = n_scans,
          stringsAsFactors = FALSE
        )
      }
    )

    pair_df <- dplyr::bind_rows(pair_results)

    if (nrow(pair_df) == 0) {
      next
    }

    pair_df <- pair_df |>
      dplyr::mutate(
        p_adj = stats::p.adjust(
          .data$p.value,
          method = "BH"
        )
      )

    result_list[[
      paste(
        region_name,
        context_name,
        sep = "__"
      )
    ]] <- pair_df
  }

  cooc_df <- dplyr::bind_rows(
    result_list
  )

  if (nrow(cooc_df) == 0) {
    return(cooc_df)
  }

  cooc_df |>
    dplyr::select(
      .data$region,
      .data$context,
      .data$celltype_1,
      .data$celltype_2,
      .data$r,
      .data$p.value,
      .data$p_adj,
      .data$n_scans
    )
}
