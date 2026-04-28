#' Contrast utilities for spatial proportion models
#'
#' Utilities for modeling inferred spatial cell-type proportions using
#' beta mixed-effects models.
#'
#' Contrasts included:
#' - Amyloid effect: Amyloid vs AmyloidFree within AD-CAA
#' - Disease effect: AD-CAA vs Control in AmyloidFree ROIs
#' - Weighted overall AD effect: pathology-weighted AD vs Control
#' - Max pathology effect: AD Amyloid vs Control AmyloidFree
#'
#' @keywords internal
NULL


# ============================================================
# Data preparation
# ============================================================

#' Prepare spatial proportion data
#'
#' @param df Long-format spatial cell-type proportion dataframe.
#' @param abundance_col Column containing relative abundance values.
#' @param disease_col Disease-status column.
#' @param pathology_col Pathology column.
#' @param region_col Region column.
#' @param scan_col Scan/sample ID column.
#' @param celltype_col Cell-type column.
#' @param eps Small value for clipping proportions into (0, 1).
#'
#' @return Cleaned dataframe.
#' @export
prepare_spatial_proportion_data <- function(df,
                                            abundance_col = "rel_abundance",
                                            disease_col = "disease_status",
                                            pathology_col = "pathology",
                                            region_col = "region",
                                            scan_col = "Scan_ID",
                                            celltype_col = "celltype",
                                            eps = 1e-6) {
  required_cols <- c(
    abundance_col,
    disease_col,
    pathology_col,
    region_col,
    scan_col,
    celltype_col
  )

  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df[[abundance_col]] <- as.numeric(df[[abundance_col]])

  df[[abundance_col]] <- pmin(
    pmax(df[[abundance_col]], eps),
    1 - eps
  )

  df[[disease_col]] <- factor(
    df[[disease_col]],
    levels = c("Control", "AD-CAA")
  )

  df[[pathology_col]] <- factor(
    df[[pathology_col]],
    levels = c("AmyloidFree", "Amyloid")
  )

  df[[region_col]] <- factor(df[[region_col]])
  df[[scan_col]] <- factor(df[[scan_col]])
  df[[celltype_col]] <- factor(df[[celltype_col]])

  return(df)
}


# ============================================================
# Internal helper
# ============================================================

#' Run beta mixed model safely for one subset
#'
#' @keywords internal
.run_beta_model <- function(tmp,
                            formula,
                            min_rows = 5,
                            min_scans = 2) {
  if (nrow(tmp) < min_rows) return(NULL)
  if (dplyr::n_distinct(tmp$Scan_ID) < min_scans) return(NULL)

  fit <- try(
    glmmTMB::glmmTMB(
      formula,
      data = tmp,
      family = glmmTMB::beta_family(link = "logit")
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) return(NULL)

  return(fit)
}


# ============================================================
# Amyloid effect
# ============================================================

#' Run amyloid effect model
#'
#' Tests Amyloid vs AmyloidFree within AD-CAA only.
#'
#' Model:
#' rel_abundance ~ pathology + (1 | Scan_ID)
#'
#' @param df Cleaned dataframe.
#' @param abundance_col Relative abundance column.
#'
#' @return List with means and contrasts.
#' @export
run_amyloid_effect <- function(df,
                               abundance_col = "rel_abundance") {
  dat_ad <- df |>
    dplyr::filter(.data$disease_status == "AD-CAA")

  mean_results <- list()
  contrast_results <- list()

  for (rg in levels(dat_ad$region)) {
    for (ct in levels(dat_ad$celltype)) {

      tmp <- dat_ad |>
        dplyr::filter(.data$region == rg, .data$celltype == ct) |>
        dplyr::filter(!is.na(.data[[abundance_col]])) |>
        droplevels()

      if (dplyr::n_distinct(tmp$pathology) < 2) next

      fit <- .run_beta_model(
        tmp = tmp,
        formula = stats::as.formula(
          paste0(abundance_col, " ~ pathology + (1 | Scan_ID)")
        )
      )

      if (is.null(fit)) next

      emm <- emmeans::emmeans(fit, ~ pathology, type = "response")

      mean_results[[paste(rg, ct, sep = "__")]] <- as.data.frame(emm) |>
        dplyr::rename(mean_abundance = response) |>
        dplyr::mutate(celltype = ct, region = rg)

      cont <- emmeans::contrast(
        emm,
        method = list(Amyloid_effect = c(-1, 1))
      )

      contrast_results[[paste(rg, ct, sep = "__")]] <-
        as.data.frame(summary(cont, infer = TRUE)) |>
        dplyr::mutate(celltype = ct, region = rg)
    }
  }

  means_df <- dplyr::bind_rows(mean_results)

  contrast_df <- dplyr::bind_rows(contrast_results) |>
    dplyr::group_by(.data$region) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup()

  list(
    means = means_df,
    contrasts = contrast_df
  )
}


# ============================================================
# Disease effect
# ============================================================

#' Run disease effect model
#'
#' Tests AD-CAA vs Control using AmyloidFree ROIs only.
#'
#' Model:
#' rel_abundance ~ disease_status + (1 | Scan_ID)
#'
#' @param df Cleaned dataframe.
#' @param abundance_col Relative abundance column.
#'
#' @return List with means and contrasts.
#' @export
run_disease_effect <- function(df,
                               abundance_col = "rel_abundance") {
  dat_dis <- df |>
    dplyr::filter(.data$pathology == "AmyloidFree")

  mean_results <- list()
  contrast_results <- list()

  for (rg in levels(dat_dis$region)) {
    for (ct in levels(dat_dis$celltype)) {

      tmp <- dat_dis |>
        dplyr::filter(.data$region == rg, .data$celltype == ct) |>
        dplyr::filter(!is.na(.data[[abundance_col]])) |>
        droplevels()

      if (dplyr::n_distinct(tmp$disease_status) < 2) next

      fit <- .run_beta_model(
        tmp = tmp,
        formula = stats::as.formula(
          paste0(abundance_col, " ~ disease_status + (1 | Scan_ID)")
        )
      )

      if (is.null(fit)) next

      emm <- emmeans::emmeans(fit, ~ disease_status, type = "response")

      mean_results[[paste(rg, ct, sep = "__")]] <- as.data.frame(emm) |>
        dplyr::rename(mean_abundance = response) |>
        dplyr::mutate(celltype = ct, region = rg)

      cont <- emmeans::contrast(
        emm,
        method = list(Disease_effect = c(-1, 1))
      )

      contrast_results[[paste(rg, ct, sep = "__")]] <-
        as.data.frame(summary(cont, infer = TRUE)) |>
        dplyr::mutate(celltype = ct, region = rg)
    }
  }

  means_df <- dplyr::bind_rows(mean_results)

  contrast_df <- dplyr::bind_rows(contrast_results) |>
    dplyr::group_by(.data$region) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup()

  list(
    means = means_df,
    contrasts = contrast_df
  )
}


# ============================================================
# Weighted overall effect
# ============================================================

#' Run weighted overall AD vs Control effect
#'
#' Defines three groups:
#' - Control_AmyloidFree
#' - AD_AmyloidFree
#' - AD_Amyloid
#'
#' Then estimates:
#' AD_overall_vs_Control =
#'   ((1 - w) * AD_AmyloidFree + w * AD_Amyloid) - Control_AmyloidFree
#'
#' where w is the region-specific proportion of AD ROIs that are Amyloid.
#'
#' @param df Cleaned dataframe.
#' @param abundance_col Relative abundance column.
#'
#' @return List with means, contrasts, and region weights.
#' @export
run_weighted_overall_effect <- function(df,
                                        abundance_col = "rel_abundance") {
  overall_dat <- df |>
    dplyr::mutate(
      group = dplyr::case_when(
        .data$disease_status == "Control" ~ "Control_AmyloidFree",
        .data$disease_status == "AD-CAA" & .data$pathology == "AmyloidFree" ~ "AD_AmyloidFree",
        .data$disease_status == "AD-CAA" & .data$pathology == "Amyloid" ~ "AD_Amyloid",
        TRUE ~ NA_character_
      ),
      group = factor(
        .data$group,
        levels = c("Control_AmyloidFree", "AD_AmyloidFree", "AD_Amyloid")
      )
    ) |>
    dplyr::filter(!is.na(.data$group))

  mean_results <- list()
  contrast_results <- list()
  weight_results <- list()

  for (rg in levels(overall_dat$region)) {

    w <- overall_dat |>
      dplyr::filter(.data$region == rg, .data$disease_status == "AD-CAA") |>
      dplyr::summarise(w = mean(.data$pathology == "Amyloid")) |>
      dplyr::pull(.data$w)

    weight_results[[rg]] <- data.frame(region = rg, weight = w)

    for (ct in levels(overall_dat$celltype)) {

      tmp <- overall_dat |>
        dplyr::filter(.data$region == rg, .data$celltype == ct) |>
        dplyr::filter(!is.na(.data[[abundance_col]]), is.finite(.data[[abundance_col]])) |>
        droplevels()

      if (dplyr::n_distinct(tmp$group) < 2) next

      fit <- .run_beta_model(
        tmp = tmp,
        formula = stats::as.formula(
          paste0(abundance_col, " ~ group + (1 | Scan_ID)")
        )
      )

      if (is.null(fit)) next

      emm_response <- try(
        emmeans::emmeans(fit, ~ group, type = "response"),
        silent = TRUE
      )

      if (!inherits(emm_response, "try-error")) {
        emm_df <- as.data.frame(summary(emm_response))

        value_col <- intersect(
          c("response", "prob", "emmean"),
          colnames(emm_df)
        )[1]

        if (!is.na(value_col)) {
          mean_results[[paste(rg, ct, sep = "__")]] <- emm_df |>
            dplyr::rename(mean_abundance = dplyr::all_of(value_col)) |>
            dplyr::mutate(celltype = ct, region = rg, weight = w)
        }
      }

      emm_link <- emmeans::emmeans(fit, ~ group)

      contr <- try(
        emmeans::contrast(
          emm_link,
          method = list(
            Disease_effect = c(-1, 1, 0),
            Amyloid_effect = c(0, -1, 1),
            AD_overall_vs_Control = c(-1, 1 - w, w)
          )
        ),
        silent = TRUE
      )

      if (inherits(contr, "try-error")) next

      contrast_results[[paste(rg, ct, sep = "__")]] <-
        as.data.frame(summary(contr, infer = TRUE)) |>
        dplyr::mutate(celltype = ct, region = rg, weight = w)
    }
  }

  means_df <- dplyr::bind_rows(mean_results)

  contrast_df <- dplyr::bind_rows(contrast_results) |>
    dplyr::group_by(.data$region, .data$contrast) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup()

  weights_df <- dplyr::bind_rows(weight_results)

  list(
    means = means_df,
    contrasts = contrast_df,
    weights = weights_df
  )
}


# ============================================================
# Max pathology effect
# ============================================================

#' Run maximal pathology effect
#'
#' Tests AD Amyloid ROIs vs Control AmyloidFree ROIs.
#'
#' Model:
#' rel_abundance ~ group2 + (1 | Scan_ID)
#'
#' @param df Cleaned dataframe.
#' @param abundance_col Relative abundance column.
#'
#' @return List with means, contrasts, and fit status.
#' @export
run_max_pathology_effect <- function(df,
                                     abundance_col = "rel_abundance") {
  comparison_dat <- df |>
    dplyr::mutate(
      group2 = dplyr::case_when(
        .data$disease_status == "Control" & .data$pathology == "AmyloidFree" ~ "Control_AmyloidFree",
        .data$disease_status == "AD-CAA" & .data$pathology == "Amyloid" ~ "AD_Amyloid",
        TRUE ~ NA_character_
      ),
      group2 = factor(
        .data$group2,
        levels = c("Control_AmyloidFree", "AD_Amyloid")
      )
    ) |>
    dplyr::filter(!is.na(.data$group2))

  mean_results <- list()
  contrast_results <- list()
  fit_status <- list()

  for (rg in levels(comparison_dat$region)) {
    for (ct in levels(comparison_dat$celltype)) {

      tmp <- comparison_dat |>
        dplyr::filter(.data$region == rg, .data$celltype == ct) |>
        dplyr::filter(!is.na(.data[[abundance_col]]), is.finite(.data[[abundance_col]])) |>
        droplevels()

      status <- data.frame(
        region = rg,
        celltype = ct,
        n_rows = nrow(tmp),
        n_scans = dplyr::n_distinct(tmp$Scan_ID),
        n_groups = dplyr::n_distinct(tmp$group2),
        fit_ok = FALSE,
        pd_hessian = NA,
        stringsAsFactors = FALSE
      )

      if (nrow(tmp) < 5 || dplyr::n_distinct(tmp$Scan_ID) < 2 ||
          dplyr::n_distinct(tmp$group2) < 2) {
        fit_status[[paste(rg, ct, sep = "__")]] <- status
        next
      }

      fit <- .run_beta_model(
        tmp = tmp,
        formula = stats::as.formula(
          paste0(abundance_col, " ~ group2 + (1 | Scan_ID)")
        )
      )

      if (is.null(fit)) {
        fit_status[[paste(rg, ct, sep = "__")]] <- status
        next
      }

      status$fit_ok <- TRUE
      status$pd_hessian <- isTRUE(fit$sdr$pdHess)

      emm_response <- try(
        emmeans::emmeans(fit, ~ group2, type = "response"),
        silent = TRUE
      )

      if (!inherits(emm_response, "try-error")) {
        emm_df <- as.data.frame(summary(emm_response))

        value_col <- intersect(
          c("response", "prob", "emmean"),
          colnames(emm_df)
        )[1]

        if (!is.na(value_col)) {
          mean_results[[paste(rg, ct, sep = "__")]] <- emm_df |>
            dplyr::rename(mean_abundance = dplyr::all_of(value_col)) |>
            dplyr::mutate(celltype = ct, region = rg)
        }
      }

      emm_link <- emmeans::emmeans(fit, ~ group2)

      contr <- try(
        emmeans::contrast(
          emm_link,
          method = list(AD_Amyloid_vs_Control = c(-1, 1))
        ),
        silent = TRUE
      )

      if (!inherits(contr, "try-error")) {
        contrast_results[[paste(rg, ct, sep = "__")]] <-
          as.data.frame(summary(contr, infer = TRUE)) |>
          dplyr::mutate(celltype = ct, region = rg)
      }

      fit_status[[paste(rg, ct, sep = "__")]] <- status
    }
  }

  means_df <- dplyr::bind_rows(mean_results)

  contrast_df <- dplyr::bind_rows(contrast_results) |>
    dplyr::group_by(.data$region) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p.value, method = "BH")) |>
    dplyr::ungroup()

  fit_status_df <- dplyr::bind_rows(fit_status)

  list(
    means = means_df,
    contrasts = contrast_df,
    fit_status = fit_status_df
  )
}


# ============================================================
# Summary formatting
# ============================================================

#' Format contrast output into unified summary table
#'
#' @param contrast_df Contrast dataframe.
#' @param contrast_type Label for contrast class.
#'
#' @return Standardized contrast summary.
#' @export
format_contrast_summary <- function(contrast_df,
                                    contrast_type) {
  if (nrow(contrast_df) == 0) {
    return(data.frame())
  }

  out <- contrast_df |>
    dplyr::mutate(
      log2_OR = dplyr::case_when(
        "odds.ratio" %in% colnames(contrast_df) ~ log2(.data$odds.ratio),
        "estimate" %in% colnames(contrast_df) ~ .data$estimate / log(2),
        TRUE ~ NA_real_
      ),
      contrast_type = contrast_type
    ) |>
    dplyr::select(
      dplyr::any_of(c(
        "celltype",
        "region",
        "contrast",
        "contrast_type",
        "log2_OR",
        "p.value",
        "p_adj",
        "weight"
      ))
    )

  return(out)
}
