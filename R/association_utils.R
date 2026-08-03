# ============================================================
# Scan-aware association utilities
# ============================================================

#' Correlation effect size with scan-aware inference
#'
#' Reports Pearson or Spearman correlation as the descriptive effect size,
#' while obtaining the p-value from a Gaussian mixed model:
#'
#' `y ~ x + (1 | scan)`
#'
#' For Spearman associations, x and y are rank-transformed before fitting.
#'
#' @param x,y Numeric vectors.
#' @param scan Scan/sample grouping vector.
#' @param method Correlation method.
#' @param min_n Minimum number of complete observations.
#' @param min_scans Minimum number of independent scans.
#'
#' @return Named list containing the correlation, mixed-model slope, p-value,
#'   sample sizes, and model status.
#' @export
scan_aware_association <- function(
    x,
    y,
    scan,
    method = c("spearman", "pearson"),
    min_n = 8L,
    min_scans = 3L
) {
  method <- match.arg(method)

  x <- as.numeric(x)
  y <- as.numeric(y)
  scan <- as.character(scan)

  complete <- is.finite(x) &
    is.finite(y) &
    !is.na(scan) &
    nzchar(scan)

  x <- x[complete]
  y <- y[complete]
  scan <- factor(scan[complete])

  n <- length(x)
  n_scans <- nlevels(
    droplevels(scan)
  )

  result <- list(
    estimate = NA_real_,
    slope = NA_real_,
    slope_se = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_,
    n = n,
    n_scans = n_scans,
    scan_aware = FALSE,
    status = "not_evaluated"
  )

  if (n < min_n) {
    result$status <- "insufficient_rois"
    return(result)
  }

  if (
    !is.finite(stats::sd(x)) ||
      stats::sd(x) <= 0
  ) {
    result$status <- "constant_expression"
    return(result)
  }

  if (
    !is.finite(stats::sd(y)) ||
      stats::sd(y) <= 0
  ) {
    result$status <- "constant_proportion"
    return(result)
  }

  result$estimate <- unname(
    suppressWarnings(
      stats::cor(
        x,
        y,
        method = method
      )
    )
  )

  if (n_scans < min_scans) {
    result$status <- "insufficient_scans"
    return(result)
  }

  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    result$status <- "glmmTMB_unavailable"
    return(result)
  }

  model_x <- x
  model_y <- y

  if (method == "spearman") {
    model_x <- rank(
      model_x,
      ties.method = "average"
    )

    model_y <- rank(
      model_y,
      ties.method = "average"
    )
  }

  model_x <- as.numeric(
    scale(model_x)
  )

  model_y <- as.numeric(
    scale(model_y)
  )

  dat <- data.frame(
    model_y = model_y,
    model_x = model_x,
    scan = scan
  )

  fit <- try(
    glmmTMB::glmmTMB(
      model_y ~ model_x + (1 | scan),
      data = dat,
      family = gaussian()
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    result$status <- "model_failed"
    return(result)
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
    result$status <- "convergence_problem"
    return(result)
  }

  coefficient_table <- summary(
    fit
  )$coefficients$cond

  if (!"model_x" %in% rownames(coefficient_table)) {
    result$status <- "coefficient_missing"
    return(result)
  }

  result$slope <- as.numeric(
    coefficient_table[
      "model_x",
      "Estimate"
    ]
  )

  result$slope_se <- as.numeric(
    coefficient_table[
      "model_x",
      "Std. Error"
    ]
  )

  result$statistic <- as.numeric(
    coefficient_table[
      "model_x",
      "z value"
    ]
  )

  result$p.value <- as.numeric(
    coefficient_table[
      "model_x",
      "Pr(>|z|)"
    ]
  )

  result$scan_aware <- TRUE
  result$status <- "ok"

  result
}
