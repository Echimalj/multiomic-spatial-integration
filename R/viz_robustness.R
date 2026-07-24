#' Figure for the leave-one-Scan-out robustness check
#'
#' Requires `R/viz_theme.R` to be sourced. Uses ggplot2.
#'
#' @keywords internal
NULL


#' Robustness caterpillar plot for one contrast
#'
#' Reads a `<contrast>_robustness.csv` and plots the full estimate against
#' the leave-one-Scan-out estimate range. Rows whose significance call changes
#' when any one scan is removed are highlighted.
#'
#' @param robustness_csv Path to a `<contrast>_robustness.csv`.
#' @param contrast_label Label for the title and filename.
#' @param output_dir Figure output directory.
#' @param top_n Number of most-variable rows to display.
#' @return The ggplot object invisibly.
#' @export
plot_robustness_caterpillar <- function(
    robustness_csv,
    contrast_label,
    output_dir,
    top_n = 30
) {
  df <- utils::read.csv(
    robustness_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(df) == 0) {
    message(
      "plot_robustness_caterpillar: empty input for ",
      contrast_label,
      "."
    )
    return(invisible(NULL))
  }

  required_cols <- c(
    "region",
    "celltype",
    "full_estimate",
    "min_loo_estimate",
    "max_loo_estimate",
    "max_abs_delta",
    "any_sig_flip"
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(df)
  )

  if (length(missing_cols) > 0) {
    stop(
      "plot_robustness_caterpillar: missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  numeric_cols <- c(
    "full_estimate",
    "min_loo_estimate",
    "max_loo_estimate",
    "max_abs_delta"
  )

  df[numeric_cols] <- lapply(
    df[numeric_cols],
    function(x) suppressWarnings(as.numeric(x))
  )

  df <- df[
    !is.na(df$full_estimate) &
      !is.na(df$min_loo_estimate) &
      !is.na(df$max_loo_estimate) &
      !is.na(df$max_abs_delta),
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    message(
      "plot_robustness_caterpillar: no complete rows for ",
      contrast_label,
      "."
    )
    return(invisible(NULL))
  }

  df <- df[
    order(-df$max_abs_delta),
    ,
    drop = FALSE
  ]

  df <- utils::head(
    df,
    top_n
  )

  df$label <- paste(
    df$region,
    df$celltype,
    sep = " | "
  )

  df <- df[
    order(df$full_estimate),
    ,
    drop = FALSE
  ]

  # Multiple rows may have the same region-celltype label.
  # Factor levels must be unique.
  df$label <- factor(
    df$label,
    levels = unique(df$label)
  )

  df$any_sig_flip <- as.logical(
    df$any_sig_flip
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(y = label)
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        xmin = min_loo_estimate,
        xmax = max_loo_estimate
      ),
      color = "grey65",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = full_estimate,
        color = any_sig_flip
      ),
      size = 2
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "FALSE" = "grey30",
        "TRUE" = "#B2182B"
      ),
      labels = c(
        "FALSE" = "stable",
        "TRUE" = "significance flips"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = paste0(
        contrast_label,
        ": leave-one-Scan-out robustness"
      ),
      subtitle = paste0(
        "Full estimate (point) vs. leave-one-out range (bar); top ",
        top_n,
        " most variable"
      ),
      x = "Contrast estimate",
      y = NULL
    ) +
    theme_analysis()

  save_figure(
    p,
    paste0(
      "robustness_",
      contrast_label
    ),
    output_dir,
    width = 9,
    height = 8
  )

  invisible(p)
}
