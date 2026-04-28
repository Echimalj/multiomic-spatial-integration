#' Plotting utilities for spatial analysis
#' @keywords internal
NULL

#' Boxplot by condition
#' @export
plot_celltype_box <- function(df,
                              celltype,
                              x = "disease_status",
                              y = "rel_abundance") {

  library(ggplot2)

  df_sub <- df[df$celltype == celltype, ]

  ggplot(df_sub, aes_string(x = x, y = y)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.4) +
    theme_classic() +
    ggtitle(celltype)
}

#' Dotplot for multiple cell types
#' @export
plot_celltype_dotplot <- function(df,
                                  group_col = "disease_status",
                                  value_col = "rel_abundance") {

  library(ggplot2)
  library(dplyr)

  df_summary <- df %>%
    group_by(celltype, .data[[group_col]]) %>%
    summarize(mean_abundance = mean(.data[[value_col]]), .groups = "drop")

  ggplot(df_summary, aes(x = celltype, y = .data[[group_col]],
                         size = mean_abundance, color = mean_abundance)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
