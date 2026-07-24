# =============================================================================
# All-cell-type marker-specificity figure
#
# Generates:
#   A. Abundance–abundance Spearman correlation
#   B. Marker-score–marker-score Spearman correlation
#   C. Raw abundance–marker Spearman correlation
#   D. Partial abundance–marker correlation
#   E. Change in self-correlation after partial adjustment
#
# Cell types are grouped by broad biological class and hierarchically clustered
# within each class. The same ordering is used across all heatmaps.
#
# Inputs:
#   results/marker_specificity_diagnostics/
#     abundance_abundance_spearman.csv
#     marker_score_marker_score_spearman.csv
#     abundance_marker_raw_spearman.csv
#     abundance_marker_partial_spearman.csv
#
# Outputs:
#   results/marker_specificity_diagnostics/figures/
#     all_celltype_marker_specificity_grouped.pdf
#     all_celltype_marker_specificity_grouped.png
#     all_celltype_marker_specificity_order.csv
#     all_celltype_marker_specificity_class_boundaries.csv
#     all_celltype_marker_specificity_diagonal.csv
# =============================================================================


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

input_dir <- "results/marker_specificity_diagnostics"
output_dir <- file.path(input_dir, "figures")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# Helper: read a correlation matrix
# -----------------------------------------------------------------------------

read_cor_matrix <- function(filename) {
  path <- file.path(input_dir, filename)

  if (!file.exists(path)) {
    stop("Required input file not found: ", path)
  }

  x <- read.csv(
    path,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  x <- as.matrix(x)
  storage.mode(x) <- "numeric"

  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("Matrix is missing row or column names: ", path)
  }

  x
}


# -----------------------------------------------------------------------------
# Load matrices
# -----------------------------------------------------------------------------

abundance_cor <- read_cor_matrix(
  "abundance_abundance_spearman.csv"
)

marker_cor <- read_cor_matrix(
  "marker_score_marker_score_spearman.csv"
)

raw_cor <- read_cor_matrix(
  "abundance_marker_raw_spearman.csv"
)

partial_cor <- read_cor_matrix(
  "abundance_marker_partial_spearman.csv"
)


# -----------------------------------------------------------------------------
# Identify shared cell types
# -----------------------------------------------------------------------------

common_celltypes <- Reduce(
  intersect,
  list(
    rownames(abundance_cor),
    colnames(abundance_cor),
    rownames(marker_cor),
    colnames(marker_cor),
    rownames(raw_cor),
    colnames(raw_cor),
    rownames(partial_cor),
    colnames(partial_cor)
  )
)

if (length(common_celltypes) == 0) {
  stop("No shared cell types were found across the four matrices.")
}

cat("Shared cell types:", length(common_celltypes), "\n")

if (length(common_celltypes) != 46) {
  warning(
    "Expected 46 shared cell types but found ",
    length(common_celltypes),
    "."
  )
}

abundance_cor <- abundance_cor[
  common_celltypes,
  common_celltypes,
  drop = FALSE
]

marker_cor <- marker_cor[
  common_celltypes,
  common_celltypes,
  drop = FALSE
]

raw_cor <- raw_cor[
  common_celltypes,
  common_celltypes,
  drop = FALSE
]

partial_cor <- partial_cor[
  common_celltypes,
  common_celltypes,
  drop = FALSE
]


# -----------------------------------------------------------------------------
# Broad biological classes
# -----------------------------------------------------------------------------

get_cell_class <- function(celltype) {
  if (grepl("^ExNeuron", celltype)) {
    return("Excitatory neurons")
  }

  if (grepl("^InhNeuron", celltype)) {
    return("Inhibitory neurons")
  }

  if (grepl("^Astrocytes", celltype)) {
    return("Astrocytes")
  }

  if (grepl("^Microglia", celltype)) {
    return("Microglia")
  }

  if (grepl("^Oligodendrocytes", celltype)) {
    return("Oligodendrocytes")
  }

  if (grepl("^OPC", celltype)) {
    return("OPCs")
  }

  if (celltype == "Endothelial") {
    return("Endothelial")
  }

  if (celltype == "Pericytes") {
    return("Pericytes")
  }

  if (celltype == "SMC") {
    return("SMC")
  }

  if (grepl("^VLMC", celltype)) {
    return("VLMC")
  }

  if (celltype == "Fibroblast") {
    return("Fibroblast")
  }

  "Other"
}


class_levels <- c(
  "Excitatory neurons",
  "Inhibitory neurons",
  "Astrocytes",
  "Microglia",
  "Oligodendrocytes",
  "OPCs",
  "Endothelial",
  "Pericytes",
  "SMC",
  "VLMC",
  "Fibroblast",
  "Other"
)

cell_class <- vapply(
  common_celltypes,
  get_cell_class,
  character(1)
)

names(cell_class) <- common_celltypes


# -----------------------------------------------------------------------------
# Colors used only for broad-class annotation strips
# -----------------------------------------------------------------------------

class_colors <- c(
  "Excitatory neurons" = "#4C78A8",
  "Inhibitory neurons" = "#F58518",
  "Astrocytes" = "#54A24B",
  "Microglia" = "#E45756",
  "Oligodendrocytes" = "#B279A2",
  "OPCs" = "#FF9DA6",
  "Endothelial" = "#9D755D",
  "Pericytes" = "#BAB0AC",
  "SMC" = "#79706E",
  "VLMC" = "#72B7B2",
  "Fibroblast" = "#D4A72C",
  "Other" = "#8C8C8C"
)


# -----------------------------------------------------------------------------
# Cluster cell types within each broad biological class
#
# The clustering distance is based on the average of:
#   abundance–abundance correlation
#   marker-score–marker-score correlation
#
# This preserves class order while grouping biologically similar subtypes.
# -----------------------------------------------------------------------------

cluster_within_class <- function(members) {
  if (length(members) <= 1) {
    return(members)
  }

  abundance_sub <- abundance_cor[
    members,
    members,
    drop = FALSE
  ]

  marker_sub <- marker_cor[
    members,
    members,
    drop = FALSE
  ]

  combined_similarity <- (
    abundance_sub + marker_sub
  ) / 2

  combined_similarity[!is.finite(combined_similarity)] <- 0
  combined_similarity <- (
    combined_similarity + t(combined_similarity)
  ) / 2

  diag(combined_similarity) <- 1

  distance_matrix <- 1 - combined_similarity
  distance_matrix[distance_matrix < 0] <- 0
  diag(distance_matrix) <- 0

  hc <- stats::hclust(
    stats::as.dist(distance_matrix),
    method = "average"
  )

  members[hc$order]
}


celltype_order <- character(0)

for (cls in class_levels) {
  members <- common_celltypes[
    cell_class[common_celltypes] == cls
  ]

  if (length(members) == 0) {
    next
  }

  celltype_order <- c(
    celltype_order,
    cluster_within_class(members)
  )
}

if (!setequal(celltype_order, common_celltypes)) {
  stop(
    "Internal ordering error: clustered order does not contain ",
    "the same cell types as the input matrices."
  )
}


# -----------------------------------------------------------------------------
# Order all matrices identically
# -----------------------------------------------------------------------------

abundance_cor <- abundance_cor[
  celltype_order,
  celltype_order,
  drop = FALSE
]

marker_cor <- marker_cor[
  celltype_order,
  celltype_order,
  drop = FALSE
]

raw_cor <- raw_cor[
  celltype_order,
  celltype_order,
  drop = FALSE
]

partial_cor <- partial_cor[
  celltype_order,
  celltype_order,
  drop = FALSE
]

ordered_classes <- unname(cell_class[celltype_order])


# -----------------------------------------------------------------------------
# Class boundaries and annotation metadata
# -----------------------------------------------------------------------------

class_runs <- rle(ordered_classes)

class_end <- cumsum(class_runs$lengths)
class_start <- c(1, head(class_end, -1) + 1)
class_mid <- (class_start + class_end) / 2
class_boundaries <- head(class_end, -1) + 0.5

class_summary <- data.frame(
  cell_class = class_runs$values,
  start = class_start,
  end = class_end,
  midpoint = class_mid,
  n_celltypes = class_runs$lengths,
  stringsAsFactors = FALSE
)

order_table <- data.frame(
  order_index = seq_along(celltype_order),
  celltype = celltype_order,
  cell_class = ordered_classes,
  stringsAsFactors = FALSE
)

write.csv(
  order_table,
  file.path(
    output_dir,
    "all_celltype_marker_specificity_order.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  class_summary,
  file.path(
    output_dir,
    "all_celltype_marker_specificity_class_boundaries.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# -----------------------------------------------------------------------------
# Diagonal summary
# -----------------------------------------------------------------------------

diagonal_summary <- data.frame(
  order_index = seq_along(celltype_order),
  celltype = celltype_order,
  cell_class = ordered_classes,
  abundance_self = diag(abundance_cor),
  marker_score_self = diag(marker_cor),
  raw_self_rho = diag(raw_cor),
  partial_self_rho = diag(partial_cor),
  stringsAsFactors = FALSE
)

diagonal_summary$self_rho_change <-
  diagonal_summary$partial_self_rho -
  diagonal_summary$raw_self_rho

write.csv(
  diagonal_summary,
  file.path(
    output_dir,
    "all_celltype_marker_specificity_diagonal.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


# -----------------------------------------------------------------------------
# Heatmap color scales
#
# Panels A/B use the full correlation range.
# Panels C/D use a narrower scale so subtype structure is visible.
# -----------------------------------------------------------------------------

palette_function <- grDevices::colorRampPalette(
  c(
    "#2166AC",
    "#67A9CF",
    "#D1E5F0",
    "#F7F7F7",
    "#FDDBC7",
    "#EF8A62",
    "#B2182B"
  )
)

heatmap_palette <- palette_function(200)

limit_ab <- 1
limit_cd <- 0.6

breaks_ab <- seq(
  -limit_ab,
  limit_ab,
  length.out = length(heatmap_palette) + 1
)

breaks_cd <- seq(
  -limit_cd,
  limit_cd,
  length.out = length(heatmap_palette) + 1
)


# -----------------------------------------------------------------------------
# Utility: clamp values to the plotting limits
# -----------------------------------------------------------------------------

clamp_matrix <- function(x, limit) {
  x[x > limit] <- limit
  x[x < -limit] <- -limit
  x
}


# -----------------------------------------------------------------------------
# Utility: draw biological-class annotation strips
# -----------------------------------------------------------------------------

draw_class_strips <- function(n) {
  strip_colors <- unname(
    class_colors[ordered_classes]
  )

  # Top strip
  graphics::rect(
    xleft = seq_len(n) - 0.5,
    ybottom = n + 0.52,
    xright = seq_len(n) + 0.5,
    ytop = n + 1.15,
    col = strip_colors,
    border = NA,
    xpd = NA
  )

  # Left strip
  graphics::rect(
    xleft = 0.5 - 1.15,
    ybottom = seq_len(n) - 0.5,
    xright = 0.5 - 0.52,
    ytop = seq_len(n) + 0.5,
    col = rev(strip_colors),
    border = NA,
    xpd = NA
  )
}


# -----------------------------------------------------------------------------
# Utility: draw a heatmap panel
# -----------------------------------------------------------------------------

draw_heatmap_panel <- function(mat,
                               title,
                               panel_label,
                               plot_limit,
                               plot_breaks,
                               show_x_labels = TRUE,
                               show_y_labels = TRUE) {
  n <- nrow(mat)

  mat <- clamp_matrix(mat, plot_limit)

  # image() places its first matrix row at the bottom, so rows are reversed.
  plot_mat <- mat[n:1, , drop = FALSE]

  graphics::par(
    mar = c(
      if (show_x_labels) 10.5 else 2,
      if (show_y_labels) 10.5 else 2,
      4.0,
      1.0
    ),
    xpd = NA
  )

  graphics::image(
    x = seq_len(n),
    y = seq_len(n),
    z = t(plot_mat),
    col = heatmap_palette,
    breaks = plot_breaks,
    axes = FALSE,
    xlab = "",
    ylab = "",
    useRaster = TRUE,
    asp = 1,
    xlim = c(0.5, n + 0.5),
    ylim = c(0.5, n + 0.5)
  )

  if (show_x_labels) {
    graphics::axis(
      side = 1,
      at = seq_len(n),
      labels = celltype_order,
      las = 2,
      tick = FALSE,
      line = -0.4,
      cex.axis = 0.42
    )
  }

  if (show_y_labels) {
    graphics::axis(
      side = 2,
      at = seq_len(n),
      labels = rev(celltype_order),
      las = 2,
      tick = FALSE,
      line = -0.4,
      cex.axis = 0.42
    )
  }

  # Major-class separators
  for (boundary in class_boundaries) {
    graphics::abline(
      v = boundary,
      lwd = 0.9,
      col = "grey20"
    )

    graphics::abline(
      h = n - boundary + 1,
      lwd = 0.9,
      col = "grey20"
    )
  }

  draw_class_strips(n)

  graphics::box(
    lwd = 0.8,
    col = "grey20"
  )

  graphics::title(
    main = title,
    font.main = 2,
    cex.main = 0.95,
    line = 1.6
  )

  graphics::mtext(
    panel_label,
    side = 3,
    adj = 0,
    line = 2.0,
    font = 2,
    cex = 1.3
  )
}


# -----------------------------------------------------------------------------
# Utility: draw heatmap legend
# -----------------------------------------------------------------------------

draw_heatmap_legend <- function(plot_limit, label) {
  graphics::par(
    mar = c(2.8, 4, 0.7, 4)
  )

  legend_values <- seq(
    -plot_limit,
    plot_limit,
    length.out = length(heatmap_palette)
  )

  graphics::image(
    x = legend_values,
    y = 1,
    z = matrix(
      legend_values,
      nrow = length(legend_values),
      ncol = 1
    ),
    col = heatmap_palette,
    breaks = seq(
      -plot_limit,
      plot_limit,
      length.out = length(heatmap_palette) + 1
    ),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )

  ticks <- c(
    -plot_limit,
    -plot_limit / 2,
    0,
    plot_limit / 2,
    plot_limit
  )

  graphics::axis(
    side = 1,
    at = ticks,
    labels = format(
      ticks,
      trim = TRUE,
      nsmall = 1
    ),
    cex.axis = 0.8
  )

  graphics::mtext(
    label,
    side = 1,
    line = 1.8,
    cex = 0.85
  )

  graphics::box()
}


# -----------------------------------------------------------------------------
# Panel E: change in self-correlation
# -----------------------------------------------------------------------------
draw_change_panel <- function() {
  plot_df <- diagonal_summary[
    order(
      diagonal_summary$self_rho_change,
      decreasing = FALSE
    ),
  ]

  n <- nrow(plot_df)
  y <- seq_len(n)

  positive_color <- "#B2182B"
  negative_color <- "#2166AC"

  bar_colors <- ifelse(
    plot_df$self_rho_change >= 0,
    positive_color,
    negative_color
  )

  graphics::par(
    mar = c(4.5, 10.5, 3.8, 1.5),
    xpd = FALSE
  )

  x_min <- min(
    -0.30,
    min(plot_df$self_rho_change, na.rm = TRUE)
  )

  x_max <- max(
    0.30,
    max(plot_df$self_rho_change, na.rm = TRUE)
  )

  graphics::plot(
    NA,
    xlim = c(x_min, x_max),
    ylim = c(0.5, n + 0.5),
    axes = FALSE,
    xlab = "",
    ylab = "",
    bty = "n"
  )

  grid_ticks <- pretty(c(x_min, x_max))

  graphics::abline(
    v = grid_ticks,
    col = "grey90",
    lwd = 0.7
  )

  graphics::abline(
    v = 0,
    col = "grey30",
    lwd = 1
  )

  graphics::rect(
    xleft = pmin(0, plot_df$self_rho_change),
    ybottom = y - 0.32,
    xright = pmax(0, plot_df$self_rho_change),
    ytop = y + 0.32,
    col = bar_colors,
    border = NA
  )

  graphics::axis(
    side = 1,
    at = grid_ticks,
    cex.axis = 0.8
  )

  graphics::axis(
    side = 2,
    at = y,
    labels = plot_df$celltype,
    las = 2,
    tick = FALSE,
    cex.axis = 0.56
  )

  graphics::mtext(
    "Change in self-correlation (partial - raw)",
    side = 1,
    line = 2.5,
    cex = 0.9
  )

  graphics::title(
    main = "Change in cell-type self-correlation after adjustment",
    font.main = 2,
    cex.main = 1,
    line = 1.4
  )

  graphics::mtext(
    "E",
    side = 3,
    adj = 0,
    line = 1.7,
    font = 2,
    cex = 1.3
  )

  # Label the five largest increases and five largest decreases.
  label_indices <- unique(c(
    head(order(plot_df$self_rho_change), 5),
    head(order(
      plot_df$self_rho_change,
      decreasing = TRUE
    ), 5)
  ))

  for (i in label_indices) {
    value <- plot_df$self_rho_change[i]

    graphics::text(
      x = value,
      y = y[i],
      labels = sprintf("%.2f", value),
      pos = if (value >= 0) 4 else 2,
      offset = 0.35,
      cex = 0.58,
      xpd = NA
    )
  }
}


# -----------------------------------------------------------------------------
# Broad-class legend
# -----------------------------------------------------------------------------

draw_class_legend <- function() {
  graphics::par(
    mar = c(0.5, 1, 0.5, 1)
  )

  graphics::plot.new()

  classes_present <- class_summary$cell_class

  labels <- paste0(
    classes_present,
    " (n=",
    class_summary$n_celltypes,
    ")"
  )

  n_cols <- 4
  n_rows <- ceiling(length(labels) / n_cols)

  for (i in seq_along(labels)) {
    col_index <- (i - 1) %% n_cols
    row_index <- floor((i - 1) / n_cols)

    x <- 0.02 + col_index * 0.245
    y <- 0.84 - row_index * (0.70 / max(1, n_rows - 1))

    graphics::points(
      x = x,
      y = y,
      pch = 15,
      cex = 1.0,
      col = class_colors[classes_present[i]]
    )

    graphics::text(
      x = x + 0.018,
      y = y,
      labels = labels[i],
      adj = c(0, 0.5),
      cex = 0.72
    )
  }
}


# -----------------------------------------------------------------------------
# Build the complete figure
# -----------------------------------------------------------------------------

draw_full_figure <- function() {
  layout_matrix <- matrix(
  c(
    1, 2,
    3, 4,
    5, 6,
    7, 7
  ),
  nrow = 4,
  byrow = TRUE
)

graphics::layout(
  layout_matrix,
  heights = c(
    1,
    1,
    0.09,
    0.65
  )
)

  draw_heatmap_panel(
    abundance_cor,
    title = "Abundance-abundance correlation",
    panel_label = "A",
    plot_limit = limit_ab,
    plot_breaks = breaks_ab
  )

  draw_heatmap_panel(
    marker_cor,
    title = "Marker-score-marker-score correlation",
    panel_label = "B",
    plot_limit = limit_ab,
    plot_breaks = breaks_ab
  )

  draw_heatmap_panel(
    raw_cor,
    title = "Raw abundance-marker correlation",
    panel_label = "C",
    plot_limit = limit_cd,
    plot_breaks = breaks_cd
  )

  draw_heatmap_panel(
    partial_cor,
    title = paste0(
      "Partial abundance-marker correlation\n",
      "adjusted for shared excitatory abundance"
    ),
    panel_label = "D",
    plot_limit = limit_cd,
    plot_breaks = breaks_cd
  )

  draw_heatmap_legend(
    plot_limit = limit_ab,
    label = "Panels A–B: Spearman correlation"
  )

  draw_heatmap_legend(
    plot_limit = limit_cd,
    label = "Panels C–D: Spearman correlation"
  )

  draw_change_panel()

}


# -----------------------------------------------------------------------------
# Save PDF
# -----------------------------------------------------------------------------

pdf_file <- file.path(
  output_dir,
  "all_celltype_marker_specificity_grouped.pdf"
)

grDevices::pdf(
  pdf_file,
  width = 20,
  height = 23,
  onefile = TRUE,
  useDingbats = FALSE
)

draw_full_figure()

grDevices::dev.off()


# -----------------------------------------------------------------------------
# Save high-resolution PNG
# -----------------------------------------------------------------------------

png_file <- file.path(
  output_dir,
  "all_celltype_marker_specificity_grouped.png"
)

grDevices::png(
  png_file,
  width = 6000,
  height = 6900,
  res = 300
)

draw_full_figure()

grDevices::dev.off()


# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------

cat("\nFigure generation completed.\n")
cat("PDF:", pdf_file, "\n")
cat("PNG:", png_file, "\n")

cat(
  "Cell-type order:",
  file.path(
    output_dir,
    "all_celltype_marker_specificity_order.csv"
  ),
  "\n"
)

cat(
  "Class boundaries:",
  file.path(
    output_dir,
    "all_celltype_marker_specificity_class_boundaries.csv"
  ),
  "\n"
)

cat(
  "Diagonal summary:",
  file.path(
    output_dir,
    "all_celltype_marker_specificity_diagonal.csv"
  ),
  "\n"
)

cat("\nSelf-correlation summary:\n")
cat(
  "Positive changes:",
  sum(
    diagonal_summary$self_rho_change > 0,
    na.rm = TRUE
  ),
  "of",
  nrow(diagonal_summary),
  "\n"
)

cat(
  "Negative changes:",
  sum(
    diagonal_summary$self_rho_change < 0,
    na.rm = TRUE
  ),
  "of",
  nrow(diagonal_summary),
  "\n"
)

cat(
  "Median change:",
  median(
    diagonal_summary$self_rho_change,
    na.rm = TRUE
  ),
  "\n"
)
