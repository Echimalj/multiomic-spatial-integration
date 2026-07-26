#' Figures for the spatial statistics analyses
#'
#' Covers:
#'   - contrast effect-size heatmaps with scan-balanced abundance barplots
#'   - region-heterogeneity tests
#'   - cell-type co-occurrence heatmaps
#'   - optional cell-type co-occurrence networks
#'
#' Requires R/viz_theme.R to be sourced first.
#'
#' @keywords internal
NULL


# ============================================================
# Internal helpers
# ============================================================

.require_columns <- function(df, required, object_name = "data frame") {
  missing_columns <- setdiff(required, colnames(df))

  if (length(missing_columns) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      ". Available columns: ",
      paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


.normalize_label <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", as.character(x)))
}


.filename_safe <- function(x) {
  x <- gsub("[^[:alnum:]_-]+", "_", as.character(x))
  gsub("_+", "_", x)
}


# Save the ComplexHeatmap spatial figure as EPS only.
# Save the ComplexHeatmap spatial figure as EPS only.
.save_complex_heatmap <- function(
    heatmap_object,
    filename,
    output_dir,
    width = 10,
    height = 10,
    annotation_legends = NULL,
    column_title = NULL
) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  eps_path <- file.path(
    output_dir,
    paste0(filename, ".eps")
  )

  grDevices::postscript(
    file = eps_path,
    width = width,
    height = height,
    horizontal = FALSE,
    onefile = FALSE,
    paper = "special",
    family = "Helvetica",
    encoding = "ISOLatin1.enc"
  )

  device_open <- TRUE

  on.exit(
    {
      if (device_open) {
        grDevices::dev.off()
      }
    },
    add = TRUE
  )

  grid::grid.newpage()

  draw_args <- list(
    object = heatmap_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE,
    padding = grid::unit(
      c(10, 8, 8, 8),
      "mm"
    )
  )

  if (!is.null(annotation_legends)) {

  if (inherits(annotation_legends, "Legends")) {
    annotation_legends <- list(annotation_legends)
  }

  draw_args$annotation_legend_list <- annotation_legends
}

  if (!is.null(column_title)) {
    draw_args$column_title <- column_title
    draw_args$column_title_gp <- grid::gpar(
      fontsize = 11,
      fontface = "bold"
    )
  }

  do.call(
    ComplexHeatmap::draw,
    draw_args
  )

  grDevices::dev.off()
  device_open <- FALSE

  message("Saved: ", eps_path)

  invisible(eps_path)
}

.prepare_contrast_heatmap_data <- function(
    df,
    effect_col = "log2_OR"
) {
  if (!"type" %in% colnames(df)) {
    if ("contrast_type" %in% colnames(df)) {
      df$type <- df$contrast_type
    } else {
      stop(
        "The combined contrast table must contain either `type` or ",
        "`contrast_type`.",
        call. = FALSE
      )
    }
  }

  .require_columns(
    df,
    c(
      "celltype",
      "region",
      "type",
      effect_col,
      "p_adj"
    ),
    "Combined spatial contrast table"
  )

  df$type <- as.character(df$type)
  df$celltype <- as.character(df$celltype)
  df$region <- as.character(df$region)

  df[[effect_col]] <- suppressWarnings(
    as.numeric(df[[effect_col]])
  )

  df$p_adj <- suppressWarnings(
    as.numeric(df$p_adj)
  )

  type_key <- .normalize_label(df$type)

  df$type[type_key == "amyloid"] <- "Amyloid"
  df$type[type_key == "disease"] <- "Disease"
  df$type[type_key %in% c("overall", "overalleffect")] <- "Overall"
  df$type[
    type_key %in% c(
      "maxpathology",
      "maximum pathology",
      "maxpathologyeffect"
    )
  ] <- "MaxPathology"

  df <- df[
    !is.na(df$celltype) &
      !is.na(df$region) &
      !is.na(df$type) &
      is.finite(df[[effect_col]]),
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    stop(
      "No valid contrast rows remain after filtering.",
      call. = FALSE
    )
  }

  duplicate_key <- paste(
    df$celltype,
    df$region,
    df$type,
    sep = "||"
  )

  if (anyDuplicated(duplicate_key)) {
    warning(
      "Duplicate celltype-region-type rows were found. ",
      "The first row for each combination will be used."
    )

    df <- df[
      !duplicated(duplicate_key),
      ,
      drop = FALSE
    ]
  }

  df
}


.calculate_scan_balanced_abundance <- function(abundance_df) {
  .require_columns(
    abundance_df,
    c(
      "Scan_ID",
      "region",
      "celltype",
      "disease_status",
      "pathology",
      "rel_abundance"
    ),
    "Cell-type abundance table"
  )

  abundance_df$rel_abundance <- suppressWarnings(
    as.numeric(abundance_df$rel_abundance)
  )

  disease_key <- .normalize_label(
    abundance_df$disease_status
  )

  pathology_key <- .normalize_label(
    abundance_df$pathology
  )

  is_control <- disease_key %in% c(
    "control",
    "ctrl"
  )

  is_ad_caa <- disease_key %in% c(
    "adcaa",
    "ad",
    "case"
  )

  is_amyloid_free <- pathology_key %in% c(
    "amyloidfree",
    "abetafree",
    "amyloidnegative",
    "amyloidneg",
    "negative",
    "none"
  )

  is_amyloid <- pathology_key %in% c(
    "amyloid",
    "amyloidpositive",
    "amyloidpos",
    "abeta",
    "positive"
  )

  abundance_df$abundance_group <- NA_character_

  abundance_df$abundance_group[
    is_control & is_amyloid_free
  ] <- "Control_AF"

  abundance_df$abundance_group[
    is_ad_caa & is_amyloid_free
  ] <- "AD_AF"

  abundance_df$abundance_group[
    is_ad_caa & is_amyloid
  ] <- "AD_Abeta"

  abundance_df <- abundance_df[
    !is.na(abundance_df$abundance_group) &
      !is.na(abundance_df$rel_abundance),
    ,
    drop = FALSE
  ]

  if (nrow(abundance_df) == 0) {
    stop(
      "No rows could be assigned to Control_AF, AD_AF, or AD_Abeta. ",
      "Check disease_status and pathology labels.",
      call. = FALSE
    )
  }

  scan_context_df <- dplyr::summarise(
    dplyr::group_by(
      abundance_df,
      Scan_ID,
      region,
      celltype,
      abundance_group
    ),
    scan_mean_abundance = mean(
      rel_abundance,
      na.rm = TRUE
    ),
    n_rois = dplyr::n(),
    .groups = "drop"
  )

  mean_context_long <- dplyr::summarise(
    dplyr::group_by(
      scan_context_df,
      region,
      celltype,
      abundance_group
    ),
    mean_pct = 100 * mean(
      scan_mean_abundance,
      na.rm = TRUE
    ),
    median_pct = 100 * stats::median(
      scan_mean_abundance,
      na.rm = TRUE
    ),
    sd_pct = 100 * stats::sd(
      scan_mean_abundance,
      na.rm = TRUE
    ),
    n_scans = dplyr::n_distinct(Scan_ID),
    .groups = "drop"
  )

  tidyr::pivot_wider(
    dplyr::select(
      mean_context_long,
      region,
      celltype,
      abundance_group,
      mean_pct
    ),
    names_from = abundance_group,
    values_from = mean_pct
  )
}


.make_abundance_annotation_matrix <- function(
    mean_context_df,
    region_name,
    celltype_order
) {
  required_groups <- c(
    "Control_AF",
    "AD_AF",
    "AD_Abeta"
  )

  region_df <- mean_context_df[
    mean_context_df$region == region_name,
    ,
    drop = FALSE
  ]

  for (group_name in required_groups) {
    if (!group_name %in% colnames(region_df)) {
      region_df[[group_name]] <- NA_real_
    }
  }

  match_index <- match(
    celltype_order,
    region_df$celltype
  )

  anno_mat <- as.matrix(
    region_df[
      match_index,
      required_groups,
      drop = FALSE
    ]
  )

  storage.mode(anno_mat) <- "numeric"
  rownames(anno_mat) <- celltype_order
  colnames(anno_mat) <- c(
    "Control AF",
    "AD AF",
    "AD Abeta"
  )

  anno_mat[is.na(anno_mat)] <- 0

  anno_mat
}


.choose_abundance_breaks <- function(max_pct) {
  if (!is.finite(max_pct) || max_pct <= 0) {
    return(c(0, 1))
  }

  candidates <- c(
    0,
    0.5,
    1,
    2,
    5,
    10,
    25,
    50,
    75,
    100
  )

  breaks <- candidates[
    candidates <= max_pct
  ]

  if (length(breaks) == 0) {
    breaks <- 0
  }

  upper_value <- signif(
    max_pct,
    digits = 2
  )

  if (
    upper_value > max(breaks) &&
      upper_value >= 1.25 * max(breaks)
  ) {
    breaks <- c(
      breaks,
      upper_value
    )
  }

  breaks <- sort(
    unique(breaks)
  )

  # Keep the axis readable: no more than five labels.
  if (length(breaks) > 5) {
    keep <- unique(
      round(
        seq(
          1,
          length(breaks),
          length.out = 5
        )
      )
    )

    breaks <- breaks[keep]
  }

  if (length(breaks) == 1) {
    breaks <- c(
      0,
      max(1, upper_value)
    )
  }

  breaks
}

.celltype_lineage <- function(celltypes){

    lineage <- rep("Other", length(celltypes))

    lineage[grepl("^ExNeuron",celltypes)] <- "Excitatory"
    lineage[grepl("^InNeuron",celltypes)] <- "Inhibitory"
    lineage[grepl("^Astro",celltypes)] <- "Astrocyte"
    lineage[grepl("^Oligo",celltypes)] <- "Oligodendrocyte"
    lineage[grepl("^OPC",celltypes)] <- "OPC"
    lineage[grepl("^Micro",celltypes)] <- "Microglia"
    lineage[grepl("^Endo",celltypes)] <- "Endothelial"
    lineage[grepl("^Pericyte",celltypes)] <- "Pericyte"
    lineage[grepl("^SMC",celltypes)] <- "SMC"
    lineage[grepl("^VLMC",celltypes)] <- "VLMC"
    lineage[grepl("^Fibroblast",celltypes)] <- "Fibroblast"

    lineage
}

# ============================================================
# Contrast heatmap with abundance barplots
# ============================================================

#' Effect-size heatmaps with scan-balanced abundance barplots
#'
#' Generates one EPS figure per region. The left annotation contains grouped
#' horizontal bars for Control amyloid-free, AD amyloid-free, and AD
#' amyloid-positive ROIs. Bar lengths use log1p-transformed percentages for
#' visualization, while the axis labels remain on the original percentage
#' scale.
#'
#' @param combined_csv Path to combined_spatial_contrast_summary.csv.
#' @param output_dir Figure output directory.
#' @param abundance_csv Path to the long-format cell-type abundance table.
#' @param padj_cutoff Largest adjusted p-value receiving one asterisk.
#' @param effect_limit Symmetric display limit for log2 odds ratios.
#' @return A named list of ComplexHeatmap objects, invisibly.
#' @export
plot_contrast_effsize_heatmap <- function(
    combined_csv,
    output_dir,
    abundance_csv =
      "results/cell_proportions/spatial_celltype_proportions_for_R.csv",
    effect_col = "log2_OR",
    effect_label = "log2 OR",
    padj_cutoff = 0.05,
    effect_limit = 2.5
) {
if (
  length(effect_col) != 1 ||
    is.na(effect_col) ||
    !nzchar(effect_col)
) {
  stop(
    "`effect_col` must be one non-empty column name.",
    call. = FALSE
  )
}

if (
  length(effect_label) != 1 ||
    is.na(effect_label) ||
    !nzchar(effect_label)
) {
  stop(
    "`effect_label` must be one non-empty display label.",
    call. = FALSE
  )
}

if (
  length(effect_limit) != 1 ||
    !is.finite(effect_limit) ||
    effect_limit <= 0
) {
  stop(
    "`effect_limit` must be one positive finite number.",
    call. = FALSE
  )
}
  required_packages <- c(
    "ComplexHeatmap",
    "circlize",
    "dplyr",
    "tidyr",
    "grid"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_packages) > 0) {
    stop(
      "plot_contrast_effsize_heatmap requires: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  contrast_df <- utils::read.csv(
    combined_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  abundance_df <- utils::read.csv(
    abundance_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  contrast_df <- .prepare_contrast_heatmap_data(
    contrast_df,
    effect_col = effect_col
  )

  mean_context_df <- .calculate_scan_balanced_abundance(
    abundance_df
  )

  desired_contrasts <- c(
    "Amyloid",
    "Disease",
	"MaxPathology",
    "Overall"
  )

  available_contrasts <- desired_contrasts[
    desired_contrasts %in% unique(contrast_df$type)
  ]

  if (length(available_contrasts) == 0) {
    stop(
      "No expected contrast types were found. Observed types: ",
      paste(sort(unique(contrast_df$type)), collapse = ", "),
      call. = FALSE
    )
  }

  contrast_df$type <- factor(
    contrast_df$type,
    levels = available_contrasts
  )

  preferred_region_order <- c(
    "Parenchyma",
    "Arteries",
    "Capillaries"
  )

  observed_regions <- unique(
    contrast_df$region
  )

  regions <- c(
    preferred_region_order[
      preferred_region_order %in% observed_regions
    ],
    setdiff(
      observed_regions,
      preferred_region_order
    )
  )

  celltype_order <- unique(
    contrast_df$celltype
  )

  effect_colors <- circlize::colorRamp2(
    c(
      -effect_limit,
      0,
      effect_limit
    ),
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )

  # ASCII-only legend labels render reliably in EPS.
  effect_legend_breaks <- c(
    -effect_limit,
    -1,
    0,
    1,
    effect_limit
  )

  effect_legend_labels <- c(
    paste0("-", effect_limit),
    "-1",
    "0",
    "1",
    as.character(effect_limit)
  )

  abundance_fill <- c(
    "Control AF" = "#0072B2",
    "AD AF" = "#009E73",
    "AD Abeta" = "#D55E00"
  )

  abundance_legend <- ComplexHeatmap::Legend(
    title = "ROI context",
    labels = names(abundance_fill),
    legend_gp = grid::gpar(
      fill = unname(abundance_fill),
      col = NA
    ),
    title_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    labels_gp = grid::gpar(
      fontsize = 8
    )
  )

  effect_legend <- ComplexHeatmap::Legend(
    title = effect_label,
    col_fun = effect_colors,
    at = effect_legend_breaks,
    labels = effect_legend_labels,
    title_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    labels_gp = grid::gpar(
      fontsize = 8
    ),
    grid_width = grid::unit(
      4,
      "mm"
    ),
    legend_height = grid::unit(
      4,
    "cm"
    )
  )
# ------------------------------------------------------------
# Custom legends
# ------------------------------------------------------------

significance_legend <- ComplexHeatmap::Legend(
  title = "Significance",
  labels = c(
    "*  FDR < 0.05",
    "**  FDR < 0.01",
    "***  FDR < 0.001"
  ),
  type = "points",
  pch = c(
    "*",
    "**",
    "***"
  ),
  size = grid::unit(
    3.5,
    "mm"
  ),
  legend_gp = grid::gpar(
    col = "grey10",
    fontsize = 9,
    fontface = "bold"
  ),
  title_gp = grid::gpar(
    fontsize = 9,
    fontface = "bold"
  ),
  labels_gp = grid::gpar(
    fontsize = 8
  ),
  gap = grid::unit(
    1.5,
    "mm"
  )
  )

  stacked_legend <- ComplexHeatmap::packLegend(
    effect_legend,
    significance_legend,
    abundance_legend,
    direction = "vertical",
    gap = grid::unit(
      5,
      "mm"
    )
  ) 

  heatmap_list <- list()

  for (region_name in regions) {
    region_df <- contrast_df[
      contrast_df$region == region_name,
      ,
      drop = FALSE
    ]

    if (nrow(region_df) == 0) {
      next
    }

    heat_wide <- tidyr::pivot_wider(
      dplyr::select(
        region_df,
        celltype,
        type,
        dplyr::all_of(effect_col)
      ),
      names_from = type,
      values_from = dplyr::all_of(effect_col)
    )

    region_celltypes <- celltype_order[
      celltype_order %in% heat_wide$celltype
    ]

    heat_wide <- heat_wide[
      match(
        region_celltypes,
        heat_wide$celltype
      ),
      ,
      drop = FALSE
    ]

    for (contrast_name in available_contrasts) {
      if (!contrast_name %in% colnames(heat_wide)) {
        heat_wide[[contrast_name]] <- NA_real_
      }
    }

    effect_mat <- as.matrix(
      heat_wide[
        ,
        available_contrasts,
        drop = FALSE
      ]
    )

    storage.mode(effect_mat) <- "numeric"
    rownames(effect_mat) <- heat_wide$celltype

    effect_mat_display <- pmax(
      pmin(
        effect_mat,
        effect_limit
      ),
      -effect_limit
    )
	
	contrast_display_labels <- c(
  	Amyloid = "Amyloid\n(AD Abeta+ vs Abeta-)",
  	Disease = "Disease\n(AD AF vs Control AF)",
  	MaxPathology = "Maximum pathology\n(AD Abeta+ vs Control AF)",
  	Overall = "Overall\n(weighted)"
	)	

	colnames(effect_mat_display) <- unname(
 	 contrast_display_labels[
 	   colnames(effect_mat_display)
 	 ]
	)

    significance_df <- region_df

    significance_df$sig <- dplyr::case_when(
      is.na(significance_df$p_adj) ~ "",
      significance_df$p_adj < 0.001 ~ "***",
      significance_df$p_adj < 0.01 ~ "**",
      significance_df$p_adj < padj_cutoff ~ "*",
      TRUE ~ ""
    )

    significance_wide <- tidyr::pivot_wider(
      dplyr::select(
        significance_df,
        celltype,
        type,
        sig
      ),
      names_from = type,
      values_from = sig,
      values_fill = ""
    )

    significance_wide <- significance_wide[
      match(
        rownames(effect_mat_display),
        significance_wide$celltype
      ),
      ,
      drop = FALSE
    ]

    for (contrast_name in available_contrasts) {
      if (!contrast_name %in% colnames(significance_wide)) {
        significance_wide[[contrast_name]] <- ""
      }
    }

    sig_mat <- as.matrix(
      significance_wide[
        ,
        available_contrasts,
        drop = FALSE
      ]
    )

    rownames(sig_mat) <- significance_wide$celltype
	colnames(sig_mat) <- colnames(effect_mat_display)

    abundance_mat <- .make_abundance_annotation_matrix(
      mean_context_df = mean_context_df,
      region_name = region_name,
      celltype_order = rownames(effect_mat_display)
    )

    abundance_max <- max(
      abundance_mat,
      na.rm = TRUE
    )

    abundance_breaks <- .choose_abundance_breaks(
      abundance_max
    )

    abundance_mat_display <- log1p(
      abundance_mat
    )

    abundance_annotation <- ComplexHeatmap::rowAnnotation(
      `Scan-balanced mean abundance (%)` =
        ComplexHeatmap::anno_barplot(
          abundance_mat_display,
          beside = TRUE,
          gp = grid::gpar(
            fill = unname(abundance_fill),
            col = NA
          ),
          bar_width = 0.82,
          border = FALSE,
          axis = TRUE,
          axis_param = list(
            side = "top",
            at = log1p(abundance_breaks),
            labels = format(
              abundance_breaks,
              trim = TRUE,
              scientific = FALSE,
              digits = 3
            ),
            labels_rot = 0,
            gp = grid::gpar(
              fontsize = 7
            )
          )
        ),
      annotation_name_side = "top",
      annotation_name_rot = 0,
      annotation_name_gp = grid::gpar(
        fontsize = 9,
        fontface = "bold"
      ),
      width = grid::unit(
        5.6,
        "cm"
      )
    )
 
  heatmap_object <- ComplexHeatmap::Heatmap(
    effect_mat_display,
    name = effect_label,
    col = effect_colors,
    na_col = "grey92",
    show_heatmap_legend = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_side = "right",
    row_names_gp = grid::gpar(
      fontsize = 8
    ),
    column_names_rot = 45,
    column_names_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    rect_gp = grid::gpar(
      col = "grey92",
      lwd = 0.6
    ),
    left_annotation = abundance_annotation,
    width = grid::unit(
      6.2,
      "cm"
    ),
    cell_fun = function(
        j,
        i,
        x,
        y,
        width,
        height,
        fill
    ) {
      sig_label <- sig_mat[i, j]

      if (
        !is.na(sig_label) &&
          nzchar(sig_label)
      ) {
        grid::grid.text(
          sig_label,
          x = x,
          y = y,
          gp = grid::gpar(
            fontsize = 9,
            fontface = "bold",
            col = "grey10"
          )
        )
      }
    }
  )

  region_filename <- paste0(
    "contrast_effsize_heatmap_",
    .filename_safe(region_name)
  )

  .save_complex_heatmap(
    heatmap_object = heatmap_object,
    filename = region_filename,
    output_dir = output_dir,
    width = 12.5,
    height = 10.5,
    annotation_legends = stacked_legend,
    column_title = NULL
  )

    heatmap_list[[region_name]] <- heatmap_object
  }

  invisible(heatmap_list)
}


# ============================================================
# Region heterogeneity
# ============================================================

#' Region-heterogeneity interaction significance plot
#'
#' @param interaction_csv Path to interaction-test results.
#' @param output_dir Figure output directory.
#' @param padj_cutoff Adjusted significance threshold.
#' @return The ggplot object invisibly.
#' @export
plot_region_interaction <- function(
    interaction_csv,
    output_dir,
    padj_cutoff = 0.05
) {
  df <- utils::read.csv(
    interaction_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  .require_columns(
    df,
    "celltype",
    "Region interaction table"
  )

  padj_col <- if (
    "p_adj" %in% colnames(df)
  ) {
    "p_adj"
  } else if (
    "p.value" %in% colnames(df)
  ) {
    "p.value"
  } else {
    stop(
      "Region interaction table contains neither p_adj nor p.value.",
      call. = FALSE
    )
  }

  df[[padj_col]] <- suppressWarnings(
    as.numeric(df[[padj_col]])
  )

  positive_p <- df[[padj_col]][
    is.finite(df[[padj_col]]) &
      df[[padj_col]] > 0
  ]

  smallest_positive <- if (
    length(positive_p) > 0
  ) {
    min(positive_p)
  } else {
    .Machine$double.xmin
  }

  plot_p <- df[[padj_col]]
  plot_p[!is.na(plot_p) & plot_p <= 0] <- smallest_positive

  df$neglog_p <- -log10(plot_p)
  df$significant <- (
    !is.na(df[[padj_col]]) &
      df[[padj_col]] < padj_cutoff
  )

  df <- df[
    order(df$neglog_p),
    ,
    drop = FALSE
  ]

  df$celltype <- factor(
    df$celltype,
    levels = unique(df$celltype)
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = neglog_p,
      y = celltype,
      fill = significant
    )
  ) +
    ggplot2::geom_col(
      width = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "grey50"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "FALSE" = "grey75",
        "TRUE" = "#B2182B"
      ),
      labels = c(
        "FALSE" = "n.s.",
        "TRUE" = paste0(
          "adjusted p < ",
          padj_cutoff
        )
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Region heterogeneity of disease/amyloid effects",
      subtitle = "Group x region interaction per cell type",
      x = expression(
        -log[10]~"(adjusted p)"
      ),
      y = NULL
    ) +
    theme_analysis()

  save_figure(
    p,
    "region_interaction_significance",
    output_dir,
    width = 8,
    height = 9
  )

  invisible(p)
}


# ============================================================
# Co-occurrence heatmap
# ============================================================

#' Cell-type co-occurrence correlation heatmap
#'
#' @param cooccurrence_csv Path to celltype_cooccurrence_by_region.csv.
#' @param output_dir Figure output directory.
#' @return The ggplot object invisibly.
#' @export
plot_cooccurrence_heatmap <- function(
    cooccurrence_csv,
    output_dir
) {
  df <- utils::read.csv(
    cooccurrence_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  .require_columns(
    df,
    c(
      "celltype_1",
      "celltype_2",
      "region",
      "r"
    ),
    "Co-occurrence table"
  )

  df$r <- suppressWarnings(
    as.numeric(df$r)
  )

  swap <- df
  swap$celltype_1 <- df$celltype_2
  swap$celltype_2 <- df$celltype_1

  df_sym <- rbind(
    df,
    swap
  )

  p <- ggplot2::ggplot(
    df_sym,
    ggplot2::aes(
      x = celltype_1,
      y = celltype_2,
      fill = r
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(
      ggplot2::vars(region)
    ) +
    scale_fill_diverging(
      name = "CLR correlation",
      limits = c(-1, 1)
    ) +
    ggplot2::labs(
      title = "Cell-type co-occurrence by region",
      subtitle = "Centered log-ratio correlation",
      x = NULL,
      y = NULL
    ) +
    theme_analysis(
      base_size = 8
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6
      ),
      axis.text.y = ggplot2::element_text(
        size = 6
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 9
	  )
    )

  save_figure(
    p,
    "cooccurrence_heatmap",
    output_dir,
    width = 15,
    height = 6.5
  )

  invisible(p)
}


# ============================================================
# Optional co-occurrence network
# ============================================================

#' Cell-type co-occurrence network of significant pairs
#'
#' @param cooccurrence_csv Path to celltype_cooccurrence_by_region.csv.
#' @param output_dir Figure output directory.
#' @param padj_cutoff Adjusted significance threshold for edges.
#' @return The combined plot invisibly, or NULL when unavailable.
#' @export
plot_cooccurrence_network <- function(
    cooccurrence_csv,
    output_dir,
    padj_cutoff = 0.05
) {
  if (
    !requireNamespace(
      "igraph",
      quietly = TRUE
    ) ||
      !requireNamespace(
        "ggraph",
        quietly = TRUE
      )
  ) {
    message(
      "plot_cooccurrence_network: igraph/ggraph not installed; ",
      "skipping network."
    )

    return(invisible(NULL))
  }

  df <- utils::read.csv(
    cooccurrence_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  .require_columns(
    df,
    c(
      "celltype_1",
      "celltype_2",
      "region",
      "r",
      "p_adj"
    ),
    "Co-occurrence network table"
  )

  df$r <- suppressWarnings(
    as.numeric(df$r)
  )

  df$p_adj <- suppressWarnings(
    as.numeric(df$p_adj)
  )

  df <- df[
    !is.na(df$p_adj) &
      df$p_adj < padj_cutoff &
      !is.na(df$r),
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    message(
      "plot_cooccurrence_network: no significant pairs at adjusted p < ",
      padj_cutoff,
      "."
    )

    return(invisible(NULL))
  }

  plots <- lapply(
    sort(unique(df$region)),
    function(region_name) {
      region_df <- df[
        df$region == region_name,
        ,
        drop = FALSE
      ]

      graph <- igraph::graph_from_data_frame(
        region_df[
          ,
          c(
            "celltype_1",
            "celltype_2",
            "r"
          )
        ],
        directed = FALSE
      )

      ggraph::ggraph(
        graph,
        layout = "fr"
      ) +
        ggraph::geom_edge_link(
          ggplot2::aes(
            edge_color = r
          ),
          edge_width = 0.6
        ) +
        ggraph::geom_node_point(
          size = 2,
          color = "grey30"
        ) +
        ggraph::geom_node_text(
          ggplot2::aes(
            label = name
          ),
          repel = TRUE,
          size = 2.5
        ) +
        ggraph::scale_edge_color_gradient2(
          low = "#2166AC",
          mid = "grey80",
          high = "#B2182B",
          midpoint = 0,
          name = "Correlation"
        ) +
        ggplot2::labs(
          title = region_name
        ) +
        ggplot2::theme_void()
    }
  )

  combined <- if (
    requireNamespace(
      "patchwork",
      quietly = TRUE
    )
  ) {
    patchwork::wrap_plots(
      plots,
      nrow = 1
    )
  } else {
    plots[[1]]
  }

  save_figure(
    combined,
    "cooccurrence_network",
    output_dir,
    width = 15,
    height = 5
  )

  invisible(combined)
}

# ============================================================
# Context-specific co-occurrence network
# ============================================================

#' Prepare context-specific co-occurrence edges
#'
#' Filters correlations by absolute effect size and optionally retains only
#' the strongest edges within each region x context panel.
#'
#' @keywords internal
.prepare_context_network_edges <- function(
    df,
    edge_threshold = 0.90,
    max_edges_per_panel = 120
) {
  .require_columns(
    df,
    c(
      "region",
      "context",
      "celltype_1",
      "celltype_2",
      "r",
      "n_scans"
    ),
    "Context-specific co-occurrence table"
  )

  df$region <- as.character(df$region)
  df$context <- as.character(df$context)
  df$celltype_1 <- as.character(df$celltype_1)
  df$celltype_2 <- as.character(df$celltype_2)
  df$r <- suppressWarnings(as.numeric(df$r))
  df$n_scans <- suppressWarnings(as.numeric(df$n_scans))

  df <- df[
    !is.na(df$region) &
      !is.na(df$context) &
      !is.na(df$celltype_1) &
      !is.na(df$celltype_2) &
      is.finite(df$r) &
      abs(df$r) >= edge_threshold,
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    return(df)
  }

  df$abs_r <- abs(df$r)

  df$correlation_sign <- ifelse(
    df$r >= 0,
    "Positive",
    "Negative"
  )

  df <- df |>
    dplyr::group_by(
      region,
      context
    ) |>
    dplyr::arrange(
      dplyr::desc(abs_r),
      .by_group = TRUE
    )

  if (
    !is.null(max_edges_per_panel) &&
      is.finite(max_edges_per_panel)
  ) {
    df <- df |>
      dplyr::slice_head(
        n = max_edges_per_panel
      )
  }

  df |>
    dplyr::ungroup()
}


#' Calculate a fixed union-graph layout for one region
#'
#' Creates a graph from all edges retained in any context for a region.
#' Every cell type is included as a vertex, including isolated vertices.
#' The resulting coordinates are reused for all context panels.
#'
#' @keywords internal
.make_region_union_layout <- function(
    region_edges,
    region_celltypes,
    layout_seed = 1234
) {
  region_celltypes <- sort(
    unique(
      as.character(region_celltypes)
    )
  )

  vertices <- data.frame(
    name = region_celltypes,
    lineage = .celltype_lineage(region_celltypes),
    stringsAsFactors = FALSE
  )

  if (nrow(region_edges) > 0) {
    union_edges <- region_edges |>
      dplyr::mutate(
        pair_1 = pmin(
          celltype_1,
          celltype_2
        ),
        pair_2 = pmax(
          celltype_1,
          celltype_2
        )
      ) |>
      dplyr::group_by(
        pair_1,
        pair_2
      ) |>
      dplyr::summarise(
        union_weight = max(
          abs_r,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) |>
      dplyr::transmute(
        from = pair_1,
        to = pair_2,
        union_weight = union_weight
      )
  } else {
    union_edges <- data.frame(
      from = character(),
      to = character(),
      union_weight = numeric(),
      stringsAsFactors = FALSE
    )
  }

  union_graph <- igraph::graph_from_data_frame(
    d = union_edges,
    directed = FALSE,
    vertices = vertices
  )

  union_degree <- igraph::degree(
    union_graph,
    mode = "all"
  )

  if (igraph::ecount(union_graph) == 0) {
    coordinates <- igraph::layout_in_circle(
      union_graph
    )
  } else {
    set.seed(layout_seed)

    coordinates <- igraph::layout_with_fr(
      union_graph,
      weights = igraph::E(union_graph)$union_weight,
      niter = 2000,
      grid = "nogrid"
    )
  }

  # Normalize each region to the same plotting range.
  normalize_coordinate <- function(x) {
    x_range <- range(
      x,
      finite = TRUE
    )

    if (
      !all(is.finite(x_range)) ||
        diff(x_range) == 0
    ) {
      return(rep(0, length(x)))
    }

    scales::rescale(
      x,
      to = c(-1, 1)
    )
  }

  layout_df <- data.frame(
    name = igraph::V(union_graph)$name,
    x = normalize_coordinate(coordinates[, 1]),
    y = normalize_coordinate(coordinates[, 2]),
    union_degree = as.numeric(
      union_degree[
        igraph::V(union_graph)$name
      ]
    ),
    stringsAsFactors = FALSE
  )

  layout_df$lineage <- .celltype_lineage(
    layout_df$name
  )

  layout_df
}


#' Plot context-specific cell-type co-occurrence networks
#'
#' Creates a 3 x 3 panel figure showing strong CLR correlations across
#' Arteries, Capillaries, and Parenchyma in Control amyloid-free,
#' AD amyloid-free, and AD amyloid-positive contexts.
#'
#' A union graph is constructed independently for each region using every
#' retained edge present in any context. Its layout is then reused across
#' all three context panels for that region, ensuring that changes between
#' panels represent changes in edges rather than changes in node placement.
#'
#' @param cooccurrence_csv Path to
#'   `celltype_cooccurrence_by_context.csv`.
#' @param output_dir Figure output directory.
#' @param edge_threshold Minimum absolute CLR correlation retained.
#' @param max_edges_per_panel Maximum number of strongest edges shown in
#'   each region x context panel. Set to `NULL` to disable the cap.
#' @param layout_seed Random seed used for the union layouts.
#' @param label_reference_context Context in which node labels are displayed.
#'   Since coordinates are fixed within regions, labels need only be shown
#'   in one row.
#'
#' @return The assembled patchwork figure, invisibly.
#' @export
plot_context_cooccurrence_network <- function(
    cooccurrence_csv,
    output_dir,
    edge_threshold = 0.90,
    max_edges_per_panel = 120,
    layout_seed = 1234,
    label_reference_context = "Control_AF"
) {
  required_packages <- c(
    "igraph",
    "ggraph",
    "patchwork",
    "scales"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing_packages) > 0) {
    message(
      "plot_context_cooccurrence_network: missing package(s): ",
      paste(missing_packages, collapse = ", "),
      ". Figure skipped."
    )

    return(invisible(NULL))
  }

  df <- utils::read.csv(
    cooccurrence_csv,
    stringsAsFactors = FALSE
  )

  .require_columns(
    df,
    c(
      "region",
      "context",
      "celltype_1",
      "celltype_2",
      "r",
      "n_scans"
    ),
    "Context-specific co-occurrence table"
  )

  region_order <- c(
    "Arteries",
    "Capillaries",
    "Parenchyma"
  )

  context_order <- c(
    "Control_AF",
    "AD_AF",
    "AD_Abeta"
  )

  context_labels <- c(
  Control_AF = "Control, amyloid-free",
  AD_AF = "AD, amyloid-free",
  AD_Abeta = "AD, amyloid"
)

  # Keep expected strata while allowing the function to report missing ones.
  df <- df[
    df$region %in% region_order &
      df$context %in% context_order,
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    message(
      "plot_context_cooccurrence_network: no recognized ",
      "region/context rows were found."
    )

    return(invisible(NULL))
  }

  # Preserve every cell type found in the full table, not only those passing
  # the edge filter. This keeps isolated nodes visible.
  region_celltypes <- lapply(
    region_order,
    function(region_name) {
      region_df <- df[
        df$region == region_name,
        ,
        drop = FALSE
      ]

      sort(
        unique(
          c(
            region_df$celltype_1,
            region_df$celltype_2
          )
        )
      )
    }
  )

  names(region_celltypes) <- region_order

  edge_df <- .prepare_context_network_edges(
    df = df,
    edge_threshold = edge_threshold,
    max_edges_per_panel = max_edges_per_panel
  )

  lineage_palette <- c(
    Excitatory = "#E69F00",
    Inhibitory = "#D55E00",
    Astrocyte = "#009E73",
    Oligodendrocyte = "#0072B2",
    OPC = "#56B4E9",
    Microglia = "#CC79A7",
    Endothelial = "#8C564B",
    Pericyte = "#9467BD",
    SMC = "#7F7F7F",
    VLMC = "#E377C2",
    Fibroblast = "#BCBD22",
    Other = "#BDBDBD"
  )

  correlation_palette <- c(
    Positive = "#B2182B",
    Negative = "#2166AC"
  )

  region_layouts <- list()

  for (region_index in seq_along(region_order)) {
    region_name <- region_order[region_index]

    region_edges <- edge_df[
      edge_df$region == region_name,
      ,
      drop = FALSE
    ]

    region_layouts[[region_name]] <-
      .make_region_union_layout(
        region_edges = region_edges,
        region_celltypes = region_celltypes[[region_name]],
        layout_seed = layout_seed + region_index
      )
  }
  # ------------------------------------------------------------
  # Shared node-size scale across every panel
  # ------------------------------------------------------------

  all_union_degrees <- unlist(
    lapply(
      region_layouts,
      function(layout_df) {
        layout_df$union_degree
      }
    ),
    use.names = FALSE
  )

  all_union_degrees <- all_union_degrees[
    is.finite(all_union_degrees)
  ]

  union_degree_max <- max(
    all_union_degrees,
    na.rm = TRUE
  )

  union_degree_breaks <- pretty(
    c(0, union_degree_max),
    n = 4
  )

  union_degree_breaks <- union_degree_breaks[
    union_degree_breaks > 0 &
      union_degree_breaks <= union_degree_max
  ]

  plot_list <- list()

  for (context_name in context_order) {
    for (region_name in region_order) {
      panel_edges <- edge_df[
        edge_df$region == region_name &
          edge_df$context == context_name,
        ,
        drop = FALSE
      ]

      layout_df <- region_layouts[[region_name]]

      vertex_df <- layout_df[
        ,
        c(
          "name",
          "lineage",
          "union_degree"
        ),
        drop = FALSE
      ]

      if (nrow(panel_edges) > 0) {
        graph_edges <- panel_edges |>
          dplyr::transmute(
            from = celltype_1,
            to = celltype_2,
            r = r,
            abs_r = abs_r,
            correlation_sign = correlation_sign
          )
      } else {
        graph_edges <- data.frame(
          from = character(),
          to = character(),
          r = numeric(),
          abs_r = numeric(),
          correlation_sign = character(),
          stringsAsFactors = FALSE
        )
      }

      panel_graph <- igraph::graph_from_data_frame(
        d = graph_edges,
        directed = FALSE,
        vertices = vertex_df
      )

      graph_vertex_order <- igraph::V(panel_graph)$name

      panel_layout <- layout_df[
        match(
          graph_vertex_order,
          layout_df$name
        ),
        ,
        drop = FALSE
      ]

      n_scans_values <- unique(
        df$n_scans[
          df$region == region_name &
            df$context == context_name
        ]
      )

      n_scans_values <- n_scans_values[
        is.finite(n_scans_values)
      ]

      n_scans <- if (length(n_scans_values) > 0) {
        n_scans_values[1]
      } else {
        NA_real_
      }

      panel_subtitle <- paste0(
        "n = ",
        ifelse(
          is.na(n_scans),
          "NA",
          n_scans
        ),
        " scans; ",
        nrow(panel_edges),
        " edges"
      )

      panel_plot <- ggraph::ggraph(
        panel_graph,
        layout = "manual",
        x = panel_layout$x,
        y = panel_layout$y
      )

      if (igraph::ecount(panel_graph) > 0) {
        panel_plot <- panel_plot +
          ggraph::geom_edge_link(
            ggplot2::aes(
              edge_colour = correlation_sign,
              edge_width = abs_r,
              edge_alpha = abs_r
            ),
            lineend = "round",
            show.legend = TRUE
          )
      }

      panel_plot <- panel_plot +
        ggraph::geom_node_point(
          ggplot2::aes(
            fill = lineage,
            size = union_degree
          ),
          shape = 21,
          colour = "grey20",
          stroke = 0.3
        )

      # Label nodes only in one reference context. Because coordinates are
      # fixed within each region, the labels apply to the panels below it.
      if (identical(
        context_name,
        label_reference_context
      )) {
        panel_plot <- panel_plot +
          ggraph::geom_node_text(
            ggplot2::aes(
              label = name
            ),
            repel = TRUE,
            size = 2.1,
            family = "sans",
            max.overlaps = Inf,
            box.padding = 0.25,
            point.padding = 0.15
          )
      }

      panel_plot <- panel_plot +
        ggraph::scale_edge_colour_manual(
          name = "CLR correlation",
          values = correlation_palette,
          drop = FALSE
        ) +
        ggraph::scale_edge_width_continuous(
          name = "|r|",
          range = c(0.25, 1.6),
          limits = c(
            edge_threshold,
            1
          )
        ) +
        ggraph::scale_edge_alpha_continuous(
          range = c(0.25, 0.85),
          limits = c(
            edge_threshold,
            1
          ),
          guide = "none"
        ) +
        ggplot2::scale_fill_manual(
          name = "Cell lineage",
          values = lineage_palette,
          drop = FALSE
        ) +
        ggplot2::scale_size_continuous(
          name = "Union degree",
          limits = c(
            0,
            union_degree_max
          ),
          breaks = union_degree_breaks,
          range = c(
            2.2,
            7
          ),
          guide = ggplot2::guide_legend(
            order = 1,
            override.aes = list(
              alpha = 1
            )
          )
        ) +
		ggplot2::guides(
 		 size = ggplot2::guide_legend(
		    order = 1
		  ),
		  fill = ggplot2::guide_legend(
		    order = 2
		  ),
		  edge_colour = ggplot2::guide_legend(
		    order = 3
		  ),
		  edge_width = ggplot2::guide_legend(
		    order = 4
		  )
		)+
        ggplot2::coord_equal(
          xlim = c(-1.35, 1.35),
          ylim = c(-1.35, 1.35),
          clip = "off"
        ) +
        ggplot2::labs(
          title = paste0(
            region_name,
            "\n",
            unname(context_labels[context_name])
          ),
          subtitle = panel_subtitle
        ) +
        ggplot2::theme_void(
          base_family = "sans"
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 9,
            face = "bold",
            hjust = 0.5,
            lineheight = 1.05,
            margin = ggplot2::margin(
              b = 2
            )
          ),
          plot.subtitle = ggplot2::element_text(
            size = 7.5,
            colour = "grey35",
            hjust = 0.5
          ),
          plot.margin = ggplot2::margin(
            8,
            8,
            8,
            8
          ),
          legend.title = ggplot2::element_text(
            size = 9,
            face = "bold"
          ),
          legend.text = ggplot2::element_text(
            size = 8
          )
        )

      plot_key <- paste(
        context_name,
        region_name,
        sep = "__"
      )

      plot_list[[plot_key]] <- panel_plot
    }
  }

  row_plots <- list()

  for (context_name in context_order) {
    context_panels <- lapply(
      region_order,
      function(region_name) {
        plot_list[[
          paste(
            context_name,
            region_name,
            sep = "__"
          )
        ]]
      }
    )

    row_plot <- patchwork::wrap_plots(
      context_panels,
      nrow = 1
    )

    row_plots[[context_name]] <- row_plot
  }

  assembled_plot <- patchwork::wrap_plots(
  row_plots,
  ncol = 1,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Cell-type co-occurrence networks by spatial context",
    subtitle = paste0(
      "Strong CLR correlations across independent scans; ",
      "edges retained at |r| >= ",
      edge_threshold,
      if (
        is.null(max_edges_per_panel)
      ) {
        ""
      } else {
        paste0(
          "; up to ",
          max_edges_per_panel,
          " strongest edges per panel"
        )
      }
    ),
    caption = paste0(
      "Node coordinates are derived from the union of retained edges ",
      "across all three contexts within each region and are held fixed ",
      "across rows. Edge color indicates correlation direction and edge ",
      "width indicates |r|. Networks are descriptive because each panel ",
      "contains only 3-4 independent scans."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 15,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10
      ),
      plot.caption = ggplot2::element_text(
        size = 8,
        colour = "grey35",
        hjust = 0
      ),
      legend.position = "right"
    )
  )


  save_figure(
    assembled_plot,
    "cooccurrence_network_by_context",
    output_dir,
    width = 16,
    height = 18
  )

  invisible(assembled_plot)
}
