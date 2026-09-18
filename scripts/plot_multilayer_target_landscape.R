#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

# ============================================================
# FIGURE P
# Continuous multi-layer target evidence landscape
#
# Each point = one evaluated AGORA target
#
# X = AGORA nomination frequency
# Y = number of spatial FDR-supported effects (0-9)
# Size = fine-cell FDR associations
# Color = lineage attribution evidence
# ============================================================

# ------------------------------------------------------------
# Input
# ------------------------------------------------------------

single_file <- paste0(
  "results/gene_panel_validation/",
  "AGORA_Single_Nomination/",
  "integrated_gene_panel_summary.csv"
)

high_file <- paste0(
  "results/gene_panel_validation/",
  "AGORA_High_Nomination/",
  "integrated_gene_panel_summary.csv"
)

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------

single <- read_csv(
  single_file,
  show_col_types = FALSE
) %>%
  mutate(
    total_nominations = 1
  )

high <- read_csv(
  high_file,
  show_col_types = FALSE
) %>%
  mutate(
    total_nominations = as.numeric(
      sub(
        "AGORA nominations: ",
        "",
        category
      )
    )
  )

x <- bind_rows(
  single,
  high
) %>%
  filter(
    !is.na(n_spatial_fdr_significant)
  )

cat("\n========================================\n")
cat("TARGET LANDSCAPE\n")
cat("========================================\n")

cat(
  "Targets plotted:",
  nrow(x),
  "\n"
)

print(
  table(
    x$total_nominations
  )
)

# ============================================================
# Fine-cell evidence
# ============================================================

x <- x %>%
  mutate(
    n_fine_fdr_total =
      coalesce(
        n_fine_fdr_significant_Amyloid,
        0
      ) +
      coalesce(
        n_fine_fdr_significant_AmyloidFree,
        0
      ),

    fine_cell_bin = case_when(
      n_fine_fdr_total == 0 ~ "0",
      n_fine_fdr_total <= 5 ~ "1-5",
      n_fine_fdr_total <= 15 ~ "6-15",
      TRUE ~ ">15"
    ),

    fine_cell_bin = factor(
      fine_cell_bin,
      levels = c(
        "0",
        "1-5",
        "6-15",
        ">15"
      )
    )
  )

# ============================================================
# Lineage evidence
# ============================================================

x <- x %>%
  mutate(
    lineage_class = case_when(
      lineage_evidence ==
        "Exploratory nominal" ~
        "Exploratory lineage association",

      TRUE ~
        "No lineage association"
    ),

    lineage_class = factor(
      lineage_class,
      levels = c(
        "No lineage association",
        "Exploratory lineage association"
      )
    )
  )

# ============================================================
# Continuous jitter
# ============================================================

set.seed(1234)

x <- x %>%
  mutate(
    x_jitter_width = case_when(
      total_nominations == 1 ~ 0.38,
      total_nominations == 3 ~ 0.22,
      total_nominations == 4 ~ 0.18,
      total_nominations == 5 ~ 0.14,
      TRUE ~ 0.15
    ),

    plot_x =
      total_nominations +
      runif(
        n(),
        min = -x_jitter_width,
        max = x_jitter_width
      ),

    plot_y =
      n_spatial_fdr_significant +
      runif(
        n(),
        min = -0.08,
        max = 0.08
      )
  )

# ============================================================
# Genes to label
# ============================================================

label_genes <- c(
  "CLU",
  "PLEC",
  "BIN1",
  "RUFY3",
  "FAM107A",
  "TAL1",
  "SYNPO",
  "STX1B"
)

# ============================================================
# Plot
# ============================================================

p <- ggplot(
  x,
  aes(
    x = plot_x,
    y = plot_y
  )
) +

  # ----------------------------------------------------------
  # High-spatial-support reference line
  # ----------------------------------------------------------

  geom_hline(
    yintercept = 3,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey45"
  ) +

  annotate(
    "text",
    x = 5.28,
    y = 3.13,
    label = "Broad spatial support",
    hjust = 1,
    size = 3.4,
    color = "grey35"
  ) +

  # ----------------------------------------------------------
  # All targets
  # ----------------------------------------------------------

  geom_point(
    aes(
      size = fine_cell_bin,
      color = lineage_class
    ),
    alpha = 0.62
  ) +

  # ----------------------------------------------------------
  # Highlight labeled targets
  # ----------------------------------------------------------

  geom_point(
    data = x %>%
      filter(
        gene %in% label_genes
      ),

    aes(
      x = plot_x,
      y = plot_y,
      size = fine_cell_bin
    ),

    shape = 21,
    fill = NA,
    color = "black",
    stroke = 0.9,
    show.legend = FALSE
  ) +

  # ----------------------------------------------------------
  # Labels
  # ----------------------------------------------------------

  ggrepel::geom_text_repel(
    data = x %>%
      filter(
        gene %in% label_genes
      ),

    aes(
      x = plot_x,
      y = plot_y,
      label = gene
    ),

    size = 4,
    fontface = "bold",

    box.padding = 0.45,
    point.padding = 0.35,

    min.segment.length = 0,

    segment.color = "grey40",
    segment.size = 0.4,

    max.overlaps = Inf,

    show.legend = FALSE
  ) +

  # ==========================================================
  # X axis
  # ==========================================================

  scale_x_continuous(
  breaks = 1:5,
  labels = 1:5,
  limits = c(0.35, 5.65),
  expand = c(0, 0)
) +

  # ==========================================================
  # Y axis
  # ==========================================================

  scale_y_continuous(
    breaks = 0:6,

    limits = c(
      -0.35,
      6.6
    )
  ) +

  # ==========================================================
  # Size
  # ==========================================================

  scale_size_manual(
    values = c(
      "0" = 1.7,
      "1-5" = 2.8,
      "6-15" = 4.0,
      ">15" = 5.3
    ),

    name =
      "Fine-cell FDR\nassociations"
  ) +

  # ==========================================================
  # Color
  # ==========================================================

  scale_color_manual(
    values = c(
      "No lineage association" =
        "#F28E82",

      "Exploratory lineage association" =
        "#26BFC7"
    ),

    name =
      "Lineage attribution"
  ) +

  # ==========================================================
  # Labels
  # ==========================================================

  labs(
    title = NULL,

    subtitle =
      "Each point represents one evaluated AGORA therapeutic target",

    x =
      "AGORA nomination frequency",

    y =
      "FDR-supported spatial effects\n(out of 9)"
  ) +

  # ==========================================================
  # Theme
  # ==========================================================

  theme_classic(
    base_size = 13
  ) +

  theme(
    axis.title =
      element_text(
        face = "bold",
        size = 13
      ),

    axis.text =
      element_text(
        size = 11
      ),

    axis.text.x =
      element_text(
        face = "bold"
      ),

    legend.position =
      "right",

    legend.title =
      element_text(
        face = "bold",
        size = 10
      ),

    legend.text =
      element_text(
        size = 9
      ),

    plot.subtitle =
      element_text(
        size = 11,
        margin = margin(
          b = 12
        )
      ),

    plot.margin =
      margin(
        t = 10,
        r = 15,
        b = 10,
        l = 10
      )
  )

print(p)

# ============================================================
# Save
# ============================================================

output_dir <- paste0(
  "results/figures/gene_panel_validation/",
  "AGORA_Nomination_Comparison"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pdf_file <- file.path(
  output_dir,
  "Figure_P_multilayer_target_continuous.pdf"
)

png_file <- file.path(
  output_dir,
  "Figure_P_multilayer_target_continuous.png"
)

ggsave(
  pdf_file,
  p,
  width = 12.5,
  height = 4.8,
  device = cairo_pdf
)

ggsave(
  png_file,
  p,
  width = 12.5,
  height = 4.8,
  dpi = 350
)

write_csv(
  x,
  file.path(
    output_dir,
    "Figure_P_multilayer_target_continuous_data.csv"
  )
)

cat(
  "\n========================================\n",
  "WROTE CONTINUOUS FIGURE P\n",
  "========================================\n",
  pdf_file,
  "\n",
  png_file,
  "\n",
  sep = ""
)
