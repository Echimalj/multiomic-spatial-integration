#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ============================================================
# FIGURE A
# AGORA nomination frequency vs human spatial support
#
# Question:
# Does repeated AGORA nomination predict whether a target
# shows at least one FDR-supported spatial perturbation?
# ============================================================

# ------------------------------------------------------------
# Input
# ------------------------------------------------------------

input_file <- paste0(
  "results/gene_panel_validation/",
  "AGORA_Nomination_Comparison/",
  "gene_level_spatial_evidence.csv"
)

x <- read_csv(
  input_file,
  show_col_types = FALSE
)

cat("\nNomination groups:\n")
print(
  table(
    x$nomination_group,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# Define binary spatial-support variable
# ------------------------------------------------------------

x <- x %>%
  mutate(
    spatial_fdr_support = n_fdr > 0,

    # Use simple ASCII labels to avoid font/rendering issues
    group_label = recode(
      nomination_group,
      "1 nomination" = "1 nomination",
      "3–5 nominations" = "3-5 nominations"
    )
  )

# ------------------------------------------------------------
# Summarize by nomination group
# ------------------------------------------------------------

support_summary <- x %>%
  group_by(group_label) %>%
  summarise(
    n_genes = n(),

    n_supported = sum(
      spatial_fdr_support,
      na.rm = TRUE
    ),

    n_not_supported = sum(
      !spatial_fdr_support,
      na.rm = TRUE
    ),

    proportion_supported =
      n_supported / n_genes,

    percent_supported =
      100 * proportion_supported,

    .groups = "drop"
  )

cat("\n========================================\n")
cat("SPATIAL SUPPORT SUMMARY\n")
cat("========================================\n")

print(
  support_summary,
  width = Inf
)

# ------------------------------------------------------------
# Fisher's exact test
#
# Important:
# Test is performed on the original groups, not plotting labels.
# ------------------------------------------------------------

contingency_table <- table(
  x$nomination_group,
  x$spatial_fdr_support
)

cat("\n========================================\n")
cat("FISHER'S EXACT TEST\n")
cat("========================================\n")

print(contingency_table)

fisher_result <- fisher.test(
  contingency_table
)

print(fisher_result)

fisher_p <- fisher_result$p.value

# ------------------------------------------------------------
# Set plotting order
# ------------------------------------------------------------

support_summary$group_label <- factor(
  support_summary$group_label,
  levels = c(
    "1 nomination",
    "3-5 nominations"
  )
)

# ------------------------------------------------------------
# Bar labels
# ------------------------------------------------------------

support_summary <- support_summary %>%
  mutate(
    bar_label = paste0(
      n_supported,
      "/",
      n_genes,
      "\n",
      sprintf(
        "%.1f%%",
        percent_supported
      )
    )
  )

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot(
  support_summary,
  aes(
    x = group_label,
    y = proportion_supported,
    fill = group_label
  )
) +

  # Bars
  geom_col(
    width = 0.52
  ) +

  # Different colors for each nomination group
  scale_fill_manual(
    values = c(
      "1 nomination" = "#4C78A8",
      "3-5 nominations" = "#E07B39"
    ),
    guide = "none"
  ) +

  # Counts and percentages
  geom_text(
    aes(
      label = bar_label
    ),
    vjust = -0.50,
    size = 4.7,
    fontface = "bold",
    lineheight = 1.05
  ) +

  # ----------------------------------------------------------
  # Statistical comparison bracket
  # ----------------------------------------------------------

  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = 0.505,
    yend = 0.505,
    linewidth = 0.6
  ) +

  annotate(
    "segment",
    x = 1,
    xend = 1,
    y = 0.485,
    yend = 0.505,
    linewidth = 0.6
  ) +

  annotate(
    "segment",
    x = 2,
    xend = 2,
    y = 0.485,
    yend = 0.505,
    linewidth = 0.6
  ) +

  annotate(
    "text",
    x = 1.5,
    y = 0.535,
    label = paste0(
      "Fisher's exact p = ",
      sprintf(
        "%.2f",
        fisher_p
      )
    ),
    size = 4
  ) +

  # ----------------------------------------------------------
  # Y axis
  # ----------------------------------------------------------

  scale_y_continuous(
    labels = function(z) {
      paste0(
        round(
          z * 100
        ),
        "%"
      )
    },

    breaks = seq(
      0,
      0.6,
      by = 0.1
    ),

    limits = c(
      0,
      0.60
    ),

    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  # ----------------------------------------------------------
  # Labels
  # ----------------------------------------------------------

  labs(
    title =
      "Spatial support provides evidence beyond nomination frequency",

    subtitle =
      "Targets with at least one FDR-supported disease, amyloid, or maximum-pathology effect",

    x = NULL,

    y =
      "Targets with spatial\nFDR evidence"
  ) +

  # ----------------------------------------------------------
  # Theme
  # ----------------------------------------------------------

  theme_classic(
    base_size = 13
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0,
      margin = margin(
        b = 4
      )
    ),

    plot.subtitle = element_text(
      size = 11,
      hjust = 0,
      margin = margin(
        b = 14
      )
    ),

    axis.title.y = element_text(
      face = "bold",
      size = 12,
      margin = margin(
        r = 10
      )
    ),

    axis.text.x = element_text(
      face = "bold",
      size = 12,
      margin = margin(
        t = 8
      )
    ),

    axis.text.y = element_text(
      size = 11
    ),

    axis.ticks.x = element_blank(),

    legend.position = "none",

    plot.margin = margin(
      t = 12,
      r = 20,
      b = 12,
      l = 12
    )
  )

print(p)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

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
  "Figure_A_nomination_vs_spatial_support.pdf"
)

png_file <- file.path(
  output_dir,
  "Figure_A_nomination_vs_spatial_support.png"
)

ggsave(
  filename = pdf_file,
  plot = p,
  width = 7.5,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = png_file,
  plot = p,
  width = 7.5,
  height = 6,
  dpi = 350
)

# ------------------------------------------------------------
# Save numerical summary
# ------------------------------------------------------------

write_csv(
  support_summary,
  file.path(
    output_dir,
    "Figure_A_nomination_vs_spatial_support_summary.csv"
  )
)

cat(
  "\n========================================\n",
  "WROTE FIGURE A\n",
  "========================================\n",
  pdf_file,
  "\n",
  png_file,
  "\n",
  sep = ""
)
