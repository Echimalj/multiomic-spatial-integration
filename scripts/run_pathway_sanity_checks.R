# ============================================================
# Pathway sanity checks
# ============================================================

library(dplyr)

source("R/enrichment_utils.R")

markers <- read.table("results/AD_CAA_cluster_markers.txt", header = TRUE)

genes <- extract_top_genes(markers)

genes <- format_gene_list(genes)

write.table(
  genes,
  file = "results/pathway_gene_list.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
