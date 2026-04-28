# ============================================================
# Spatial statistical modeling (beta mixed-effects)
# ============================================================

library(glmmTMB)
library(emmeans)
library(dplyr)

source("R/contrast_utils.R")

df <- read.csv("results/cell_proportions/spatial_celltype_proportions_for_R.csv")

df <- prepare_contrasts(df)

# Example: model per celltype
run_model <- function(df, celltype) {

  df_sub <- df[df$celltype == celltype, ]

  model <- glmmTMB(
    rel_abundance ~ disease_status + (1 | ROI),
    data = df_sub,
    family = beta_family(link = "logit")
  )

  res <- emmeans(model, pairwise ~ disease_status)

  return(res)
}

results <- lapply(unique(df$celltype), function(ct) {
  run_model(df, ct)
})

names(results) <- unique(df$celltype)
