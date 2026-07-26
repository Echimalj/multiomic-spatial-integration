###############################################################################
# test_absolute_abundance_utils.R
#
# Regression tests for R/absolute_abundance_utils.R
###############################################################################

library(glmmTMB)
library(emmeans)

source("R/contrast_utils.R")
source("R/absolute_abundance_utils.R")

cat("=========================================================\n")
cat("Generating synthetic dataset\n")
cat("=========================================================\n")

n_scans <- 12
rois_per_scan <- 8

test_df <- expand.grid(
  Scan_ID = paste0("Scan", seq_len(n_scans)),
  ROI = seq_len(rois_per_scan),
  celltype = "Astrocytes",
  KEEP.OUT.ATTRS = FALSE
)

test_df$disease_status <- rep(
  rep(c("Control", "AD-CAA"), each = rois_per_scan),
  length.out = nrow(test_df)
)

test_df$pathology <- ifelse(
  test_df$disease_status == "Control",
  "AmyloidFree",
  rep(c("AmyloidFree", "Amyloid"), length.out = nrow(test_df))
)

test_df$region <- rep(
  c("Frontal", "Temporal"),
  length.out = nrow(test_df)
)

test_df$size_factor_mor <- exp(
  rnorm(nrow(test_df), mean = 0, sd = 0.25)
)

scan_effect <- rnorm(n_scans, mean = 0, sd = 0.15)
names(scan_effect) <- paste0("Scan", seq_len(n_scans))

true_log_fc <- log(2)

test_df$abs_abundance <- exp(
  log(test_df$size_factor_mor) +
    true_log_fc *
    (test_df$disease_status == "AD-CAA") +
    scan_effect[test_df$Scan_ID] +
    rnorm(nrow(test_df), sd = 0.08)
)

###############################################################################
# Helper
###############################################################################

fit_model <- function(df){

  prepared <- prepare_absolute_abundance_data(
    df = df,
    abundance_col = "abs_abundance",
    offset_col = "size_factor_mor"
  )

  glmmTMB(
    log_abs_abundance ~ disease_status +
      offset(log_offset) +
      (1 | Scan_ID),
    data = prepared,
    family = gaussian()
  )

}

get_contrast <- function(fit){

  res <- contrast(
    emmeans(fit, ~ disease_status),
    method = list(
      AD_vs_Control = c(-1,1)
    )
  )

  as.data.frame(res)$estimate[[1]]

}

###############################################################################
# TEST 1
###############################################################################

cat("\n=========================================================\n")
cat("TEST 1 : Known effect recovery\n")
cat("=========================================================\n")

fit <- fit_model(test_df)

estimate <- get_contrast(fit)

cat("Expected log effect :", true_log_fc, "\n")
cat("Estimated log effect:", estimate, "\n")
cat("Expected fold change:", exp(true_log_fc), "\n")
cat("Estimated fold change:", exp(estimate), "\n")

stopifnot(abs(estimate - true_log_fc) < 0.15)

cat("PASS\n")


###############################################################################
# TEST 2
###############################################################################

cat("\n=========================================================\n")
cat("TEST 2 : Offset shift invariance\n")
cat("=========================================================\n")

prepared_original <- prepare_absolute_abundance_data(
  test_df,
  abundance_col="abs_abundance",
  offset_col="size_factor_mor"
)

prepared_shifted <- prepared_original
prepared_shifted$log_offset <- prepared_shifted$log_offset + 5

fit_original <- glmmTMB(
  log_abs_abundance ~ disease_status +
    offset(log_offset) +
    (1|Scan_ID),
  data=prepared_original,
  family=gaussian()
)

fit_shifted <- glmmTMB(
  log_abs_abundance ~ disease_status +
    offset(log_offset) +
    (1|Scan_ID),
  data=prepared_shifted,
  family=gaussian()
)

contrast_original <- get_contrast(fit_original)
contrast_shifted <- get_contrast(fit_shifted)

difference <- abs(contrast_original - contrast_shifted)

cat("Original contrast :", contrast_original, "\n")
cat("Shifted contrast  :", contrast_shifted, "\n")
cat("Absolute difference:", difference, "\n")

stopifnot(difference < 1e-6)

cat("PASS\n")


###############################################################################
# TEST 3
###############################################################################

cat("\n=========================================================\n")
cat("TEST 3 : Global abundance scaling invariance\n")
cat("=========================================================\n")

scaled_df <- test_df
scaled_df$abs_abundance <- scaled_df$abs_abundance * 100

fit_scaled <- fit_model(scaled_df)

contrast_scaled <- get_contrast(fit_scaled)

difference_scaled <- abs(contrast_original - contrast_scaled)

cat("Original contrast :", contrast_original, "\n")
cat("Scaled contrast   :", contrast_scaled, "\n")
cat("Absolute difference:", difference_scaled, "\n")

stopifnot(difference_scaled < 1e-6)

cat("PASS\n")


###############################################################################
# TEST 4
###############################################################################

cat("\n=========================================================\n")
cat("TEST 4 : Negative abundance handling\n")
cat("=========================================================\n")

bad_df <- test_df
bad_df$abs_abundance[1] <- -1

negative_result <- try(

  prepare_absolute_abundance_data(
    bad_df,
    abundance_col="abs_abundance",
    offset_col="size_factor_mor"
  ),

  silent=TRUE

)

if(inherits(negative_result,"try-error")){

  cat("PASS : negative abundance rejected.\n")

}else{

  cat("WARNING : negative abundance accepted.\n")
  cat("This should probably be rejected explicitly.\n")

}


###############################################################################
# TEST 5
###############################################################################

cat("\n=========================================================\n")
cat("TEST 5 : Zero offset handling\n")
cat("=========================================================\n")

bad_df2 <- test_df
bad_df2$size_factor_mor[1] <- 0

zero_result <- try(

  prepare_absolute_abundance_data(
    bad_df2,
    abundance_col="abs_abundance",
    offset_col="size_factor_mor"
  ),

  silent=TRUE

)

if(inherits(zero_result,"try-error")){

  cat("PASS : zero offset rejected.\n")

}else{

  cat("WARNING : zero offset accepted.\n")
  cat("This should probably be rejected explicitly.\n")

}
###############################################################################
# TEST 6
###############################################################################

cat("\n=========================================================\n")
cat("TEST 6 : Condition-specific scaling\n")
cat("=========================================================\n")

scaled_condition <- test_df

## Artificially double ALL AD abundances
scaled_condition$abs_abundance[
  scaled_condition$disease_status == "AD-CAA"
] <-
scaled_condition$abs_abundance[
  scaled_condition$disease_status == "AD-CAA"
] * 2

fit_condition_scaled <- fit_model(scaled_condition)

condition_scaled_contrast <- get_contrast(fit_condition_scaled)

expected_contrast <- contrast_original + log(2)

cat("Original contrast          :", contrast_original, "\n")
cat("Condition-scaled contrast  :", condition_scaled_contrast, "\n")
cat("Expected contrast          :", expected_contrast, "\n")

difference <- abs(
  condition_scaled_contrast - expected_contrast
)

cat("Absolute difference        :", difference, "\n")

stopifnot(difference < 1e-6)

cat("PASS\n")

###############################################################################
# SUMMARY
###############################################################################

cat("\n=========================================================\n")
cat("All regression tests completed.\n")
cat("=========================================================\n")
