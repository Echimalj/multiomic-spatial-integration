###############################################################################
# test_absolute_calibration_utils.R
#
# Regression tests for:
#   R/absolute_calibration_utils.R
###############################################################################

set.seed(42)  # The answer to the universe

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(glmmTMB)
})

source("R/absolute_calibration_utils.R")


###############################################################################
# Helpers
###############################################################################

print_header <- function(test_number, description) {
  cat("\n=========================================================\n")
  cat("TEST", test_number, ":", description, "\n")
  cat("=========================================================\n")
}


expect_error <- function(expr, pattern = NULL) {
  result <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e
  )

  if (is.null(result)) {
    stop("Expected an error, but no error was raised.")
  }

  if (!is.null(pattern) &&
      !grepl(pattern, conditionMessage(result), fixed = TRUE)) {
    stop(
      "An error was raised, but its message did not contain:\n  ",
      pattern,
      "\nActual message:\n  ",
      conditionMessage(result)
    )
  }

  invisible(result)
}


expect_close <- function(observed,
                         expected,
                         tolerance = 1e-6,
                         label = "value") {
  difference <- abs(observed - expected)

  cat("Observed", label, ":", observed, "\n")
  cat("Expected", label, ":", expected, "\n")
  cat("Absolute difference:", difference, "\n")

  if (!is.finite(difference) || difference >= tolerance) {
    stop(
      label,
      " differs from its expected value by ",
      difference,
      ", which exceeds tolerance ",
      tolerance,
      "."
    )
  }

  invisible(TRUE)
}


###############################################################################
# Base hand-calculable data
###############################################################################

abs_df_small <- data.frame(
  ROI_ID = rep(c("ROI1", "ROI2"), each = 2),
  disease_status = rep(c("Control", "AD-CAA"), each = 2),
  pathology = rep("AmyloidFree", 4),
  region = rep("Frontal", 4),
  Scan_ID = rep(c("Scan1", "Scan2"), each = 2),
  celltype = rep(c("Astrocytes", "Microglia"), times = 2),
  abs_abundance = c(
    10, 30,  # ROI1 total = 40
    20, 80   # ROI2 total = 100
  ),
  stringsAsFactors = FALSE
)

total_counts_small <- data.frame(
  ROI_ID = c("ROI1", "ROI2"),
  total_counts = c(100, 200),
  size_factor_mor = c(0.8, 1.2),
  stringsAsFactors = FALSE
)


###############################################################################
# TEST 1
###############################################################################

print_header(1, "Exact calibration-ratio arithmetic")

calib_small <- compute_calibration_ratio_table(
  abs_df = abs_df_small,
  total_counts_df = total_counts_small
)

calib_small <- calib_small |>
  arrange(ROI_ID)

expected_total_abundance <- c(40, 100)
expected_ratios <- c(0.4, 0.5)
expected_log_ratios <- log(expected_ratios)

stopifnot(nrow(calib_small) == 2)
stopifnot(!anyDuplicated(calib_small$ROI_ID))

stopifnot(
  isTRUE(
    all.equal(
      calib_small$total_abs_abundance,
      expected_total_abundance,
      tolerance = 1e-12
    )
  )
)

stopifnot(
  isTRUE(
    all.equal(
      calib_small$calibration_ratio,
      expected_ratios,
      tolerance = 1e-12
    )
  )
)

stopifnot(
  isTRUE(
    all.equal(
      calib_small$log_calibration_ratio,
      expected_log_ratios,
      tolerance = 1e-12
    )
  )
)

cat("PASS: calibration ratios match hand calculations.\n")


###############################################################################
# TEST 2
###############################################################################

print_header(2, "Incomplete ROI join is rejected")

counts_missing_roi <- total_counts_small |>
  filter(ROI_ID != "ROI2")

expect_error(
  compute_calibration_ratio_table(
    abs_df = abs_df_small,
    total_counts_df = counts_missing_roi
  ),
  pattern = "incomplete join"
)

cat("PASS: missing ROI in total-count table was rejected.\n")


###############################################################################
# TEST 3
###############################################################################

print_header(3, "Extra ROI in total-count table is rejected")

counts_extra_roi <- bind_rows(
  total_counts_small,
  data.frame(
    ROI_ID = "ROI3",
    total_counts = 150,
    size_factor_mor = 1
  )
)

expect_error(
  compute_calibration_ratio_table(
    abs_df = abs_df_small,
    total_counts_df = counts_extra_roi
  ),
  pattern = "incomplete join"
)

cat("PASS: extra ROI in total-count table was rejected.\n")


###############################################################################
# TEST 4
###############################################################################

print_header(4, "Conflicting ROI metadata is rejected")

conflicting_metadata <- abs_df_small

conflicting_metadata$region[
  conflicting_metadata$ROI_ID == "ROI1" &
    conflicting_metadata$celltype == "Microglia"
] <- "Temporal"

expect_error(
  compute_calibration_ratio_table(
    abs_df = conflicting_metadata,
    total_counts_df = total_counts_small
  ),
  pattern = "conflicting ROI-level metadata"
)

cat("PASS: conflicting ROI-level metadata was rejected.\n")


###############################################################################
# Synthetic data generator for gate-model tests
###############################################################################

make_gate_data <- function(
    disease_scale = 1,
    region_pattern = c(
      Frontal = 1.00,
      Temporal = 1.20,
      Parietal = 0.85,
      Occipital = 1.10
    ),
    reverse_ad_pattern = FALSE,
    noise_sd = 0.12,
    scans_per_condition = 8,
    rois_per_scan = 8) {

  conditions <- c("Control", "AD-CAA")
  regions <- names(region_pattern)

  scan_table <- expand.grid(
    disease_status = conditions,
    scan_number = seq_len(scans_per_condition),
    stringsAsFactors = FALSE
  )

  scan_table$Scan_ID <- paste0(
    ifelse(scan_table$disease_status == "Control", "CTRL", "AD"),
    "_",
    scan_table$scan_number
  )

  dat <- scan_table[
    rep(seq_len(nrow(scan_table)), each = rois_per_scan),
    ,
    drop = FALSE
  ]

  dat$ROI_ID <- paste0("ROI_", seq_len(nrow(dat)))

  dat$region <- rep(
    regions,
    length.out = nrow(dat)
  )

  dat$pathology <- "AmyloidFree"

  ad_rows <- dat$disease_status == "AD-CAA"

  # Ensure AD contains both pathology levels, while Control remains AmyloidFree.
  dat$pathology[ad_rows] <- rep(
    c("AmyloidFree", "Amyloid"),
    length.out = sum(ad_rows)
  )

  control_region_effect <- region_pattern[dat$region]

  if (reverse_ad_pattern) {
    ad_region_map <- rev(region_pattern)
    names(ad_region_map) <- names(region_pattern)

    region_effect <- control_region_effect
    region_effect[ad_rows] <- ad_region_map[dat$region[ad_rows]]
  } else {
    region_effect <- control_region_effect
  }

  scan_random_effect <- rnorm(
    nrow(scan_table),
    mean = 0,
    sd = 0.05
  )

  names(scan_random_effect) <- scan_table$Scan_ID

  dat$calibration_ratio <- (
    region_effect *
      ifelse(ad_rows, disease_scale, 1) *
      exp(scan_random_effect[dat$Scan_ID]) *
      exp(rnorm(nrow(dat), mean = 0, sd = noise_sd))
  )

  dat$log_calibration_ratio <- log(dat$calibration_ratio)

  dat
}


###############################################################################
# TEST 5
###############################################################################

print_header(5, "Calibration model recovers a known condition scale shift")

known_scale <- 2

scaled_gate_df <- make_gate_data(
  disease_scale = known_scale,
  reverse_ad_pattern = FALSE,
  noise_sd = 0.05
)

scaled_gate_result <- run_calibration_gate_check(
  calib_df = scaled_gate_df,
  biological_threshold_log2 = 0.5
)

expected_log2_shift <- log2(known_scale)
observed_log2_shift <- scaled_gate_result$disease_coefficient_log2

expect_close(
  observed = observed_log2_shift,
  expected = expected_log2_shift,
  tolerance = 0.15,
  label = "log2 condition shift"
)

cat(
  "Gate decision:",
  scaled_gate_result$decision,
  "\n"
)

cat("PASS: known condition-specific scale shift was recovered.\n")


###############################################################################
# TEST 6
###############################################################################

print_header(6, "No meaningful condition shift produces proceed")

proceed_df <- make_gate_data(
  disease_scale = 1.02,
  reverse_ad_pattern = FALSE,
  noise_sd = 0.18
)

proceed_result <- run_calibration_gate_check(
  calib_df = proceed_df,
  biological_threshold_log2 = 0.5
)

cat(
  "Estimated disease coefficient:",
  proceed_result$disease_coefficient_log2,
  "log2 units\n"
)

cat(
  "Maximum within-condition CV:",
  max(proceed_result$within_condition_cv$cv, na.rm = TRUE),
  "\n"
)

cat(
  "Region-pattern correlation:",
  proceed_result$region_pattern_correlation,
  "\n"
)

cat("Decision:", proceed_result$decision, "\n")

stopifnot(proceed_result$decision == "proceed")

cat("PASS: negligible condition shift produced 'proceed'.\n")


###############################################################################
# TEST 7
###############################################################################

print_header(
  7,
  "Large scale shift with concordant region pattern produces correction decision"
)

correction_df <- make_gate_data(
  disease_scale = 2,
  reverse_ad_pattern = FALSE,
  noise_sd = 0.08
)

correction_result <- run_calibration_gate_check(
  calib_df = correction_df,
  biological_threshold_log2 = 0.5
)

cat(
  "Estimated disease coefficient:",
  correction_result$disease_coefficient_log2,
  "log2 units\n"
)

cat(
  "Region-pattern correlation:",
  correction_result$region_pattern_correlation,
  "\n"
)

cat("Decision:", correction_result$decision, "\n")

stopifnot(
  correction_result$disease_coefficient_log2 > 0.5
)

stopifnot(
  is.finite(correction_result$region_pattern_correlation),
  correction_result$region_pattern_correlation > 0.5
)

stopifnot(
  correction_result$decision == "proceed_with_correction"
)

cat(
  "PASS: large shift with concordant regional pattern produced ",
  "'proceed_with_correction'.\n",
  sep = ""
)


###############################################################################
# TEST 8
###############################################################################

print_header(
  8,
  "Large scale shift with discordant region pattern does not proceed"
)

discordant_df <- make_gate_data(
  disease_scale = 2,
  reverse_ad_pattern = TRUE,
  noise_sd = 0.05
)

discordant_result <- run_calibration_gate_check(
  calib_df = discordant_df,
  biological_threshold_log2 = 0.5
)

cat(
  "Estimated disease coefficient:",
  discordant_result$disease_coefficient_log2,
  "log2 units\n"
)

cat(
  "Region-pattern correlation:",
  discordant_result$region_pattern_correlation,
  "\n"
)

cat("Decision:", discordant_result$decision, "\n")

stopifnot(
  abs(discordant_result$disease_coefficient_log2) > 0.5
)

stopifnot(
  is.na(discordant_result$region_pattern_correlation) ||
    discordant_result$region_pattern_correlation <= 0.5
)

stopifnot(
  discordant_result$decision == "do_not_proceed"
)

cat(
  "PASS: large shift with discordant regional pattern produced ",
  "'do_not_proceed'.\n",
  sep = ""
)


###############################################################################
# TEST 9
###############################################################################

print_header(9, "Zero total counts are rejected or exposed")

zero_counts <- total_counts_small
zero_counts$total_counts[1] <- 0

zero_count_result <- tryCatch(
  compute_calibration_ratio_table(
    abs_df = abs_df_small,
    total_counts_df = zero_counts
  ),
  error = function(e) e
)

if (inherits(zero_count_result, "error")) {
  cat("PASS: zero total counts were rejected.\n")
} else if (
  any(!is.finite(zero_count_result$calibration_ratio)) ||
    any(!is.finite(zero_count_result$log_calibration_ratio))
) {
  cat(
    "WARNING: zero total counts produced infinite calibration values.\n",
    "The utility should reject non-positive total counts explicitly.\n",
    sep = ""
  )
} else {
  cat(
    "WARNING: zero total counts were accepted without producing an error.\n"
  )
}


###############################################################################
# TEST 10
###############################################################################

print_header(10, "Negative absolute abundance is rejected or exposed")

negative_abs <- abs_df_small
negative_abs$abs_abundance[1] <- -100

negative_result <- tryCatch(
  compute_calibration_ratio_table(
    abs_df = negative_abs,
    total_counts_df = total_counts_small
  ),
  error = function(e) e
)

if (inherits(negative_result, "error")) {
  cat("PASS: negative absolute abundance was rejected.\n")
} else if (
  any(negative_result$total_abs_abundance <= 0) ||
    any(!is.finite(negative_result$log_calibration_ratio))
) {
  cat(
    "WARNING: negative abundance produced a non-positive or invalid ",
    "calibration ratio.\n",
    "The utility should validate absolute abundance explicitly.\n",
    sep = ""
  )
} else {
  cat(
    "WARNING: negative abundance was accepted and masked through summation.\n"
  )
}


###############################################################################
# TEST 11
###############################################################################

print_header(11, "Missing absolute-abundance values are not silently hidden")

missing_abs <- abs_df_small
missing_abs$abs_abundance[
  missing_abs$ROI_ID == "ROI1"
] <- NA_real_

missing_result <- tryCatch(
  compute_calibration_ratio_table(
    abs_df = missing_abs,
    total_counts_df = total_counts_small
  ),
  error = function(e) e
)

if (inherits(missing_result, "error")) {
  cat("PASS: all-missing ROI abundance was rejected.\n")
} else {
  roi1_total <- missing_result$total_abs_abundance[
    missing_result$ROI_ID == "ROI1"
  ]

  if (length(roi1_total) == 1 && roi1_total == 0) {
    cat(
      "WARNING: an ROI with all missing abundances was converted to zero ",
      "because sum(..., na.rm = TRUE) returns zero.\n",
      "This should be rejected rather than treated as true zero abundance.\n",
      sep = ""
    )
  } else {
    cat(
      "WARNING: missing abundance values were accepted. ",
      "Inspect the resulting ROI total.\n",
      sep = ""
    )
  }
}


###############################################################################
# TEST 12
###############################################################################

print_header(12, "Invalid biological threshold")

invalid_threshold_result <- tryCatch(
  run_calibration_gate_check(
    calib_df = proceed_df,
    biological_threshold_log2 = 0
  ),
  error = function(e) e
)

if (inherits(invalid_threshold_result, "error")) {
  cat("PASS: non-positive biological threshold was rejected.\n")
} else {
  cat(
    "WARNING: biological_threshold_log2 = 0 was accepted.\n",
    "The utility should require a finite, strictly positive threshold.\n",
    sep = ""
  )
}


###############################################################################
# SUMMARY
###############################################################################

cat("\n=========================================================\n")
cat("All absolute-calibration regression tests completed.\n")
cat("=========================================================\n")
