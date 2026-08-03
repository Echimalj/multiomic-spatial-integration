# ============================================================
# AGORA analysis tiers
#
# The full AGORA Wall of Targets contains hundreds of nominated
# genes. Running the complete gene-panel workflow on all targets
# would generate:
#
# - thousands of spatial contrasts,
# - thousands of lineage-interaction models,
# - and tens of thousands of fine-cell-type mixed models.
#
# For the primary AGORA analysis, we therefore define a focused
# Tier 1 panel containing targets with at least three independent
# nominations.
#
# This threshold:
#
# - prioritizes repeatedly nominated AD targets;
# - reduces the multiple-testing burden;
# - keeps the analysis comparable in size to the Hajjar panel;
# - and preserves all targets with the strongest cross-team support.
#
# In the current AGORA download, Tier 1 contains 38 targets (03 Aug2026) :
#
# - 4 targets with 5 nominations;
# - 10 targets with 4 nominations;
# - 24 targets with 3 nominations.
#
# All 38 Tier 1 targets are evaluable in the current GeoMx and
# inferred-reference datasets.
# ============================================================


# ============================================================
# AGORA nominated-target utilities
# ============================================================

#' Download nominated Alzheimer's disease targets from AGORA
#'
#' Queries the public AGORA comparison-tools API page by page and returns
#' one row per nominated target. Nested metadata fields are collapsed into
#' semicolon-delimited strings so the output can be written to CSV.
#'
#' @param page_size Number of records requested per API page.
#' @param max_retries Number of attempts per page.
#' @param retry_wait Seconds to wait between retries.
#' @param verbose Print download progress.
#'
#' @return Data frame containing AGORA target metadata and the standard
#'   gene-panel columns required by the validation pipeline.
#' @export
download_agora_targets <- function(
    page_size = 100L,
    max_retries = 3L,
    retry_wait = 2,
    verbose = TRUE
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "download_agora_targets: package `jsonlite` is required.",
      call. = FALSE
    )
  }

  page_size <- as.integer(page_size)
  max_retries <- as.integer(max_retries)

  if (!is.finite(page_size) || page_size < 1L) {
    stop(
      "download_agora_targets: `page_size` must be a positive integer.",
      call. = FALSE
    )
  }

  if (!is.finite(max_retries) || max_retries < 1L) {
    stop(
      "download_agora_targets: `max_retries` must be a positive integer.",
      call. = FALSE
    )
  }

  api_base <- paste0(
    "https://agora.adknowledgeportal.org/",
    "api/v1/comparison-tools/targets"
  )

  build_url <- function(page_number) {
    paste0(
      api_base,
      "?itemFilterType=exclude",
      "&pageNumber=", page_number,
      "&pageSize=", page_size,
      "&sortFields=total_nominations",
      "&sortFields=hgnc_symbol",
      "&sortOrders=-1",
      "&sortOrders=1"
    )
  }

  read_page <- function(page_number) {
    request_url <- build_url(page_number)

    last_error <- NULL

    for (attempt in seq_len(max_retries)) {
      response <- try(
        jsonlite::fromJSON(
          request_url,
          simplifyVector = FALSE
        ),
        silent = TRUE
      )

      if (!inherits(response, "try-error")) {
        return(response)
      }

      last_error <- response

      if (attempt < max_retries) {
        if (verbose) {
          message(
            "AGORA page ",
            page_number,
            " failed on attempt ",
            attempt,
            "; retrying..."
          )
        }

        Sys.sleep(retry_wait)
      }
    }

    stop(
      "download_agora_targets: failed to retrieve page ",
      page_number,
      " after ",
      max_retries,
      " attempts. Last error: ",
      as.character(last_error),
      call. = FALSE
    )
  }

  collapse_field <- function(x) {
    if (is.null(x) || length(x) == 0L) {
      return(NA_character_)
    }

    values <- unique(
      trimws(
        as.character(
          unlist(x, use.names = FALSE)
        )
      )
    )

    values <- values[
      !is.na(values) &
        nzchar(values)
    ]

    if (length(values) == 0L) {
      NA_character_
    } else {
      paste(values, collapse = "; ")
    }
  }

  target_to_row <- function(target) {
    data.frame(
      agora_id = if (is.null(target[["_id"]])) {
        NA_character_
      } else {
        as.character(target[["_id"]])
      },
      ensembl_gene_id = if (is.null(target[["ensembl_gene_id"]])) {
        NA_character_
      } else {
        as.character(target[["ensembl_gene_id"]])
      },
      gene = if (is.null(target[["hgnc_symbol"]])) {
        NA_character_
      } else {
        as.character(target[["hgnc_symbol"]])
      },
      total_nominations = if (is.null(target[["total_nominations"]])) {
        NA_integer_
      } else {
        as.integer(target[["total_nominations"]])
      },
      initial_nomination = if (is.null(target[["initial_nomination"]])) {
        NA_integer_
      } else {
        as.integer(target[["initial_nomination"]])
      },
      nominating_teams = collapse_field(
        target[["nominating_teams"]]
      ),
      cohort_studies = collapse_field(
        target[["cohort_studies"]]
      ),
      input_data = collapse_field(
        target[["input_data"]]
      ),
      programs = collapse_field(
        target[["programs"]]
      ),
      pharos_class = if (is.null(target[["pharos_class"]])) {
        NA_character_
      } else {
        as.character(target[["pharos_class"]])
      },
      stringsAsFactors = FALSE
    )
  }

  if (verbose) {
    message("Requesting AGORA page 0...")
  }

  first_page <- read_page(0L)

  if (
    is.null(first_page[["page"]]) ||
      is.null(first_page[["page"]][["totalPages"]]) ||
      is.null(first_page[["page"]][["totalElements"]])
  ) {
    stop(
      "download_agora_targets: API response did not contain expected ",
      "pagination metadata.",
      call. = FALSE
    )
  }

  total_pages <- as.integer(
    first_page[["page"]][["totalPages"]]
  )

  expected_total <- as.integer(
    first_page[["page"]][["totalElements"]]
  )

  if (verbose) {
    message(
      "AGORA reports ",
      expected_total,
      " targets across ",
      total_pages,
      " page(s)."
    )
  }

  page_numbers <- seq.int(
    from = 0L,
    to = total_pages - 1L
  )

  page_responses <- vector(
    "list",
    length(page_numbers)
  )

  page_responses[[1]] <- first_page

  if (length(page_numbers) > 1L) {
    for (i in seq.int(
      from = 2L,
      to = length(page_numbers)
    )) {
      page_number <- page_numbers[[i]]
  
      if (verbose) {
        message(
          "Requesting AGORA page ",
          page_number,
          " of ",
          total_pages - 1L,
          "..."
        )
      }

      page_responses[[i]] <- read_page(
        page_number
      )
    }
  }

  target_records <- unlist(
    lapply(
      page_responses,
      function(page_response) {
        page_response[["nominatedTargets"]]
      }
    ),
    recursive = FALSE
  )

  if (length(target_records) == 0L) {
    stop(
      "download_agora_targets: no nominated targets were returned.",
      call. = FALSE
    )
  }

record_ids <- vapply(
  target_records,
  function(target) {
    value <- target[["_id"]]

    if (is.null(value)) {
      NA_character_
    } else {
      as.character(value)
    }
  },
  character(1)
)

n_duplicate_records <- sum(
  duplicated(record_ids) &
    !is.na(record_ids)
)

if (n_duplicate_records > 0L) {
  warning(
    "download_agora_targets: API pagination returned ",
    n_duplicate_records,
    " duplicated target record(s).",
    call. = FALSE
  )
}

  targets <- do.call(
    rbind,
    lapply(
      target_records,
      target_to_row
    )
  )

  rownames(targets) <- NULL

  targets$gene <- trimws(
    as.character(targets$gene)
  )

  targets <- targets[
    !is.na(targets$gene) &
      nzchar(targets$gene),
    ,
    drop = FALSE
  ]

  duplicated_symbols <- unique(
    targets$gene[
      duplicated(targets$gene)
    ]
  )

  if (length(duplicated_symbols) > 0L) {
    warning(
      "download_agora_targets: removing duplicated HGNC symbols: ",
      paste(
        utils::head(duplicated_symbols, 10),
        collapse = ", "
      ),
      call. = FALSE
    )

    targets <- targets[
      !duplicated(targets$gene),
      ,
      drop = FALSE
    ]
  }

  targets$panel <- "AGORA_Wall_of_Targets"
  targets$category <- "AGORA nominated target"
  targets$headline <- FALSE
  targets$source <- "AGORA AD Knowledge Portal"

  standard_columns <- c(
    "gene",
    "panel",
    "category",
    "headline",
    "source"
  )

  metadata_columns <- setdiff(
    names(targets),
    standard_columns
  )

  targets <- targets[
    ,
    c(
      standard_columns,
      metadata_columns
    ),
    drop = FALSE
  ]

  targets <- targets[
    order(
      -targets$total_nominations,
      targets$gene,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  rownames(targets) <- NULL

  if (nrow(targets) != expected_total) {
    warning(
      "download_agora_targets: API reported ",
      expected_total,
      " targets, but ",
      nrow(targets),
      " unique valid HGNC symbols were retained.",
      call. = FALSE
    )
  }

  attr(targets, "downloaded_at") <- format(
    Sys.time(),
    tz = "UTC",
    usetz = TRUE
  )

  attr(targets, "api_total_elements") <- expected_total

  targets
}

#' Create a high-nomination AGORA target panel
#'
#' Defines the primary AGORA analysis panel by retaining nominated
#' targets supported by at least `minimum_nominations` independent
#' nominations.
#'
#' The default threshold of three nominations is intended to prioritize
#' repeatedly nominated AD targets while keeping the downstream spatial,
#' lineage, and fine-cell-type analyses computationally and statistically
#' manageable.
#'
#' @param agora_targets Data frame returned by `download_agora_targets()`
#'   or loaded from `agora_wall_of_targets.csv`.
#' @param minimum_nominations Minimum number of AGORA nominations required.
#' @param panel_name Panel label written to the standard `panel` column.
#' @param headline_threshold Nomination count used to mark headline targets.
#'
#' @return Panel-compatible data frame containing the retained AGORA targets
#'   and all original AGORA metadata.
#' @export
create_agora_tier1_panel <- function(
    agora_targets,
    minimum_nominations = 3L,
    panel_name = "AGORA_High_Nomination",
    headline_threshold = 4L
) {
  if (!is.data.frame(agora_targets)) {
    stop(
      "create_agora_tier1_panel: `agora_targets` must be a data frame.",
      call. = FALSE
    )
  }

  required_cols <- c(
    "gene",
    "total_nominations"
  )

  missing_cols <- setdiff(
    required_cols,
    names(agora_targets)
  )

  if (length(missing_cols) > 0L) {
    stop(
      "create_agora_tier1_panel: missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  minimum_nominations <- as.integer(
    minimum_nominations
  )

  headline_threshold <- as.integer(
    headline_threshold
  )

  if (
    !is.finite(minimum_nominations) ||
      minimum_nominations < 1L
  ) {
    stop(
      "create_agora_tier1_panel: `minimum_nominations` must be a ",
      "positive integer.",
      call. = FALSE
    )
  }

  if (
    !is.finite(headline_threshold) ||
      headline_threshold < minimum_nominations
  ) {
    stop(
      "create_agora_tier1_panel: `headline_threshold` must be greater ",
      "than or equal to `minimum_nominations`.",
      call. = FALSE
    )
  }

  targets <- agora_targets

  targets$gene <- trimws(
    as.character(targets$gene)
  )

  targets$total_nominations <- suppressWarnings(
    as.integer(targets$total_nominations)
  )

  invalid_gene <- is.na(targets$gene) |
    !nzchar(targets$gene)

  if (any(invalid_gene)) {
    warning(
      "create_agora_tier1_panel: removing ",
      sum(invalid_gene),
      " row(s) with missing or empty gene symbols.",
      call. = FALSE
    )

    targets <- targets[
      !invalid_gene,
      ,
      drop = FALSE
    ]
  }

  tier1 <- targets[
    is.finite(targets$total_nominations) &
      targets$total_nominations >= minimum_nominations,
    ,
    drop = FALSE
  ]

  if (nrow(tier1) == 0L) {
    stop(
      "create_agora_tier1_panel: no targets met the nomination threshold ",
      "of ",
      minimum_nominations,
      ".",
      call. = FALSE
    )
  }

  duplicated_genes <- unique(
    tier1$gene[
      duplicated(tier1$gene)
    ]
  )

  if (length(duplicated_genes) > 0L) {
    warning(
      "create_agora_tier1_panel: removing duplicated gene symbols: ",
      paste(
        utils::head(duplicated_genes, 10),
        collapse = ", "
      ),
      call. = FALSE
    )

    tier1 <- tier1[
      !duplicated(tier1$gene),
      ,
      drop = FALSE
    ]
  }

  tier1$panel <- panel_name

  tier1$category <- paste0(
    "AGORA nominations: ",
    tier1$total_nominations
  )

  tier1$headline <-
    tier1$total_nominations >= headline_threshold

  if (!"source" %in% names(tier1)) {
    tier1$source <- "AGORA AD Knowledge Portal"
  }

  standard_columns <- c(
    "gene",
    "panel",
    "category",
    "headline",
    "source"
  )

  metadata_columns <- setdiff(
    names(tier1),
    standard_columns
  )

  tier1 <- tier1[
    ,
    c(
      standard_columns,
      metadata_columns
    ),
    drop = FALSE
  ]

  tier1 <- tier1[
    order(
      -tier1$total_nominations,
      tier1$gene,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  rownames(tier1) <- NULL

  message(
    "Created AGORA Tier 1 panel with ",
    nrow(tier1),
    " target(s) having at least ",
    minimum_nominations,
    " nominations."
  )

  nomination_distribution <- table(
    tier1$total_nominations,
    useNA = "ifany"
  )

  message(
    "Tier 1 nomination distribution: ",
    paste(
      names(nomination_distribution),
      as.integer(nomination_distribution),
      sep = "=",
      collapse = ", "
    ),
    "."
  )

  tier1
}
