#' Contrast utilities for spatial proportion models
#' @keywords internal
NULL

#' Prepare contrast variables
#'
#' @param df Long-format dataframe
#' @param disease_col Column with disease status
#' @param pathology_col Column with pathology
#' @param region_col Column with region
#'
#' @return Dataframe with factors
#' @export
prepare_contrasts <- function(df,
                              disease_col = "disease_status",
                              pathology_col = "pathology",
                              region_col = "region") {

  df[[disease_col]] <- factor(df[[disease_col]])
  df[[pathology_col]] <- factor(df[[pathology_col]])
  df[[region_col]] <- factor(df[[region_col]])

  return(df)
}

#' Create weighted amyloid burden (optional)
#'
#' @param df Dataframe
#' @param amyloid_col Column with amyloid score
#'
#' @return Updated dataframe
#' @export
add_weighted_amyloid <- function(df,
                                 amyloid_col = "amyloid_score") {

  if (!amyloid_col %in% colnames(df)) {
    warning("Amyloid column not found")
    return(df)
  }

  df$amyloid_scaled <- scale(df[[amyloid_col]])

  return(df)
}

