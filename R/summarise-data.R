#' Summarise continuous variables
#'
#' @param data Data set
#' (e.g. produced by `prepare_proteomics_analysis_dataset()`)
#' @param cols Names of columns containing continuous variables to be
#' summarised.
#' @param by Names of columns that contain categorical participant
#' characteristics to produce subgroup summaries.
#'
#' @export
summarise_continuous <- function(data,
                                 cols,
                                 by = NULL) {

  calc_stats <- function(x) {
    quantiles <- as.numeric(stats::quantile(x, probs = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99), na.rm = TRUE))

    if (sum(!is.na(x)) < 1) {
      return(
        data.frame(
          n_nonmissing = sum(!is.na(x)),
          n_missing = sum(is.na(x)),
          percent_missing = 100 * sum(is.na(x)) / length(x)
        )
      )
    }

    data.frame(
      n_nonmissing = sum(!is.na(x)),
      n_missing = sum(is.na(x)),
      percent_missing = 100 * sum(is.na(x)) / length(x),
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      p1 = quantiles[1],
      p10 = quantiles[2],
      p25 = quantiles[3],
      p50 = quantiles[4],
      p75 = quantiles[5],
      p90 = quantiles[6],
      p99 = quantiles[7],
      max = max(x, na.rm = TRUE)
    )
  }

  get_summaries <- function(data) {
    purrr::list_rbind(purrr::map(rlang::set_names(cols), ~ calc_stats(data[[.]])),
                      names_to = "var")
  }

  result <- get_summaries(data)

  summary_statistics_by_subgroup <- NULL
  if (!is.null(by)) {
    grouped_data <- purrr::list_rbind(
      purrr::map(
        by,
        ~ tidyr::nest(dplyr::group_by(data, dplyr::across(tidyselect::all_of(.))))
      )
    )

    summary_statistics_by_subgroup <- grouped_data %>%
      dplyr::mutate(res = purrr::map(data, get_summaries)) %>%
      tidyr::unnest("res") %>%
      dplyr::select(-"data")

    result <- dplyr::bind_rows(result, summary_statistics_by_subgroup) %>%
      dplyr::relocate(tidyselect::all_of(by), .after = 1)  %>%
      dplyr::arrange(match(.data$var, cols))
  }

  return(result)
}
