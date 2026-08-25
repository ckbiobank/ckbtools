#' Lexis expansion of data by age-at-risk
#'
#' @param data A data frame
#' @param col.id Name of column containing unique participant IDs
#' @param col.dob Name of column containing date of birth (Must be of class `Date`)
#' @param col.entry_date Name of column containing entry date (Must be of class `Date`)
#' @param col.exit_date Name of column containing exit date (Must be of class `Date`)
#' @param col.status Name of column containing status
#' @param col.keep Names of additional columns to include in the result (e.g. exposure and covariate variables)
#' @param ages Vector of ages (years) to create age-at-risk strata. Any participant follow-up outside these ages is dropped. (Must be a numeric vector with distinct, increasing values. Default: c(0, 100))
#' @param age.labels Vector of name for age-at-risk strata (Default: `paste0(ages[-length(ages)], "-", ages[-1])`)
#' @param null.status Value to use in status column for intervals before the final interval (Default: 0)
#'
#' @return Data frame expanded by age-at-risk
#'
#' @export

expand_age_at_risk <- function(data,
                               col.id         = "csid",
                               col.dob        = "dob_anon",
                               col.entry_date = "study_date",
                               col.exit_date  = "endpoint_date",
                               col.status     = "endpoint",
                               col.keep       = NULL,
                               ages           = c(0, 100),
                               age.labels     = paste0(ages[-length(ages)], "-", ages[-1]),
                               null.status    = 0) {

  # check arguments
  column_names <- c(col.id, col.dob, col.entry_date, col.exit_date, col.status)

  ## missing columns
  if (length(setdiff(column_names, names(data))) != 0) {
    cli::cli_abort(c("{.field {setdiff(column_names, names(data))}} do{?es/} not exist in data.",
                     "i" = "Use arguments {.arg col.id}, {.arg col.dob}, {.arg col.entry_date}, {.arg col.exit_date} and {.arg col.status} to set names of columns."))
  }

  if (length(setdiff(col.keep, names(data))) != 0) {
    cli::cli_abort(c("{.field {setdiff(col.keep, names(data))}} do{?es/} not exist in data."))
  }

  ## format of date columns
  if (!inherits(data[[col.dob]], "Date") ||
      !inherits(data[[col.entry_date]], "Date") ||
      !inherits(data[[col.exit_date]], "Date")) {
    cli::cli_abort("{.field {col.dob}}, {.field {col.entry_date}} and {.field {col.exit_date}} must be dates.")
  }

  ## missing values
  cols_with_missing <- column_names[1:4][sapply(data[column_names[1:4]], function(col) any(is.na(col)))]
  if (length(cols_with_missing) > 0) {
    cli::cli_abort("{.field {cols_with_missing}} contain{?s/} missing values.")
  }

  if (any(is.na(data[[col.status]]))) {
    cli::cli_warn("{.field {col.status}} contains missing values, these observations will be dropped.")
  }

  ## non-integer dates
  if (any(floor(as.numeric(data[[col.dob]])) != as.numeric(data[[col.dob]])) ||
      any(floor(as.numeric(data[[col.entry_date]])) != as.numeric(data[[col.entry_date]])) ||
      any(floor(as.numeric(data[[col.exit_date]])) != as.numeric(data[[col.exit_date]]))) {
    cli::cli_warn("Data contains non-integer dates.")
  }

  ## ordering of dates
  if (any(data[[col.entry_date]] > data[[col.exit_date]])) {
    cli::cli_abort("{.field {col.exit_date}} must not be before {.field {col.entry_date}}.")
  }

  if (any(data[[col.dob]] > data[[col.entry_date]])) {
    cli::cli_abort("{.field {col.entry_date}} must not be before {.field {col.dob}}.")
  }

  ## ages vector
  if (!is.numeric(ages) || !identical(ages, sort(unique(ages)))) {
    cli::cli_abort("{.arg ages} must be a numeric vector with distinct, increasing values.")
  }

  # end of checking arguments

  # remove grouping
  if (dplyr::is_grouped_df(data)) {
    cli::cli_warn("Removed grouping by {.field {dplyr::group_vars(data)}}")
    data <- dplyr::ungroup(data)
  }


  expanded_data <- data %>%

    # keep only needed columns, to reduce memory requirements
    dplyr::select(tidyselect::all_of(c(col.id, col.dob, col.entry_date, col.exit_date, col.status))) %>%

    # remove incomplete
    na.omit() %>%

    # create row for each participant and each age group
    tidyr::crossing(tibble::tibble(agegrp_start = ages[-length(ages)],
                                   agegrp_end = ages[-1],
                                   agegrp = age.labels)) %>%

    # start date of each interval is birthday (i.e. calendar years after dob)
    dplyr::mutate(start_int = lubridate::add_with_rollback(.data[[col.dob]],
                                                           lubridate::years(.data$agegrp_start),
                                                           roll_to_first = TRUE)) %>%

    # remove intervals that start after censoring date
    dplyr::filter(.data$start_int <= .data[[col.exit_date]]) %>%

    # end date of each interval is day before next interval
    dplyr::mutate(end_int = lubridate::add_with_rollback(.data[[col.dob]],
                                                         lubridate::years(.data$agegrp_end),
                                                         roll_to_first = TRUE) - 1) %>%

    # remove intervals that end before entering study
    dplyr::filter(.data$end_int >= .data[[col.entry_date]]) %>%

    # correct start and end date for intervals partially during study
    dplyr::mutate(start_int = pmax(.data[[col.entry_date]], .data$start_int),
                  end_int   = pmin(.data[[col.exit_date]], .data$end_int)) %>%


    # endpoint can only be 1 for the last interval
    dplyr::mutate(status = dplyr::if_else(.data$end_int == .data[[col.exit_date]],
                                          .data[[col.status]],
                                          .env$null.status),
                  # calculate days / years from start of study
                  ## end intervals at + 0.95
                  t_start_days   = lubridate::interval(.data[[col.entry_date]], .data$start_int) %/% lubridate::days(1),
                  t_end_days     = (lubridate::interval(.data[[col.entry_date]], .data$end_int) %/% lubridate::days(1)) + 0.95,
                  time_in        = .data$t_start_days / 365.25,
                  time_out       = .data$t_end_days / 365.25) %>%

    dplyr::left_join(dplyr::select(data, tidyselect::all_of(c(col.id, col.keep))), by = col.id)

  return(expanded_data)
}
