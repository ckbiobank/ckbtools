#' Calculate rates and their standard errors using various methods
#'
#' @param data A data frame.
#' @param per_pyar Crude and standardised incidence rates are given per `per_pyar` person-years at risk. (Default: 1000)
#' @param method Method to calculate confidence intervals, one of "norm" or "gamma". Default "norm" is Normal approximation.
#'
#' @return A data frame
#'
#' @keywords internal
#'
calc_rates <- function(data, per_pyar, method = c("norm", "gamma")) {
  method <- rlang::arg_match(method)

  if (method == "norm") {
    data %>%
      dplyr::summarise(crude_nevents = sum(.data$nevents),
                       crude_pyar    = sum(.data$pyar),
                       crude_ir      = sum(.data$nevents) / sum(.data$pyar) * per_pyar,
                       std_ir        = sum(.data$nevents / .data$pyar * .data$weight) * per_pyar,
                       var           = sum((.data$weight * .data$weight) * .data$nevents / .data$pyar ^ 2)) %>%
      dplyr::mutate(std_ir_se = sqrt(.data$var) * per_pyar) %>%
      dplyr::mutate(std_ir_lci = .data$std_ir - 1.96 * .data$std_ir_se,
                    std_ir_uci = .data$std_ir + 1.96 * .data$std_ir_se) %>%
      dplyr::select(-"var")
  } else if (method == "gamma") {
    q <- asht::wspoissonTest(data$nevents, data$weight / data$pyar, mult = per_pyar)
    tibble::tibble(crude_nevents = sum(data$nevents),
                   crude_pyar    = sum(data$pyar),
                   crude_ir      = sum(data$nevents) / sum(data$pyar) * per_pyar,
                   std_ir        = as.numeric(q$estimate),
                   std_ir_lci    = as.numeric(q$conf.int[[1]]),
                   std_ir_uci    = as.numeric(q$conf.int[[2]]))
  }
}



#' Calculate standardised incidence rates
#'
#'
#' @details
#' # Methods
#' See `vignette("incidence_rates")` for details.
#'
#' @inheritParams calc_rates
#' @param col.dob Name of column containing date of birth. (Must be of class `Date`.)
#' @param col.entry_date Name of column containing entry date. (Must be of class `Date`.)
#' @param col.exit_date Name of column containing exit date. (Must be of class `Date`.)
#' @param col.status Name of column containing status. This should be 1 if the event has occurred at the exit date, or 0 if it has not.
#' @param ages Vector of ages (years) to create age-at-risk strata. (Must be a numeric vector with distinct, increasing values. Default: `c(0, 100)`)
#' @param group Name of column containing groups to separately calculate incidence rates.
#' @param adjust_for Vector of names of columns containing adjustment variables. Include "agegrp" to adjust for age-at-risk strata.
#' @param weights A data frame with columns that match the `adjust_for` variables, and one row for each combination of levels of those variables, and a column named ‘std_pyar’ with the weights (person-years at risk for that stratum in the reference/standard population).
#' @return A named list:
#'  - incidence_rates_in_all_strata
#'  - standardised_incidence_rates
#'
#' @export
incidence_rates <- function(data,
                            col.dob        = "dob_anon",
                            col.entry_date = "study_date",
                            col.exit_date  = "endpoint_date",
                            col.status     = "endpoint",
                            ages           = c(0, 100),
                            group          = NULL,
                            adjust_for     = NULL,
                            per_pyar       = 1000,
                            weights        = NULL,
                            method         = c("norm", "gamma")) {

  method <- rlang::arg_match(method)

  # check arguments
  group_without_agegrp <-  group[group != "agegrp"]
  adjust_for_without_agegrp <- adjust_for[adjust_for != "agegrp"]
  column_names <- c(col.dob, col.entry_date, col.exit_date, col.status, group_without_agegrp, adjust_for_without_agegrp)

  ## reserved names
  reserved_names <- intersect(c("calendar", "age", "agegrp"), names(data))
  if (length(reserved_names) > 0) {
    cli::cli_abort("Column name{?s} {.field {reserved_names}} cannot be used in {.arg data}.")
  }

  ## missing columns
  if (length(setdiff(column_names, names(data))) != 0) {
    cli::cli_abort(c("{.field {setdiff(column_names, names(data))}} do{?es/} not exist in data.",
                     "i" = "Use arguments {.arg col.dob}, {.arg col.entry_date}, {.arg col.exit_date}, {.arg col.status}, {.arg group} and {.arg adjust_for} to set names of columns."))
  }

  ## format of date columns
  if (!inherits(data[[col.dob]], "Date") ||
      !inherits(data[[col.entry_date]], "Date") ||
      !inherits(data[[col.exit_date]], "Date")) {
    cli::cli_abort("{.field {col.dob}}, {.field {col.entry_date}} and {.field {col.exit_date}} must be dates.")
  }

  ## format of status column
  if (!is.numeric(data[[col.status]])) {
    cli::cli_abort("{.field {col.status}} must be numeric.")
  }

  ## format of group columns
  for (x in group_without_agegrp){
    if (!(is.factor(data[[x]]) || is.character(data[[x]]))) {
      if (!is.numeric(data[[x]])) {
        cli::cli_abort(c("{.field {x}} is a {tolower(class(data[[x]]))} column.",
                         i = "{.arg group} columns should be factor or character vectors (or numeric with few levels)."))
      }
      n_unique <- length(unique(data[[x]]))
      if (n_unique > 10) {
        cli::cli_abort(c("{.field {x}} is numeric and has {n_unique} level{?s}.",
                         i = "{.arg group} columns should be factor or character vectors (or numeric with few levels)."))
      }
      cli::cli_warn(c("{.field {x}} is numeric and has {n_unique} level{?s}, so treating as a group variable.",
                      i = "Change {.field {x}} to a factor to suppress this warning."))
    }
  }

  ## format of adjust_for columns
  for (x in adjust_for_without_agegrp){
    if (!(is.factor(data[[x]]) || is.character(data[[x]]))) {
      if (!is.numeric(data[[x]])) {
        cli::cli_abort(c("{.field {x}} is a {tolower(class(data[[x]]))} column.",
                         i = "{.arg adjust_for} columns should be factor or character vectors (or numeric with few levels)."))
      }
      n_unique <- length(unique(data[[x]]))
      if (n_unique > 10) {
        cli::cli_abort(c("{.field {x}} is numeric and has {n_unique} level{?s}.",
                         i = "{.arg adjust_for} columns should be factor or character vectors (or numeric with few levels)."))
      }
      cli::cli_warn(c("{.field {x}} is numeric and has {n_unique} level{?s}, so treating as a group variable.",
                      i = "Change {.field {x}} to a factor to suppress this warning."))
    }
  }

  ## missing values
  cols_with_missing <- column_names[sapply(data[column_names], function(col) any(is.na(col)))]
  if (length(cols_with_missing) > 0) {
    cli::cli_abort("{.field {cols_with_missing}} contain{?s/} missing values.")
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
    cli::cli_abort("{.arg ages} must a numeric vector with distinct, increasing values.")
  }

  ## asht package
  if (method == "gamma" && !requireNamespace("asht", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg asht} must be installed for the gamma method for confidence intervals.")
  }

  # end of checking arguments


  if (is.null(adjust_for)) {
    data$all <- 1
    adjust_for <- "all"
  }

  data_in_years <- data %>%
    dplyr::mutate(age_entry  = Epi::cal.yr(.data[[col.entry_date]]) - Epi::cal.yr(.data[[col.dob]]),
                  age_exit   = Epi::cal.yr(.data[[col.exit_date]]) - Epi::cal.yr(.data[[col.dob]]),
                  status     = .data[[col.status]],
                  entry_date = Epi::cal.yr(.data[[col.entry_date]]),
                  exit_date  = Epi::cal.yr(.data[[col.exit_date]]))

  lexis_data <- Epi::Lexis(
    entry        = list("calendar" = entry_date, "age" = age_entry),
    exit         = list("calendar" = exit_date),
    entry.status = 0,
    exit.status  = status,
    data         = data_in_years
  )

  expanded <- Epi::splitLexis(lexis_data, breaks = ages, time.scale = "age") %>%
    dplyr::mutate(agegrp = cut(.data$age,
                               breaks = ages,
                               right  = FALSE)) %>%
    dplyr::filter(!is.na(.data$agegrp))


  # calculate PYAR, no. of events, and incidence rate in groups formed by
  # group column and all strata columns
  ir_in_all_strata <- expanded %>%
    dplyr::group_by(!!!rlang::syms(c(group, adjust_for))) %>%
    dplyr::summarise(pyar    = sum(.data$lex.dur),
                     nevents = sum(.data$lex.Xst),
                     .groups = "drop") %>%
    dplyr::mutate(incidencerate = .data$nevents / .data$pyar)

  # calculate PYAR in groups formed by strata columns (standard population = all individuals)
  # or load from data set
  if (is.null(weights)) {
    pyar_across_strata <- ir_in_all_strata %>%
      dplyr::group_by(!!!rlang::syms(adjust_for)) %>%
      dplyr::summarise(std_pyar = sum(.data$pyar),
                       .groups = "drop")
  } else {
    pyar_across_strata <- weights
  }

  total_std_pyar <- sum(pyar_across_strata$std_pyar)

  # calculate standardised incidence rates in groups formed by group column
  ir_in_all_strata <- pyar_across_strata %>%
    tidyr::crossing(ir_in_all_strata[, group]) %>%
    dplyr::full_join(ir_in_all_strata, by = c(group, adjust_for)) %>%
    dplyr::mutate(pyar = dplyr::coalesce(.data$pyar, 0))

  if (any(ir_in_all_strata$pyar == 0)) {
    cli::cli_warn("There are strata with no time at risk. Returning only incidence rates within strata.")
    return(list(incidence_rates_in_all_strata = ir_in_all_strata))
  }

  ir_in_all_strata <- ir_in_all_strata %>%
    dplyr::mutate(weight = .data$std_pyar / total_std_pyar,
                  nevents = dplyr::coalesce(.data$nevents, 0)) %>%
    dplyr::arrange(dplyr::across(tidyselect::all_of(group)))

  if (any(ir_in_all_strata$nevents < 1)) {
    cli::cli_warn("There are strata with zero events.")
  } else  if (min(ir_in_all_strata$nevents) < 20) {
    cli::cli_warn("There are strata with fewer than 20 events.")
  }


  if (is.null(group)) {
    std_ir <- calc_rates(ir_in_all_strata, per_pyar, method)
  } else {
    std_ir <- ir_in_all_strata %>%
      dplyr::group_by(!!!rlang::syms(group)) %>%
      tidyr::nest() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(res = list(calc_rates(data, per_pyar, method))) %>%
      tidyr::unnest(cols = "res") %>%
      dplyr::ungroup() %>%
      dplyr::select(-"data")
  }

  return(list(incidence_rates_in_all_strata = ir_in_all_strata,
              standardised_incidence_rates = std_ir))
}
