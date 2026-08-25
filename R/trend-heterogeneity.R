#' Test for trend
#'
#' A chi-squared test for trend (on 1 degree of freedom) given a vector of independent estimates and their standard errors.
#'
#'
#' @details
#' See `vignette("trend_heterogeneity")` for details.
#'
#' @param beta vector of independent estimates
#' @param se vector of their standard errors
#'
#' @return A named list:
#'  - test statistic (chi-squared test statistic)
#'  - p (p-value)
#'
#' @export

trend <- function(beta, se) {
  # check arguments
  if (!is.numeric(beta) || !is.numeric(se)) {
    cli::cli_abort("{.var beta} and {.var se} must be numeric vectors.")
  }
  if (length(beta) != length(se)) {
    cli::cli_abort("{.var beta} and {.var se} must have same length.")
  }
  if (length(beta) < 2) {
    cli::cli_abort("At least two estimates are required.")
  }
  if (any(is.na(beta)) || any(is.na(se))) {
    cli::cli_abort("{.var beta} and {.var se} must not contain missing values.")
  }
  if (!all(se > 0)) {
    cli::cli_abort("Standard errors must be greater than zero.")
  }

  z <- seq_along(beta)
  w <- 1 / (se^2)
  test <- (sum(w * beta * (z - (sum(z * w) / sum(w))))^2) / sum(w * ((z - sum(z * w) / sum(w))^2))
  p <- stats::pchisq(test, df = 1, lower = FALSE)
  return(list("test statistic" = test, "p" = p))
}




#' Test for heterogeneity
#'
#' A chi-squared test (on \eqn{k - 1} degrees of freedom) for heterogeneity given a vector of \eqn{k} independent estimates and their standard errors.
#'
#' @details
#' See `vignette("trend_heterogeneity")` for details.
#'
#' @param beta vector of independent estimates
#' @param se vector of their standard errors
#'
#' @return A named list:
#'  - test statistic (chi-squared test statistic)
#'  - degrees of freedom (number of degrees of freedom used for the test, equal to `length(beta)-1`)
#'  - p (p-value)
#'
#' @export

heterogeneity <- function(beta, se) {
  # check arguments
  if (!is.numeric(beta) || !is.numeric(se)) {
    cli::cli_abort("{.var beta} and {.var se} must be numeric vectors.")
  }
  if (length(beta) != length(se)) {
    cli::cli_abort("{.var beta} and {.var se} must have same length.")
  }
  if (length(beta) < 2) {
    cli::cli_abort("At least two estimates are required.")
  }
  if (any(is.na(beta)) || any(is.na(se))) {
    cli::cli_abort("{.var beta} and {.var se} must not contain missing values.")
  }
  if (!all(se > 0)) {
    cli::cli_abort("Standard errors must be greater than zero.")
  }

  # degrees of freedom
  df <- length(beta) - 1

  # expected_beta: inverse variance weighted average of betas
  expected_beta <- sum(beta / se^2) / sum(1 / se^2)

  heterogeneity_test_statistic <- sum(((beta - expected_beta) / se)^2)

  p <- stats::pchisq(heterogeneity_test_statistic, df = df, lower = FALSE)

  return(list("test statistic" = heterogeneity_test_statistic,
              "degrees of freedom" = df,
              "p" = p))
}
