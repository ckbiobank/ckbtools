#' Reduce the size of model objects
#'
#' @keywords internal
reduce_model_object <- function(x) {
  x$data              <- NULL
  x$model             <- NULL
  x$linear.predictors <- NULL
  x$y                 <- NULL
  x$prior.weights     <- NULL
  x$weights           <- NULL
  x$residuals         <- NULL
  x$fitted.values     <- NULL
  x$family            <- NULL
  x$qr$qr             <- NULL
  attr(x$terms, ".Environment")   <- NULL
  attr(x$formula, ".Environment") <- NULL
  return(x)
}


#' Convert error/warning objects to character
#'
#' @keywords internal
to_character_vector <- function(obj) {
  if (is.null(obj) || length(obj) == 0) {
    return(NA_character_)
  }
  if (identical(as.character(obj), "")) {
    return(NA_character_)
  }
  return(paste(as.character(obj), collapse = "; "))
}


#' Create a tibble from the output of linear_model_safely(),
#' logistic_model_safely(), cox_model_safely()
#'
#' @keywords internal
build_nested_tibble_from_safely_output <- function(z) {
  result <- tibble::tibble(error = NA,
                           warnings = NA,
                           messages = NA,
                           output = NA,
                           model = list(z$result))
  ## add error, warnings, messages and output
  result$error <- z$error
  result$warnings <- z$warnings
  result$messages <- z$messages
  result$output <- z$output
  return(result)
}

#' @keywords internal
build_tibble_from_safely_output <- function(z) {
  result <- if (is.null(z$result)) data.frame(term = NA) else z$result
  result <- tibble::tibble(error = NA,
                           warnings = NA,
                           messages = NA,
                           output = NA,
                           result)
  ## add error, warnings, messages and output
  ## to every row of the result data frame
  result$error <- z$error
  result$warnings <- z$warnings
  result$messages <- z$messages
  result$output <- z$output
  return(result)
}


#' Find terms in a model that relate to exposure
#'
#' @keywords internal
find_exposure_terms <- function(data, model, exposure, exposure_term) {

  xlev <- model$xlevels[[exposure]]
  if (!is.null(xlev) && !is.ordered(model$model[[exposure]])) {
    return(paste0(exposure, model$xlevels[[exposure]]))
  }

  model_matrix <- stats::model.matrix(stats::delete.response(model$terms), data)
  assign <- attr(model_matrix, "assign")
  term_labels <- attr(model$terms, "term.labels")
  term_index <- which(
    vapply(
      term_labels,
      function(x) exposure %in% all.vars(stats::as.formula(paste("~", x))),
      logical(1)
    )
  )
  return(colnames(model_matrix)[assign %in% term_index])
}


#' organise coefs data frame
#' either keep only exposure terms, or move them to the top
#' @keywords internal
organise_coefs <- function(return, coefs, exposure_terms) {
  coefs <- coefs[order(match(coefs$term,
                             exposure_terms,
                             nomatch = length(exposure_terms) + 1)), ]
  if (return == "exposure") {
    coefs <- coefs[coefs$term %in% exposure_terms, ]
  }
  return(coefs)
}


#' Epi::float as a tibble
#'
#' @keywords internal
get_epi_float <- function(model, exposure) {
  far <- Epi::float(model)
  result <- tibble::tibble(
    term = paste0(exposure, names(far$coef)),
    float_estimate = as.numeric(far$coef),
    float_se = sqrt(far$var),
    limits = NA_real_,
    row.names = NULL
  )
  result$limits[1:2] <- far$limits
  return(result)
}


#' qvcalc::qvcalc::float as a tibble
#'
#' @keywords internal
get_qvcalc <- function(model, exposure, exposure_terms) {
  qv <- qvcalc::qvcalc(model, factorname = exposure)
  result <- tibble::tibble(
    term = exposure_terms,
    qvcalc_estimate = as.numeric(qv$qvframe$estimate),
    quasi_se = as.numeric(qv$qvframe$quasiSE),
    worstErrors = NA_real_,
    row.names = NULL
  )
  result$worstErrors[1:2] <- qvcalc::worstErrors(qv)
  return(result)
}
