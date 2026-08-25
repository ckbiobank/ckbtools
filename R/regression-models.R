
#' Estimate the association between an exposure and a continuous outcome
#'
#' This function is used to estimate the (adjusted) linear association
#' between an exposure variable and a continuous outcome variable using
#' a linear regression model.
#' It is a wrapper for the \link[stats]{lm} function for convenience
#' and provides results from some additional analysis steps if requested.
#'
#' @details
#' # Additional results
#' The `additional` argument is used to add results from additional analyses.
#' It must be a character vector containing one or more of:
#' \itemize{
#'  \item "anova": ANOVA F-test comparing the full model to one without the
#'  exposure. Output will have additional columns "F" (F test
#'  statistics), "Df" (degrees of freedom), "p_F" (p-value).
#'  \item "float": Standard errors for the exposure calculated by \link{float}.
#'  Output will have additional columns "float_estimate", "float_se" and "limits".
#'  \item "qvcalc": Standard errors for the exposure calculated by
#'  \link{qvcalc}. Output will have additional columns "qvcalc_estimate",
#'  "quasi_se" and "worstErrors".
#'  \item "emmeans": Estimated marginal means for the exposure calculated using
#'  \link[emmeans]{emmeans} with proportional weights for factor variables
#'  that are averaged over. Output will have additional columns "emmean",
#'  "emmean_se" and "emmean_quasi_se".
#'  \item "" (default): Nothing.
#'  }
#' These will only be included when relevant e.g. standard errors from
#' \link{float} when the exposure has more than two levels.
#' This only has an effect if `return` is "coefs" or "exposure".
#'
#' See \link{float} and \link{qvcalc} for details about these results. Note that
#' "limits" and "worstErrors" do not correspond to particular groups but
#' accuracy other all contrasts.
#'
#' @param data A data frame.
#' @param exposure Name of column containing exposure variable.
#' @param outcome Name of column containing outcome variable.
#' @param adjust_for Character vector of terms to include in the model formula.
#' @param exclude Names of columns that identify observations that should be
#' excluded.
#' @param exposure_term The term(s) to include in the model for the exposure.
#' This can be a character vector or a function that take the exposure name as
#' the argument and returns a character vector. (Default: exposure)
#' @param additional Extra results to include.
#' See "Additional results" for details.
#' @param return
#' What should the function return:
#' \itemize{
#'  \item "model" (default): The model object (object of class "lm").
#'  \item "reduced-model": The model object with some large components removed.
#'  \item "coefs": A data frame of all estimated coefficients and other key
#'  model statistics.
#'  \item "exposure": A data frame of estimated coefficient(s) for the exposure
#'  and other key model statistics.
#' }
#' @param return_fn
#' If return="model", name of function to apply to model before it's returned.
#' (Default: identity)
#' @param nuisance If "emmeans" is requested via the `additional` argument, use
#' this to specify predictors to omit from the reference grid. See the help for
#' \link[emmeans]{emmeans} and \link[emmeans]{ref_grid}.
#' @param ... Arguments passed to `stats::lm()`
#'
#' @returns Either a model object, reduced model object or data frame of
#' coefficients.
#'
#' @examples
#' linear_model(ckbtools_participant_data,
#'              exposure = "age_at_study_date",
#'              outcome = "sbp_mean",
#'              adjust_for = "sex",
#'              return = "coefs")
#'
#' # factor exposure, calculate estimated marginal means
#' linear_model(ckbtools_participant_data,
#'              exposure = "smoking",
#'              outcome = "sbp_mean",
#'              adjust_for = c("age_at_study_date", "sex"),
#'              additional = "emmeans",
#'              return = "coefs")
#'
#' # return the model object
#' linear_model(ckbtools_participant_data,
#'              exposure = "bmi_grp",
#'              outcome = "sbp_mean",
#'              adjust_for = c("age_at_study_date", "sex"),
#'              return = "model")
#'
#' @export
linear_model <- function(
  data,
  exposure,
  outcome,
  adjust_for    = NULL,
  exclude       = NULL,
  exposure_term = exposure,
  additional    = "",
  return        = c("model", "reduced-model", "coefs", "exposure"),
  return_fn     = identity,
  nuisance      = character(0),
  ...
) {

  # check arguments
  return <- rlang::arg_match(return)
  additional <- rlang::arg_match(additional,
                                 values = c("anova",
                                            "float",
                                            "qvcalc",
                                            "emmeans",
                                            ""),
                                 multiple = TRUE)

  # exclude observations
  if (!is.null(exclude)) {
    data <- data %>%
      dplyr::filter(!dplyr::if_any(tidyselect::any_of(exclude), ~ .x %in% TRUE))  # treats NA as FALSE
  }

  # linear regression model
  if (is.function(exposure_term)) {
    exposure_term <- exposure_term(exposure)
  }
  model_formula <- stats::reformulate(c(exposure_term, adjust_for), outcome)
  model <- stats::lm(formula = model_formula, data = data, ...)

  # return model
  if (return == "model") {
    return(return_fn(model))
  }

  # return reduced model
  if (return == "reduced-model") {
    return(reduce_model_object(model))
  }

  # get results from model
  coefs <- broom::tidy(model)
  model_frame <- stats::model.frame(model)
  exposure_terms <- find_exposure_terms(data, model, exposure, exposure_term)

  # Count observations per group
  if (!is.null(exposure) &&
      length(model$xlevels[[exposure]]) > 1) {
    counts <- table(model_frame[[exposure]])
    coefs <- dplyr::full_join(coefs,
                              data.frame(term = paste0(exposure, names(counts)),
                                         n_group = as.numeric(counts)),
                              by = "term")
  }

  # Epi::float standard errors, if model has intercept
  if ("float" %in% additional &&
      !is.null(exposure) &&
      !is.ordered(model$model[[exposure]]) &&
      length(model$xlevels[[exposure]]) > 2  &&
      attr(stats::terms(model), "intercept") == 1) {
    coefs <- dplyr::full_join(coefs,
                              get_epi_float(model, exposure),
                              by = "term")
  }

  # qvcalc::qvcalc standard errors
  if ("qvcalc" %in% additional &&
      !is.null(exposure) &&
      !is.ordered(model$model[[exposure]]) &&
      length(model$xlevels[[exposure]]) > 2) {
    coefs <- dplyr::full_join(coefs,
                              get_qvcalc(model, exposure, exposure_terms),
                              by = "term")
  }

  # Estimates marginal means
  if ("emmeans" %in% additional &&
      !is.null(exposure) &&
      !is.ordered(model$model[[exposure]]) &&
      length(model$xlevels[[exposure]]) > 1) {

    if (!requireNamespace("emmeans", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg emmeans} must be installed for emmeans.")
    }

    mm <- emmeans::emmeans(model,
                           specs = exposure,
                           weights = "proportional",
                           nuisance = nuisance)
    summary_mm <- summary(mm)

    marginal_means <- tibble::tibble(
      term = paste0(exposure, summary_mm[, 1]),
      emmean = summary_mm$emmean,
      emmean_se = summary_mm$SE,
      emmean_quasi_se = NA_real_,
      row.names = NULL
    )

    if (length(model$xlevels[[exposure]]) > 2) {
      qvcalc_mm <- qvcalc(stats::vcov(mm),
                          estimates = summary_mm$emmean)
      marginal_means$emmean_quasi_se <- as.numeric(qvcalc_mm$qvframe$quasiSE)
    }


    coefs <- dplyr::full_join(coefs, marginal_means, by = "term")
  }

  # organise coefs
  coefs <- organise_coefs(return, coefs, exposure_terms)

  # either keep only exposure terms, or move them to the top
  if (return == "exposure") {
    coefs <- coefs[coefs$term %in% exposure_terms, ]
  } else {
    coefs <- coefs[order(match(coefs$term,
                               exposure_terms,
                               nomatch = length(exposure_terms) + 1)), ]
  }

  # ANOVA
  if ("anova" %in% additional) {
    null_model_formula <- stats::reformulate(adjust_for %||% "1",
                                             outcome)
    null_model <- stats::lm(null_model_formula,
                            data = data[!is.na(data[[exposure]]), ],
                            ...)
    anova_result <- stats::anova(null_model, model)
    coefs[, c("F", "Res.Df", "Df", "p_F")] <- NA
    coefs$F[[1]] <- anova_result$F[[2]]
    coefs$Res.Df[[1]] <- anova_result$Res.Df[[2]]
    coefs$Df[[1]] <- anova_result$Df[[2]]
    coefs$p_F[[1]] <- anova_result$`Pr(>F)`[[2]]
  }

  # Number of observations and model formula
  coefs[, c("n", "formula")] <- NA
  coefs$n[[1]] <- nrow(model_frame)
  coefs <- dplyr::relocate(coefs, tidyselect::any_of("n_group"), .after = "n")
  coefs$formula[[1]] <- deparse1(model_formula)

  if ("anova" %in% additional) {
    coefs[, "null_model_formula"] <- NA
    coefs$null_model_formula[[1]] <- deparse1(null_model_formula)
  }

  # return coefficients
  return(coefs)
}


#' Estimate the association between an exposure and a binary outcome
#'
#' This function is used to estimate the (adjusted) association (log odds ratio)
#' between an exposure variable and a binary outcome variable using
#' a logistic regression model.
#' It is a wrapper for the \link[stats]{glm} function with `family = "binomial"`
#' for convenience and provides results from some additional analysis steps if
#' requested.
#'
#' @details
#' # Additional results
#' The `additional` argument is used to add results from additional analyses.
#' It must be a character vector containing one or more of:
#' \itemize{
#'  \item "lrt": Likehood ratio test comparing the full model to one without the
#'  exposure. Output will have additional columns "Chisq" (chi-square test
#'  statistics), "Df" (degrees of freedom), "p_Chisq" (p-value).
#'  \item "float": Standard errors for the exposure calculated by \link{float}
#'  or by \link[Epi]{ftrend} if the model has no intercept term.
#'  Output will have additional columns "float_estimate", "float_se" and
#'  "limits".
#'  \item "qvcalc": Standard errors for the exposure calculated by \link{qvcalc}.
#'  Output will have additional columns "qvcalc_estimate", "quasi_se" and
#'  "worstErrors".
#'  \item "" (default): Nothing.
#'  }
#' These will only be included when relevant e.g. standard errors from
#' \link{float} when the exposure has more than two levels.
#' This only has an effect if `return` is "coefs" or "exposure".
#'
#' See \link{float} and \link{qvcalc} for details about these results. Note that
#' "limits" and "worstErrors" do not correspond to particular groups but
#' accuracy other all contrasts.
#'
#' @param data A data frame.
#' @param exposure Name of column containing exposure variable.
#' @param outcome Name of the column containing binary outcome variable.
#' @param adjust_for Character vector of terms to include in the model formula.
#' @param exclude Names of columns that identify observations that are to be
#' excluded.
#' @param exclude_controls Names of columns that identify observations that are
#' to be exluded when outcome is equal to zero.
#' @param exposure_term The term(s) to include in the model for the exposure.
#' This can be a character vector or a function that take the exposure name as
#' the argument and returns a character vector. (Default: exposure)
#' @param additional Extra results to include.
#' See "Additional results" for details.
#' @param return
#' What should the function return:
#' \itemize{
#'  \item "model" (default): The model object (object of class "glm").
#'  \item "reduced-model": The model object with some large components removed.
#'  \item "coefs": A data frame of all estimated coefficients and other key
#'  model statistics.
#'  \item "exposure": A data frame of estimated coefficient(s) for the exposure
#'  and other key model statistics.
#' }
#' @param return_fn
#' If return="model", name of function to apply to model before it's returned.
#' (Default: identity)
#' @param ... Arguments passed to `stats::glm()`
#'
#' @returns Either a model object, reduced model object or data frame of
#' coefficients.
#'
#' @examples
#' logistic_model(ckbtools_participant_data,
#'                exposure = "bmi_calc",
#'                outcome = "has_diabetes",
#'                adjust_for = c("age_at_study_date", "sex"),
#'                return = "coefs")
#'
#' logistic_model(ckbtools_participant_data,
#'                exposure = "bmi_grp",
#'                outcome = "has_diabetes",
#'                adjust_for = c("age_at_study_date", "sex"),
#'                return = "coefs")
#'
#' @export
logistic_model <- function(
  data,
  exposure,
  outcome,
  adjust_for       = NULL,
  exclude          = NULL,
  exclude_controls = NULL,
  exposure_term    = exposure,
  additional       = "",
  return           = c("model", "reduced-model", "coefs", "exposure"),
  return_fn     = identity,
  ...
) {

  # check arguments
  return <- rlang::arg_match(return)
  additional <- rlang::arg_match(additional,
                                 values = c("lrt", "float", "qvcalc", "logLik", ""),
                                 multiple = TRUE)

  # exclude observations
  if (!is.null(exclude)) {
    data <- data %>%
      dplyr::filter(!dplyr::if_any(tidyselect::any_of(exclude), ~ .x %in% TRUE))  # treats NA as FALSE
  }

  # exclude observations from controls only
  if (!is.null(exclude_controls) && any(exclude_controls %in% names(data))) {
    data <- data %>%
      dplyr::filter(.data[[outcome]] | !dplyr::if_any(tidyselect::any_of(exclude_controls), ~ .x %in% TRUE))
  }


  # logistic regression model
  if (is.function(exposure_term)) {
    exposure_term <- exposure_term(exposure)
  }
  model_formula <- stats::reformulate(c(exposure_term, adjust_for),
                                      as.name(outcome))
  model <- stats::glm(formula = model_formula,
                      family = "binomial",
                      data = data,
                      ...)

  # return model
  if (return == "model") {
    return(return_fn(model))
  }

  # return reduced model
  if (return == "reduced-model") {
    return(reduce_model_object(model))
  }

  # get results from model
  coefs <- broom::tidy(model)
  model_frame <- stats::model.frame(model)
  exposure_terms <- find_exposure_terms(data, model, exposure, exposure_term)

  # Count observations per group
  if (!is.null(exposure) &&
      length(model$xlevels[[exposure]]) > 1) {
    counts <- table(model_frame[[exposure]])
    counts_cases <- table(model_frame[model_frame[[outcome]] == 1, ][[exposure]])
    counts_controls <- table(model_frame[model_frame[[outcome]] == 0, ][[exposure]])
    coefs <- dplyr::full_join(coefs,
                              data.frame(term = paste0(exposure, names(counts)),
                                         n_group = as.numeric(counts),
                                         ncases_group = as.numeric(counts_cases),
                                         ncontrols_group = as.numeric(counts_controls)),
                              by = "term")
  }

  # Epi::float standard errors
  # or Epi::ftrend standard errors if model has no intercept
  if ("float" %in% additional &&
      !is.null(exposure) &&
      !is.ordered(model$model[[exposure]]) &&
      length(model$xlevels[[exposure]]) > 2) {

    if (attr(stats::terms(model), "intercept") == 0) {
      far <- Epi::ftrend(model)
      flt <- tibble::tibble(
        term = paste0(exposure, names(far$coef)),
        float_estimate = far$coef,
        float_se = sqrt(diag(far$vcov)),
        row.names = NULL
      )
    } else {
      flt <- get_epi_float(model, exposure)
    }
    coefs <- dplyr::full_join(coefs, flt, by = "term")
  }

  # qvcalc::qvcalc standard errors
  if ("qvcalc" %in% additional &&
      !is.null(exposure) &&
      !is.ordered(model$model[[exposure]]) &&
      length(model$xlevels[[exposure]]) > 2) {
    coefs <- dplyr::full_join(coefs,
                              get_qvcalc(model, exposure, exposure_terms),
                              by = "term")
  }

  # organise coefs
  coefs <- organise_coefs(return, coefs, exposure_terms)

  # LRT
  if ("lrt" %in% additional) {
    null_model_formula <- stats::reformulate(adjust_for %||% "1",
                                             as.name(outcome))
    null_model <- stats::glm(null_model_formula,
                             family = "binomial",
                             data = data[!is.na(data[[exposure]]), ],
                             ...)
    anova_result <- stats::anova(null_model, model, test = "LRT")
    coefs[, c("Chisq", "Df", "p_Chisq")] <- NA
    coefs$Chisq[[1]] <- anova_result$Deviance[[2]]
    coefs$Df[[1]] <- anova_result$Df[[2]]
    coefs$p_Chisq[[1]] <- anova_result$`Pr(>Chi)`[[2]]
  }

  # Number of observations, 'events' and model formula
  coefs[, c("n", "ncases", "ncontrols", "formula")] <- NA
  coefs$n[[1]] <- nrow(model_frame)
  coefs$ncases[[1]] <- sum(model_frame[[outcome]] == 1)
  coefs$ncontrols[[1]] <- sum(model_frame[[outcome]] == 0)
  coefs <- dplyr::relocate(coefs,
                           tidyselect::any_of(c("n_group",
                                                "ncases_group",
                                                "ncontrols_group")),
                           .after = "ncontrols")
  coefs$formula[[1]] <- deparse1(model_formula)

  if ("lrt" %in% additional) {
    coefs[, "null_model_formula"] <- NA
    coefs$null_model_formula[[1]] <- deparse1(null_model_formula)
  }

  if ("logLik" %in% additional) {
    coefs[, c("logLik", "logLik_df")] <- NA
    coefs$logLik[[1]] <- stats::logLik(model)[[1]]
    coefs$logLik_df[[1]] <- attr(stats::logLik(model), "df")
  }

  # return coefficients
  return(coefs)
}

#' Estimate the association between an exposure and incident disease
#'
#' This function is used to estimate the (adjusted) association
#' between an exposure variable and an incident disease endpoint
#' using a Cox proportional hazards model.
#' It is a wrapper for the \link[survival]{coxph} function for convenience
#' and provides results from some additional analysis steps if requested.
#'
#' @details
#'
#' # Time scale
#' For analysis of incident disease using time-on-study as the time scale,
#' `outcome`, `outcome_date` and `start_date` are used to define follow up
#' time (right censored data).
#'
#' For other anlyses, the `time_in` and `time_out` arguments can be used.
#' If `time_in` and `time_out` are given, then these are starting and ending
#' time for interval survival data.
#' If `time_out` (but not `time_in`) is given, then this is the follow up time
#' (right censored data).
#'
#' # Additional results
#' The `additional` argument is used to add results from additional analyses.
#' It must be a character vector containing one or more of:
#' \itemize{
#'  \item "lrt": Likehood ratio test comparing the full model to one without the
#'  exposure. Output will have additional columns "Chisq" (chi-square test
#'  statistics), "Df" (degrees of freedom), "p_Chisq" (p-value).
#'  \item "float": Standard errors for the exposure calculated by \link{float}.
#'  Output will have additional columns "float_estimate", "float_se" and
#'  "limits".
#'  \item "qvcalc": Standard errors for the exposure calculated by \link{qvcalc}.
#'  Output will have additional columns "qvcalc_estimate", "quasi_se" and
#'  "worstErrors".
#'  \item "" (default): Nothing.
#'  }
#' These will only be included when relevant e.g. standard errors from
#' \link{float} when the exposure has more than two levels.
#' This only has an effect if `return` is "coefs" or "exposure".
#'
#' See \link{float} and \link{qvcalc} for details about these results. Note that
#' "limits" and "worstErrors" do not correspond to particular groups but
#' accuracy other all contrasts.
#'
#' @param data A data frame.
#' @param exposure Name of column containing exposure variable.
#' @param outcome Name of the column containing status indicator
#' (normally 0 = no disease, 1 = disease).
#' @param outcome_date Name of column containing date of event or censoring
#' date. (Default: `paste0(outcome, "_date)`)
#' @param time_in Name of column containing starting time for interval data.
#' @param time_out Name of column containing follow up time (or ending time for
#' interval data).
#' @param start_date Name of column containing date of start of follow up.
#' (Default: "study_date")
#' @param adjust_for Character vector of terms to include in the model formula.
#' @param exclude Names of columns that identify observations that should be
#' excluded.
#' @param exposure_term The term(s) to include in the model for the exposure.
#' This can be a character vector or a function that take the exposure name as
#' the argument and returns a character vector. (Default: exposure)
#' @param cluster Name of column that clusters the observations, for the
#' purposes of a robust variance.
#' @param init_values A named vector of initial values for model coefficients.
#' @param additional Extra results to include.
#' See "Additional results" for details.
#' @param return
#' What should the function return:
#' \itemize{
#'  \item "model" (default): The model object (see ?survival::coxph.object).
#'  \item "reduced-model": The model object with some large components removed.
#'  \item "coefs": A data frame of all estimated coefficients and other key
#'  model statistics.
#'  \item "exposure": A data frame of estimated coefficient(s) for the exposure
#'  and other key model statistics.
#' }
#' @param return_fn
#' If return="model", name of function to apply to model before it's returned.
#' (Default: identity)
#' @param ... Arguments passed to `survival::coxph()`
#'
#' @returns Either a model object, reduced model object or data frame of
#' coefficients.
#'
#' @examples
#' cox_model(ckbtools_participant_data,
#'           exposure = "bmi_calc",
#'           outcome = "all_cause_mortality",
#'           adjust_for = c("age_at_study_date", "strata(sex)"),
#'           return = "coefs")
#'
#' cox_model(ckbtools_participant_data,
#'           exposure = "bmi_grp",
#'           outcome = "all_cause_mortality",
#'           adjust_for = c("age_at_study_date", "strata(sex)"),
#'           return = "coefs")
#'
#' @export
cox_model <- function(
  data,
  exposure,
  outcome,
  outcome_date  = paste0(outcome, "_date"),
  start_date    = "study_date",
  time_in       = NULL,
  time_out      = NULL,
  adjust_for    = NULL,
  exclude       = NULL,
  exposure_term = exposure,
  cluster       = NULL,
  init_values   = c(0),
  additional    = "",
  return        = c("model", "reduced-model", "coefs", "exposure"),
  return_fn     = identity,
  ...
) {

  # check arguments
  return <- rlang::arg_match(return)
  additional <- rlang::arg_match(additional,
                                 values = c("lrt", "float", "qvcalc", "logLik", ""),
                                 multiple = TRUE)

  # exclude observations
  if (!is.null(exclude)) {
    data <- data %>%
      dplyr::filter(!dplyr::if_any(tidyselect::any_of(exclude), ~ .x %in% TRUE))  # treats NA as FALSE
  }

  # model formula
  if (is.function(exposure_term)) {
    exposure_term <- exposure_term(exposure)
  }
  if (!is.null(time_out)) {
    if (!is.null(time_in)) {
      model_response <- glue::glue("Surv({time_in}, {time_out}, {outcome})")
    } else {
      model_response <- glue::glue("Surv({time_out}, {outcome})")
    }
  } else {
    data$study_time <- data[[outcome_date]] - data[[start_date]]
    model_response <- paste0("Surv(study_time, ", outcome, ")")
  }
  model_formula <- stats::reformulate(
    c(exposure_term, adjust_for),
    model_response
  )

  # create init values
  draft_model <- survival::coxph(formula = model_formula,
                                 data = data,
                                 y = FALSE,
                                 iter.max = 0,
                                 cluster = if (is.null(cluster)) NULL else eval(as.name(cluster)),
                                 ...)
  init_values_aligned <- init_values[match(names(draft_model$coefficients),
                                           names(init_values))]
  init_values_aligned[is.na(init_values_aligned)] <- 0

  # Cox regression model
  model <- survival::coxph(model_formula,
                           data = data,
                           init = init_values_aligned,
                           y = FALSE,
                           iter.max = 40,
                           cluster = if (is.null(cluster)) NULL else eval(as.name(cluster)),
                           ...)

  # return model
  if (return == "model") {
    return(return_fn(model))
  }

  # return reduced model
  if (return == "reduced-model") {
    return(reduce_model_object(model))
  }

  # get results from model
  model_summary <- summary(model)
  model_frame <- stats::model.frame(model)
  if (is.null(stats::coef(model))) cli::cli_abort("No coefficients to return.")
  coefs <- broom::tidy(model)
  exposure_terms <- NULL
  if (!is.null(exposure)) {
    exposure_terms <- find_exposure_terms(data, model, exposure, exposure_term)
  }

  # Count observations per group
  if (!is.null(exposure) &&
      length(model$xlevels[[exposure]]) > 1) {
    counts <- table(model_frame[[exposure]])
    counts_events <- table(model_frame[model_frame[, 1][, "status"] == 1, ][[exposure]])
    coefs <- dplyr::full_join(coefs,
                              data.frame(term = paste0(exposure, names(counts)),
                                         n_group = as.numeric(counts),
                                         nevent_group = as.numeric(counts_events)),
                              by = "term")
  }

  # Epi::float standard errors
  if ("float" %in% additional &&
      !is.null(exposure) &&
      !is.ordered(model$model[[exposure]]) &&
      length(model$xlevels[[exposure]]) > 2) {
    coefs <- dplyr::full_join(coefs, get_epi_float(model, exposure), by = "term")
  }

  # qvcalc::qvcalc standard errors
  if ("qvcalc" %in% additional &&
      !is.null(exposure) &&
      !is.ordered(model$model[[exposure]]) &&
      length(model$xlevels[[exposure]]) > 2) {
    coefs <- dplyr::full_join(coefs,
                              get_qvcalc(model, exposure, exposure_terms),
                              by = "term")
  }

  # organise coefs
  coefs <- organise_coefs(return, coefs, exposure_terms)

  # LRT
  if ("lrt" %in% additional) {
    coefs[, c("Chisq", "Df", "p_Chisq")] <- NA

    if (!is.null(exposure)) {
      null_model_formula <- stats::reformulate(
        adjust_for %||% "1",
        model_response
      )
      draft_model <- survival::coxph(null_model_formula,
                                     data = data[!is.na(data[[exposure]]), ],
                                     y = FALSE,
                                     iter.max = 0,
                                     ...)
      init_values_aligned <- init_values[match(names(draft_model$coefficients),
                                               names(init_values))]
      init_values_aligned[is.na(init_values_aligned)] <- 0
      null_model <- survival::coxph(null_model_formula,
                                    data = data[!is.na(data[[exposure]]), ],
                                    init = init_values_aligned,
                                    y = FALSE,
                                    iter.max = 40,
                                    ...)
      anova_result <- stats::anova(model, null_model)
      coefs$Chisq[[1]] <- anova_result$Chisq[[2]]
      coefs$Df[[1]] <- anova_result$Df[[2]]
      coefs$p_Chisq[[1]] <- anova_result$`Pr(>|Chi|)`[[2]]
    }
  }

  # Number of observations, events and model formula
  coefs[, c("n", "nevent", "formula")] <- NA
  coefs$n[[1]] <- as.integer(model_summary$n)
  coefs$nevent[[1]] <- as.integer(model_summary$nevent)
  coefs <- dplyr::relocate(coefs,
                           tidyselect::any_of(c("n_group", "nevent_group")),
                           .after = "nevent")
  coefs$formula[[1]] <- deparse1(model_formula)

  if ("lrt" %in% additional) {
    coefs[, "null_model_formula"] <- NA
    coefs$null_model_formula[[1]] <- deparse1(null_model_formula)
  }

  if ("logLik" %in% additional) {
    coefs[, c("logLik", "logLik_df")] <- NA
    coefs$logLik[[1]] <- stats::logLik(model)[[1]]
    coefs$logLik_df[[1]] <- attr(stats::logLik(model), "df")
  }

  # return coefficients
  return(coefs)
}


#' Safe versions of regression model function
#'
#' These functions are wrappers around the regression model functions that
#' will capture side-effects such as errors and warnings and return a list.
#' This is particularly useful when using the model functions iteratively.
#'
#' @param ... Arguments are passed to the model function.
#'
#' @returns A named list:
#' \itemize{
#'  \item result: What is returned by the model function.
#'  \item errors: Errors produced by the model function as a character string.
#'  \item warnings: Warnings produced by the model function as a character
#'  string.
#'  \item messages: Messages produced by the model function as a character
#'  string.
#'  \item output: Output produced by the model function as a character string.
#' }
#'
#' @seealso [linear_model()], [logistic_model()], [cox_model()]
#'
#' @name safemodels
NULL


#' @rdname safemodels
#' @export
linear_model_safely <- function(...) {
  x <- purrr::quietly(purrr::safely(linear_model))(...)
  result <- x[["result"]][["result"]]
  error <- to_character_vector(x[["result"]][["error"]])
  warnings <- to_character_vector(x[["warnings"]])
  messages <- to_character_vector(x[["messages"]])
  output <- to_character_vector(x[["output"]])

  result_list <- tibble::lst(result, error, warnings, messages, output)
  return(result_list)
}


#' @rdname safemodels
#' @export
logistic_model_safely <- function(...) {
  x <- purrr::quietly(purrr::safely(logistic_model))(...)
  result <- x[["result"]][["result"]]
  error <- to_character_vector(x[["result"]][["error"]])
  warnings <- to_character_vector(x[["warnings"]])
  messages <- to_character_vector(x[["messages"]])
  output <- to_character_vector(x[["output"]])

  result_list <- tibble::lst(result, error, warnings, messages, output)
  return(result_list)
}


#' @rdname safemodels
#' @export
cox_model_safely <- function(...) {
  x <- purrr::quietly(purrr::safely(cox_model))(...)
  result <- x[["result"]][["result"]]
  error <- to_character_vector(x[["result"]][["error"]])
  warnings <- to_character_vector(x[["warnings"]])
  messages <- to_character_vector(x[["messages"]])
  output <- to_character_vector(x[["output"]])

  result_list <- tibble::lst(result, error, warnings, messages, output)
  return(result_list)
}

#' Multiple linear, logistic or Cox proportional hazards regression models
#'
#'
#' @param data A data frame.
#' @param model Type of model to use:
#' \itemize{
#'  \item "linear": Linear regression using `linear_model_safely()`.
#'  \item "logistic": Logistic regression using `logistic_model_safely()`.
#'  \item "cox": Cox proportional hazards model using `cox_model_safely()`.
#' }
#' @param exposures Names of columns containing exposure variables.
#' @param outcomes Names of columns containing outcome variables.
#' @param outcome_dates Names of columns containing outcome dates (for Cox model only).
#' @param start_date Name of column containing date of start of follow up.
#' (Default: "study_date")
#' @param time_in Names of columns containing starting time for interval data (for Cox model only).
#' @param time_out Names of columns containing follow up time (or ending time for
#' interval data) (for Cox model only).
#' @param age_at_risk A list of arguments for `age_at_risk()` to
#' expand data by age-at-risk for each outcome before fitting models for that
#' outcome. If this is used, "strata(agegrp)" is automatically added to each
#' element of `adjust_for`.
#' @param by Name of column (or multiple columns) to create subgroups.
#' @param exclude Names of columns that identify observations that should be
#' excluded.
#' @param exclude_by_exposure A named list of character vectors of column names.
#' When the name of an element matches the exposure, the column names will be
#' added to `exclude`.
#' @param exclude_by_outcome A named list of character vectors of column names.
#' When the name of an element matches the outcome, the column names will be
#' added to `exclude`.
#' @param exclude_controls Names of columns that identify observations that are
#' to be exluded when outcome is equal to zero.
#' @param exclude_controls_by_outcome A named list of character vectors of
#' column names. When the name of an element matches the outcome, the column
#' names will be added to `exclude_controls`.
#' @param adjust_for A named list of character vectors of terms to include in
#' model formulas.
#' @param remove_adjust_for A named list of character vectors of terms.
#' When the name of an element matches either the outcome or exposure, terms
#' given in the character vector will be removed from every element of
#' `adjust_for`.
#' @param add_adjust_for A named list of character vectors of terms.
#' When the name of an element matches either the outcome or exposure, terms
#' given in the character vector will be added to every element of `adjust_for`.
#' @param cluster Name of column that clusters the observations, for the
#' purposes of a robust variance.
#' @param init_values For Cox models, initial values for model coefficients.
#' This must be a data frame with columns "outcome", "subgroup", "adjustment",
#' "term" and "estimate", where term and estimate are the coefficient names and
#' initial values.
#' @param return
#' What should the model function return:
#' \itemize{
#'  \item "exposure" (default): A data frame of estimated coefficient(s) for the exposure
#'  and other key model statistics.
#'  \item "coefs": A data frame of all estimated coefficients and other key
#'  model statistics.
#'  \item "model": The model object.
#'  \item "reduced-model": The model object with some large components removed.
#' }
#' Note that R model objects can be very large (they can contain the model data,
#' residuals, fitted values, etc.) so it is not recommended to set
#' `return="model"` for a large number of models.
#' @param verbose "progress" (default) to show progress bars, "updates" for
#' messages without progress bars, or "quiet" to hide.
#' @param dryrun Return only the data frame of models to be run.
#' Columns (except 'subgroup', 'adjustment' and 'outcome_name') are the
#' arguments that will be used in the model function. (Default: FALSE)
#' @param ... Arguments passed to the model function.
#'
#' @returns A "tibble" data frame. If `return` is "model" or "reduced-model"
#' then it will contain a list-column called "model".
#'
#'
#' @seealso [linear_model()], [logistic_model()], [cox_model()],
#' [linear_model_safely()], [logistic_model_safely()], [cox_model_safely()]
#'
#'
#' @section Iteration:
#'
#' This function provides a convenient way to iterate over models using
#' [linear_model()], [logistic_model()] or [cox_model()]. First a data frame
#' defining all models is created. If the number of exposures is greater than or
#' equal to the number of outcomes (or Cox models with age-at-risk expansion are
#' used) then the models are split by outcome, otherwise they are split by
#' exposure. Within each set of models, [purrr::pmap()] is used to iterate over
#' the fitting of models. Results are then combined into a single data frame.
#'
#' For iteration of other models or analyses, you must write your own code. The
#' Iteration chapter of the first edition of “R for Data Science” provides an
#' introduction to iteration using R: <https://r4ds.had.co.nz/iteration.html>
#'
#'
#' @examples
#'
#' # linear regression models
#' multiple_models(
#'   ckbtools_participant_data,
#'   model = "linear",
#'   exposures = c("age_at_study_date", "sex", "bmi_calc", "bmi_grp"),
#'   outcomes = c("sbp_mean")
#' )
#'
#' # logistic regression models
#' multiple_models(
#'   ckbtools_participant_data,
#'   model = "logistic",
#'   exposures = c("age_at_study_date", "sex", "bmi_calc", "bmi_grp"),
#'   outcomes = c("has_diabetes")
#' )
#'
#' # Cox proportional hazards models
#' multiple_models(
#'   ckbtools_participant_data,
#'   model = "cox",
#'   exposures = c("age_at_study_date", "sex", "bmi_calc", "bmi_grp"),
#'   outcomes = c("all_cause_mortality", "ihd")
#' )
#'
#'
#' @export
multiple_models <- function(
  data,
  model,
  exposures,
  outcomes,
  outcome_dates               = paste0(outcomes, "_date"),
  start_date                  = "study_date",
  time_in                     = NULL,
  time_out                    = NULL,
  age_at_risk                 = NULL,
  by                          = NULL,
  exclude                     = NULL,
  exclude_controls            = NULL,
  exclude_by_exposure         = NULL,
  exclude_by_outcome          = NULL,
  exclude_controls_by_outcome = NULL,
  adjust_for                  = list(none = NULL),
  remove_adjust_for           = NULL,
  add_adjust_for              = NULL,
  cluster                     = NULL,
  init_values                 = NULL,
  return                      = c("coefs", "exposure", "reduced-model", "model"),
  verbose                     = if (rlang::is_interactive()) "progress" else "updates",
  dryrun                      = FALSE,
  ...
) {
  # check arguments
  verbose <- rlang::arg_match(verbose, c("progress", "updates", "quiet"))
  model <- rlang::arg_match(model, c("linear", "logistic", "cox"))
  return <- rlang::arg_match(return)
  dots <- rlang::list2(...)

  age_at_risk_used <- model == "cox" && !is.null(age_at_risk)
  if (age_at_risk_used) {
    adjust_for <- lapply(adjust_for, function(x) c(x, "strata(agegrp)"))
  }

  # keep only relevant columns in data
  data <- data %>%
    dplyr::select(dplyr::any_of(
      unique(c(exposures,
               outcomes,
               outcome_dates,
               start_date,
               time_in,
               time_out,
               unlist(age_at_risk[c("col.id", "col.dob", "col.keep")]),
               "csid", "dob_anon",  ## defaults from expand_age_at_risk()
               by,
               cluster,
               exclude,
               exclude_controls,
               unlist(c(exclude_by_exposure[exposures],
                        exclude_by_outcome[outcomes],
                        exclude_controls_by_outcome[outcomes])),
               if (length(unlist(adjust_for)) > 0) all.vars(stats::reformulate(unlist(adjust_for))),
               if (length(unlist(add_adjust_for)) > 0) all.vars(stats::reformulate(unlist(add_adjust_for)))))
    ))

  # split data
  if (!is.null(by)) {
    split_data <- split(data,
                        interaction(data[by], drop = FALSE, sep = "_"),
                        sep = "_")
    subgroup_names <- names(split_data)
  } else {
    split_data <- list(all = data)
    subgroup_names <- "all"
  }

  # remove original data
  rm(data)

  # create model function
  model_fn <- switch(
    model,
    linear   = linear_model_safely,
    logistic = logistic_model_safely,
    cox      = cox_model_safely,
    stop("Unsupported model: ", model)
  )

  # create data frame of all models
  all_models_grid <- tidyr::crossing(
    exposure   = exposures,
    outcome    = outcomes,
    subgroup   = subgroup_names,
    adjustment = names(adjust_for)
  ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(exclude = list(c(exclude,
                                   exclude_by_outcome[[.data$outcome]],
                                   exclude_by_exposure[[.data$exposure]])),
                  adjust_for = list(setdiff(c(adjust_for[[.data$adjustment]],
                                              add_adjust_for[[.data$outcome]],
                                              add_adjust_for[[.data$exposure]]),
                                            c(remove_adjust_for[[.data$outcome]],
                                              remove_adjust_for[[.data$exposure]])))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(outcome_name = .data$outcome, .before = "outcome") %>%
    dplyr::arrange(match(.data$exposure, .env$exposures),
                   match(.data$outcome, .env$outcomes),
                   match(.data$subgroup, .env$subgroup_names),
                   match(.data$adjustment, names(.env$adjust_for))) %>%
    dplyr::mutate(model_id = seq_len(dplyr::n()), .before = "exposure")

  # for logistic model, add exclude_control
  if (model == "logistic") {
    all_models_grid <- all_models_grid %>%
      dplyr::rowwise() %>%
      dplyr::mutate(exclude_controls = list(c(exclude_controls,
                                              exclude_controls_by_outcome[[.data$outcome]]))) %>%
      dplyr::ungroup()
  }

  # for cox model, add outcome_date, start_date, time_in, time_out, cluster
  if (model == "cox") {
    all_models_grid <- all_models_grid %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        outcome_date = outcome_dates[[match(.data$outcome, outcomes)]],
        start_date   = start_date,
        time_in      = time_in[[match(.data$outcome, outcomes)]],
        time_out     = time_out[[match(.data$outcome, outcomes)]],
        cluster      = cluster
      ) %>%
      dplyr::ungroup()

    ## if init_values are supplied, then add init_values column
    if (!is.null(init_values)) {
      init_values <- init_values %>%
        dplyr::select("outcome", "subgroup", "adjustment", "term", "estimate") %>%
        dplyr::group_by(.data$outcome, .data$subgroup, .data$adjustment) %>%
        dplyr::summarise(init_values = list(stats::setNames(.data$estimate, .data$term)))

      all_models_grid <- all_models_grid %>%
        dplyr::left_join(init_values, by = c("outcome", "subgroup", "adjustment")) %>%
        dplyr::mutate(init_values = dplyr::coalesce(init_values, list(0)))
    }
  }

  # prepare for age at risk expansion
  if (age_at_risk_used) {
    all_models_grid$outcome <- "status"
    all_models_grid$time_in <- "time_in"
    all_models_grid$time_out <- "time_out"
    if (length(unlist(c(adjust_for, exclude_by_exposure, exclude_by_outcome))) > 0) {
      all_adjust_for <- setdiff(all.vars(stats::reformulate(unlist(c(adjust_for, exclude_by_exposure, exclude_by_outcome)))), "agegrp")
    } else {
      all_adjust_for <- NULL
    }
  }

  ## if no age-at-risk expansion, data is ready for models
  data_for_models <- split_data

  # if dry run, return models data frame
  if (dryrun) {
    return(all_models_grid)
  }

  # determine if outcome or exposure is outer loop
  outcome_is_outer_loop <- (length(outcomes) <= length(exposures)) ||
    (model == "cox" && age_at_risk_used)

  if (outcome_is_outer_loop) {
    outer_iteration <- outcomes
  } else {
    outer_iteration <- exposures
  }

  results <- stats::setNames(vector("list", length(outer_iteration)), outer_iteration)

  # start outer loop
  for (outer_item in outer_iteration) {

    # create data frame of models in this iteration of the outer loop
    if (outcome_is_outer_loop) {
      models_grid <- all_models_grid %>%
        dplyr::filter(.data$outcome_name == .env$outer_item)
    } else {
      models_grid <- all_models_grid %>%
        dplyr::filter(.data$exposure == .env$outer_item)
    }

    # age-at-risk expansion
    if (age_at_risk_used && outcome_is_outer_loop) {

      ## remove existing data
      rm("data_for_models")

      if (identical(verbose, "updates")) {
        cli::cli_alert_success("Expanding data by age-at-risk for {outer_item}")
      }

      if (identical(verbose, "progress")) {
        cli::cli_progress_step("Expanding data by age-at-risk for {outer_item}")
      }
      data_for_models <- lapply(seq_along(split_data), function(x) {

        update_args <- list(col.entry_date = start_date,
                            col.exit_date = outcome_dates[[match(outer_item, outcomes)]],
                            col.status = outer_item)
        age_at_risk$col.keep <- unique(c(age_at_risk$col.keep, all_adjust_for, exposures)) ## size problem if many exposures?
        do.call(expand_age_at_risk,
                c(utils::modifyList(age_at_risk, update_args), data = list(split_data[[x]])))
      })
      names(data_for_models) <- names(split_data)
      if (identical(verbose, "progress")) {
        cli::cli_progress_done()
      }
    }
    ## end age-at-risk expansion

    # set up progress bar
    prog_bar_start <- glue::glue("Fitting {nrow(models_grid)} models for {outer_item} ")
    prog_bar <- list(format = paste0(prog_bar_start,
                                     "({cli::pb_current}/{cli::pb_total}) ",
                                     "{cli::pb_bar} {cli::pb_percent} | ",
                                     "ETA: {cli::pb_eta}"),
                     format_done = paste0("{cli::col_green(cli::symbol$tick)} ",
                                          prog_bar_start,
                                          "{cli::col_grey(sprintf('[%s]', cli::pb_elapsed))}"),
                     clear = FALSE,
                     show_after = 0)
    if (!identical(verbose, "progress")) {
      prog_bar <- FALSE
    }
    if (identical(verbose, "updates")) {
      cli::cli_alert_success(prog_bar_start)
    }

    # iterate over models
    results[[outer_item]] <- models_grid %>%
      dplyr::mutate(result =
                      purrr::pmap(
                        dplyr::select(models_grid, -"model_id", -"adjustment", -"outcome_name"),
                        function(..., subgroup) {
                          rlang::exec(model_fn, ..., data = data_for_models[[subgroup]], return = return, !!!dots)
                        },
                        .progress = prog_bar
                      )) %>%
      dplyr::select("model_id", outcome = "outcome_name", "exposure", "subgroup", "adjustment", "result")
  }

  # combine results from each loop, unnest results and organise
  if (return %in% c("reduced-model", "model")) {
    results <- dplyr::bind_rows(results) %>%
      dplyr::mutate(res = purrr::map(.data$result, build_nested_tibble_from_safely_output)) %>%
      dplyr::select(-"result") %>%
      tidyr::unnest("res") %>%
      dplyr::arrange(.data$model_id)
  } else {
    results <- dplyr::bind_rows(results) %>%
      dplyr::mutate(res = purrr::map(.data$result, build_tibble_from_safely_output)) %>%
      dplyr::select(-"result") %>%
      tidyr::unnest("res") %>%
      dplyr::arrange(.data$model_id)
  }

  return(results)
}
