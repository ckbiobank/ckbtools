
test_that("linear_model with continuous exposure", {
  model <- lm(sbp_mean ~ bmi_calc + age_at_study_date_x100 + is_female,
              data = ckbtools_participant_data)
  ## coefficient for exposure is moved first
  estimates <- as.numeric(coef(model)[c(2, 1, 3, 4)])
  std.errors <- as.numeric(sqrt(diag(vcov(model)))[c(2, 1, 3, 4)])
  p.values <- as.numeric(summary(model)$coefficients[c(2, 1, 3, 4), "Pr(>|t|)"])
  expected_result <- tibble::tibble(
    term      = c("bmi_calc", "(Intercept)", "age_at_study_date_x100", "is_female"),
    estimate  = estimates,
    std.error = std.errors,
    statistic = estimates / std.errors,
    p.value   = p.values,
    n         = c(nrow(ckbtools_participant_data), NA, NA, NA),
    formula   = c("sbp_mean ~ bmi_calc + age_at_study_date_x100 + is_female", NA, NA, NA)
  )

  actual_result <- linear_model(
    data       = ckbtools_participant_data,
    exposure   = "bmi_calc",
    outcome    = "sbp_mean",
    adjust_for = c("age_at_study_date_x100", "is_female"),
    return     = "coefs"
  )

  expect_equal(actual_result, expected_result)
})

test_that("linear_model with categorical exposure", {
  model <- lm(sbp_mean ~ smoking + age_at_study_date_x100 + is_female,
              data = ckbtools_participant_data)
  model_frame <- stats::model.frame(model)
  estimates <- c(NA_real_, as.numeric(coef(model)[2:4]))
  std.errors <- c(NA_real_, as.numeric(sqrt(diag(vcov(model)))[c(2:4)]))
  p.values <- c(NA_real_, as.numeric(summary(model)$coefficients[2:4, "Pr(>|t|)"]))

  # float
  flt <- Epi::float(model)

  # emmeans marginal means
  marginal_means <- emmeans::emmeans(model, specs = "smoking", weights = "proportional")
  summary_marginal_means <- summary(marginal_means)
  qvcalc_marginal_means <- qvcalc::qvcalc(stats::vcov(marginal_means))

  # qvcalc
  qvcalc_result <- qvcalc::qvcalc(model, "smoking")

  # ANOVA (full model vs model without exposure)
  null_model <- lm(sbp_mean ~ age_at_study_date_x100 + is_female,
                   data = ckbtools_participant_data[!is.na(ckbtools_participant_data[["smoking"]]), ])
  anova_result <- stats::anova(null_model, model)

  # expected result
  expected_result <- tibble::tibble(
    term               = paste0("smoking", levels(ckbtools_participant_data$smoking)),
    estimate           = estimates,
    std.error          = std.errors,
    statistic          = estimates / std.errors,
    p.value            = p.values,
    float_estimate     = as.numeric(flt$coef),
    float_se           = sqrt(flt$var),
    limits             = c(flt$limits, NA, NA),
    qvcalc_estimate    = as.numeric(qvcalc_result$qvframe$estimate[, 1]),
    quasi_se           = qvcalc_result$qvframe$quasiSE,
    worstErrors        = c(qvcalc::worstErrors(qvcalc_result), NA, NA),
    emmean             = summary_marginal_means$emmean,
    emmean_se          = summary_marginal_means$SE,
    emmean_quasi_se    = as.numeric(qvcalc_marginal_means$qvframe$quasiSE),
    F                  = c(anova_result$F[[2]], NA, NA, NA),
    Res.Df             = c(anova_result$Res.Df[[2]], NA, NA, NA),
    Df                 = c(anova_result$Df[[2]], NA, NA, NA),
    p_F                = c(anova_result$`Pr(>F)`[[2]], NA, NA, NA),
    n                  = c(nrow(model_frame), NA, NA, NA),
    n_group            = as.numeric(table(model_frame[["smoking"]])),
    formula            = c("sbp_mean ~ smoking + age_at_study_date_x100 + is_female", NA, NA, NA),
    null_model_formula = c("sbp_mean ~ age_at_study_date_x100 + is_female", NA, NA, NA),
  )

  actual_result <- linear_model(
    data       = ckbtools_participant_data,
    exposure   = "smoking",
    outcome    = "sbp_mean",
    adjust_for = c("age_at_study_date_x100", "is_female"),
    additional = c("float", "anova", "emmeans", "qvcalc"),
    return     = "exposure"
  )

  expect_equal(actual_result, expected_result)
})


test_that("linear_model with exclusion", {
  participant_data_mod <- ckbtools_participant_data
  participant_data_mod$exclude_row <- FALSE
  participant_data_mod$exclude_row[1:10] <- TRUE
  participant_data_mod$exclude_row[11] <- NA   # NA should be treated as FALSE by linear_model()

  model <- lm(sbp_mean ~ bmi_calc + age_at_study_date_x100 + is_female,
              data = participant_data_mod[!participant_data_mod$exclude_row %in% TRUE, ])
  n_obs <- nrow(model.frame(model))

  result <- linear_model(
    data       = participant_data_mod,
    exposure   = "bmi_calc",
    outcome    = "sbp_mean",
    adjust_for = c("age_at_study_date_x100", "is_female"),
    exclude    = "exclude_row",
    return     = "exposure"
  )

  expect_equal(c(result$estimate[[1]], result$n[[1]]),
               c(coef(model)[[2]], n_obs))
})

test_that("linear_model with squared exposure term", {
  model <- lm(sbp_mean ~ bmi_calc + I(bmi_calc^2) +
                age_at_study_date_x100 + is_female,
              data = ckbtools_participant_data)
  estimates <- as.numeric(coef(model)[c(2, 3, 1, 4, 5)])
  std.errors <- as.numeric(sqrt(diag(vcov(model)))[c(2, 3, 1, 4, 5)])
  p.values <- as.numeric(summary(model)$coefficients[c(2, 3, 1, 4, 5), "Pr(>|t|)"])
  expected_result <- tibble::tibble(
    term      = c("bmi_calc", "I(bmi_calc^2)",
                  "(Intercept)", "age_at_study_date_x100", "is_female"),
    estimate  = estimates,
    std.error = std.errors,
    statistic = estimates / std.errors,
    p.value   = p.values,
    n         = NA_integer_,
    formula   = NA_character_
  )
  expected_result$n[[1]] <- nrow(ckbtools_participant_data)
  expected_result$formula[[1]] <- "sbp_mean ~ bmi_calc + I(bmi_calc^2) + age_at_study_date_x100 + is_female"

  actual_result <- linear_model(
    data          = ckbtools_participant_data,
    exposure      = "bmi_calc",
    exposure_term = c("bmi_calc", "I(bmi_calc^2)"),
    outcome       = "sbp_mean",
    adjust_for    = c("age_at_study_date_x100", "is_female"),
    return        = "coefs"
  )

  expect_equal(actual_result, expected_result)
})


test_that("logistic_model with continuous exposure", {

  model <- glm(
    has_diabetes ~ bmi_calc + age_at_study_date + sex + sbp_mean,
    family = "binomial",
    data = ckbtools_participant_data
  )
  ## coefficient for exposure is moved first
  estimates <- as.numeric(coef(model)[c(2, 1, 3, 4, 5)])
  std.errors <- as.numeric(sqrt(diag(vcov(model)))[c(2, 1, 3, 4, 5)])

  expected_result <- tibble::tibble(
    term      = c("bmi_calc", "(Intercept)",
                  "age_at_study_date", "sexFemale", "sbp_mean"),
    estimate  = estimates,
    std.error = std.errors,
    statistic = estimates / std.errors,
    p.value   = 2 * (1 - pnorm(abs(estimates / std.errors))),
    n         = NA_integer_,
    ncases    = NA_integer_,
    ncontrols = NA_integer_,
    formula   = NA_character_
  )

  model_frame <- stats::model.frame(model)
  expected_result$n[[1]] <- nrow(model_frame)
  expected_result$ncases[[1]] <- sum(model_frame$has_diabetes == 1)
  expected_result$ncontrols[[1]] <- sum(model_frame$has_diabetes == 0)
  expected_result$formula[[1]] <- "has_diabetes ~ bmi_calc + age_at_study_date + sex + sbp_mean"

  actual_result <- logistic_model(
    data       = ckbtools_participant_data,
    exposure   = "bmi_calc",
    outcome    = "has_diabetes",
    adjust_for = c("age_at_study_date", "sex", "sbp_mean"),
    return     = "coefs"
  )

  expect_equal(actual_result, expected_result)
})


test_that("logistic_model with categorical exposure", {

  model <- glm(
    has_diabetes ~ smoking + age_at_study_date + sex,
    family = binomial,
    data = ckbtools_participant_data
  )
  model_frame <- stats::model.frame(model)
  ## coefficient for exposure is moved first
  estimates <- c(NA_real_, as.numeric(coef(model)[c(2, 3, 4, 1, 5, 6)]))
  std.errors <- c(NA_real_, as.numeric(sqrt(diag(vcov(model)))[c(2, 3, 4, 1, 5, 6)]))
  p.values <- c(NA_real_, as.numeric(summary(model)$coefficients[c(2, 3, 4, 1, 5, 6), "Pr(>|z|)"]))

  null_model <- glm(
    has_diabetes ~ age_at_study_date + sex,
    family = binomial,
    data = ckbtools_participant_data
  )
  a <- anova(null_model, model, test = "LRT")

  expected_result <- tibble::tibble(
    term      = c(paste0("smoking", levels(ckbtools_participant_data$smoking)),
                  "(Intercept)", "age_at_study_date", "sexFemale"),
    estimate  = estimates,
    std.error = std.errors,
    statistic = estimates / std.errors,
    p.value   = 2 * (1 - pnorm(abs(estimates / std.errors))),
    Chisq     = NA_real_,
    Df        = NA_real_,
    p_Chisq   = NA_real_,
    n         = NA_integer_,
    ncases    = NA_integer_,
    ncontrols = NA_integer_,
    n_group         = NA_integer_,
    ncases_group    = NA_integer_,
    ncontrols_group = NA_integer_,
    formula   = NA_character_,
    null_model_formula = NA_character_
  )

  expected_result$Chisq[[1]] <- a$Deviance[[2]]
  expected_result$Df[[1]] <- a$Df[[2]]
  expected_result$p_Chisq[[1]] <- a$`Pr(>Chi)`[[2]]
  expected_result$n[[1]] <- nrow(model_frame)
  expected_result$ncases[[1]] <- sum(model_frame$has_diabetes == 1, na.rm = TRUE)
  expected_result$ncontrols[[1]] <- sum(model_frame$has_diabetes == 0, na.rm = TRUE)
  expected_result$n_group[1:4]   <- as.numeric(table(model_frame$smoking))
  expected_result$ncases_group[1:4] <- as.numeric(table(model_frame$smoking[model_frame$has_diabetes == 1]))
  expected_result$ncontrols_group[1:4] <- as.numeric(table(model_frame$smoking[model_frame$has_diabetes == 0]))
  expected_result$formula[[1]]   <- "has_diabetes ~ smoking + age_at_study_date + sex"
  expected_result$null_model_formula[[1]]   <- "has_diabetes ~ age_at_study_date + sex"

  actual_result <- logistic_model(
    data       = ckbtools_participant_data,
    exposure   = "smoking",
    outcome    = "has_diabetes",
    adjust_for = c("age_at_study_date", "sex"),
    additional = "lrt",
    return     = "coefs"
  )

  expect_equal(actual_result, expected_result)
})


test_that("logistic_model with exclusion", {
  participant_data_mod <- ckbtools_participant_data
  participant_data_mod$exclude_row <- FALSE
  participant_data_mod$exclude_row[1:100] <- TRUE
  participant_data_mod$exclude_row[101:120] <- NA

  ### exclude in controls only
  model <- glm(
    has_diabetes ~ bmi_calc + age_at_study_date + sex,
    family = "binomial",
    data = participant_data_mod[is.na(participant_data_mod$exclude_row) |
                                  !participant_data_mod$exclude_row, ]
  )
  ## coefficient for exposure is moved first
  estimates <- as.numeric(coef(model)[c(2, 1, 3, 4)])
  std.errors <- as.numeric(sqrt(diag(vcov(model)))[c(2, 1, 3, 4)])

  expected_result <- tibble::tibble(
    term      = c("bmi_calc", "(Intercept)",
                  "age_at_study_date", "sexFemale"),
    estimate  = estimates,
    std.error = std.errors,
    statistic = estimates / std.errors,
    p.value   = 2 * (1 - pnorm(abs(estimates / std.errors))),
    n         = NA_integer_,
    ncases    = NA_integer_,
    ncontrols = NA_integer_,
    formula   = NA_character_
  )

  model_frame <- stats::model.frame(model)
  expected_result$n[[1]] <- nrow(model_frame)
  expected_result$ncases[[1]] <- sum(model_frame$has_diabetes == 1)
  expected_result$ncontrols[[1]] <- sum(model_frame$has_diabetes == 0)
  expected_result$formula[[1]] <- "has_diabetes ~ bmi_calc + age_at_study_date + sex"

  actual_result <- logistic_model(
    data       = participant_data_mod,
    exposure   = "bmi_calc",
    outcome    = "has_diabetes",
    adjust_for = c("age_at_study_date", "sex"),
    exclude    = "exclude_row",
    return     = "coefs"
  )

  expect_equal(actual_result, expected_result)


  ### exclude in controls only
  model <- glm(
    has_diabetes ~ bmi_calc + age_at_study_date + sex,
    family = "binomial",
    data = participant_data_mod[participant_data_mod$has_diabetes == 1 |
                                  is.na(participant_data_mod$exclude_row) |
                                  !participant_data_mod$exclude_row, ]
  )
  ## coefficient for exposure is moved first
  estimates <- as.numeric(coef(model)[c(2, 1, 3, 4)])
  std.errors <- as.numeric(sqrt(diag(vcov(model)))[c(2, 1, 3, 4)])

  expected_result <- tibble::tibble(
    term      = c("bmi_calc", "(Intercept)",
                  "age_at_study_date", "sexFemale"),
    estimate  = estimates,
    std.error = std.errors,
    statistic = estimates / std.errors,
    p.value   = 2 * (1 - pnorm(abs(estimates / std.errors))),
    n         = NA_integer_,
    ncases    = NA_integer_,
    ncontrols = NA_integer_,
    formula   = NA_character_
  )

  model_frame <- stats::model.frame(model)
  expected_result$n[[1]] <- nrow(model_frame)
  expected_result$ncases[[1]] <- sum(model_frame$has_diabetes == 1)
  expected_result$ncontrols[[1]] <- sum(model_frame$has_diabetes == 0)
  expected_result$formula[[1]] <- "has_diabetes ~ bmi_calc + age_at_study_date + sex"

  actual_result <- logistic_model(
    data       = participant_data_mod,
    exposure   = "bmi_calc",
    outcome    = "has_diabetes",
    adjust_for = c("age_at_study_date", "sex"),
    exclude_controls = "exclude_row",
    return     = "coefs"
  )

  expect_equal(actual_result, expected_result)
})


test_that("cox_model with continuous exposure", {

  ckbtools_participant_data$study_time <- ckbtools_participant_data$all_cause_mortality_date - ckbtools_participant_data$study_date


  model <- survival::coxph(
    Surv(study_time, all_cause_mortality) ~ bmi_calc + age_at_study_date + sex + strata(region),
    data = ckbtools_participant_data
  )

  estimates <- as.numeric(coef(model))
  std.errors <- as.numeric(sqrt(diag(vcov(model))))

  expected_result <- tibble::tibble(
    term      = c("bmi_calc", "age_at_study_date", "sexFemale"),
    estimate  = estimates,
    std.error = std.errors,
    statistic = estimates / std.errors,
    p.value   = 2 * (1 - pnorm(abs(estimates / std.errors))),
    n         = NA_integer_,
    nevent    = NA_integer_,
    formula   = NA_character_
  )

  model_frame <- stats::model.frame(model)
  expected_result$n[[1]] <- nrow(model_frame)
  expected_result$nevent[[1]] <- sum(model_frame[, 1][, 2])
  expected_result$formula[[1]] <- "Surv(study_time, all_cause_mortality) ~ bmi_calc + age_at_study_date + sex + strata(region)"

  actual_result <- cox_model(
    data       = ckbtools_participant_data,
    exposure   = "bmi_calc",
    outcome    = "all_cause_mortality",
    adjust_for = c("age_at_study_date", "sex", "strata(region)"),
    return     = "coefs"
  )

  expect_equal(actual_result, expected_result)
})


test_that("cox_model with categorical exposure", {

  ckbtools_participant_data$study_time <- ckbtools_participant_data$ihd_date - ckbtools_participant_data$study_date


  model <- survival::coxph(
    Surv(study_time, ihd) ~ smoking + age_at_study_date + sex + strata(region),
    data = ckbtools_participant_data
  )

  null_model <- survival::coxph(
    Surv(study_time, ihd) ~ age_at_study_date + sex + strata(region),
    data = ckbtools_participant_data
  )
  lrt <- stats::anova(null_model, model)

  estimates <- c(NA_real_, as.numeric(coef(model)))
  std.errors <- c(NA_real_, as.numeric(sqrt(diag(vcov(model)))))

  expected_result <- tibble::tibble(
    term      = c(paste0("smoking", levels(ckbtools_participant_data$smoking)),
                  "age_at_study_date", "sexFemale"),
    estimate  = estimates,
    std.error = std.errors,
    statistic = estimates / std.errors,
    p.value   = 2 * (1 - pnorm(abs(estimates / std.errors))),
    Chisq     = NA_real_,
    Df        = NA_integer_,
    p_Chisq   = NA_real_,
    n         = NA_integer_,
    nevent    = NA_integer_,
    n_group   = NA_integer_,
    nevent_group = NA_integer_,
    formula   = NA_character_,
    null_model_formula = NA_character_
  )

  model_frame <- stats::model.frame(model)
  expected_result$n[[1]] <- nrow(model_frame)
  expected_result$nevent[[1]] <- sum(model_frame[, 1][, 2])
  expected_result$n_group[1:4]   <- as.numeric(table(model_frame$smoking))
  expected_result$nevent_group[1:4] <- as.numeric(table(model_frame$smoking[model_frame[, 1][, 2] == 1]))
  expected_result$Chisq[[1]] <- lrt$Chisq[[2]]
  expected_result$Df[[1]] <- lrt$Df[[2]]
  expected_result$p_Chisq[[1]] <- lrt$`Pr(>|Chi|)`[[2]]
  expected_result$formula[[1]] <- "Surv(study_time, ihd) ~ smoking + age_at_study_date + sex + strata(region)"
  expected_result$null_model_formula[[1]] <- "Surv(study_time, ihd) ~ age_at_study_date + sex + strata(region)"

  actual_result <- cox_model(
    data       = ckbtools_participant_data,
    exposure   = "smoking",
    outcome    = "ihd",
    adjust_for = c("age_at_study_date", "sex", "strata(region)"),
    return     = "coefs",
    additional = "lrt"
  )

  expect_equal(actual_result, expected_result)
})

test_that("return reduced model", {

  actual_result <- linear_model(
    data       = ckbtools_participant_data,
    exposure   = "bmi_calc",
    outcome    = "sbp_mean",
    adjust_for = c("age_at_study_date_x100", "is_female"),
    return     = "reduced-model"
  )
  expect_null(actual_result$data)
  expect_null(actual_result$model)
  expect_null(actual_result$y)
  expect_null(actual_result$residuals)
  expect_null(attr(actual_result$terms, ".Environment"))
  expect_null(attr(actual_result$formula, ".Environment"))
})


test_that("capture error in multiple_models", {

  actual_result <- multiple_models(
    data      = ckbtools_participant_data,
    model     = "linear",
    exposures = "bmi_clac", # mis-spelled
    outcomes  = "sbp_mean",
    verbose   = "quiet"
  )
  expect_equal(actual_result$error[[1]],
               "Error in eval(predvars, data, env): object 'bmi_clac' not found\n")
})
