
test_that("multiple linear models dry run", {

  exposures <- c("sex", "age_at_baseline_group", "education_4_groups", "smoking", "has_diabetes")
  outcomes <- c("sbp_mean", "bmi_calc")

  result <- multiple_models(
    data = ckbtools_participant_data,
    model = "linear",
    exposures = exposures,
    outcomes = outcomes,
    verbose = "quiet",
    return = "model",
    dryrun = TRUE
  )

  expected_result <- tidyr::crossing(
    exposure = exposures,
    outcome_name = outcomes
  ) %>%
    dplyr::arrange(match(.data$exposure, .env$exposures),
                   match(.data$outcome_name, .env$outcomes)) %>%
    dplyr::mutate(model_id = seq_len(dplyr::n()), .before = "exposure") %>%
    dplyr::mutate(outcome = outcome_name,
                  subgroup = "all",
                  adjustment = "none",
                  exclude = rep(list(NULL), dplyr::n()),
                  adjust_for = rep(list(NULL), dplyr::n()))

  expect_identical(result, expected_result)

})

test_that("multiple linear models", {
  result <- multiple_models(
    data = ckbtools_participant_data,
    model = "linear",
    exposures = "sex",
    outcomes = c("sbp_mean", "bmi_calc"),
    verbose = "quiet",
    return = "model"
  ) %>%
    dplyr::mutate(model = purrr:::map(model, broom::tidy))

  expected_result <- tibble::tibble(
    model_id = 1:2,
    outcome = c("sbp_mean", "bmi_calc"),
    exposure = "sex",
    subgroup = "all",
    adjustment = "none",
    error = NA_character_,
    warnings = NA_character_,
    messages = NA_character_,
    output = NA_character_,
    model = list(
      broom::tidy(lm(sbp_mean ~ sex, data = ckbtools_participant_data)),
      broom::tidy(lm(bmi_calc ~ sex, data = ckbtools_participant_data))
      )
  )

  expect_identical(result, expected_result)

})

test_that("multiple linear models: snapshot", {
  result <- multiple_models(
    data = ckbtools_participant_data,
    model = "linear",
    exposures = c("sex", "age_at_baseline_group", "education_4_groups", "smoking", "has_diabetes"),
    outcomes = c("sbp_mean", "bmi_calc"),
    verbose = "quiet"
  )
  temp_file <- tempfile(fileext = ".csv")
  write.csv(result, file = temp_file)
  expect_snapshot_file(temp_file, "multiple_models1.csv")
})

test_that("multiple logistic models: snapshot", {
  result <- multiple_models(
    data = ckbtools_participant_data,
    model = "logistic",
    exposures = c("age_at_baseline_group", "education_4_groups", "smoking", "has_diabetes"),
    outcomes = c("ihd", "all_cause_mortality"),
    by = "sex",
    verbose = "quiet"
  )
  temp_file <- tempfile(fileext = ".csv")
  write.csv(result, file = temp_file)
  expect_snapshot_file(temp_file, "multiple_models2.csv")
})

test_that("multiple cox models: snapshot", {
  result <- multiple_models(
    data = ckbtools_participant_data,
    model = "cox",
    exposures = c("age_at_baseline_group", "education_4_groups", "smoking", "has_diabetes"),
    outcomes = c("ihd", "all_cause_mortality"),
    verbose = "quiet"
  )
  temp_file <- tempfile(fileext = ".csv")
  write.csv(result, file = temp_file)
  expect_snapshot_file(temp_file, "multiple_models3.csv")
})

test_that("multiple cox models - age at risk: snapshot", {
  result <- multiple_models(
    data = ckbtools_participant_data,
    model = "cox",
    exposures = c("education_4_groups", "smoking", "has_diabetes"),
    outcomes = c("ihd", "all_cause_mortality"),
    age_at_risk = list(ages = seq(40, 90, 5)),
    verbose = "quiet"
  )
  temp_file <- tempfile(fileext = ".csv")
  write.csv(result, file = temp_file)
  expect_snapshot_file(temp_file, "multiple_models4.csv")
})
