


test_that("incidence_rates: error if columns called age and calendar in data", {
  ckbtools_participant_data$age <- runif(40, 80, 5000)
  ckbtools_participant_data$calendar <- 0
  ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc,
                                           breaks = c(0, 25, 30, Inf),
                                           right  = FALSE,
                                           labels = c("<25", "25-29.9", "30+"))
  expect_error(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ep0002_date",
                    col.status     = "ep0002"),
    "Column names calendar and age cannot be used in `data`."
  )
})

test_that("incidence_rates: error if missing column", {
  participant_data <- tibble::tibble(
    csid = c(1, 2),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    incidence_rates(participant_data),
    "dob_anon does not exist in data."
  )
})

test_that("incidence_rates: error if date columns are not valid", {
  participant_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = c("1960-01-01", "1970-01-01"),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    incidence_rates(participant_data),
    "dob.anon, study_date and endpoint_date must be dates."
  )
})

test_that("incidence_rates: error if date columns are not valid", {
  participant_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = c("1960-01-01", "1970-01-01"),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    incidence_rates(participant_data),
    "dob.anon, study_date and endpoint_date must be dates."
  )
})

test_that("incidence_rates: error if date columns are not valid", {
  participant_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("1960-01-01", "1970-01-01")),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c("dead", "alive")
  )

  expect_error(
    incidence_rates(participant_data),
    "endpoint must be numeric."
  )
})




test_that("incidence_rates: errors/warnings for group and adjust_for columns", {
  ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc,
                                           breaks = c(0, 25, 30, Inf),
                                           right  = FALSE,
                                           labels = c("<25", "25-29.9", "30+"))

  expect_error(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    group         = "bmi_calc"),
    "bmi_calc is numeric and has 236 levels."
  )
  expect_error(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    group      = "bmi_grp",
                    adjust_for = "bmi_calc"),
    "bmi_calc is numeric and has 236 levels."
  )

  expect_error(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    group         = "ihd_date"),
    "ihd_date is a date column."
  )

  expect_error(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    group          = "bmi_grp",
                    adjust_for     = "ihd_date"),
    "ihd_date is a date column."
  )

  expect_warning(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    group          = "is_female"),
    "is_female is numeric and has 2 levels, so treating as a group variable."
  )
  expect_warning(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    group      = NULL,
                    adjust_for = "is_female"),
    "is_female is numeric and has 2 levels, so treating as a group variable."
  )

})



test_that("incidence_rates: error if any missing data", {
  ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc,
                                           breaks = c(0, 25, 30, Inf),
                                           right  = FALSE)
  ckbtools_participant_data$waist_grp <- cut(ckbtools_participant_data$waist_mm,
                                             breaks = c(0, 800, 900),
                                             right  = FALSE,
                                             labels = c("A", "B"))

  expect_error(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    group         = "waist_grp"),
    "waist_grp contains missing values."
  )

  expect_error(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    group          = "bmi_grp",
                    adjust_for = "waist_grp"),
    "waist_grp contains missing values."
  )
})




test_that("incidence_rates: ages must be a numeric ", {
  ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc,
                                           breaks = c(0, 25, 30, Inf),
                                           right  = FALSE,
                                           labels = c("<25", "25-29.9", "30+"))

  expect_error(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    ages           = c(50, 30, 60)),
    "`ages` must a numeric vector with distinct, increasing values."
  )

  expect_error(
    incidence_rates(ckbtools_participant_data,
                    col.dob        = "dob_anon",
                    col.entry_date = "study_date",
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    ages           = c("thirty", "forty", "fifty"),
                    group         = "bmi_grp"),
    "`ages` must a numeric vector with distinct, increasing values."
  )
})


test_that("incidence_rates: error if entry date variable is later than exit date variable", {
  participant_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("1960-01-01", "1970-01-01")),
    study_date = as.Date(c("2010-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2000-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    incidence_rates(participant_data),
    "endpoint_date must not be before study_date."
  )
})


test_that("incidence_rates: error if DOB variable is later than entry date variable", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("2000-01-01", "1970-01-01")),
    study_date = as.Date(c("1960-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    incidence_rates(input_data),
    "study_date must not be before dob_anon."
  )
})


test_that("incidence_rates: correct output if no standardisation", {
  ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc,
                                           breaks = c(0, 25, 30, Inf),
                                           right  = FALSE,
                                           labels = c("<25", "25-29.9", "30+"))

  output <- incidence_rates(ckbtools_participant_data,
                            col.exit_date  = "ihd_date",
                            col.status     = "ihd",
                            group         = "bmi_grp")

  expect_equal(output$incidence_rates_in_all_strata$incidencerate,
               output$standardised_incidence_rates$crude_ir / 1000)
  expect_equal(output$standardised_incidence_rates$crude_ir,
               output$standardised_incidence_rates$std_ir)
})




test_that("incidence_rates: consistent totals", {
  ckbtools_participant_data$met_grp <- cut(ckbtools_participant_data$met,
                                           breaks = c(0, 10, 20, 30, Inf),
                                           right  = FALSE)
  ckbtools_participant_data$region_is_urban <- factor(ckbtools_participant_data$region_is_urban)

  output <- incidence_rates(ckbtools_participant_data,
                            col.exit_date = "ihd_date",
                            col.status    = "ihd",
                            ages          = seq(40, 80, 10),
                            group         = "met_grp",
                            adjust_for    = "region_is_urban",
                            method        = "norm")

  expect_equal(sum(output$incidence_rates_in_all_strata$pyar),
               sum(output$standardised_incidence_rates$crude_pyar))
  expect_equal(sum(output$incidence_rates_in_all_strata$nevents),
               sum(output$standardised_incidence_rates$crude_nevents))
})


test_that("incidence_rates: warning if no time at risk in strata", {
  ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc,
                                           breaks = c(0, 25, 30, Inf),
                                           right  = FALSE,
                                           labels = c("<25", "25-29.9", "30+"))
  ckbtools_participant_data$is_female <- factor(ckbtools_participant_data$is_female)
  ckbtools_participant_data <- ckbtools_participant_data %>%
    dplyr::filter(!(is_female == "0" &
                      region == "Gansu" &
                      bmi_grp == "30+"))

  expect_warning(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date = "ihd_date",
                    col.status    = "ihd",
                    group         = "bmi_grp",
                    adjust_for    = c("is_female", "region")),
    "There are strata with no time at risk. Returning only incidence rates within strata."
  )
})


test_that("incidence_rates: warning if strata with zero or < 20 events", {
  ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc,
                                           breaks = c(0, 25, 30, Inf),
                                           right  = FALSE,
                                           labels = c("<25", "25-29.9", "30+"))
  ckbtools_participant_data$is_female <- factor(ckbtools_participant_data$is_female)

  expect_warning(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "ihd_date",
                    col.status     = "ihd",
                    group          = "bmi_grp",
                    adjust_for     = "region"),
    "There are strata with fewer than 20 events."
  )

  expect_warning(
    incidence_rates(ckbtools_participant_data,
                    col.exit_date  = "respiratory_diseases_date",
                    col.status     = "respiratory_diseases",
                    group          = "bmi_grp",
                    adjust_for     = c("is_female", "region")),
    "There are strata with zero events."
  )
})


test_that("incidence rates: example output", {
  ckbtools_participant_data$bmi_grp <- cut(ckbtools_participant_data$bmi_calc,
                                           breaks = c(0, 25, 30, Inf),
                                           right  = FALSE,
                                           labels = c("<25", "25-29.9", "30+"))

  suppressWarnings(
    output <- incidence_rates(ckbtools_participant_data,
                              col.exit_date  = "ihd_date",
                              col.status     = "ihd",
                              group          = "bmi_grp",
                              ages            = seq(40, 80, 10),
                              adjust_for     = "is_female")
  )

  suppressWarnings(
    output_gamma <- incidence_rates(ckbtools_participant_data,
                                    col.exit_date  = "ihd_date",
                                    col.status     = "ihd",
                                    group          = "bmi_grp",
                                    ages            = seq(40, 80, 10),
                                    adjust_for     = "is_female",
                                    method         = "gamma")
  )

  suppressWarnings(
    output_external_weights <- incidence_rates(
      ckbtools_participant_data,
      col.exit_date  = "ihd_date",
      col.status     = "ihd",
      group          = "bmi_grp",
      ages            = seq(40, 80, 10),
      adjust_for     = "is_female",
      weights        = data.frame(is_female = c(0, 1),
                                  std_pyar = c(1, 1))
    )
  )

  expect_equal(output$standardised_incidence_rates$std_ir,
               output_gamma$standardised_incidence_rates$std_ir)

  temp_file <- tempfile(fileext = ".csv")
  write.csv(output, file = temp_file)
  expect_snapshot_file(temp_file, "incidence_rate_output.csv")

  temp_file <- tempfile(fileext = ".csv")
  write.csv(output_gamma, file = temp_file)
  expect_snapshot_file(temp_file, "incidence_rate_output_gamma.csv")

  temp_file <- tempfile(fileext = ".csv")
  write.csv(output_external_weights, file = temp_file)
  expect_snapshot_file(temp_file, "incidence_rate_output_external_weights.csv")
})
