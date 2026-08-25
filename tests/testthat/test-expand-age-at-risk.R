test_that("expand_age_at_risk: no expansion", {
  output <- expand_age_at_risk(ckbtools_participant_data,
                               col.exit_date = "ihd_date",
                               col.status   = "ihd")

  expect_equal(output$csid, ckbtools_participant_data$csid)
  expect_equal(output$dob_anon, ckbtools_participant_data$dob_anon)
  expect_equal(output$study_date, ckbtools_participant_data$study_date)
  expect_equal(output$start_int, ckbtools_participant_data$study_date)
  expect_equal(output$end_int, ckbtools_participant_data$ihd_date)
  expect_equal(output$status, ckbtools_participant_data$ihd)
  expect_equal(output$t_start_days, rep(0, nrow(ckbtools_participant_data)))
  expect_equal(output$t_end_days,
               as.numeric(ckbtools_participant_data$ihd_date - ckbtools_participant_data$study_date) + 0.95)
})


test_that("expand_age_at_risk: error if columns are not present", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob = as.Date(c("1960-01-01", "1970-01-01")),
    start_date = as.Date(c("2000-01-01", "2005-01-01")),
    end_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    expand_age_at_risk(input_data),
    "dob_anon, study_date, and endpoint_date do not exist in data."
  )

  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob = as.Date(c("1960-01-01", "1970-01-01")),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    expand_age_at_risk(input_data),
    "dob_anon does not exist in data."
  )
})


test_that("expand_age_at_risk: error if date columns are not valid", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = c("1960-01-01", "1970-01-01"),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    expand_age_at_risk(input_data),
    "dob.anon, study_date and endpoint_date must be dates."
  )
})


test_that("expand_age_at_risk: warning for non-integer dates", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("1960-01-01", "1970-01-01")) + 0.5,
    study_date = as.Date(c("2000-01-01", "2005-01-01")) + 0.25,
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")) + 0.75,
    endpoint = c(1, 0)
  )

  expect_warning(
    expand_age_at_risk(input_data),
    "Data contains non-integer dates."
  )
})


test_that("expand_age_at_risk: no warning for integer dates only", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("1960-01-01", "1970-01-01")),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_silent(
    expand_age_at_risk(input_data)
  )
})


test_that("expand_age_at_risk: warning if status column contains missing values", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("1960-01-01", "1970-01-01")),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, NA)  # Includes missing value
  )

  expect_warning(
    expand_age_at_risk(input_data),
    "endpoint contains missing values, these observations will be dropped."
  )
})

test_that("expand_age_at_risk: error if id, dob, entry_date or exit_date column contains missing values", {
  input_data <- tibble::tibble(
    csid          = c(NA, 2),
    dob_anon      = as.Date(c("1960-01-01", "1970-01-01")),
    study_date    = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint      = c(1, 0)
  )
  expect_error(
    expand_age_at_risk(input_data),
    "csid contains missing values."
  )

  input_data <- tibble::tibble(
    csid          = c(1, 2),
    dob_anon      = as.Date(c(NA, "1970-01-01")),
    study_date    = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c(NA, "2015-01-01")),
    endpoint      = c(1, 0)
  )
  expect_error(
    expand_age_at_risk(input_data),
    "dob_anon and endpoint_date contain missing values."
  )

  input_data <- tibble::tibble(
    csid          = c(1, 2),
    dob_anon      = as.Date(c("1960-01-01", "1970-01-01")),
    study_date    = as.Date(c(NA, "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint      = c(1, 0)
  )
  expect_error(
    expand_age_at_risk(input_data),
    "study_date contains missing values."
  )

  input_data <- tibble::tibble(
    csid          = c(1, 2),
    dob_anon      = as.Date(c("1960-01-01", "1970-01-01")),
    study_date    = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c(NA, "2015-01-01")),
    endpoint      = c(1, 0)
  )
  expect_error(
    expand_age_at_risk(input_data),
    "endpoint_date contains missing values."
  )
})


test_that("expand_age_at_risk: ages must be a numeric ", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("1960-01-01", "1970-01-01")),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )
  expect_error(
    expand_age_at_risk(input_data, ages = c(50, 30, 60)),
    "`ages` must be a numeric vector with distinct, increasing values."
  )
  expect_error(
    expand_age_at_risk(input_data, ages = c("thirty", "forty", "fifty")),
    "`ages` must be a numeric vector with distinct, increasing values."
  )
})


test_that("expand_age_at_risk: error if DOB variable is later than entry date variable", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("2000-01-01", "1970-01-01")),
    study_date = as.Date(c("1960-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    expand_age_at_risk(input_data),
    "study_date must not be before dob_anon."
  )
})



test_that("expand_age_at_risk: error if entry date variable is later than exit date variable", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("1960-01-01", "1970-01-01")),
    study_date = as.Date(c("2010-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2000-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  expect_error(
    expand_age_at_risk(input_data),
    "endpoint_date must not be before study_date."
  )
})




test_that("expand_age_at_risk: example output using package data", {
  participant_data_long <- expand_age_at_risk(data          = ckbtools_participant_data,
                                              col.exit_date = "ihd_date",
                                              col.status    = "ihd",
                                              ages          = seq(from = 40, to = 80, by = 5))

  temp_file <- tempfile(fileext = ".csv")
  write.csv(participant_data_long[1:100, ], file = temp_file)
  expect_snapshot_file(temp_file, "expand_age_at_risk1.csv")

  participant_data_long <- expand_age_at_risk(data          = ckbtools_participant_data,
                                              col.exit_date = "all_cause_mortality_date",
                                              col.status    = "all_cause_mortality",
                                              ages          = seq(from = 40, to = 80, by = 5))

  temp_file <- tempfile(fileext = ".csv")
  write.csv(participant_data_long[1:100, ], file = temp_file)
  expect_snapshot_file(temp_file, "expand_age_at_risk2.csv")
})



test_that("expand_age_at_risk: example output", {
  participant_data <- data.frame(csid = 1,
                                 dob_anon = as.Date("1940-02-15"),
                                 study_date = as.Date("2005-06-10"),
                                 endpoint_date = as.Date("2010-11-02"),
                                 endpoint = 1)
  output <- expand_age_at_risk(participant_data,
                               ages = seq(30, 80, by = 5))

  expect_equal(nrow(output), 2)
  expect_equal(output$agegrp, c("65-70", "70-75"))
  expect_equal(output$status, c(0, 1))
  expect_equal(output$t_start_days, c(0, 1711))
  expect_equal(output$t_end_days, c(1710.95, 1971.95))
})



test_that("expand_age_at_risk: example output", {
  participant_data <- data.frame(
    csid = 1,
    dob_anon = as.Date("1975-06-15"),
    study_date = as.Date("2005-01-01"),
    endpoint_date = as.Date("2017-12-31"),
    endpoint = 0
  )

  output <- expand_age_at_risk(participant_data,
                               ages = seq(30, 80, by = 5))
  expect_equal(nrow(output), 3)
  expect_equal(output$agegrp, c("30-35", "35-40", "40-45"))
  expect_equal(output$status, c(0, 0, 0))
  expect_equal(output$t_start_days, c(165, 1991, 3817))
  expect_equal(output$t_end_days, c(1990.95, 3816.95, 4747.95))
})




test_that("expand_age_at_risk: output date columns are valid Date type and output endpoint is valid numeric type", {
  input_data <- tibble::tibble(
    csid = c(1, 2),
    dob_anon = as.Date(c("1960-01-01", "1970-01-01")),
    study_date = as.Date(c("2000-01-01", "2005-01-01")),
    endpoint_date = as.Date(c("2010-01-01", "2015-01-01")),
    endpoint = c(1, 0)
  )

  output <- expand_age_at_risk(input_data, ages = seq(30, 80, by = 5))

  expect_true(all(sapply(output[, c("dob_anon", "study_date", "endpoint_date")], lubridate::is.Date)),
              expect_true(is.numeric, output$endpoint))
})
