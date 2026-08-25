
test_that("summarise_continuous", {
  result <- summarise_continuous(
    data = ckbtools_participant_data,
    cols = c("sbp_mean", "bmi_calc", "waist_mm", "met")
  )

  expect_equal(
    result$mean,
    unname(apply(ckbtools_participant_data[, c("sbp_mean", "bmi_calc", "waist_mm", "met")],
                 2,
                 mean,
                 na.rm = TRUE))
  )

  expect_equal(
    result$sd,
    unname(apply(ckbtools_participant_data[, c("sbp_mean", "bmi_calc", "waist_mm", "met")],
                 2,
                 sd,
                 na.rm = TRUE))
  )

  expect_equal(
    result$min,
    unname(apply(ckbtools_participant_data[, c("sbp_mean", "bmi_calc", "waist_mm", "met")],
                 2,
                 min,
                 na.rm = TRUE))
  )

  expect_equal(
    result$max,
    unname(apply(ckbtools_participant_data[, c("sbp_mean", "bmi_calc", "waist_mm", "met")],
                 2,
                 max,
                 na.rm = TRUE))
  )

  expect_equal(
    result$p50,
    unname(apply(ckbtools_participant_data[, c("sbp_mean", "bmi_calc", "waist_mm", "met")],
                 2,
                 median,
                 na.rm = TRUE))
  )
})
