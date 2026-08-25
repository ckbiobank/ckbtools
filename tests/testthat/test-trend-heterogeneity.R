

test_that("trend: correct output for increasing trend", {
  result <- trend(c(1.2, 1.4, 1.8), c(0.2, 0.3, 0.4))

  # Manually calculate expected values
  w <- 1 / c(0.2, 0.3, 0.4)^2
  z <- seq_along(c(1.2, 1.4, 1.8))
  z_centered <- z - (sum(z * w) / sum(w))
  numerator <- sum(w * c(1.2, 1.4, 1.8) * z_centered)
  denominator <- sum(w * z_centered^2)
  expected_test_statistic <- (numerator^2) / denominator
  expected_p <- pchisq(expected_test_statistic, df = 1, lower.tail = FALSE)

  expect_equal(round(result$`test statistic`, 6), round(expected_test_statistic, 6))
  expect_equal(round(result$p, 6), round(expected_p, 6))
})


test_that("trend: correct output for decreasing trend", {
  result <- trend(c(1.8, 1.4, 1.2), c(0.2, 0.3, 0.4))

  w <- 1 / c(0.2, 0.3, 0.4)^2
  z <- seq_along(c(1.8, 1.4, 1.2))
  z_centered <- z - (sum(z * w) / sum(w))
  numerator <- sum(w * c(1.8, 1.4, 1.2) * z_centered)
  denominator <- sum(w * z_centered^2)
  expected_test_statistic <- (numerator^2) / denominator
  expected_p <- pchisq(expected_test_statistic, df = 1, lower.tail = FALSE)

  expect_equal(round(result$`test statistic`, 6), round(expected_test_statistic, 6))
  expect_equal(round(result$p, 6), round(expected_p, 6))
})


test_that("trend: correct output for flat trend", {
  result <- trend(c(1, 1, 1), c(0.2, 0.3, 0.4))

  w <- 1 / c(0.2, 0.3, 0.4)^2
  z <- seq_along(c(1, 1, 1))
  z_centered <- z - (sum(z * w) / sum(w))
  numerator <- sum(w * c(1, 1, 1) * z_centered)
  denominator <- sum(w * z_centered^2)
  expected_test_statistic <- (numerator^2) / denominator
  expected_p <- pchisq(expected_test_statistic, df = 1, lower.tail = FALSE)

  expect_equal(round(result$`test statistic`, 6), round(expected_test_statistic, 6))
  expect_equal(round(result$p, 6), round(expected_p, 6))
})


test_that("heterogeneity: correct result for perfectly homogeneous data", {
  result <- heterogeneity(c(1.0, 1.0, 1.0), c(0.2, 0.2, 0.2))
  expect_equal(result$`test statistic`, 0)
  expect_equal(result$`degrees of freedom`, 2)
  expect_equal(result$p, 1)
})


test_that("trend/heterogeneity: error for mismatched lengths", {
  expect_error(trend(c(1.2, 1.4, 1.8), c(0.2, 0.3)),
               "`beta` and `se` must have same length.")
  expect_error(heterogeneity(c(1.2, 1.4, 1.8), c(0.2, 0.3)),
               "`beta` and `se` must have same length.")
})


test_that("trend/heterogeneity: error for fewer than two points", {
  expect_error(trend(c(1.2), c(0.2)),
               "At least two estimates are required.")
  expect_error(heterogeneity(c(1.2), c(0.2)),
               "At least two estimates are required.")
})


test_that("trend/heterogeneity: error for non-numeric values", {
  expect_error(
    trend(c("a", "b", "c"), c(0.2, 0.3, 0.4)),
    "`beta` and `se` must be numeric vectors."
  )
  expect_error(
    heterogeneity(c(1.2, "a", 1.8), c(0.2, 0.3, 0.4)),
    "`beta` and `se` must be numeric vectors."
  )
  expect_error(
    heterogeneity(c(1.2, 1.4, 1.8), c(0.2, "b", 0.4)),
    "`beta` and `se` must be numeric vectors."
  )
})


test_that("trend/heterogeneity: error for missing values in beta", {
  expect_error(trend(c(1.2, NA, 1.8), c(0.2, 0.3, 0.4)),
               "`beta` and `se` must not contain missing values.")
  expect_error(heterogeneity(c(1.2, NA, 1.8), c(0.2, 0.3, 0.4)),
               "`beta` and `se` must not contain missing values.")
})


test_that("trend: error for zero or negative values in se", {
  expect_error(trend(c(1.2, 1.4, 1.8), c(0.2, 0, 0.4)),
               "Standard errors must be greater than zero")
  expect_error(trend(c(1.2, 1.4, 1.8), c(0.2, -0.3, 0.4)),
               "Standard errors must be greater than zero")
  expect_error(heterogeneity(c(1.2, 1.4, 1.8), c(0.2, 0, 0.4)),
               "Standard errors must be greater than zero")
  expect_error(heterogeneity(c(1.2, 1.4, 1.8), c(0.2, -0.3, 0.4)),
               "Standard errors must be greater than zero")
})
