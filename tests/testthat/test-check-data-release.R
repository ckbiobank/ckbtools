
test_that("check_data_release enforces option & version", {
  withr::local_options(list(cli.num_colors = 0))  # no ANSI colors

  withr::local_options(list(ckb.data.release = NULL))
  expect_error(check_data_release(), "Use `options\\(\\)` to set .*ckb\\.data\\.release")

  withr::local_options(list(ckb.data.release = "20.01"))
  expect_error(check_data_release(), "This function is not validated for CKB data release ")

  withr::local_options(list(ckb.data.release = "19.02"))
  expect_invisible(check_data_release())
})
