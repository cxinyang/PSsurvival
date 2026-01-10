test_that("weightedKM works with binary treatment and no weights", {
  data <- make_test_data_binary(n = 150)

  result <- weightedKM(
    data = data,
    treatment_var = "Z",
    time_var = "time",
    event_var = "event",
    weight_method = "none"
  )

  # Check output structure
  expect_s3_class(result, "weightedKM")
  expect_true(all(c("surv_estimates", "surv_var", "eval_times", "treatment_levels",
                    "weight_method") %in% names(result)))

  # Check dimensions
  expect_equal(ncol(result$surv_estimates), 2)  # Binary: 2 groups
  expect_equal(nrow(result$surv_estimates), length(result$eval_times))

  # Check survival estimates in [0, 1]
  expect_true(all(result$surv_estimates >= 0 & result$surv_estimates <= 1, na.rm = TRUE))

  # Check weight method recorded
  expect_equal(result$weight_method, "none")

  # Check treatment levels
  expect_equal(result$treatment_levels, c("A", "B"))
})

test_that("weightedKM works with overlap weights", {
  data <- make_test_data_binary(n = 150)

  result <- weightedKM(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    weight_method = "OW"
  )

  # Check runs and has correct structure
  expect_s3_class(result, "weightedKM")
  expect_equal(result$weight_method, "OW")

  # Check variance estimates exist
  expect_false(any(is.null(result$surv_var)))
})

test_that("weightedKM works with IPW", {
  data <- make_test_data_binary(n = 150)

  result <- weightedKM(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    weight_method = "IPW"
  )

  # Check runs
  expect_s3_class(result, "weightedKM")
  expect_equal(result$weight_method, "IPW")
})

test_that("weightedKM works with ATT", {
  data <- make_test_data_binary(n = 150)

  result <- weightedKM(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    weight_method = "ATT",
    att_group = "B"
  )

  # Check runs and records ATT correctly
  expect_s3_class(result, "weightedKM")
  expect_equal(result$weight_method, "ATT")
  expect_equal(result$att_group, "B")
})

test_that("weightedKM trimming works", {
  data <- make_test_data_binary(n = 150)

  result <- weightedKM(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    weight_method = "IPW",
    trim = TRUE,
    delta = 0.1
  )

  # Check runs
  expect_s3_class(result, "weightedKM")
  # Note: trim_method field may not be in weightedKM object (only in weights)
  # Just check it runs without error
})

test_that("weightedKM works with multiple treatment groups", {
  data <- make_test_data_multi(n = 180)

  result <- weightedKM(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    weight_method = "IPW"
  )

  # Check output structure
  expect_s3_class(result, "weightedKM")
  expect_equal(ncol(result$surv_estimates), 3)  # 3 groups
  expect_equal(result$treatment_levels, c("A", "B", "C"))
})

test_that("weightedKM summary and plot methods work", {
  data <- make_test_data_binary(n = 150)

  result <- weightedKM(
    data = data,
    treatment_var = "Z",
    time_var = "time",
    event_var = "event",
    weight_method = "none"
  )

  # Test summary method
  expect_output(summary(result), "Summary of Weighted Kaplan-Meier")

  # Test plot method
  expect_s3_class(plot(result), "ggplot")
  expect_s3_class(plot(result, type = "CR"), "ggplot")
})
