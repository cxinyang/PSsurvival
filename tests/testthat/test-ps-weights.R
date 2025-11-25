# Tests for propensity score estimation and weight calculation
# These tests call internal functions directly

test_that("PS estimates are in [0, 1] for binary treatment", {
  data <- make_test_data_binary(n = 150)

  ps_result <- estimate_ps(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    ps_control = list()
  )

  # Check all PS values in [0, 1]
  expect_true(all(ps_result$ps_matrix >= 0))
  expect_true(all(ps_result$ps_matrix <= 1))

  # Check observed PS in [0, 1]
  expect_true(all(ps_result$ps >= 0))
  expect_true(all(ps_result$ps <= 1))
})

test_that("PS rows sum to 1 for multiple treatment groups", {
  data <- make_test_data_multi(n = 180)

  ps_result <- estimate_ps(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    ps_control = list()
  )

  # Check rows sum to 1 (ignore names from rowSums)
  row_sums <- rowSums(ps_result$ps_matrix)
  expect_equal(as.numeric(row_sums), rep(1, nrow(data)), tolerance = 1e-6)
})

test_that("Binary treatment uses glm for PS model", {
  data <- make_test_data_binary(n = 150)

  ps_result <- estimate_ps(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    ps_control = list()
  )

  # Check model object is glm
  expect_s3_class(ps_result$ps_model, "glm")
})

test_that("Multiple treatment uses multinom for PS model", {
  data <- make_test_data_multi(n = 180)

  ps_result <- estimate_ps(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    ps_control = list()
  )

  # Check model object is multinom (from nnet package)
  expect_s3_class(ps_result$ps_model, "multinom")
})

test_that("Overlap weights are in [0, 1]", {
  data <- make_test_data_binary(n = 150)

  # First estimate PS
  ps_result <- estimate_ps(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    ps_control = list()
  )

  # Then estimate overlap weights
  weight_result <- estimate_weights(
    ps_result = ps_result,
    data = data,
    treatment_var = "Z",
    estimand = "overlap",
    att_group = NULL,
    trim = NULL,
    delta = NULL,
    alpha = NULL
  )

  # Check all weights in [0, 1]
  expect_true(all(weight_result$weights >= 0))
  expect_true(all(weight_result$weights <= 1))
})

test_that("ATT weights: treated group has weight = 1", {
  data <- make_test_data_binary(n = 150)

  # First estimate PS
  ps_result <- estimate_ps(
    data = data,
    treatment_var = "Z",
    ps_formula = Z ~ X1 + X2 + B1,
    ps_control = list()
  )

  # Estimate ATT weights targeting group B
  weight_result <- estimate_weights(
    ps_result = ps_result,
    data = data,
    treatment_var = "Z",
    estimand = "ATT",
    att_group = "B",
    trim = NULL,
    delta = NULL,
    alpha = NULL
  )

  # Check weights for treated group (B) are all 1
  is_treated <- data$Z == "B"
  expect_equal(weight_result$weights[is_treated], rep(1, sum(is_treated)), tolerance = 1e-6)

  # Check att_group recorded
  expect_equal(weight_result$att_group, "B")
})
