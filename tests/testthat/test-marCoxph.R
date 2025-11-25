test_that("marCoxph works with binary treatment and overlap weights", {
  data <- make_test_data_binary(n = 150)

  result <- marCoxph(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    reference_level = "A",
    estimand = "overlap",
    variance_method = "robust"
  )

  # Check output structure
  expect_s3_class(result, "marCoxph")
  expect_true(all(c("coxph_fitted", "logHR_est", "logHR_se_robust",
                    "treatment_var", "treatment_levels", "reference_level",
                    "n_levels", "estimand") %in% names(result)))

  # Check logHR estimates
  expect_length(result$logHR_est, 1)  # Binary: 1 HR (B vs A)
  expect_named(result$logHR_est, "Z:B")

  # Check SE present
  expect_length(result$logHR_se_robust, 1)
  expect_true(result$logHR_se_robust > 0)

  # Check estimand recorded
  expect_equal(result$estimand, "overlap")

  # Check reference level
  expect_equal(result$reference_level, "A")
})

test_that("marCoxph works with binary treatment and ATE (IPTW)", {
  data <- make_test_data_binary(n = 150)

  result <- marCoxph(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    reference_level = "A",
    estimand = "ATE",
    variance_method = "robust"
  )

  # Check runs and has correct structure
  expect_s3_class(result, "marCoxph")
  expect_equal(result$estimand, "ATE")
})

test_that("marCoxph works with binary treatment and ATT", {
  data <- make_test_data_binary(n = 150)

  result <- marCoxph(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    reference_level = "A",
    estimand = "ATT",
    att_group = "B",
    variance_method = "robust"
  )

  # Check runs and records ATT correctly
  expect_s3_class(result, "marCoxph")
  expect_equal(result$estimand, "ATT")
  expect_equal(result$att_group, "B")
})

test_that("marCoxph works with multiple treatment groups", {
  data <- make_test_data_multi(n = 180)

  result <- marCoxph(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    reference_level = "A",
    estimand = "ATE",
    variance_method = "robust"
  )

  # Check output structure
  expect_s3_class(result, "marCoxph")
  expect_equal(result$n_levels, 3)
  expect_equal(result$treatment_levels, c("A", "B", "C"))

  # Check number of hazard ratios = n_levels - 1
  expect_length(result$logHR_est, 2)  # B vs A, C vs A

  # Check names formatted correctly
  expect_named(result$logHR_est, c("Z:B", "Z:C"))
})

test_that("marCoxph trimming works (symmetric and asymmetric)", {
  data <- make_test_data_binary(n = 150)

  # Symmetric trimming
  result_sym <- marCoxph(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    reference_level = "A",
    estimand = "ATE",
    trim = "symmetric",
    delta = 0.1,
    variance_method = "robust"
  )

  # Check runs
  expect_s3_class(result_sym, "marCoxph")

  # Check sample size reduced due to trimming
  expect_true(sum(result_sym$n_coxph_fitted) < 150)

  # Asymmetric trimming
  result_asym <- marCoxph(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    reference_level = "A",
    estimand = "ATE",
    trim = "asymmetric",
    alpha = 0.05,
    variance_method = "robust"
  )

  # Check runs and sample size reduced
  expect_s3_class(result_asym, "marCoxph")
  expect_true(sum(result_asym$n_coxph_fitted) < 150)
})

test_that("marCoxph bootstrap variance works with full and strata sampling", {
  data <- make_test_data_binary(n = 150)

  # Bootstrap with full sampling
  result_full <- marCoxph(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    reference_level = "A",
    estimand = "overlap",
    variance_method = "bootstrap",
    boot_level = "full",
    B = 20,  # Small for speed
    seed = 999
  )

  # Check bootstrap SE present
  expect_false(is.null(result_full$logHR_se_bootstrap))
  expect_length(result_full$logHR_se_bootstrap, 1)
  expect_true(result_full$logHR_se_bootstrap > 0)

  # Check boot_result structure
  expect_false(is.null(result_full$boot_result))
  expect_equal(result_full$boot_result$B, 20)
  expect_equal(result_full$variance_method, "bootstrap-full")

  # Bootstrap with strata sampling
  result_strata <- marCoxph(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    reference_level = "A",
    estimand = "overlap",
    variance_method = "bootstrap",
    boot_level = "strata",
    B = 20,
    seed = 888
  )

  # Check bootstrap SE present
  expect_false(is.null(result_strata$logHR_se_bootstrap))
  expect_equal(result_strata$variance_method, "bootstrap-strata")
})

test_that("marCoxph summary method works", {
  data <- make_test_data_binary(n = 150)

  result <- marCoxph(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    time_var = "time",
    event_var = "event",
    reference_level = "A",
    estimand = "overlap",
    variance_method = "robust"
  )

  # Test summary method (prints style)
  expect_output(summary(result), "Marginal Cox Model Summary")

  # Test summary method (returns style)
  summ <- summary(result, style = "returns")
  expect_type(summ, "list")
  expect_true(all(c("logHR", "HR", "SE", "logHR_CI_lower", "logHR_CI_upper",
                    "HR_CI_lower", "HR_CI_upper") %in% names(summ)))

  # Check exponentiated HR
  expect_equal(summ$HR, exp(summ$logHR), tolerance = 1e-6)
})
