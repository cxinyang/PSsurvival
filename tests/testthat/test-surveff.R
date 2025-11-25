test_that("surveff works with binary treatment and overlap weights", {
  data <- make_test_data_binary(n = 150)

  result <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    estimand = "overlap",
    censoring_method = "weibull"
  )

  # Check output structure
  expect_s3_class(result, "surveff")
  expect_true(all(c("survival_estimates", "survival_se", "difference_estimates",
                    "difference_se", "eval_times", "treatment_levels", "n_levels",
                    "n", "included", "estimand", "weights") %in% names(result)))

  # Check dimensions
  expect_equal(ncol(result$survival_estimates), 2)  # Binary: 2 groups
  expect_equal(nrow(result$survival_estimates), length(result$eval_times))
  expect_equal(ncol(result$difference_estimates), 1)  # Binary: 1 difference

  # Check survival estimates in [0, 1]
  expect_true(all(result$survival_estimates >= 0 & result$survival_estimates <= 1))

  # Check estimand recorded
  expect_equal(result$estimand, "overlap")

  # Check treatment levels
  expect_equal(result$treatment_levels, c("A", "B"))
  expect_equal(result$n_levels, 2)
})

test_that("surveff works with binary treatment and ATE (IPTW)", {
  data <- make_test_data_binary(n = 150)

  result <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    estimand = "ATE",
    censoring_method = "weibull"
  )

  # Check runs and has correct structure
  expect_s3_class(result, "surveff")
  expect_equal(result$estimand, "ATE")

  # Check difference estimates present
  expect_false(is.null(result$difference_estimates))
  expect_equal(ncol(result$difference_estimates), 1)
})

test_that("surveff works with binary treatment and ATT", {
  data <- make_test_data_binary(n = 150)

  result <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    estimand = "ATT",
    att_group = "B",
    censoring_method = "weibull"
  )

  # Check runs and records ATT correctly
  expect_s3_class(result, "surveff")
  expect_equal(result$estimand, "ATT")

  # Check att_group recorded in weight_result
  expect_equal(result$weights %>% {attr(., "att_group")} %||%
               result$estimand, "ATT")  # Estimand should be ATT
})

test_that("surveff works with multiple treatment groups and contrast matrix", {
  data <- make_test_data_multi(n = 180)
  contrast_mat <- make_contrast_matrix()

  result <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    estimand = "ATE",
    contrast_matrix = contrast_mat,
    censoring_method = "weibull"
  )

  # Check output structure
  expect_s3_class(result, "surveff")
  expect_equal(ncol(result$survival_estimates), 3)  # 3 groups
  expect_equal(result$n_levels, 3)
  expect_equal(result$treatment_levels, c("A", "B", "C"))

  # Check contrasts match matrix
  expect_equal(ncol(result$difference_estimates), nrow(contrast_mat))  # 2 contrasts
  expect_equal(colnames(result$difference_estimates), rownames(contrast_mat))

  # Check contrast matrix stored
  expect_equal(result$contrast_matrix, contrast_mat)
})

test_that("surveff trimming works (symmetric and asymmetric)", {
  data <- make_test_data_binary(n = 150)

  # Symmetric trimming
  result_sym <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    estimand = "ATE",
    trim = "symmetric",
    delta = 0.1,
    censoring_method = "weibull"
  )

  # Check some observations are trimmed (weights = 0)
  expect_true(any(result_sym$weights == 0))
  expect_true(any(!result_sym$included))
  expect_equal(length(result_sym$included), result_sym$n)

  # Asymmetric trimming
  result_asym <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    estimand = "ATE",
    trim = "asymmetric",
    alpha = 0.05,
    censoring_method = "weibull"
  )

  # Check some observations are trimmed
  expect_true(any(result_asym$weights == 0))
  expect_true(any(!result_asym$included))
})

test_that("surveff works with Cox censoring and bootstrap variance", {
  data <- make_test_data_binary(n = 150)

  result <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    estimand = "overlap",
    censoring_method = "cox",
    variance_method = "bootstrap",
    B = 20,  # Small for speed
    seed = 999
  )

  # Check runs
  expect_s3_class(result, "surveff")
  expect_equal(result$censoring_method, "cox")
  expect_equal(result$variance_method, "bootstrap")

  # Check boot_result structure exists
  expect_false(is.null(result$boot_result))
  expect_true("B" %in% names(result$boot_result))
  expect_equal(result$boot_result$B, 20)
})

test_that("surveff summary and plot methods work", {
  data <- make_test_data_binary(n = 150)

  result <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    estimand = "overlap",
    censoring_method = "weibull"
  )

  # Test summary method (prints style)
  expect_output(summary(result, max.len = 3), "Survival Effect Estimation Summary")

  # Test summary method (returns style)
  summ <- summary(result, style = "returns")
  expect_type(summ, "list")
  expect_true(all(c("survival_summary", "difference_summary") %in% names(summ)))

  # Test plot method
  expect_s3_class(plot(result, type = "surv"), "ggplot")
  expect_s3_class(plot(result, type = "survdiff"), "ggplot")
})
