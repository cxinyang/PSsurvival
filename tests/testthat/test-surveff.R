test_that("surveff works with binary treatment and overlap weights", {
  data <- make_test_data_binary(n = 150)

  result <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    weight_method = "OW",
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

  # Check estimand recorded (internal mapping: OW -> overlap)
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
    weight_method = "IPW",
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
    weight_method = "ATT",
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
    weight_method = "IPW",
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

test_that("surveff trimming works (symmetric)", {
  data <- make_test_data_binary(n = 150)

  # Symmetric trimming
  result_sym <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    weight_method = "IPW",
    trim = TRUE,
    delta = 0.1,
    censoring_method = "weibull"
  )

  # Check some observations are trimmed (weights = 0)
  expect_true(any(result_sym$weights == 0))
  expect_true(any(!result_sym$included))
  expect_equal(length(result_sym$included), result_sym$n)
})

test_that("surveff works with Cox censoring and bootstrap variance", {
  data <- make_test_data_binary(n = 150)

  result <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    weight_method = "OW",
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
    weight_method = "OW",
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

test_that("tidy.surveff and df.residual.surveff return broom-style output", {
  data <- make_test_data_binary(n = 150)

  result <- surveff(
    data = data,
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    eval_times = c(0.5, 1, 1.5),
    weight_method = "OW",
    censoring_method = "weibull"
  )

  td <- tidy.surveff(result)
  expect_s3_class(td, "data.frame")
  expect_true(all(c("term", "parameter", "estimate", "std.error",
                    "statistic", "p.value") %in% names(td)))

  # 2 groups + 1 contrast, each at 3 eval_times = 9 rows
  n_groups <- ncol(result$survival_estimates)
  n_contr <- ncol(result$difference_estimates)
  n_times <- length(result$eval_times)
  expect_equal(nrow(td), (n_groups + n_contr) * n_times)

  # Terms are unique (required for pool() to group correctly across imputations)
  expect_equal(length(unique(td$term)), nrow(td))

  # Survival-row estimates match the source matrix (natural probability scale)
  surv_rows <- td[td$parameter == "survival", ]
  expect_equal(surv_rows$estimate,
               as.numeric(result$survival_estimates))
  expect_equal(surv_rows$std.error,
               as.numeric(result$survival_se))

  # Difference-row estimates match the source matrix
  diff_rows <- td[td$parameter == "difference", ]
  expect_equal(diff_rows$estimate, as.numeric(result$difference_estimates))

  # conf.int columns
  td_ci <- tidy.surveff(result, conf.int = TRUE)
  expect_true(all(c("conf.low", "conf.high") %in% names(td_ci)))

  # tidy() must absorb the extra args mice injects (effects, parametric, dfcom)
  expect_equal(
    tidy.surveff(result, effects = "fixed", parametric = TRUE, dfcom = Inf),
    td
  )

  # df.residual is Inf (large-sample pooling for curve estimators)
  expect_identical(df.residual.surveff(result), Inf)

  # Dispatch via the generics/stats generics (as used by mice::pool())
  expect_equal(generics::tidy(result), td)
  expect_identical(stats::df.residual(result), Inf)
})
