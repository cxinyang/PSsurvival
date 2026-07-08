# Integration tests for the mice multiple-imputation workflow.
# The tidy() / df.residual() methods themselves are unit-tested in
# test-marCoxph.R and test-surveff.R; here we check the full
# with() -> pool() pipeline, so these are skipped when mice is absent.

test_that("mice::pool() reproduces Rubin's rules for marCoxph", {
  skip_if_not_installed("mice")

  d <- make_test_data_binary(n = 300, seed = 1)
  d$X1[1:30] <- NA
  imp <- suppressWarnings(mice::mice(d, m = 5, printFlag = FALSE, seed = 1))

  fit <- with(imp, marCoxph(
    data = data.frame(Z, X1, X2, B1, time, event),
    ps_formula = Z ~ X1 + X2 + B1, time_var = "time", event_var = "event",
    reference_level = "A", weight_method = "OW", variance_method = "robust"))
  pooled <- summary(mice::pool(fit))
  expect_equal(nrow(pooled), 1L)

  # Hand-computed Rubin (1987) rules must match pool() exactly
  tds <- lapply(fit$analyses, tidy)
  Q <- vapply(tds, function(z) z$estimate[1], numeric(1))
  U <- vapply(tds, function(z) z$std.error[1]^2, numeric(1))
  m <- length(Q)
  expect_equal(pooled$estimate, mean(Q))
  expect_equal(pooled$std.error, sqrt(mean(U) + (1 + 1 / m) * var(Q)))
})

test_that("mice::pool() pools every (group, time) parameter for surveff", {
  skip_if_not_installed("mice")

  d <- make_test_data_binary(n = 300, seed = 2)
  d$X1[1:30] <- NA
  imp <- suppressWarnings(mice::mice(d, m = 5, printFlag = FALSE, seed = 2))

  fit <- with(imp, surveff(
    data = data.frame(Z, X1, X2, B1, time, event),
    ps_formula = Z ~ X1 + X2 + B1,
    censoring_formula = survival::Surv(time, event) ~ X1,
    eval_times = c(0.5, 1), weight_method = "OW", censoring_method = "weibull"))
  # pool() emits an expected "large sample assumed" note (df = Inf for curves)
  pooled <- suppressWarnings(summary(mice::pool(fit)))

  # 2 groups + 1 contrast, each at 2 times = 6 pooled parameters
  expect_equal(nrow(pooled), 6L)
  expect_true(all(is.finite(pooled$estimate)))
  # pooled survival probabilities remain in [0, 1]
  surv <- pooled$estimate[grepl("^survival:", pooled$term)]
  expect_true(all(surv >= 0 & surv <= 1))
})
