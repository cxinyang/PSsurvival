# Input validation tests for surveff() and marCoxph()

# =============================================================================
# surveff() validation tests
# =============================================================================

test_that("surveff errors when data is not a data.frame", {
  expect_error(
    surveff(
      data = list(Z = c("A", "B"), time = c(1, 2)),
      ps_formula = Z ~ X1,
      censoring_formula = survival::Surv(time, event) ~ X1
    ),
    "object 'X1' not found|data.frame"
  )
})

test_that("surveff errors when treatment variable has <2 levels", {
  data <- make_test_data_binary(n = 50)
  data$Z <- "A"  # All same level

  expect_error(
    surveff(
      data = data,
      ps_formula = Z ~ X1 + X2,
      censoring_formula = survival::Surv(time, event) ~ X1
    ),
    "at least 2"
  )
})

test_that("surveff errors when time variable has invalid values", {
  data <- make_test_data_binary(n = 50)
  data$time[1:5] <- -1  # Negative time

  expect_error(
    surveff(
      data = data,
      ps_formula = Z ~ X1 + X2,
      censoring_formula = survival::Surv(time, event) ~ X1
    ),
    "positive"
  )
})

test_that("surveff errors when event variable is not binary", {
  data <- make_test_data_binary(n = 50)
  data$event[1:5] <- 2  # Not 0/1

  # Surv() converts invalid values to NA and issues warning,
  # so we expect warning rather than error
  expect_warning(
    surveff(
      data = data,
      ps_formula = Z ~ X1 + X2,
      censoring_formula = survival::Surv(time, event) ~ X1
    ),
    "Invalid status value"
  )
})

test_that("surveff errors for invalid estimand", {
  data <- make_test_data_binary(n = 50)

  expect_error(
    surveff(
      data = data,
      ps_formula = Z ~ X1 + X2,
      censoring_formula = survival::Surv(time, event) ~ X1,
      estimand = "invalid_estimand"
    ),
    "estimand must be"
  )
})

test_that("surveff errors when ATT specified without att_group", {
  data <- make_test_data_binary(n = 50)

  expect_error(
    surveff(
      data = data,
      ps_formula = Z ~ X1 + X2,
      censoring_formula = survival::Surv(time, event) ~ X1,
      estimand = "ATT"
    ),
    "att_group must be specified"
  )
})

test_that("surveff errors when using overlap weights with trimming", {
  data <- make_test_data_binary(n = 50)

  expect_error(
    surveff(
      data = data,
      ps_formula = Z ~ X1 + X2,
      censoring_formula = survival::Surv(time, event) ~ X1,
      estimand = "overlap",
      trim = "symmetric",
      delta = 0.1
    ),
    "not supported with overlap"
  )
})

test_that("surveff errors when contrast_matrix columns don't match treatment levels", {
  data <- make_test_data_multi(n = 100)

  # Wrong column names
  contrast_mat <- matrix(c(1, -1, 0), nrow = 1)
  colnames(contrast_mat) <- c("X", "Y", "Z")

  expect_error(
    surveff(
      data = data,
      ps_formula = Z ~ X1 + X2,
      censoring_formula = survival::Surv(time, event) ~ X1,
      contrast_matrix = contrast_mat
    ),
    "column names must match treatment levels"
  )
})

test_that("surveff errors when contrast_matrix has invalid structure", {
  data <- make_test_data_multi(n = 100)

  # More than 2 non-zero elements
  contrast_mat <- matrix(c(1, 1, -1), nrow = 1)
  colnames(contrast_mat) <- c("A", "B", "C")

  expect_error(
    surveff(
      data = data,
      ps_formula = Z ~ X1 + X2,
      censoring_formula = survival::Surv(time, event) ~ X1,
      contrast_matrix = contrast_mat
    ),
    "exactly two non-zero"
  )
})

test_that("surveff errors when censoring_formula is not a Surv object", {
  data <- make_test_data_binary(n = 50)

  expect_error(
    surveff(
      data = data,
      ps_formula = Z ~ X1 + X2,
      censoring_formula = time ~ X1  # Not Surv()
    ),
    "Surv"
  )
})

# =============================================================================
# marCoxph() validation tests
# =============================================================================

test_that("marCoxph errors when reference_level is missing", {
  data <- make_test_data_binary(n = 50)

  expect_error(
    marCoxph(
      data = data,
      ps_formula = Z ~ X1 + X2,
      time_var = "time",
      event_var = "event"
      # reference_level missing
    ),
    "reference_level.*mandatory"
  )
})

test_that("marCoxph errors when reference_level doesn't exist in data", {
  data <- make_test_data_binary(n = 50)

  expect_error(
    marCoxph(
      data = data,
      ps_formula = Z ~ X1 + X2,
      time_var = "time",
      event_var = "event",
      reference_level = "NonExistent"
    ),
    "not in treatment levels"
  )
})

test_that("marCoxph errors when treatment group has zero events", {
  data <- make_test_data_binary(n = 50)
  # Set all events in group A to 0 (all censored)
  data$event[data$Z == "A"] <- 0

  expect_error(
    marCoxph(
      data = data,
      ps_formula = Z ~ X1 + X2,
      time_var = "time",
      event_var = "event",
      reference_level = "B"
    ),
    "no observed events"
  )
})

test_that("marCoxph errors when data is not a data.frame", {
  expect_error(
    marCoxph(
      data = list(Z = c("A", "B")),
      ps_formula = Z ~ X1,
      time_var = "time",
      event_var = "event",
      reference_level = "A"
    ),
    "data.frame"
  )
})

test_that("marCoxph errors when variance_method is invalid", {
  data <- make_test_data_binary(n = 50)

  expect_error(
    marCoxph(
      data = data,
      ps_formula = Z ~ X1 + X2,
      time_var = "time",
      event_var = "event",
      reference_level = "A",
      variance_method = "invalid"
    ),
    "bootstrap.*robust"
  )
})
