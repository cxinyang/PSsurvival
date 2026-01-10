# Test data fixtures for unit tests
# These functions create small, reproducible datasets for testing

#' Create binary treatment test dataset
#'
#' @param n Sample size
#' @param seed Random seed for reproducibility
#' @return Data frame with binary treatment (A/B), covariates, and survival outcome
make_test_data_binary <- function(n = 150, seed = 12345) {
  set.seed(seed)

  # Covariates
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  B1 <- rbinom(n, 1, 0.5)

  # Treatment assignment (with confounding)
  ps <- plogis(0.5 * X1 - 0.3 * X2 + 0.2 * B1)
  Z_numeric <- rbinom(n, 1, ps)
  Z <- c("A", "B")[Z_numeric + 1]

  # Survival time (Weibull, treatment effect present)
  lambda <- exp(-0.5 * X1 + 0.3 * Z_numeric)
  time_event <- rweibull(n, shape = 1.5, scale = 1 / lambda)

  # Censoring time (uniform, ~30% censoring)
  censor_time <- runif(n, 0, 20)

  # Observed outcome
  event <- as.numeric(time_event < censor_time)
  time <- pmin(time_event, censor_time)

  data.frame(
    Z = Z,
    X1 = X1,
    X2 = X2,
    B1 = B1,
    time = time,
    event = event,
    stringsAsFactors = FALSE
  )
}

#' Create multiple treatment test dataset
#'
#' @param n Sample size
#' @param seed Random seed for reproducibility
#' @return Data frame with 3 treatment groups (A/B/C), covariates, and survival outcome
make_test_data_multi <- function(n = 180, seed = 54321) {
  set.seed(seed)

  # Covariates
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  B1 <- rbinom(n, 1, 0.5)

  # Treatment assignment (multinomial)
  # Log-odds for groups B and C vs A
  logit_B <- 0.4 * X1 - 0.2 * X2 + 0.1 * B1
  logit_C <- -0.3 * X1 + 0.5 * X2 - 0.2 * B1

  # Multinomial probabilities
  exp_B <- exp(logit_B)
  exp_C <- exp(logit_C)
  denom <- 1 + exp_B + exp_C

  prob_A <- 1 / denom
  prob_B <- exp_B / denom
  prob_C <- exp_C / denom

  # Sample treatment
  Z_numeric <- apply(cbind(prob_A, prob_B, prob_C), 1, function(p) {
    sample(0:2, 1, prob = p)
  })
  Z <- c("A", "B", "C")[Z_numeric + 1]

  # Survival time (treatment effects)
  lambda <- exp(-0.4 * X1 + 0.2 * (Z_numeric == 1) + 0.4 * (Z_numeric == 2))
  time_event <- rweibull(n, shape = 1.5, scale = 1 / lambda)

  # Censoring time (~30% censoring)
  censor_time <- runif(n, 0, 18)

  # Observed outcome
  event <- as.numeric(time_event < censor_time)
  time <- pmin(time_event, censor_time)

  data.frame(
    Z = Z,
    X1 = X1,
    X2 = X2,
    B1 = B1,
    time = time,
    event = event,
    stringsAsFactors = FALSE
  )
}

#' Create contrast matrix for multiple treatment comparisons
#'
#' @return Matrix with contrasts for A vs B, A vs C
make_contrast_matrix <- function() {
  contrast_mat <- matrix(
    c(1, -1,  0,   # A vs B
      1,  0, -1),  # A vs C
    nrow = 2, byrow = TRUE
  )
  colnames(contrast_mat) <- c("A", "B", "C")
  rownames(contrast_mat) <- c("A vs B", "A vs C")
  contrast_mat
}
