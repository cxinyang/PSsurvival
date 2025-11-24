#' Survival Effect Estimation with Weibull Censoring Scores
#'
#' @description
#' Implements Estimator I (AJE 2022, Eq. 2) combining propensity score weighting
#' and inverse probability of censoring weighting to estimate counterfactual
#' survival functions and treatment effects.


# Helper Functions ---------------------------------------------------------

#' Compute Etau Normalization Constant
#'
#' @param estimand Estimand type ("ATE", "ATT", or for overlap).
#' @param method Method ("IPW" or "overlap").
#' @param ps_matrix Propensity score matrix [n x J].
#' @param n_levels Number of treatment levels.
#' @param att_group Target group for ATT (if applicable).
#' @param treatment_levels Vector of treatment levels.
#'
#' @return Numeric scalar Etau.
#'
#' @keywords internal
compute_etau <- function(estimand, method, ps_matrix, n_levels, att_group = NULL,
                        treatment_levels = NULL) {
  if (estimand == "ATE") {
    # ATE: h = 1
    return(1)
  } else if (estimand == "ATT") {
    # ATT: h = e_target (propensity score of target group)
    att_col <- which(treatment_levels == att_group)
    return(mean(ps_matrix[, att_col]))
  } else if (method == "overlap") {
    # Overlap weighting
    if (n_levels == 2) {
      # Binary: h = e*(1-e)
      return(mean(ps_matrix[, 1] * ps_matrix[, 2]))
    } else {
      # Multiple groups: h = 1/(sum 1/e_l)
      return(mean(1 / rowSums(1 / ps_matrix)))
    }
  } else {
    # Default fallback
    return(1)
  }
}


# Survival Function Estimation ---------------------------------------------

#' Estimate Counterfactual Survival Functions Using Weibull Censoring Scores
#'
#' @param data Data frame.
#' @param time_var Name of time variable.
#' @param event_var Name of event indicator (1 = event, 0 = censored).
#' @param treatment_var Name of treatment variable.
#' @param eval_times Numeric vector of time points. If NULL, uses all unique event times.
#' @param weight_result Output from \code{estimate_weights()}.
#' @param censoring_formula Formula for censoring score model.
#' @param censoring_control Control parameters for \code{survreg()}.
#'   Default \code{list(maxiter = 350)}.
#'
#' @return List containing:
#'   \item{survival_matrix}{Matrix [time x J] of survival function estimates S^(j)(t).}
#'   \item{eval_times}{Time points evaluated.}
#'   \item{treatment_levels}{Treatment level values.}
#'   \item{n_levels}{Number of treatment levels.}
#'   \item{weights_by_group}{List of weight vectors by treatment group.}
#'   \item{censoring_scores_by_group}{List of censoring score vectors by group.}
#'   \item{method, estimand}{Weighting method and target estimand.}
#'   \item{Etau}{Normalization constant.}
#'   \item{censoring_result, ps_result, weight_result}{Input objects.}
#'   \item{design_matrices}{List with W (PS model) and X (censoring model).}
#'   \item{Z_matrix}{Binary indicator matrix [n x J] for treatment groups.}
#'   \item{time_vec, event_vec}{Original time and event vectors.}
#'
#' @details
#' Estimates counterfactual survival function for each treatment group j:
#' \deqn{S^{(j)}_w(t) = 1 - \frac{\sum_i w_i I(A_i=j) \delta_i I(T_i \leq t) / K_c^{(j)}(T_i, X_i)}{\sum_i w_i I(A_i=j)}}
#'
#' @keywords internal
surv_weibull <- function(data, time_var, event_var, treatment_var,
                         eval_times = NULL,
                         weight_result,
                         censoring_formula,
                         censoring_control = list(maxiter = 350)) {

  # Extract PS result from weight_result
  ps_result <- weight_result$ps_result

  # Extract basic quantities
  n <- nrow(data)
  treatment <- data[[treatment_var]]

  # Extract time and event from data
  # time_var and event_var are simple column names (validated to exist in data)
  time_vec <- data[[time_var]]
  event_vec <- data[[event_var]]

  weights <- weight_result$weights
  method <- weight_result$method
  estimand <- weight_result$estimand
  n_levels <- weight_result$n_levels
  treatment_levels <- weight_result$treatment_levels

  # Handle trimming: exclude observations where weights == 0
  trimmed <- (weights == 0)
  if (any(trimmed)) {
    message("Excluding ", sum(trimmed), " trimmed observations (weights == 0) from calculations.")
    data <- data[!trimmed, , drop = FALSE]
    treatment <- treatment[!trimmed]
    time_vec <- time_vec[!trimmed]
    event_vec <- event_vec[!trimmed]
    weights <- weights[!trimmed]
    n <- sum(!trimmed)
  }

  # Set eval_times
  if (is.null(eval_times)) {
    eval_times <- sort(unique(time_vec[event_vec == 1]))
  }
  n_times <- length(eval_times)

  # Estimate censoring scores
  censoring_result <- estimate_censoring_score_weibull(
    data = data,
    time_var = time_var,
    treatment_var = treatment_var,
    formula = censoring_formula,
    control = censoring_control
  )

  # Construct design matrices (using trimmed data)
  if (ps_result$n_levels == 2) {
    ps_formula <- ps_result$ps_model$formula
  } else {
    ps_formula <- stats::formula(ps_result$ps_model)
  }
  W_rhs <- stats::as.formula(paste("~", as.character(ps_formula)[3]))
  W <- stats::model.matrix(W_rhs, data = data)

  # Handle censoring design matrix
  # When censoring_formula is NULL (no censoring adjustment), X is not used for variance
  # Set X to NULL to signal this to var_surv_weibull_analytical
  if (is.null(censoring_formula)) {
    X <- NULL
  } else {
    cens_formula_char <- as.character(censoring_formula)
    X_rhs <- stats::as.formula(paste("~", cens_formula_char[3]))
    X <- stats::model.matrix(X_rhs, data = data)
  }

  # Extract ps vectors and matrices from ps_result
  # Note: If trimming occurred, ps_result already contains refitted PS on trimmed data
  # So ps_result$ps and ps_result$ps_matrix already have the correct dimensions (n, not n_original)
  ps <- ps_result$ps
  ps_matrix <- ps_result$ps_matrix

  # Create binary indicator matrix Z [n x J]
  Z_matrix <- matrix(0, nrow = n, ncol = n_levels)
  colnames(Z_matrix) <- as.character(treatment_levels)
  for (j in 1:n_levels) {
    Z_matrix[, j] <- as.numeric(treatment == treatment_levels[j])
  }

  # Compute Etau (normalization constant for tilting function h)
  Etau <- compute_etau(
    estimand = estimand,
    method = method,
    ps_matrix = ps_matrix,
    n_levels = n_levels,
    att_group = weight_result$att_group,
    treatment_levels = treatment_levels
  )

  # Organize weights and censoring scores by group
  weights_by_group <- vector("list", n_levels)
  names(weights_by_group) <- as.character(treatment_levels)
  censoring_scores_by_group <- vector("list", n_levels)
  names(censoring_scores_by_group) <- as.character(treatment_levels)

  for (j in 1:n_levels) {
    weights_by_group[[j]] <- weights * Z_matrix[, j]
    censoring_scores_by_group[[j]] <- censoring_result$censoring_matrix[, j]
  }

  # Estimate survival functions for each group
  survival_matrix <- matrix(NA, nrow = n_times, ncol = n_levels)
  colnames(survival_matrix) <- as.character(treatment_levels)

  for (j in 1:n_levels) {
    w_j <- weights_by_group[[j]]
    Kc_j <- censoring_scores_by_group[[j]]

    for (i in 1:n_times) {
      u <- eval_times[i]
      numerator <- sum(w_j * event_vec * (time_vec <= u) / Kc_j)
      denominator <- sum(w_j)
      survival_matrix[i, j] <- 1 - numerator / denominator
    }
  }

  # Validate survival estimates: truncate to [0, 1]
  # If < 0, set to 0; if > 1, set to 1
  for (j in 1:n_levels) {
    for (i in 1:n_times) {
      if (!is.na(survival_matrix[i, j])) {
        if (survival_matrix[i, j] < 0) {
          survival_matrix[i, j] <- 0
        } else if (survival_matrix[i, j] > 1) {
          survival_matrix[i, j] <- 1
        }
      }
    }
  }

  # Return results
  # Have some redundancy to reduce nesting
  result <- list(
    survival_matrix = survival_matrix,
    eval_times = eval_times,
    treatment_levels = treatment_levels,
    n_levels = n_levels,
    weights_by_group = weights_by_group,
    censoring_scores_by_group = censoring_scores_by_group,
    method = method,
    estimand = estimand,
    Etau = Etau,
    censoring_result = censoring_result,
    ps_result = ps_result,
    weight_result = weight_result,
    design_matrices = list(W = W, X = X),
    Z_matrix = Z_matrix,
    time_vec = time_vec,
    event_vec = event_vec,
    ps = ps,  # Trimmed ps vector
    ps_matrix = ps_matrix,  # Trimmed ps_matrix
    trimmed = trimmed  # Trimming indicator (original length)
  )

  return(result)
}


# Variance Estimation Function ---------------------------------------------

#' Compute Analytical M-Estimation Variance for Binary Treatment Survival Functions
#'
#' @description
#' Computes analytical variance estimates using M-estimation for binary treatment.
#' Calculates variances for S^(0)(t), S^(1)(t), and their difference S^(1)(t) - S^(0)(t).
#'
#' @param surv_result Output from \code{surv_weibull()} with binary treatment (2 levels).
#'
#' @return List containing:
#'   \item{var_matrix}{Matrix [time x 3] of variances: [var(S0), var(S1), var(S1-S0)].}
#'   \item{se_matrix}{Matrix [time x 3] of standard errors: [se(S0), se(S1), se(S1-S0)].}
#'   \item{influence_components}{List of Itheta and Igamma arrays for delta variance.}
#'   \item{Etau}{Normalization constant.}
#'   \item{n}{Sample size after trimming.}
#'
#' @details
#' Implements M-estimation variance for binary treatment survival functions.
#' For each group j: \deqn{I_j = \frac{1}{E_\tau}(I_{\tau,j} + I_{\theta_j} + I_{\gamma_j} + I_{\beta,j})}
#' \deqn{Var(S^{(j)}) = \sum I_j^2 / n^2}
#'
#' For the difference: \deqn{I_{diff} = \frac{1}{E_\tau}(I_{\tau,diff} + I_{\theta_1} - I_{\theta_0} + I_{\gamma_1} - I_{\gamma_0} + I_{\beta,diff})}
#' \deqn{Var(S^{(1)} - S^{(0)}) = \sum I_{diff}^2 / n^2}
#'
#' @keywords internal
var_surv_weibull_analytical <- function(surv_result) {

  # Check binary treatment
  if (surv_result$n_levels != 2) {
    stop("var_surv_weibull_analytical() only supports binary treatment (2 levels). ",
         "Found ", surv_result$n_levels, " levels.", call. = FALSE)
  }

  # Extract quantities
  n <- nrow(surv_result$design_matrices$W)
  n_times <- length(surv_result$eval_times)
  treatment_levels <- surv_result$treatment_levels
  eval_times <- surv_result$eval_times
  time_vec <- surv_result$time_vec
  event_vec <- surv_result$event_vec
  Z_matrix <- surv_result$Z_matrix
  W <- surv_result$design_matrices$W
  X <- surv_result$design_matrices$X
  Etau <- surv_result$Etau
  survival_matrix <- surv_result$survival_matrix
  weights_by_group <- surv_result$weights_by_group
  censoring_scores_by_group <- surv_result$censoring_scores_by_group
  censoring_result <- surv_result$censoring_result
  estimand <- surv_result$estimand
  method <- surv_result$method

  # Check if no censoring adjustment (X is NULL when censoring_formula was NULL)
  no_censoring_adjustment <- is.null(X)

  # Extract ps_matrix (trimmed) from surv_result
  # ps_matrix[, 1] = P(Z = 0 | W), ps_matrix[, 2] = P(Z = 1 | W)
  ps_matrix <- surv_result$ps_matrix
  ps1 <- ps_matrix[, 2]  # P(Z = 1 | W) for all individuals

  # Initialize variance matrices: columns for S0, S1, and (S1-S0)
  var_matrix <- matrix(NA, nrow = n_times, ncol = 3)
  colnames(var_matrix) <- c(as.character(treatment_levels), "diff")

  # Compute expected information matrices (E)
  Ebeta <- crossprod(sqrt(ps1 * (1 - ps1)) * W) / n

  Z_0 <- Z_matrix[, 1]
  Z_1 <- Z_matrix[, 2]

  # Censoring model E matrices for both groups (only needed when censoring adjustment is used)
  if (!no_censoring_adjustment) {
    theta_0 <- censoring_result$parameters[[as.character(treatment_levels[1])]]$theta
    gamma_0 <- censoring_result$parameters[[as.character(treatment_levels[1])]]$gamma
    theta_1 <- censoring_result$parameters[[as.character(treatment_levels[2])]]$theta
    gamma_1 <- censoring_result$parameters[[as.character(treatment_levels[2])]]$gamma

    E_theta0 <- crossprod(sqrt(Z_0 * (time_vec^gamma_0) * exp(c(X %*% theta_0))) * X) / n
    E_gamma0 <- (1/n) * sum(Z_0 * ((1 - event_vec) / (gamma_0^2) +
                              (time_vec^gamma_0) * (log(time_vec)^2) * exp(c(X %*% theta_0))))

    E_theta1 <- crossprod(sqrt(Z_1 * (time_vec^gamma_1) * exp(c(X %*% theta_1))) * X) / n
    E_gamma1 <- (1/n) * sum(Z_1 * ((1 - event_vec) / (gamma_1^2) +
                              (time_vec^gamma_1) * (log(time_vec)^2) * exp(c(X %*% theta_1))))
  }

  # Pre-compute Hbeta for both groups at all times
  Hbeta_0_array <- matrix(NA, nrow = n_times, ncol = ncol(W))
  Hbeta_1_array <- matrix(NA, nrow = n_times, ncol = ncol(W))

  w_0 <- weights_by_group[[1]]
  w_1 <- weights_by_group[[2]]
  Kc_0 <- censoring_scores_by_group[[1]]
  Kc_1 <- censoring_scores_by_group[[2]]

  for (i in 1:n_times) {
    u <- eval_times[i]
    tau_0 <- survival_matrix[i, 1]
    tau_1 <- survival_matrix[i, 2]

    # Compute dw_i^(j)/dbeta based on weighting scheme
    # ps1 = P(Z=1|X), (1-ps1) = P(Z=0|X)
    # dps1/dbeta = ps1*(1-ps1)*W

    if (estimand == "ATE") {
      # ATE (IPW): w^(1) = 1/ps1, w^(0) = 1/(1-ps1)
      # Group 1 (treated): d(1/ps1)/dbeta = -(1-ps1)/ps1 * W
      dw_dbeta_1 <- -(1 - ps1) / ps1 * W
      # Group 0 (control): d(1/(1-ps1))/dbeta = ps1/(1-ps1) * W
      dw_dbeta_0 <- ps1 / (1 - ps1) * W
    } else if (estimand == "ATT") {
      att_group <- surv_result$weight_result$att_group
      att_group_idx <- which(treatment_levels == att_group)

      if (att_group_idx == 2) {
        # Target is group 1: w^(1) = 1, w^(0) = ps1/(1-ps1)
        dw_dbeta_1 <- 0 * W
        dw_dbeta_0 <- ps1 / (1 - ps1) * W
      } else {
        # Target is group 0: w^(1) = (1-ps1)/ps1, w^(0) = 1
        dw_dbeta_1 <- -(1 - ps1) / ps1 * W
        dw_dbeta_0 <- 0 * W
      }
    } else if (method == "overlap") {
      # Overlap: w^(1) = 1-ps1, w^(0) = ps1
      # Group 1: d(1-ps1)/dbeta = -ps1*(1-ps1)*W
      dw_dbeta_1 <- -ps1 * (1 - ps1) * W
      # Group 0: d(ps1)/dbeta = ps1*(1-ps1)*W
      dw_dbeta_0 <- ps1 * (1 - ps1) * W
    } else {
      stop("Unknown weighting scheme for binary treatment", call. = FALSE)
    }

    # Compute Hbeta for both groups
    Hbeta_0_array[i, ] <- (1/n) * t(Z_0 * (event_vec * (time_vec <= u) / Kc_0 - (1 - tau_0))) %*% dw_dbeta_0
    Hbeta_1_array[i, ] <- (1/n) * t(Z_1 * (event_vec * (time_vec <= u) / Kc_1 - (1 - tau_1))) %*% dw_dbeta_1
  }

  # Store influence function components for delta variance computation
  # (only used when censoring adjustment is applied)
  Itheta_0_array <- matrix(0, nrow = n, ncol = n_times)
  Igamma_0_array <- matrix(0, nrow = n, ncol = n_times)
  Itheta_1_array <- matrix(0, nrow = n, ncol = n_times)
  Igamma_1_array <- matrix(0, nrow = n, ncol = n_times)

  # Compute variance for each group at each time
  for (i in 1:n_times) {
    u <- eval_times[i]
    tau_0 <- survival_matrix[i, 1]
    tau_1 <- survival_matrix[i, 2]

    # Itau: Compute separately for each group (not as difference)
    Itau0 <- w_0 * (event_vec * (time_vec <= u) / Kc_0 - (1 - tau_0))
    Itau1 <- w_1 * (event_vec * (time_vec <= u) / Kc_1 - (1 - tau_1))

    # Ibeta: Compute separately for each group
    Ibeta0 <- ((Z_1 - ps1) * W) %*% t(Hbeta_0_array[i, ] %*% solve(Ebeta, tol = 1e-25))
    Ibeta1 <- ((Z_1 - ps1) * W) %*% t(Hbeta_1_array[i, ] %*% solve(Ebeta, tol = 1e-25))

    if (no_censoring_adjustment) {
      # No censoring model: influence function only has Itau and Ibeta components
      I_0 <- (1 / Etau) * (c(Itau0) + c(Ibeta0))
      I_1 <- (1 / Etau) * (c(Itau1) + c(Ibeta1))

      var_matrix[i, 1] <- sum(I_0^2) / (n^2)
      var_matrix[i, 2] <- sum(I_1^2) / (n^2)

      # Variance for difference
      Itau_diff <- Itau1 - Itau0
      Ibeta_diff <- ((Z_1 - ps1) * W) %*% t((Hbeta_1_array[i, ] - Hbeta_0_array[i, ]) %*% solve(Ebeta, tol = 1e-25))
      I_diff <- (1 / Etau) * (c(Itau_diff) + c(Ibeta_diff))

      var_matrix[i, 3] <- sum(I_diff^2) / (n^2)

    } else {
      # With censoring model: compute censoring-related influence components

      # Derivatives (H) for group 0
      Htheta_0 <- (1/n) * t((w_0 * event_vec * (time_vec <= u)) / Kc_0 *
                             (time_vec^gamma_0) * exp(c(X %*% theta_0))) %*% X
      Hgamma_0 <- (1/n) * sum((w_0 * event_vec * (time_vec <= u)) / Kc_0 *
                               (time_vec^gamma_0) * log(time_vec) * exp(c(X %*% theta_0)))

      # Influence functions (I) for group 0
      Itheta_0 <- (Z_0 * (((1 - event_vec) - (time_vec^gamma_0) * exp(c(X %*% theta_0))) * X)) %*%
                  t(Htheta_0 %*% solve(E_theta0, tol = 1e-25))
      Igamma_0 <- Hgamma_0 * (1 / E_gamma0) * Z_0 *
                  ((1 - event_vec) * (1/gamma_0 + log(time_vec)) -
                   (time_vec^gamma_0) * log(time_vec) * exp(c(X %*% theta_0)))

      Itheta_0_array[, i] <- c(Itheta_0)
      Igamma_0_array[, i] <- c(Igamma_0)

      # Derivatives (H) for group 1
      Htheta_1 <- (1/n) * t((w_1 * event_vec * (time_vec <= u)) / Kc_1 *
                             (time_vec^gamma_1) * exp(c(X %*% theta_1))) %*% X
      Hgamma_1 <- (1/n) * sum((w_1 * event_vec * (time_vec <= u)) / Kc_1 *
                               (time_vec^gamma_1) * log(time_vec) * exp(c(X %*% theta_1)))

      # Influence functions (I) for group 1
      Itheta_1 <- (Z_1 * (((1 - event_vec) - (time_vec^gamma_1) * exp(c(X %*% theta_1))) * X)) %*%
                  t(Htheta_1 %*% solve(E_theta1, tol = 1e-25))
      Igamma_1 <- Hgamma_1 * (1 / E_gamma1) * Z_1 *
                  ((1 - event_vec) * (1/gamma_1 + log(time_vec)) -
                   (time_vec^gamma_1) * log(time_vec) * exp(c(X %*% theta_1)))

      Itheta_1_array[, i] <- c(Itheta_1)
      Igamma_1_array[, i] <- c(Igamma_1)

      # Combined influence function for each group
      I_0 <- (1 / Etau) * (c(Itau0) + c(Itheta_0) + c(Igamma_0) + c(Ibeta0))
      I_1 <- (1 / Etau) * (c(Itau1) + c(Itheta_1) + c(Igamma_1) + c(Ibeta1))

      var_matrix[i, 1] <- sum(I_0^2) / (n^2)
      var_matrix[i, 2] <- sum(I_1^2) / (n^2)

      # Compute variance for the difference (S1 - S0)
      Itau_diff <- Itau1 - Itau0
      Ibeta_diff <- ((Z_1 - ps1) * W) %*% t((Hbeta_1_array[i, ] - Hbeta_0_array[i, ]) %*% solve(Ebeta, tol = 1e-25))
      I_diff <- (1 / Etau) * (c(Itau_diff) + c(Itheta_1) - c(Itheta_0) + c(Igamma_1) - c(Igamma_0) + c(Ibeta_diff))

      var_matrix[i, 3] <- sum(I_diff^2) / (n^2)
    }
  }

  se_matrix <- sqrt(var_matrix)

  return(list(
    var_matrix = var_matrix,
    se_matrix = se_matrix,
    influence_components = list(
      "1" = list(Itheta = Itheta_0_array, Igamma = Igamma_0_array),
      "2" = list(Itheta = Itheta_1_array, Igamma = Igamma_1_array)
    ),
    Etau = Etau,
    n = n
  ))
}


# Survival Effect Estimation Wrapper --------------------------------------

#' Wrapper for Weibull Survival Effect Estimation with Variance
#'
#' @description
#' High-level wrapper combining propensity score estimation, weighting (with optional
#' trimming), survival function estimation, and variance estimation (analytical or bootstrap).
#'
#' @param data Data frame (possibly processed by data validation).
#' @param treatment_var Name of treatment variable.
#' @param ps_formula Formula for propensity score model.
#' @param time_var Name of time variable.
#' @param event_var Name of event indicator (1 = event, 0 = censored).
#' @param censoring_formula Formula for censoring score model.
#' @param eval_times Numeric vector of time points. If NULL, uses all unique event times.
#' @param estimand Target estimand: "ATE", "ATT", or "overlap". Default "ATE".
#' @param att_group Target group for ATT (required if estimand = "ATT").
#' @param trim Trimming method: "symmetric" or "asymmetric". Default NULL (no trimming).
#' @param delta Threshold for symmetric trimming. Required if trim = "symmetric".
#' @param alpha Percentile for asymmetric trimming. Required if trim = "asymmetric".
#' @param contrast_matrix Optional matrix for treatment differences. Each row represents
#'   one contrast with exactly two non-zero elements: -1 and 1. Column names must
#'   match treatment levels. For binary treatment, ignored (always estimates S1-S0).
#'   For >2 groups, required to estimate differences (otherwise returns NULL).
#' @param variance_method Variance estimation method: "analytical" (binary only) or
#'   "bootstrap". Default "analytical" for binary, "bootstrap" for multiple groups.
#' @param B Number of bootstrap iterations. Default 100. Ignored if variance_method = "analytical".
#' @param parallel Logical. Use parallel bootstrap via \code{mclapply}. Default FALSE.
#' @param mc.cores Number of cores for parallel bootstrap. Default 2.
#' @param seed Random seed for bootstrap reproducibility. Default NULL.
#' @param censoring_control Control parameters for \code{survreg()}. Default \code{list(maxiter = 350)}.
#' @param ps_control Control parameters for PS model. Default \code{list()}.
#' @param boot_level Bootstrap sampling level: "full" (default) or "strata".
#'   "full" resamples from entire dataset (observational studies). "strata"
#'   resamples within treatment groups preserving group sizes (RCTs). Only
#'   used if variance_method = "bootstrap".
#'
#' @return List containing:
#'   \item{survival_estimates}{Matrix [time x J] of survival function estimates.}
#'   \item{survival_variances}{Matrix [time x J] of survival function variances.}
#'   \item{difference_estimates}{Matrix [time x K] of treatment differences (NULL for >2 groups without contrast_matrix).}
#'   \item{difference_variances}{Matrix [time x K] of difference variances (NULL for >2 groups without contrast_matrix).}
#'   \item{eval_times}{Time points evaluated.}
#'   \item{treatment_levels}{Treatment level values.}
#'   \item{n_levels}{Number of treatment levels.}
#'   \item{contrast_matrix}{Contrast matrix used (NULL if not applicable).}
#'   \item{surv_result}{Output from \code{surv_weibull()}.}
#'   \item{variance_method}{Method used for variance estimation.}
#'   \item{boot_result}{Bootstrap results (NULL if variance_method = "analytical").}
#'
#' @details
#' This function implements the complete estimation workflow:
#'
#' **Without trimming:**
#' 1. Estimate propensity scores on full data
#' 2. Estimate weights from PS
#' 3. Estimate censoring scores on full data
#' 4. Estimate survival functions
#' 5. Estimate differences (if applicable)
#' 6. Estimate variances
#'
#' **With trimming:**
#' 1. Estimate propensity scores on full data (PS_full)
#' 2. Use PS_full to identify observations to trim
#' 3. Re-estimate propensity scores on trimmed data (PS_trimmed)
#' 4. Estimate weights from PS_trimmed
#' 5. Estimate censoring scores on trimmed data
#' 6. Estimate survival functions on trimmed data
#' 7. Estimate differences (if applicable)
#' 8. Estimate variances
#'
#' **Treatment differences:**
#' - Binary (2 groups): Always estimates S1 - S0 (second level minus first level),
#'   ignoring contrast_matrix even if provided.
#' - Multiple groups (>2): Requires contrast_matrix to estimate differences.
#'   Returns NULL for difference_* if contrast_matrix not provided.
#'
#' **Variance estimation for differences:**
#' - Binary with analytical: Uses influence function approach with proper
#'   correlation structure (NOT sum of variances).
#' - Binary with bootstrap: Differences columns in each bootstrap sample,
#'   calculates sample variance across B trials.
#' - Multiple groups: Always uses bootstrap method.
#'
#' @keywords internal
surveff_weibull <- function(data,
                            treatment_var,
                            ps_formula,
                            time_var,
                            event_var,
                            censoring_formula,
                            eval_times = NULL,
                            estimand = "ATE",
                            att_group = NULL,
                            trim = NULL,
                            delta = NULL,
                            alpha = NULL,
                            contrast_matrix = NULL,
                            variance_method = NULL,
                            B = 100,
                            parallel = FALSE,
                            mc.cores = 2,
                            seed = NULL,
                            censoring_control = list(maxiter = 350),
                            ps_control = list(),
                            boot_level = "full") {

  # Step 1: Estimate propensity scores on full data
  ps_result_full <- estimate_ps(data, treatment_var, ps_formula, ps_control)

  n_levels <- ps_result_full$n_levels
  treatment_levels <- ps_result_full$treatment_levels

  # Step 2: Estimate weights (potentially trimmed)
  weight_result <- estimate_weights(
    ps_result = ps_result_full,
    data = data,
    treatment_var = treatment_var,
    estimand = estimand,
    att_group = att_group,
    trim = trim,
    delta = delta,
    alpha = alpha
  )

  # Step 3: Estimate survival functions
  # Note: weight_result$ps_result contains the correct PS (refitted after trimming if applicable)
  surv_result <- surv_weibull(
    data = data,
    time_var = time_var,
    event_var = event_var,
    treatment_var = treatment_var,
    eval_times = eval_times,
    weight_result = weight_result,
    censoring_formula = censoring_formula,
    censoring_control = censoring_control
  )

  n_times <- length(surv_result$eval_times)

  # Step 4: Determine variance method
  if (is.null(variance_method)) {
    variance_method <- if (n_levels == 2) "analytical" else "bootstrap"
  }

  # Validate variance method
  if (variance_method == "analytical" && n_levels != 2) {
    stop("Analytical variance estimation only supports binary treatment. ",
         "Use variance_method = 'bootstrap' for multiple groups.", call. = FALSE)
  }

  # Step 5: Estimate variances for survival functions
  if (variance_method == "analytical") {
    var_result <- var_surv_weibull_analytical(surv_result)
    survival_variances <- var_result$var_matrix[, 1:n_levels, drop = FALSE]
    boot_result <- NULL
  } else {
    # Bootstrap variance
    boot_result <- var_surv_weibull_bootstrap(
      data = data,
      surv_result = surv_result,
      treatment_var = treatment_var,
      time_var = time_var,
      event_var = event_var,
      censoring_formula = censoring_formula,
      censoring_control = censoring_control,
      B = B,
      parallel = parallel,
      mc.cores = mc.cores,
      seed = seed,
      boot_level = boot_level
    )
    survival_variances <- boot_result$var_matrix
  }

  # Set variance to NA where point estimate is NA
  for (j in 1:n_levels) {
    for (i in 1:n_times) {
      if (is.na(surv_result$survival_matrix[i, j])) {
        survival_variances[i, j] <- NA
      }
    }
  }

  # Step 6: Estimate treatment differences and associated variance
  difference_estimates <- NULL
  difference_variances <- NULL

  if (n_levels == 2) {
    # Binary treatment: always estimate S1 - S0 (ignore contrast_matrix)
    diff_values <- surv_result$survival_matrix[, 2] - surv_result$survival_matrix[, 1]
    # Set to NA if either group is NA
    for (i in 1:n_times) {
      if (is.na(surv_result$survival_matrix[i, 1]) || is.na(surv_result$survival_matrix[i, 2])) {
        diff_values[i] <- NA
      }
    }
    difference_estimates <- matrix(diff_values, ncol = 1)
    colnames(difference_estimates) <- paste0(treatment_levels[2], " vs ", treatment_levels[1])

    # Variance for difference
    if (variance_method == "analytical") {
      # Use analytical variance (column 3 of var_matrix from var_surv_weibull_analytical)
      difference_variances <- matrix(var_result$var_matrix[, 3], ncol = 1)
      colnames(difference_variances) <- paste0(treatment_levels[2], " vs ", treatment_levels[1])
    } else {
      # Bootstrap: difference columns in each sample, calculate variance
      diff_array <- array(NA, dim = c(n_times, B))
      for (b in 1:B) {
        boot_mat <- boot_result$boot_samples[[b]]
        if (!is.null(boot_mat) && is.matrix(boot_mat) && ncol(boot_mat) == 2) {
          diff_array[, b] <- boot_mat[, 2] - boot_mat[, 1]
        }
      }
      diff_var <- apply(diff_array, 1, function(x) {
        vals_valid <- x[!is.na(x)]
        if (length(vals_valid) > 1) stats::var(vals_valid) else NA
      })
      difference_variances <- matrix(diff_var, ncol = 1)
      colnames(difference_variances) <- paste0(treatment_levels[2], " vs ", treatment_levels[1])
    }

    # Set difference variance to NA where difference estimate is NA
    for (i in 1:n_times) {
      if (is.na(difference_estimates[i, 1])) {
        difference_variances[i, 1] <- NA
      }
    }

  } else if (n_levels > 2) {
    # Multiple groups: require contrast_matrix
    if (!is.null(contrast_matrix)) {
      # Validate contrast_matrix
      if (!is.matrix(contrast_matrix)) {
        stop("contrast_matrix must be a matrix", call. = FALSE)
      }
      if (ncol(contrast_matrix) != n_levels) {
        stop("contrast_matrix must have ", n_levels, " columns matching number of treatment groups", call. = FALSE)
      }
      if (!all(colnames(contrast_matrix) %in% as.character(treatment_levels))) {
        stop("contrast_matrix column names must match treatment levels", call. = FALSE)
      }

      # Reorder contrast_matrix columns to match treatment_levels order
      # This ensures alignment when multiplying survival_matrix %*% c_vec
      # treatment_levels are sorted alphabetically (from estimate_ps)
      contrast_matrix <- contrast_matrix[, as.character(treatment_levels), drop = FALSE]

      # Check each row has exactly two non-zero elements: -1 and 1
      for (i in 1:nrow(contrast_matrix)) {
        nonzero <- contrast_matrix[i, ][contrast_matrix[i, ] != 0]
        if (length(nonzero) != 2 || !all(sort(nonzero) == c(-1, 1))) {
          stop("Each row of contrast_matrix must contain exactly two non-zero elements: -1 and 1",
               call. = FALSE)
        }
      }

      # Estimate differences using bootstrap
      n_contrasts <- nrow(contrast_matrix)
      difference_estimates <- matrix(NA, nrow = n_times, ncol = n_contrasts)
      difference_variances <- matrix(NA, nrow = n_times, ncol = n_contrasts)

      contrast_names <- character(n_contrasts)
      for (k in 1:n_contrasts) {
        c_vec <- contrast_matrix[k, ]
        group1_idx <- which(c_vec == -1)
        group2_idx <- which(c_vec == 1)
        contrast_names[k] <- paste0(treatment_levels[group2_idx], " vs ", treatment_levels[group1_idx])

        # Point estimate - check for NA in involved groups
        group_indices <- which(c_vec != 0)  # Groups involved in this contrast
        for (t_idx in 1:n_times) {
          involved_values <- surv_result$survival_matrix[t_idx, group_indices]
          if (any(is.na(involved_values))) {
            difference_estimates[t_idx, k] <- NA
          } else {
            # Only sum over involved groups to avoid NaN contamination
            difference_estimates[t_idx, k] <- sum(surv_result$survival_matrix[t_idx, group_indices] * c_vec[group_indices])
          }
        }

        # Bootstrap variance
        diff_array <- array(NA, dim = c(n_times, B))
        # Identify which groups are involved in this contrast
        groups_involved <- as.character(treatment_levels[group_indices])

        for (b in 1:B) {
          boot_mat <- boot_result$boot_samples[[b]]
          # Check that boot_mat has all groups involved in this specific contrast
          if (!is.null(boot_mat) && is.matrix(boot_mat) &&
              all(groups_involved %in% colnames(boot_mat))) {
            # Only use columns involved in this contrast to avoid NA contamination
            boot_mat_involved <- boot_mat[, groups_involved, drop = FALSE]
            c_vec_involved <- c_vec[group_indices]
            # Calculate difference using only involved groups
            diff_array[, b] <- as.numeric(boot_mat_involved %*% c_vec_involved)
          }
        }
        diff_var <- apply(diff_array, 1, function(x) {
          vals_valid <- x[!is.na(x)]
          if (length(vals_valid) > 1) stats::var(vals_valid) else NA
        })
        difference_variances[, k] <- diff_var
      }

      # Set difference variance to NA where difference estimate is NA
      for (k in 1:n_contrasts) {
        for (i in 1:n_times) {
          if (is.na(difference_estimates[i, k])) {
            difference_variances[i, k] <- NA
          }
        }
      }

      colnames(difference_estimates) <- contrast_names
      colnames(difference_variances) <- contrast_names
    }
    # else: contrast_matrix is NULL, so difference_estimates and difference_variances remain NULL
  }

  # Return results
  return(list(
    survival_estimates = surv_result$survival_matrix,
    survival_variances = survival_variances,
    difference_estimates = difference_estimates,
    difference_variances = difference_variances,
    eval_times = surv_result$eval_times,
    treatment_levels = treatment_levels,
    n_levels = n_levels,
    contrast_matrix = contrast_matrix,
    surv_result = surv_result,
    variance_method = variance_method,
    boot_result = boot_result
  ))
}
