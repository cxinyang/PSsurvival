#' Survival Effect Estimation with Cox Censoring Scores
#'
#' @description
#' Implements Estimator I (AJE 2022, Eq. 2) combining propensity score weighting
#' and inverse probability of censoring weighting to estimate counterfactual
#' survival functions and treatment effects. Variance estimation uses bootstrap only.

# Survival Function Estimation ---------------------------------------------

#' Estimate Counterfactual Survival Functions Using Cox Censoring Scores
#'
#' @param data Data frame.
#' @param time_var Name of time variable.
#' @param event_var Name of event indicator (1 = event, 0 = censored).
#' @param treatment_var Name of treatment variable.
#' @param eval_times Numeric vector of time points. If NULL, uses all unique event times.
#' @param weight_result Output from \code{estimate_weights()}.
#' @param censoring_formula Formula for censoring score model.
#' @param censoring_control Control parameters for \code{coxph()}.
#'   Default \code{list()}.
#' @param ties Tie handling method for Cox model. Default "efron".
#'
#' @return List containing:
#'   \item{survival_matrix}{Matrix [time x J] of survival function estimates S^(j)(t).}
#'   \item{eval_times}{Time points evaluated.}
#'   \item{treatment_levels}{Treatment level values.}
#'   \item{n_levels}{Number of treatment levels.}
#'   \item{weights_by_group}{List of weight vectors by treatment group.}
#'   \item{censoring_scores_by_group}{List of censoring score vectors by group.}
#'   \item{method, estimand}{Weighting method and target estimand.}
#'   \item{censoring_result, ps_result, weight_result}{Input objects.}
#'   \item{Z_matrix}{Binary indicator matrix [n x J] for treatment groups.}
#'   \item{time_vec, event_vec}{Original time and event vectors.}
#'
#' @details
#' Estimates counterfactual survival function for each treatment group j:
#' \deqn{S^{(j)}_w(t) = 1 - \frac{\sum_i w_i I(A_i=j) \delta_i I(T_i \leq t) / K_c^{(j)}(T_i, X_i)}{\sum_i w_i I(A_i=j)}}
#'
#' Censoring scores are estimated using Cox proportional hazards models fit separately
#' within each treatment group. Variance estimation for Cox-based survival functions
#' is performed using bootstrap only (no analytical variance available).
#'
#' @keywords internal
surv_cox <- function(data, time_var, event_var, treatment_var,
                     eval_times = NULL,
                     weight_result,
                     censoring_formula,
                     censoring_control = list(),
                     ties = "efron") {
  
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
  censoring_result <- estimate_censoring_score_cox(
    data = data,
    time_var = time_var,
    treatment_var = treatment_var,
    formula = censoring_formula,
    control = censoring_control,
    ties = ties
  )

  # Create binary indicator matrix Z [n x J]
  Z_matrix <- matrix(0, nrow = n, ncol = n_levels)
  colnames(Z_matrix) <- as.character(treatment_levels)
  for (j in 1:n_levels) {
    Z_matrix[, j] <- as.numeric(treatment == treatment_levels[j])
  }
  
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
  # Note: treatment_levels, n_levels, method, estimand are included at top-level
  # for convenience even though they exist in nested objects (weight_result,
  # ps_result, censoring_result). This avoids requiring downstream code to
  # access deeply nested values repeatedly.
  result <- list(
    survival_matrix = survival_matrix,
    eval_times = eval_times,
    treatment_levels = treatment_levels,
    n_levels = n_levels,
    weights_by_group = weights_by_group,
    censoring_scores_by_group = censoring_scores_by_group,
    method = method,
    estimand = estimand,
    censoring_result = censoring_result,
    ps_result = ps_result,
    weight_result = weight_result,
    Z_matrix = Z_matrix,
    time_vec = time_vec,
    event_vec = event_vec
  )

  return(result)
}


# Survival Effect Estimation Wrapper --------------------------------------

#' Wrapper for Cox Survival Effect Estimation with Variance
#'
#' @description
#' High-level wrapper combining propensity score estimation, weighting (with optional
#' trimming), survival function estimation, and variance estimation using bootstrap.
#' Unlike Weibull-based estimation, Cox censoring scores do not support analytical
#' variance; bootstrap is always used.
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
#' @param B Number of bootstrap iterations. Default 100.
#' @param parallel Logical. Use parallel bootstrap via \code{mclapply}. Default FALSE.
#' @param mc.cores Number of cores for parallel bootstrap. Default 2.
#' @param seed Random seed for bootstrap reproducibility. Default NULL.
#' @param censoring_control Control parameters for \code{coxph()}. Default \code{list()}.
#' @param ties Tie handling method for Cox model. Default "efron".
#' @param ps_control Control parameters for PS model. Default \code{list()}.
#' @param boot_level Bootstrap sampling level: "full" (default) or "strata".
#'   "full" resamples from entire dataset (standard for observational studies). "strata"
#'   resamples within each treatment group preserving group sizes (useful when treatment assignment 
#'   follows a stratified or fixed-ratio design). Only used if \code{variance_method = "bootstrap"}.
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
#'   \item{surv_result}{Output from \code{surv_cox()}.}
#'   \item{variance_method}{Always "bootstrap" for Cox-based estimation.}
#'   \item{boot_result}{Bootstrap results.}
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
#' 6. Estimate variances using bootstrap
#'
#' **With trimming:**
#' 1. Estimate propensity scores on full data (PS_full)
#' 2. Use PS_full to identify observations to trim
#' 3. Re-estimate propensity scores on trimmed data (PS_trimmed)
#' 4. Estimate weights from PS_trimmed
#' 5. Estimate censoring scores on trimmed data
#' 6. Estimate survival functions on trimmed data
#' 7. Estimate differences (if applicable)
#' 8. Estimate variances using bootstrap
#'
#' **Treatment differences:**
#' - Binary (2 groups): Always estimates S1 - S0 (second level minus first level),
#'   ignoring contrast_matrix even if provided.
#' - Multiple groups (>2): Requires contrast_matrix to estimate differences.
#'   Returns NULL for difference_* if contrast_matrix not provided.
#'
#' **Variance estimation for differences:**
#' - Binary: Differences columns in each bootstrap sample, calculates sample variance across B trials.
#' - Multiple groups: Always uses bootstrap method.
#'
#' @keywords internal
surveff_cox <- function(data,
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
                        B = 100,
                        parallel = FALSE,
                        mc.cores = 2,
                        seed = NULL,
                        censoring_control = list(),
                        ties = "efron",
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
  surv_result <- surv_cox(
    data = data,
    time_var = time_var,
    event_var = event_var,
    treatment_var = treatment_var,
    eval_times = eval_times,
    weight_result = weight_result,
    censoring_formula = censoring_formula,
    censoring_control = censoring_control,
    ties = ties
  )

  n_times <- length(surv_result$eval_times)

  # Step 4: Estimate variances for survival functions (always bootstrap for Cox)
  variance_method <- "bootstrap"

  boot_result <- var_surv_cox_bootstrap(
    data = data,
    surv_result = surv_result,
    treatment_var = treatment_var,
    time_var = time_var,
    event_var = event_var,
    censoring_formula = censoring_formula,
    censoring_control = censoring_control,
    ties = ties,
    B = B,
    parallel = parallel,
    mc.cores = mc.cores,
    seed = seed,
    boot_level = boot_level
  )
  survival_variances <- boot_result$var_matrix

  # Set variance to NA where point estimate is NA
  for (j in 1:n_levels) {
    for (i in 1:n_times) {
      if (is.na(surv_result$survival_matrix[i, j])) {
        survival_variances[i, j] <- NA
      }
    }
  }

  # Step 5: Estimate treatment differences and associated variance
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

    # Bootstrap variance for difference
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
