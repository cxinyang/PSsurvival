# Bootstrap Variance Estimation --------------------------------------------

#' Generate Bootstrap Sample Indices
#'
#' @description
#' Helper function to generate bootstrap sample indices for either full sample
#' resampling (observational studies) or stratified resampling within treatment
#' groups (randomized controlled trials with fixed allocation).
#'
#' @param data Data frame to bootstrap from.
#' @param treatment_var Name of treatment variable.
#' @param boot_level Bootstrap sampling level: "full" (default) samples from
#'   entire dataset, "strata" samples within each treatment group preserving
#'   group sizes.
#'
#' @return Integer vector of bootstrap indices of length \code{nrow(data)}.
#'
#' @keywords internal
generate_boot_indices <- function(data, treatment_var, boot_level = "full") {
  n <- nrow(data)

  if (boot_level == "full") {
    # Standard bootstrap: sample n observations from full dataset
    return(sample(1:n, size = n, replace = TRUE))

  } else if (boot_level == "strata") {
    # Stratified bootstrap: sample within each treatment group
    treatment <- data[[treatment_var]]
    treatment_levels <- sort(unique(treatment))

    boot_indices <- integer(n)
    idx <- 1
    for (grp in treatment_levels) {
      grp_indices <- which(treatment == grp)
      n_grp <- length(grp_indices)
      # Sample n_grp observations with replacement from group grp
      sampled <- sample(grp_indices, size = n_grp, replace = TRUE)
      boot_indices[idx:(idx + n_grp - 1)] <- sampled
      idx <- idx + n_grp
    }
    return(boot_indices)

  } else {
    stop("boot_level must be 'full' or 'strata'", call. = FALSE)
  }
}


#' Bootstrap Variance Estimation for Weibull Survival Functions
#'
#' @description
#' Estimates bootstrap variance for counterfactual survival functions from
#' \code{surv_weibull()}. Supports binary and multiple treatment groups.
#'
#' @param data Data frame used in original \code{surv_weibull()} call.
#' @param surv_result Output from \code{surv_weibull()}.
#' @param treatment_var Name of treatment variable.
#' @param time_var Name of time variable.
#' @param event_var Name of event indicator variable.
#' @param censoring_formula Formula for censoring score model.
#' @param censoring_control Control parameters for \code{survreg()}.
#'   Default \code{list(maxiter = 350)}.
#' @param B Number of bootstrap iterations. Default 100.
#' @param parallel Logical. If TRUE, use parallel computation via \code{mclapply}.
#'   Default FALSE.
#' @param mc.cores Number of cores for parallel computation. Default 2.
#' @param seed Optional random seed for reproducibility. Ensures identical results
#'   across runs and between sequential and parallel execution. Default NULL.
#' @param boot_level Bootstrap sampling level: "full" (default) or "strata".
#'   "full" resamples from entire dataset (observational studies). "strata"
#'   resamples within treatment groups preserving group sizes (RCTs).
#'
#' @return List containing:
#'   \item{var_matrix}{Matrix [time x group] of bootstrap variances.}
#'   \item{se_matrix}{Matrix [time x group] of bootstrap standard errors.}
#'   \item{boot_samples}{List of length B with survival matrices from each iteration.}
#'   \item{boot_allocation}{Matrix [B x group] of sample sizes for each group in each bootstrap iteration.}
#'   \item{n_success_by_group}{Matrix [time x group] of successful iterations.}
#'   \item{n_failed_by_group}{Matrix [time x group] of failed iterations.}
#'   \item{B}{Total bootstrap iterations.}
#'
#' @details
#' Each bootstrap iteration resamples observations with replacement and re-estimates
#' propensity scores, weights, and survival functions using the same specifications
#' as the original analysis. Treatment groups may fail independently if not sampled
#' or if models fail to converge.
#'
#' @keywords internal
var_surv_weibull_bootstrap <- function(data, surv_result, treatment_var, time_var, event_var,
                                       censoring_formula,
                                       censoring_control = list(maxiter = 350),
                                       B = 100,
                                       parallel = FALSE,
                                       mc.cores = 2,
                                       seed = NULL,
                                       boot_level = "full") {

  # Extract specifications
  eval_times <- surv_result$eval_times
  n_times <- length(eval_times)
  treatment_levels <- surv_result$treatment_levels
  n_levels <- surv_result$n_levels
  estimand <- surv_result$estimand
  weight_result <- surv_result$weight_result

  if (surv_result$ps_result$n_levels == 2) {
    ps_formula <- surv_result$ps_result$ps_model$formula
  } else {
    ps_formula <- stats::formula(surv_result$ps_result$ps_model)
  }

  trim_method <- weight_result$trim_method
  delta <- weight_result$delta
  alpha <- weight_result$alpha
  att_group <- weight_result$att_group

  n <- nrow(data)

  # Generate iteration-specific seeds
  if (!is.null(seed)) {
    set.seed(seed)
    iteration_seeds <- sample.int(.Machine$integer.max, B)
  } else {
    iteration_seeds <- rep(NA, B)
  }

  # Bootstrap iteration
  boot_iteration <- function(b) {
    suppressWarnings(suppressMessages({
      if (!is.na(iteration_seeds[b])) {
        set.seed(iteration_seeds[b])
      }

      boot_indices <- generate_boot_indices(data, treatment_var, boot_level)
      boot_data <- data[boot_indices, , drop = FALSE]

      # Count sample sizes per treatment group
      boot_treatment <- boot_data[[treatment_var]]
      group_sizes <- sapply(treatment_levels, function(g) sum(boot_treatment == g))
      names(group_sizes) <- as.character(treatment_levels)

      result_matrix <- matrix(NA, nrow = n_times, ncol = n_levels)
      colnames(result_matrix) <- as.character(treatment_levels)

      ps_result_boot <- tryCatch({
        estimate_ps(boot_data, treatment_var, ps_formula, ps_control = list())
      }, error = function(e) NULL)

      if (is.null(ps_result_boot)) {
        return(list(survival = result_matrix, group_sizes = group_sizes))
      }

      weight_result_boot <- tryCatch({
        if (is.null(trim_method)) {
          estimate_weights(ps_result_boot, boot_data, treatment_var, estimand, att_group = att_group)
        } else if (trim_method == "symmetric") {
          estimate_weights(ps_result_boot, boot_data, treatment_var, estimand,
                          att_group = att_group, trim = "symmetric", delta = delta)
        } else if (trim_method == "asymmetric") {
          estimate_weights(ps_result_boot, boot_data, treatment_var, estimand,
                          att_group = att_group, trim = "asymmetric", alpha = alpha)
        }
      }, error = function(e) NULL)

      if (is.null(weight_result_boot)) {
        return(list(survival = result_matrix, group_sizes = group_sizes))
      }

      surv_result_boot <- tryCatch({
        surv_weibull(boot_data, time_var, event_var, treatment_var, eval_times,
                    weight_result_boot, censoring_formula, censoring_control)
      }, error = function(e) NULL)

      if (!is.null(surv_result_boot)) {
        result_matrix <- surv_result_boot$survival_matrix
      }

      return(list(survival = result_matrix, group_sizes = group_sizes))
    }))
  }

  # Run bootstrap
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package 'parallel' is required for parallel bootstrap. ",
           "Please install it or set parallel = FALSE.", call. = FALSE)
    }
    #message("Running ", B, " bootstrap iterations in parallel using ", mc.cores, " cores...")
    boot_results_list <- parallel::mclapply(1:B, boot_iteration, mc.cores = mc.cores)
  } else {
    #message("Running ", B, " bootstrap iterations sequentially...")
    boot_results_list <- lapply(1:B, boot_iteration)
  }

  # Extract survival matrices and group sizes separately
  boot_samples <- lapply(boot_results_list, function(x) x$survival)

  # Create allocation matrix (B rows x n_levels columns)
  boot_allocation <- matrix(NA, nrow = B, ncol = n_levels)
  colnames(boot_allocation) <- as.character(treatment_levels)
  for (b in 1:B) {
    if (!is.null(boot_results_list[[b]]$group_sizes)) {
      boot_allocation[b, ] <- boot_results_list[[b]]$group_sizes
    }
  }

  # Aggregate results
  boot_array <- array(NA, dim = c(n_times, n_levels, B))
  dimnames(boot_array)[[2]] <- as.character(treatment_levels)

  for (b in 1:B) {
    boot_mat <- boot_samples[[b]]
    if (!is.null(boot_mat) && is.matrix(boot_mat) && !is.null(colnames(boot_mat))) {
      for (grp in as.character(treatment_levels)) {
        if (grp %in% colnames(boot_mat)) {
          boot_array[, grp, b] <- boot_mat[, grp]
        }
      }
    }
  }

  # Calculate variance and success counts
  var_matrix <- matrix(NA, nrow = n_times, ncol = n_levels)
  colnames(var_matrix) <- as.character(treatment_levels)

  n_success_by_group <- matrix(NA, nrow = n_times, ncol = n_levels)
  colnames(n_success_by_group) <- as.character(treatment_levels)

  for (j in 1:n_levels) {
    grp <- as.character(treatment_levels[j])
    for (i in 1:n_times) {
      vals <- boot_array[i, grp, ]
      vals_valid <- vals[!is.na(vals)]
      if (length(vals_valid) > 1) {
        var_matrix[i, grp] <- stats::var(vals_valid)
      }
      n_success_by_group[i, grp] <- sum(!is.na(vals))
    }
  }

  n_failed_by_group <- B - n_success_by_group
  se_matrix <- sqrt(var_matrix)

  return(list(
    var_matrix = var_matrix,
    se_matrix = se_matrix,
    boot_samples = boot_samples,
    boot_allocation = boot_allocation,
    n_success_by_group = n_success_by_group,
    n_failed_by_group = n_failed_by_group,
    B = B
  ))
}


#' Bootstrap Variance Estimation for Cox Survival Functions
#'
#' @description
#' Estimates bootstrap variance for counterfactual survival functions from
#' \code{surv_cox()}. Supports binary and multiple treatment groups.
#'
#' @param data Data frame used in original \code{surv_cox()} call.
#' @param surv_result Output from \code{surv_cox()}.
#' @param treatment_var Name of treatment variable.
#' @param time_var Name of time variable.
#' @param event_var Name of event indicator variable.
#' @param censoring_formula Formula for censoring score model.
#' @param censoring_control Control parameters for \code{coxph()}.
#'   Default \code{list()}.
#' @param ties Tie handling method for Cox model. Default "efron".
#' @param B Number of bootstrap iterations. Default 100.
#' @param parallel Logical. If TRUE, use parallel computation via \code{mclapply}.
#'   Default FALSE.
#' @param mc.cores Number of cores for parallel computation. Default 2.
#' @param seed Optional random seed for reproducibility. Ensures identical results
#'   across runs and between sequential and parallel execution. Default NULL.
#' @param boot_level Bootstrap sampling level: "full" (default) samples from
#'   entire dataset, "strata" samples within each treatment group preserving
#'   group sizes.
#'
#' @return List containing:
#'   \item{var_matrix}{Matrix [time x group] of bootstrap variances.}
#'   \item{se_matrix}{Matrix [time x group] of bootstrap standard errors.}
#'   \item{boot_samples}{List of length B with survival matrices from each iteration.}
#'   \item{boot_allocation}{Matrix [B x group] of sample sizes for each group in each bootstrap iteration.}
#'   \item{n_success_by_group}{Matrix [time x group] of successful iterations.}
#'   \item{n_failed_by_group}{Matrix [time x group] of failed iterations.}
#'   \item{B}{Total bootstrap iterations.}
#'
#' @details
#' Each bootstrap iteration resamples observations with replacement and re-estimates
#' propensity scores, weights, and survival functions using the same specifications
#' as the original analysis. Treatment groups may fail independently if not sampled
#' or if models fail to converge.
#'
#' @keywords internal
var_surv_cox_bootstrap <- function(data, surv_result, treatment_var, time_var, event_var,
                                   censoring_formula,
                                   censoring_control = list(),
                                   ties = "efron",
                                   B = 100,
                                   parallel = FALSE,
                                   mc.cores = 2,
                                   seed = NULL,
                                   boot_level = "full") {

  # Extract specifications
  eval_times <- surv_result$eval_times
  n_times <- length(eval_times)
  treatment_levels <- surv_result$treatment_levels
  n_levels <- surv_result$n_levels
  estimand <- surv_result$estimand
  weight_result <- surv_result$weight_result

  if (surv_result$ps_result$n_levels == 2) {
    ps_formula <- surv_result$ps_result$ps_model$formula
  } else {
    ps_formula <- stats::formula(surv_result$ps_result$ps_model)
  }

  trim_method <- weight_result$trim_method
  delta <- weight_result$delta
  alpha <- weight_result$alpha
  att_group <- weight_result$att_group

  n <- nrow(data)

  # Generate iteration-specific seeds
  if (!is.null(seed)) {
    set.seed(seed)
    iteration_seeds <- sample.int(.Machine$integer.max, B)
  } else {
    iteration_seeds <- rep(NA, B)
  }

  # Bootstrap iteration
  boot_iteration <- function(b) {
    suppressWarnings(suppressMessages({
      if (!is.na(iteration_seeds[b])) {
        set.seed(iteration_seeds[b])
      }

      boot_indices <- generate_boot_indices(data, treatment_var, boot_level)
      boot_data <- data[boot_indices, , drop = FALSE]

      # Count sample sizes per treatment group
      boot_treatment <- boot_data[[treatment_var]]
      group_sizes <- sapply(treatment_levels, function(g) sum(boot_treatment == g))
      names(group_sizes) <- as.character(treatment_levels)

      result_matrix <- matrix(NA, nrow = n_times, ncol = n_levels)
      colnames(result_matrix) <- as.character(treatment_levels)

      ps_result_boot <- tryCatch({
        estimate_ps(boot_data, treatment_var, ps_formula, ps_control = list())
      }, error = function(e) NULL)

      if (is.null(ps_result_boot)) {
        return(list(survival = result_matrix, group_sizes = group_sizes))
      }

      weight_result_boot <- tryCatch({
        if (is.null(trim_method)) {
          estimate_weights(ps_result_boot, boot_data, treatment_var, estimand, att_group = att_group)
        } else if (trim_method == "symmetric") {
          estimate_weights(ps_result_boot, boot_data, treatment_var, estimand,
                          att_group = att_group, trim = "symmetric", delta = delta)
        } else if (trim_method == "asymmetric") {
          estimate_weights(ps_result_boot, boot_data, treatment_var, estimand,
                          att_group = att_group, trim = "asymmetric", alpha = alpha)
        }
      }, error = function(e) NULL)

      if (is.null(weight_result_boot)) {
        return(list(survival = result_matrix, group_sizes = group_sizes))
      }

      surv_result_boot <- tryCatch({
        surv_cox(boot_data, time_var, event_var, treatment_var, eval_times,
                weight_result_boot, censoring_formula, censoring_control, ties)
      }, error = function(e) NULL)

      if (!is.null(surv_result_boot)) {
        result_matrix <- surv_result_boot$survival_matrix
      }

      return(list(survival = result_matrix, group_sizes = group_sizes))
    }))
  }

  # Run bootstrap
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("Package 'parallel' is required for parallel bootstrap. ",
           "Please install it or set parallel = FALSE.", call. = FALSE)
    }
    #message("Running ", B, " bootstrap iterations in parallel using ", mc.cores, " cores...")
    boot_results_list <- parallel::mclapply(1:B, boot_iteration, mc.cores = mc.cores)
  } else {
    #message("Running ", B, " bootstrap iterations sequentially...")
    boot_results_list <- lapply(1:B, boot_iteration)
  }

  # Extract survival matrices and group sizes separately
  boot_samples <- lapply(boot_results_list, function(x) x$survival)

  # Create allocation matrix (B rows x n_levels columns)
  boot_allocation <- matrix(NA, nrow = B, ncol = n_levels)
  colnames(boot_allocation) <- as.character(treatment_levels)
  for (b in 1:B) {
    if (!is.null(boot_results_list[[b]]$group_sizes)) {
      boot_allocation[b, ] <- boot_results_list[[b]]$group_sizes
    }
  }

  # Aggregate results
  boot_array <- array(NA, dim = c(n_times, n_levels, B))
  dimnames(boot_array)[[2]] <- as.character(treatment_levels)

  for (b in 1:B) {
    boot_mat <- boot_samples[[b]]
    if (!is.null(boot_mat) && is.matrix(boot_mat) && !is.null(colnames(boot_mat))) {
      for (grp in as.character(treatment_levels)) {
        if (grp %in% colnames(boot_mat)) {
          boot_array[, grp, b] <- boot_mat[, grp]
        }
      }
    }
  }

  # Calculate variance and success counts
  var_matrix <- matrix(NA, nrow = n_times, ncol = n_levels)
  colnames(var_matrix) <- as.character(treatment_levels)

  n_success_by_group <- matrix(NA, nrow = n_times, ncol = n_levels)
  colnames(n_success_by_group) <- as.character(treatment_levels)

  for (j in 1:n_levels) {
    grp <- as.character(treatment_levels[j])
    for (i in 1:n_times) {
      vals <- boot_array[i, grp, ]
      vals_valid <- vals[!is.na(vals)]
      if (length(vals_valid) > 1) {
        var_matrix[i, grp] <- stats::var(vals_valid)
      }
      n_success_by_group[i, grp] <- sum(!is.na(vals))
    }
  }

  n_failed_by_group <- B - n_success_by_group
  se_matrix <- sqrt(var_matrix)

  return(list(
    var_matrix = var_matrix,
    se_matrix = se_matrix,
    boot_samples = boot_samples,
    boot_allocation = boot_allocation,
    n_success_by_group = n_success_by_group,
    n_failed_by_group = n_failed_by_group,
    B = B
  ))
}
