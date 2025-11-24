#' Bootstrap Variance Estimation for Marginal Cox Model
#'
#' @description
#' Estimates variance of marginal hazard ratios via bootstrap resampling.

# Bootstrap Variance Estimation ----------------------------------------------

#' Bootstrap Variance for Marginal Cox Model
#'
#' @description
#' Performs bootstrap resampling to estimate variance of log hazard ratios
#' from weighted marginal Cox model. Supports full (unstratified) and
#' stratified bootstrap by treatment group.
#'
#' @param data A data.frame containing the complete-case analysis data.
#' @param treatment_var Character string specifying the treatment variable name.
#' @param time_var Character string specifying the time variable name.
#' @param event_var Character string specifying the event variable name.
#' @param ps_formula A formula object for the propensity score model.
#' @param treatment_levels Vector of treatment levels (from main fit_marginal_cox).
#' @param reference_level Reference treatment level (from main fit_marginal_cox).
#' @param estimand Character string: "ATE", "ATT", or "overlap".
#' @param att_group For ATT, which group to target. NULL otherwise.
#' @param trim Trimming method: NULL, "symmetric", or "asymmetric".
#' @param delta Symmetric trimming threshold (NULL uses defaults).
#' @param alpha Asymmetric trimming threshold (NULL uses defaults).
#' @param boot_level Bootstrap sampling level: "full" (default) or "strata".
#'   "full" resamples from entire dataset (observational studies). "strata"
#'   resamples within treatment groups preserving group sizes (RCTs).
#' @param B Integer number of bootstrap iterations. Default 100.
#' @param robust Logical. Use robust variance in Cox model? Default TRUE.
#' @param parallel Logical. If TRUE, use parallel computation via mclapply.
#'   Default FALSE.
#' @param mc.cores Integer number of cores for parallel processing. Default 2.
#' @param seed Optional random seed for reproducibility. Default NULL.
#'
#' @return A list containing:
#'   \item{boot_samples}{List of length B with hr_estimates from each iteration.}
#'   \item{boot_allocation}{Matrix (B x n_levels) of group sample sizes per trial.}
#'   \item{n_used_boot}{Matrix (B x n_levels) of sample sizes used in Cox model
#'     per trial (after trimming).}
#'   \item{events_used_boot}{Matrix (B x n_levels) of event counts used in Cox
#'     model per trial (after trimming).}
#'   \item{n_success_by_group}{Named integer vector of successful estimates per group
#'     (non-NA across B trials).}
#'   \item{B}{Number of bootstrap iterations.}
#'   \item{boot_level}{Bootstrap method used.}
#'   \item{treatment_levels}{Treatment levels used.}
#'   \item{reference_level}{Reference level used.}
#'
#' @details
#' \strong{Bootstrap Workflow:}
#' For each bootstrap iteration:
#' \enumerate{
#'   \item Resample data with replacement (full or stratified by treatment)
#'   \item Estimate propensity scores on bootstrap sample
#'   \item Calculate weights (with optional trimming)
#'   \item Fit marginal Cox model using fit_marginal_cox with functionality="boot"
#'   \item Record estimates, sample sizes, and event counts
#' }
#'
#' \strong{Parallel Processing:}
#' Uses parallel::mclapply for parallel bootstrap. Set ncores > 1 to enable.
#' Note: mclapply uses forking (not available on Windows).
#'
#' \strong{Error Handling:}
#' Uses fit_marginal_cox(..., functionality="boot") which returns NA for failed
#' estimates instead of throwing errors. This ensures all B trials complete.
#'
#' @keywords internal
var_marginalcox_bootstrap <- function(data, treatment_var, time_var, event_var,
                                      ps_formula, treatment_levels, reference_level,
                                      estimand = "ATE", att_group = NULL,
                                      trim = NULL, delta = NULL, alpha = NULL,
                                      boot_level = "full", B = 100,
                                      robust = TRUE, parallel = FALSE, mc.cores = 2,
                                      seed = NULL) {

  # Basic dimensions
  n <- nrow(data)
  n_levels <- length(treatment_levels)
  non_ref_levels <- setdiff(treatment_levels, reference_level)
  n_coef <- length(non_ref_levels)

  # Generate iteration-specific seeds
  if (!is.null(seed)) {
    set.seed(seed)
    iteration_seeds <- sample.int(.Machine$integer.max, B)
  } else {
    iteration_seeds <- rep(NA, B)
  }

  # Helper: return structure for failed bootstrap trials
  failed_boot <- function(group_sizes, boot_indices) {
    list(
      hr_estimates = stats::setNames(rep(NA_real_, n_coef), paste0("trt", non_ref_levels)),
      group_sizes = group_sizes,
      n_used = stats::setNames(rep(0, n_levels), as.character(treatment_levels)),
      events_used = stats::setNames(rep(0, n_levels), as.character(treatment_levels)),
      boot_indices = boot_indices
    )
  }

  # Bootstrap function for single iteration
  boot_iteration <- function(b) {
    if (!is.na(iteration_seeds[b])) {
      set.seed(iteration_seeds[b])
    }

    # Resample data
    if (boot_level == "full") {
      boot_indices <- sample(1:n, n, replace = TRUE)
    } else {  # strata
      treatment_original <- data[[treatment_var]]
      boot_indices <- unlist(lapply(treatment_levels, function(lev) {
        idx_lev <- which(treatment_original == lev)
        if (length(idx_lev) > 1) {
          sample(idx_lev, length(idx_lev), replace = TRUE)
          } else {idx_lev}
      }))
    }

    boot_data <- data[boot_indices, , drop = FALSE]
    boot_treatment <- boot_data[[treatment_var]]

    # Count group allocation (before trimming)
    group_sizes <- stats::setNames(
      sapply(treatment_levels, function(lv) sum(boot_treatment == lv)),
      as.character(treatment_levels)
    )

    # Estimate propensity scores
    ps_result_boot <- tryCatch(
      estimate_ps(boot_data, treatment_var, ps_formula, ps_control = list(trace = FALSE)),
      error = function(e) NULL
    )
    if (is.null(ps_result_boot)) return(failed_boot(group_sizes, boot_indices))

    # Calculate weights
    weights_result_boot <- tryCatch(
      estimate_weights(ps_result_boot, boot_data, treatment_var, estimand, att_group, trim, delta, alpha),
      error = function(e) NULL
    )
    if (is.null(weights_result_boot)) return(failed_boot(group_sizes, boot_indices))

    # Fit marginal Cox model
    fit_boot <- fit_marginal_cox(
      data = boot_data,
      treatment_var = treatment_var,
      time_var = time_var,
      event_var = event_var,
      weights = weights_result_boot$weights,
      treatment_levels = treatment_levels,
      reference_level = reference_level,
      robust = robust,
      functionality = "boot"
    )

    list(
      hr_estimates = fit_boot$hr_estimates,
      group_sizes = group_sizes,
      n_used = fit_boot$n_per_group_used,
      events_used = fit_boot$events_per_group_used,
      boot_indices = boot_indices
    )
  }

  # Run bootstrap iterations (suppress messages and warnings)
  boot_results_list <- suppressWarnings(suppressMessages({
    if (parallel) {
      if (!requireNamespace("parallel", quietly = TRUE)) {
        stop("Package 'parallel' is required for parallel bootstrap. ",
             "Please install it or set parallel = FALSE.", call. = FALSE)
      }
      parallel::mclapply(1:B, boot_iteration, mc.cores = mc.cores)
    } else {
      lapply(1:B, boot_iteration)
    }
  }))

  # Extract hr_estimates as boot_samples (list of vectors)
  boot_samples <- lapply(boot_results_list, function(x) x$hr_estimates)

  # Collapse results into matrices
  boot_allocation <- t(sapply(boot_results_list, function(x) x$group_sizes))
  n_used_boot <- t(sapply(boot_results_list, function(x) x$n_used))
  events_used_boot <- t(sapply(boot_results_list, function(x) x$events_used))
  boot_indices_matrix <- matrix(unlist(lapply(boot_results_list, function(x) x$boot_indices)),
                                 nrow = n, ncol = B)

  # Count successes (non-NA estimates per group)
  hr_boot_matrix <- do.call(rbind, boot_samples)
  n_success_by_group <- colSums(!is.na(hr_boot_matrix))

  # Return results
  list(
    boot_samples = boot_samples,
    boot_allocation = boot_allocation,
    n_used_boot = n_used_boot,
    events_used_boot = events_used_boot,
    boot_indices_matrix = boot_indices_matrix,
    n_success_by_group = n_success_by_group,
    B = B,
    boot_level = boot_level,
    treatment_levels = treatment_levels,
    reference_level = reference_level
  )
}
