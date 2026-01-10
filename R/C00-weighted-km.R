# Declare non-standard evaluation (NSE) variables used in ggplot2::aes()
# These column names are evaluated within data frame context in plot.weightedKM()
utils::globalVariables(c("Time", "Estimate", "Strata", "CI_lower", "CI_upper", "Y", "Label", "Type"))

#' Weighted Kaplan-Meier Estimation with Propensity Score Weights
#'
#' @description
#' Computes weighted Kaplan-Meier survival or cumulative incidence curves using
#' propensity score weights. Supports multiple treatment groups with various
#' weighting schemes (IPW, OW, or ATT) and optional trimming. Special case
#' \code{weight_method = "none"} provides classical (unweighted) Kaplan-Meier.
#'
#' @param data Data frame containing treatment, survival outcome, and covariates.
#' @param treatment_var Character string specifying the name of the treatment variable.
#' @param ps_formula Formula for propensity score model: \code{treatment ~ covariates}.
#'   Required unless \code{weight_method = "none"}.
#' @param time_var Character string specifying the time-to-event variable name.
#' @param event_var Character string specifying the event indicator variable name.
#'   Should be coded as 1=event, 0=censored.
#' @param weight_method Weighting method: "IPW" (inverse probability weighting),
#'   "OW" (overlap weighting), "ATT" (average treatment effect on the treated),
#'   or "none" (unweighted). Default "IPW".
#' @param att_group Target group for ATT. Required if \code{weight_method = "ATT"}.
#' @param trim Logical. To perform symmetric propensity score trimming? Default FALSE.
#'   If TRUE, symmetric trimming is applied (Crump extension for multiple treatments).
#'   See \code{\link{estimate_weights}} for trimming details. Ignored if
#'   \code{weight_method = "none"} or \code{weight_method = "OW"}. Asymmetric
#'   trimming is no longer supported due to poor statistical performance.
#' @param delta Threshold for symmetric trimming in \eqn{(0, 1/J]}, where \eqn{J} is the number
#'   of treatment levels. Default NULL uses recommended values: 0.1 for binary
#'   treatment, 0.067 for 3 groups, \eqn{1/(2J)} for \eqn{J \ge 4} (Yoshida et al., 2019).
#'   Used only if \code{trim = TRUE}.
#' @param ps_control Control parameters for propensity score model. Default \code{list()}.
#'   Ignored if \code{weight_method = "none"}.
#'
#' @details
#' **Weighting Methods:**
#'
#' The \code{weight_method} parameter specifies the target population for causal
#' inference:
#'
#' \itemize{
#'   \item \strong{IPW (Inverse Probability Weighting)}: Observations are weighted
#'     by the inverse probability of their observed treatment, \eqn{w_i = 1/e_j(X_i)}
#'     where j is the observed treatment group. Inference targets the combined
#'     population.
#'
#'   \item \strong{OW (Overlap Weighting)}: Observations are weighted by overlap
#'     weights, which extends to multiple treatment groups (Li et al., 2018; 
#'     Li and Li, 2019). Inference targets the population at clinical equipoise 
#'     (overlap population).
#'
#'   \item \strong{ATT (Average Treatment Effect on the Treated)}: IPW weights
#'     tilted toward a specified target group. Observations in the target group
#'     receive weight 1, others receive \eqn{w_i = e_{\text{target}}(X_i) / e_j(X_i)}.
#'     Inference targets the specified treatment group population.
#'
#'   \item \strong{none}: No weighting applied (all weights = 1). Produces classical
#'     unweighted Kaplan-Meier estimates stratified by treatment.
#' }
#'
#' @return Object of class "weightedKM" containing:
#'   \item{eval_times}{Numeric vector of all unique event times.}
#'   \item{surv_estimates}{Matrix [n_times x n_groups] of survival probability estimates.}
#'   \item{surv_var}{Matrix [n_times x n_groups] of variances.}
#'   \item{n_risk}{Matrix [n_times x n_groups] of weighted number at risk.}
#'   \item{n_event}{Matrix [n_times x n_groups] of weighted number of events.}
#'   \item{n_acc_event}{Matrix [n_times x n_groups] of cumulative weighted events up to each time.}
#'   \item{treatment_levels}{Sorted unique treatment values.}
#'   \item{weight_method}{Weighting method used.}
#'   \item{att_group}{Target group for ATT (NULL if not applicable).}
#'   \item{trim}{Logical indicating whether trimming was applied.}
#'   \item{delta}{Symmetric trimming threshold (NULL if trim = FALSE).}
#'   \item{n}{Number of complete cases used in analysis.}
#'   \item{ps_result}{Propensity score estimation results (NULL if \code{weight_method = "none"}).}
#'   \item{weight_result}{Weight estimation results (NULL if \code{weight_method = "none"}).}
#'   \item{weights}{Vector of weights used in estimation (all 1s if \code{weight_method = "none"}).}
#'
#' @examples
#' \donttest{
#' # Example 1: Classical (unweighted) Kaplan-Meier for binary treatment
#' data(simdata_bin)
#' result1 <- weightedKM(
#'   data = simdata_bin,
#'   treatment_var = "Z",
#'   time_var = "time",
#'   event_var = "event",
#'   weight_method = "none"
#' )
#' plot(result1)
#'
#' # Example 2: Overlap-weighted Kaplan-Meier
#' result2 <- weightedKM(
#'   data = simdata_bin,
#'   treatment_var = "Z",
#'   ps_formula = Z ~ X1 + X2 + X3 + B1 + B2,
#'   time_var = "time",
#'   event_var = "event",
#'   weight_method = "OW"
#' )
#' summary(result2)
#'
#' # Example 3: IPW-weighted Kaplan-Meier for multiple treatments
#' data(simdata_multi)
#' result3 <- weightedKM(
#'   data = simdata_multi,
#'   treatment_var = "Z",
#'   ps_formula = Z ~ X1 + X2 + X3 + B1 + B2,
#'   time_var = "time",
#'   event_var = "event",
#'   weight_method = "IPW"
#' )
#' plot(result3)
#'
#' # Example 4: ATT with symmetric trimming
#' result4 <- weightedKM(
#'   data = simdata_multi,
#'   treatment_var = "Z",
#'   ps_formula = Z ~ X1 + X2 + X3 + B1 + B2,
#'   time_var = "time",
#'   event_var = "event",
#'   weight_method = "ATT",
#'   att_group = "A",
#'   trim = TRUE,
#'   delta = 0.1
#' )
#' summary(result4)
#' }
#'
#' @references
#' Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via
#' propensity score weighting. \emph{Journal of the American Statistical Association},
#' 113(521), 390-400.
#'
#' Li, F., & Li, F. (2019). Propensity score weighting for causal inference with
#' multiple treatments. \emph{The Annals of Applied Statistics}, 13(4), 2389-2415.
#'
#' Yoshida, K., et al. (2019). Multinomial extension of propensity score trimming
#' methods: A simulation study. \emph{American Journal of Epidemiology}, 188(3),
#' 609-616.
#'
#' @export
weightedKM <- function(data,
                       treatment_var,
                       time_var,
                       event_var,
                       ps_formula = NULL,
                       weight_method = "IPW",
                       att_group = NULL,
                       trim = FALSE,
                       delta = NULL,
                       ps_control = list()) {

  # ============================================================================
  # STEP 1: Input Validation
  # ============================================================================

  # Basic parameter validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  if (!weight_method %in% c("IPW", "OW", "ATT", "none")) {
    stop("'weight_method' must be 'IPW', 'OW', 'ATT', or 'none'", call. = FALSE)
  }
  if (!is.list(ps_control)) {
    stop("'ps_control' must be a list.", call. = FALSE)
  }

  # Check treatment_var
  if (!is.character(treatment_var) || length(treatment_var) != 1 || !treatment_var %in% names(data)) {
    stop("'treatment_var' must be a single character string and exist in data.", call. = FALSE)
  }

  # Check time_var and event_var
  if (!is.character(time_var) || length(time_var) != 1 || !time_var %in% names(data)) {
    stop("'time_var' must be a single character string and exist in data.", call. = FALSE)
  }
  if (!is.character(event_var) || length(event_var) != 1 || !event_var %in% names(data)) {
    stop("'event_var' must be a single character string and exist in data.", call. = FALSE)
  }

  # Validate weight_method-specific requirements
  if (weight_method == "ATT" && is.null(att_group)) {
    stop("'att_group' is required when weight_method = 'ATT'.", call. = FALSE)
  }
  if (!is.logical(trim)) {
    stop("'trim' must be logical (TRUE or FALSE).", call. = FALSE)
  }
  if (weight_method == "OW" && trim) {
    stop("Trimming is not supported with overlap weights.", call. = FALSE)
  }
  if (weight_method == "none" && trim) {
    warning("'trim' is ignored when weight_method = 'none'.", call. = FALSE)
  }
  if (weight_method != "none" && is.null(ps_formula)) {
    stop("'ps_formula' is required when weight_method is '", weight_method, "'.", call. = FALSE)
  }
  if (trim && !is.null(delta)) {
    if (!is.numeric(delta) || delta <= 0 || delta >= 0.5) {
      stop("'delta' must be between 0 and 0.5 when trim = TRUE, or NULL to use recommended defaults.", call. = FALSE)
    }
  }

  # ============================================================================
  # STEP 2: Create Complete-Case Dataset
  # ============================================================================

  if (weight_method == "none") {
    # No PS modeling: only require treatment, time, event to be complete
    all_complete <- !is.na(data[[treatment_var]]) & !is.na(data[[time_var]]) & !is.na(data[[event_var]])

    if (sum(all_complete) == 0) {
      stop("No complete cases available for analysis.", call. = FALSE)
    }

    data_clean <- data[all_complete, , drop = FALSE]

  } else {
    # PS modeling required: validate ps_formula and require complete covariates
    # Validate ps_formula (check treatment variable matches)
    ps_validated <- validate_ps_formula(ps_formula)
    ps_treatment_var <- ps_validated$treatment_var

    if (ps_treatment_var != treatment_var) {
      stop("Treatment variable in 'ps_formula' ('", ps_treatment_var,
           "') does not match 'treatment_var' ('", treatment_var, "').", call. = FALSE)
    }

    # Create complete-case dataset (ps_formula covariates + time + event)
    ps_mf <- tryCatch(
      stats::model.frame(ps_formula, data = data, na.action = stats::na.pass),
      error = function(e) stop("Error in 'ps_formula': ", conditionMessage(e), call. = FALSE)
    )
    all_complete <- stats::complete.cases(ps_mf) & !is.na(data[[time_var]]) & !is.na(data[[event_var]])

    if (sum(all_complete) == 0) {
      stop("No complete cases available for analysis.", call. = FALSE)
    }

    data_clean <- data[all_complete, , drop = FALSE]
  }

  # ============================================================================
  # STEP 3: Validate Data Variables (Complete Cases)
  # ============================================================================

  treatment <- data_clean[[treatment_var]]
  time <- data_clean[[time_var]]
  event <- data_clean[[event_var]]

  treatment_levels <- sort(unique(treatment))

  # For PS-based methods, require at least 2 treatment levels
  if (weight_method != "none" && length(treatment_levels) < 2) {
    stop("Treatment variable must have at least 2 levels in complete cases for propensity score weighting.",
         call. = FALSE)
  }

  if (!is.numeric(time) || any(time <= 0)) {
    stop("'time_var' must be numeric and positive.", call. = FALSE)
  }

  if (!all(event %in% c(0, 1))) {
    stop("'event_var' must be binary (0/1).", call. = FALSE)
  }

  # Validate att_group if specified
  if (weight_method == "ATT" && !att_group %in% treatment_levels) {
    stop("'att_group' '", att_group, "' not in treatment levels: ",
         paste(treatment_levels, collapse = ", "), call. = FALSE)
  }

  n_complete <- nrow(data_clean)

  # ============================================================================
  # STEP 4: Estimate Propensity Scores and Weights
  # ============================================================================

  if (weight_method == "none") {
    # No PS modeling: assign unit weights
    ps_result <- NULL
    weight_result <- NULL
    weights <- rep(1, n_complete)

  } else {
    # Estimate propensity scores
    ps_result <- estimate_ps(
      data = data_clean,
      treatment_var = treatment_var,
      ps_formula = ps_formula,
      ps_control = ps_control
    )

    # Map weight_method to estimand for estimate_weights()
    estimand <- switch(weight_method,
                       "IPW" = "ATE",
                       "OW" = "overlap",
                       "ATT" = "ATT")

    # Estimate weights
    weight_result <- estimate_weights(
      ps_result = ps_result,
      data = data_clean,
      treatment_var = treatment_var,
      estimand = estimand,
      att_group = att_group,
      trim = if (trim) "symmetric" else NULL,
      delta = if (trim) delta else NULL,
      alpha = NULL
    )

    weights <- weight_result$weights
  }

  # ============================================================================
  # STEP 5: Compute Weighted Kaplan-Meier
  # ============================================================================

  km_result <- estimate_weighted_km(
    data = data_clean,
    time_var = time_var,
    event_var = event_var,
    treatment_var = treatment_var,
    weights = weights,
    treatment_levels = treatment_levels
  )

  # ============================================================================
  # STEP 6: Prepare Output
  # ============================================================================

  result <- list(
    eval_times = km_result$eval_times,
    surv_estimates = km_result$surv_estimates,
    surv_var = km_result$surv_var,
    n_risk = km_result$n_risk,
    n_event = km_result$n_event,
    n_acc_event = km_result$n_acc_event,
    treatment_levels = treatment_levels,
    weight_method = weight_method,
    att_group = att_group,
    trim = trim,
    delta = if (trim) delta else NULL,
    n = n_complete,
    ps_result = ps_result,
    weight_result = weight_result,
    weights = weights
  )

  class(result) <- "weightedKM"
  return(result)
}


# Print Method ----------------------------------------------------------------

#' Print Method for Weighted Kaplan-Meier Estimates
#'
#' @description
#' Prints a summary of weighted Kaplan-Meier survival estimates.
#'
#' @param x An object of class "weightedKM" from \code{weightedKM()}.
#' @param print.digits Number of decimal places for printed output. Default 3.
#' @param print.rows Number of rows to print for each treatment group. Default 10.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @examples
#' \donttest{
#' data(simdata_bin)
#' result <- weightedKM(
#'   data = simdata_bin,
#'   treatment_var = "Z",
#'   time_var = "time",
#'   event_var = "event",
#'   weight_method = "none"
#' )
#' print(result)
#' }
#'
#' @export
print.weightedKM <- function(x, print.digits = 3, print.rows = 10, ...) {

  # Input validation
  if (!is.numeric(print.digits) || print.digits < 0) {
    stop("'print.digits' must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(print.rows) || print.rows < 1) {
    stop("'print.rows' must be a positive integer.", call. = FALSE)
  }

  # Header
  cat("\n")
  cat("===============================================================================\n")
  cat("Weighted Kaplan-Meier Survival Estimates\n")
  cat("===============================================================================\n\n")

  # Summary information
  cat("Weight Method:", x$weight_method, "\n")
  if (!is.null(x$att_group)) {
    cat("ATT Group:", x$att_group, "\n")
  }
  if (!is.null(x$trim) && x$trim) {
    cat("Trimming: Symmetric (delta =", x$delta, ")\n")
  }
  cat("Treatment Levels:", paste(x$treatment_levels, collapse = ", "), "\n")
  cat("Number of Complete Cases:", x$n, "\n")
  cat("Number of Unique Event Times:", length(x$eval_times), "\n")
  cat("\n")

  # Print survival estimates for each group
  eval_times <- x$eval_times
  n_times <- length(eval_times)

  for (grp in x$treatment_levels) {
    cat("-------------------------------------------------------------------------------\n")
    cat("Treatment Group:", grp, "\n")
    cat("-------------------------------------------------------------------------------\n")

    # Extract data for this group
    surv_est <- x$surv_estimates[, grp]
    se_est <- sqrt(x$surv_var[, grp])

    # Create matrix for printing
    print_matrix <- matrix(
      c(eval_times, surv_est, se_est),
      ncol = 3,
      dimnames = list(NULL, c("time", "survival", "se"))
    )

    # Determine rows to print
    rows_to_print <- min(print.rows, n_times)

    # Round and print
    print_matrix_rounded <- round(print_matrix[1:rows_to_print, , drop = FALSE], print.digits)
    print(print_matrix_rounded)

    if (n_times > print.rows) {
      cat("... (", n_times - print.rows, " more rows)\n", sep = "")
    }
    cat("\n")
  }

  cat("===============================================================================\n")
  cat("Use summary() for confidence intervals. Use plot() to visualize.\n")
  cat("===============================================================================\n\n")

  invisible(x)
}


#' Summary Method for Weighted Kaplan-Meier Estimates
#'
#' @description
#' Generates summary tables of weighted Kaplan-Meier survival or cumulative risk
#' estimates with confidence intervals for each treatment group.
#'
#' @param object An object of class "weightedKM" from \code{weightedKM()}.
#' @param type Type of estimate to summarize: "Kaplan-Meier" (survival probabilities,
#'   default) or "CR" (cumulative risk, aka. cumulative incidence = 1 - survival).
#' @param conf_type Type of confidence interval: "plain", "log", or
#'   "log-log" (default). See \code{?plot.weightedKM} for details.
#' @param conf_level Confidence level for intervals. Default 0.95.
#' @param print.digits Number of decimal places for printed output. Default 3.
#' @param print.rows Number of rows to print for each treatment group. Default 10.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list with one element per treatment group. Each element is a matrix
#'   with columns:
#'   \itemize{
#'     \item \code{time}: Evaluation time points
#'     \item \code{estimate}: Survival probability or cumulative risk
#'     \item \code{se}: Standard error
#'     \item \code{CI_lower}: Lower confidence bound
#'     \item \code{CI_upper}: Upper confidence bound
#'   }
#'   The list is returned invisibly with full precision. Printed output is
#'   rounded to \code{print.digits} decimal places and shows first \code{print.rows}
#'   rows per group.
#'
#' @details
#' This method provides tabular summaries of weighted Kaplan-Meier estimates with
#' confidence intervals. It uses the same CI calculation methods as
#' \code{plot.weightedKM()}.
#'
#' When \code{type = "CR"}, the function transforms survival estimates to
#' cumulative risk \eqn{1 - S} and calculates confidence intervals on that scale.
#'
#' The returned list contains full-precision matrices that can be used for
#' further analysis. The printed output is rounded for readability.
#'
#' @export
summary.weightedKM <- function(object,
                               type = "Kaplan-Meier",
                               conf_type = "log-log",
                               conf_level = 0.95,
                               print.digits = 3,
                               print.rows = 10,
                               ...) {
  
  # ============================================================================
  # STEP 1: Input Validation
  # ============================================================================
  
  if (!type %in% c("Kaplan-Meier", "CR")) {
    stop("'type' must be 'Kaplan-Meier' or 'CR'.", call. = FALSE)
  }
  if (!conf_type %in% c("plain", "log", "log-log")) {
    stop("'conf_type' must be 'plain', 'log', or 'log-log'.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(print.digits) || print.digits < 0) {
    stop("'print.digits' must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(print.rows) || print.rows < 1) {
    stop("'print.rows' must be a positive integer.", call. = FALSE)
  }
  
  # ============================================================================
  # STEP 2: Extract and Transform Data
  # ============================================================================
  
  eval_times <- object$eval_times
  treatment_levels <- object$treatment_levels
  
  # Transform to cumulative risk if requested
  if (type == "CR") {
    # S stores cumulative risk (1 - survival)
    S_matrix <- 1 - object$surv_estimates
    quantity_name <- "Cumulative Risk"
    S_at_zero <- 0  # CR starts at 0
  } else {
    # S stores survival
    S_matrix <- object$surv_estimates
    quantity_name <- "Survival Probability"
    S_at_zero <- 1  # Survival starts at 1
  }
  
  # Variance is the same for both S and 1-S
  var_matrix <- object$surv_var
  SE_matrix <- sqrt(var_matrix)
  
  # ============================================================================
  # STEP 2.5: Prepend time = 0 row
  # ============================================================================
  
  # Summary should start at time = 0 with S = S_at_zero
  # At time 0: variance = 0, SE = 0, CI is degenerate [S_at_zero, S_at_zero]
  
  n_cols <- ncol(S_matrix)
  
  # Prepend time = 0
  eval_times <- c(0, eval_times)
  
  # Prepend S = S_at_zero for all groups
  S_matrix <- rbind(
    matrix(S_at_zero, nrow = 1, ncol = n_cols, dimnames = list(NULL, colnames(S_matrix))),
    S_matrix
  )
  
  # Prepend variance = 0 and SE = 0 for all groups
  var_matrix <- rbind(
    matrix(0, nrow = 1, ncol = n_cols, dimnames = list(NULL, colnames(var_matrix))),
    var_matrix
  )
  SE_matrix <- rbind(
    matrix(0, nrow = 1, ncol = n_cols, dimnames = list(NULL, colnames(SE_matrix))),
    SE_matrix
  )
  
  # ============================================================================
  # STEP 3: Calculate Confidence Intervals
  # ============================================================================
  
  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  
  # Initialize CI matrices
  n_times <- length(eval_times)
  CI_lower_matrix <- matrix(NA, nrow = n_times, ncol = ncol(S_matrix))
  CI_upper_matrix <- matrix(NA, nrow = n_times, ncol = ncol(S_matrix))
  colnames(CI_lower_matrix) <- colnames(S_matrix)
  colnames(CI_upper_matrix) <- colnames(S_matrix)
  
  for (grp in treatment_levels) {
    S <- S_matrix[, grp]
    SE <- SE_matrix[, grp]
    
    if (conf_type == "plain") {
      # Plain CI: S ± z × SE
      CI_lower <- S - z_alpha * SE
      CI_upper <- S + z_alpha * SE
      
      # Truncate to [0, 1]
      CI_lower <- pmax(CI_lower, 0)
      CI_upper <- pmin(CI_upper, 1)
      
    } else if (conf_type == "log") {
      # Log transformation: CI = S^(exp(±z × SE/S))
      # Avoid division by zero or log of zero
      valid <- S > 0 & !is.na(S) & !is.na(SE)
      
      CI_lower <- rep(NA, length(S))
      CI_upper <- rep(NA, length(S))
      
      CI_lower[valid] <- S[valid] ^ exp(-z_alpha * SE[valid] / S[valid])
      CI_upper[valid] <- S[valid] ^ exp(z_alpha * SE[valid] / S[valid])
      
    } else if (conf_type == "log-log") {
      # Log-log transformation: CI = S^(exp(±z × SE / (S × |log(S)|)))
      # Avoid division by zero or log of zero/one
      valid <- S > 0 & S < 1 & !is.na(S) & !is.na(SE)
      
      CI_lower <- rep(NA, length(S))
      CI_upper <- rep(NA, length(S))
      
      denominator <- S[valid] * abs(log(S[valid]))
      CI_lower[valid] <- S[valid] ^ exp(-z_alpha * SE[valid] / denominator)
      CI_upper[valid] <- S[valid] ^ exp(z_alpha * SE[valid] / denominator)
    }
    
    # Override CI at time = 0 to degenerate interval [S_at_zero, S_at_zero]
    # At time 0, SE = 0, so CI should be degenerate regardless of transformation
    CI_lower[1] <- S_at_zero
    CI_upper[1] <- S_at_zero
    
    CI_lower_matrix[, grp] <- CI_lower
    CI_upper_matrix[, grp] <- CI_upper
  }
  
  # ============================================================================
  # STEP 4: Build Summary List
  # ============================================================================
  
  summary_list <- lapply(treatment_levels, function(grp) {
    matrix(
      c(eval_times,
        S_matrix[, grp],
        SE_matrix[, grp],
        CI_lower_matrix[, grp],
        CI_upper_matrix[, grp]),
      ncol = 5,
      dimnames = list(
        NULL,
        c("time", "estimate", "se", "CI_lower", "CI_upper")
      )
    )
  })
  names(summary_list) <- treatment_levels
  
  # ============================================================================
  # STEP 5: Print Summary
  # ============================================================================
  
  cat("\n")
  cat("===============================================================================\n")
  cat("Summary of Weighted Kaplan-Meier Estimates\n")
  cat("===============================================================================\n\n")
  cat("Type:", quantity_name, "\n")
  cat("Confidence Level:", paste0(conf_level * 100, "%"), "\n")
  cat("Confidence Type:", conf_type, "\n")
  cat("Weight Method:", object$weight_method, "\n")
  if (!is.null(object$att_group)) {
    cat("ATT Group:", object$att_group, "\n")
  }
  if (!is.null(object$trim_method)) {
    cat("Trimming:", object$trim_method, "\n")
  }
  cat("\n")
  
  for (grp in treatment_levels) {
    cat("-------------------------------------------------------------------------------\n")
    cat("Treatment Group:", grp, "\n")
    cat("-------------------------------------------------------------------------------\n")
    
    summ_grp <- summary_list[[grp]]
    n_rows <- nrow(summ_grp)
    
    # Determine rows to print
    rows_to_print <- min(print.rows, n_rows)
    
    # Round for printing
    summ_print <- round(summ_grp[1:rows_to_print, , drop = FALSE], print.digits)
    
    # Print with row numbers
    print(summ_print)
    
    if (n_rows > print.rows) {
      cat("... (", n_rows - print.rows, " more rows)\n", sep = "")
    }
    cat("\n")
  }
  
  cat("===============================================================================\n")
  cat("Note: Full-precision results returned invisibly. Use print.rows to show more.\n")
  cat("===============================================================================\n\n")
  
  # ============================================================================
  # STEP 6: Return Results Invisibly
  # ============================================================================
  
  invisible(summary_list)
}


#' Plot Method for Weighted Kaplan-Meier Estimates
#'
#' @param x An object of class "weightedKM" from \code{weightedKM()}.
#' @param type Type of curve to plot: "Kaplan-Meier" (survival probabilities,
#'   default) or "CR" (cumulative risk, aka. cumulative incidence = 1 - survival).
#' @param include_CI Logical. Include confidence interval ribbons? Default TRUE.
#' @param conf_type Type of confidence interval: "plain", "log", or
#'   "log-log" (default). See Details.
#' @param conf_level Confidence level for intervals. Default 0.95.
#' @param max_time Maximum time to display on x-axis. Default is maximum observed
#'   event time.
#' @param strata_to_plot Character vector of treatment levels to plot. Default
#'   plots all groups.
#' @param strata_colors Character vector of colors for each stratum in
#'   \code{strata_to_plot}. Must match length. Default uses ggplot2 colors.
#' @param curve_width Line width for survival curves. Default 1.
#' @param CI_alpha Transparency level for confidence interval ribbons (0-1).
#'   Default 0.3.
#' @param legend_position Position of legend: "right" or "bottom". Default "right".
#' @param legend_title Title for legend. Default "Treatment".
#' @param plot_title Main plot title. Default depends on \code{type}.
#' @param xlab X-axis label. Default "Time".
#' @param ylab Y-axis label. Default depends on \code{type}.
#' @param xlim Numeric vector of length 2 specifying x-axis limits. Default \code{c(0, max_time)}.
#' @param ylim Numeric vector of length 2 specifying y-axis limits. Default \code{c(0, 1)}.
#' @param risk_table Logical. Add risk table below the plot? Default FALSE.
#'   Requires \code{cowplot} package if TRUE.
#' @param risk_table_stats Character vector specifying statistics to display in risk table.
#'   Options: "n.risk" (number at risk), "n.acc.event" (cumulative events).
#'   Default \code{c("n.risk", "n.acc.event")}.
#' @param risk_table_height Numeric. Relative height of risk table (0-1). Default 0.25.
#' @param risk_table_breaks Numeric vector of time points for risk table columns.
#'   If NULL (default), automatically determined based on \code{max_time}.
#' @param risk_table_fontsize Numeric. Font size for risk table text. Default 3.5.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot2 object if \code{risk_table = FALSE}, or a combined plot
#'   (cowplot object) if \code{risk_table = TRUE}.
#'
#' @details
#' When \code{type = "CR"}, the function plots \eqn{1 - S(t)} representing the
#' probability of experiencing the event by time t. Variance is the same as
#' for survival, but confidence intervals are calculated on the CR scale.
#'
#' @export
plot.weightedKM <- function(x,
                            type = "Kaplan-Meier",
                            include_CI = TRUE,
                            conf_type = "log-log",
                            conf_level = 0.95,
                            max_time = NULL,
                            strata_to_plot = NULL,
                            strata_colors = NULL,
                            curve_width = 1,
                            CI_alpha = 0.3,
                            legend_position = "right",
                            legend_title = NULL,
                            plot_title = NULL,
                            xlab = "Time",
                            ylab = NULL,
                            xlim = NULL,
                            ylim = NULL,
                            risk_table = FALSE,
                            risk_table_stats = c("n.risk", "n.acc.event"),
                            risk_table_height = 0.25,
                            risk_table_breaks = NULL,
                            risk_table_fontsize = 3.5,
                            ...) {
  
  # ============================================================================
  # STEP 1: Check Package Availability
  # ============================================================================
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.", call. = FALSE)
  }
  
  if (risk_table && !requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package 'cowplot' is required for risk tables. Please install it or set risk_table = FALSE.",
         call. = FALSE)
  }
  
  # ============================================================================
  # STEP 2: Input Validation
  # ============================================================================
  
  if (!type %in% c("Kaplan-Meier", "CR")) {
    stop("'type' must be 'Kaplan-Meier' or 'CR'.", call. = FALSE)
  }
  if (!conf_type %in% c("plain", "log", "log-log")) {
    stop("'conf_type' must be 'plain', 'log', or 'log-log'.", call. = FALSE)
  }
  if (!legend_position %in% c("right", "bottom")) {
    stop("'legend_position' must be 'right' or 'bottom'.", call. = FALSE)
  }
  if (!is.logical(include_CI)) {
    stop("'include_CI' must be logical.", call. = FALSE)
  }
  if (!is.numeric(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be between 0 and 1.", call. = FALSE)
  }
  if (CI_alpha < 0 || CI_alpha > 1) {
    stop("'CI_alpha' must be between 0 and 1.", call. = FALSE)
  }
  if (!is.logical(risk_table)) {
    stop("'risk_table' must be logical.", call. = FALSE)
  }
  if (risk_table) {
    valid_stats <- c("n.risk", "n.acc.event")
    if (!all(risk_table_stats %in% valid_stats)) {
      stop("'risk_table_stats' must be a subset of: ", paste(valid_stats, collapse = ", "),
           call. = FALSE)
    }
    if (!is.numeric(risk_table_height) || risk_table_height <= 0 || risk_table_height >= 1) {
      stop("'risk_table_height' must be between 0 and 1.", call. = FALSE)
    }
  }
  
  # ============================================================================
  # STEP 3: Determine max_time
  # ============================================================================
  
  if (is.null(max_time)) {
    max_time <- max(x$eval_times)
  } else {
    if (!is.numeric(max_time) || length(max_time) != 1 || max_time <= 0) {
      stop("'max_time' must be a single positive number.", call. = FALSE)
    }
  }
  
  # Filter eval_times to max_time
  time_idx <- x$eval_times <= max_time
  if (sum(time_idx) == 0) {
    stop("No evaluation times available within max_time = ", max_time, call. = FALSE)
  }
  
  # ============================================================================
  # STEP 4: Determine Strata to Plot
  # ============================================================================
  
  all_strata <- as.character(x$treatment_levels)
  if (is.null(strata_to_plot)) {
    strata_to_plot <- all_strata
  } else {
    strata_to_plot <- as.character(strata_to_plot)
    if (!all(strata_to_plot %in% all_strata)) {
      stop("'strata_to_plot' must be a subset of treatment levels: ",
           paste(all_strata, collapse = ", "), call. = FALSE)
    }
  }
  
  # Validate strata_colors
  if (!is.null(strata_colors)) {
    if (length(strata_colors) != length(strata_to_plot)) {
      stop("'strata_colors' must have length ", length(strata_to_plot),
           " (matching strata_to_plot).", call. = FALSE)
    }
  }
  
  # ============================================================================
  # STEP 5: Extract and Transform Data
  # ============================================================================
  
  eval_times <- x$eval_times[time_idx]
  
  # Transform to cumulative risk if requested
  if (type == "CR") {
    # S stores cumulative risk (1 - survival)
    S_matrix <- 1 - x$surv_estimates[time_idx, , drop = FALSE]
    default_ylab <- "Cumulative Risk"
    default_plot_title <- "Weighted Accumulated Risk Curves"
    S_at_zero <- 0  # CR starts at 0
  } else {
    # S stores survival
    S_matrix <- x$surv_estimates[time_idx, , drop = FALSE]
    default_ylab <- "Survival Probability"
    default_plot_title <- "Weighted Kaplan-Meier Curves"
    S_at_zero <- 1  # Survival starts at 1
  }
  
  # Set ylab if not provided
  if (is.null(ylab)) ylab <- default_ylab
  
  # Variance is the same for both S and 1-S
  var_matrix <- x$surv_var[time_idx, , drop = FALSE]
  SE_matrix <- sqrt(var_matrix)
  
  # ============================================================================
  # STEP 5.5: Prepend time = 0 row
  # ============================================================================
  
  # KM curves should start at (0, 1) for survival or (0, 0) for CR
  # At time 0: variance = 0, SE = 0, CI is degenerate [S_at_zero, S_at_zero]
  
  n_times_orig <- length(eval_times)
  n_cols <- ncol(S_matrix)
  
  # Prepend time = 0
  eval_times <- c(0, eval_times)
  
  # Prepend S = S_at_zero for all groups
  S_matrix <- rbind(
    matrix(S_at_zero, nrow = 1, ncol = n_cols, dimnames = list(NULL, colnames(S_matrix))),
    S_matrix
  )
  
  # Prepend variance = 0 and SE = 0 for all groups
  var_matrix <- rbind(
    matrix(0, nrow = 1, ncol = n_cols, dimnames = list(NULL, colnames(var_matrix))),
    var_matrix
  )
  SE_matrix <- rbind(
    matrix(0, nrow = 1, ncol = n_cols, dimnames = list(NULL, colnames(SE_matrix))),
    SE_matrix
  )
  
  # ============================================================================
  # STEP 6: Calculate Confidence Intervals
  # ============================================================================
  
  z_alpha <- stats::qnorm(1 - (1 - conf_level) / 2)
  
  # Initialize CI matrices
  n_times <- length(eval_times)
  n_groups <- length(strata_to_plot)
  CI_lower_matrix <- matrix(NA, nrow = n_times, ncol = ncol(S_matrix))
  CI_upper_matrix <- matrix(NA, nrow = n_times, ncol = ncol(S_matrix))
  colnames(CI_lower_matrix) <- colnames(S_matrix)
  colnames(CI_upper_matrix) <- colnames(S_matrix)
  
  for (grp in strata_to_plot) {
    S <- S_matrix[, grp]
    SE <- SE_matrix[, grp]
    
    if (conf_type == "plain") {
      # Plain CI: S ± z × SE
      CI_lower <- S - z_alpha * SE
      CI_upper <- S + z_alpha * SE
      
      # Truncate to [0, 1]
      CI_lower <- pmax(CI_lower, 0)
      CI_upper <- pmin(CI_upper, 1)
      
    } else if (conf_type == "log") {
      # Log transformation: CI = S^(exp(±z × SE/S))
      # Avoid division by zero or log of zero
      valid <- S > 0 & !is.na(S) & !is.na(SE)
      
      CI_lower <- rep(NA, length(S))
      CI_upper <- rep(NA, length(S))
      
      CI_lower[valid] <- S[valid] ^ exp(-z_alpha * SE[valid] / S[valid])
      CI_upper[valid] <- S[valid] ^ exp(z_alpha * SE[valid] / S[valid])
      
    } else if (conf_type == "log-log") {
      # Log-log transformation: CI = S^(exp(±z × SE / (S × |log(S)|)))
      # Avoid division by zero or log of zero/one
      valid <- S > 0 & S < 1 & !is.na(S) & !is.na(SE)
      
      CI_lower <- rep(NA, length(S))
      CI_upper <- rep(NA, length(S))
      
      denominator <- S[valid] * abs(log(S[valid]))
      CI_lower[valid] <- S[valid] ^ exp(-z_alpha * SE[valid] / denominator)
      CI_upper[valid] <- S[valid] ^ exp(z_alpha * SE[valid] / denominator)
    }
    
    # Override CI at time = 0 to degenerate interval [S_at_zero, S_at_zero]
    # At time 0, SE = 0, so CI should be degenerate regardless of transformation
    CI_lower[1] <- S_at_zero
    CI_upper[1] <- S_at_zero
    
    CI_lower_matrix[, grp] <- CI_lower
    CI_upper_matrix[, grp] <- CI_upper
  }
  
  # ============================================================================
  # STEP 7: Build Plot Data Frame
  # ============================================================================
  
  plot_data <- do.call(rbind, lapply(strata_to_plot, function(grp) {
    data.frame(
      Time = eval_times,
      Estimate = S_matrix[, grp],
      CI_lower = CI_lower_matrix[, grp],
      CI_upper = CI_upper_matrix[, grp],
      Strata = grp,
      stringsAsFactors = FALSE
    )
  }))
  
  # Set factor levels for consistent ordering
  plot_data$Strata <- factor(plot_data$Strata, levels = strata_to_plot)
  
  # ============================================================================
  # STEP 8: Set Default Titles
  # ============================================================================
  
  if (is.null(plot_title)) plot_title <- default_plot_title
  if (is.null(legend_title)) legend_title <- "Treatment"
  
  # ============================================================================
  # STEP 9: Create ggplot
  # ============================================================================
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Time, y = Estimate,
                                               color = Strata, fill = Strata))
  
  # Add CI ribbon if requested
  if (include_CI) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = CI_lower, ymax = CI_upper),
      alpha = CI_alpha, color = NA
    )
  }
  
  # Add curves
  p <- p + ggplot2::geom_line(linewidth = curve_width)
  
  # Set axis limits and labels
  if (is.null(xlim)) xlim <- c(0, max_time)
  if (is.null(ylim)) ylim <- c(0, 1)
  
  p <- p + ggplot2::scale_x_continuous(limits = xlim)
  p <- p + ggplot2::scale_y_continuous(limits = ylim)
  
  p <- p + ggplot2::labs(
    title = plot_title,
    x = xlab,
    y = ylab,
    color = legend_title,
    fill = legend_title
  )
  
  # Apply custom colors if provided
  if (!is.null(strata_colors)) {
    p <- p + ggplot2::scale_color_manual(values = strata_colors)
    p <- p + ggplot2::scale_fill_manual(values = strata_colors)
  }
  
  # Apply theme and legend position
  p <- p + ggplot2::theme_minimal()
  if (legend_position == "bottom") {
    p <- p + ggplot2::theme(legend.position = "bottom",
                            plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    p <- p + ggplot2::theme(legend.position = "right",
                            plot.title = ggplot2::element_text(hjust = 0.5))
  }
  
  # ============================================================================
  # STEP 10: Add Risk Table if Requested
  # ============================================================================
  
  if (!risk_table) {
    return(p)
  }
  
  # Determine time breaks for risk table
  if (is.null(risk_table_breaks)) {
    # Create nice breaks: always start at 0, then evenly spaced
    n_breaks <- 4  # Total number of breaks (including 0)
    if (max_time <= 10) {
      # For small max_time, use intervals of 2 or 5
      interval <- ifelse(max_time <= 5, 2, 5)
    } else {
      # Round to nearest 10, 20, 30, etc.
      interval <- ceiling((max_time / (n_breaks - 1)) / 10) * 10
    }
    risk_table_breaks <- seq(0, max_time, by = interval)
    # Remove breaks beyond max_time
    risk_table_breaks <- risk_table_breaks[risk_table_breaks <= max_time]
    # Ensure we have at least 0 and max_time
    if (!max_time %in% risk_table_breaks) {
      risk_table_breaks <- c(risk_table_breaks, max_time)
    }
  }
  
  # Extract n_risk and n_acc_event at break times
  # Use right-continuous step function (standard for KM)
  # Note: Use original event times from x, not the modified eval_times (which has time=0 prepended)
  original_times <- x$eval_times[time_idx]  # Same filtering as main plot
  
  risk_table_data <- do.call(rbind, lapply(strata_to_plot, function(grp) {
    rows_list <- lapply(risk_table_breaks, function(t_break) {
      # Find the last event time <= t_break
      idx <- which(original_times <= t_break)
      if (length(idx) == 0) {
        # Before first event (including t_break = 0): use first event time as proxy
        # n_risk[1, grp] represents all at risk before first event
        n_risk_val <- x$n_risk[1, grp]
        n_acc_event_val <- 0
      } else {
        idx <- max(idx)
        n_risk_val <- x$n_risk[idx, grp]
        n_acc_event_val <- x$n_acc_event[idx, grp]
      }
      
      # Prepare row
      row_data <- data.frame(
        Strata = grp,
        Time = t_break,
        stringsAsFactors = FALSE
      )
      
      # Add requested statistics
      if ("n.risk" %in% risk_table_stats) {
        row_data$n.risk <- round(n_risk_val, 0)
      }
      if ("n.acc.event" %in% risk_table_stats) {
        row_data$n.acc.event <- round(n_acc_event_val, 0)
      }
      
      row_data
    })
    do.call(rbind, rows_list)
  }))
  
  # Create combined label for statistics
  if (all(c("n.risk", "n.acc.event") %in% risk_table_stats)) {
    risk_table_data$Label <- paste0(risk_table_data$n.risk, " (", risk_table_data$n.acc.event, ")")
  } else if ("n.risk" %in% risk_table_stats) {
    risk_table_data$Label <- as.character(risk_table_data$n.risk)
  } else {
    risk_table_data$Label <- as.character(risk_table_data$n.acc.event)
  }
  
  # Create table with 2 rows per group:
  # Row 1 (top): Group name
  # Row 2 (bottom): Statistics at each time point
  # Order matches legend: A, B, C, D from top to bottom
  
  n_groups <- length(strata_to_plot)
  
  table_plot_data <- do.call(rbind, lapply(seq_along(strata_to_plot), function(i) {
    grp <- strata_to_plot[i]
    grp_data <- risk_table_data[risk_table_data$Strata == grp, ]
    
    # Calculate Y positions (higher Y = visually on top)
    # For group i: name at higher Y, stats at lower Y
    y_name <- (n_groups - i + 1) * 2
    y_stats <- (n_groups - i + 1) * 2 - 1
    
    # Row 1: Group name (placed at leftmost position)
    name_row <- data.frame(
      Time = xlim[1],  # Leftmost position
      Y = y_name,
      Label = as.character(grp),
      Type = "name",
      Strata = grp,
      stringsAsFactors = FALSE
    )
    
    # Row 2: Statistics at each time point
    stats_rows <- data.frame(
      Time = grp_data$Time,
      Y = y_stats,
      Label = grp_data$Label,
      Type = "stats",
      Strata = grp,
      stringsAsFactors = FALSE
    )
    
    rbind(name_row, stats_rows)
  }))
  
  # Determine y-axis label based on statistics shown
  if (all(c("n.risk", "n.acc.event") %in% risk_table_stats)) {
    y_axis_label <- "At Risk (Events)"
  } else if ("n.risk" %in% risk_table_stats) {
    y_axis_label <- "At Risk"
  } else {
    y_axis_label <- "Cumulative Events"
  }
  
  # Create risk table plot with 2-row structure per group
  # All text is black regardless of strata_colors
  p_table <- ggplot2::ggplot(table_plot_data, ggplot2::aes(x = Time, y = Y)) +
    # Add group names (bold, black, left-aligned)
    ggplot2::geom_text(
      data = subset(table_plot_data, Type == "name"),
      ggplot2::aes(label = Label),
      size = risk_table_fontsize,
      fontface = "bold",
      hjust = 0,  # Left-aligned
      vjust = 0.5,
      color = "black"
    ) +
    # Add statistics (centered at each time point)
    ggplot2::geom_text(
      data = subset(table_plot_data, Type == "stats"),
      ggplot2::aes(label = Label),
      size = risk_table_fontsize,
      hjust = 0.5,  # Centered
      vjust = 0.5,
      color = "black"
    ) +
    ggplot2::scale_x_continuous(limits = xlim, breaks = risk_table_breaks) +
    ggplot2::scale_y_continuous(limits = c(0, n_groups * 2 + 1)) +
    ggplot2::labs(x = "", y = y_axis_label) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold", angle = 90, margin = ggplot2::margin(r = 10)),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  # Combine plots using cowplot
  combined_plot <- cowplot::plot_grid(
    p, p_table,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(1 - risk_table_height, risk_table_height)
  )
  
  return(combined_plot)
}