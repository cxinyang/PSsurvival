#' Propensity Score Weighting for PSsurvival Package
#'
#' @description
#' Functions for calculating propensity score weights with optional trimming
#' for binary and multiple treatment groups.


# Internal Weight Calculation Functions ----------------------------------

#' Calculate Inverse Probability Weights (IPW)
#'
#' @description
#' Calculates inverse probability weights for binary or multiple treatment groups.
#' For ATE (Average Treatment Effect), weights are 1/e_j(X) for each treatment j.
#' For ATT (Average Treatment Effect on the Treated), weights target a specific
#' treatment group.
#'
#' @param ps_result A list returned by \code{estimate_ps()}, containing:
#'   \itemize{
#'     \item ps_matrix: Matrix of propensity scores (n x J) for multiple treatments,
#'       or NULL for binary treatment
#'     \item ps: Vector of propensity scores for observed treatment
#'     \item n_levels: Number of treatment levels
#'     \item treatment_levels: Vector of treatment level values
#'   }
#' @param data A data.frame containing the treatment variable.
#' @param treatment_var A character string specifying the name of the treatment
#'   variable in \code{data}.
#' @param estimand Character string specifying the estimand. One of "ATE" (default)
#'   or "ATT" (Average Treatment Effect on the Treated).
#' @param att_group For ATT estimation, specifies which treatment group to target.
#'   This is MANDATORY when estimand = "ATT".
#'
#' @return A numeric vector of IPW weights with length equal to nrow(data).
#'
#' @details
#' For ATE with multiple treatments, weights are:
#' \deqn{w_j(X_i) = \frac{1}{e_j(X_i)} \text{ for individual i in treatment j}}
#'
#' For ATT targeting treatment k, weights are:
#' \deqn{w_j(X_i) = \begin{cases}
#'   1 & \text{if } j = k \\
#'   \frac{e_k(X_i)}{e_j(X_i)} & \text{if } j \neq k
#' \end{cases}}
#'
#' @keywords internal
calculate_ipw_weights <- function(ps_result, data, treatment_var,
                                   estimand = "ATE", att_group = NULL) {

  # Extract treatment variable
  treatment <- data[[treatment_var]]
  n <- length(treatment)
  n_levels <- ps_result$n_levels
  treatment_levels <- ps_result$treatment_levels

  # Input validation
  if (!estimand %in% c("ATE", "ATT")) {
    stop("estimand must be either 'ATE' or 'ATT'", call. = FALSE)
  }

  if (estimand == "ATT") {
    if (is.null(att_group)) {
      stop("att_group must be specified when estimand = 'ATT'.\n",
           "  Available treatment levels: ", paste(treatment_levels, collapse = ", "),
           call. = FALSE)
    }
    if (!att_group %in% treatment_levels) {
      stop("att_group '", att_group, "' is not one of the treatment levels.\n",
           "  Available treatment levels: ", paste(treatment_levels, collapse = ", "),
           call. = FALSE)
    }
  }

  # Binary treatment case
  if (n_levels == 2) {

    if (estimand == "ATE") {
      # ATE weights: 1/P(Z=observed|X)
      weights <- 1 / ps_result$ps

    } else {  # ATT
      is_att_group <- (treatment == att_group)
      # P(Z = att_group | X): ps if in att_group, 1-ps otherwise
      prob_att_group <- ifelse(is_att_group, ps_result$ps, 1 - ps_result$ps)
      # ATT weights: 1 for att_group, prob_att_group/ps_result$ps for others
      weights <- ifelse(is_att_group, 1, prob_att_group / ps_result$ps)
    }

  } else {
    # Multiple treatment case
    ps_matrix <- ps_result$ps_matrix

    if (estimand == "ATE") {
      # ATE weights: 1/e_j(X) for observed treatment j
      weights <- 1 / ps_result$ps

    } else {  # ATT
      # Get column index for target group
      att_col <- which(treatment_levels == att_group)

      # ATT weights: 1 for att_group, e_att(X)/e_j(X) for others
      weights <- numeric(n)
      is_att_group <- (treatment == att_group)

      weights[is_att_group] <- 1
      weights[!is_att_group] <- ps_matrix[!is_att_group, att_col] / ps_result$ps[!is_att_group]
    }
  }

  return(weights)
}


#' Calculate Overlap Weights
#'
#' @description
#' Calculates overlap weights using the Li & Li (2019) formula for binary or
#' multiple treatment groups. Overlap weights target the population with the
#' most equipoise and are bounded between 0 and 1.
#'
#' @param ps_result A list returned by \code{estimate_ps()}, containing:
#'   \itemize{
#'     \item ps_matrix: Matrix of propensity scores (n x J) for multiple treatments,
#'       or NULL for binary treatment
#'     \item ps: Vector of propensity scores for observed treatment
#'     \item n_levels: Number of treatment levels
#'     \item treatment_levels: Vector of treatment level values
#'   }
#' @param data A data.frame containing the treatment variable.
#' @param treatment_var A character string specifying the name of the treatment
#'   variable in \code{data}.
#'
#' @return A numeric vector of overlap weights with length equal to nrow(data).
#'
#' @details
#' For multiple treatments, overlap weights are calculated as:
#' \deqn{w_j(X_i) = \frac{1/e_j(X_i)}{\sum_{l=1}^{J} 1/e_l(X_i)}}
#'
#' For binary treatment (J=2), this reduces to:
#' \deqn{w(X) = 1 - P(Z = observed | X)}
#' which equals the probability of receiving the opposite treatment.
#'
#' @references
#' Li, F., & Li, F. (2019). Propensity score weighting for causal inference with
#' multiple treatments. \emph{The Annals of Applied Statistics}, 13(4), 2389-2415.
#'
#' @keywords internal
calculate_overlap_weights <- function(ps_result, data, treatment_var) {

  n <- length(data[[treatment_var]])
  n_levels <- ps_result$n_levels

  # Binary treatment case
  if (n_levels == 2) {
    # Overlap weight = probability of opposite treatment
    # ps_result$ps = P(Z = observed | X)
    # 1 - ps_result$ps = P(Z = opposite | X)
    weights <- 1 - ps_result$ps

  } else {
    # Multiple treatment case
    ps_matrix <- ps_result$ps_matrix

    # Calculate overlap weights: (1/e_j) / sum(1/e_l)
    # For numerical stability: 1 / (e_j * sum(1/e_l))
    sum_inv_ps <- rowSums(1 / ps_matrix)
    weights <- 1 / (ps_result$ps * sum_inv_ps)
  }

  return(weights)
}


# Trimming Functions ------------------------------------------------------

#' Symmetric Propensity Score Trimming (Crump Extension)
#'
#' @description
#' Performs symmetric trimming based on minimum propensity score across all
#' treatment groups. Implements the Crump extension to multiple treatments as
#' described in Yoshida et al. (2019).
#'
#' @param ps_result A list returned by \code{estimate_ps()}.
#' @param data A data.frame containing the treatment variable.
#' @param treatment_var A character string specifying the name of the treatment
#'   variable in \code{data}.
#' @param delta Numeric trimming threshold in (0, 1/J] where J is the number
#'   of treatment levels. Default is NULL, which uses recommended values:
#'   0.1 for binary treatment, 0.067 for 3 groups, 1/(2*J) for J >= 4.
#'
#' @return A logical vector of length n, where TRUE indicates the observation
#'   should be kept and FALSE indicates it should be trimmed.
#'
#' @details
#' The symmetric trimming rule retains observation i if:
#' \deqn{\min_j\{e_{ji}\} \geq \delta}
#'
#' For binary treatment (J=2), this reduces to: \eqn{e(X) \in [\delta, 1-\delta]}.
#'
#' @references
#' Yoshida, K., et al. (2019). Multinomial extension of propensity score trimming
#' methods: A simulation study. \emph{American Journal of Epidemiology}, 188(3),
#' 609-616.
#'
#' @keywords internal
trim_weights_symmetric <- function(ps_result, data, treatment_var, delta = NULL) {

  n_levels <- ps_result$n_levels
  treatment <- data[[treatment_var]]
  treatment_levels <- ps_result$treatment_levels

  # Set default delta based on number of groups
  if (is.null(delta)) {
    if (n_levels == 2) {
      delta <- 0.10
    } else if (n_levels == 3) {
      delta <- 0.067
    } else {
      delta <- 1 / (2 * n_levels)
    }
    message("Using default symmetric trimming threshold delta = ",
            round(delta, 4), " for ", n_levels, " treatment groups")
  }

  # Validate delta
  max_delta <- 1 / n_levels
  if (delta <= 0 || delta > max_delta) {
    stop("delta must be in (0, ", round(max_delta, 4), "] for ", n_levels,
         " treatment groups. Got delta = ", delta, call. = FALSE)
  }

  # Binary treatment case
  if (n_levels == 2) {
    ps <- ps_result$ps_matrix[, 2]
    keep <- (ps >= delta) & (ps <= (1 - delta))

  } else {
    # Multiple treatment case
    ps_matrix <- ps_result$ps_matrix
    min_ps <- apply(ps_matrix, 1, min)
    keep <- (min_ps >= delta)
  }

  return(keep)
}


#' Asymmetric Propensity Score Trimming (Sturmer Extension)
#'
#' @description
#' Performs asymmetric (percentile-based) trimming using within-group percentile
#' thresholds. Implements the Sturmer extension to multiple treatments as
#' described in Yoshida et al. (2019).
#'
#' @param ps_result A list returned by \code{estimate_ps()}.
#' @param data A data.frame containing the treatment variable.
#' @param treatment_var A character string specifying the name of the treatment
#'   variable in \code{data}.
#' @param alpha Numeric percentile threshold in (0, 0.5). Default is NULL,
#'   which uses recommended values: 0.05 for binary treatment, 0.033 for 3 groups,
#'   1/(2*J) for J >= 4.
#'
#' @return A logical vector of length n, where TRUE indicates the observation
#'   should be kept and FALSE indicates it should be trimmed.
#'
#' @details
#' The asymmetric trimming rule retains observation i if:
#' \deqn{e_{ji} \geq F^{-1}_{e_{ji}|A_i=j}(\alpha|j) \text{ for all } j}
#'
#' where \eqn{F^{-1}_{e_{ji}|A_i=j}(\alpha|j)} is the \eqn{\alpha}-percentile of
#' propensity scores \eqn{e_{ji}} among individuals who actually received
#' treatment j.
#'
#' @references
#' Yoshida, K., et al. (2019). Multinomial extension of propensity score trimming
#' methods: A simulation study. \emph{American Journal of Epidemiology}, 188(3),
#' 609-616.
#'
#' @keywords internal
trim_weights_asymmetric <- function(ps_result, data, treatment_var, alpha = NULL) {

  n_levels <- ps_result$n_levels
  treatment <- data[[treatment_var]]
  treatment_levels <- ps_result$treatment_levels
  n <- length(treatment)

  # Set default alpha based on number of groups
  if (is.null(alpha)) {
    if (n_levels == 2) {
      alpha <- 0.05
    } else if (n_levels == 3) {
      alpha <- 0.033
    } else {
      alpha <- 1 / (2 * n_levels)
    }
    message("Using default asymmetric trimming threshold alpha = ",
            round(alpha, 4), " for ", n_levels, " treatment groups")
  }

  # Validate alpha
  if (alpha <= 0 || alpha >= 0.5) {
    stop("alpha must be in (0, 0.5). Got alpha = ", alpha, call. = FALSE)
  }

  # Binary treatment case
  if (n_levels == 2) {
    # Extract P(Z=level2|X) for all observations from ps_matrix
    ps <- ps_result$ps_matrix[, 2]

    # Separate by observed treatment group
    is_level2 <- (treatment == treatment_levels[2])
    ps0 <- ps[!is_level2]
    ps1 <- ps[is_level2]

    # Overlap bounds
    lps <- max(min(ps0), min(ps1))
    ups <- min(max(ps0), max(ps1))

    # Percentile thresholds
    alpha0 <- as.numeric(stats::quantile(ps0, probs = 1 - alpha, na.rm = TRUE))
    alpha1 <- as.numeric(stats::quantile(ps1, probs = alpha, na.rm = TRUE))

    # Trimming rule
    keep <- rep(NA, n)
    keep[!is_level2] <- ((ps0 >= alpha1) & (ps0 <= alpha0) & (ps0 >= lps) & (ps0 <= ups))
    keep[is_level2] <- ((ps1 >= alpha1) & (ps1 <= alpha0) & (ps1 >= lps) & (ps1 <= ups))

  } else {
    # Multiple treatment case
    ps_matrix <- ps_result$ps_matrix

    # Calculate threshold for each treatment group
    thresholds <- numeric(n_levels)
    for (j in 1:n_levels) {
      # Get propensity scores e_ji for treatment j, among those who received j
      in_group_j <- (treatment == treatment_levels[j])
      ps_j_in_group <- ps_matrix[in_group_j, j]

      # Calculate alpha-percentile
      thresholds[j] <- stats::quantile(ps_j_in_group, probs = alpha, na.rm = TRUE)
    }

    # For each observation, check if all propensity scores >= their thresholds
    keep <- rep(TRUE, n)
    for (j in 1:n_levels) {
      keep <- keep & (ps_matrix[, j] >= thresholds[j])
    }
  }

  return(keep)
}


# Main User-Facing Function -----------------------------------------------

#' Estimate Propensity Score Weights
#'
#' @description
#' Calculates propensity score weights for causal inference with optional trimming.
#' Supports ATE, ATT, and overlap population estimands for binary and multiple
#' treatment groups.
#'
#' @param ps_result A list returned by \code{estimate_ps()}, containing the
#'   fitted propensity score model and estimated propensity scores.
#' @param data A data.frame containing the treatment variable (same data used
#'   in \code{estimate_ps()}).
#' @param treatment_var A character string specifying the name of the treatment
#'   variable in \code{data}.
#' @param estimand Character string specifying the target population. One of:
#'   \itemize{
#'     \item "ATE": Average Treatment Effect (default). Uses IPW method.
#'     \item "ATT": Average Treatment Effect on the Treated. Uses IPW method.
#'     \item "overlap": Overlap population (Li & Li, 2019). Uses overlap weighting.
#'   }
#' @param att_group For ATT estimation, specifies which treatment group to target.
#'   This is MANDATORY when estimand = "ATT". Ignored for other estimands.
#' @param trim Character string specifying the trimming method, or NULL for no
#'   trimming (default). Options: "symmetric" (Crump extension) or "asymmetric"
#'   (Sturmer extension). Trimming is NOT supported with overlap estimand.
#' @param delta Trimming threshold for symmetric trimming in (0, 1/J], where J is
#'   the number of treatment levels. If NULL (default), uses recommended values
#'   from Yoshida et al. (2019). Ignored unless trim = "symmetric".
#' @param alpha Percentile threshold for asymmetric trimming in (0, 0.5). If NULL
#'   (default), uses recommended values from Yoshida et al. (2019). Ignored unless
#'   trim = "asymmetric".
#'
#' @return A list containing:
#'   \item{weights}{Numeric vector of weights (length = nrow(data)).}
#'   \item{trim_summary}{Data frame with trimming summary by treatment group.}
#'   \item{ess}{Named numeric vector of effective sample size by treatment group.}
#'   \item{method}{Character string: "IPW" for ATE/ATT, "overlap" for overlap.}
#'   \item{estimand}{Character string of estimand used.}
#'   \item{att_group}{Target group for ATT (NULL if not applicable).}
#'   \item{trim_method}{Character string of trimming method (NULL if no trimming).}
#'   \item{delta}{Numeric trimming threshold used for symmetric trimming (NULL if not applicable).}
#'   \item{alpha}{Numeric percentile threshold used for asymmetric trimming (NULL if not applicable).}
#'   \item{n_levels}{Number of treatment levels.}
#'   \item{treatment_levels}{Vector of treatment level values.}
#'   \item{ps_result}{PS result object (refitted after trimming if trimming was applied).}
#'
#' @details
#' \strong{Trimming Workflow:}
#' When trimming is requested, the function: (1) identifies observations to trim
#' using PS from full data, (2) re-estimates PS on trimmed data, (3) calculates
#' weights from re-estimated PS. This ensures trimming uses the original covariate
#' distribution while weights reflect the overlapping population.
#'
#' Overlap weights do not support trimming (already bounded in [0,1]).
#'
#' @examples
#' \donttest{
#' # Example 1: Overlap weighting for binary treatment
#' data(simdata_bin)
#' ps_bin <- estimate_ps(
#'   data = simdata_bin,
#'   treatment_var = "Z",
#'   ps_formula = Z ~ X1 + X2 + X3 + B1 + B2
#' )
#' weights_ow <- estimate_weights(
#'   ps_result = ps_bin,
#'   data = simdata_bin,
#'   treatment_var = "Z",
#'   estimand = "overlap"
#' )
#' summary(weights_ow$weights)
#'
#' # Example 2: ATT with multiple treatments
#' data(simdata_multi)
#' ps_multi <- estimate_ps(
#'   data = simdata_multi,
#'   treatment_var = "Z",
#'   ps_formula = Z ~ X1 + X2 + X3 + B1 + B2
#' )
#' weights_att <- estimate_weights(
#'   ps_result = ps_multi,
#'   data = simdata_multi,
#'   treatment_var = "Z",
#'   estimand = "ATT",
#'   att_group = "C"
#' )
#' summary(weights_att$weights)
#' }
#'
#' @references
#' Li, F., & Li, F. (2019). Propensity score weighting for causal inference with
#' multiple treatments. \emph{The Annals of Applied Statistics}, 13(4), 2389-2415.
#'
#' Yoshida, K., et al. (2019). Multinomial extension of propensity score trimming
#' methods: A simulation study. \emph{American Journal of Epidemiology}, 188(3),
#' 609-616.
#'
#' Crump, R. K., et al. (2009). Dealing with limited overlap in estimation of
#' average treatment effects. \emph{Biometrika}, 96(1), 187-199.
#'
#' @export
estimate_weights <- function(ps_result, data, treatment_var,
                             estimand = "ATE",
                             att_group = NULL,
                             trim = NULL,
                             delta = NULL,
                             alpha = NULL) {

  # Input validation
  if (!estimand %in% c("ATE", "ATT", "overlap")) {
    stop("estimand must be 'ATE', 'ATT', or 'overlap'", call. = FALSE)
  }

  if (!is.null(trim) && !trim %in% c("symmetric", "asymmetric")) {
    stop("trim must be NULL, 'symmetric', or 'asymmetric'", call. = FALSE)
  }

  # Overlap weights do not support trimming
  if (estimand == "overlap" && !is.null(trim)) {
    stop("Trimming is not supported with overlap weights.\n",
         "  Overlap weights have downweighted the tail populations and been bounded [0,1].\n",
         "  Use estimand = 'ATE' or 'ATT' with trim if needed.",
         call. = FALSE)
  }

  # Determine method from estimand
  method <- if (estimand == "overlap") "overlap" else "IPW"

  # Extract basic info
  treatment <- data[[treatment_var]]
  n_original <- length(treatment)
  n_levels <- ps_result$n_levels
  treatment_levels <- ps_result$treatment_levels

  # Trimming workflow
  delta_used <- NULL
  alpha_used <- NULL

  if (!is.null(trim)) {
    # Step 1: Identify observations to trim using original PS
    if (trim == "symmetric") {
      # Set default delta if not specified
      if (is.null(delta)) {
        if (n_levels == 2) {
          delta <- 0.10
        } else if (n_levels == 3) {
          delta <- 0.067
        } else {
          delta <- 1 / (2 * n_levels)
        }
      }
      delta_used <- delta
      keep <- trim_weights_symmetric(ps_result, data, treatment_var, delta)
    } else {
      # Asymmetric trimming
      # Set default alpha if not specified
      if (is.null(alpha)) {
        if (n_levels == 2) {
          alpha <- 0.05
        } else if (n_levels == 3) {
          alpha <- 0.033
        } else {
          alpha <- 1 / (2 * n_levels)
        }
      }
      alpha_used <- alpha
      keep <- trim_weights_asymmetric(ps_result, data, treatment_var, alpha)
    }

    # Step 2: Subset data
    data_trimmed <- data[keep, , drop = FALSE]

    # Step 3: Re-estimate PS on trimmed data
    ps_formula_obj <- if (n_levels == 2) {
      ps_result$ps_model$formula
    } else {
      stats::formula(ps_result$ps_model)
    }

    ps_result_final <- tryCatch(
      estimate_ps(
        data = data_trimmed,
        treatment_var = treatment_var,
        ps_formula = ps_formula_obj,
        ps_control = list(trace = FALSE)
      ),
      error = function(e) {
        stop("Propensity score model cannot be refitted after trimming; ",
             "error from estimate_ps: ", conditionMessage(e), call. = FALSE)
      }
    )

    data_final <- data_trimmed
  } else {
    # No trimming
    ps_result_final <- ps_result
    data_final <- data
    keep <- rep(TRUE, n_original)
  }

  # Calculate weights
  if (estimand == "overlap") {
    weights_trimmed <- calculate_overlap_weights(ps_result_final, data_final, treatment_var)
    att_group_out <- NULL
  } else {
    weights_trimmed <- calculate_ipw_weights(
      ps_result_final, data_final, treatment_var, estimand, att_group
    )
    att_group_out <- att_group
  }

  # Expand weights back to original data size (trimmed obs get weight 0)
  weights_full <- numeric(n_original)
  weights_full[keep] <- weights_trimmed
  weights_full[!keep] <- 0

  # Calculate trimming summary
  treatment_original <- data[[treatment_var]]
  n_per_group_original <- table(treatment_original)[as.character(treatment_levels)]
  n_per_group_trimmed <- table(treatment_original[!keep])[as.character(treatment_levels)]
  n_per_group_trimmed[is.na(n_per_group_trimmed)] <- 0

  trim_summary <- data.frame(
    treatment = treatment_levels,
    n_original = as.numeric(n_per_group_original),
    n_after_trim = as.numeric(n_per_group_original) - as.numeric(n_per_group_trimmed),
    n_trimmed = as.numeric(n_per_group_trimmed),
    pct_trimmed = as.numeric(n_per_group_trimmed) / as.numeric(n_per_group_original) * 100,
    stringsAsFactors = FALSE
  )

  # Calculate effective sample size (ESS) by treatment group
  # Only for non-trimmed observations
  ess <- numeric(n_levels)
  names(ess) <- as.character(treatment_levels)
  for (j in 1:n_levels) {
    in_group_kept <- (treatment_original == treatment_levels[j]) & keep
    w_group <- weights_full[in_group_kept]
    if (sum(w_group > 0) > 0) {
      ess[j] <- sum(w_group)^2 / sum(w_group^2)
    } else {
      ess[j] <- 0
    }
  }

  # Return results
  result <- list(
    weights = weights_full,
    trim_summary = trim_summary,
    ess = ess,
    method = method,
    estimand = estimand,
    att_group = att_group_out,
    trim_method = trim,
    delta = delta_used,
    alpha = alpha_used,
    n_levels = n_levels,
    treatment_levels = treatment_levels,
    ps_result = ps_result_final
  )

  return(result)
}
