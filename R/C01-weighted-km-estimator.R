#' Weighted Kaplan-Meier Estimation with Classic Greenwood Variance
#'
#' @description
#' Core estimation function for computing weighted Kaplan-Meier survival curves
#' with propensity score weights. Handles multiple treatment groups simultaneously
#' and uses the classic weighted Greenwood formula for variance estimation.

# Core Estimation Function -----------------------------------------------

#' Estimate Weighted Kaplan-Meier Curves for All Treatment Groups
#'
#' @description
#' Computes weighted Kaplan-Meier survival estimates and variances for all
#' treatment groups using propensity score weights. Uses classic weighted
#' Greenwood formula: \eqn{Var[S(t)] = [S(t)]^2 \sum (D_l / (R_l (R_l - D_l)))}.
#'
#' @param data A data.frame containing the complete-case analysis data.
#' @param time_var A character string specifying the name of the time variable.
#' @param event_var A character string specifying the name of the event variable.
#'   Should be coded as 1 = event, 0 = censored.
#' @param treatment_var A character string specifying the name of the treatment
#'   variable in \code{data}.
#' @param weights A numeric vector of propensity score weights with length equal
#'   to nrow(data). Each observation has one weight corresponding to its observed
#'   treatment group. For ATE: \eqn{w_i = 1/e_j(X_i)} where j is observed treatment.
#' @param treatment_levels A vector of unique treatment values (sorted). Should
#'   match the levels from \code{estimate_ps()}.
#'
#' @return A list containing:
#'   \item{eval_times}{Numeric vector of all unique event times where survival
#'     is estimated.}
#'   \item{surv_estimates}{Matrix [n_times x n_groups] of survival estimates.
#'     Column names are treatment levels.}
#'   \item{surv_var}{Matrix [n_times x n_groups] of variances for survival.}
#'   \item{n_risk}{Matrix [n_times x n_groups] of weighted number at risk (R).}
#'   \item{n_event}{Matrix [n_times x n_groups] of weighted number of events (D).}
#'   \item{n_acc_event}{Matrix [n_times x n_groups] of cumulative weighted events up to each time.}
#'   \item{treatment_levels}{Treatment levels (column names for matrices).}
#'   \item{n_levels}{Number of treatment groups.}
#'
#' @details
#' **Weighted Kaplan-Meier Formula:**
#'
#' For treatment group j, at each event time \eqn{t_l}:
#' \deqn{R_l = \sum_{i: T_i \ge t_l, Z_i = j} w_{i,j}}
#' \deqn{D_l = \sum_{i: T_i = t_l, \delta_i = 1, Z_i = j} w_{i,j}}
#' \deqn{\hat{S}^w_j(t) = \prod_{t_l \le t} \left(1 - \frac{D_l}{R_l}\right)}
#'
#' where \eqn{R_l} is the weighted number at risk and \eqn{D_l} is the weighted number of events.
#' Ties between events and censorings are handled using the Breslow method.
#'
#' **Classic Weighted Greenwood Variance:**
#'
#' \deqn{Var[\hat{S}^w_j(t)] = [\hat{S}^w_j(t)]^2 \sum_{t_l \le t} \frac{D_l}{R_l (R_l - D_l)}}
#'
#' This is the standard weighted extension of Greenwood's formula. When all
#' weights equal 1, reduces to classical Greenwood's formula.
#'
#' **Weight Structure:**
#'
#' The weight vector has length nrow(data). Each observation i in treatment group j
#' has weight \eqn{w_i} based on its propensity score for group j. For ATE estimation,
#' \eqn{w_i = 1/e_j(X_i)}. When computing weighted KM for group j, only observations
#' with \eqn{Z_i = j} and their corresponding weights are used.
#'
#' **Handling Edge Cases:**
#'
#' - If weighted at-risk count \eqn{R_l = 0} at time t, survival remains constant
#'   after t (last observation censored).
#' - Variance is undefined when \eqn{R_l - D_l \le 0}; set to NA for that time point.
#'
#' @keywords internal
estimate_weighted_km <- function(data,
                                  time_var,
                                  event_var,
                                  treatment_var,
                                  weights,
                                  treatment_levels) {

  # ============================================================================
  # Input Validation
  # ============================================================================

  if (length(weights) != nrow(data)) {
    stop("Length of 'weights' (", length(weights), ") must equal nrow(data) (",
         nrow(data), ").", call. = FALSE)
  }

  n_levels <- length(treatment_levels)
  treatment <- data[[treatment_var]]
  time <- data[[time_var]]
  event <- data[[event_var]]

  # ============================================================================
  # Determine Evaluation Times (all unique event times)
  # ============================================================================

  # Use all unique event times across all groups
  eval_times <- sort(unique(time[event == 1]))
  if (length(eval_times) == 0) {
    stop("No events observed in data. Cannot estimate survival curves.",
         call. = FALSE)
  }

  n_times <- length(eval_times)

  # ============================================================================
  # Initialize Output Matrices
  # ============================================================================

  surv_estimates <- matrix(NA, nrow = n_times, ncol = n_levels)
  surv_var <- matrix(NA, nrow = n_times, ncol = n_levels)
  n_risk <- matrix(NA, nrow = n_times, ncol = n_levels)
  n_event <- matrix(NA, nrow = n_times, ncol = n_levels)

  colnames(surv_estimates) <- as.character(treatment_levels)
  colnames(surv_var) <- as.character(treatment_levels)
  colnames(n_risk) <- as.character(treatment_levels)
  colnames(n_event) <- as.character(treatment_levels)

  # ============================================================================
  # Compute Weighted KM for Each Treatment Group
  # ============================================================================

  for (j in 1:n_levels) {
    grp <- treatment_levels[j]

    # Subset to current treatment group
    in_group <- (treatment == grp)
    n_in_group <- sum(in_group)

    if (n_in_group == 0) {
      # No observations in this group (e.g., all trimmed)
      # Leave as NA
      next
    }

    time_grp <- time[in_group]
    event_grp <- event[in_group]
    weights_grp <- weights[in_group]

    # Compute weighted KM using product-limit formula
    # Track cumulative survival and variance components
    surv_current <- 1  # S(t) starts at 1
    var_sum <- 0       # Sum for Greenwood variance

    for (t_idx in 1:n_times) {
      t <- eval_times[t_idx]

      # Weighted number at risk (R_l): sum of weights for individuals still at risk
      R_w <- sum(weights_grp[time_grp >= t])
      n_risk[t_idx, j] <- R_w

      # Weighted number of events (D_l) at exactly time t
      D_w <- sum(weights_grp[time_grp == t & event_grp == 1])
      n_event[t_idx, j] <- D_w

      # Update survival estimate if there are events at time t
      if (D_w > 0 && R_w > 0) {
        # Survival decrement at time t
        surv_current <- surv_current * (1 - D_w / R_w)

        # Add Greenwood variance component: D_w / (R_w * (R_w - D_w))
        # Check denominator to avoid division by zero or negative values
        if (R_w > D_w) {
          var_sum <- var_sum + D_w / (R_w * (R_w - D_w))
        } else {
          # R_w = D_w means all at-risk individuals had event at time t
          # Variance becomes undefined; set to NA
          var_sum <- NA
        }
      }

      # Store results for this time point
      surv_estimates[t_idx, j] <- surv_current

      # Greenwood variance: [S(t)]^2 * sum(variance components)
      if (!is.na(var_sum)) {
        surv_var[t_idx, j] <- (surv_current^2) * var_sum
      } else {
        surv_var[t_idx, j] <- NA
      }
    }
  }

  # ============================================================================
  # Compute Cumulative Events
  # ============================================================================

  # Cumulative weighted events up to each time point
  n_acc_event <- apply(n_event, 2, cumsum)
  colnames(n_acc_event) <- colnames(n_event)

  # ============================================================================
  # Return Results
  # ============================================================================

  return(list(
    eval_times = eval_times,
    surv_estimates = surv_estimates,
    surv_var = surv_var,
    n_risk = n_risk,
    n_event = n_event,
    n_acc_event = n_acc_event,
    treatment_levels = treatment_levels,
    n_levels = n_levels
  ))
}
