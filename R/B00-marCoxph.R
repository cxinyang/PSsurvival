#' Marginal Cox Model with Propensity Score Weighting
#'
#' @description
#' Main user interface for estimating marginal hazard ratios using propensity
#' score weighting. Supports binary and multiple treatment groups with various
#' weighting schemes (ATE, ATT, overlap) and optional trimming. Variance can be
#' estimated via bootstrap or robust sandwich estimator.
#'
#' @param data Data frame containing treatment, survival outcome, and covariates.
#' @param ps_formula Formula for propensity score model: \code{treatment ~ covariates}.
#' @param time_var Character string specifying the time-to-event variable name.
#' @param event_var Character string specifying the event indicator variable name.
#'   Should be coded as 1=event, 0=censored.
#' @param reference_level Treatment level to use as reference in Cox model. MANDATORY.
#'   Must be one of the treatment levels.
#' @param estimand Target estimand: "ATE" (average treatment effect), "ATT" (average
#'   treatment effect on the treated), or "overlap" (overlap weighting). Default "ATE".
#' @param att_group Target group for ATT. Required if \code{estimand = "ATT"}.
#' @param trim Trimming method: "symmetric" or "asymmetric". Default NULL (no trimming).
#' @param delta Threshold for symmetric trimming (e.g., 0.1). Required if \code{trim = "symmetric"}.
#' @param alpha Percentile for asymmetric trimming (e.g., 0.05). Required if \code{trim = "asymmetric"}.
#' @param variance_method Variance estimation method: "bootstrap" (default) or "robust".
#'   "bootstrap" resamples the entire analysis pipeline. "robust" uses the sandwich
#'   variance estimator from \code{coxph()} without bootstrap.
#' @param boot_level Bootstrap sampling level: "full" (default) or "strata".
#'   "full" resamples from entire dataset (standard for observational studies). "strata"
#'   resamples within each treatment group preserving group sizes (useful when treatment assignment 
#'   follows a stratified or fixed-ratio design). Only used if \code{variance_method = "bootstrap"}.
#' @param B Number of bootstrap iterations. Default 100. Used only if \code{variance_method = "bootstrap"}.
#' @param parallel Logical. Use parallel bootstrap computation? Default FALSE.
#' @param mc.cores Number of cores for parallel bootstrap. Default 2.
#' @param seed Random seed for bootstrap reproducibility. Default NULL.
#' @param ps_control Control parameters for propensity score model. Default \code{list()}.
#' @param robust Logical. Use robust (sandwich) variance in Cox model fitting?
#'   Default TRUE. When TRUE, \code{coxph()} is called with \code{robust = TRUE}.
#'
#' @return Object of class "marCoxph" containing:
#'   \item{coxph_fitted}{Fitted \code{coxph} model object.}
#'   \item{logHR_est}{Named vector of estimated log hazard ratios. Names are formatted
#'     as "treatment_var:level" (e.g., "Z:B" for treatment Z, level B vs reference).}
#'   \item{logHR_se_robust}{Named vector of robust standard errors from \code{coxph}.}
#'   \item{logHR_se_bootstrap}{Named vector of bootstrap standard errors. NULL if
#'     \code{variance_method = "robust"}.}
#'   \item{n_coxph_fitted}{Named vector of sample sizes per treatment group used in
#'     Cox model fitting (after trimming).}
#'   \item{events_coxph_fitted}{Named vector of event counts per treatment group used
#'     in Cox model fitting (after trimming).}
#'   \item{variance_method}{Variance method used: "bootstrap-full", "bootstrap-strata",
#'     or "robust".}
#'   \item{estimand}{Target estimand used.}
#'   \item{att_group}{Target group for ATT (NULL if not applicable).}
#'   \item{trim_method}{Trimming method (NULL if no trimming).}
#'   \item{delta}{Symmetric trimming threshold (NULL if not applicable).}
#'   \item{alpha}{Asymmetric trimming threshold (NULL if not applicable).}
#'   \item{treatment_var}{Name of treatment variable.}
#'   \item{treatment_levels}{Sorted unique treatment values.}
#'   \item{reference_level}{Reference level used in Cox model.}
#'   \item{n_levels}{Number of treatment groups.}
#'   \item{n}{Number of complete cases used in analysis.}
#'   \item{ps_result}{Propensity score estimation results.}
#'   \item{weight_result}{Weight estimation results.}
#'   \item{boot_result}{Bootstrap results (NULL if \code{variance_method = "robust"}).
#'     Contains: boot_samples, boot_allocation, n_success_by_group, B.}
#'
#' @details
#' **Analysis Workflow:**
#' 1. Extract treatment variable from \code{ps_formula}.
#' 2. Estimate propensity scores using multinomial logistic regression (or logistic
#'    for binary treatment).
#' 3. Calculate propensity score weights based on \code{estimand} and optional \code{trim}.
#' 4. Fit marginal Cox model \code{Surv(time, event) ~ treatment} with weights.
#' 5. Estimate variance via bootstrap (resampling full pipeline) or robust sandwich
#'    estimator.
#'
#' **Variance Estimation:**
#' - \code{bootstrap}: Resamples data (full or stratified), re-estimates PS and weights,
#'   re-fits Cox model. Provides bootstrap SE for log hazard ratios.
#' - \code{robust}: Uses robust sandwich variance from \code{coxph()} directly. No
#'   bootstrap performed (faster but may be less accurate with extreme weights).
#'
#' **Trimming:**
#' - Symmetric: Crump extension for multiple treatments (Yoshida et al., 2019).
#' - Asymmetric: Sturmer extension for multiple treatments (Yoshida et al., 2019).
#' - Not supported with overlap weights (already bounded [0,1]).
#'
#' @examples
#' \donttest{
#' # Example 1: Binary treatment with overlap weighting
#' data(simdata_bin)
#' result1 <- marCoxph(
#'   data = simdata_bin,
#'   ps_formula = Z ~ X1 + X2 + X3 + B1 + B2,
#'   time_var = "time",
#'   event_var = "event",
#'   reference_level = "A",
#'   estimand = "overlap"
#' )
#' summary(result1)
#'
#' # Example 2: Multiple treatments with ATT and robust variance
#' data(simdata_multi)
#' result2 <- marCoxph(
#'   data = simdata_multi,
#'   ps_formula = Z ~ X1 + X2 + X3 + B1 + B2,
#'   time_var = "time",
#'   event_var = "event",
#'   reference_level = "C",
#'   estimand = "ATT",
#'   att_group = "C",
#'   variance_method = "robust"
#' )
#' summary(result2)
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
#' @export
marCoxph <- function(data,
                     ps_formula,
                     time_var,
                     event_var,
                     reference_level,
                     estimand = "ATE",
                     att_group = NULL,
                     trim = NULL,
                     delta = NULL,
                     alpha = NULL,
                     variance_method = "bootstrap",
                     boot_level = "full",
                     B = 100,
                     parallel = FALSE,
                     mc.cores = 2,
                     seed = NULL,
                     ps_control = list(),
                     robust = TRUE) {

  # ============================================================================
  # STEP 1: Input Validation and Data Cleaning
  # ============================================================================

  # Basic parameter validation
  if (missing(reference_level) || is.null(reference_level)) {
    stop("'reference_level' is mandatory and must be specified.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }
  if (!variance_method %in% c("bootstrap", "robust")) {
    stop("'variance_method' must be 'bootstrap' or 'robust'", call. = FALSE)
  }
  if (!estimand %in% c("ATE", "ATT", "overlap")) {
    stop("'estimand' must be 'ATE', 'ATT', or 'overlap'", call. = FALSE)
  }
  if (estimand == "overlap" && !is.null(trim)) {
    stop("Trimming is not supported with overlap weights.\n",
         "  Use estimand = 'ATE' or 'ATT' with trim if needed.", call. = FALSE)
  }
  if (!is.null(trim) && !trim %in% c("symmetric", "asymmetric")) {
    stop("'trim' must be NULL, 'symmetric', or 'asymmetric'", call. = FALSE)
  }
  if (!is.list(ps_control)) {
    stop("'ps_control' must be a list.", call. = FALSE)
  }

  # Validate and extract treatment variable from ps_formula
  ps_validated <- validate_ps_formula(ps_formula)
  treatment_var <- ps_validated$treatment_var

  # Check time_var and event_var
  if (!is.character(time_var) || length(time_var) != 1 || !time_var %in% names(data)) {
    stop("'time_var' must be a single character string and exist in data.", call. = FALSE)
  }
  if (!is.character(event_var) || length(event_var) != 1 || !event_var %in% names(data)) {
    stop("'event_var' must be a single character string and exist in data.", call. = FALSE)
  }

  # Create complete-case dataset
  ps_mf <- tryCatch(
    stats::model.frame(ps_formula, data = data, na.action = stats::na.pass),
    error = function(e) stop("Error in 'ps_formula': ", conditionMessage(e), call. = FALSE)
  )
  all_complete <- stats::complete.cases(ps_mf) & !is.na(data[[time_var]]) & !is.na(data[[event_var]])

  if (sum(all_complete) == 0) {
    stop("No complete cases available for analysis.", call. = FALSE)
  }

  data_clean <- data[all_complete, , drop = FALSE]

  # Validate data variables (complete cases)
  treatment <- data_clean[[treatment_var]]
  time <- data_clean[[time_var]]
  event <- data_clean[[event_var]]

  treatment_levels_unique <- sort(unique(treatment))
  if (length(treatment_levels_unique) < 2) {
    stop("Treatment variable must have at least 2 levels in complete cases.", call. = FALSE)
  }
  if (!reference_level %in% treatment_levels_unique) {
    stop("'reference_level' '", reference_level, "' not in treatment levels: ",
         paste(treatment_levels_unique, collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(time) || any(time <= 0)) {
    stop("'time_var' must be numeric and positive.", call. = FALSE)
  }

  if (!all(event %in% c(0, 1))) {
    stop("'event_var' must be binary (0/1).", call. = FALSE)
  }

  # Check events per treatment group
  for (trt_level in treatment_levels_unique) {
    if (sum(event[treatment == trt_level]) == 0) {
      stop("Treatment group '", trt_level, "' has no observed events.", call. = FALSE)
    }
  }

  # Validate bootstrap parameters if needed
  if (variance_method == "bootstrap") {
    if (!is.numeric(B) || B < 1) stop("'B' must be a positive integer.", call. = FALSE)
    if (!is.logical(parallel)) stop("'parallel' must be logical.", call. = FALSE)
    if (!is.numeric(mc.cores) || mc.cores < 1) stop("'mc.cores' must be positive.", call. = FALSE)
    if (!boot_level %in% c("full", "strata")) stop("'boot_level' must be 'full' or 'strata'.", call. = FALSE)
  }

  n_complete <- nrow(data_clean)

  # ============================================================================
  # STEP 2: Statistical Analysis (using data_clean)
  # ============================================================================

  # Estimate propensity scores
  ps_result <- estimate_ps(
    data = data_clean,
    treatment_var = treatment_var,
    ps_formula = ps_formula,
    ps_control = ps_control
  )

  # Extract treatment levels (should match treatment_levels_unique from validation)
  treatment_levels <- ps_result$treatment_levels

  # Estimate weights
  weight_result <- estimate_weights(
    ps_result = ps_result,
    data = data_clean,
    treatment_var = treatment_var,
    estimand = estimand,
    att_group = att_group,
    trim = trim,
    delta = delta,
    alpha = alpha
  )

  # Fit marginal Cox model
  cox_fit <- fit_marginal_cox(
    data = data_clean,
    treatment_var = treatment_var,
    time_var = time_var,
    event_var = event_var,
    weights = weight_result$weights,
    treatment_levels = treatment_levels,
    reference_level = reference_level,
    robust = robust,
    functionality = "main"
  )

  # Format logHR names: "treatment_var:level"
  non_ref_levels <- setdiff(treatment_levels, reference_level)
  logHR_names <- paste0(treatment_var, ":", non_ref_levels)
  names(cox_fit$hr_estimates) <- logHR_names
  names(cox_fit$hr_se_robust) <- logHR_names

  # ============================================================================
  # STEP 3: Variance Estimation
  # ============================================================================

  boot_result <- NULL
  logHR_se_bootstrap <- NULL

  if (variance_method == "bootstrap") {
    boot_result_raw <- var_marginalcox_bootstrap(
      data = data_clean,
      treatment_var = treatment_var,
      time_var = time_var,
      event_var = event_var,
      ps_formula = ps_formula,
      treatment_levels = treatment_levels,
      reference_level = reference_level,
      estimand = estimand,
      att_group = att_group,
      trim = trim,
      delta = delta,
      alpha = alpha,
      boot_level = boot_level,
      B = B,
      robust = robust,
      parallel = parallel,
      mc.cores = mc.cores,
      seed = seed
    )

    # Calculate bootstrap SE
    hr_boot_matrix <- do.call(rbind, boot_result_raw$boot_samples)
    logHR_se_bootstrap <- apply(hr_boot_matrix, 2, function(x) stats::sd(x, na.rm = TRUE))
    names(logHR_se_bootstrap) <- logHR_names

    # Store essential bootstrap results (exclude boot_indices_matrix)
    boot_result <- list(
      boot_samples = boot_result_raw$boot_samples,
      boot_allocation = boot_result_raw$boot_allocation,
      n_success_by_group = boot_result_raw$n_success_by_group,
      B = boot_result_raw$B
    )

    variance_method_full <- paste0("bootstrap-", boot_level)
  } else {
    variance_method_full <- "robust"
  }

  # Format output
  output <- list(
    coxph_fitted = cox_fit$cox_model,
    logHR_est = cox_fit$hr_estimates,
    logHR_se_robust = cox_fit$hr_se_robust,
    logHR_se_bootstrap = logHR_se_bootstrap,
    n_coxph_fitted = cox_fit$n_per_group_used,
    events_coxph_fitted = cox_fit$events_per_group_used,
    variance_method = variance_method_full,
    estimand = estimand,
    att_group = weight_result$att_group,
    trim_method = weight_result$trim_method,
    delta = weight_result$delta,
    alpha = weight_result$alpha,
    treatment_var = treatment_var,
    treatment_levels = treatment_levels,
    reference_level = reference_level,
    n_levels = cox_fit$n_levels,
    n = n_complete,
    ps_result = ps_result,
    weight_result = weight_result,
    boot_result = boot_result
  )

  class(output) <- "marCoxph"
  return(output)
}


#' Print Method for marCoxph Objects
#'
#' @param x A \code{marCoxph} object.
#' @param max.len Maximum number of treatment comparisons to print. Default 10.
#' @param round.digits Number of digits for rounding displayed values. Default 4.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @export
print.marCoxph <- function(x, max.len = 10, round.digits = 4, ...) {
  cat("\n=== Marginal Cox Proportional Hazards Model ===\n\n")
  cat("Treatment variable:", x$treatment_var, "\n")
  cat("Treatment levels:", paste(x$treatment_levels, collapse = ", "), "\n")
  cat("Reference level:", x$reference_level, "\n")
  cat("Number of groups:", x$n_levels, "\n")
  cat("Sample size:", x$n, "complete cases\n")
  cat("Estimand:", x$estimand, "\n")
  if (!is.null(x$att_group)) {
    cat("ATT target group:", x$att_group, "\n")
  }
  cat("Variance method:", x$variance_method, "\n")
  if (!is.null(x$boot_result)) {
    cat("Bootstrap iterations:", x$boot_result$B, "\n")
  }

  if (!is.null(x$trim_method)) {
    cat("Trimming:", x$trim_method)
    if (!is.null(x$delta)) cat(" (delta = ", round(x$delta, 4), ")", sep = "")
    if (!is.null(x$alpha)) cat(" (alpha = ", round(x$alpha, 4), ")", sep = "")
    cat("\n")
  } else {
    cat("Trimming: None\n")
  }

  cat("\nSample sizes used in Cox model:\n")
  n_df <- data.frame(
    treatment = names(x$n_coxph_fitted),
    n = as.numeric(x$n_coxph_fitted),
    events = as.numeric(x$events_coxph_fitted)
  )
  print(n_df, row.names = FALSE)

  # Determine number of comparisons to show
  n_comparisons <- length(x$logHR_est)
  n_show <- min(n_comparisons, max.len)
  show_note <- n_comparisons > max.len

  cat("\nLog hazard ratio estimates:\n")
  est_df <- data.frame(
    treatment = names(x$logHR_est)[1:n_show],
    logHR = round(x$logHR_est[1:n_show], round.digits),
    SE_robust = round(x$logHR_se_robust[1:n_show], round.digits)
  )
  if (!is.null(x$logHR_se_bootstrap)) {
    est_df$SE_bootstrap <- round(x$logHR_se_bootstrap[1:n_show], round.digits)
  }
  print(est_df, row.names = FALSE)
  if (show_note) cat("  ... (", n_comparisons - n_show, " more comparisons not shown)\n", sep = "")

  cat("\nNote: Use summary() for confidence intervals and exponentiated hazard ratios.\n")

  invisible(x)
}


#' Summary Method for marCoxph Objects
#'
#' @param object A \code{marCoxph} object.
#' @param conf_level Confidence level for intervals. Default 0.95.
#' @param round.digits Number of digits for rounding displayed values. Default 4.
#'   Only used if \code{style = "prints"}.
#' @param style Output style: "prints" (print formatted tables) or "returns"
#'   (return vectors). Default "prints".
#' @param ... Additional arguments (ignored).
#'
#' @return If \code{style = "prints"}, returns invisibly. If \code{style = "returns"},
#'   returns a list with:
#'   \item{logHR}{Named vector of log hazard ratio estimates.}
#'   \item{logHR_CI_lower}{Named vector of lower CI bounds on log scale.}
#'   \item{logHR_CI_upper}{Named vector of upper CI bounds on log scale.}
#'   \item{SE}{Named vector of standard errors on log scale (from variance_method).}
#'   \item{HR}{Named vector of hazard ratio estimates (original scale).}
#'   \item{HR_CI_lower}{Named vector of lower CI bounds on original scale.}
#'   \item{HR_CI_upper}{Named vector of upper CI bounds on original scale.}
#'   \item{variance_method}{Variance method used.}
#'   \item{conf_level}{Confidence level used.}
#'   \item{n_per_group}{Named vector of sample sizes per group in Cox model.}
#'   \item{events_per_group}{Named vector of event counts per group in Cox model.}
#'
#' @details
#' Confidence intervals are Wald-type intervals calculated as:
#' \itemize{
#'   \item Log scale: logHR ± z_crit * SE
#'   \item Original scale: exp(logHR ± z_crit * SE)
#' }
#'
#' The SE used depends on \code{variance_method} from the original \code{marCoxph} call:
#' \itemize{
#'   \item "robust": Uses \code{logHR_se_robust} from sandwich estimator.
#'   \item "bootstrap-full" or "bootstrap-strata": Uses \code{logHR_se_bootstrap}.
#' }
#'
#' @export
summary.marCoxph <- function(object, conf_level = 0.95, round.digits = 4,
                             style = "prints", ...) {
  if (!style %in% c("prints", "returns")) {
    stop("Invalid 'style' argument. Must be 'prints' or 'returns'.", call. = FALSE)
  }

  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  ci_pct <- round(conf_level * 100)

  # Determine which SE to use based on variance_method
  if (grepl("^bootstrap", object$variance_method)) {
    se_use <- object$logHR_se_bootstrap
    if (is.null(se_use)) {
      stop("Bootstrap SE not available. This should not happen.", call. = FALSE)
    }
  } else {  # robust
    se_use <- object$logHR_se_robust
  }

  # Calculate CIs on log scale
  logHR_lower <- object$logHR_est - z_crit * se_use
  logHR_upper <- object$logHR_est + z_crit * se_use

  # Calculate estimates and CIs on original scale (HR)
  HR <- exp(object$logHR_est)
  HR_lower <- exp(logHR_lower)
  HR_upper <- exp(logHR_upper)

  if (style == "returns") {
    return(list(
      logHR = object$logHR_est,
      logHR_CI_lower = logHR_lower,
      logHR_CI_upper = logHR_upper,
      SE = se_use,
      HR = HR,
      HR_CI_lower = HR_lower,
      HR_CI_upper = HR_upper,
      variance_method = object$variance_method,
      conf_level = conf_level,
      n_per_group = object$n_coxph_fitted,
      events_per_group = object$events_coxph_fitted
    ))
  }

  # style == "prints"
  cat("\n=== Marginal Cox Model Summary ===\n\n")
  cat("Treatment variable:", object$treatment_var, "\n")
  cat("Treatment levels:", paste(object$treatment_levels, collapse = ", "), "\n")
  cat("Reference group:", object$reference_level, "\n")
  cat("Number of groups:", object$n_levels, "\n")
  cat("Sample size:", object$n, "complete cases\n")
  cat("Estimand:", object$estimand, "\n")
  if (!is.null(object$att_group)) {
    cat("ATT target group:", object$att_group, "\n")
  }
  cat("Variance method:", object$variance_method, "\n")
  if (!is.null(object$boot_result)) {
    cat("Bootstrap iterations:", object$boot_result$B, "\n")
  }
  if (!is.null(object$trim_method)) {
    cat("Trimming:", object$trim_method)
    if (!is.null(object$delta)) cat(" (delta = ", round(object$delta, 4), ")", sep = "")
    if (!is.null(object$alpha)) cat(" (alpha = ", round(object$alpha, 4), ")", sep = "")
    cat("\n")
  }
  cat("Confidence level:", conf_level, "\n\n")

  # Sample sizes per group
  cat("--- Sample Sizes in Cox Model ---\n\n")
  n_df <- data.frame(
    treatment = names(object$n_coxph_fitted),
    n = as.numeric(object$n_coxph_fitted),
    events = as.numeric(object$events_coxph_fitted)
  )
  print(n_df, row.names = FALSE)
  cat("\n")

  # Log HR estimates with CIs
  cat("--- Log Hazard Ratios (log scale) ---\n\n")
  logHR_df <- data.frame(
    treatment = names(object$logHR_est),
    logHR = round(object$logHR_est, round.digits),
    SE = round(se_use, round.digits),
    CI.lower = round(logHR_lower, round.digits),
    CI.upper = round(logHR_upper, round.digits),
    row.names = NULL
  )
  colnames(logHR_df)[4:5] <- paste0(ci_pct, "%CI.", c("lower", "upper"))
  print(logHR_df, row.names = FALSE)
  cat("\n")

  # HR estimates with CIs (original scale)
  cat("--- Hazard Ratios (original scale) ---\n\n")
  HR_df <- data.frame(
    treatment = names(object$logHR_est),
    HR = round(HR, round.digits),
    CI.lower = round(HR_lower, round.digits),
    CI.upper = round(HR_upper, round.digits),
    row.names = NULL
  )
  colnames(HR_df)[3:4] <- paste0(ci_pct, "%CI.", c("lower", "upper"))
  print(HR_df, row.names = FALSE)
  cat("\n")

  invisible(object)
}
