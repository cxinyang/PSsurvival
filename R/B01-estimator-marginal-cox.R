#' Marginal Cox Model Estimation with Propensity Score Weighting
#'
#' @description
#' Fits a weighted marginal Cox proportional hazards model to estimate
#' marginal hazard ratios between treatment groups using propensity score weights.

# Core Estimation Function -----------------------------------------------

#' Fit Weighted Marginal Cox Model
#'
#' @description
#' Fits a marginal Cox model `Surv(time, event) ~ treatment` with propensity
#' score weights to estimate marginal hazard ratios. Automatically handles
#' zero weights (from trimming) by subsetting data before fitting.
#'
#' @param data A data.frame containing the complete-case analysis data.
#' @param treatment_var A character string specifying the name of the treatment
#'   variable in \code{data}.
#' @param time_var A character string specifying the name of the time variable
#'   in \code{data}.
#' @param event_var A character string specifying the name of the event variable
#'   in \code{data}. Should be coded as 1 = event, 0 = censored.
#' @param weights A numeric vector of propensity score weights with length equal
#'   to nrow(data). Returned from \code{estimate_weights()}. May contain zeros
#'   for trimmed observations (these will be excluded before fitting).
#' @param treatment_levels A vector of unique treatment values (sorted). Should
#'   match the levels from \code{estimate_ps()}.
#' @param reference_level Which treatment level to use as reference in the Cox
#'   model. MANDATORY parameter.
#' @param robust Logical. Use robust (sandwich) variance estimator? Default TRUE.
#'   When TRUE, uses \code{coxph(..., robust = TRUE)}.
#' @param functionality Character string indicating purpose: "main" for main
#'   point estimation or "boot" for bootstrap. Default "main". Controls error
#'   handling behavior when groups are missing or have no events.
#'
#' @return A list containing:
#'   \item{cox_model}{The fitted coxph model object. NULL if fitting failed in bootstrap mode.}
#'   \item{hr_estimates}{Named numeric vector of log hazard ratios (coefficients).
#'     Length = n_levels - 1. Names indicate which group is compared to reference
#'     (e.g., "trtB" means group B vs reference). Contains NA for groups that are
#'     missing in data or have no events.}
#'   \item{hr_se_robust}{Named numeric vector of robust standard errors for log HRs.
#'     Same length and names as hr_estimates. NA for failed groups.}
#'   \item{reference_level}{The treatment level used as reference.}
#'   \item{treatment_levels}{Vector of all treatment levels.}
#'   \item{n_levels}{Number of treatment levels.}
#'   \item{n_per_group_original}{Named numeric vector of sample sizes per group
#'     before trimming (from original data).}
#'   \item{n_per_group_used}{Named numeric vector of sample sizes per group
#'     actually used in Cox model (after excluding zero weights).}
#'   \item{events_per_group_original}{Named numeric vector of event counts per
#'     group before trimming.}
#'   \item{events_per_group_used}{Named numeric vector of event counts per group
#'     in fitted Cox model.}
#'
#' @details
#' \strong{Functionality Modes:}
#'
#' **"main" mode (point estimation):**
#' - If reference group is missing or has no events after trimming: throws error
#' - If non-reference group is missing: sets its HR and SE to NA, continues
#' - If non-reference group has no events: sets its HR and SE to NA if coxph fails
#' - Does NOT suppress warnings/messages from coxph
#'
#' **"boot" mode (bootstrap):**
#' - If reference group is missing or has no events: returns hr_estimates with all NA, no error
#' - If non-reference group is missing: sets its HR and SE to NA
#' - Does NOT suppress warnings (this is handled in bootstrap wrapper function)
#'
#' \strong{Model Formula:}
#' Fits \code{Surv(time, event) ~ treatment} where treatment is a factor with
#' k levels. Cox regression automatically creates k-1 dummy variables relative
#' to the reference level.
#'
#' \strong{Zero Weights:}
#' coxph does not accept zero or negative weights. Observations with weight <= 0
#' are excluded before model fitting. This handles trimming automatically.
#'
#' \strong{Coefficients:}
#' Coefficients represent log hazard ratios. To get hazard ratios, use
#' \code{exp(hr_estimates)}. A positive coefficient means higher hazard
#' (worse survival) compared to reference group.
#'
#' \strong{Robust Variance:}
#' When \code{robust = TRUE}, the sandwich variance estimator accounts for
#' weighting and provides more conservative standard errors. This is recommended
#' for propensity score weighted analyses.
#'
#' @keywords internal
fit_marginal_cox <- function(data, treatment_var, time_var, event_var, weights,
                              treatment_levels, reference_level,
                              robust = TRUE, functionality = "main") {

  # Validate essential inputs (only what's not covered in U01-data-validation)
  if (length(weights) != nrow(data)) {
    stop("Length of 'weights' (", length(weights), ") must equal nrow(data) (",
         nrow(data), ").", call. = FALSE)
  }

  if (!functionality %in% c("main", "boot")) {
    stop("'functionality' must be 'main' or 'boot'. Got: ", functionality, call. = FALSE)
  }

  if (missing(reference_level) || is.null(reference_level)) {
    stop("'reference_level' is mandatory and must be specified.", call. = FALSE)
  }

  if (!reference_level %in% treatment_levels) {
    stop("reference_level '", reference_level, "' is not in treatment_levels: ",
         paste(treatment_levels, collapse = ", "), call. = FALSE)
  }

  n_levels <- length(treatment_levels)

  # Calculate sample sizes and event counts BEFORE trimming
  treatment_original <- data[[treatment_var]]
  event_original <- data[[event_var]]

  n_per_group_original <- stats::setNames(
    sapply(treatment_levels, function(lv) sum(treatment_original == lv)),
    as.character(treatment_levels)
  )

  events_per_group_original <- stats::setNames(
    sapply(treatment_levels, function(lv) sum(treatment_original == lv & event_original == 1)),
    as.character(treatment_levels)
  )

  # Subset to non-zero weights (coxph requires weights > 0)
  keep <- weights > 0
  n_used_total <- sum(keep)

  # Initialize return structure for hr_estimates (all NA by default)
  non_ref_levels <- setdiff(treatment_levels, reference_level)
  hr_names <- paste0("trt", non_ref_levels)
  hr_estimates <- stats::setNames(rep(NA_real_, length(non_ref_levels)), hr_names)
  hr_se_robust <- stats::setNames(rep(NA_real_, length(non_ref_levels)), hr_names)

  # Calculate sample sizes and events AFTER trimming
  if (n_used_total > 0) {
    data_subset <- data[keep, , drop = FALSE]
    weights_subset <- weights[keep]
    treatment_subset_raw <- data_subset[[treatment_var]]
    time_subset <- data_subset[[time_var]]
    event_subset <- data_subset[[event_var]]

    n_per_group_used <- stats::setNames(
      sapply(treatment_levels, function(lv) sum(treatment_subset_raw == lv)),
      as.character(treatment_levels)
    )

    events_per_group_used <- stats::setNames(
      sapply(treatment_levels, function(lv) sum(treatment_subset_raw == lv & event_subset == 1)),
      as.character(treatment_levels)
    )
  } else {
    n_per_group_used <- stats::setNames(rep(0, n_levels), as.character(treatment_levels))
    events_per_group_used <- stats::setNames(rep(0, n_levels), as.character(treatment_levels))
  }

  # Helper function for failed fit returns
  return_failed_fit <- function() {
    list(
      cox_model = NULL,
      hr_estimates = hr_estimates,
      hr_se_robust = hr_se_robust,
      reference_level = reference_level,
      treatment_levels = treatment_levels,
      n_levels = n_levels,
      n_per_group_original = n_per_group_original,
      n_per_group_used = n_per_group_used,
      events_per_group_original = events_per_group_original,
      events_per_group_used = events_per_group_used
    )
  }

  # Check reference group conditions
  ref_n <- n_per_group_used[as.character(reference_level)]
  ref_events <- events_per_group_used[as.character(reference_level)]

  if (ref_n == 0 || ref_events == 0) {
    if (functionality == "main") {
      stop("Reference group '", reference_level, "' has ",
           ifelse(ref_n == 0, "no observations", "no events"),
           " after trimming. Cannot fit marginal Cox model.", call. = FALSE)
    } else {
      return(return_failed_fit())
    }
  }

  # Check if there are ANY non-zero weights
  if (n_used_total == 0) {
    if (functionality == "main") {
      stop("All weights are zero. Cannot fit Cox model.", call. = FALSE)
    } else {
      return(return_failed_fit())
    }
  }

  # Prepare treatment variable as factor with reference level first
  treatment_subset <- factor(treatment_subset_raw, levels = c(
    reference_level,
    setdiff(treatment_levels, reference_level)
  ))

  # Fit marginal Cox model: Surv(time, event) ~ treatment
  surv_obj <- survival::Surv(time = time_subset, event = event_subset)
  cox_formula <- surv_obj ~ treatment_subset

  cox_model <- tryCatch({
    survival::coxph(
      formula = cox_formula,
      data = data_subset,
      weights = weights_subset,
      robust = robust,
      model = TRUE
    )
  }, error = function(e) {
    if (functionality == "main") {
      stop("Cox model fitting failed: ", conditionMessage(e), call. = FALSE)
    } else {
      return(NULL)
    }
  })

  if (is.null(cox_model)) {
    return(return_failed_fit())
  }

  # Extract coefficients and robust SEs (coefficients are log hazard ratios)
  coef_from_model <- stats::coef(cox_model)
  cox_summary <- summary(cox_model)

  if (robust) {
    if ("robust se" %in% colnames(cox_summary$coefficients)) {
      se_from_model <- cox_summary$coefficients[, "robust se"]
    } else {
      if (functionality == "main") {
        warning("Robust SE not found in Cox model summary. Using regular SE.", call. = FALSE)
      }
      se_from_model <- cox_summary$coefficients[, "se(coef)"]
    }
  } else {
    se_from_model <- cox_summary$coefficients[, "se(coef)"]
  }

  # Map coefficients to hr_estimates structure
  # coef_from_model has names like "treatment_subsetB" -> extract "B"
  estimated_levels <- gsub("^treatment_subset", "", names(coef_from_model))

  for (i in seq_along(estimated_levels)) {
    hr_name <- paste0("trt", estimated_levels[i])
    if (hr_name %in% names(hr_estimates)) {
      hr_estimates[hr_name] <- coef_from_model[i]
      hr_se_robust[hr_name] <- se_from_model[i]
    }
  }

  # Return results
  return(list(
    cox_model = cox_model,
    hr_estimates = hr_estimates,
    hr_se_robust = hr_se_robust,
    reference_level = reference_level,
    treatment_levels = treatment_levels,
    n_levels = n_levels,
    n_per_group_original = n_per_group_original,
    n_per_group_used = n_per_group_used,
    events_per_group_original = events_per_group_original,
    events_per_group_used = events_per_group_used
  ))
}
