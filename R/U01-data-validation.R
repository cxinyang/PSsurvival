#' Data Validation Functions for PSsurvival Package
#'
#' @keywords internal


# Formula Validation and Parsing -------------------------------------------

#' Validate Propensity Score Formula
#'
#' Checks that the propensity score formula is valid and extracts the
#' treatment variable name from the left-hand side.
#'
#' @param ps_formula A formula object for the propensity score model.
#'   Should be of the form `treatment ~ covariates`.
#'
#' @return A list containing:
#'   \item{formula}{The validated formula object}
#'   \item{treatment_var}{Character string of the treatment variable name}
#'
#' @keywords internal
validate_ps_formula <- function(ps_formula) {
  # Check if formula object
  if (!inherits(ps_formula, "formula")) {
    stop("'ps_formula' must be a formula object.", call. = FALSE)
  }

  # Check if formula has both LHS and RHS
  if (length(ps_formula) != 3) {
    stop("'ps_formula' must have both left-hand side (treatment) and right-hand side (covariates).",
         call. = FALSE)
  }

  # Extract treatment variable name from LHS
  treatment_var <- as.character(ps_formula[[2]])

  if (length(treatment_var) == 0 || treatment_var == "") {
    stop("Left-hand side of 'ps_formula' must specify the treatment variable.",
         call. = FALSE)
  }

  # Return validated formula and extracted variable name
  list(
    formula = ps_formula,
    treatment_var = treatment_var
  )
}


#' Validate Censoring Formula
#'
#' Checks that the censoring formula is valid, uses Surv() notation, and
#' extracts the time and event variable names.
#'
#' @param censor_formula A formula object for the censoring model.
#'   Should be of the form `Surv(time, event) ~ covariates`.
#'   Use `Surv(time, event) ~ 0` to indicate no censoring adjustment
#'   (all censoring scores set to 1).
#'
#' @return A list containing:
#'   \item{formula}{The validated formula object, or NULL if RHS is 0
#'     (no censoring adjustment)}
#'   \item{time_var}{Character string of the time variable name}
#'   \item{event_var}{Character string of the event variable name}
#'
#' @keywords internal
validate_censor_formula <- function(censor_formula) {
  # Check if formula object
  if (!inherits(censor_formula, "formula")) {
    stop("'censor_formula' must be a formula object.", call. = FALSE)
  }

  # Check if formula has both LHS and RHS
  if (length(censor_formula) != 3) {
    stop("'censor_formula' must have both left-hand side (Surv object) and right-hand side (covariates).",
         call. = FALSE)
  }

  # Extract LHS
  lhs <- censor_formula[[2]]

  # Check if LHS is a call to Surv()
  # Handle both Surv() and survival::Surv() forms
  if (!is.call(lhs)) {
    stop("Left-hand side of 'censor_formula' must use Surv() notation, e.g., Surv(time, event).",
         call. = FALSE)
  }

  # Get the function name - handle namespace-qualified calls like survival::Surv
  fn_name <- lhs[[1]]
  if (is.call(fn_name) && identical(fn_name[[1]], as.name("::"))) {
    # The package called must be "Survival"
    if (fn_name[[2]] != "survival") {
      stop("Must use package 'survival' for Surv() formula object.")
    } else {
      # Namespace-qualified: survival::Surv -> extract "Surv"
      fn_char <- as.character(fn_name[[3]])
      }
  } else {
    fn_char <- as.character(fn_name)
  }

  if (fn_char != "Surv") {
    stop("Left-hand side of 'censor_formula' must use Surv() notation, e.g., Surv(time, event).",
         call. = FALSE)
  }

  # Surv() should have at least 2 arguments: time and event
  if (length(lhs) < 3) {
    stop("Surv() must have at least two arguments: time and event.",
         call. = FALSE)
  }

  # Extract time and event variables
  # 2nd element of call is time (1st is function name "Surv")
  time_expr <- lhs[[2]]
  # 3rd element of call is event
  event_expr <- lhs[[3]]

  # Extract time variable name (should be a simple name)
  if (is.name(time_expr)) {
    time_var <- as.character(time_expr)
  } else {
    # If it's an expression, deparse it
    time_var <- deparse(time_expr)
  }

  # Extract event variable name
  # Handle the common case of I(1 - Event) or I(Event)
  if (is.call(event_expr) && as.character(event_expr[[1]]) == "I") {
    # Expression is wrapped in I()
    inner_expr <- event_expr[[2]]

    # Check if inner expression is a binary operation like 1 - Event
    if (is.call(inner_expr) && as.character(inner_expr[[1]]) %in% c("-", "+")) {
      # Binary operation: extract the variable (not the constant)
      left_part <- inner_expr[[2]]
      right_part <- inner_expr[[3]]

      if (is.name(right_part)) {
        event_var <- as.character(right_part)
      } else if (is.name(left_part)) {
        event_var <- as.character(left_part)
      } else {
        # Can't determine variable name, use full expression
        event_var <- deparse(event_expr)
      }
    } else if (is.name(inner_expr)) {
      # I(Event) - just a variable wrapped in I()
      event_var <- as.character(inner_expr)
    } else {
      # Complex expression, use full deparse
      event_var <- deparse(event_expr)
    }
  } else if (is.name(event_expr)) {
    # Simple variable name
    event_var <- as.character(event_expr)
  } else {
    # Other expression, deparse it
    event_var <- deparse(event_expr)
  }


  # Check if RHS is 0 (no censoring adjustment)
  # When RHS is 0, return formula = NULL to signal no censoring model
  rhs <- censor_formula[[3]]
  no_censoring_adjustment <- is.numeric(rhs) && length(rhs) == 1 && rhs == 0

  # Return validated formula and extracted variable names
  list(
    formula = if (no_censoring_adjustment) NULL else censor_formula,
    time_var = time_var,
    event_var = event_var
  )
}


#' Check if Variables Exist in Data
#'
#' Verifies that all variables referenced in formulas exist in the dataset.
#'
#' @param data A data.frame containing the analysis data.
#' @param treatment_var Character string of treatment variable name.
#' @param ps_formula Formula object for propensity score model.
#' @param censor_formula Formula object for censoring model.
#'
#' @return Invisible NULL if all checks pass; otherwise throws an error.
#'
#' @keywords internal
check_variables_exist <- function(data, treatment_var, ps_formula, censor_formula) {
  # Get all variable names from data
  data_vars <- names(data)

  # Check treatment variable
  if (!treatment_var %in% data_vars) {
    stop("Treatment variable '", treatment_var, "' not found in data.",
         call. = FALSE)
  }

  # Try to create model frames to check if all formula variables exist
  # This will catch any missing variables in the formulas

  # Check PS formula variables
  tryCatch(
    {
      stats::model.frame(ps_formula, data = data, na.action = stats::na.pass)
    },
    error = function(e) {
      stop("Error in 'ps_formula': ", conditionMessage(e), call. = FALSE)
    }
  )

  # Check censoring formula variables
  tryCatch(
    {
      stats::model.frame(censor_formula, data = data, na.action = stats::na.pass)
    },
    error = function(e) {
      stop("Error in 'censor_formula': ", conditionMessage(e), call. = FALSE)
    }
  )

  invisible(NULL)
}



# Variable-Specific Validation ---------------------------------------------

#' Validate Data Variables
#'
#' Performs all variable-specific validation checks on treatment, time, event,
#' and covariates. Only stops execution if critical errors are found.
#'
#' @param data A data.frame containing the analysis data.
#' @param treatment_var Character string of treatment variable name.
#' @param ps_formula Formula object for propensity score model.
#' @param censor_formula Formula object for censoring model.
#'
#' @return Invisible NULL if all checks pass; otherwise throws an error.
#'
#' @keywords internal
validate_data_variables <- function(data, treatment_var, ps_formula, censor_formula) {

  # Check treatment variable
  treatment <- data[[treatment_var]]
  treatment_levels <- unique(treatment[!is.na(treatment)])
  n_levels <- length(treatment_levels)

  if (n_levels < 2) {
    stop("Treatment variable '", treatment_var, "' must have at least 2 levels. ",
         "Found ", n_levels, " level(s).", call. = FALSE)
  }

  # Check time and event variables via model frame
  mf <- stats::model.frame(censor_formula, data = data, na.action = stats::na.pass)
  surv_obj <- stats::model.response(mf)

  # Extract time
  time <- surv_obj[, "time"]
  if (!is.numeric(time)) {
    stop("Time variable must be numeric. Found type: ", class(time)[1], call. = FALSE)
  }

  # Check for non-positive time (among non-NA)
  time_complete <- time[!is.na(time)]
  if (any(time_complete <= 0)) {
    n_nonpositive <- sum(time_complete <= 0)
    stop("Time variable must be positive. Found ", n_nonpositive,
         " non-positive value(s).", call. = FALSE)
  }

  # Extract event
  event <- surv_obj[, "status"]
  event_complete <- event[!is.na(event)]
  unique_values <- sort(unique(event_complete))
  is_binary <- length(unique_values) == 2 && all(unique_values %in% c(0, 1))

  if (!is_binary) {
    stop("Event variable must be binary (0/1). Found values: ",
         paste(unique_values, collapse = ", "), call. = FALSE)
  }

  # Check covariates - verify model matrices can be constructed
  tryCatch(
    {
      ps_mf <- stats::model.frame(ps_formula, data = data, na.action = stats::na.pass)
      ps_mm <- stats::model.matrix(ps_formula, data = ps_mf)
    },
    error = function(e) {
      stop("Error constructing model matrix for 'ps_formula': ",
           conditionMessage(e), call. = FALSE)
    }
  )

  tryCatch(
    {
      censor_mf <- stats::model.frame(censor_formula, data = data, na.action = stats::na.pass)
      censor_mm <- stats::model.matrix(censor_formula, data = censor_mf)
    },
    error = function(e) {
      stop("Error constructing model matrix for 'censor_formula': ",
           conditionMessage(e), call. = FALSE)
    }
  )

  invisible(NULL)
}


# Data Structure Validation ------------------------------------------------

#' Check Data Structure
#'
#' Validates that the data structure is appropriate for model fitting.
#' Only stops if models cannot be fit numerically.
#'
#' @param data A data.frame containing the analysis data.
#' @param treatment_var Character string of treatment variable name.
#' @param ps_formula Formula object for propensity score model.
#' @param censor_formula Formula object for censoring model.
#'
#' @return Invisible NULL if checks pass; otherwise throws an error.
#'
#' @keywords internal
check_data_structure <- function(data, treatment_var, ps_formula, censor_formula) {

  # Check if data is a data frame (or coercible)
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame.", call. = FALSE)
  }

  # Get complete cases (removing NAs across all variables in formulas)
  ps_mf <- stats::model.frame(ps_formula, data = data, na.action = stats::na.pass)
  censor_mf <- stats::model.frame(censor_formula, data = data, na.action = stats::na.pass)

  # Combine to find rows with complete data for analysis
  ps_complete <- stats::complete.cases(ps_mf)
  censor_complete <- stats::complete.cases(censor_mf)
  all_complete <- ps_complete & censor_complete

  n_complete <- sum(all_complete)

  if (n_complete == 0) {
    stop("No complete cases available for analysis after removing missing values.",
         call. = FALSE)
  }

  # Check if enough observations in each treatment group (complete cases only)
  treatment_complete <- data[[treatment_var]][all_complete]
  treatment_table <- table(treatment_complete)

  if (any(treatment_table == 0)) {
    stop("After removing missing values, at least one treatment group has no observations.",
         call. = FALSE)
  }

  # Check if sufficient observations relative to parameters in PS model
  # This ensures the PS model can be fit
  ps_mm <- stats::model.matrix(ps_formula, data = ps_mf[all_complete, , drop = FALSE])
  n_ps_params <- ncol(ps_mm)

  if (n_complete <= n_ps_params) {
    stop("Insufficient complete observations (n = ", n_complete,
         ") relative to PS model parameters (p = ", n_ps_params,
         "). Need at least n > p.", call. = FALSE)
  }

  # Check if sufficient observations in each treatment group for censoring models
  # Each treatment group needs enough observations to fit the censoring model
  censor_mm <- stats::model.matrix(censor_formula, data = censor_mf[all_complete, , drop = FALSE])
  n_censor_params <- ncol(censor_mm)

  min_group_size <- min(treatment_table)
  if (min_group_size < n_censor_params) {
    warning("Smallest treatment group has ", min_group_size, " complete observations, ",
         "but censoring model has ", n_censor_params, " parameters. ",
         "Censoring model may not be uniquely identified.", call. = FALSE)
  }

  # Check if there are any events in each treatment group (for survival models)
  surv_obj <- stats::model.response(censor_mf[all_complete, , drop = FALSE])
  events <- surv_obj[, "status"]

  for (trt_level in unique(treatment_complete)) {
    n_events_in_group <- sum(events[treatment_complete == trt_level])
    if (n_events_in_group == 0) {
      stop("Treatment group '", trt_level, "' has no observed events. ",
           "Cannot fit censoring model for this group.", call. = FALSE)
    }
  }

  invisible(NULL)
}


# Method-Specific Validation -----------------------------------------------

#' Validate Method Arguments
#'
#' Validates weight_method, censor_method, trimming parameters, time_points,
#' and control lists. Synthesizes all method-specific checks.
#'
#' @param weight_method Character string specifying weighting method.
#' @param censor_method Character string specifying censoring method.
#' @param trim_alpha Numeric, symmetric trimming threshold.
#' @param trim_q Numeric, asymmetric trimming quantile.
#' @param time_points Numeric vector or NULL.
#' @param conf_level Numeric, confidence level.
#' @param ps_control List of PS model control parameters.
#' @param censor_control List of censoring model control parameters.
#' @param bootstrap_control List of bootstrap control parameters or NULL.
#' @param max_time Numeric, maximum observed time (to check time_points range).
#'
#' @return Invisible NULL if checks pass; otherwise throws an error.
#'
#' @keywords internal
validate_method_arguments <- function(weight_method, censor_method, trim_alpha, trim_q,
                                      time_points, conf_level, ps_control, censor_control,
                                      bootstrap_control, max_time) {

  # Validate weight_method
  valid_weight_methods <- c("OW", "IPTW", "Symmetric", "Asymmetric")
  weight_method <- match.arg(weight_method, valid_weight_methods)

  # Validate trim_alpha for Symmetric method
  if (weight_method == "Symmetric") {
    if (!is.numeric(trim_alpha) || length(trim_alpha) != 1) {
      stop("'trim_alpha' must be a single numeric value.", call. = FALSE)
    }
    if (trim_alpha <= 0 || trim_alpha >= 0.5) {
      stop("'trim_alpha' must be in the range (0, 0.5). Found: ", trim_alpha, call. = FALSE)
    }
  }

  # Validate trim_q for Asymmetric method
  if (weight_method == "Asymmetric") {
    if (!is.numeric(trim_q) || length(trim_q) != 1) {
      stop("'trim_q' must be a single numeric value.", call. = FALSE)
    }
    if (trim_q <= 0 || trim_q >= 0.5) {
      stop("'trim_q' must be in the range (0, 0.5). Found: ", trim_q, call. = FALSE)
    }
  }

  # Validate censor_method
  valid_censor_methods <- c("weibull", "cox")
  censor_method <- match.arg(censor_method, valid_censor_methods)

  # Validate time_points
  if (!is.null(time_points)) {
    if (!is.numeric(time_points)) {
      stop("'time_points' must be a numeric vector or NULL.", call. = FALSE)
    }
    if (any(time_points < 0)) {
      stop("'time_points' must be non-negative.", call. = FALSE)
    }

    # Check if time_points exceed max observed time (warning only)
    if (any(time_points > max_time)) {
      warning("Some 'time_points' exceed the maximum observed time (",
              round(max_time, 2), "). Extrapolation may be unreliable.",
              call. = FALSE)
    }
  }

  # Validate conf_level
  if (!is.numeric(conf_level) || length(conf_level) != 1) {
    stop("'conf_level' must be a single numeric value.", call. = FALSE)
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be in the range (0, 1). Found: ", conf_level, call. = FALSE)
  }

  # Validate control lists
  if (!is.list(ps_control)) {
    stop("'ps_control' must be a list.", call. = FALSE)
  }
  if (!is.list(censor_control)) {
    stop("'censor_control' must be a list.", call. = FALSE)
  }

  # Validate bootstrap_control if Cox method is used
  if (censor_method == "cox" && !is.null(bootstrap_control)) {
    if (!is.list(bootstrap_control)) {
      stop("'bootstrap_control' must be a list or NULL.", call. = FALSE)
    }

    # Check bootstrap-specific parameters
    if (!is.null(bootstrap_control$n_boot)) {
      if (!is.numeric(bootstrap_control$n_boot) || bootstrap_control$n_boot < 1) {
        stop("'bootstrap_control$n_boot' must be a positive integer.", call. = FALSE)
      }
    }

    if (!is.null(bootstrap_control$parallel)) {
      if (!is.logical(bootstrap_control$parallel)) {
        stop("'bootstrap_control$parallel' must be logical (TRUE/FALSE).", call. = FALSE)
      }
    }

    if (!is.null(bootstrap_control$n_cores)) {
      if (!is.numeric(bootstrap_control$n_cores) || bootstrap_control$n_cores < 1) {
        stop("'bootstrap_control$n_cores' must be a positive integer.", call. = FALSE)
      }
    }

    if (!is.null(bootstrap_control$seed)) {
      if (!is.numeric(bootstrap_control$seed)) {
        stop("'bootstrap_control$seed' must be numeric.", call. = FALSE)
      }
    }
  }

  invisible(NULL)
}


# Umbrella Validation Function ---------------------------------------------

#' Validate All Inputs for PSsurvdiff
#'
#' Umbrella function that calls all validation functions and returns the
#' cleaned dataset with complete cases ready for model fitting.
#'
#' @param data A data.frame containing the analysis data.
#' @param ps_formula Formula object for propensity score model.
#' @param censor_formula Formula object for censoring model.
#' @param weight_method Character string specifying weighting method.
#' @param censor_method Character string specifying censoring method.
#' @param trim_alpha Numeric, symmetric trimming threshold.
#' @param trim_q Numeric, asymmetric trimming quantile.
#' @param time_points Numeric vector or NULL.
#' @param conf_level Numeric, confidence level.
#' @param ps_control List of PS model control parameters.
#' @param censor_control List of censoring model control parameters.
#' @param bootstrap_control List of bootstrap control parameters or NULL.
#'
#' @return A list containing:
#'   \item{data_clean}{Data frame with complete cases only}
#'   \item{treatment_var}{Character string of treatment variable name}
#'   \item{time_var}{Character string of time variable name (possibly an expression)}
#'   \item{event_var}{Character string of event variable name (possibly an expression)}
#'   \item{censor_formula}{The validated censoring formula, or NULL if no censoring
#'     adjustment was requested (i.e., original formula had RHS of 0)}
#'   \item{n_complete}{Integer, number of complete cases used in analysis}
#'
#' @keywords internal
validate_PSsurvdiff_inputs <- function(data, ps_formula, censor_formula,
                                       weight_method, censor_method,
                                       trim_alpha, trim_q, time_points, conf_level,
                                       ps_control, censor_control, bootstrap_control) {

  # Step 1: Validate formulas and extract variable names
  ps_validated <- validate_ps_formula(ps_formula)
  censor_validated <- validate_censor_formula(censor_formula)

  treatment_var <- ps_validated$treatment_var
  time_var <- censor_validated$time_var
  event_var <- censor_validated$event_var

  # Get the validated censoring formula (NULL if RHS was 0, meaning no censoring adjustment)
  censor_formula_validated <- censor_validated$formula

  # Create a temp formula for validation when no censoring adjustment is requested

  # This temp formula only checks time/event variables exist and have valid values
  # It uses ~1 (intercept only) since there are no covariates to check
  if (is.null(censor_formula_validated)) {
    censor_formula_for_validation <- stats::as.formula(
      paste0("Surv(", time_var, ", ", event_var, ") ~ 1")
    )
  } else {
    censor_formula_for_validation <- censor_formula_validated
  }

  # Step 2: Check if variables exist in data
  check_variables_exist(data, treatment_var, ps_formula, censor_formula_for_validation)

  # Step 3: Create complete-case dataset
  # Get model frames for both formulas to identify complete cases
  ps_mf <- stats::model.frame(ps_formula, data = data, na.action = stats::na.pass)
  censor_mf <- stats::model.frame(censor_formula_for_validation, data = data, na.action = stats::na.pass)

  # Identify rows with complete data across both formulas
  ps_complete <- stats::complete.cases(ps_mf)
  censor_complete <- stats::complete.cases(censor_mf)
  all_complete <- ps_complete & censor_complete

  # Extract complete cases from original data
  data_clean <- data[all_complete, , drop = FALSE]
  n_complete <- nrow(data_clean)

  # Check if any complete cases exist
  if (n_complete == 0) {
    stop("No complete cases available for analysis after removing missing values.",
         call. = FALSE)
  }

  # Step 4: Validate complete-case data and extract max_time
  # Validate variable types and values using complete-case data
  # Use censor_formula_for_validation to ensure proper validation even when no censoring adjustment
  validate_data_variables(data_clean, treatment_var, ps_formula, censor_formula_for_validation)

  # Validate data structure (sample sizes, events, etc.) using complete-case data
  check_data_structure(data_clean, treatment_var, ps_formula, censor_formula_for_validation)

  # Extract max_time from complete-case data
  censor_mf_clean <- stats::model.frame(censor_formula_for_validation, data = data_clean, na.action = stats::na.pass)
  surv_obj_clean <- stats::model.response(censor_mf_clean)
  max_time <- max(surv_obj_clean[, "time"])

  # Step 5: Validate method-specific arguments
  validate_method_arguments(weight_method, censor_method, trim_alpha, trim_q,
                           time_points, conf_level, ps_control, censor_control,
                           bootstrap_control, max_time)

  # Return cleaned data and variable names
  # censor_formula is NULL when no censoring adjustment requested (RHS was 0)
  list(
    data_clean = data_clean,
    treatment_var = treatment_var,
    time_var = time_var,
    event_var = event_var,
    censor_formula = censor_formula_validated,
    n_complete = n_complete
  )
}