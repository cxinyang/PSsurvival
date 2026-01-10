#' Censoring Score Estimation
#'
#' @description
#' Estimate censoring scores P(C >= T | X) using Weibull or Cox models fit
#' separately within each treatment group.

# Weibull Censoring Score Estimation -----------------------------------------

#' Estimate Censoring Scores Using Weibull Regression
#'
#' @param data Data frame.
#' @param time_var Name of time variable.
#' @param treatment_var Name of treatment variable.
#' @param formula Censoring model formula. Use \code{Surv(time, censor_indicator) ~ X1 + X2}
#'   where \code{censor_indicator = 1} indicates censoring. If event is coded
#'   canonically (event=1, censored=0), use \code{I(1-event)}. Otherwise, use
#'   the appropriate transformation. Treatment is automatically removed if included.
#' @param control Control parameters for \code{survreg()}. Default
#'   \code{list(maxiter = 350)}.
#'
#' @return List with class "censoring_score_weibull":
#'   \item{censoring_models}{Fitted \code{survreg} objects by treatment level.}
#'   \item{censoring_scores}{P(C >= T_i | Z_i, X_i) for observed treatment.}
#'   \item{censoring_matrix}{(n x J) matrix of P(C >= T_i | Z=j, X_i).}
#'   \item{n_levels}{Number of treatment levels.}
#'   \item{treatment_levels}{Sorted treatment values.}
#'   \item{model_type}{"weibull".}
#'   \item{parameters}{Transformed parameters (\code{theta}, \code{gamma}, \code{coef},
#'     \code{vcov}, \code{scale}) by treatment level.}
#'   \item{linear_predictors_matrix}{(n x J) matrix of linear predictors.}
#'
#' @details
#' Fits Weibull models within each treatment group. Censoring scores computed as:
#' \deqn{K_c^{(j)}(t, X) = \exp(-\exp(X'\theta_j) \cdot t^{\gamma_j})}
#' where \eqn{\theta_j = -\beta_j/\sigma_j}, \eqn{\gamma_j = 1/\sigma_j}.
#'
#' @keywords internal
estimate_censoring_score_weibull <- function(data, time_var, treatment_var,
                                             formula,
                                             control = list(maxiter = 350)) {

  # Extract treatment variable and get levels
  treatment <- data[[treatment_var]]
  treatment_levels <- sort(unique(treatment))
  n_levels <- length(treatment_levels)
  n <- nrow(data)

  # Handle no censoring adjustment case (formula = NULL)
  # When formula is NULL, all censoring scores are 1 (no IPCW adjustment)
  if (is.null(formula)) {
    censoring_matrix <- matrix(1, nrow = n, ncol = n_levels)
    colnames(censoring_matrix) <- as.character(treatment_levels)

    result <- list(
      censoring_models = NULL,
      censoring_scores = rep(1, n),
      censoring_matrix = censoring_matrix,
      n_levels = n_levels,
      treatment_levels = treatment_levels,
      model_type = "weibull",
      parameters = NULL,
      linear_predictors_matrix = NULL,
      no_censoring_adjustment = TRUE
    )

    class(result) <- "censoring_score_weibull"
    return(result)
  }

  # Check if treatment is in formula and remove if present
  formula_vars <- all.vars(formula)
  if (treatment_var %in% formula_vars) {
    message("Treatment variable '", treatment_var, "' found in formula. ",
            "Removing it as censoring models are fit within each treatment group separately.")
    # Remove treatment from RHS
    formula_rhs <- as.character(formula)[3]
    formula_rhs_clean <- gsub(paste0("\\s*\\+?\\s*", treatment_var, "\\s*\\+?\\s*"), " ", formula_rhs)
    formula_rhs_clean <- gsub("^\\s*\\+\\s*|\\s*\\+\\s*$", "", formula_rhs_clean)  # Clean leading/trailing +
    formula <- stats::as.formula(paste(as.character(formula)[2], "~", formula_rhs_clean))
  }

  # Initialize storage
  censoring_models <- vector("list", n_levels)
  names(censoring_models) <- as.character(treatment_levels)

  parameters <- vector("list", n_levels)
  names(parameters) <- as.character(treatment_levels)

  censoring_matrix <- matrix(NA, nrow = n, ncol = n_levels)
  colnames(censoring_matrix) <- as.character(treatment_levels)

  linear_predictors_matrix <- matrix(NA, nrow = n, ncol = n_levels)
  colnames(linear_predictors_matrix) <- as.character(treatment_levels)

  # Extract observed time for all individuals
  time_vec <- data[[time_var]]

  # Extract RHS formula for design matrix construction (remove Surv() LHS)
  rhs_formula <- stats::as.formula(paste("~", as.character(formula)[3]))

  # Fit models within each treatment group
  for (j in 1:n_levels) {
    level_j <- treatment_levels[j]

    # Subset data for this treatment group
    data_j <- data[treatment == level_j, , drop = FALSE]

    # Fit Weibull model
    model_j <- survival::survreg(
      formula = formula,
      data = data_j,
      dist = "weibull",
      control = control
    )

    # Store model
    censoring_models[[as.character(level_j)]] <- model_j

    # Extract parameters for this group
    scale_j <- model_j$scale
    coef_j <- stats::coef(model_j)
    theta_j <- -coef_j / scale_j
    gamma_j <- 1 / scale_j
    vcov_j <- stats::vcov(model_j)

    # Store parameters
    parameters[[as.character(level_j)]] <- list(
      theta = theta_j,
      gamma = gamma_j,
      coef = coef_j,
      vcov = vcov_j,
      scale = scale_j
    )

    # Construct design matrix for ALL individuals (not just group j)
    # This allows counterfactual censoring score calculation
    X_all <- stats::model.matrix(rhs_formula, data = data)

    # Calculate linear predictors for all individuals under treatment j
    linear_pred_j <- c(X_all %*% theta_j)
    linear_predictors_matrix[, j] <- linear_pred_j

    # Calculate censoring scores for all individuals at their observed times
    # K_c^{(j)}(T_i, X_i) = exp(-exp(X_i' theta_j) * T_i^gamma_j)
    censoring_matrix[, j] <- exp(-exp(linear_pred_j) * time_vec^gamma_j)
  }

  # Extract censoring scores for observed treatment
  treatment_indices <- match(treatment, treatment_levels)
  censoring_scores <- censoring_matrix[cbind(1:n, treatment_indices)]

  # Return results
  result <- list(
    censoring_models = censoring_models,
    censoring_scores = censoring_scores,
    censoring_matrix = censoring_matrix,
    n_levels = n_levels,
    treatment_levels = treatment_levels,
    model_type = "weibull",
    parameters = parameters,
    linear_predictors_matrix = linear_predictors_matrix
  )

  class(result) <- "censoring_score_weibull"
  return(result)
}


# Cox Censoring Score Estimation ---------------------------------------------

#' Estimate Censoring Scores Using Cox Regression
#'
#' @param data Data frame.
#' @param time_var Name of time variable.
#' @param treatment_var Name of treatment variable.
#' @param formula Censoring model formula. Use \code{Surv(time, censor_indicator) ~ X1 + X2}
#'   where \code{censor_indicator = 1} indicates censoring. If event is coded
#'   canonically (event=1, censored=0), use \code{I(1-event)}. Otherwise, use
#'   the appropriate transformation. Treatment is automatically removed if included.
#' @param control Control parameters for \code{coxph()}. Default \code{list()}.
#' @param ties Tie handling method. Default "efron".
#'
#' @return List with class "censoring_score_cox":
#'   \item{censoring_models}{Fitted \code{coxph} objects by treatment level.}
#'   \item{censoring_scores}{P(C >= T_i | Z_i, X_i) for observed treatment.}
#'   \item{censoring_matrix}{(n x J) matrix of P(C >= T_i | Z=j, X_i).}
#'   \item{n_levels}{Number of treatment levels.}
#'   \item{treatment_levels}{Sorted treatment values.}
#'   \item{model_type}{"cox".}
#'   \item{baseline_hazards}{Baseline cumulative hazards by treatment level.}
#'   \item{coef_list}{Coefficient vectors by treatment level.}
#'   \item{vcov_list}{Variance-covariance matrices by treatment level.}
#'   \item{linear_predictors_matrix}{(n x J) matrix of linear predictors.}
#'
#' @details
#' Fits Cox models within each treatment group. Censoring scores computed as:
#' \deqn{K_c^{(j)}(t, X) = \exp(-H_0^{(j)}(t) \cdot \exp(\beta_j' X))}
#' where \eqn{H_0^{(j)}(t)} is cumulative baseline hazard. Baseline hazards
#' evaluated at nearest time point for each individual.
#'
#' @keywords internal
estimate_censoring_score_cox <- function(data, time_var, treatment_var,
                                         formula,
                                         control = list(), ties = "efron") {

  # Extract treatment variable and get levels
  treatment <- data[[treatment_var]]
  treatment_levels <- sort(unique(treatment))
  n_levels <- length(treatment_levels)
  n <- nrow(data)

  # Handle no censoring adjustment case (formula = NULL)
  # When formula is NULL, all censoring scores are 1 (no IPCW adjustment)
  if (is.null(formula)) {
    censoring_matrix <- matrix(1, nrow = n, ncol = n_levels)
    colnames(censoring_matrix) <- as.character(treatment_levels)

    result <- list(
      censoring_models = NULL,
      censoring_scores = rep(1, n),
      censoring_matrix = censoring_matrix,
      n_levels = n_levels,
      treatment_levels = treatment_levels,
      model_type = "cox",
      baseline_hazards = NULL,
      coef_list = NULL,
      vcov_list = NULL,
      linear_predictors_matrix = NULL,
      no_censoring_adjustment = TRUE
    )

    class(result) <- "censoring_score_cox"
    return(result)
  }

  # Check if treatment is in formula and remove if present
  formula_vars <- all.vars(formula)
  if (treatment_var %in% formula_vars) {
    message("Treatment variable '", treatment_var, "' found in formula. ",
            "Removing it as censoring models are fit within each treatment group separately.")
    # Remove treatment from RHS
    formula_rhs <- as.character(formula)[3]
    formula_rhs_clean <- gsub(paste0("\\s*\\+?\\s*", treatment_var, "\\s*\\+?\\s*"), " ", formula_rhs)
    formula_rhs_clean <- gsub("^\\s*\\+\\s*|\\s*\\+\\s*$", "", formula_rhs_clean)  # Clean leading/trailing +
    formula <- stats::as.formula(paste(as.character(formula)[2], "~", formula_rhs_clean))
  }

  # Initialize storage
  censoring_models <- vector("list", n_levels)
  names(censoring_models) <- as.character(treatment_levels)

  baseline_hazards <- vector("list", n_levels)
  names(baseline_hazards) <- as.character(treatment_levels)

  coef_list <- vector("list", n_levels)
  names(coef_list) <- as.character(treatment_levels)

  vcov_list <- vector("list", n_levels)
  names(vcov_list) <- as.character(treatment_levels)

  censoring_matrix <- matrix(NA, nrow = n, ncol = n_levels)
  colnames(censoring_matrix) <- as.character(treatment_levels)

  linear_predictors_matrix <- matrix(NA, nrow = n, ncol = n_levels)
  colnames(linear_predictors_matrix) <- as.character(treatment_levels)

  # Extract observed time for all individuals
  time_vec <- data[[time_var]]

  # Extract RHS formula for design matrix construction (remove Surv() LHS)
  rhs_formula <- stats::as.formula(paste("~", as.character(formula)[3]))

  # Fit models within each treatment group
  for (j in 1:n_levels) {
    level_j <- treatment_levels[j]

    # Subset data for this treatment group
    data_j <- data[treatment == level_j, , drop = FALSE]

    # Fit Cox model
    # Note: model=TRUE stores the model frame, needed for basehaz() to work correctly
    model_j <- survival::coxph(
      formula = formula,
      data = data_j,
      ties = ties,
      control = control,
      model = TRUE
    )

    # Store model
    censoring_models[[as.character(level_j)]] <- model_j

    # Extract coefficients and vcov
    coef_j <- stats::coef(model_j)
    vcov_j <- stats::vcov(model_j)

    coef_list[[as.character(level_j)]] <- coef_j
    vcov_list[[as.character(level_j)]] <- vcov_j

    # Extract baseline hazard with error handling
    # Note: Use centered=TRUE to match archived CensorScoreFun.R
    # The baseline hazard is centered at mean covariate values
    baseline_hazard_j <- tryCatch(
      {
        survival::basehaz(model_j, centered = TRUE)
      },
      error = function(e) {
        stop("Failed to extract baseline hazard for treatment group '", level_j, "'. ",
             "Error: ", conditionMessage(e), ". ",
             "This may occur if there are insufficient events in this group.",
             call. = FALSE)
      }
    )

    # Store baseline hazard
    baseline_hazards[[as.character(level_j)]] <- baseline_hazard_j

    # Construct design matrix for ALL individuals
    # Note: Cox models don't have intercepts, so we need to remove it from the design matrix
    X_all <- stats::model.matrix(rhs_formula, data = data)
    # Remove intercept column if present (Cox models don't use intercepts)
    if ("(Intercept)" %in% colnames(X_all)) {
      X_all <- X_all[, -1, drop = FALSE]
    }

    # Calculate linear predictors for all individuals under treatment j
    # NOTE: When using centered=TRUE in basehaz(), we must use centered linear predictors
    # Centering: LP_centered = X * beta - X_mean * beta, where X_mean is from training data (group j)
    linear_pred_j_uncentered <- c(X_all %*% coef_j)
    centering_constant <- sum(model_j$means * coef_j)
    linear_pred_j <- linear_pred_j_uncentered - centering_constant
    linear_predictors_matrix[, j] <- linear_pred_j

    # Calculate censoring scores for all individuals at their observed times
    # K_c^{(j)}(T_i, X_i) = exp(-H_0^{(j)}(T_i) * exp(beta_j' (X_i - X_mean_j)))

    # For each individual's time T_i, find closest time in baseline hazard
    H0_at_Ti <- sapply(time_vec, function(t) {
      idx <- which.min(abs(baseline_hazard_j$time - t))[1]
      baseline_hazard_j$hazard[idx]
    })

    # Calculate censoring scores
    censoring_matrix[, j] <- exp(-H0_at_Ti * exp(linear_pred_j))
  }

  # Extract censoring scores for observed treatment
  treatment_indices <- match(treatment, treatment_levels)
  censoring_scores <- censoring_matrix[cbind(1:n, treatment_indices)]

  # Return results
  result <- list(
    censoring_models = censoring_models,
    censoring_scores = censoring_scores,
    censoring_matrix = censoring_matrix,
    n_levels = n_levels,
    treatment_levels = treatment_levels,
    model_type = "cox",
    baseline_hazards = baseline_hazards,
    coef_list = coef_list,
    vcov_list = vcov_list,
    linear_predictors_matrix = linear_predictors_matrix
  )

  class(result) <- "censoring_score_cox"
  return(result)
}
