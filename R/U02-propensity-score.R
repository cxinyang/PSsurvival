#' Propensity Score Estimation for PSsurvival Package
#'
#' @description
#' Functions for estimating propensity scores for binary and multiple treatment groups.


# Propensity Score Estimation ---------------------------------------------

#' Estimate Propensity Scores
#'
#' Fits a propensity score model and extracts propensity scores for binary or
#' multiple treatment groups. For binary treatments, uses binomial logistic
#' regression. For multiple treatments (>2 levels), uses multinomial logistic
#' regression to estimate generalized propensity scores.
#'
#' @param data A data.frame containing the analysis data (typically the cleaned
#'   data with complete cases).
#' @param treatment_var A character string specifying the name of the treatment
#'   variable in \code{data}. Can be numeric, character, or factor with any
#'   coding (e.g., 0/1, 1/2, "Control"/"Treated"). Function assumes treatment
#'   has been validated for 2 or more levels.
#' @param ps_formula A formula object for the propensity score model, of the
#'   form \code{treatment ~ covariates}.
#' @param ps_control An optional list of control parameters to pass to the
#'   model fitting function (\code{glm} for binary treatment or
#'   \code{nnet::multinom} for multiple treatments). Default is an empty list.
#'
#' @return A list with the following components:
#'   \item{ps_model}{The fitted propensity score model object (class \code{glm}
#'     for binary treatment or \code{multinom} for multiple treatments).}
#'   \item{ps}{A numeric vector of propensity scores representing the probability
#'     of receiving the actual treatment each individual received. Length equals
#'     the number of rows in \code{data}.}
#'   \item{ps_matrix}{A numeric matrix of dimension n × K where n is the number
#'     of observations and K is the number of treatment levels. Each row contains
#'     the predicted probabilities for all treatment levels. Column names correspond
#'     to treatment levels.}
#'   \item{n_levels}{An integer indicating the number of treatment levels.}
#'   \item{treatment_levels}{A vector of unique treatment values sorted by
#'     \code{sort()}: numerically for numeric, alphabetically for character,
#'     by factor level order for factor.}
#'
#' @details
#' \strong{Propensity Score Definition:}
#' Returns P(Z = observed | X) for each individual, not P(Z=1|X) for all (as in
#' Rosenbaum & Rubin 1983). This definition enables direct use in IPW and extends
#' naturally to multiple treatments.
#'
#' \strong{Binary Treatments (2 levels):}
#' Fits binomial logistic regression via \code{glm()}. Treatment is factorized
#' with levels sorted by \code{sort()}: numerically for numeric, alphabetically
#' for character, by factor level order for factor. Returns P(Z = observed | X).
#'
#' \strong{Multiple Treatments (>2 levels):}
#' Fits multinomial logistic regression via \code{nnet::multinom()}. Returns
#' P(Z = observed | X) for each individual from the generalized PS matrix.
#'
#' \strong{Control Parameters (\code{ps_control}):}
#' \itemize{
#'   \item Binary: \code{glm.control()} parameters (default: \code{epsilon=1e-08, maxit=25})
#'   \item Multiple: \code{multinom()} parameters (default: \code{MaxNWts=10000, maxit=100, trace=FALSE})
#' }
#'
#' @examples
#' \donttest{
#' # Example 1: Binary treatment
#' data(simdata_bin)
#' ps_bin <- estimate_ps(
#'   data = simdata_bin,
#'   treatment_var = "Z",
#'   ps_formula = Z ~ X1 + X2 + X3 + B1 + B2
#' )
#' summary(ps_bin$ps)
#' table(simdata_bin$Z)
#'
#' # Example 2: Multiple treatments
#' data(simdata_multi)
#' ps_multi <- estimate_ps(
#'   data = simdata_multi,
#'   treatment_var = "Z",
#'   ps_formula = Z ~ X1 + X2 + X3 + B1 + B2
#' )
#' head(ps_multi$ps_matrix)
#' }
#'
#' @export
estimate_ps <- function(data, treatment_var, ps_formula, ps_control = list()) {

  # Extract treatment variable
  treatment <- data[[treatment_var]]

  # Get treatment levels
  treatment_levels <- sort(unique(treatment))
  n_levels <- length(treatment_levels)

  # Check if treatment has at least 2 levels (should already be validated, but double-check)
  if (n_levels < 2) {
    stop("Treatment variable must have at least 2 levels. Found ", n_levels, " level(s).",
         call. = FALSE)
  }

  # Binary treatment case (2 levels)
  if (n_levels == 2) {

    # Fit binomial logistic regression
    # Merge ps_control with default glm.control() settings
    default_control <- list(epsilon = 1e-08, maxit = 25, trace = FALSE)
    glm_control <- utils::modifyList(default_control, ps_control)

    # For binary treatment, need to ensure treatment is a factor for glm
    # Character and non-0/1 numeric treatments need to be converted
    data_glm <- data
    data_glm[[treatment_var]] <- factor(data_glm[[treatment_var]], levels = treatment_levels)

    ps_model <- stats::glm(
      formula = ps_formula,
      data = data_glm,
      family = stats::binomial(link = "logit"),
      control = do.call(stats::glm.control, glm_control)
    )

    # Extract propensity scores
    # glm fitted.values give P(Y = second factor level)
    # treatment_levels[2] is the second level after sort(unique(treatment))
    ps_fitted <- ps_model$fitted.values

    # Create ps_matrix [n x 2] for consistency with multiple treatment case
    ps_matrix <- cbind(1 - ps_fitted, ps_fitted)
    colnames(ps_matrix) <- as.character(treatment_levels)

    # Return P(Treatment = observed | X) for each individual
    second_level_indicator <- as.numeric(treatment == treatment_levels[2])
    ps <- ifelse(second_level_indicator == 1, ps_fitted, 1 - ps_fitted)

    # Return results
    return(list(
      ps_model = ps_model,
      ps = ps,
      ps_matrix = ps_matrix,
      n_levels = n_levels,
      treatment_levels = treatment_levels
    ))

  } else {
    # Multiple treatment case (>2 levels)

    # Check if nnet package is available
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("Package 'nnet' is required for multiple treatment propensity score estimation. ",
           "Please install it with: install.packages('nnet')",
           call. = FALSE)
    }

    # Fit multinomial logistic regression
    # Merge ps_control with default multinom() settings
    default_control <- list(MaxNWts = 10000, maxit = 100, trace = FALSE)
    multinom_control <- utils::modifyList(default_control, ps_control)

    # Need to convert treatment to factor for multinom
    data_multinom <- data
    data_multinom[[treatment_var]] <- factor(data_multinom[[treatment_var]])

    ps_model <- nnet::multinom(
      formula = ps_formula,
      data = data_multinom,
      MaxNWts = multinom_control$MaxNWts,
      maxit = multinom_control$maxit,
      trace = multinom_control$trace
    )

    # Extract predicted probabilities for all treatment levels
    # This returns a matrix of dimension n × K
    ps_matrix <- stats::predict(ps_model, newdata = data_multinom, type = "probs")

    # If only 3 levels, multinom sometimes returns a vector instead of matrix
    # for one of the dimensions, so ensure it's always a matrix
    if (!is.matrix(ps_matrix)) {
      ps_matrix <- matrix(ps_matrix, ncol = n_levels, byrow = FALSE)
      colnames(ps_matrix) <- levels(data_multinom[[treatment_var]])
    }

    # Extract the propensity score for the actual treatment received
    # For each individual i with observed treatment t, extract P(Treatment_i = t | X_i)
    # from ps_matrix. Use match() to find column indices, then matrix indexing.
    treatment_indices <- match(treatment, treatment_levels)
    ps <- ps_matrix[cbind(seq_along(treatment_indices), treatment_indices)]

    # Return results
    return(list(
      ps_model = ps_model,
      ps = ps,
      ps_matrix = ps_matrix,
      n_levels = n_levels,
      treatment_levels = treatment_levels
    ))
  }
}
