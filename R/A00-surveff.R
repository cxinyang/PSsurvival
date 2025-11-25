# Declare non-standard evaluation (NSE) variables used in ggplot2::aes()
# These column names are evaluated within data frame context in plot.surveff()
utils::globalVariables(c("Time", "Estimate", "Strata", "CI_lower", "CI_upper"))

#' Survival Effect Estimation with Propensity Score Weighting
#'
#' @description
#' Main user interface for estimating counterfactual survival functions and
#' treatment effects using propensity score weighting and inverse probability
#' of censoring weighting. Supports binary and multiple treatment groups with
#' various weighting schemes (ATE, ATT, overlap) and optional trimming.
#'
#' @param data Data frame containing treatment, outcome, and covariates.
#' @param ps_formula Formula for propensity score model: \code{treatment ~ covariates}.
#' @param censoring_formula Formula for censoring model: \code{Surv(time, event) ~ covariates}.
#'   Event should be coded as 1=event, 0=censored. Use \code{I(1-event)} if reversed.
#' @param eval_times Numeric vector of time points for evaluation. If NULL (default),
#'   uses all unique event times.
#' @param estimand Target estimand: "ATE" (average treatment effect), "ATT" (average
#'   treatment effect on the treated), or "overlap" (overlap weighting). Default "ATE".
#' @param att_group Target group for ATT. Required if \code{estimand = "ATT"}.
#' @param trim Trimming method: "symmetric" or "asymmetric". Default NULL (no trimming).
#' @param delta Threshold for symmetric trimming (e.g., 0.1). Required if \code{trim = "symmetric"}.
#' @param alpha Percentile for asymmetric trimming (e.g., 0.05). Required if \code{trim = "asymmetric"}.
#' @param contrast_matrix Optional matrix for treatment differences in multiple group settings.
#'   Each row defines one contrast with exactly two non-zero elements: -1 and 1.
#'   Column names must match treatment levels. For binary treatment, always estimates
#'   second level minus first level (S1 - S0), ignoring this parameter.
#' @param censoring_method Method for censoring score estimation: "weibull" or "cox".
#'   Default "weibull".
#' @param variance_method Variance estimation method: "analytical" (binary treatment with
#'   Weibull censoring only) or "bootstrap". Default "analytical" for binary Weibull,
#'   "bootstrap" otherwise. Cox censoring always uses bootstrap.
#' @param B Number of bootstrap iterations. Default 100. Used only if \code{variance_method = "bootstrap"}.
#' @param parallel Logical. Use parallel bootstrap computation? Default FALSE.
#' @param mc.cores Number of cores for parallel bootstrap. Default 2.
#' @param seed Random seed for bootstrap reproducibility. Default NULL.
#' @param censoring_control Control parameters passed to censoring model fitting function.
#'   For Weibull: passed to \code{survreg()}, default \code{list(maxiter = 350)}.
#'   For Cox: passed to \code{coxph()}, default \code{list()}.
#' @param ties Tie handling method for Cox models. Default "efron". Ignored for Weibull.
#' @param ps_control Control parameters for propensity score model. Default \code{list()}.
#' @param boot_level Bootstrap sampling level: "full" (default) or "strata".
#'   "full" resamples from entire dataset (standard for observational studies). "strata"
#'   resamples within each treatment group preserving group sizes (useful when treatment assignment 
#'   follows a stratified or fixed-ratio design). Only used if \code{variance_method = "bootstrap"}.
#'
#' @return List containing:
#'   \item{survival_estimates}{Matrix [time x J] of survival function estimates for each group.}
#'   \item{survival_se}{Matrix [time x J] of standard errors for survival functions.}
#'   \item{difference_estimates}{Matrix [time x K] of treatment effect estimates.
#'     For binary treatment: single column with S1-S0. For multiple groups: contrasts
#'     from \code{contrast_matrix}, or NULL if not provided.}
#'   \item{difference_se}{Matrix [time x K] of standard errors for treatment effects.}
#'   \item{eval_times}{Time points evaluated.}
#'   \item{treatment_levels}{Sorted unique treatment values.}
#'   \item{n_levels}{Number of treatment groups.}
#'   \item{n}{Sample size (complete cases after data validation).}
#'   \item{included}{Logical vector [n] indicating inclusion in analysis. TRUE = included,
#'     FALSE = excluded due to trimming.}
#'   \item{estimand}{Estimand used.}
#'   \item{censoring_method}{Censoring method used.}
#'   \item{variance_method}{Variance method used.}
#'   \item{contrast_matrix}{Contrast matrix used (NULL if not applicable).}
#'   \item{ps_model}{Fitted propensity score model (glm or multinom object).}
#'   \item{censoring_models}{Named list of fitted censoring models by treatment group.}
#'   \item{weights}{Numeric vector [n] of final weights (0 for trimmed observations).}
#'   \item{trim_summary}{Data frame with trimming summary by treatment group.}
#'   \item{ess}{Named numeric vector of effective sample size by treatment group.}
#'   \item{boot_result}{Bootstrap results (NULL if analytical variance used).}
#'
#' @details
#' **Variance Estimation:**
#' - Analytical: Binary treatment with Weibull censoring only (M-estimation).
#' - Bootstrap: All settings (resamples entire pipeline).
#' - Cox: Always uses bootstrap.
#'
#' **Treatment Effects:**
#' - Binary: S1 - S0 (second level minus first).
#' - Multiple groups: Requires \code{contrast_matrix} for pairwise comparisons.
#'
#' @examples
#' \dontrun{
#' # Binary treatment with overlap weights
#' result <- surveff(
#'   data = mydata,
#'   ps_formula = trt ~ age + sex,
#'   censoring_formula = Surv(time, event) ~ age,
#'   estimand = "overlap",
#'   censoring_method = "weibull"
#' )
#'
#' # Multiple groups with ATE and symmetric trimming
#' result <- surveff(
#'   data = mydata,
#'   ps_formula = group ~ age + sex + comorbidity,
#'   censoring_formula = Surv(time, event) ~ age + comorbidity,
#'   estimand = "ATE",
#'   trim = "symmetric",
#'   delta = 0.1,
#'   censoring_method = "cox",
#'   variance_method = "bootstrap",
#'   B = 500,
#'   parallel = TRUE,
#'   mc.cores = 4
#' )
#' }
#'
#' @export
surveff <- function(data,
                    ps_formula,
                    censoring_formula,
                    eval_times = NULL,
                    estimand = "ATE",
                    att_group = NULL,
                    trim = NULL,
                    delta = NULL,
                    alpha = NULL,
                    contrast_matrix = NULL,
                    censoring_method = "weibull",
                    variance_method = NULL,
                    B = 100,
                    parallel = FALSE,
                    mc.cores = 2,
                    seed = NULL,
                    censoring_control = NULL,
                    ties = "efron",
                    ps_control = list(),
                    boot_level = "full") {

  # Step 1: Map surveff parameters to validation function parameters
  # Check for invalid combinations
  if (estimand == "overlap" && !is.null(trim)) {
    stop("Trimming is not supported with overlap weights.\n",
         "  Overlap weights have downweighted tail populations and been bounded [0,1].\n",
         "  Use estimand = 'ATE' or 'ATT' with trim if needed.",
         call. = FALSE)
  }

  # Map estimand + trim to weight_method
  if (estimand == "overlap") {
    weight_method <- "OW"
  } else if (is.null(trim)) {
    weight_method <- "IPTW"
  } else if (trim == "symmetric") {
    weight_method <- "Symmetric"
  } else if (trim == "asymmetric") {
    weight_method <- "Asymmetric"
  } else {
    stop("Invalid trim value: '", trim, "'. Must be 'symmetric', 'asymmetric', or NULL.",
         call. = FALSE)
  }

  # Set default censoring_control based on method
  if (is.null(censoring_control)) {
    censoring_control <- if (censoring_method == "weibull") list(maxiter = 350) else list()
  }

  # Construct bootstrap_control from individual parameters
  bootstrap_control <- list(
    n_boot = B,
    parallel = parallel,
    n_cores = mc.cores,
    seed = seed
  )

  # Step 2: Validate all inputs and get complete-case dataset
  validated <- validate_PSsurvdiff_inputs(
    data = data,
    ps_formula = ps_formula,
    censor_formula = censoring_formula,
    weight_method = weight_method,
    censor_method = censoring_method,
    trim_alpha = delta,
    trim_q = alpha,
    time_points = eval_times,
    conf_level = 0.95,  # Default for validation; not used in estimation
    ps_control = ps_control,
    censor_control = censoring_control,
    bootstrap_control = bootstrap_control
  )

  # Extract validated components
  data <- validated$data_clean
  treatment_var <- validated$treatment_var
  time_var <- validated$time_var
  event_var <- validated$event_var
  # Use validated censoring formula (NULL when ~0 was specified for no censoring adjustment)
  censoring_formula_validated <- validated$censor_formula

  # Route to appropriate implementation
  if (censoring_method == "weibull") {
    result <- surveff_weibull(
      data = data,
      treatment_var = treatment_var,
      ps_formula = ps_formula,
      time_var = time_var,
      event_var = event_var,
      censoring_formula = censoring_formula_validated,
      eval_times = eval_times,
      estimand = estimand,
      att_group = att_group,
      trim = trim,
      delta = delta,
      alpha = alpha,
      contrast_matrix = contrast_matrix,
      variance_method = variance_method,
      B = B,
      parallel = parallel,
      mc.cores = mc.cores,
      seed = seed,
      censoring_control = censoring_control,
      ps_control = ps_control,
      boot_level = boot_level
    )
  } else if (censoring_method == "cox") {
    # Handle variance_method validation for Cox
    if (!is.null(variance_method) && variance_method == "analytical") {
      message("Analytical variance is not available for Cox censoring. Using bootstrap instead.")
      variance_method <- "bootstrap"
    }

    result <- surveff_cox(
      data = data,
      treatment_var = treatment_var,
      ps_formula = ps_formula,
      time_var = time_var,
      event_var = event_var,
      censoring_formula = censoring_formula_validated,
      eval_times = eval_times,
      estimand = estimand,
      att_group = att_group,
      trim = trim,
      delta = delta,
      alpha = alpha,
      contrast_matrix = contrast_matrix,
      B = B,
      parallel = parallel,
      mc.cores = mc.cores,
      seed = seed,
      censoring_control = censoring_control,
      ties = ties,
      ps_control = ps_control,
      boot_level = boot_level
    )
  } else {
    stop("Invalid censoring_method: '", censoring_method, "'. ",
         "Must be 'weibull' or 'cox'.", call. = FALSE)
  }

  # Format output - extract and flatten nested structures
  surv_res <- result$surv_result
  weight_res <- surv_res$weight_result
  ps_res <- surv_res$ps_result
  cens_res <- surv_res$censoring_result

  # Create included vector (TRUE = included after trimming, FALSE = trimmed)
  n_original <- nrow(data)
  included <- weight_res$weights > 0

  output <- list(
    # Primary results
    survival_estimates = result$survival_estimates,
    survival_se = sqrt(result$survival_variances),
    difference_estimates = result$difference_estimates,
    difference_se = if (!is.null(result$difference_variances)) sqrt(result$difference_variances) else NULL,

    # Study design info
    eval_times = result$eval_times,
    treatment_levels = result$treatment_levels,
    n_levels = result$n_levels,
    n = n_original,
    included = included,

    # Method specifications
    estimand = surv_res$estimand,
    censoring_method = censoring_method,
    variance_method = result$variance_method,
    contrast_matrix = result$contrast_matrix,

    # Model objects (flattened from nested structures)
    ps_model = ps_res$ps_model,
    censoring_models = cens_res$censoring_models,

    # Weights summary (flattened from nested structures)
    weights = weight_res$weights,
    trim_summary = weight_res$trim_summary,
    ess = weight_res$ess,
    ps_matrix = ps_res$ps_matrix,
    ps_obs = ps_res$ps,
    
    # censoring summary
    censoring_scores = cens_res$censoring_scores,
    
    # Bootstrap results
    boot_result = result$boot_result
  )

  class(output) <- "surveff"
  return(output)
}


#' Print Method for surveff Objects
#'
#' @param x A \code{surveff} object.
#' @param max.len Maximum number of rows (time points) to print. Default 6.
#' @param round.digits Number of digits for rounding displayed values. Default 4.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.surveff <- function(x, max.len = 6, round.digits = 4, ...) {
  cat("\n=== Survival Effect Estimation ===\n\n")
  cat("Treatment groups:", paste(x$treatment_levels, collapse = ", "), "\n")
  cat("Number of groups:", x$n_levels, "\n")
  cat("Estimand:", x$estimand, "\n")
  cat("Censoring method:", x$censoring_method, "\n")
  cat("Variance method:", x$variance_method, "\n")
  cat("Evaluation times:", length(x$eval_times), "time points\n\n")

  # Determine rows to display
  n_times <- nrow(x$survival_estimates)
  rows_to_show <- if (n_times <= max.len) 1:n_times else 1:max.len
  show_note <- n_times > max.len

  cat("Survival function estimates:\n")
  print(round(x$survival_estimates[rows_to_show, , drop = FALSE], round.digits))
  if (show_note) cat("  ... (", n_times - max.len, " more rows not shown)\n", sep = "")

  if (!is.null(x$difference_estimates)) {
    cat("\nTreatment effect estimates:\n")
    print(round(x$difference_estimates[rows_to_show, , drop = FALSE], round.digits))
    if (show_note) cat("  ... (", n_times - max.len, " more rows not shown)\n", sep = "")
  } else {
    cat("\nTreatment effects: Not estimated (provide contrast_matrix for multiple groups)\n")
  }

  invisible(x)
}


#' Summary Method for surveff Objects
#'
#' @param object A \code{surveff} object.
#' @param conf_level Confidence level for intervals. Default 0.95.
#' @param max.len Maximum number of rows (time points) to print. Default 6.
#'   Only used if \code{style = "prints"}.
#' @param round.digits Number of digits for rounding displayed values. Default 4.
#'   Only used if \code{style = "prints"}.
#' @param style Output style: "prints" (print formatted tables) or "returns"
#'   (return list of matrices). Default "prints".
#' @param ... Additional arguments (ignored).
#'
#' @return If \code{style = "prints"}, returns invisibly. If \code{style = "returns"},
#'   returns a list with:
#'   \item{survival_summary}{List of matrices, one per treatment group, with columns:
#'     Time, Estimate, SE, CI.lower, CI.upper}
#'   \item{difference_summary}{List of matrices, one per contrast, with same columns.
#'     NULL if no contrasts estimated.}
#'
#' @export
summary.surveff <- function(object, conf_level = 0.95, max.len = 6, round.digits = 4,
                            style = "prints", ...) {
  if (!style %in% c("prints", "returns")) {
    stop("Invalid 'style' argument. Must be 'prints' or 'returns'.", call. = FALSE)
  }

  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  ci_pct <- round(conf_level * 100)

  # Helper to build result table with CIs
  build_table <- function(est, se, truncate_01 = FALSE, round_vals = FALSE) {
    lower <- if (truncate_01) pmax(0, est - z_crit * se) else est - z_crit * se
    upper <- if (truncate_01) pmin(1, est + z_crit * se) else est + z_crit * se

    if (round_vals) {
      tbl <- data.frame(Time = object$eval_times, Estimate = round(est, round.digits),
                        SE = round(se, round.digits), CI.lower = round(lower, round.digits),
                        CI.upper = round(upper, round.digits), row.names = NULL)
    } else {
      tbl <- data.frame(Time = object$eval_times, Estimate = est, SE = se,
                        CI.lower = lower, CI.upper = upper, row.names = NULL)
    }
    colnames(tbl)[4:5] <- paste0(ci_pct, "%CI.", c("lower", "upper"))
    tbl
  }

  if (style == "returns") {
    # Build survival summary
    survival_list <- lapply(1:object$n_levels, function(j) {
      grp <- as.character(object$treatment_levels[j])
      build_table(object$survival_estimates[, grp], object$survival_se[, grp], truncate_01 = TRUE)
    })
    names(survival_list) <- as.character(object$treatment_levels)

    # Build difference summary
    difference_list <- NULL
    if (!is.null(object$difference_estimates)) {
      difference_list <- lapply(1:ncol(object$difference_estimates), function(k) {
        build_table(object$difference_estimates[, k], object$difference_se[, k])
      })
      names(difference_list) <- colnames(object$difference_estimates)
    }

    return(list(survival_summary = survival_list, difference_summary = difference_list))
  }

  # style == "prints"
  n_times <- nrow(object$survival_estimates)
  rows_to_show <- if (n_times <= max.len) 1:n_times else 1:max.len
  show_note <- n_times > max.len

  cat("\n=== Survival Effect Estimation Summary ===\n\n")
  cat("Treatment groups:", paste(object$treatment_levels, collapse = ", "), "\n")
  cat("Number of groups:", object$n_levels, "\n")
  cat("Sample size:", object$n, "complete cases\n")
  cat("Estimand:", object$estimand, "\n")
  cat("Censoring method:", object$censoring_method, "\n")
  cat("Variance method:", object$variance_method, "\n")
  if (object$variance_method == "bootstrap") cat("Bootstrap iterations:", object$boot_result$B, "\n")
  cat("Evaluation times:", length(object$eval_times), "time points\n")
  cat("Confidence level:", conf_level, "\n\n")

  # Print survival estimates
  cat("--- Survival Function Estimates ---\n\n")
  for (j in 1:object$n_levels) {
    grp <- as.character(object$treatment_levels[j])
    cat("Group", grp, ":\n")
    tbl <- build_table(object$survival_estimates[, grp], object$survival_se[, grp],
                       truncate_01 = TRUE, round_vals = TRUE)
    print(tbl[rows_to_show, , drop = FALSE], row.names = FALSE)
    if (show_note) cat("  ... (", n_times - max.len, " more rows not shown)\n", sep = "")
    cat("\n")
  }

  # Print difference estimates
  if (!is.null(object$difference_estimates)) {
    cat("--- Treatment Effect Estimates ---\n\n")
    for (k in 1:ncol(object$difference_estimates)) {
      cat(colnames(object$difference_estimates)[k], ":\n")
      tbl <- build_table(object$difference_estimates[, k], object$difference_se[, k],
                         round_vals = TRUE)
      print(tbl[rows_to_show, , drop = FALSE], row.names = FALSE)
      if (show_note) cat("  ... (", n_times - max.len, " more rows not shown)\n", sep = "")
      cat("\n")
    }
  } else {
    cat("--- Treatment Effects ---\n")
    cat("Not estimated (provide contrast_matrix for multiple groups)\n\n")
  }

  invisible(object)
}


#' Plot Method for surveff Objects
#'
#' @param x A \code{surveff} object.
#' @param type Type of plot: "surv" for survival curves or "survdiff" for
#'   treatment effect curves. Default "surv".
#' @param max_time Maximum time to display on x-axis. If NULL, uses max(eval_times).
#' @param strata_to_plot Vector of strata to plot. For \code{type = "surv"}, must be
#'   subset of treatment_levels. For \code{type = "survdiff"}, must be subset of
#'   contrast names (column names of difference_estimates). If NULL, plots all available strata.
#' @param strata_colors Vector of color names/codes for strata. Length must match
#'   strata_to_plot. Order matches strata order. If NULL, uses ggplot2 default colors.
#' @param conf_level Confidence level for confidence intervals. Default 0.95.
#' @param include_CI Logical. Include confidence interval ribbons? Default TRUE.
#' @param curve_width Line width for survival/difference curves. Default 1.
#' @param CI_alpha Transparency level for CI ribbons (0-1). Default 0.3.
#' @param legend_position Position of legend: "right" or "bottom". Default "right".
#' @param legend_title Title for legend. If NULL, uses "Treatment" for type="surv"
#'   or "Comparison" for type="survdiff".
#' @param plot_title Plot title. If NULL, uses default title based on type.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot2 object.
#'
#' @details
#' Creates publication-ready plots of survival curves or treatment effects over time.
#'
#' For \code{type = "surv"}: Plots estimated survival functions with optional
#' confidence intervals. Y-axis ranges from 0 to 1.
#'
#' For \code{type = "survdiff"}: Plots estimated treatment effects (survival differences)
#' with optional confidence intervals. Y-axis is not constrained to [0,1].
#'
#' @export
plot.surveff <- function(x, type = "surv", max_time = NULL, strata_to_plot = NULL,
                         strata_colors = NULL, conf_level = 0.95, include_CI = TRUE,
                         curve_width = 1, CI_alpha = 0.3, legend_position = "right",
                         legend_title = NULL, plot_title = NULL, ...) {

  # Check ggplot2 availability
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.", call. = FALSE)
  }

  # Validate inputs
  if (!type %in% c("surv", "survdiff")) {
    stop("'type' must be 'surv' or 'survdiff'.", call. = FALSE)
  }
  if (!legend_position %in% c("right", "bottom")) {
    stop("'legend_position' must be 'right' or 'bottom'.", call. = FALSE)
  }
  if (!is.logical(include_CI)) {
    stop("'include_CI' must be logical.", call. = FALSE)
  }
  if (CI_alpha < 0 || CI_alpha > 1) {
    stop("'CI_alpha' must be between 0 and 1.", call. = FALSE)
  }
  if (type == "survdiff" && is.null(x$difference_estimates)) {
    stop("No treatment effect estimates available. Provide contrast_matrix to surveff() for multiple groups.",
         call. = FALSE)
  }

  # Determine max_time
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

  # Get summary with CIs (used for both types)
  summ <- summary(x, conf_level = conf_level, style = "returns")

  # Determine strata and prepare data based on type
  if (type == "surv") {
    # Determine strata (treatment groups)
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

    # Build data frame for survival curves
    plot_data <- do.call(rbind, lapply(strata_to_plot, function(grp) {
      summ_grp <- summ$survival_summary[[grp]]
      data.frame(
        Time = summ_grp$Time[time_idx],
        Estimate = summ_grp$Estimate[time_idx],
        CI_lower = summ_grp[time_idx, 4],
        CI_upper = summ_grp[time_idx, 5],
        Strata = grp,
        stringsAsFactors = FALSE
      )
    }))

    ylab <- "Survival Probability"
    ylim <- c(0, 1)
    default_legend_title <- "Treatment"
    default_plot_title <- "Estimated Survival Curves by Groups"

  } else {  # survdiff
    # Determine strata (contrasts)
    all_strata <- colnames(x$difference_estimates)
    if (is.null(strata_to_plot)) {
      strata_to_plot <- all_strata
    } else {
      strata_to_plot <- as.character(strata_to_plot)
      if (!all(strata_to_plot %in% all_strata)) {
        stop("'strata_to_plot' must be a subset of contrast names: ",
             paste(all_strata, collapse = ", "), call. = FALSE)
      }
    }

    # Build data frame for treatment effect curves
    plot_data <- do.call(rbind, lapply(strata_to_plot, function(contrast) {
      summ_contrast <- summ$difference_summary[[contrast]]
      data.frame(
        Time = summ_contrast$Time[time_idx],
        Estimate = summ_contrast$Estimate[time_idx],
        CI_lower = summ_contrast[time_idx, 4],
        CI_upper = summ_contrast[time_idx, 5],
        Strata = contrast,
        stringsAsFactors = FALSE
      )
    }))

    ylab <- "Survival Difference"
    ylim <- NULL
    default_legend_title <- "Comparison"
    default_plot_title <- "Estimated Treatment Effect Curves"
  }

  # Validate strata_colors
  if (!is.null(strata_colors)) {
    if (length(strata_colors) != length(strata_to_plot)) {
      stop("'strata_colors' must have length ", length(strata_to_plot), " (matching strata_to_plot).",
           call. = FALSE)
    }
  }

  # Set factor levels for consistent ordering
  plot_data$Strata <- factor(plot_data$Strata, levels = strata_to_plot)

  # Set default titles if not provided
  if (is.null(plot_title)) plot_title <- default_plot_title
  if (is.null(legend_title)) legend_title <- default_legend_title

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Time, y = Estimate, color = Strata, fill = Strata))

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
  p <- p + ggplot2::scale_x_continuous(limits = c(0, max_time), expand = c(0, 0))
  if (!is.null(ylim)) {
    p <- p + ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0))
  }

  p <- p + ggplot2::labs(
    title = plot_title,
    x = "Time",
    y = ylab,
    color = legend_title,
    fill = legend_title
  )

  # Apply custom colors if provided
  if (!is.null(strata_colors)) {
    p <- p + ggplot2::scale_color_manual(values = strata_colors) +
      ggplot2::scale_fill_manual(values = strata_colors)
  }

  # Theme with light grey grid, centered title, no top/right borders
  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = legend_position,
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_line(color = "grey95", linewidth = 0.2),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.5)
    )

  return(p)
}
