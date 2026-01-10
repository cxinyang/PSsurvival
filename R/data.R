#' Simulated Survival Data with Binary Treatment
#'
#' A simulated dataset for demonstrating propensity score weighting methods
#' in survival analysis with a binary treatment.
#'
#' @format A data frame with 1000 observations and 8 variables:
#' \describe{
#'   \item{X1}{Continuous covariate (standard normal).}
#'   \item{X2}{Continuous covariate (standard normal).}
#'   \item{X3}{Continuous covariate (standard normal).}
#'   \item{B1}{Binary covariate (0/1).}
#'   \item{B2}{Binary covariate (0/1).}
#'   \item{Z}{Treatment group: "A" or "B". Distribution is approximately 40:60.}
#'   \item{time}{Observed follow-up time (event or censoring), range 0-20.}
#'   \item{event}{Event indicator: 1 = event observed, 0 = censored.}
#' }
#'
#' @details
#' The data were generated with the following characteristics:
#' \itemize{
#'   \item Treatment assignment depends on X1, X2, and B1 via logistic model.
#'   \item Survival times follow Weibull distributions with group-specific scales
#'     (group A has better survival than group B).
#'   \item Censoring times follow an exponential distribution depending on X1 and B1.
#'   \item Administrative censoring occurs at time 20.
#'   \item Overall censoring rate is approximately 30%.
#' }
#'
#' @examples
#' data(simdata_bin)
#' head(simdata_bin)
#' table(simdata_bin$Z)
#'
#' @seealso \code{\link{simdata_multi}} for a dataset with 4 treatment groups.
"simdata_bin"


#' Simulated Survival Data with Multiple Treatments
#'
#' A simulated dataset for demonstrating propensity score weighting methods
#' in survival analysis with four treatment groups.
#'
#' @format A data frame with 1000 observations and 8 variables:
#' \describe{
#'   \item{X1}{Continuous covariate (standard normal).}
#'   \item{X2}{Continuous covariate (standard normal).}
#'   \item{X3}{Continuous covariate (standard normal).}
#'   \item{B1}{Binary covariate (0/1).}
#'   \item{B2}{Binary covariate (0/1).}
#'   \item{Z}{Treatment group: "A", "B", "C", or "D". Distribution is approximately
#'     20:20:20:35.}
#'   \item{time}{Observed follow-up time (event or censoring), range 0-20.}
#'   \item{event}{Event indicator: 1 = event observed, 0 = censored.}
#' }
#'
#' @details
#' The data were generated with the following characteristics:
#' \itemize{
#'   \item Treatment assignment depends on X1, X2, X3, B1, and B2 via multinomial
#'     logistic model.
#'   \item Survival times follow Weibull distributions with group-specific scales.
#'     Survival ordering (best to worst): C > A > B > D.
#'   \item Censoring times follow an exponential distribution depending on X1 and B1.
#'   \item Administrative censoring occurs at time 20.
#'   \item Overall censoring rate is approximately 30%.
#' }
#'
#' @examples
#' data(simdata_multi)
#' head(simdata_multi)
#' table(simdata_multi$Z)
#'
#' @seealso \code{\link{simdata_bin}} for a dataset with binary treatment.
"simdata_multi"
