# PSsurvival 0.2.0

## Breaking Changes

* **API Unification**: Changed `estimand` parameter to `weight_method` in `surveff()` and `marCoxph()` for consistency across all package functions
  - `estimand = "ATE"` → `weight_method = "IPW"`
  - `estimand = "overlap"` → `weight_method = "OW"`
  - `estimand = "ATT"` → `weight_method = "ATT"` (name unchanged, but parameter name changed)

* Changed `trim` parameter from character to logical in `surveff()`, `marCoxph()`, and `weightedKM()`
  - `trim = "symmetric"` → `trim = TRUE`
  - `trim = NULL` → `trim = FALSE`

* Removed `alpha` parameter and asymmetric trimming support from all functions due to poor statistical performance in practice
  - Only symmetric trimming (Crump extension for multiple treatments) is now supported

* Default `delta` is now `NULL` (automatic selection based on number of treatment groups) instead of a fixed value

## New Features

* **New function**: `weightedKM()` for weighted Kaplan-Meier estimation with propensity score weights
  - Supports IPW, overlap weights, ATT, and classical unweighted Kaplan-Meier
  - Uses classic weighted Greenwood variance formula
  - Handles binary and multiple treatment groups

* **New S3 methods** for `weightedKM` objects:
  - `plot.weightedKM()`: Visualize weighted Kaplan-Meier or cumulative risk curves with confidence intervals
  - `summary.weightedKM()`: Tabular summaries with confidence intervals
  - `print.weightedKM()`: Formatted console output

* **Risk table support** in `plot.weightedKM()`:
  - Display number at risk and cumulative events below survival/risk curves
  - Customizable statistics, breaks, height, and font size
  - Automatic group name positioning and formatting

* **Enhanced visualization options**:
  - Three confidence interval types: plain, log, log-log (default)
  - Cumulative risk plots (1 - survival)
  - Customizable colors, widths, transparency, and legends

## Documentation

* Unified weighting method descriptions across all main functions (`surveff()`, `marCoxph()`, `weightedKM()`)
* Added comprehensive `@details` sections explaining IPW, OW, and ATT weighting methods
* Updated all examples to use new API with `weight_method` parameter

## Internal

* Moved `ggplot2` and `cowplot` from Suggests to Imports (required for plotting functions)
* Parameter mapping layer in `surveff()` and `marCoxph()` to transform user-facing API to internal implementation

## Migration Guide for Existing Users

If upgrading from version 0.1.0:

* Replace `estimand = "ATE"` with `weight_method = "IPW"`
* Replace `estimand = "overlap"` with `weight_method = "OW"`
* Replace `estimand = "ATT"` with `weight_method = "ATT"`
* Replace `trim = "symmetric", delta = 0.1` with `trim = TRUE, delta = 0.1`
* Replace `trim = NULL` with `trim = FALSE` (or simply omit, as FALSE is default)
* Remove `alpha` parameter entirely (asymmetric trimming no longer supported)
* Set `delta = NULL` to use automatic defaults (0.1 for binary, 0.067 for 3 groups, 1/(2J) for J ≥ 4)

Example migration:
```r
# Old (v0.1.0)
result <- surveff(
  data = mydata,
  ps_formula = Z ~ X1 + X2,
  censoring_formula = Surv(time, event) ~ X1,
  estimand = "overlap",
  censoring_method = "weibull"
)

# New (v0.2.0)
result <- surveff(
  data = mydata,
  ps_formula = Z ~ X1 + X2,
  censoring_formula = Surv(time, event) ~ X1,
  weight_method = "OW",
  censoring_method = "weibull"
)
```

# PSsurvival 0.1.0

* Initial CRAN release

* Documentation improvements for CRAN resubmission:
  - Added `@return` tags to print methods
  - Removed examples from internal functions
  - Updated examples to use `\donttest{}` with executable code
  - Enhanced output documentation for all exported functions

* Implements propensity score weighting methods for estimating counterfactual
  survival functions and marginal hazard ratios in observational studies with
  time-to-event outcomes

* Main functions:
  - `surveff()`: Estimates counterfactual survival curves and survival
    differences over time
  - `marCoxph()`: Estimates marginal hazard ratios via weighted Cox regression

* Supports binary and multiple (>2) treatment groups

* Weighting methods:
  - Inverse probability weighting (ATE)
  - Treated/target group weighting (ATT)
  - Overlap weighting

* Propensity score trimming:
  - Symmetric trimming (Crump et al.)
  - Asymmetric trimming (Sturmer et al.)

* Censoring adjustment methods:
  - Weibull accelerated failure time models
  - Cox proportional hazards models

* Variance estimation:
  - Analytical M-estimation (binary treatment with Weibull censoring)
  - Bootstrap (all configurations, with full or stratified resampling)

* Comprehensive documentation including vignette with examples

* Based on methods from Cheng et al. (2022) <doi:10.1093/aje/kwac043> and
  Li & Li (2019) <doi:10.1214/19-AOAS1282>
