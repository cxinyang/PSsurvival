# PSsurvival

<!-- badges: start -->
[![R-CMD-check](https://github.com/cxinyang/PSsurvival/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cxinyang/PSsurvival/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Propensity score methods for survival analysis.

## Overview

PSsurvival implements propensity score methods for observational studies with time-to-event outcomes. The package provides three main functions:

**Counterfactual survival functions** (`surveff`): Estimates group-specific survival curves and survival differences over time, adjusting for confounding via propensity score weighting and for censoring via inverse probability of censoring weighting.

**Marginal hazard ratios** (`marCoxph`): Fits weighted marginal Cox proportional hazards models to estimate marginal hazard ratios between treatment groups.

**Weighted Kaplan-Meier curves** (`weightedKM`): Estimates weighted Kaplan-Meier (KM) and cumulative risk (CR) curves with propensity score weights.

All functions support:

- Binary and multiple (>2) treatment groups
- Weighting schemes: inverse probability weighting (IPW), overlap weighting (OW), and average treatment effect on the treated (ATT)
- Propensity score trimming: symmetric trimming (Crump extension for multiple treatments)
- Variance estimation: analytical M-estimation (binary treatment with Weibull censoring) or bootstrap

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("cxinyang/PSsurvival")

# Install from R CRAN
install.package("PSsurvival")
```

## Usage

```r
library(PSsurvival)

# Counterfactual survival curves with overlap weighting
result <- surveff(
  data = mydata,
  ps_formula = treatment ~ X1 + X2 + X3,
  censoring_formula = Surv(time, event) ~ X1 + X2,
  weight_method = "OW",
  censoring_method = "weibull"
)

summary(result)
plot(result)

# Marginal hazard ratio with IPW
hr_result <- marCoxph(
  data = mydata,
  ps_formula = treatment ~ X1 + X2 + X3,
  time_var = "time",
  event_var = "event",
  reference_level = "control",
  weight_method = "IPW"
)

summary(hr_result)

# Weighted Kaplan-Meier with risk table
km_result <- weightedKM(
  data = mydata,
  treatment_var = "treatment",
  ps_formula = treatment ~ X1 + X2 + X3,
  time_var = "time",
  event_var = "event",
  weight_method = "OW"
)

plot(km_result, risk_table = TRUE)
summary(km_result)
```

## Details

**Propensity score estimation**: Uses logistic regression for binary treatments and multinomial logistic regression (via `nnet::multinom`) for multiple treatments.

**Censoring adjustment** (`surveff` only): Models the censoring distribution within each treatment group using either Weibull accelerated failure time models or Cox proportional hazards models.

**Variance estimation**: For binary treatments with Weibull censoring, analytical variance based on M-estimation theory is available. Bootstrap variance (resampling the full estimation pipeline) is supported for all configurations.

## References
Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via propensity score weighting. *Journal of the American Statistical Association*, 113(521), 390-400.

Li, F., & Li, F. (2019). Propensity score weighting for causal inference with multiple treatments. *The Annals of Applied Statistics*, 13(4), 2389-2415.

Cheng, C., Li, F., Thomas, L. E., & Li, F. (2022). Addressing extreme propensity scores in estimating counterfactual survival functions via the overlap weights. *American Journal of Epidemiology*, 191(6), 1140-1151.
