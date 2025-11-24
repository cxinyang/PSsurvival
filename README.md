# PSsurvival

<!-- badges: start -->
[![R-CMD-check](https://github.com/cxinyang/PSsurvival/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cxinyang/PSsurvival/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Propensity score weighting methods for causal survival analysis with binary and multiple treatments.

## Overview

PSsurvival implements inverse probability weighting approaches for estimating causal effects in observational studies with time-to-event outcomes. The package provides two main analyses:

**Counterfactual survival functions** (`surveff`): Estimates group-specific survival curves and survival differences over time, adjusting for confounding via propensity score weighting and for censoring via inverse probability of censoring weighting.

**Marginal hazard ratios** (`marCoxph`): Fits weighted Cox proportional hazards models to estimate marginal hazard ratios between treatment groups.

Both functions support:

- Binary and multiple (>2) treatment groups
- Weighting schemes: inverse probability weighting (ATE), treated/target group weighting (ATT), and overlap weighting
- Propensity score trimming: symmetric (Crump) and asymmetric (Sturmer) methods
- Variance estimation: analytical M-estimation (binary treatment with Weibull censoring) or bootstrap

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("cxinyang/PSsurvival")
```

## Usage

```r
library(PSsurvival)

# Counterfactual survival curves with overlap weighting
result <- surveff(
  data = mydata,
  ps_formula = treatment ~ X1 + X2 + X3,
  censoring_formula = Surv(time, event) ~ X1 + X2,
  estimand = "overlap",
  censoring_method = "weibull"
)

summary(result)
plot(result)

# Marginal hazard ratio with ATE weighting
hr_result <- marCoxph(
  data = mydata,
  ps_formula = treatment ~ X1 + X2 + X3,
  time_var = "time",
  event_var = "event",
  reference_level = "control",
  estimand = "ATE"
)

summary(hr_result)
```

## Details

**Propensity score estimation**: Uses logistic regression for binary treatments and multinomial logistic regression (via `nnet::multinom`) for multiple treatments.

**Censoring adjustment** (`surveff` only): Models the censoring distribution within each treatment group using either Weibull accelerated failure time models or Cox proportional hazards models.

**Variance estimation**: For binary treatments with Weibull censoring, analytical variance based on M-estimation theory is available. Bootstrap variance (resampling the full estimation pipeline) is supported for all configurations.

## References
Cheng, C., Li, F., Thomas, L. E., & Li, F. (2022). Addressing extreme propensity scores in estimating counterfactual survival functions via the overlap weights. *American Journal of Epidemiology*, 191(6), 1140-1151.

Li, F., & Li, F. (2019). Propensity score weighting for causal inference with multiple treatments. *The Annals of Applied Statistics*, 13(4), 2389-2415.
