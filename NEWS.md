# PSsurvival 0.1.0

* Initial CRAN release (under review)

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
