# CRAN Resubmission Comments - PSsurvival 0.1.0

## Resubmission Date
2025-12-02

## Response to Reviewer Feedback

Thank you for reviewing our package. We have addressed all issues raised:

### 1. Missing \value Documentation for Print Methods

**Issue:** `print.surveff()` and `print.marCoxph()` lacked `\value{}` tags in documentation.

**Fixed:** Added `@return` tags to both functions (R/A00-surveff.R:312, R/B00-marCoxph.R:377). Both now document: "Invisibly returns the input object \code{x}."

### 2. Examples for Unexported Functions

**Issue:** Two internal functions (`estimate_censoring_score_weibull()`, `estimate_censoring_score_cox()`) included `\dontrun{}` examples.

**Fixed:** Removed `@examples` sections from both functions (R/U04-censoring-scores.R). All 25 internal functions now have no examples, as required for `@keywords internal` functions.

### 3. \dontrun{} vs \donttest{} Usage

**Issue:** Examples used `\dontrun{}` instead of executable `\donttest{}`.

**Fixed:** All examples for exported functions (`surveff()`, `marCoxph()`, `estimate_ps()`, `estimate_weights()`) now use `\donttest{}` and run with package data (`simdata_bin`, `simdata_multi`). Examples are executable and demonstrate realistic usage.

### 4. Output Documentation Completeness

**Issue:** Reviewer requested clear documentation of output structure, class, and meaning for all exported functions.

**Verified:** All 9 exported functions have comprehensive `\value{}` sections explaining:
- Structure: List components, dimensions, data types
- Class: Object classes (e.g., "marCoxph", "surveff", ggplot2)
- Meaning: Interpretation of each output component

Spot-checked examples:
- `surveff()`: 15 components documented with types and dimensions
- `marCoxph()`: 18 components with class name explicitly stated
- `estimate_ps()`: 5 components with matrix dimensions specified
- `estimate_weights()`: 12 components with detailed explanations

## Test Environments

Tested on:
- Local: macOS Sequoia 15.5, R 4.5.1
- GitHub Actions (5 configurations)

## R CMD check Results

**Current status:** 0 errors | 0 warnings | 0 notes

All documentation issues resolved. The package passes `devtools::check()` with no issues.

## Additional Changes

- Removed `PSweight_estimate` class and its print/summary methods for API consistency with `estimate_ps()` (which is also exported without S3 methods)
- All internal functions verified to have `@keywords internal` and no examples
- All exported functions with examples use `\donttest{}` with executable code
- Changed the title to "Propensity Score Methods for Survival Analysis" for succinctness

## Notes

- This is a resubmission addressing reviewer feedback
- 99 unit tests passing (testthat)
- Comprehensive documentation verified for all functions
- Examples use package datasets and demonstrate intended workflows
