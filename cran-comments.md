# CRAN Submission Comments - PSsurvival 0.1.0

## Submission Date
2025-01-24

## Test Environments
- Local: macOS Sequoia 15.5, R 4.5.1 (0 errors, 0 warnings, 0 notes)
- WinBuilder R-release 4.5.2 (0 errors, 0 warnings, 1 note)
- WinBuilder R-devel (0 errors, 0 warnings, 1 note)
- GitHub Actions:
  - macOS-latest, R-release
  - Windows-latest, R-release
  - Ubuntu-latest, R-devel
  - Ubuntu-latest, R-release
  - Ubuntu-latest, R-oldrel-1

## R CMD check results
0 errors | 0 warnings | 1 note

The single NOTE is:
```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Chengxin Yang <chengxin.yang@duke.edu>'

New submission

Possibly misspelled words in DESCRIPTION:
  ATT (23:12)
  Cheng (26:22)
  Crump (23:71)
  Sturmer (24:21)
  al (26:31)
  estimands (23:40)
  et (26:28)
```

These are all correct:
- ATT: Acronym for "Average Treatment effect on the Treated"
- Cheng, Crump, Sturmer: Author names from references
- estimands: Valid statistical term (plural of estimand)
- et al: Standard Latin abbreviation in citations

## Submission Comments Sent to CRAN
This is a new submission.

The package implements propensity score weighting methods for causal
survival analysis with time-to-event outcomes.

Test results:
- Local macOS: 0 errors, 0 warnings, 0 notes
- WinBuilder R-release: 0 errors, 0 warnings, 1 note
- WinBuilder R-devel: 0 errors, 0 warnings, 1 note
- GitHub Actions (5 platforms): All passing

The single NOTE is for "New submission" and flagged technical terms
(ATT, estimands) and author names (Cheng, Crump, Sturmer) in the
DESCRIPTION, which are all correct.

## Notes
- First CRAN submission
- 99 passing unit tests using testthat
- Comprehensive vignette included
- All examples run successfully
