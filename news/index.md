# Changelog

## pcvr 1.4.0

Added `stat_growthss` to make ggplot layers of models from `growthSS`.

## pcvr 1.3.1

CRAN release: 2025-08-18

Cran release

## pcvr 1.3.0.4

Simplified logic in `frem` to always use an interaction term between
design variables.

## pcvr 1.3.0.3

Adding option to return exact parameters per each subject simulated in
`growthSim`.

## pcvr 1.3.0.2

Bug fixes for growthPlot with brms models that have only 1 group.

Explicit error handling for conjugate samples of 1L.

Update to field capacity calculations in articles.

## pcvr 1.3.0.1

Added logic to joyplots to handle continuous variables on y axis (time).

## pcvr 1.3.0

Added 4 and 5 parameter logistic curves to growthSS

Changed `bw.*` functions to `pcv.*` prefixes to reduce bellwether system
artefacts.

## pcvr 1.2.0.1

Allowing `brms` models from `growthSS` to specify subset of parameters
to estimate per group.

## pcvr 1.2.0

CRAN release: 2025-04-16

Added `conjugate` S3 method for outputs from conjugate class. Changed
internals for several conjugate distributions. Added Bayes factors to
conjugate.

## pcvr 1.1.1.0

CRAN release: 2024-11-06

Fixing inelegant failure on `cran` when building the bellwether vignette
without wifi and resubmitting to `cran`.

## pcvr 1.1.0.1

Adding arguments to `mvSim` to allow for easier simulation of
longitudinal multi-value traits from the growth model options supported
in `growthSim`.

## pcvr 1.1.0.0

Allowing for multiple grouping variables to be used in `growthSS` and
downstream functions.

## pcvr 1.0.0.6

Fixing `growthSS` behavior with non-integer time options.

## pcvr 1.0.0.5

Added S3 class for `growthSS` and `mvSS` output (`pcvrss`) with
print/summary methods.

Bug fixes for brms models using splines without grouping variables.

Bug fixes for brms model plotting without grouping variables.

## pcvr 1.0.0.4

Added error handling for examples that read data from github in case
they run in a session without an internet connection.

## pcvr 1.0.0.3

Changes to `frem` to allow for using datetimes and for cases where a
single timepoint cannot be modeled due to singular grouping.

## pcvr 1.0.0.2

Change plotting methods to avoid printing internal labels.

Fixed problems where `growthSim` could simulate changepoint data
incorrectly for different changepoints between groups and where the
`Stan` `inv_logit` function may not be available for model prediction.

## pcvr 1.0.0.1

New `mvSS` function for simplified interface to modeling multi-value
traits via `growthSS`, either longitudinally or at one timepoint.

## pcvr 1.0.0

CRAN release: 2024-09-05

- Initial CRAN submission.
