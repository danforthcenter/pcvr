# pcvr 1.0.0.5

Added S3 class for `growthSS` and `mvSS` output (`pcvrss`) with print/summary methods.

# pcvr 1.0.0.4

Added error handling for examples that read data from github in case they run in a session without an
internet connection.

# pcvr 1.0.0.3

Changes to `frem` to allow for using datetimes and for cases where a single timepoint cannot be modeled
due to singular grouping.

# pcvr 1.0.0.2

Change plotting methods to avoid printing internal labels.

Fixed problems where `growthSim` could simulate changepoint data incorrectly for different changepoints between groups and where the `Stan` `inv_logit` function may not be available for model prediction.

# pcvr 1.0.0.1

New `mvSS` function for simplified interface to modeling multi-value traits via `growthSS`, either longitudinally or at one timepoint.

# pcvr 1.0.0

* Initial CRAN submission.
