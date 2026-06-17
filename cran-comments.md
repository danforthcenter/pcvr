# pcvr 1.4.1

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 1 notes ✖

❯ checking for future file timestamps ... NOTE
  unable to verify current time

## Notes

Resubmitting due to cran testing error probably related to either not installing
suggested packages or the change from `mc-stan.org` to `stan-dev.r--universe.dev`
for some additional packages and for development features
(ggplot2 geom support, more distributions, etc).
Previously a NOTE related to the failure specified:

```
❯ checking package dependencies ... NOTE
Package suggested but not available for checking: ‘cmdstanr’
Suggests or Enhances not in mainstream repositories:
  cmdstanr
Availability using Additional_repositories specification:
  cmdstanr   yes   https://stan-dev.r-universe.dev
```

But the `cmdstanr` suggestion is available from the new additional repository specified in DESCRIPTION.

No notes seem critical.

Names (PlantCV, Kruschke) and words (phenotyping) in DESCRIPTION are not misspelled.
The R CMD check on MacOs via github actions yields a NOTE about installed package size as well.

Additionally there are some examples wrapped in `\donttest` that may take several minutes to run.

## Test environments

```
- {os: macos-latest,   r: 'release'}
- {os: windows-latest, r: 'release'}
- {os: windows-latest, r: 'devel'}
- {os: ubuntu-latest,   r: 'devel'}
- {os: ubuntu-latest,   r: 'release'}
- {os: ubuntu-latest,   r: 'oldrel-1'}
```
