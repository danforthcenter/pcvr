# pcvr 1.0.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

❯ checking CRAN incoming feasibility ... [4s/6s] NOTE
  Maintainer: ‘Josh Sumner <jsumner@danforthcenter.org>’
  
  Possibly misspelled words in DESCRIPTION:
    Kruschke (16:2, 17:2, 17:46)
    Phenotyping (3:14)
    phenotyping (12:44)

Suggests or Enhances not in mainstream repositories:
  cmdstanr
Availability using Additional_repositories specification:
  cmdstanr   yes   https://mc-stan.org/r-packages/
❯ checking package namespace information ... OK
❯ checking package dependencies ... NOTE
Package suggested but not available for checking: ‘cmdstanr’

## Notes

No notes seem critical.
Names (PlantCV, Kruschke) and words (phenotyping) in DESCRIPTION are not misspelled.
The `cmdstanr` suggestion is available from the additional repository specified in DESCRIPTION.
The R CMD check on MacOs via github actions yields a NOTE about installed package size as well.

```
* checking installed package size ... NOTE
  installed size is  6.4Mb
  sub-directories of 1Mb or more:
    doc   4.1Mb
    R     2.0Mb
```

Additionally there are some examples wrapped in `\donttest` that may take several minutes to run.

* This is a new release.

## Test environments

```
- {os: macos-latest,   r: 'release'}
- {os: windows-latest, r: 'release'}
- {os: windows-latest, r: 'devel'}
- {os: ubuntu-latest,   r: 'devel'}
- {os: ubuntu-latest,   r: 'release'}
- {os: ubuntu-latest,   r: 'oldrel-1'}
```
