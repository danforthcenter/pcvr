# pcvr 1.0.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 1 note ✖

❯ checking CRAN incoming feasibility ... [4s/21s] NOTE
  Maintainer: ‘Josh Sumner <jsumner@danforthcenter.org>’
  
  New submission
  
  Possibly misspelled words in DESCRIPTION:
    Phenotyping (3:14)
    PlantCV (10:179)
    phenotyping (10:44)
  
  Suggests or Enhances not in mainstream repositories:
    cmdstanr
  Availability using Additional_repositories specification:
    cmdstanr   yes   https://mc-stan.org/r-packages/

## Notes

No notes seem critical.
Names (PlantCV) and words (phenotyping) in DESCRIPTION are not misspelled.
The `cmdstanr` suggestion is available from the additional repository specified in DESCRIPTION.
The R CMD check on MacOs via github actions yields a NOTE about installed package size as well.

```
* checking installed package size ... NOTE
  installed size is  6.4Mb
  sub-directories of 1Mb or more:
    doc   4.1Mb
    R     2.0Mb
```

* This is a new release.

## Test environments

```
- {os: macos-latest,   r: 'release'}
- {os: windows-latest, r: 'release'}
- {os: ubuntu-latest,   r: 'devel'}
- {os: ubuntu-latest,   r: 'release'}
- {os: ubuntu-latest,   r: 'oldrel-1'}
```
