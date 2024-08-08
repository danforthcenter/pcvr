# pcvr 1.0.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

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

❯ checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found


## Notes

No notes seem critical.
Names (PlantCV) and words (phenotyping) in DESCRIPTION are not misspelled.
The `cmdstanr` suggestion is available from the additional repository speficied in DESCRIPOTION.
The HTML version of manual works on other systems.

* This is a new release.

## Test environments

```
- {os: macos-latest,   r: 'release'}
- {os: windows-latest, r: 'release'}
- {os: ubuntu-latest,   r: 'devel'}
- {os: ubuntu-latest,   r: 'release'}
- {os: ubuntu-latest,   r: 'oldrel-1'}
```
