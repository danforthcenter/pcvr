# Read in lemnatech watering data from metadata.json files

Read in lemnatech watering data from metadata.json files

## Usage

``` r
pcv.water(file = NULL, envKey = "environment")
```

## Arguments

- file:

  Path to a json file of lemnatech metadata.

- envKey:

  Character string representing the json key for environment data. By
  default this is set to "environment". Currently there are no
  situations where this makes sense to change.

## Value

A data frame containing the bellwether watering data

## Examples

``` r
tryCatch(
  {
    w <- pcv.water("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/metadata.json")
  },
  error = function(e) {
    message(e)
  }
)
#> Using the first watering time, 2023-04-13 23:28:17.58, as beginning of experiment to assign DAS
```
