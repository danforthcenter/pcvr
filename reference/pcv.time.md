# Time conversion and plotting for bellwether data

Time conversion and plotting for bellwether data

## Usage

``` r
pcv.time(
  df = NULL,
  mode = c("DAS", "DAP", "DAE"),
  plantingDelay = NULL,
  phenotype = NULL,
  cutoff = 1,
  timeCol = "timestamp",
  group = "Barcodes",
  plot = TRUE,
  format = "%Y-%m-%d %H:%M:%S",
  traitCol = "trait",
  valueCol = "value",
  index = NULL,
  digits = 0
)
```

## Arguments

- df:

  Data frame to use, this can be in wide or long format.

- mode:

  One of "DAS", "DAP" or "DAE" (Days After Planting and Days After
  Emergence). Defaults to adding all columns. Note that if timeCol is
  not numeric then DAS is always returned.

- plantingDelay:

  If \`mode\` includes "DAP" then \`plantingDelay\` is used to adjust
  "DAS"

- phenotype:

  If \`mode\` includes "DAE" then this is the phenotype used to classify
  emergence.

- cutoff:

  If \`mode\` includes "DAE" then this value is used to classify
  emergence. Defaults to 1, meaning an image with a value of 1 or more
  for \`phenotype\` has "emerged".

- timeCol:

  Column of input time values, defaults to "timestamp". If this is not
  numeric then it is assumed to be a timestamp in the format of the
  format argument.

- group:

  Grouping variables to specify unique plants as a character vector.
  This defaults to "Barcodes". These taken together should identify a
  unique plant across time, although often "angle" or "rotation" should
  be added.

- plot:

  Logical, should plots of the new time variables be printed?

- format:

  An R POSIXct format, defaults to lemnatech standard format. This is
  only used if timeCol is not a numeric.

- traitCol:

  Column with phenotype names, defaults to "trait". This should
  generally not need to be changed from the default. If this and
  valueCol are present in colnames(df) then the data is assumed to be in
  long format.

- valueCol:

  Column with phenotype values, defaults to "value". This should
  generally not need to be changed from the default.

- index:

  Optionally a time to use as the beginning of the experiment. This may
  be useful if you have multiple datasets or you are adding data from
  pcv.water and plants were watered before being imaged or if you want
  to index days off of midnight. This defaults to NULL but will take any
  value coercible to POSIXct by `as.POSIXct(... , tz="UTC")` such as
  "2020-01-01 18:30:00"

- digits:

  Number of digits to round DAS to if timeCol is not numeric, defaults
  to 0.

## Value

The input dataframe with new numeric columns for different ways of
describing time in the experiment. If plot is TRUE then a ggplot is also
returned as part of a list.

## Examples

``` r
# \donttest{
f <- "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv"
tryCatch(
  {
    sv <- read.pcv(
      f,
      mode = "wide", reader = "fread"
    )
    sv$genotype = substr(sv$barcode, 3, 5)
    sv$genotype = ifelse(sv$genotype == "002", "B73",
      ifelse(sv$genotype == "003", "W605S",
        ifelse(sv$genotype == "004", "MM", "Mo17")
      )
    )
    sv$fertilizer = substr(sv$barcode, 8, 8)
    sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
      ifelse(sv$fertilizer == "B", "50", "0")
    )
    sv <- pcv.time(sv,
      plantingDelay = 0, phenotype = "area_pixels", cutoff = 10,
      timeCol = "timestamp", group = c("barcode", "rotation"), plot = FALSE
    )

    svl <- read.pcv(
      f,
      mode = "long", reader = "fread"
    )
    svl$genotype = substr(svl$barcode, 3, 5)
    svl$genotype = ifelse(svl$genotype == "002", "B73",
      ifelse(svl$genotype == "003", "W605S",
        ifelse(svl$genotype == "004", "MM", "Mo17")
      )
    )
    svl$fertilizer = substr(svl$barcode, 8, 8)
    svl$fertilizer = ifelse(svl$fertilizer == "A", "100",
      ifelse(svl$fertilizer == "B", "50", "0")
    )
    svl <- pcv.time(svl,
      plantingDelay = 0, phenotype = "area_pixels", cutoff = 10, timeCol = "timestamp",
      group = c("barcode", "rotation"), plot = FALSE
    )
  },
  error = function(e) {
    message(e)
  }
)
# }
```
