# Reduce phenotypes in longitudinal data to cumulative sums of phenotypes.

Often in bellwether experiments we are curious about the effect of some
treatment vs control. For certain routes in analysing the data this
requires considering phenotypes as relative differences compared to a
control.

## Usage

``` r
cumulativePheno(
  df,
  phenotypes = NULL,
  group = "barcode",
  timeCol = "DAS",
  traitCol = "trait",
  valueCol = "value"
)
```

## Arguments

- df:

  Dataframe to use, this can be in long or wide format.

- phenotypes:

  A character vector of column names for the phenotypes that should be
  compared against control.

- group:

  A character vector of column names that identify groups in the data.
  Defaults to "barcode". These groups will be calibrated separately,
  with the exception of the group that identifies a control within the
  greater hierarchy.

- timeCol:

  Column name to use for time data.

- traitCol:

  Column with phenotype names, defaults to "trait". This should
  generally not need to be changed from the default. If this and
  valueCol are present in colnames(df) then the data is assumed to be in
  long format.

- valueCol:

  Column with phenotype values, defaults to "value". This should
  generally not need to be changed from the default.

## Value

A dataframe with cumulative sum columns added for specified phenotypes

## Examples

``` r
# \donttest{
f <- "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv"
tryCatch(
  {
    sv <- read.pcv(
      f,
      reader = "fread"
    )
    sv$genotype <- substr(sv$barcode, 3, 5)
    sv$genotype <- ifelse(sv$genotype == "002", "B73",
      ifelse(sv$genotype == "003", "W605S",
        ifelse(sv$genotype == "004", "MM", "Mo17")
      )
    )
    sv$fertilizer <- substr(sv$barcode, 8, 8)
    sv$fertilizer <- ifelse(sv$fertilizer == "A", "100",
      ifelse(sv$fertilizer == "B", "50", "0")
    )

    sv <- pcv.time(sv,
      plantingDelay = 0, phenotype = "area_pixels", cutoff = 10,
      timeCol = "timestamp", group = c("barcode", "rotation"), plot = TRUE
    )$data
    sv <- pcv.outliers(sv,
      phenotype = "area_pixels", group = c("DAS", "genotype", "fertilizer"),
      plotgroup = c("barcode", "rotation")
    )$data
    phenotypes <- colnames(sv)[19:35]
    phenoForm <- paste0("cbind(", paste0(phenotypes, collapse = ", "), ")")
    groupForm <- "DAS+DAP+barcode+genotype+fertilizer"
    form <- as.formula(paste0(phenoForm, "~", groupForm))
    sv <- aggregate(form, data = sv, mean, na.rm = TRUE)
    pixels_per_cmsq <- 42.5^2 # pixel per cm^2
    sv$area_cm2 <- sv$area_pixels / pixels_per_cmsq
    sv$height_cm <- sv$height_pixels / 42.5
    df <- sv
    phenotypes <- c("area_cm2", "height_cm")
    group <- c("barcode")
    timeCol <- "DAS"
    df <- cumulativePheno(df, phenotypes, group, timeCol)


    sv_l <- read.pcv(
      f,
      mode = "long", reader = "fread"
    )
    sv_l$genotype <- substr(sv_l$barcode, 3, 5)
    sv_l$genotype <- ifelse(sv_l$genotype == "002", "B73",
      ifelse(sv_l$genotype == "003", "W605S",
        ifelse(sv_l$genotype == "004", "MM", "Mo17")
      )
    )
    sv_l$fertilizer <- substr(sv_l$barcode, 8, 8)
    sv_l$fertilizer <- ifelse(sv_l$fertilizer == "A", "100",
      ifelse(sv_l$fertilizer == "B", "50", "0")
    )
    sv_l <- pcv.time(sv_l,
      plantingDelay = 0, phenotype = "area_pixels", cutoff = 10,
      timeCol = "timestamp", group = c("barcode", "rotation")
    )$data
    sv_l <- cumulativePheno(sv_l,
      phenotypes = c("area_pixels", "height_pixels"),
      group = c("barcode", "rotation"), timeCol = "DAS"
    )
  },
  error = function(e) {
    message(e)
  }
)
#> Warning: 14 groupings had all observations removed
# }
```
