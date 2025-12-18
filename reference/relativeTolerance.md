# Calculate relative tolerance of some phenotype(s) relative to control

Often in bellwether experiments we are curious about the effect of some
treatment vs control. For certain routes in analysing the data this
requires considering phenotypes as relative differences compared to a
control. Note that the `conjugate` function can also be useful in
considering the relative tolerance to stress between groups and that
growth models are another suggested way to test relative tolerance
questions.

## Usage

``` r
relativeTolerance(
  df,
  phenotypes = NULL,
  grouping = NULL,
  control = NULL,
  controlGroup = NULL,
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

- grouping:

  A character vector of column names that identify groups in the data.
  These groups will be calibrated separately, with the exception of the
  group that identifies a control within the greater hierarchy. Note
  that for levels of grouping where the control group does not exist the
  output will be NA.

- control:

  A column name for the variable to be used to select the control
  observations. If left NULL (the default) then this will be taken as
  the first string in the group argument.

- controlGroup:

  The level of the control variable to compare groups against.

- traitCol:

  Column with phenotype names, defaults to "trait". This should
  generally not need to be changed from the default. If this and
  valueCol are present in colnames(df) then the data is assumed to be in
  long format.

- valueCol:

  Column with phenotype values, defaults to "value". This should
  generally not need to be changed from the default.

## Value

A dataframe with relative tolerance columns added.

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
      plantingDelay = 0, phenotype = "area_pixels",
      cutoff = 10, timeCol = "timestamp", group = c("barcode", "rotation"), plot = FALSE
    )
    phenotypes <- colnames(sv)[19:35]
    phenoForm <- paste0("cbind(", paste0(phenotypes, collapse = ", "), ")")
    groupForm <- "DAS+DAP+barcode+genotype+fertilizer"
    form <- as.formula(paste0(phenoForm, "~", groupForm))
    sv <- aggregate(form, data = sv, mean, na.rm = TRUE)
    sv <- pcv.outliers(sv,
      phenotype = "area_pixels",
      group = c("DAS", "genotype", "fertilizer"),
      plotgroup = c("barcode")
    )$data

    pixels_per_cmsq <- 42.5^2 # pixel per cm^2
    sv$area_cm2 <- sv$area_pixels / pixels_per_cmsq
    sv$height_cm <- sv$height_pixels / 42.5

    df <- sv
    phenotypes <- c("area_cm2", "height_cm")
    grouping <- c("fertilizer", "genotype", "DAS")
    controlGroup <- "100"
    control <- "fertilizer"

    rt <- relativeTolerance(df, phenotypes, grouping, control, controlGroup)
    head(rt)
    sapply(rt, function(c) sum(is.na(c)))
  },
  error = function(e) {
    message(e)
  }
)
#> Warning: 14 groupings had all observations removed
#> fertilizer   genotype        DAS  phenotype     mu_rel     se_rel     mu_trt 
#>          0          0          0          0         20         42          0 
#>     se_trt mu_control se_control 
#>         20         20         32 
# }
```
