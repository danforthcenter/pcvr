# Run Partial Least Squares Regression on spectral data

Partial Least Squares Regression (plsr) is often used to analyze
spectral data.

## Usage

``` r
pcv.plsr(df, resps = NULL, spectra = NULL, train = 0.8, cv = 10, ...)
```

## Arguments

- df:

  Data frame containing metadata and spectral histogram data

- resps:

  Vector of response variables.

- spectra:

  Either one column name (in the case of long data) or a set of columns
  in the case of wide data. If a single character string is provided and
  it is not one of the column names then it is taken to be a pattern
  that will match some set of column names in the data to use (see
  examples).

- train:

  Proportion of data to use as training data.

- cv:

  Number of cross validation iterations.

- ...:

  Further arguments passed to caret::train.

## Value

a list of lists each with model performance, prediction target, model,
plot, N components, and variable influence on projection components for
each response variable.

## Details

Note that columns that sum to 0 in the training or test data will be
removed. This function also uses the 'pls' method from the pls package.

## Examples

``` r
if (rlang::is_installed("pls")) {
  dists <- list(
    rlnorm = list(meanlog = log(40), sdlog = 0.5),
    rlnorm = list(meanlog = log(60), sdlog = 0.35)
  )
  mv <- mvSim(
    dists = dists, n_samples = 100, counts = 1000,
    min_bin = 1, max_bin = 180, wide = TRUE
  )
  sv <- growthSim("logistic",
    n = 5, t = 20,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  d <- cbind(sv, mv[, -1])
  # note that this requires the "pls" package to be installed.
  x <- pcv.plsr(df = d, resps = "y", spectra = grepl("^sim_", colnames(d)))
}
```
