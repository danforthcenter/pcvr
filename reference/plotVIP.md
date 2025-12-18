# Plot Variable Influence on Projection

This function is used to visualize variable influence on projection
(vip) from a plsr model.

## Usage

``` r
plotVIP(plsrObject, i = 1, mean = FALSE, removePattern = ".*_")
```

## Arguments

- plsrObject:

  Output from pcv.plsr

- i:

  An index from the plsrObject to use if the plsrObject contains models
  for several outcomes. Can be a name or a position. Defaults to 1.

- mean:

  Logical, should the mean be plotted (TRUE) or should the components be
  shown individually (FALSE, the default).

- removePattern:

  A pattern to remove to make the wavelength column into a numeric.

## Value

A ggplot showing variable influence on projection

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
  x <- pcv.plsr(df = d, resps = "y", spectra = grepl("^sim_", colnames(d)))
  plotVIP(x)
}
```
