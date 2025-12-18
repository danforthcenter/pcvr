# Function to visualize common `quantreg::nlrq` growth models.

Models fit using
[growthSS](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
inputs by
[fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)
(and similar models made through other means) can be visualized easily
using this function. This will generally be called by `growthPlot`.

## Usage

``` r
nlrqPlot(
  fit,
  form,
  df = NULL,
  groups = NULL,
  timeRange = NULL,
  facetGroups = TRUE,
  groupFill = FALSE,
  virMaps = c("plasma")
)
```

## Arguments

- fit:

  A model fit, or list of model fits, returned by `fitGrowth` with
  type="nlrq".

- form:

  A formula similar to that in `growthSS` inputs (or the `pcvrForm` part
  of the output) specifying the outcome, predictor, and grouping
  structure of the data as `outcome ~ predictor|individual/group`. If
  the individual and group are specified then the observed growth lines
  are plotted.

- df:

  A dataframe to use in plotting observed growth curves on top of the
  model. This must be supplied for nlrq models.

- groups:

  An optional set of groups to keep in the plot. Defaults to NULL in
  which case all groups in the model are plotted.

- timeRange:

  An optional range of times to use. This can be used to view
  predictions for future data if the available data has not reached some
  point (such as asymptotic size).

- facetGroups:

  logical, should groups be separated in facets? Defaults to TRUE.

- groupFill:

  logical, should groups have different colors? Defaults to FALSE. If
  TRUE then viridis colormaps are used in the order of virMaps

- virMaps:

  order of viridis maps to use. Will be recycled to necessary length.
  Defaults to "plasma", but will generally be informed by growthPlot's
  default.

## Value

Returns a ggplot showing an nlrq model's quantiles and optionally the
individual growth lines.

## Examples

``` r
simdf <- growthSim("logistic",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
)
ss <- growthSS(
  model = "logistic", form = y ~ time | id / group,
  tau = c(0.5, 0.9), df = simdf, start = NULL, type = "nlrq"
)
#> Individual is not used with type = 'nlrq'.
fit <- fitGrowth(ss)
nlrqPlot(fit, form = ss$pcvrForm, df = ss$df, groups = "a", timeRange = 1:20)
#> Warning: Removed 20 rows containing missing values or values outside the scale range
#> (`geom_line()`).

nlrqPlot(fit, form = ss$pcvrForm, df = ss$df, groupFill = TRUE, virMaps = c("plasma", "viridis"))
#> Warning: Removed 25 rows containing missing values or values outside the scale range
#> (`geom_line()`).
#> Warning: Removed 25 rows containing missing values or values outside the scale range
#> (`geom_line()`).


ss <- growthSS(
  model = "logistic", form = y ~ time,
  tau = c(0.5, 0.9), df = simdf, start = NULL, type = "nlrq"
)
fit <- fitGrowth(ss)
nlrqPlot(fit, form = ss$pcvrForm, df = ss$df)
#> Warning: Removed 25 rows containing missing values or values outside the scale range
#> (`geom_line()`).

```
