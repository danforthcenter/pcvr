# Function to visualize common `nlme::nlme` growth models.

Models fit using
[growthSS](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
inputs by
[fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)
(and similar models made through other means) can be visualized easily
using this function. This will generally be called by `growthPlot`.

## Usage

``` r
nlmePlot(
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

  A model fit returned by `fitGrowth` with type="nlme".

- form:

  A formula similar to that in `growthSS` inputs (or the `pcvrForm` part
  of the output) specifying the outcome, predictor, and grouping
  structure of the data as `outcome ~ predictor|individual/group`

- df:

  A dataframe to use in plotting observed growth curves on top of the
  model. This must be supplied for nlme models.

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
  TRUE then viridis colormaps are used in the order of virMaps.

- virMaps:

  order of viridis maps to use. Will be recycled to necessary length.
  Defaults to "plasma", but will generally be informed by growthPlot's
  default.

## Value

Returns a ggplot showing an nlme model's credible intervals and
optionally the individual growth lines.

## Examples

``` r
simdf <- growthSim("logistic",
  n = 10, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
)

ss <- growthSS(
  model = "logistic", form = y ~ time | id / group, sigma = "none",
  df = simdf, start = NULL, type = "nlme"
)

fit <- fitGrowth(ss)

nlmePlot(fit, form = ss$pcvrForm, groups = NULL, df = ss$df, timeRange = NULL)

nlmePlot(fit, form = ss$pcvrForm, groups = "a", df = ss$df, timeRange = 1:10, groupFill = TRUE)

```
