# Check priors used in ease of use brms functions

Check priors used in ease of use brms functions

## Usage

``` r
plotPrior(priors, type = "density", n = 200, t = 25)
```

## Arguments

- priors:

  A named list of means for prior distributions. This takes the same
  input as the prior argument of
  [`growthSS`](https://danforthcenter.github.io/pcvr/reference/growthSS.md).
  Alternatively, if given the output of growthSS this will preform a
  prior predictive check and return a plot from
  [`growthPlot`](https://danforthcenter.github.io/pcvr/reference/growthPlot.md)
  of that check ignoring all other arguments. Note that all priors must
  be proper in that case (non-flat) and the fit is likely to be strange
  looking due to how thick tailed the default priors from
  [`growthSS`](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
  are.

- type:

  Either "density", the default, or a model as would be specified in
  `growthSS` or `growthSim` such as "logistic", "gompertz",
  "monomolecular", "exponential", "linear", "power law", "double
  logistic", or "double gompertz". If this is a model type then n draws
  from the prior will be simulated as growth trendlines and densities
  will be plotted on margins for some distributions.

- n:

  Numeric, if type is a model then how many draws from the prior should
  be simulated?

- t:

  Numeric, time passed to growthSim. Defaults to 25 (the growthSim
  default).

## Value

A named list of plots showing prior distributions that `growthSS` would
use, optionally with a plot of simulated growth curves using draws from
those priors.

## See also

[barg](https://danforthcenter.github.io/pcvr/reference/barg.md) for
Bayesian model reporting metrics,
[growthSim](https://danforthcenter.github.io/pcvr/reference/growthSim.md)
for simulating data using similar specification.

## Examples

``` r
set.seed(123)
priors <- list("A" = c(100, 130), "B" = c(10, 8), "C" = c(0.2, 0.1))
plotPrior(priors)
#> $A

#> 
#> $B

#> 
#> $C

#> 

plotPrior(priors, "gompertz")[[1]]

```
