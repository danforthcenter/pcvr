# `pcvr`: An R Package for Plant Phenotyping Statistics ![pcvr Logo](reference/figures/pcvr_logo.png)

`pcvr` provides R functions for use with `PlantCV` output or other
phenotype data with the goal of lowering the barrier to entry for
Bayesian statistics and non-linear modeling.

## Installation

The release version of `pcvr` can be installed from `CRAN`

``` r
install.packages("pcvr")
library(pcvr)
```

Alternatively the development version of `pcvr` can be installed using
remotes/devtools `install_github` as shown below. Note that the default
behavior in devtools/remotes is to only install true dependencies. Some
functions in pcvr use specific packages that would otherwise not be
needed for most work, notably the `brms` modeling functions. To install
suggested packages (see DESCRIPTION file) add `dependencies=T` to the
`install_github` function call.

``` r
devtools::install_github("danforthcenter/pcvr")
library(pcvr)
```

## Vignettes

See the `Vignettes` tab above for several example workflows for common
plant phenotyping tasks.

## Tutorials

See the `Quarto Tutorials` tab above for links to Quarto presentations
on github that go more in depth about the rationale behind several
`pcvr` functions/options.

## Function Reference

Functions are separated by broad goal/type of data they use under the
`Functions` tab above.

## Getting started

Please see the `bellwether` vignette (named for the high throughput
phenotyping facility at the Donald Danforth Plant Science Center) for a
general introduction to `pcvr`.

``` r
vignette("bellwether", package="pcvr")
# or 
browseVignettes("pcvr")
```

## Feedback

Please report bugs and make feature requests with issues the [github
page](https://github.com/danforthcenter/pcvr).
