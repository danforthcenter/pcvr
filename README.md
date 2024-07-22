## pcvr
<!-- badges: start -->
  [![R-CMD-check](https://github.com/danforthcenter/pcvr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/danforthcenter/pcvr/actions/workflows/R-CMD-check.yaml)
  [![Coverage
Status](https://codecov.io/github/danforthcenter/pcvr/coverage.svg?branch=master)](https://app.codecov.io/github/danforthcenter/pcvr)
  <!-- badges: end -->

R functions for use with plantCV output

## Installation

`pcvr` can be installed using remotes/devtools `install_github` as shown below.
Note that the default behavior in devtools/remotes is to only install true dependencies. Several functions in pcvr use specific packages that would otherwise not be needed for most work, notably the brms supporting functions. To install suggested packages (see DESCRIPTION file) add dependencies=T to the `install_github` function call.

```
devtools::install_github("danforthcenter/pcvr", build_vignettes=T) # to install without building vignettes set the build_vignettes argument to F or exclude it.
library(pcvr)
```

## Usage

Please see the `bellwether` vignette for a general introduction to `pcvr`. In the future we expect to have more vignettes and to specialize each to some degree.

```
vignette("bellwether", package="pcvr")
# or 
browseVignettes("pcvr")
```


## Feedback

Please report bugs and make feature requests with issues in this repository.
