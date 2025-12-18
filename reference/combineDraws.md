# Combine Draws From brms Models

Helper function for binding draws from several `brms` models to make a
data.frame for use with
[`brms::hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html).
This will also check that the draws are comparable using basic model
metrics.

## Usage

``` r
combineDraws(..., names = NULL, message = TRUE)
```

## Arguments

- ...:

  Some number of brmsfit objects and/or dataframes of draws (should
  generally be the same type of model fit to different data)

- names:

  Optional vector of names for the models/data.frames. The default
  (NULL) will use the object names of arguments to `...` as a vector of
  names to uniquely identify columns in the output data frame. See
  details for using this with `do.call` on a list of models.

- message:

  Logical, should messages about possible problems be printed? Default
  is TRUE. This will warn if models may not have converged, if there are
  different numbers of draws in the objects, or if models have different
  formulations.

## Value

A data.frame of posterior draws, labeled to show which object they come
from.

Returns a dataframe of posterior draws.

## Details

If you fit models as part of a loop/apply function and end up with a
list of models it may be helpful to call this function on the list. In
that case object names are not parsed well from the list by default so
passing the `names` argument is helpful and can be done as
`do.call(combineDraws, c(fits, list(names = names(fits))))`.

## Examples

``` r
# note that this example will fit several models using Stan and may run slowly.
# \donttest{
simdf <- growthSim("logistic",
  n = 20, t = 25,
  params = list(
    "A" = c(200, 160, 220, 200, 140, 300),
    "B" = c(13, 11, 10, 9, 16, 12),
    "C" = c(3, 3.5, 3.2, 2.8, 3.3, 2.5)
  )
)
ss_ab <- growthSS(
  model = "logistic", form = y ~ time | id / group,
  sigma = "logistic", df = simdf[simdf$group %in% c("a", "b"), ],
  start = list(
    "A" = 130, "B" = 12, "C" = 3,
    "sigmaA" = 15, "sigmaB" = 10, "sigmaC" = 3
  ), type = "brms"
)

ss_cd <- growthSS(
  model = "logistic", form = y ~ time | id / group,
  sigma = "logistic", df = simdf[simdf$group %in% c("c", "d"), ],
  start = list(
    "A" = 130, "B" = 12, "C" = 3,
    "sigmaA" = 15, "sigmaB" = 10, "sigmaC" = 3
  ), type = "brms"
)

ss_ef <- growthSS(
  model = "logistic", form = y ~ time | id / group,
  sigma = "logistic", df = simdf[simdf$group %in% c("e", "f"), ],
  start = list(
    "A" = 130, "B" = 12, "C" = 3,
    "sigmaA" = 15, "sigmaB" = 10, "sigmaC" = 3
  ), type = "brms"
)
ss_ef2 <- growthSS(
  model = "gompertz", form = y ~ time | id / group,
  sigma = "logistic", df = simdf[simdf$group %in% c("e", "f"), ],
  start = list(
    "A" = 130, "B" = 12, "C" = 3,
    "sigmaA" = 15, "sigmaB" = 10, "sigmaC" = 3
  ), type = "brms"
)


fit_ab <- fitGrowth(ss_ab, chains = 1, cores = 1, iter = 1000)
#> Start sampling
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#> Intercept_nu
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:   1 / 1000 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is -nan, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is -nan, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[501] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Degrees of freedom parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is -nan, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration: 100 / 1000 [ 10%]  (Warmup) 
#> Chain 1 Iteration: 200 / 1000 [ 20%]  (Warmup) 
#> Chain 1 Iteration: 300 / 1000 [ 30%]  (Warmup) 
#> Chain 1 Iteration: 400 / 1000 [ 40%]  (Warmup) 
#> Chain 1 Iteration: 500 / 1000 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 501 / 1000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 600 / 1000 [ 60%]  (Sampling) 
#> Chain 1 Iteration: 700 / 1000 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 800 / 1000 [ 80%]  (Sampling) 
#> Chain 1 Iteration: 900 / 1000 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1000 [100%]  (Sampling) 
#> Chain 1 finished in 9.7 seconds.
fit_ab2 <- fitGrowth(ss_ab, chains = 1, cores = 1, iter = 1200)
#> Start sampling
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#> Intercept_nu
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:    1 / 1200 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[501] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[501] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[501] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[2] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Degrees of freedom parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is -nan, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[501] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration:  100 / 1200 [  8%]  (Warmup) 
#> Chain 1 Iteration:  200 / 1200 [ 16%]  (Warmup) 
#> Chain 1 Iteration:  300 / 1200 [ 25%]  (Warmup) 
#> Chain 1 Iteration:  400 / 1200 [ 33%]  (Warmup) 
#> Chain 1 Iteration:  500 / 1200 [ 41%]  (Warmup) 
#> Chain 1 Iteration:  600 / 1200 [ 50%]  (Warmup) 
#> Chain 1 Iteration:  601 / 1200 [ 50%]  (Sampling) 
#> Chain 1 Iteration:  700 / 1200 [ 58%]  (Sampling) 
#> Chain 1 Iteration:  800 / 1200 [ 66%]  (Sampling) 
#> Chain 1 Iteration:  900 / 1200 [ 75%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1200 [ 83%]  (Sampling) 
#> Chain 1 Iteration: 1100 / 1200 [ 91%]  (Sampling) 
#> Chain 1 Iteration: 1200 / 1200 [100%]  (Sampling) 
#> Chain 1 finished in 10.0 seconds.
fit_cd <- fitGrowth(ss_cd, chains = 1, cores = 1, iter = 1000)
#> Start sampling
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#> Intercept_nu
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:   1 / 1000 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is -nan, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is -nan, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[501] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[501] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Degrees of freedom parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration: 100 / 1000 [ 10%]  (Warmup) 
#> Chain 1 Iteration: 200 / 1000 [ 20%]  (Warmup) 
#> Chain 1 Iteration: 300 / 1000 [ 30%]  (Warmup) 
#> Chain 1 Iteration: 400 / 1000 [ 40%]  (Warmup) 
#> Chain 1 Iteration: 500 / 1000 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 501 / 1000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 600 / 1000 [ 60%]  (Sampling) 
#> Chain 1 Iteration: 700 / 1000 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 800 / 1000 [ 80%]  (Sampling) 
#> Chain 1 Iteration: 900 / 1000 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1000 [100%]  (Sampling) 
#> Chain 1 finished in 8.5 seconds.
fit_ef <- fitGrowth(ss_ef, chains = 1, cores = 1, iter = 1000)
#> Start sampling
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#> Intercept_nu
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:   1 / 1000 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is -nan, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is -nan, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[501] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Location parameter[1] is -nan, but must be finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is -nan, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[501] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff229b39e9.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration: 100 / 1000 [ 10%]  (Warmup) 
#> Chain 1 Iteration: 200 / 1000 [ 20%]  (Warmup) 
#> Chain 1 Iteration: 300 / 1000 [ 30%]  (Warmup) 
#> Chain 1 Iteration: 400 / 1000 [ 40%]  (Warmup) 
#> Chain 1 Iteration: 500 / 1000 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 501 / 1000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 600 / 1000 [ 60%]  (Sampling) 
#> Chain 1 Iteration: 700 / 1000 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 800 / 1000 [ 80%]  (Sampling) 
#> Chain 1 Iteration: 900 / 1000 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1000 [100%]  (Sampling) 
#> Chain 1 finished in 8.6 seconds.
fit_ef2 <- fitGrowth(ss_ef2, chains = 1, cores = 1, iter = 1000)
#> Start sampling
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#> Intercept_nu
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:   1 / 1000 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff543cba8e.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff543cba8e.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff543cba8e.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff543cba8e.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Location parameter[1] is -nan, but must be finite! (in '/tmp/RtmpfB9aQg/model-1eff543cba8e.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff543cba8e.stan', line 120, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration: 100 / 1000 [ 10%]  (Warmup) 
#> Chain 1 Iteration: 200 / 1000 [ 20%]  (Warmup) 
#> Chain 1 Iteration: 300 / 1000 [ 30%]  (Warmup) 
#> Chain 1 Iteration: 400 / 1000 [ 40%]  (Warmup) 
#> Chain 1 Iteration: 500 / 1000 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 501 / 1000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 600 / 1000 [ 60%]  (Sampling) 
#> Chain 1 Iteration: 700 / 1000 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 800 / 1000 [ 80%]  (Sampling) 
#> Chain 1 Iteration: 900 / 1000 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1000 [100%]  (Sampling) 
#> Chain 1 finished in 21.5 seconds.

x <- combineDraws(fit_ab, fit_cd, fit_ef)
draws_ef <- as.data.frame(fit_ef)
draws_ef <- draws_ef[, grepl("^b_", colnames(draws_ef))]
x2 <- combineDraws(fit_ab2, fit_cd, draws_ef)
#> fit_cd has fewer than 600 draws and will be padded with 100 NAs
#> draws_ef has fewer than 600 draws and will be padded with 100 NAs
x3 <- combineDraws(fit_ab, fit_cd, fit_ef2)
#> Some of these models have different growth formulas, consider if this is what you want.
#> fit_ab: y~A/(1 + exp((B - time)/C)), fit_cd: y~A/(1 + exp((B - time)/C)), fit_ef2: y~A * exp(-B * exp(-C * time))
# }
```
