# Function to visualize brms survival models specified using growthSS.

Models fit using
[growthSS](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
inputs by
[fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)
(and similar models made through other means) can be visualized easily
using this function. This will generally be called by `growthPlot`.

## Usage

``` r
brmSurvPlot(
  fit,
  form,
  df = NULL,
  groups = NULL,
  timeRange = NULL,
  facetGroups = TRUE
)
```

## Arguments

- fit:

  A brmsfit object, similar to those fit with
  [`growthSS`](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
  outputs.

- form:

  A formula similar to that in `growthSS` inputs specifying the outcome,
  predictor, and grouping structure of the data as
  `outcome ~ predictor|individual/group`.

- df:

  An optional dataframe to use in plotting observed growth curves on top
  of the model.

- groups:

  An optional set of groups to keep in the plot. Defaults to NULL in
  which case all groups in the model are plotted.

- timeRange:

  An optional range of times to use. This can be used to view
  predictions for future data if the available data has not reached some
  point (such as asymptotic size), although prediction using splines
  outside of the observed range is not necessarily reliable.

- facetGroups:

  logical, should groups be separated in facets? Defaults to TRUE.

## Value

Returns a ggplot showing a brms model's credible intervals and
optionally the individual growth lines.

## Examples

``` r
# \donttest{
set.seed(123)
df <- growthSim("exponential",
  n = 20, t = 50,
  params = list("A" = c(1, 1), "B" = c(0.15, 0.1))
)
ss1 <- growthSS(
  model = "survival weibull", form = y > 100 ~ time | id / group,
  df = df, start = c(0, 5)
)
#> Prior is numeric, replicating to 2 length 2 elements (mu, sd) and assuming order a, b
fit1 <- fitGrowth(ss1, iter = 600, cores = 2, chains = 2, backend = "cmdstanr")
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 Iteration:   1 / 600 [  0%]  (Warmup) 
#> Chain 1 Iteration: 100 / 600 [ 16%]  (Warmup) 
#> Chain 1 Iteration: 200 / 600 [ 33%]  (Warmup) 
#> Chain 1 Iteration: 300 / 600 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 301 / 600 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 400 / 600 [ 66%]  (Sampling) 
#> Chain 1 Iteration: 500 / 600 [ 83%]  (Sampling) 
#> Chain 1 Iteration: 600 / 600 [100%]  (Sampling) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 44, column 2 to column 43)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 44, column 2 to column 43)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 44, column 2 to column 43)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 44, column 2 to column 43)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 2 Iteration:   1 / 600 [  0%]  (Warmup) 
#> Chain 2 Iteration: 100 / 600 [ 16%]  (Warmup) 
#> Chain 2 Iteration: 200 / 600 [ 33%]  (Warmup) 
#> Chain 2 Iteration: 300 / 600 [ 50%]  (Warmup) 
#> Chain 2 Iteration: 301 / 600 [ 50%]  (Sampling) 
#> Chain 2 Iteration: 400 / 600 [ 66%]  (Sampling) 
#> Chain 2 Iteration: 500 / 600 [ 83%]  (Sampling) 
#> Chain 2 Iteration: 600 / 600 [100%]  (Sampling) 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 44, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 44, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 44, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: gamma_lpdf: Random variable is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 44, column 2 to column 43)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: weibull_lpdf: Scale parameter[1] is 0, but must be positive finite! (in '/tmp/RtmpVMOOl9/model-2eda2879f52a.stan', line 54, column 4 to column 104)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 1 finished in 0.0 seconds.
#> Chain 2 finished in 0.0 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.0 seconds.
#> Total execution time: 0.2 seconds.
#> 
brmSurvPlot(fit1, form = ss1$pcvrForm, df = ss1$df)


# note that using the cumulative hazard to calculate survival is likely to underestimate
# survival in these plots if events do not start immediately.
ss2 <- growthSS(
  model = "survival binomial", form = y > 100 ~ time | id / group,
  df = df, start = c(-4, 3)
)
#> Prior is numeric, replicating to 2 length 2 elements (mu, sd) and assuming order a, b
#> Priors and parameters are not the same length. Output will assume that priors are for groups and are in order: a, b
fit2 <- fitGrowth(ss2, iter = 600, cores = 2, chains = 2, backend = "cmdstanr")
#> Start sampling
#> Running MCMC with 2 parallel chains...
#> 
#> Chain 1 Iteration:   1 / 600 [  0%]  (Warmup) 
#> Chain 1 Iteration: 100 / 600 [ 16%]  (Warmup) 
#> Chain 1 Iteration: 200 / 600 [ 33%]  (Warmup) 
#> Chain 1 Iteration: 300 / 600 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 301 / 600 [ 50%]  (Sampling) 
#> Chain 2 Iteration:   1 / 600 [  0%]  (Warmup) 
#> Chain 2 Iteration: 100 / 600 [ 16%]  (Warmup) 
#> Chain 2 Iteration: 200 / 600 [ 33%]  (Warmup) 
#> Chain 2 Iteration: 300 / 600 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 400 / 600 [ 66%]  (Sampling) 
#> Chain 1 Iteration: 500 / 600 [ 83%]  (Sampling) 
#> Chain 1 Iteration: 600 / 600 [100%]  (Sampling) 
#> Chain 2 Iteration: 301 / 600 [ 50%]  (Sampling) 
#> Chain 2 Iteration: 400 / 600 [ 66%]  (Sampling) 
#> Chain 2 Iteration: 500 / 600 [ 83%]  (Sampling) 
#> Chain 2 Iteration: 600 / 600 [100%]  (Sampling) 
#> Chain 1 finished in 0.2 seconds.
#> Chain 2 finished in 0.2 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.2 seconds.
#> Total execution time: 0.3 seconds.
#> 
brmSurvPlot(fit2, form = ss2$pcvrForm, df = ss2$df)

# }
```
