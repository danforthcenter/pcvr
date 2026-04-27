# Function to help fulfill elements of the Bayesian Analysis Reporting Guidelines.

The Bayesian Analysis Reporting Guidelines were put forward by Kruschke
(https://www.nature.com/articles/s41562-021-01177-7) to aide in
reproducibility and documentation of Bayesian statistical analyses that
are sometimes unfamiliar to reviewers or scientists. The purpose of this
function is to summarize goodness of fit metrics from one or more
Bayesian models made by
[growthSS](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
and
[fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md).
See details for explanations of those metrics and the output.

## Usage

``` r
barg(fit, ss = NULL, priors = NULL)
```

## Arguments

- fit:

  A conjugate object, brmsfit object, or a list of brmsfit objects in
  the case that you split models to run on subsets of the data for
  computational simplicity.

- ss:

  The growthSS output used to specify the model. If fit is a list then
  this can either be one growthSS list in which case the priors are
  assumed to be the same for each model or it can be a list of the same
  length as fit. Note that the only parts of this which are used are the
  `call$start` which is expected to be a call, `pcvrForm`, and `df` list
  elements, so if you have a list of brmsfit objects and no ss object
  you can specify a stand-in list. This can also be left NULL (the
  default) and posterior predictive plots and prior predictive plots
  will not be made.

- priors:

  A list of priors similar to how they are specified in conjugate but
  named for the distribution you plan to use, see details and examples.

## Value

A named list containing Rhat, ESS, NEFF, and Trace/Prior/Posterior
Predictive plots. See details for interpretation.

## Details

The majority of the Bayesian Analysis and Reporting Guidelines are
geared towards statistical methods that use MCMC or other numeric
approximations. For those cases (here meaning brms models fit by
`fitGrowth` and `growthSS`) the output will contain:

- **General**: This includes chain number, length, and total divergent
  transitions per model. Divergent transitions are a marker that the
  MCMC had something go wrong. Conceptually it may be helpful to think
  about rolling a marble over a 3D curve then having the marble suddenly
  jolt in an unexpected direction, something happened that suggests a
  problem/misunderstood surface. In practice you want extremely few
  (ideally no) divergences. If you do have divergences then consider
  specifying more control parameters (see brms::brm or examples for
  [fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)).
  If the problem persists then the model may need to be simplified. For
  more information on MCMC and divergence see the stan manual
  (https://mc-stan.org/docs/2_19/reference-manual/divergent-transitions).

- **ESS**: ESS stands for Effective Sample Size and is a goodness of fit
  metric that approximates the number of independent replicates that
  would equate to the same amount of information as the (autocorrelated)
  MCMC iterations. ESS of 1000+ is often considered as a pretty stable
  value, but more is better. Still, 100 per chain may be plenty
  depending on your applications and the inference you wish to do. One
  of the benefits to using lots of chains and/or longer chains is that
  you will get more complete information and that benefit will be shown
  by a larger ESS. This is separated into "bulk" and "tail" to represent
  the middle and tails of the posterior distribution, since those can
  sometimes have very different sampling behavior. A summary and the
  total values are returned, with the summary being useful if several
  models are included in a list for fit argument

- **Rhat**: Rhat is a measure of "chain mixture". It compares the
  between vs within chain values to assess how well the chains mixed. If
  chains did not mix well then Rhat will be greater than 1, with 1.05
  being a broadly agreed upon cutoff to signify a problem. Running
  longer chains should result in lower Rhat values. The default in brms
  is to run 4 chains, partially to ensure that there is a good chance to
  check that the chains mixed well via Rhat. A summary and the total
  values are returned, with the summary being useful if several models
  are included in a list for fit argument

- **NEFF**: NEFF is the NEFF ratio (Effective Sample Size over Total
  MCMC Sample Size). Values greater than 0.5 are generally considered
  good, but there is a consensus that lower can be fine down to about
  0.1. A summary and the total values are returned, with the summary
  being useful if several models are included in a list for fit argument

- **mcmcTrace**: A plot of each model's MCMC traces. Ideally these
  should be very mixed and stationary. For more options for visualizing
  MCMC diagnostics see
  [`bayesplot::mcmc_trace`](https://mc-stan.org/bayesplot/reference/MCMC-traces.html).

- **priorPredictive**: A plot of data simulated from the prior using
  [plotPrior](https://danforthcenter.github.io/pcvr/reference/plotPrior.md).
  This should generate data that is biologically plausible for your
  situation, but it will probably be much more variable than your data.
  That is the effect of the mildly informative thick tailed lognormal
  priors. If you specified non-default style priors then this currently
  will not work.

- **posteriorPredictive**: A plot of each model's posterior predictive
  interval over time. This is the same as plots returned from
  [growthPlot](https://danforthcenter.github.io/pcvr/reference/growthPlot.md)
  and shows 1-99 coming to a mean yellow trend line. These should
  encompass the overwhelming majority of your data and ideally match the
  variance pattern that you see in your data. If parts of the predicted
  interval are biologically impossible (area below 0, percentage about
  100 model should be reconsidered.

For analytic solutions (ie, the `conjugate` class) there are fewer
elements.

- **priorSensitivity**: Patchwork of prior sensitivity plots showing the
  distribution of posterior probabilities, any interpretation changes
  from those tests, and the random priors that were used. This is only
  returned if the `priors` argument is specified (see below).

- **posteriorPredictive**: Plot of posterior predictive distributions
  similar to that from a non-longitudinal `fitGrowth` model fit with
  brms.

- **Summary**: The summary of the `conjugate` object.

Priors here are specified using a named list. For instance, to use 100
normal priors with means between 5 and 20 and standard deviations
between 5 and 10 the prior argument would be
`list("rnorm" = list("mean" = c(5, 20), "sd" = c(5, 10), "n" = 100)))`.
The priors that are used in sensitivity analysis are drawn randomly from
within the ranges specified by the provided list. If you are unsure what
random-generation function to use then check the
[conjugate](https://danforthcenter.github.io/pcvr/reference/conjugate.md)
docs where the distributions are listed for each method in the details
section.

## See also

[plotPrior](https://danforthcenter.github.io/pcvr/reference/plotPrior.md)
for visual prior predictive checks.

## Examples

``` r
# \donttest{
simdf <- growthSim("logistic",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
)
ss <- growthSS(
  model = "logistic", form = y ~ time | id / group, sigma = "logistic",
  df = simdf, start = list(
    "A" = 130, "B" = 12, "C" = 3,
    "sigmaA" = 20, "sigmaB" = 10, "sigmaC" = 2
  ), type = "brms"
)
fit_test <- fitGrowth(ss,
  iter = 600, cores = 1, chains = 1, backend = "cmdstanr",
  sample_prior = "only" # only sampling from prior for speed
)
#> Start sampling
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#> Intercept_nu
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:   1 / 600 [  0%]  (Warmup) 
#> Chain 1 Iteration: 100 / 600 [ 16%]  (Warmup) 
#> Chain 1 Iteration: 200 / 600 [ 33%]  (Warmup) 
#> Chain 1 Iteration: 300 / 600 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 301 / 600 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 400 / 600 [ 66%]  (Sampling) 
#> Chain 1 Iteration: 500 / 600 [ 83%]  (Sampling) 
#> Chain 1 Iteration: 600 / 600 [100%]  (Sampling) 
#> Chain 1 finished in 0.0 seconds.
#> Loading required namespace: rstan
b <- barg(fit_test, ss)
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
fit_2 <- fit_test
fit_list <- list(fit_test, fit_2)
x <- barg(fit_list, list(ss, ss))
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.
#> Warning: The ESS has been capped to avoid unstable estimates.

x <- conjugate(
  s1 = rnorm(10, 10, 1), s2 = rnorm(10, 13, 1.5), method = "t",
  priors = list(
    list(mu = 10, sd = 2),
    list(mu = 10, sd = 2)
  ),
  rope_range = c(-8, 8), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "unequal",
  bayes_factor = c(50, 55)
)
b <- barg(x, priors = list("rnorm" = list("n" = 10, "mean" = c(5, 20), "sd" = c(5, 10))))
# }
```
