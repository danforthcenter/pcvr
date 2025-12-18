# Ease of use growth model helper function.

Output from this should be passed to
[fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)
to fit the specified model.

## Usage

``` r
growthSS(
  model,
  form,
  sigma = NULL,
  df,
  start = NULL,
  pars = NULL,
  type = "brms",
  tau = 0.5,
  hierarchy = NULL
)
```

## Arguments

- model:

  The name of a model as a character string. Supported options are
  c("logistic", "logistic4", "logistic5", "gompertz", "weibull",
  "frechet", "gumbel", "monomolecular", "exponential", "linear", "power
  law", "bragg", "lorentz", "beta", "double logistic", "double
  gompertz", "gam", "int"), with "int" representing an intercept only
  model which is only used in brms (and is expected to only be used in
  threshold models or to model homoskedasticity). Note that the dose
  response curves (bragg, lorentz, and beta) may be difficult to fit
  using the `nlme` backend but should work well using other options. See
  [`growthSim`](https://danforthcenter.github.io/pcvr/reference/growthSim.md)
  for examples of each type of single parameterized growth curve ("gam"
  is not supported in `growthSim`). You can also specify decay models by
  including the "decay" keyword with the model name. Note that using
  "decay" is only necessary for the brms backend since otherwise the
  priors are strictly positive. In brms models the entire formula is
  negated for decay models so that lognormal priors can still be used
  when at least some coefficients would be negative. Additionally, the
  "int\_" prefix may be added to a model name to specify that an
  intercept should be included. By default these models are assumed to
  have intercepts at 0, which is often fine. If you include an intercept
  in a brms model then you would specify the prior as you would for an
  "A", "B", or "C" parameter but as "I". By default growthSS will make
  student T priors for intercept parameters in the same way that it will
  for estimated changepoints (see below). With type="brms" you can also
  specify segmented models by combining model names with a plus sign
  such as "linear + linear". In a segmented model the names for
  parameters do not follow the normal "A", "B", "C" notation, instead
  they are named for the type of model, the position in the formula,
  then for the parameter of that model. There will also be parameters to
  represent the time when growth switches from one model to another
  called either "changepointX" or "fixedChangePointX". All
  "changePointX" terms are estimated as parameters of the model.
  "fixedChangePointX" parameters are not estimated and are kept as the
  numeric value given in the priors, this is useful if your experiment
  has an intervention at a set time which you expect to change the
  growth process acutely. For the "linear + linear" example this would
  yield parameters "linear1A", "changePoint1" (or "fixedChangePoint1"),
  and "linear2A". A "linear + gompertz" model would have "linear1A",
  "changePoint1", "gompertz2A", "gompertz2B", and "gompertz2C" for
  parameters. Note that double sigmoid models are not supported as parts
  of segmented models and gams can currently only be included as the
  last part of a segmented model. When using a changepoint model it may
  be worth using segments that are simpler to fit (gompertz instead of
  EVD options, for instance). Currently "homo" and "int" are treated the
  same and "spline" and "gam" are interchangeable. Time-to-event models
  can be specified using the "survival" keyword, see details for an
  explanation of the changes that entails. Similarly, using the brms
  backend response distributions (see
  [`brms::brmsfamily`](https://paulbuerkner.com/brms/reference/brmsfamily.html))
  can be specified in the model as "family: model" so that a model of
  logistic increasing counts may be written as
  `model = "poisson: logistic"`.

- form:

  A formula describing the model. The left hand side should only be the
  outcome variable (phenotype), and a cutoff if you are making a
  survival model (see details). The right hand side needs at least the x
  variable (which should be numeric and is typically time). Grouping is
  also described in this formula using roughly lme4 style syntax,with
  formulas like `y~time|individual/group` to show that predictors should
  vary by `group` and autocorrelation between `individual:group`
  interactions should be modeled. Note that autocorrelation is only
  modeled with the "brms" backend in this way. "nlme" requires random
  effects and correlations to use the same grouping, so autocorrelation
  using the "nlme" backend works at the group level, so will slightly
  underestimate the autocorrelation at the individual level. If group
  has only one level or is not included then it will be ignored in
  formulas for growth and variance (this may be the case if you split
  data before fitting models to be able to run more smaller models each
  more quickly). To include multiple grouping variables they should be
  separated with "+" as in
  `y~time|individual/groupingVariable1 + groupingVariable2`. For some
  backends multiple grouping variables will be combined into a single
  factor of their interaction. Hierarchical models can be specified for
  the brms backend as `y~time+other_covariate|individual/group` in which
  case the parameters of the main growth model will themselves be
  estimated by models as specified in the `hierarchy` argument. For
  instance, if normally "A" had an intercept for each `group`, now it
  would be predicted as `A ~ AI + AA * covariate` where AI and AA now
  have an intercept for each `group`. Note that if you specify a
  hierarchical model then priors are required for AI and AA in the
  previous example.

- sigma:

  Other models for distributional parameters. This argument is only used
  with "brms" and "nlme" models and is handled differently for each.
  When type="brms" this can be supplied as a model or as a list of
  models. It is turned into a formula (or list of formulas) with an
  entry corresponding to each distributional parameter (after the mean)
  of the growth model family. If no family was specified
  (`model="logistic"` for instance) then the student T distribution is
  used, with additional distributional parameters sigma and nu. To check
  the naming of distributional parameters in each response family use
  `brms::brmsfamily("family")$dpars`. The supported options are the same
  as the model options (including threshold models). For distributional
  parameters that do not have a formula specified they will be modeled
  as intercept only (not by group). Parameter names are the same as
  those in the main model but with the distributional parameter name as
  a prefix. Additionally, if a linear model is used for sigma then it
  can be modeled with or without a prior, if a prior is specified
  ("sigmaA") then a non-linear formula is used and the "sigmaA"
  parameter will be included in the output instead of the default
  "sigma" term. In the rare case that you wish to model the mean and the
  3rd distributional parameter but not the 2nd then
  `sigma = list("not_estimated", "model")` would allow for that. When
  type ="nlme" the options are more limited to c("none", "power",
  "exp"), corresponding to using
  [`nlme::varIdent`](https://rdrr.io/pkg/nlme/man/varIdent.html),
  [`nlme::varPower`](https://rdrr.io/pkg/nlme/man/varPower.html), or
  [`nlme::varExp`](https://rdrr.io/pkg/nlme/man/varExp.html)
  respectively where "power" is the default.

- df:

  A dataframe to use. Must contain all the variables listed in the
  formula. Note that rows with NA or infinite values in x, y, or
  hierarchical predictors are removed.

- start:

  An optional named list of starting values OR means for prior
  distributions. If this is not provided then starting values are picked
  with [`stats::selfStart`](https://rdrr.io/r/stats/selfStart.html).
  When type = "brms" these should be provided and are treated as the
  means of lognormal priors for all growth model parameters and
  T_5(mu, 3) priors for changepoint parameters. This is done because the
  values are strictly positive and the lognormal distribution is easily
  interpreted. The changepoint priors are T distributions for symmetry,
  5 DF having been chosen for heavy but not unmanageable tails. If this
  argument is not provided then priors are made using brms::get_prior.
  Those priors are unlikely to be suitable and a different set of priors
  will need to be made for the model using
  [`brms::set_prior`](https://paulbuerkner.com/brms/reference/set_prior.html)
  for good convergence. When specifying starting values/prior means
  think of this as being similar to the `params` argument in
  `growthSim`. Names should correspond to parameter names from the
  `model` argument. A numeric vector can also be used, but specifying
  names is best practice for clarity. Additionally, due to a limitation
  in `brms` currently lower bounds cannot be set for priors for specific
  groups. If priors include multiple groups
  (`start = list(A = c(10,15), ...)`) then you will see warnings after
  the model is fit about not having specified a lower bound explicitly.
  Those warnings can safely be ignored and will be addressed if the
  necessary features are added to `brms`. See details for guidance.

- pars:

  Optionally specify which parameters should change by group. By default
  all parameters are modeled by group.

- type:

  Type of model to fit, options are "brms", "nlrq", "nlme", "nls", and
  "mgcv". Note that the "mgcv" option only supports "gam" models.
  Survival models can use the "survreg" model type (this will be called
  if any non-brms/flexsurv type is given) or the "flexsurv" model type
  which requires the flexsurv package to be installed. Note that for
  non-brms models variables in the model will be labeled by the factor
  level of the group, not necessarily by the group name This is done for
  ease of use with different modeling functions, the levels are
  alphabetically sorted and can be checked using:
  `table(ss$df$group, ss$df$group_numericLabel)`.

- tau:

  A vector of quantiles to fit for nlrq models.

- hierarchy:

  Optionally a list of model parameters that should themselves by
  modeled by another predictor variable. This is only used with the brms
  backend.

## Value

A named list of elements to make it easier to fit non linear growth
models with several R packages.

For `brms` models the output contains:

`formula`: A
[`brms::bf`](https://paulbuerkner.com/brms/reference/brmsformula.html)
formula specifying the growth model, autocorrelation, variance submodel,
and models for each variable in the growth model. `prior`: A
brmsprior/data.frame object. `initfun`: A function to randomly
initialize chains using a random draw from a gamma distribution
(confines initial values to positive and makes correct number of initial
values for chains and groups). `df` The data input, with dummy variables
added if needed and a column to link groups to their factor levels.
`family` The model family, currently this will always be "student".
`pcvrForm` The form argument unchanged. This is returned so that it can
be used later on in model visualization. Often it may be a good idea to
save the output of this function with the fit model, so having this can
be useful later on.

For [`quantreg::nlrq`](https://rdrr.io/pkg/quantreg/man/nlrq.html)
models the output contains:

`formula`: An `nls` style formula specifying the growth model with
groups if specified. `taus`: The quantiles to be fit `start`: The
starting values, typically these will be generated from the growth model
and your data in a similar way as shown in
[`stats::selfStart`](https://rdrr.io/r/stats/selfStart.html) models.
`df` The input data for the model. `pcvrForm` The form argument
unchanged.

For `nls` models the output is the same as for
[`quantreg::nlrq`](https://rdrr.io/pkg/quantreg/man/nlrq.html) models
but without `taus` returned.

For [`nlme::nlme`](https://rdrr.io/pkg/nlme/man/nlme.html) models the
output contains:

`formula`: An list of `nlme` style formulas specifying the model, fixed
and random effects, random effect grouping, and variance model
(weights). `start`: The starting values, typically these will be
generated from the growth model and your data in a similar way as shown
in [`stats::selfStart`](https://rdrr.io/r/stats/selfStart.html) models.
`df` The input data for the model. `pcvrForm` The form argument
unchanged.

For all models the type and model are also returned for simplicity
downstream.

## Details

Default priors are not provided, but these can serve as starting points
for each distribution. You are encouraged to use `growthSim` to consider
what kind of trendlines result from changes to your prior and for
interpretation of each parameter. The
[plotPrior](https://danforthcenter.github.io/pcvr/reference/plotPrior.md)
function can be used to do prior predictive checks. You should not
looking back and forth at your data trying to match your observed growth
exactly with a prior distribution, rather this should be informed by an
understanding of the plants you are using and expectations based on
previous research. For the "double" models the parameter interpretation
is the same as for their non-double counterparts except that there are A
and A2, etc. It is strongly recommended to familiarize yourself with the
double sigmoid distributions using growthSim before attempting to model
one. Additionally, those distributions are intended for use with long
delays in an experiment, think stress recovery experiments, not for
minor hiccups in plant growth.

- **Logistic**: `list('A' = 130, 'B' = 12, 'C' = 3)`

- **Logistic4**: `list('A' = 130, 'B' = 12, 'C' = 3, 'D' = 0)`

- **Logistic5**: `list('A' = 130, 'B' = 12, 'C' = 3, 'D' = 0, 'E' = 1)`

- **Gompertz**: `list('A' = 130, 'B' = 12, 'C' = 1.25)`

- **Weibull**: `list('A' = 130, 'B' = 2, 'C' = 2)`

- **Frechet**: `list('A' = 130, 'B' = 5, 'C' = 6)`

- **Gumbel**: `list('A' = 130, 'B' = 6, 'C' = 4)`

- **Double Logistic**:
  `list('A' = 130, 'B' = 12, 'C' = 3, 'A2' = 200, 'B2' = 25, 'C2' = 1)`

- **Double Gompertz**:
  `list('A' = 130, 'B' = 12, 'C' = 0.25, 'A2' = 220, 'B2' = 30, 'C2' = 0.1)`

- **Monomolecular**: `list('A' = 130, 'B' = 2)`

- **Exponential**: `list('A' = 15, 'B' = 0.1)`

- **Linear**: `list('A' = 1)`

- **Power Law**: `list('A' = 13, 'B' = 2)`

See details below about parameterization for each model option.

- **Logistic**: \`A / (1 + exp( (B-x)/C) )\` Where A is the asymptote, B
  is the inflection point, C is the growth rate.

- **Logistic4**: \`D + (A - D) / (1 + exp((B - x) / C))\` Where A is the
  asymptote (maximum), B is the inflection point, C is the growth rate,
  and D is the lower asymptote (minimum, if this is 0 then the model
  converges to 3 parameter logistic).

- **Logistic5**: \`D + ((A - D) / (1 + exp((B - x) / C)) ^ E)\` Where A
  is the asymptote (maximum), B is the inflection point, C is the growth
  rate, and D is the lower asymptote (minimum), and E is the asymmetry
  factor (1 is symmetric and converges to 4 parameter logistic).

- **Gompertz**: \`A \* exp(-B \* exp(-C\*x))\` Where A is the asymptote,
  B is the inflection point, C is the growth rate.

- **Weibull**: \`A \* (1-exp(-(x/C)^B))\` Where A is the asymptote, B is
  the weibull shape parameter, C is the weibull scale parameter.

- **Frechet**: \`A \* exp(-((x-0)/C)^(-B))\` Where A is the asymptote, B
  is the frechet shape parameter, C is the frechet scale parameter. Note
  that the location parameter (conventionally m) is 0 in these models
  for simplicity but is still included in the formula.

- **Gumbel**: \`A \* exp(-exp(-(x-B)/C))\` Where A is the asymptote, B
  is the inflection point (location), C is the growth rate (scale).

- **Double Logistic**: \`A / (1+exp((B-x)/C)) + ((A2-A)
  /(1+exp((B2-x)/C2)))\` Where A is the asymptote, B is the inflection
  point, C is the growth rate, A2 is the second asymptote, B2 is the
  second inflection point, and C2 is the second growth rate.

- **Double Gompertz**: \`A \* exp(-B \* exp(-C\*x)) + ((A2-A) \* exp(-B2
  \* exp(-C2\*(x-B))))\` Where A is the asymptote, B is the inflection
  point, C is the growth rate, A2 is the second asymptote, B2 is the
  second inflection point, and C2 is the second growth rate.

- **Monomolecular**: \`A-A \* exp(-B \* x)\` Where A is the asymptote
  and B is the growth rate.

- **Exponential**: \`A \* exp(B \* x)\` Where A is the scale parameter
  and B is the growth rate.

- **Linear**: \`A \* x\` Where A is the growth rate.

- **Power Law**: \`A \* x^(B)\` Where A is the scale parameter and B is
  the growth rate.

- **Bragg**: \`A \* exp(-B \* (x - C) ^ 2)\` This models minima and
  maxima as a dose-response curve where A is the max response, B is the
  "precision" or slope at inflection, and C is the x position of the max
  response.

- **Lorentz**: \`A / (1 + B \* (x - C) ^ 2)\` This models minima and
  maxima as a dose-response curve where A is the max response, B is the
  "precision" or slope at inflection, and C is the x position of the max
  response. Generally Bragg is preferred to Lorentz for dose-response
  curves.

- **Beta**: \`A \* (((x - D) / (C - D)) \* ((E - x) / (E - C)) ^
  ((E - C) / (C - D))) ^ B\` This models minima and maxima as a
  dose-response curve where A is the Maximum Value, B is a
  shape/concavity exponent similar to the sum of alpha and beta in a
  Beta distribution, C is the position of maximum value, D is the
  minimum position where distribution \> 0, E is the maximum position
  where distribution \> 0. This is a difficult model to fit but can
  model non-symmetric dose-response relationships which may sometimes be
  worth the extra effort.

Note that for these distributions parameters do not exist in a vacuum.
Changing one will make the others look different in the resulting data.
The `growthSim` function can be helpful in familiarizing further with
these distributions.

Using the `brms` backend the `sigma` argument optionally specifies a sub
model to account for heteroskedasticity. Any combination of models
(except for decay models) can be specified in the `sigma` term. If you
need variance to raise and lower then a gam/spline is the most
appropriate option.

Using the `brms` backend a model with lots of parameters may be
difficult to estimate if there are lots of groups. If you have very many
levels of your "group" variable in a complex model then consider fitting
models to subsets of the "group" variable and using
[combineDraws](https://danforthcenter.github.io/pcvr/reference/combineDraws.md)
to make a data.frame for hypothesis testing.

Limits on the Y variable can be specified in the `brms` backend. This
should generally be unnecessary and will make the model slower to fit
and potentially more difficult to set priors on. If you do have a
limited phenotype (besides the normal positive constraint for growth
models) then this may be helpful, one situation may be canopy coverage
percentage which is naturally bounded at an upper and lower limit. To
specify these limits add square brackets to the Y term with upper and
lower limits such as `"y[0,100] ~ time|id/group"`. Other "Additional
response information" such as resp_weights or standard errors can be
specified using the `brms` backend, with those options documented fully
in the
[`brms::brmsformula`](https://paulbuerkner.com/brms/reference/brmsformula.html)
details.

There are also three supported submodel options for `nlme` models, but a
`varFunc` object can also be supplied, see
[`?nlme::varClasses`](https://rdrr.io/pkg/nlme/man/varClasses.html).

- **none**: `varIdent(1|group)`, which models a constant variance
  separately for each group.

- **power**: `varPower(x|group)`, which models variance as a power of x
  per group.

- **exp**: `varExp(x|group)`, which models variance as an exponent of x
  per group.

Survival models can be fit using the "survival" keyword in the model
specification. Using the "brms" backend (type argument) you can specify
"weibull" (the default) or "binomial" for the distribution to use in
that model so that the final model string would be "survival binomial"
or "survival weibull" which is equivalent to "survival". Time to event
data is very different than standard phenotype data, so the formula
argument should include a cutoff for the Y variable to count as an
"event". For example, if you were checking germination using area and
wanted to use 50 pixels as a germinated plant your formula would be
`area > 50 ~ time|id/group`. Internally the input dataframe will be
converted to time-to-event data based on that formula. Alternatively you
can make your own time to event data and supply that to growthSS. In
that case your data should have columns called "n_events" (number of
individuals experiencing the event at this time) and "n_eligible"
(number of individuals who had not experienced the event at least up to
this time) for the binomial model family OR "event" (binary 1,0 for
TRUE, FALSE) for the Weibull model family. Note that since these are
linear models using different model families the priors are handled
differently. For survival models the default priors are weak
regularizing priors (Normal(0,5)) on all parameters. If you wish to
specify your own priors you can supply them as brmsprior objects or as a
list such as `priors = list("group1" = c(0,3), "group2" = c(0,1))` where
the order of values is Mu, Sigma. Any non-brms backend will instead use
[`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html) to
fit the model unless the "flexsurv" type is specified. Distributions
will be passed to `survreg` where options are "weibull", "exponential",
"gaussian", "logistic","lognormal" and "loglogistic" if type = "survreg"
or to
[`flexsurv::flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
if type = "flexsurv" where options are "gengamma", "gengamma.orig",
"genf", "genf.orig", "weibull", "gamma", "exp", "llogis", "lnorm",
"gompertz", "exponential", and "lognormal". In `flexsurvreg`
distributional modeling is supported and additional formula can be
passed as a list to the sigma argument of growthSS in the same way as to
the anc argument of
[`flexsurv::flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).
Further additional arguments should be supplied via `fitGrowth` if
desired.

## See also

[fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)
for fitting the model specified by this list and
[mvSS](https://danforthcenter.github.io/pcvr/reference/mvSS.md) for the
multi-value trait equivalent.

## Examples

``` r

simdf <- growthSim("logistic",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
)
ss <- growthSS(
  model = "logistic", form = y ~ time | id / group,
  sigma = "spline", df = simdf,
  start = list("A" = 130, "B" = 12, "C" = 3), type = "brms"
)
lapply(ss, class)
#> $formula
#> [1] "brmsformula" "bform"      
#> 
#> $prior
#> [1] "brmsprior"  "data.frame"
#> 
#> $initfun
#> [1] "function"
#> 
#> $df
#> [1] "data.frame"
#> 
#> $family
#> [1] "character"
#> 
#> $pcvrForm
#> [1] "formula"
#> 
#> $type
#> [1] "character"
#> 
#> $model
#> [1] "character"
#> 
#> $call
#> [1] "call"
#> 
ss$initfun()
#> $b_A
#> [1] 0.3400957 0.3126611
#> 
#> $b_B
#> [1] 0.8370048 2.9719342
#> 
#> $b_C
#> [1] 0.2303612 1.3446416
#> 
# the next step would typically be compiling/fitting the model
# here we use very few chains and very few iterations for speed, but more of both is better.
# \donttest{
fit_test <- fitGrowth(ss,
  iter = 500, cores = 1, chains = 1, backend = "cmdstanr",
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)
#> Start sampling
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#> Intercept_sigma, bs_sigma, zs_sigma_1_1, sds_sigma_1, zs_sigma_2_1, sds_sigma_2, Intercept_nu
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> Running MCMC with 1 chain...
#> 
#> Chain 1 Iteration:   1 / 500 [  0%]  (Warmup) 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff4fc9f62c.stan', line 131, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Scale parameter[1] is inf, but must be positive finite! (in '/tmp/RtmpfB9aQg/model-1eff4fc9f62c.stan', line 131, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: student_t_lpdf: Location parameter[1] is -nan, but must be finite! (in '/tmp/RtmpfB9aQg/model-1eff4fc9f62c.stan', line 131, column 4 to column 48)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Iteration: 100 / 500 [ 20%]  (Warmup) 
#> Chain 1 Iteration: 200 / 500 [ 40%]  (Warmup) 
#> Chain 1 Iteration: 251 / 500 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 350 / 500 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 450 / 500 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 500 / 500 [100%]  (Sampling) 
#> Chain 1 finished in 88.9 seconds.
# }


# formulas and priors will look different if there is only one group in the data

ex <- growthSim("linear", n = 20, t = 25, params = list("A" = 2))
ex_ss <- growthSS(
  model = "linear", form = y ~ time | id / group, sigma = "spline",
  df = ex, start = list("A" = 1), type = "brms"
)

ex_ss$prior # no coef level grouping for priors
#>                    prior     class coef group resp  dpar nlpar lb   ub tag
#>     student_t(3, 0, 2.5) Intercept                 sigma                  
#>         normal(2.7, 0.8) Intercept                    nu                  
#>  lognormal(log(1), 0.25)         b                           A  0 <NA>    
#>   source
#>  default
#>  default
#>     user
ex_ss$formula # intercept only model for A
#> y ~ A * time 
#> autocor ~ tructure(list(), class = "formula", .Environment = <environment>)
#> sigma ~ s(time)
#> nu ~ 1
#> A ~ 1

ex2 <- growthSim("linear", n = 20, t = 25, params = list("A" = c(2, 2.5)))
ex2_ss <- growthSS(
  model = "linear", form = y ~ time | id / group, sigma = "spline",
  df = ex2, start = list("A" = 1), type = "brms"
)
ex2_ss$prior # has coef level grouping for priors
#>                    prior     class coef group resp  dpar nlpar lb   ub tag
#>     student_t(3, 0, 2.5) Intercept                 sigma                  
#>         normal(2.7, 0.8) Intercept                    nu                  
#>  lognormal(log(1), 0.25)         b                           A  0 <NA>    
#>   source
#>  default
#>  default
#>     user
ex2_ss$formula # specifies an A intercept for each group and splines by group for sigma
#> y ~ A * time 
#> autocor ~ tructure(list(), class = "formula", .Environment = <environment>)
#> sigma ~ s(time, by = group)
#> nu ~ 1
#> A ~ 0 + group
```
