# Bayesian testing using conjugate priors and method of moments for single or multi value traits.

Function to perform bayesian tests and ROPE comparisons using single or
multi value traits with several distributions.

## Usage

``` r
conjugate(
  s1 = NULL,
  s2 = NULL,
  method = c("t", "gaussian", "beta", "binomial", "lognormal", "lognormal2", "poisson",
    "negbin", "vonmises", "vonmises2", "uniform", "pareto", "gamma", "bernoulli",
    "exponential", "bivariate_uniform", "bivariate_gaussian", "bivariate_lognormal"),
  priors = NULL,
  plot = NULL,
  rope_range = NULL,
  rope_ci = 0.89,
  cred.int.level = 0.89,
  hypothesis = "equal",
  bayes_factor = NULL,
  support = NULL
)
```

## Arguments

- s1:

  A data.frame or matrix of multi value traits or a vector of single
  value traits. If a multi value trait is used then column names should
  include a number representing the "bin". Alternatively for
  distributions other than "binomial" (which requires list data with
  "successes" and "trials" as numeric vectors in the list, see examples)
  this can be a formula specifying `outcome ~ group` where group has
  exactly 2 levels. If using wide MV trait data then the formula should
  specify column positions ~ grouping such as `1:180 ~ group`. This
  sample is shown in red if plotted.

- s2:

  An optional second sample, or if s1 is a formula then this should be a
  dataframe. This sample is shown in blue if plotted.

- method:

  The distribution/method to use. Currently "t", "gaussian", "beta",
  "binomial", "lognormal", "lognormal2", "poisson", "negbin" (negative
  binomial), "uniform", "pareto", "gamma", "bernoulli", "exponential",
  "vonmises", and "vonmises2" are supported. The count (binomial,
  poisson and negative binomial), bernoulli, exponential, and pareto
  distributions are only implemented for single value traits due to
  their updating and/or the nature of the input data. The "t" and
  "gaussian" methods both use a T distribution with "t" testing for a
  difference of means and "gaussian" testing for a difference in the
  distributions (similar to a Z test). Both Von Mises options are for
  use with circular data (for instance hue values when the circular
  quality of the data is relevant). Note that non-circular distributions
  can be compared to each other. This should only be done with caution
  and may not be supported in all downstream functions. There are also 3
  bivariate conjugate priors that are supported for use with single
  value data. Those are "bivariate_uniform", "bivariate_gaussian" and
  "bivariate_lognormal".

- priors:

  Prior distributions described as a list of lists. If this is a single
  list then it will be duplicated for the second sample, which is
  generally a good idea if both samples use the same distribution
  (method). Elements in the inner lists should be named for the
  parameter they represent (see examples). These names vary by method
  (see details). By default this is NULL and weak priors (generally
  jeffrey's priors) are used. The `posterior` part of output can also be
  recycled as a new prior if Bayesian updating is appropriate for your
  use.

- plot:

  deprecated, use `plot` method instead.

- rope_range:

  Optional vector specifying a region of practical equivalence. This
  interval is considered practically equivalent to no effect.
  Kruschke (2018) suggests c(-0.1, 0.1) as a broadly reasonable ROPE for
  standardized parameters. That range could also be rescaled by a
  standard deviation/magnitude for non-standardized parameters, but
  ultimately this should be informed by your setting and scientific
  question. See Kruschke (2018) for details on ROPE and other Bayesian
  methods to aide decision-making
  [doi:10.1177/2515245918771304](https://doi.org/10.1177/2515245918771304)
  and [doi:10.1037/a0029146](https://doi.org/10.1037/a0029146) .

- rope_ci:

  The credible interval probability to use for ROPE. Defaults to 0.89.

- cred.int.level:

  The credible interval probability to use in computing HDI for samples,
  defaults to 0.89.

- hypothesis:

  Direction of a hypothesis if two samples are provided. Options are
  "unequal", "equal", "greater", and "lesser", read as "sample1 greater
  than sample2".

- bayes_factor:

  Optional point or interval to evaluate bayes factors on. Note that
  this generally only makes sense to use if you have informative priors
  where the change in odds between prior and posterior is meaningful
  about the data. If this is non-NULL then columns of bayes factors are
  added to the summary output. Note these are only implemented for
  univariate distributions.

- support:

  Deprecated

## Value

A
[conjugate-class](https://danforthcenter.github.io/pcvr/reference/conjugate-class.md)
object with slots including:

- **summary**: A data frame containing HDI/HDE values for each sample
  and the ROPE as well as posterior probability of the hypothesis and
  ROPE test (if specified). The HDE is the "Highest Density Estimate" of
  the posterior, that is the tallest part of the probability density
  function. The HDI is the Highest Density Interval, which is an
  interval that contains X% of the posterior distribution, so
  `cred.int.level = 0.8` corresponds to an HDI that includes 80 percent
  of the posterior probability. Bayes factors are calculated as
  posterior/prior for each sample.

- **posterior**: A list of updated parameters in the same format as the
  prior for the given method. If desired this does allow for Bayesian
  updating.

- **prior**: The prior in a list with the same format as the posterior.

- **plot**: A ggplot showing the distribution of samples and optionally
  the distribution of differences/ROPE.

- **plot_parameters**: Parameters used in making a plot of the data.
  Contains support range and posterior recoded to use a density
  function.

- **data**: Data from s1 and s2 arguments.

- **call**: The function call.

## Details

Prior distributions default to be weakly informative and in some cases
you may wish to change them.

- **"t", "gaussian", and "lognormal":**
  `priors = list(mu = 0, sd = 10)`, where mu is the mean and sd is
  standard deviation. For the lognormal method these describe the normal
  distribution of the mean parameter for lognormal data and are on the
  log scale.

- **"beta", "bernoulli", and "binomial":**
  `priors = list(a = 0.5, b = 0.5)`, where a and b are shape parameters
  of the beta distribution. Note that for the binomial distribution this
  is used as the prior for success probability P, which is assumed to be
  beta distributed as in a beta-binomial distribution.

- **"lognormal2":** `priors = list(a = 1, b = 1) `, where a and b are
  the shape and scale parameters of the gamma distribution of lognormal
  data's precision parameter (using the alternative mu, precision
  parameterization).

- **"gamma":**
  `priors = list(shape = 0.5, scale = 0.5, known_shape = 1)`, where
  shape and scale are the respective parameters of the gamma distributed
  rate (inverse of scale) parameter of gamma distributed data.

- **"poisson" and "exponential":** `priors = list(a = 0.5,b = 0.5)`,
  where a and b are shape and rate parameters of the gamma distribution.

- **"negbin":** `priors = list(r = 10, a = 0.5, b = 0.5)`, where r is
  the r parameter of the negative binomial distribution (representing
  the number of successes required) and where a and b are shape
  parameters of the beta distribution. Note that the r value is not
  updated. The conjugate beta prior is only valid when r is fixed and
  known, which is a limitation for this method.

- **"uniform":** `list(scale = 0.5, location = 0.5)`, where scale is the
  scale parameter of the pareto distributed upper boundary and location
  is the location parameter of the pareto distributed upper boundary.
  Note that different sources will use different terminology for these
  parameters. These names were chosen for consistency with the
  `extraDistr` implementation of the pareto distribution. On Wikipedia
  the parameters are called shape and scale, corresponding to
  extraDistr's scale and location respectively, which can be confusing.
  Note that the lower boundary of the uniform is assumed to be 0.

- **"pareto":** `list(a = 1, b = 1, known_location = min(data))`, where
  a and b are the shape and scale parameters of the gamma distribution
  of the pareto distribution's scale parameter. In this case location is
  assumed to be constant and known, which is less of a limitation than
  knowing r for the negative binomial method since location will
  generally be right around/just under the minimum of the sample data.
  Note that the pareto method is only implemented currently for single
  value traits since one of the statistics needed to update the gamma
  distribution here is the product of the data and we do not currently
  have a method to calculate a similar sufficient statistic from multi
  value traits.

- **"vonmises":**
  `list(mu = 0, kappa = 0.5, boundary = c(-pi, pi), known_kappa = 1, n = 1)`,
  where mu is the direction of the circular distribution (the mean),
  kappa is the precision of the mean, boundary is a vector including the
  two values that are the where the circular data "wraps" around the
  circle, known_kappa is the fixed value of precision for the total
  distribution, and n is the number of prior observations. This Von
  Mises option updates the conjugate prior for the mean direction, which
  is itself Von-Mises distributed. This in some ways is analogous to the
  T method, but assuming a fixed variance when the mean is updated. Note
  that due to how the rescaling works larger circular boundaries can be
  slow to plot.

- **"vonmises2":**
  `priors = list(mu = 0, kappa = 0.5, boundary = c(-pi, pi), n = 1)`,
  where mu and kappa are mean direction and precision of the von mises
  distribution, boundary is a vector including the two values that are
  the where the circular data "wraps" around the circle, and n is the
  number of prior observations. This Von-Mises implementation does not
  assume constant variance and instead uses MLE to estimate kappa from
  the data and updates the kappa prior as a weighted average of the data
  and the prior. The mu parameter is then updated per Von-Mises
  conjugacy.

- **"bivariate_uniform":**
  `list(location_l = 1, location_u = 2, scale = 1)`, where scale is the
  shared scale parameter of the pareto distributed upper and lower
  boundaries and location l and u are the location parameters for the
  Lower (l) and Upper (u) boundaries of the uniform distribution. Note
  this uses the same terminology for the pareto distribution's
  parameters as the "uniform" method.

- **"bivariate_gaussian" and "bivariate_lognormal":**
  `list(mu = 0, sd = 10, a = 1, b = 1)`, where mu and sd are the mean
  and standard deviation of the Normal distribution of the data's mean
  and a and b are the shape and scale of the gamma distribution on
  precision. Note that internally this uses the Mu and Precision
  parameterization of the normal distribution and those are the
  parameters shown in the plot and tested, but priors use Mu and SD for
  the normal distribution of the mean.

## See also

[barg](https://danforthcenter.github.io/pcvr/reference/barg.md) for
additional reporting

## Examples

``` r
mv_ln <- mvSim(
  dists = list(
    rlnorm = list(meanlog = log(130), sdlog = log(1.2)),
    rlnorm = list(meanlog = log(100), sdlog = log(1.3))
  ),
  n_samples = 30
)

# lognormal mv
ln_mv_ex <- conjugate(
  s1 = mv_ln[1:30, -1], s2 = mv_ln[31:60, -1], method = "lognormal",
  priors = list(mu = 5, sd = 2),
  rope_range = c(-40, 40), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal", support = NULL
)
#> Warning: support argument is deprecated.

# lognormal sv
ln_sv_ex <- conjugate(
  s1 = rlnorm(100, log(130), log(1.3)), s2 = rlnorm(100, log(100), log(1.6)),
  method = "lognormal",
  priors = list(mu = 5, sd = 2),
  rope_range = NULL, rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal", support = NULL
)
#> Warning: support argument is deprecated.

# Z test mv example

mv_gauss <- mvSim(
  dists = list(
    rnorm = list(mean = 50, sd = 10),
    rnorm = list(mean = 60, sd = 12)
  ),
  n_samples = 30
)

gauss_mv_ex <- conjugate(
  s1 = mv_gauss[1:30, -1], s2 = mv_gauss[31:60, -1], method = "gaussian",
  priors = list(mu = 30, sd = 10),
  rope_range = c(-25, 25), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal", support = NULL
)
#> Warning: support argument is deprecated.

# T test sv example with two different priors

gaussianMeans_sv_ex <- conjugate(
  s1 = rnorm(10, 50, 10), s2 = rnorm(10, 60, 12), method = "t",
  priors = list(list(mu = 40, sd = 10), list(mu = 45, sd = 8)),
  rope_range = c(-5, 8), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal", support = NULL
)
#> Warning: support argument is deprecated.

# beta mv example

set.seed(123)
mv_beta <- mvSim(
  dists = list(
    rbeta = list(shape1 = 5, shape2 = 8),
    rbeta = list(shape1 = 10, shape2 = 10)
  ),
  n_samples = c(30, 20)
)

beta_mv_ex <- conjugate(
  s1 = mv_beta[1:30, -1], s2 = mv_beta[31:50, -1], method = "beta",
  priors = list(a = 0.5, b = 0.5),
  rope_range = c(-0.1, 0.1), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal",
  bayes_factor = 0.5 # note this may not be reasonable with these priors
)

# beta sv example

beta_sv_ex <- conjugate(
  s1 = rbeta(20, 5, 5), s2 = rbeta(20, 8, 5), method = "beta",
  priors = list(a = 0.5, b = 0.5),
  rope_range = c(-0.1, 0.1), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal",
  bayes_factor = c(0.5, 0.75) # note this may not be reasonable with these priors
)

# binomial sv example
# note that specifying trials = 20 would also work
# and the number of trials will be recycled to the length of successes

binomial_sv_ex <- conjugate(
  s1 = list(successes = c(15, 14, 16, 11), trials = c(20, 20, 20, 20)),
  s2 = list(successes = c(7, 8, 10, 5), trials = c(20, 20, 20, 20)), method = "binomial",
  priors = list(a = 0.5, b = 0.5),
  rope_range = c(-0.1, 0.1), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)

# poisson sv example

poisson_sv_ex <- conjugate(
  s1 = rpois(20, 10), s2 = rpois(20, 8), method = "poisson",
  priors = list(a = 0.5, b = 0.5),
  rope_range = c(-1, 1), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)

# negative binomial sv example
# knowing r (required number of successes) is an important caveat for this method.
# in the current implementation we suggest using the poisson method for data such as leaf counts

negbin_sv_ex <- conjugate(
  s1 = rnbinom(20, 10, 0.5), s2 = rnbinom(20, 10, 0.25), method = "negbin",
  priors = list(r = 10, a = 0.5, b = 0.5),
  rope_range = c(-1, 1), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)

# von mises mv example

mv_gauss <- mvSim(
  dists = list(
    rnorm = list(mean = 50, sd = 10),
    rnorm = list(mean = 60, sd = 12)
  ),
  n_samples = c(30, 40)
)
vm1_ex <- conjugate(
  s1 = mv_gauss[1:30, -1],
  s2 = mv_gauss[31:70, -1],
  method = "vonmises",
  priors = list(mu = 45, kappa = 1, boundary = c(0, 180), known_kappa = 1, n = 1),
  rope_range = c(-1, 1), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)

# von mises 2 sv example
vm2_ex <- conjugate(
  s1 = brms::rvon_mises(10, 2, 2),
  s2 = brms::rvon_mises(15, 3, 3),
  method = "vonmises2",
  priors = list(mu = 0, kappa = 0.5, boundary = c(-pi, pi), n = 1),
  cred.int.level = 0.95
)
```
