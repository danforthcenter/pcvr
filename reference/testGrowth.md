# Hypothesis testing for [fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md) models.

Hypothesis testing for
[fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)
models.

## Usage

``` r
testGrowth(ss = NULL, fit, test = "A")
```

## Arguments

- ss:

  A list output from
  [growthSS](https://danforthcenter.github.io/pcvr/reference/growthSS.md).
  This is not required for nls, nlme, and brms models if `test` is given
  in
  [`brms::hypothesis`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html)
  style as a written statement.

- fit:

  A model (or list of nlrq models) output from
  [fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md).
  For brms models this can also be a data.frame of draws.

- test:

  A description of the hypothesis to test. This can take two main forms,
  either the parameter names to vary before comparing a nested model
  ("A", "B", "C") using an anova or a hypothesis test/list of hypothesis
  tests written as character strings. The latter method is not
  implemented for `nlrq` models. If this is a vector of parameters to
  test in the model then they should be parameters which vary by group
  in your original model and that you want to test against a null model
  where they do not vary by group. Alternatively for nlrq models this
  can be a comparison of model terms written as
  `"group_X|tau|par - group_Y|tau|par"`, which uses a fat tailed T
  distribution to make comparisons on the means of each quantile
  estimate. For GAMs these tests compare the model with splines either
  by group or interacting with group to a model that ignores the
  grouping in the data. If this is a list of hypothesis tests then they
  should describe tests similar to "A.group1 - A.group2\*1.1" and can be
  thought of as contrasts. For brms models the "test" argument is passed
  to brms::hypothesis, which has extensive documentation and is very
  flexible. Note that for survreg the
  [`survival::survdiff`](https://rdrr.io/pkg/survival/man/survdiff.html)
  function is used so fewer hypothesis testing options are available and
  flexsurv models are tested using contrasts via
  [`flexsurv::standsurv`](http://chjackson.github.io/flexsurv-dev/reference/standsurv.md).

## Value

A list containing an anova object comparing non-linear growth models and
the null model.

## Details

For nls and nlme models an anova is run and returned as part of a list
along with the null model. For nlrq models several assumptions are made
and a likelihood ratio test for each tau is run and returned as a list.

## See also

[growthSS](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
and
[fitGrowth](https://danforthcenter.github.io/pcvr/reference/fitGrowth.md)
for making compatible models,
[growthPlot](https://danforthcenter.github.io/pcvr/reference/growthPlot.md)
for hypothesis testing on compatible models.

## Examples

``` r
set.seed(123)
simdf <- growthSim("logistic",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
)
ss <- suppressMessages(growthSS(
  model = "logistic", form = y ~ time | id / group,
  df = simdf, type = "nlrq"
))
fit <- fitGrowth(ss)
testGrowth(ss, fit, "A")
#> $`0.5`
#> Model 1: y ~ A[group]/(1 + exp((B[group] - time)/C[group]))
#> Model 2: y ~ A/(1 + exp((B[group] - time)/C[group]))
#>   #Df  LogLik Df  Chisq Pr(>Chisq)    
#> 1   6 -3899.4                         
#> 2   5 -4028.2 -1 257.71  < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
testGrowth(ss, fit, "a|0.5|A > b|0.5|A")
#>      Value Std. Error   t value Pr(>|t|) par group prob.greater
#> 1 202.1163   2.817899  71.72588        0   A     a    0.9998036
#> 4 166.5326   1.252554 132.95448        0   A     b           NA

ss2 <- suppressMessages(growthSS(
  model = "logistic", form = y ~ time | id / group,
  df = simdf, type = "nls"
))
fit2 <- fitGrowth(ss2)
testGrowth(ss2, fit2, "A")$anova
#> Analysis of Variance Table
#> 
#> Model 1: y ~ A/(1 + exp((B[group] - time)/C[group]))
#> Model 2: y ~ A[group]/(1 + exp((B[group] - time)/C[group]))
#>   Res.Df Res.Sum Sq Df Sum Sq F value    Pr(>F)    
#> 1    995     204775                                
#> 2    994     171686  1  33089  191.57 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
coef(fit2) # check options for contrast testing
#>        A1        A2        B1        B2        C1        C2 
#> 201.68366 165.22719  12.94484  10.80995   3.10079   3.44904 
testGrowth(ss2, fit2, "A1 - A2*1.1")
#>          Form Estimate       SE  t-value      p-value
#> 1 A1 - A2*1.1 19.93375 2.523292 7.899899 7.368006e-15
```
