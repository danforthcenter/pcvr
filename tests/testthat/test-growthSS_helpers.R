library(testthat)
library(pcvr)

test_that("GrowthSS Helpers for double logistic Data work", {
  set.seed(123)
  df <- growthSim(
    "logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  expect_warning(
    ss <- suppressMessages(growthSS(
      model = "double logistic", form = y ~ time | id / group,
      df = df, type = "nls"
    ))
  )
  expect_type(ss, "list")
  expect_warning(
    ss <- suppressMessages(growthSS(
      model = "double logistic", form = y ~ time | id / group,
      df = df, type = "nlrq", tau = c(0.45, 0.55)
    ))
  )
  expect_type(ss, "list")
  expect_warning(
    ss <- suppressMessages(growthSS(
      model = "double logistic", form = y ~ time | id / group,
      df = df, type = "nlme"
    ))
  )
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "double logistic", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3, "A2" = 100, "B2" = 10, "C2" = 3),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for double gompertz Data work", {
  set.seed(123)
  df <- growthSim(
    "logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  expect_warning(
    ss <- suppressMessages(growthSS(
      model = "double gompertz", form = y ~ time | id / group,
      df = df, type = "nls"
    ))
  )
  expect_type(ss, "list")
  expect_warning(
    ss <- suppressMessages(growthSS(
      model = "double gompertz", form = y ~ time | id / group,
      df = df, type = "nlrq", tau = c(0.45, 0.55)
    ))
  )
  expect_type(ss, "list")
  expect_warning(
    ss <- suppressMessages(growthSS(
      model = "double gompertz", form = y ~ time | id / group,
      df = df, type = "nlme"
    ))
  )
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "double gompertz", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3, "A2" = 100, "B2" = 10, "C2" = 3),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for Logistic Data work", {
  set.seed(123)
  df <- growthSim("logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | id / group,
    df = df, type = "nlrq", tau = c(0.45, 0.55)
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for 4 Parameter Logistic Data work", {
  set.seed(123)
  df <- growthSim("logistic4",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5), "D" = c(5, 20))
  )
  ss <- suppressMessages(growthSS(
    model = "logistic4", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logistic4", form = y ~ time | id / group,
    df = df, type = "nlrq", tau = c(0.45, 0.55)
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logistic4", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logistic4", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3, "D" = 5),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for 5 Parameter Logistic Data work", {
  set.seed(123)
  df <- growthSim("logistic5",
    n = 20, t = 25,
    params = list(
      "A" = c(200, 160), "B" = c(13, 11),
      "C" = c(3, 3.5), "D" = c(5, 20),
      "E" = c(1, 1.5)
    )
  )
  ss <- suppressMessages(growthSS(
    model = "logistic5", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logistic5", form = y ~ time | id / group,
    df = df, type = "nlrq", tau = c(0.45, 0.55)
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logistic5", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logistic5", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3, "D" = 5, "E" = 1),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for Gompertz Data work", {
  set.seed(123)
  df <- growthSim("gompertz",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(0.25, 0.5))
  )
  ss <- suppressMessages(growthSS(
    model = "gompertz", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gompertz", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gompertz", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gompertz", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for Weibull Data work", {
  set.seed(123)
  df <- growthSim("weibull",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  ss <- suppressMessages(growthSS(
    model = "weibull", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "weibull", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "weibull", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "weibull", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for Gumbel Data work", {
  set.seed(123)
  df <- growthSim("gumbel",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  ss <- suppressMessages(growthSS(
    model = "gumbel", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gumbel", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gumbel", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gumbel", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for Frechet Data work", {
  set.seed(123)
  df <- growthSim("frechet",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  ss <- suppressMessages(growthSS(
    model = "frechet", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "frechet", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "frechet", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "frechet", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for Monomolecular Data work", {
  set.seed(123)
  df <- growthSim("monomolecular",
    n = 20, t = 25,
    params = list("A" = c(100, 100), "B" = c(0.5, 0.5))
  )
  ss <- suppressMessages(growthSS(
    model = "monomolecular", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "monomolecular", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "monomolecular", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "monomolecular", form = y ~ time | id / group,
    start = list("A" = 15, "B" = 1),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for Power Law Data work", {
  set.seed(123)
  df <- growthSim("power law",
    n = 20, t = 25,
    params = list("A" = c(10, 10), "B" = c(0.75, 0.75))
  )
  ss <- suppressMessages(growthSS(
    model = "power law", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "power law", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "power law", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "power law", form = y ~ time | id / group,
    start = list("A" = 15, "B" = 1),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for Exponential Data work", {
  set.seed(123)
  df <- growthSim("exponential",
    n = 20, t = 25,
    params = list("A" = c(15, 15), "B" = c(0.01, 0.01))
  )
  ss <- suppressMessages(growthSS(
    model = "exponential", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "exponential", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "exponential", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "exponential", form = y ~ time | id / group,
    start = list("A" = 15, "B" = 1),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for Logarithmic Data work", {
  set.seed(123)
  df <- growthSim("logarithmic",
    n = 20, t = 25,
    params = list("A" = c(15, 15))
  )
  ss <- suppressMessages(growthSS(
    model = "logarithmic", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logarithmic", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logarithmic", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "logarithmic", form = y ~ time | id / group,
    start = list("A" = 15),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for linear Data work", {
  set.seed(123)
  df <- growthSim("linear",
    n = 20, t = 25,
    params = list("A" = c(15, 15))
  )
  ss <- suppressMessages(growthSS(
    model = "linear", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "linear", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "linear", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "linear", form = y ~ time | id / group,
    start = list("A" = 15),
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for gams work", {
  set.seed(123)
  df <- growthSim("linear",
    n = 20, t = 25,
    params = list("A" = c(15, 15))
  )
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group,
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group,
    df = df, type = "mgcv"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time,
    df = df, type = "mgcv"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for bragg DRMs run", {
  set.seed(123)
  df <- growthSim("bragg",
    n = 20, t = 100,
    list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
  )
  ss <- suppressMessages(growthSS(
    model = "bragg", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "bragg", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "bragg", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "bragg", form = y ~ time | id / group,
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for lorentz DRMs run", {
  set.seed(123)
  df <- growthSim("lorentz",
    n = 20, t = 100,
    list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
  )
  ss <- suppressMessages(growthSS(
    model = "lorentz", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "lorentz", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "lorentz", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "lorentz", form = y ~ time | id / group,
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that("GrowthSS Helpers for beta DRMs run", {
  set.seed(123)
  df <- growthSim("beta",
    n = 20, t = 100,
    list("A" = 10, "B" = 1.2, "C" = 15, "D" = 8, "E" = 19)
  )
  ss <- suppressMessages(growthSS(
    model = "beta", form = y ~ time | id / group,
    df = df, type = "nls"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "beta", form = y ~ time | id / group,
    df = df, type = "nlrq"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "beta", form = y ~ time | id / group,
    df = df, type = "nlme"
  ))
  expect_type(ss, "list")
  ss <- suppressMessages(growthSS(
    model = "beta", form = y ~ time | id / group,
    df = df, type = "brms"
  ))
  expect_type(ss, "list")
})

test_that(".brmSS messages about complex models", {
  set.seed(123)
  df <- growthSim(
    "gompertz",
    n = 20, t = 25,
    params = list(
      "A" = rnorm(20, 200, 20),
      "B" = rnorm(20, 12, 2),
      "C" = rnorm(20, 0.4, 0.1)
    )
  )
  expect_message(growthSS(
    model = "gompertz", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3),
    df = df, type = "brms"
  ))
})

test_that(".brmSS assembles decay model", {
  set.seed(123)
  df <- growthSim(
    "gompertz",
    n = 20, t = 25,
    params = list(
      "A" = rnorm(2, 200, 20),
      "B" = rnorm(2, 12, 2),
      "C" = rnorm(2, 0.4, 0.1)
    )
  )
  ss <- growthSS(
    model = "decay gompertz", form = y ~ time | id / group,
    start = list("A" = 100, "B" = 10, "C" = 3),
    df = df, type = "brms"
  )
  expect_true(grepl("^-", as.character(ss$formula$formula)[3]))
})

test_that(".brmSS warns about ambiguous hierarchy", {
  set.seed(123)
  df <- growthSim(
    "gompertz",
    n = 20, t = 25,
    params = list(
      "A" = rnorm(2, 200, 20),
      "B" = rnorm(2, 12, 2),
      "C" = rnorm(2, 0.4, 0.1)
    )
  )
  df$covar <- rnorm(nrow(df))
  expect_warning(growthSS(
    model = "gompertz", form = y ~ time + covar | id / group,
    start = list("AI" = 100, "BI" = 10, "CI" = 3),
    df = df, type = "brms"
  ))
})

test_that(".brmSS handles a truncated response", {
  set.seed(123)
  df <- growthSim(
    "gompertz",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(0.25, 0.5))
  )
  ss <- suppressMessages(growthSS(
    model = "gompertz", form = y[0, 100] ~ time | group,
    start = list("A" = 100, "B" = 10, "C" = 3),
    df = df, type = "brms"
  ))
  expect_true(grepl("trunc", as.character(ss$formula$formula[2])))
})

test_that(".decayChngptForm flips a formula returned by another chngpt function", {
  starter <- .linearChngptForm("x_var",
    position = 1, dpar = NULL,
    priors = list("linear1A" = 1, "changepoint" = 5, "linear2A" = 10)
  )
  final <- .decayChngptForm(starter)
  expect_equal(final$form, "-linear1A * x_var")
})

conditions <- data.frame(
  x = "x_var",
  dpar = rep(c("", "sigma"), times = 4),
  chngpt1 = rep(c("fixedChangePoint1", "changePoint1"), each = 4),
  chngpt2 = rep(c("fixedChangePoint2", "changePoint2"), each = 4)[c(7:8, 3:6, 1:2)],
  chngpt3 = rep(c("changePoint3", "fixedChangePoint3"), each = 4)
)

test_that(".intChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list("int1" = 1, "c" = 1, "int2" = 1, "c2" = 2, "int3" = 1, "c3" = 3, "int4" = 1)
    names(prior)[c(2, 4, 6)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .intChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that(".gamChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list("int1" = 1, "c" = 1, "int2" = 1, "c2" = 2, "int3" = 1, "c3" = 3, "int4" = 1)
    names(prior)[c(2, 4, 6)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    iter <- .gamChngptForm(conditions[i, "x"],
      position = 3,
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params", "splineVar"))
  }
  prior <- list("int1" = 1, "c" = 1, "int2" = 1, "c2" = 2, "int3" = 1, "c3" = 3, "int4" = 1)
  names(prior)[c(2, 4, 6)] <- paste0(conditions[1, "dpar"], c(conditions[1, 3:5]))
  expect_error(
    suppressMessages(
      .gamChngptForm(
        gam_conditions[1, "x"],
        position = 1,
        dpar = gam_conditions[1, "dpar"],
        priors = prior
      )
    )
  )
})

test_that(".linearChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "linear1" = 1, "c" = 1, "linear2" = 1, "c2" = 2,
      "linear3" = 1, "c3" = 3, "linear4" = 1
    )
    names(prior)[c(2, 4, 6)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .linearChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that("logarithmicChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "logarithmic1" = 1, "c" = 1, "logarithmic2" = 1, "c2" = 2,
      "logarithmic3" = 1, "c3" = 3, "logarithmic4" = 1
    )
    names(prior)[c(2, 4, 6)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .logarithmicChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that(".exponentialChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "exponential1A" = 1, "exponential1B" = 0.1,
      "c" = 1,
      "exponential2A" = 1, "exponential2B" = 0.1,
      "c2" = 2,
      "exponential3A" = 1, "exponential3B" = 0.1,
      "c3" = 3,
      "exponential4A" = 1, "exponential4B" = 0.1
    )
    names(prior)[c(3, 6, 9)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .exponentialChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that(".monomolecularChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "monomolecular1A" = 1, "monomolecular1B" = 0.1,
      "c" = 1,
      "monomolecular2A" = 1, "monomolecular2B" = 0.1,
      "c2" = 2,
      "monomolecular3A" = 1, "monomolecular3B" = 0.1,
      "c3" = 3,
      "monomolecular4A" = 1, "monomolecular4B" = 0.1
    )
    names(prior)[c(3, 6, 9)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .monomolecularChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that(".powerlawChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "powerlaw1A" = 1, "powerlaw1B" = 0.1,
      "c" = 1,
      "powerlaw2A" = 1, "powerlaw2B" = 0.1,
      "c2" = 2,
      "powerlaw3A" = 1, "powerlaw3B" = 0.1,
      "c3" = 3,
      "powerlaw4A" = 1, "powerlaw4B" = 0.1
    )
    names(prior)[c(3, 6, 9)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .powerlawChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that(".frechetChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "frechet1A" = 1, "frechet1B" = 1, "frechet1C" = 1,
      "c" = 1,
      "frechet2A" = 1, "frechet2B" = 1, "frechet2C" = 1,
      "c2" = 2,
      "frechet3A" = 1, "frechet3B" = 1, "frechet3C" = 1,
      "c3" = 3,
      "frechet4A" = 1, "frechet4B" = 1, "frechet4C" = 1
    )
    names(prior)[c(4, 8, 12)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .frechetChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that(".gompertzChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "gompertz1A" = 1, "gompertz1B" = 1, "gompertz1C" = 1,
      "c" = 1,
      "gompertz2A" = 1, "gompertz2B" = 1, "gompertz2C" = 1,
      "c2" = 2,
      "gompertz3A" = 1, "gompertz3B" = 1, "gompertz3C" = 1,
      "c3" = 3,
      "gompertz4A" = 1, "gompertz4B" = 1, "gompertz4C" = 1
    )
    names(prior)[c(4, 8, 12)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .gompertzChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that(".gumbelChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "gumbel1A" = 1, "gumbel1B" = 1, "gumbel1C" = 1,
      "c" = 1,
      "gumbel2A" = 1, "gumbel2B" = 1, "gumbel2C" = 1,
      "c2" = 2,
      "gumbel3A" = 1, "gumbel3B" = 1, "gumbel3C" = 1,
      "c3" = 3,
      "gumbel4A" = 1, "gumbel4B" = 1, "gumbel4C" = 1
    )
    names(prior)[c(4, 8, 12)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .gumbelChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that(".logisticChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "logistic1A" = 1, "logistic1B" = 1, "logistic1C" = 1,
      "c" = 1,
      "logistic2A" = 1, "logistic2B" = 1, "logistic2C" = 1,
      "c2" = 2,
      "logistic3A" = 1, "logistic3B" = 1, "logistic3C" = 1,
      "c3" = 3,
      "logistic4A" = 1, "logistic4B" = 1, "logistic4C" = 1
    )
    names(prior)[c(4, 8, 12)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .logisticChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that(".weibullChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "weibull1A" = 1, "weibull1B" = 1, "weibull1C" = 1,
      "c" = 1,
      "weibull2A" = 1, "weibull2B" = 1, "weibull2C" = 1,
      "c2" = 2,
      "weibull3A" = 1, "weibull3B" = 1, "weibull3C" = 1,
      "c3" = 3,
      "weibull4A" = 1, "weibull4B" = 1, "weibull4C" = 1
    )
    names(prior)[c(4, 8, 12)] <- paste0(conditions[i, "dpar"], c(conditions[i, 3:5]))
    for (ii in 1:3) {
      iter <- .weibullChngptForm(conditions[i, "x"],
        position = ii,
        dpar = conditions[i, "dpar"],
        priors = prior
      )
      expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
    }
  }
})

test_that("Error generation works", {
  prior <- list(
    "weibull1A" = 1, "weibull1B" = 1, "weibull1C" = 1,
    "weibull2A" = 1, "weibull2B" = 1, "weibull2C" = 1
  )
  expect_error(suppressMessages(.brmsChangePointHelper(
    model = "weibull + weibull", x = "x", y = "y",
    group = "group", dpar = FALSE,
    nTimes = 25, useGroup = TRUE, priors = prior, int = FALSE
  )))
})
