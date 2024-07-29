library(testthat)
library(pcvr)

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

test_that(".decayChngptForm flips a formula returned by another chngpt function", {
  starter <- .linearChngptForm("x_var",
    position = 1, dpar = NULL,
    priors = list("linear1A" = 1, "changepoint" = 5, "linear2A" = 10)
  )
  final <- .decayChngptForm(starter)
  expect_equal(final$form, "-linear1A * x_var")
})

conditions <- data.frame(
  x = "x_var", position = rep(c(1:2), each = 4),
  dpar = rep(c("", "sigma"), times = 4),
  chngpt = rep(c("changepoint", "fixedChangePoint", "fixedChangePoint", "changepoint"),
    times = 2
  )
)

test_that(".intChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list("int1" = 1, "c" = 5, "int2" = 10)
    names(prior)[2] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .intChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that(".gamChngptForm assembles a formula in all conditions", {
  gam_conditions <- conditions[conditions$position == 2, ]
  for (i in seq_len(nrow(gam_conditions))) {
    prior <- list("other1" = 1, "c" = 5, "int2" = 10)
    names(prior)[2] <- paste0(
      gam_conditions[i, "dpar"], gam_conditions[i, "chngpt"],
      gam_conditions[i, "position"]
    )
    iter <- .gamChngptForm(gam_conditions[i, "x"],
      position = gam_conditions[i, "position"],
      dpar = gam_conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params", "splineVar")) # gam terms
  }
  err_condition <- conditions[conditions$position == 1, ][1, ]
  prior <- list("other1" = 1, "c" = 5, "int2" = 10)
  names(prior)[2] <- paste0(
    err_condition[i, "dpar"], err_condition[i, "chngpt"],
    err_condition[i, "position"]
  )
  expect_error(
    .gamChngptForm(err_condition[i, "x"],
      position = err_condition[i, "position"],
      dpar = err_condition[i, "dpar"],
      priors = prior
    )
  )
})

test_that(".linearChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list("linear1A" = 1, "c" = 5, "linear2A" = 10)
    names(prior)[2] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .linearChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that("logarithmicChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list("logarithmic1A" = 1, "c" = 5, "logarithmic2A" = 10)
    names(prior)[2] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .logarithmicChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that(".exponentialChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "exponential1A" = 1, "exponential1B" = 0.1,
      "c" = 5,
      "exponential2A" = 2, "exponential2B" = 0.2
    )
    names(prior)[3] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .exponentialChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that(".monomolecularChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "monomolecular1A" = 1, "monomolecular1B" = 0.1,
      "changepoint" = 5,
      "monomolecular2A" = 2, "monomolecular2B" = 0.2
    )
    names(prior)[3] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .monomolecularChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that(".powerlawChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "powerlaw1A" = 1, "powerlaw1B" = 0.1,
      "changepoint" = 5,
      "powerlaw2A" = 2, "powerlaw2B" = 0.2
    )
    names(prior)[3] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .powerlawChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that(".frechetChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "frechet1A" = 1, "frechet1B" = 1, "frechet1C" = 1,
      "changepoint" = 5,
      "frechet2A" = 2, "frechet2B" = 2, "frechet2C" = 2
    )
    names(prior)[4] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .frechetChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that(".gompertzChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "gompertz1A" = 1, "gompertz1B" = 1, "gompertz1C" = 1,
      "changepoint" = 5,
      "gompertz2A" = 2, "gompertz2B" = 2, "gompertz2C" = 2
    )
    names(prior)[4] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .gompertzChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that(".gumbelChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "gumbel1A" = 1, "gumbel1B" = 1, "gumbel1C" = 1,
      "changepoint" = 5,
      "gumbel2A" = 2, "gumbel2B" = 2, "gumbel2C" = 2
    )
    names(prior)[4] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .gumbelChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that(".logisticChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "logistic1A" = 1, "logistic1B" = 1, "logistic1C" = 1,
      "changepoint" = 5,
      "logistic2A" = 2, "logistic2B" = 2, "logistic2C" = 2
    )
    names(prior)[4] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .logisticChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})

test_that(".weibullChngptForm assembles a formula in all conditions", {
  for (i in seq_len(nrow(conditions))) {
    prior <- list(
      "weibull1A" = 1, "weibull1B" = 1, "weibull1C" = 1,
      "changepoint" = 5,
      "weibull2A" = 2, "weibull2B" = 2, "weibull2C" = 2
    )
    names(prior)[4] <- paste0(conditions[i, "dpar"], conditions[i, "chngpt"], conditions[i, "position"])
    iter <- .weibullChngptForm(conditions[i, "x"],
      position = conditions[i, "position"],
      dpar = conditions[i, "dpar"],
      priors = prior
    )
    expect_equal(names(iter), c("form", "cp", "cpInt", "params"))
  }
})
