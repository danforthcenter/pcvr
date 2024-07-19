library(testthat)
library(pcvr)

all_pcvr_functions <- unclass(lsf.str(envir = asNamespace("pcvr"), all = T))
hidden_functions <- all_pcvr_functions[grepl("^[.]", all_pcvr_functions)]
length(hidden_functions) # 210 lol
relevant <- hidden_functions[grepl("_form_|formula|ChngptForm", hidden_functions)]
length(relevant) # 66
relevant

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
    df = df, type = "nlrq"
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
})

test_that("GrowthSS Helpers for bragg DRMs work", {
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

test_that("GrowthSS Helpers for lorentz DRMs work", {
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

test_that("GrowthSS Helpers for beta DRMs work", {
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

