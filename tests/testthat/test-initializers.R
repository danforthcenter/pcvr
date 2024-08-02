if (!interactive()) {
  pdf(NULL)
}
library(testthat)
library(pcvr)

test_that("double sigmoid models warn about starting values", {
  df <- data.frame()
  expect_warning(
    ss <- suppressMessages(
      .nlrqSS(
        model = "double logistic", form = y ~ time | id / group,
        tau = 0.5, df = df, start = NULL
      )
    )
  )
})

test_that("init logistic works", {
  df <- growthSim("logistic",
    n = 20, t = 25,
    params = list("A" = c(180), "B" = c(10), "C" = c(3))
  )
  expect_equal(names(.initlogistic(df, "time", "y", FALSE)), c("A", "B", "C"))
  expect_equal(names(.initlogistic(df, "time", "y", TRUE)), c("I", "A", "B", "C"))
  expect_error(.initlogistic(df[1:3, ], "time", "y", FALSE))
})

test_that("init gompertz works", {
  df <- growthSim("gompertz",
    n = 20, t = 25,
    params = list("A" = c(180), "B" = c(10), "C" = c(0.25))
  )
  expect_equal(names(.initgompertz(df, "time", "y", FALSE)), c("A", "B", "C"))
  expect_equal(names(.initgompertz(df, "time", "y", TRUE)), c("I", "A", "B", "C"))
  expect_error(.initgompertz(df[1:3, ], "time", "y", FALSE))
})

test_that("init frechet works", {
  df <- growthSim("frechet",
    n = 20, t = 25,
    params = list("A" = c(180), "B" = c(3), "C" = c(3))
  )
  expect_equal(names(.initfrechet(df, "time", "y", FALSE)), c("A", "B", "C"))
  expect_equal(names(.initfrechet(df, "time", "y", TRUE)), c("I", "A", "B", "C"))
  expect_error(.initfrechet(df[1:3, ], "time", "y", FALSE))
})

test_that("init gumbel works", {
  df <- growthSim("gumbel",
    n = 20, t = 25,
    params = list("A" = c(180), "B" = c(3), "C" = c(3))
  )
  expect_equal(names(.initgumbel(df, "time", "y", FALSE)), c("A", "B", "C"))
  expect_equal(names(.initgumbel(df, "time", "y", TRUE)), c("I", "A", "B", "C"))
  expect_error(.initgumbel(df[1:3, ], "time", "y", FALSE))
})

test_that("init weibull works", {
  df <- growthSim("weibull",
    n = 20, t = 25,
    params = list("A" = c(180), "B" = c(3), "C" = c(3))
  )
  expect_equal(names(.initweibull(df, "time", "y", FALSE)), c("A", "B", "C"))
  expect_equal(names(.initweibull(df, "time", "y", TRUE)), c("I", "A", "B", "C"))
  expect_error(.initweibull(df[1:3, ], "time", "y", FALSE))
})

test_that("init monomolecular works", {
  df <- growthSim("monomolecular",
    n = 20, t = 25,
    params = list("A" = c(100), "B" = c(0.2))
  )
  expect_equal(names(.initmonomolecular(df, "time", "y", FALSE)), c("A", "B"))
  expect_equal(names(.initmonomolecular(df, "time", "y", TRUE)), c("I", "A", "B"))
  expect_error(.initmonomolecular(df[1:3, ], "time", "y", FALSE))
})

test_that("init power law works", {
  df <- growthSim("power law",
    n = 20, t = 25,
    params = list("A" = c(1), "B" = c(1))
  )
  expect_equal(names(.initpowerlaw(df, "time", "y", FALSE)), c("A", "B"))
  expect_equal(names(.initpowerlaw(df, "time", "y", TRUE)), c("I", "A", "B"))
  expect_error(.initpowerlaw(df[1:2, ], "time", "y", FALSE))
})

test_that("init exponential works", {
  df <- growthSim("exponential",
    n = 20, t = 25,
    params = list("A" = c(1), "B" = c(0.2))
  )
  expect_equal(names(.initexponential(df, "time", "y", FALSE)), c("A", "B"))
  expect_equal(names(.initexponential(df, "time", "y", TRUE)), c("I", "A", "B"))
  expect_error(.initexponential(df[1:2, ], "time", "y", FALSE))
})

test_that("init logarithmic works", {
  set.seed(123)
  df <- growthSim("logarithmic",
    n = 20, t = 25,
    params = list("A" = c(15, 15))
  )
  expect_equal(names(.initlogarithmic(df, "time", "y", FALSE)), c("A"))
  expect_equal(names(.initlogarithmic(df, "time", "y", TRUE)), c("I", "A"))
  expect_error(.initlogarithmic(df[1, ], "time", "y", FALSE))
})

test_that("init linear works", {
  set.seed(123)
  df <- growthSim("linear",
    n = 20, t = 25,
    params = list("A" = c(15, 15))
  )
  expect_equal(names(.initlinear(df, "time", "y", FALSE)), c("A"))
  expect_equal(names(.initlinear(df, "time", "y", TRUE)), c("I", "A"))
  expect_error(.initlinear(df[1, ], "time", "y", FALSE))
})

test_that("init bragg works", {
  set.seed(123)
  df <- growthSim("bragg",
    n = 20, t = 100,
    list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
  )
  expect_equal(names(.initbragg(df, "time", "y", FALSE)), c("B", "A", "C"))
  expect_equal(names(.initbragg(df, "time", "y", TRUE)), c("I", "B", "A", "C"))
})

test_that("init lorentz works", {
  set.seed(123)
  df <- growthSim("lorentz",
    n = 20, t = 25,
    params = list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
  )
  expect_equal(names(.initlorentz(df, "time", "y", FALSE)), c("B", "A", "C"))
  expect_equal(names(.initlorentz(df, "time", "y", TRUE)), c("I", "B", "A", "C"))
})

test_that("init beta works", {
  set.seed(123)
  df <- growthSim("beta",
    n = 20, t = 25,
    params = list("A" = 10, "B" = 1.2, "C" = 15, "D" = 8, "E" = 19)
  )
  expect_equal(names(.initbeta(df, "time", "y", FALSE)), c("A", "B", "C", "D", "E"))
  expect_equal(names(.initbeta(df, "time", "y", TRUE)), c("I", "A", "B", "C", "D", "E"))
})

test_that("double logistic form works", {
  out <- .nlrq_form_doublelogistic(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_doublelogistic(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C", "A2", "B2", "C2"))
})

test_that("double gompertz form works", {
  out <- .nlrq_form_doublegompertz(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_doublegompertz(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C", "A2", "B2", "C2"))
})

test_that("logistic form works", {
  out <- .nlrq_form_logistic(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_logistic(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C"))
})

test_that("gompertz form works", {
  out <- .nlrq_form_gompertz(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_gompertz(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C"))
})

test_that("weibull form works", {
  out <- .nlrq_form_weibull(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_weibull(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C"))
})

test_that("frechet form works", {
  out <- .nlrq_form_frechet(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_frechet(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C"))
})

test_that("gumbel form works", {
  out <- .nlrq_form_gumbel(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_gumbel(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C"))
})

test_that("monomolecular form works", {
  out <- .nlrq_form_monomolecular(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_monomolecular(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B"))
})

test_that("exponential form works", {
  out <- .nlrq_form_exponential(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_exponential(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B"))
})

test_that("power law form works", {
  out <- .nlrq_form_powerlaw(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_powerlaw(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B"))
})

test_that("logarithmic form works", {
  out <- .nlrq_form_logarithmic(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_logarithmic(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "I", int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I"))
  out <- .nlrq_form_logarithmic(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = "I", int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I"))
})

test_that("linear form works", {
  out <- .nlrq_form_linear(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_linear(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "I", int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I"))
  out <- .nlrq_form_linear(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = "I", int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I"))
})

test_that("gam form works", {
  out <- .nlrq_form_gam(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_null(out$pars)
  out <- .nlrq_form_gam(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_null(out$pars)
})

test_that("bragg form works", {
  out <- .nlrq_form_bragg(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_bragg(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C"))
})

test_that("lorentz form works", {
  out <- .nlrq_form_lorentz(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_lorentz(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C"))
})

test_that("beta form works", {
  out <- .nlrq_form_beta(
    x = "x", y = "y", USEGROUP = TRUE,
    group = "group", pars = "A", int = FALSE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("A"))
  out <- .nlrq_form_beta(
    x = "x", y = "y", USEGROUP = FALSE,
    group = NULL, pars = NULL, int = TRUE
  )
  expect_equal(names(out), c("formula", "pars"))
  expect_equal(out$pars, c("I", "A", "B", "C", "D", "E"))
})
