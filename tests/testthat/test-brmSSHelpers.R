library(testthat)
library(pcvr)

test_that("prior specification in brmSSHelpers works for all inputs", {
  set.seed(123)
  simdf <- growthSim(
    "logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  ss <- growthSS(
    model = "logistic", form = y ~ time | id / group, sigma = "gam",
    start = NULL,
    df = simdf, type = "brms"
  )
  expect_s3_class(ss$prior, "brmsprior")
  ss <- growthSS(
    model = "logistic", form = y ~ time | id / group, sigma = "gam",
    start = brms::set_prior("normal(0, 10)"),
    df = simdf, type = "brms"
  )
  expect_s3_class(ss$prior, "brmsprior")
  expect_warning(
    ss <- growthSS(
      model = "logistic", form = y ~ time | id / group, sigma = "gam",
      start = c(200, 10, 3),
      df = simdf, type = "brms"
    )
  )
  expect_s3_class(ss$prior, "brmsprior")
  expect_error(
    ss <- growthSS(
      model = "logistic", form = y ~ time | id / group, sigma = "gam",
      start = c(200, 10),
      df = simdf, type = "brms"
    )
  )
  ss <- growthSS(
    model = "logistic", form = y ~ time | id / group, sigma = "gam",
    list("A" = c(130, 120), "B" = c(11, 10), "C" = c(3, 5)),
    df = simdf, type = "brms"
  )
  expect_s3_class(ss$prior, "brmsprior")
  expect_warning(
    ss <- growthSS(
      model = "logistic", form = y ~ time | id / group, sigma = "gam",
      list(c(130, 120), c(11, 10), c(3, 5)),
      df = simdf, type = "brms"
    )
  )
  expect_s3_class(ss$prior, "brmsprior")
  expect_warning(
    out <- .makePriors(
      priors = list("A" = 100, "B" = 1),
      pars = c("A", "B", "C"), df = data.frame(group = "a"), group = "group", USEGROUP = TRUE,
      sigma = NULL, family = "student", formula = y ~ x
    )
  )
  expect_s3_class(out, "brmsprior")
})

test_that("sigma helper works", {
  test1 <- .sigmaHelper(
    sigma = y ~ x,
    dpars = c("sigma", "nu"),
    family = "student",
    models = NULL
  )
  expect_equal(names(test1), c("sigma", "nu"))
  expect_error(
    .sigmaHelper(
      sigma = c("linear", "linear", "linear"),
      dpars = c("sigma", "nu"),
      family = "student"
    )
  )
})

test_that("brms form logistic works", {
  out <- .brms_form_logistic(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_logistic(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_logistic(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_logistic(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form logistic4 works", {
  out <- .brms_form_logistic4(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_logistic4(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_logistic4(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_logistic4(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form logistic5 works", {
  out <- .brms_form_logistic5(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_logistic5(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_logistic5(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_logistic5(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form gompertz works", {
  out <- .brms_form_gompertz(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_gompertz(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_gompertz(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_gompertz(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form gumbel works", {
  out <- .brms_form_gumbel(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_gumbel(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_gumbel(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_gumbel(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form frechet works", {
  out <- .brms_form_frechet(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_frechet(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_frechet(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_frechet(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form weibull works", {
  out <- .brms_form_weibull(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_weibull(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_weibull(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_weibull(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form double logistic works", {
  out <- .brms_form_doublelogistic(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_doublelogistic(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_doublelogistic(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_doublelogistic(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form double gompertz works", {
  out <- .brms_form_doublegompertz(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_doublegompertz(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_doublegompertz(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_doublegompertz(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form monomolecular works", {
  out <- .brms_form_monomolecular(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_monomolecular(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_monomolecular(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_monomolecular(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form exponential works", {
  out <- .brms_form_exponential(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_exponential(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_exponential(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_exponential(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form power law works", {
  out <- .brms_form_powerlaw(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_powerlaw(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_powerlaw(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_powerlaw(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form linear works", {
  out <- .brms_form_linear(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_linear(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_linear(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_linear(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
  out3 <- .brms_form_linear(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = list("sigmaA" = 1), int = FALSE
  )
  expect_equal(names(out3), c("form", "pars"))
})

test_that("brms form logarithmic works", {
  out <- .brms_form_logarithmic(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_logarithmic(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_logarithmic(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_logarithmic(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form gam works", {
  out <- .brms_form_gam(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = 25, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_gam(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = 5, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_gam(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = 25, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_gam(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = 5, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form bragg works", {
  out <- .brms_form_bragg(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_bragg(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_bragg(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_bragg(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form lorentz works", {
  out <- .brms_form_lorentz(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_lorentz(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_lorentz(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_lorentz(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("brms form beta works", {
  out <- .brms_form_beta(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out), c("form", "pars"))
  out <- .brms_form_beta(
    x = "x", y = "y", group = "group",
    dpar = FALSE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out), c("form", "pars"))
  out2 <- .brms_form_beta(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = FALSE
  )
  expect_equal(names(out2), c("form", "pars"))
  out2 <- .brms_form_beta(
    x = "x", y = "sigma", group = "group",
    dpar = TRUE, nTimes = NULL, useGroup = TRUE,
    prior = NULL, int = TRUE
  )
  expect_equal(names(out2), c("form", "pars"))
})

test_that("decay works with intercept model", {
  formList <- list("form" = y ~ I + A * x)
  out <- .brms_form_decay(formList, int = TRUE)
  expect_equal(as.character(out$form), as.character(y ~ I - (A * x)))
})
