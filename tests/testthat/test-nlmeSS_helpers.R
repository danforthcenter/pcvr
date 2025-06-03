library(testthat)
library(pcvr)

test_that("nlme double logistic formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_doublelogistic(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_doublelogistic(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_doublelogistic(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme double gompertz formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_doublegompertz(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_doublegompertz(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_doublegompertz(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme logistic formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_logistic(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_logistic(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_logistic(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme logistic4 formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_logistic4(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_logistic4(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_logistic4(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme logistic5 formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_logistic5(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_logistic5(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_logistic5(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})


test_that("nlme gompertz formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_gompertz(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_gompertz(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_gompertz(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme frechet formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_frechet(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_frechet(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_frechet(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme gumbel formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_gumbel(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_gumbel(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_gumbel(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme weibull formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_weibull(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_weibull(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_weibull(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme monomolecular formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_monomolecular(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_monomolecular(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_monomolecular(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme exponential formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_exponential(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_exponential(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_exponential(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme power law formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_powerlaw(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_powerlaw(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_powerlaw(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme linear formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_linear(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_linear(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_linear(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme logarithmic formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_logarithmic(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_logarithmic(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_logarithmic(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme bragg formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_bragg(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_bragg(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_bragg(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme lorentz formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_lorentz(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_lorentz(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_lorentz(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme beta formula", {
  # x, y, group, individual, matched_sigma, pars, int
  out <- .nlme_form_beta(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_beta(
    x = "x", y = "y", group = "group",
    individual = "id", matched_sigma = "power",
    pars = "A", int = TRUE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
  out <- .nlme_form_beta(
    x = "x", y = "y", group = "dummyGroup",
    individual = "id", matched_sigma = "power",
    pars = NULL, int = FALSE
  )
  expect_equal(names(out), c("model", "random", "fixed", "groups", "weights", "cor_form", "pars"))
})

test_that("nlme sigma options work", {
  out <- .nlme_sigma_form(matched_sigma = "int", x = "x", group = "group")
  expect_s3_class(out, "varFunc")
  out <- .nlme_sigma_form(matched_sigma = "power", x = "x", group = "group")
  expect_s3_class(out, "varFunc")
  out <- .nlme_sigma_form(matched_sigma = "exp", x = "x", group = "group")
  expect_s3_class(out, "varFunc")
  out <- .nlme_sigma_form(matched_sigma = nlme::varFunc(~x), x = "x", group = "group")
  expect_s3_class(out, "varFunc")
})
