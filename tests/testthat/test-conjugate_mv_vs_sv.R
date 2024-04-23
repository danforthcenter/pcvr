if (!interactive()) pdf(NULL)

test_that("conjugate vonmises method is consistent for SV and MV", {
  set.seed(123)
  method <- "vonmises"
  prior <- list(mu = 0, kappa = 1, known_kappa = 1, boundary = c(0, 180), n = 1)
  generating <- list(
    s1 = list(f = "rnorm", n = 20, mean = 20, sd = 5),
    s2 = list(f = "rnorm", n = 20, mean = 170, sd = 5)
    )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  mu1 <- as.numeric(out$posteriors[c(1,3), "mu"])
  mu2 <- as.numeric(out$posteriors[c(2,4), "mu"])
  expect_equal(
    abs(diff(mu1)/mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.1 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2)/mean(mu2)),
    0,
    tolerance = 0.1
  )
  k1 <- as.numeric(out$posteriors[c(1,3), "kappa"])
  k2 <- as.numeric(out$posteriors[c(2,4), "kappa"])
  expect_equal(
    abs(diff(k1)/mean(k1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.5 # but these have more variability
  )
  expect_equal(
    abs(diff(k2)/mean(k2)),
    0,
    tolerance = 0.5
  )
})

test_that("conjugate vonmises2 method is consistent for SV and MV", {
  set.seed(123)
  method <- "vonmises2"
  prior <- list(mu = 0, kappa = 1, known_kappa = 1, boundary = c(0, 180), n = 1)
  generating <- list(
    s1 = list(f = "rnorm", n = 20, mean = 20, sd = 5),
    s2 = list(f = "rnorm", n = 20, mean = 170, sd = 5)
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  mu1 <- as.numeric(out$posteriors[c(1,3), "mu"])
  mu2 <- as.numeric(out$posteriors[c(2,4), "mu"])
  expect_equal(
    abs(diff(mu1)/mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.1 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2)/mean(mu2)),
    0,
    tolerance = 0.1
  )
  k1 <- as.numeric(out$posteriors[c(1,3), "kappa"])
  k2 <- as.numeric(out$posteriors[c(2,4), "kappa"])
  expect_equal(
    abs(diff(k1)/mean(k1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.5 # but these have more variability
  )
  expect_equal(
    abs(diff(k2)/mean(k2)),
    0,
    tolerance = 0.5
  )
})

test_that("conjugate T method is consistent for SV and MV", {
  set.seed(324)
  method <- "t"
  prior <- list(mu=c(50,50), n=c(1,1), s2=c(100,100))
  generating <- list(
    s1 = list(f = "rnorm", n = 10, mean = 50, sd = 10),
    s2 = list(f = "rnorm", n = 10, mean = 60, sd = 15)
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  mu1 <- as.numeric(out$posteriors[c(1,3), "mu"])
  mu2 <- as.numeric(out$posteriors[c(2,4), "mu"])
  expect_equal(
    abs(diff(mu1)/mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.1 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2)/mean(mu2)),
    0,
    tolerance = 0.1
  )
  v1 <- as.numeric(out$posteriors[c(1,3), "s2"])
  v2 <- as.numeric(out$posteriors[c(2,4), "s2"])
  expect_equal(
    abs(diff(v1)/mean(v1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.25 # but these have more variability
  )
  expect_equal(
    abs(diff(v2)/mean(v2)),
    0,
    tolerance = 0.25
  )
})

test_that("conjugate Gaussian method is consistent for SV and MV", {
  set.seed(324)
  method <- "gaussian"
  prior <- list(mu=c(50,50), n=c(1,1), s2=c(100,100))
  generating <- list(
    s1 = list(f = "rnorm", n = 10, mean = 50, sd = 10),
    s2 = list(f = "rnorm", n = 10, mean = 60, sd = 15)
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  mu1 <- as.numeric(out$posteriors[c(1,3), "mu"])
  mu2 <- as.numeric(out$posteriors[c(2,4), "mu"])
  expect_equal(
    abs(diff(mu1)/mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.1 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2)/mean(mu2)),
    0,
    tolerance = 0.1
  )
  v1 <- as.numeric(out$posteriors[c(1,3), "s2"])
  v2 <- as.numeric(out$posteriors[c(2,4), "s2"])
  expect_equal(
    abs(diff(v1)/mean(v1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.25 # but these have more variability
  )
  expect_equal(
    abs(diff(v2)/mean(v2)),
    0,
    tolerance = 0.25
  )
})


test_that("conjugate Beta method is consistent for SV and MV", {
  # MV traits are stronger and this has very loose tolerances. I'd like to tighten it up, but it is
  # not obvious how much the extra information in MV traits should change this vs how much it does
  # change this. These tests are one step towards making that very stable, but it is not totally done.
  set.seed(324)
  method <- "beta"
  prior <- list( a=c(0.5, 0.5), b=c(0.5, 0.5) )
  generating <- list(
    s1 = list(f = "rbeta", n = 50, shape1 = 8, shape2 = 4),
    s2 = list(f = "rbeta", n = 50, shape1 = 7, shape2 = 6)
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  out$posteriors
  out$summaries
  a1 <- as.numeric(out$posteriors[c(1,3), "a"])
  a2 <- as.numeric(out$posteriors[c(2,4), "a"])
  expect_equal(
    abs(diff(a1)/mean(a1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.5 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(a2)/mean(a2)),
    0,
    tolerance = 0.5
  )
  b1 <- as.numeric(out$posteriors[c(1,3), "b"])
  b2 <- as.numeric(out$posteriors[c(2,4), "b"])
  expect_equal(
    abs(diff(b1)/mean(b1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.5
  )
  expect_equal(
    abs(diff(b2)/mean(b2)),
    0,
    tolerance = 0.5
  )
})

test_that("conjugate lognormal method is consistent for SV and MV", {
  set.seed(324)
  method <- "lognormal"
  prior <- list( mu_log=c(log(10),log(10)),n=c(1,1), sigma_log=c(log(3),log(3)) )
  generating <- list(
    s1 = list(f = "rlnorm", n = 20, meanlog = 3.5, sdlog = log(3)),
    s2 = list(f = "rlnorm", n = 20, meanlog = 3, sdlog = log(4))
  )
  out <- .conjugate.mv.sv.testing(method, prior, generating)
  out$posteriors
  out$summaries
  mu1 <- as.numeric(out$posteriors[c(1,3), "mu_log"])
  mu2 <- as.numeric(out$posteriors[c(2,4), "mu_log"])
  expect_equal(
    abs(diff(mu1)/mean(mu1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.15 # but give it some wiggle room per randomness
  )
  expect_equal(
    abs(diff(mu2)/mean(mu2)),
    0,
    tolerance = 0.15
  )
  v1 <- as.numeric(out$posteriors[c(1,3), "sigma_log"])
  v2 <- as.numeric(out$posteriors[c(2,4), "sigma_log"])
  expect_equal(
    abs(diff(v1)/mean(v1)),
    0, # diff should be 0, scaled or not
    tolerance = 0.25 # but these have more variability
  )
  expect_equal(
    abs(diff(v2)/mean(v2)),
    0,
    tolerance = 0.25
  )
})


