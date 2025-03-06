if (!interactive()) pdf(NULL)

test_that("conjugate multi value T works", {
  set.seed(123)
  mv_gauss <- mvSim(
    dists = list(
      rnorm = list(mean = 100, sd = 15),
      rnorm = list(mean = 70, sd = 10)
    ),
    n_samples = c(15, 20)
  )
  mv_gauss$group <- rep(c("a", "b"), times = c(15, 20))
  out <- conjugate(
    2:181 ~ group, mv_gauss,
    method = "t",
    priors = NULL,
    plot = TRUE, rope_range = c(-5, 5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )
  expect_equal(out$summary$post.prob, 0.02475074, tolerance = 1e-6)
  expect_true(out$summary$rope_prob < 1e-5)
  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))
})

test_that("conjugate multi value gaussian works", {
  set.seed(123)
  mv_gauss <- mvSim(
    dists = list(
      rnorm = list(mean = 50, sd = 10),
      rnorm = list(mean = 60, sd = 12)
    ),
    n_samples = c(30, 40)
  )

  out <- conjugate(
    s1 = mv_gauss[1:30, -1], s2 = mv_gauss[31:70, -1], method = "gaussian",
    priors = NULL,
    plot = TRUE, rope_range = c(-5, 5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )

  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))
})

test_that("conjugate multi value beta works", {
  set.seed(123)
  mv_beta <- mvSim(
    dists = list(
      rbeta = list(shape1 = 5, shape2 = 8),
      rbeta = list(shape1 = 10, shape2 = 10)
    ),
    n_samples = c(10, 10)
  )

  out <- conjugate(
    s1 = mv_beta[1:10, -1], s2 = mv_beta[11:20, -1], method = "beta",
    priors = NULL,
    plot = TRUE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )

  expect_equal(out$summary$post.prob, 0.1575291, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.4059094, tolerance = 0.0001)
  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))
  colnames(mv_beta)[ncol(mv_beta)] <- "sim_500"
  expect_error(conjugate(s1 = mv_beta[1:10, -1], s2 = mv_beta[11:20, -1], method = "beta"))
})

test_that("conjugate multi value lognormal works", {
  set.seed(123)
  mv_ln <- mvSim(
    dists = list(
      rlnorm = list(meanlog = log(50), sdlog = log(1.7)),
      rlnorm = list(meanlog = log(30), sdlog = log(2.1))
    ),
    n_samples = 30
  )
  out <- conjugate(
    s1 = mv_ln[1:30, -1], s2 = mv_ln[31:60, -1], method = "lognormal",
    priors = NULL,
    plot = TRUE, rope_range = c(-1, 1), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )
  expect_equal(out$summary$post.prob, 0.09558071, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 1, tolerance = 0.0001)
  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))
})

test_that("conjugate multi value lognormal2 works", {
  set.seed(123)
  mv_ln <- mvSim(
    dists = list(
      rlnorm = list(meanlog = log(50), sdlog = log(1.7)),
      rlnorm = list(meanlog = log(30), sdlog = log(2.1))
    ),
    n_samples = 30
  )
  out <- conjugate(
    s1 = mv_ln[1:30, -1], s2 = mv_ln[31:60, -1], method = "lognormal2",
    priors = NULL,
    plot = TRUE, rope_range = c(-1, 1), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )
  expect_equal(out$summary$post.prob, 0.3054448, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.007976632, tolerance = 0.0001)
  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))
})

test_that("conjugate multi value pareto works", {
  set.seed(123)
  library(extraDistr)
  mv <- mvSim(
    dists = list(
      rpareto = list(a = 1, b = 1),
      rpareto = list(a = 1, b = 1)
    ),
    n_samples = 30, counts = 100
  )
  out <- conjugate(
    s1 = mv[1:30, -1], s2 = mv[31:60, -1], method = "pareto",
    priors = NULL,
    plot = TRUE, rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )
  expect_equal(out$summary$post.prob, 0.9978403, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.001348163, tolerance = 1e-6)
  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))
})

test_that("conjugate multi value uniform works", {
  set.seed(123)
  mvu <- mvSim(
    dists = list(
      runif = list(min = 0, max = 100),
      runif = list(min = 0, max = 90)
    ),
    n_samples = 30
  )
  out <- conjugate(
    s1 = mvu[1:30, -1], s2 = mvu[31:60, -1], method = "uniform",
    priors = NULL,
    plot = TRUE, rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )
  expect_equal(out$summary$post.prob, 0.03885159, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0, tolerance = 1e-6)
  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))
})

test_that("conjugate multi value von mises (1) works", {
  set.seed(123)
  mv <- mvSim(
    dists = list(
      rnorm = list(mean = 50, sd = 10),
      rnorm = list(mean = 40, sd = 8)
    ),
    n_samples = 30, counts = 1000
  )
  out <- conjugate(
    s1 = mv[1:30, -1], s2 = mv[31:60, -1], method = "vonmises",
    priors = list(boundary = c(0, 100)),
    plot = TRUE, rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )
  expect_equal(out$summary$post.prob, 0.562518, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.0217953, tolerance = 1e-6)
  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))
  expect_error(conjugate(
    s1 = mv[1:30, -1], s2 = mv[31:60, -1], method = "vonmises",
    priors = NULL,
    plot = FALSE, rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  ))
})

test_that("conjugate multi value von mises (2) works", {
  set.seed(123)
  mv <- mvSim(
    dists = list(
      rnorm = list(mean = 50, sd = 10),
      rnorm = list(mean = 40, sd = 8)
    ),
    n_samples = 30, counts = 1000
  )
  out <- conjugate(
    s1 = mv[1:30, -1], s2 = mv[31:60, -1], method = "vonmises2",
    priors = list(boundary = c(0, 100)),
    plot = TRUE, rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )
  expect_equal(out$summary$post.prob, 0.5684183, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.02471632, tolerance = 1e-6)
  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))
  expect_error(conjugate(
    s1 = mv[1:30, -1], s2 = mv[31:60, -1], method = "vonmises2",
    priors = NULL,
    plot = FALSE, rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  ))
})

test_that("conjugate multi value lognormal vs gaussian", {
  set.seed(123)
  mv_ln <- mvSim(
    dists = list(
      rlnorm = list(meanlog = log(50), sdlog = log(1.7)),
      rnorm = list(mean = 5, sd = 1)
    ),
    n_samples = 30
  )
  out <- conjugate(
    s1 = mv_ln[1:30, -1], s2 = mv_ln[31:60, -1], method = c("lognormal", "gaussian"),
    priors = list(
      list(mu = 5, sd = 4),
      list(mu = 5, n = 1, s2 = 4)
    ),
    plot = FALSE, rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )

  expect_equal(out$summary$post.prob, 0.167973, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 0.28, tolerance = 0.01)

  expect_equal(unlist(lapply(out$posterior, function(p) {
    return(names(p))
  })), c("mu", "sd", "lognormal_sigma", "mu", "n", "s2"))

  expect_equal(names(out), c("summary", "posterior", "prior"))
})

test_that("bivariate conjugate multi value uniform works", {
  set.seed(123)
  mv <- mvSim(
    dists = list(
      runif = list(min = 0, max = 125),
      runif = list(min = 0, max = 150)
    ),
    n_samples = 10,
    counts = 1000,
    min_bin = 1,
    max_bin = 180,
    wide = TRUE
  )
  out <- conjugate(
    s1 = mv[1:10, -1],
    s2 = mv[11:20, -1],
    method = "bivariate_uniform", priors = NULL,
    plot = TRUE, rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal"
  )
  expect_s3_class(out$plot, "ggplot")
  expect_equal(nrow(out$summary), 2)
  expect_equal(length(out$posterior), 2)
  expect_equal(names(out$posterior[[1]]), c("scale", "location_l", "location_u"))
  expect_equal(names(out), c("summary", "posterior", "prior", "plot"))

  set.seed(123)
  mv <- mvSim(
    dists = list(
      runif = list(min = -100, max = -10),
      runif = list(min = -90, max = -15)
    ),
    n_samples = 10,
    counts = 1000,
    min_bin = -100,
    max_bin = 0,
    wide = TRUE
  )
  out2 <- conjugate(
    s1 = mv[1:10, -1],
    s2 = mv[11:20, -1],
    method = "bivariate_uniform", priors = list(location_l = -50, location_u = -45, scale = 1),
    plot = TRUE, rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal"
  )
  expect_s3_class(out2$plot, "ggplot")
  expect_equal(nrow(out2$summary), 2)
  expect_equal(length(out2$posterior), 2)
  expect_equal(names(out2$posterior[[1]]), c("scale", "location_l", "location_u"))
  expect_equal(names(out2), c("summary", "posterior", "prior", "plot"))
})
