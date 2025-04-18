if (!interactive()) pdf(NULL)

test_that("conjugate single value T works", {
  s1 <- c(
    43.8008289810423, 44.6084228775479, 68.9524219823026, 77.442231894233,
    45.2302703709121, 53.8005757403944, 33.8292993277826, 59.7018653972819,
    57.8433869136181, 52.8224097917938
  )
  s2 <- c(
    84.7860854952772, 53.38097452501, 52.352235256613, 49.2369049504088,
    72.7625716991815, 62.6982283802374, 61.2347595388326, 45.298878516913,
    39.6312400911458, 66.9134811003628
  )
  set.seed(123)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "t",
    priors = list(mu = 40, sd = 10),
    rope_range = c(-8, 8), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "unequal",
    bayes_factor = c(50, 55)
  )
  expect_equal(out$summary$post.prob, 0.4099283, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.793057, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
  b <- barg(out, priors = list("rnorm" = list("mean" = c(5, 20), "sd" = c(5, 10))))
  p <- plot(out)
  expect_s3_class(p, "ggplot")
  expect_equal(names(b), c("priorSensitivity", "posteriorPredictive", "Summary"))
  df <- data.frame(value = c(s1, s2), group = rep(c("a", "b"), each = 10))
  out2 <- conjugate(
    value ~ group, df,
    method = "t",
    priors = NULL,
    rope_range = c(-8, 8), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "lesser"
  )
  expect_equal(out2$summary$post.prob, 0.3017588, tolerance = 1e-6)
})

test_that("conjugate single value gaussian works", {
  s1 <- c(
    43.8008289810423, 44.6084228775479, 68.9524219823026, 77.442231894233,
    45.2302703709121, 53.8005757403944, 33.8292993277826, 59.7018653972819,
    57.8433869136181, 52.8224097917938
  )
  s2 <- c(
    62.2654459133762, 66.6863571733485, 61.2951438574251, 62.0014980341704,
    44.0772327229333, 56.169510174076, 71.1378538738675, 55.7547954794673,
    52.4202653287144, 63.3091644583334, 49.263640809148, 63.2460598779059,
    60.3804997092304, 25.1210401427447, 42.6563192857856
  )
  set.seed(123)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "gaussian",
    priors = NULL,
    rope_range = c(-10, 10), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = c(50, 55)
  )
  expect_s3_class(out, "conjugate")
})

test_that("conjugate single value beta works", {
  s1 <- c(
    0.703503634559075, 0.511635134690101, 0.632190260438834, 0.353931888757961,
    0.221420395118864, 0.513410161844891, 0.34422446087923, 0.377469570817219,
    0.553479714127415, 0.610722823397796, 0.642879912798542, 0.276682168393891,
    0.469471347132478, 0.444690242008423, 0.21701860450406, 0.362069754559641,
    0.324136767421681, 0.776072763733466, 0.678539925827321, 0.230328895808406
  )
  s2 <- c(
    0.609522707465377, 0.485478489613084, 0.38771181119892, 0.605383942021798,
    0.657651014793296, 0.650853576978649, 0.556595652577583, 0.719495963439006,
    0.695092908314064, 0.769699311988927, 0.665722925603491, 0.651953500947315,
    0.523204344509775, 0.736962689122744, 0.607863983829372, 0.703483180709229,
    0.418954761865872, 0.556556033829364, 0.617343726804053, 0.522669623038004
  )
  set.seed(123)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "beta",
    priors = NULL,
    rope_range = c(-0.1, 0.1), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = c(0.5)
  )
  b <- barg(out)
  expect_equal(names(b), c("posteriorPredictive", "Summary"))
  expect_equal(out$summary$post.prob, 0.02229246, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.1351534, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
  expect_error(conjugate(s1 = c(s1, -0.1), s2 = c(s2, 1.1), method = "beta"))
})

test_that("conjugate single value lognormal works", {
  set.seed(123)
  s1 <- rlnorm(100, log(130), log(1.3))
  s2 <- rlnorm(100, log(100), log(2))
  expect_warning(
    out <- conjugate(
      s1 = s1, s2 = s2,
      method = "lognormal", priors = NULL, plot = TRUE,
      rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
      hypothesis = "equal", bayes_factor = 5
    )
  )
  expect_equal(out$summary$post.prob, 0.5980101, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.1113358, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
})

test_that("conjugate single value lognormal2 works", {
  set.seed(123)
  s1 <- rlnorm(100, log(130), log(1.3))
  s2 <- rlnorm(100, log(100), log(2))
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = "lognormal2", priors = NULL,
    rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal", bayes_factor = 125
  )
  p <- plot(out)
  expect_s3_class(p, "ggplot")
  expect_equal(out$summary$post.prob, 1.069935e-09, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
})

test_that("conjugate single value poisson works", {
  set.seed(123)
  s1 <- rpois(20, 10)
  s2 <- rpois(20, 8)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "poisson",
    priors = NULL,
    rope_range = c(-1, 1), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = 9
  )
  expect_equal(out$summary$post.prob, 0.09622298, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.05594877, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
  expect_error(
    conjugate(
      s1 = c(s1, 1.5), s2 = s2, method = "poisson",
      priors = NULL,
      rope_range = c(-1, 1), rope_ci = 0.89,
      cred.int.level = 0.89, hypothesis = "equal"
    )
  )
})

test_that("conjugate single value negative binomial works", {
  set.seed(123)
  s1 <- rnbinom(20, 10, 0.5)
  s2 <- rnbinom(20, 10, 0.25)
  suppressWarnings(
    out <- conjugate(
      s1 = s1, s2 = s2, method = "negbin",
      priors = NULL,
      rope_range = c(-0.5, 0.5), rope_ci = 0.89,
      cred.int.level = 0.89, hypothesis = "equal",
      bayes_factor = 10
    )
  )
  expect_equal(out$summary$post.prob, 6.569111e-09, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 1, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
  expect_error(
    conjugate(
      c(0.5, 0.1, 1, 1.1),
      method = "negbin"
    )
  )
})

test_that("conjugate single value binomial works", {
  set.seed(123)
  s1 <- list(successes = c(15, 14, 16, 11), trials = 20)
  s2 <- list(successes = c(10, 9, 12, 10), trials = 20)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "binomial",
    priors = NULL,
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = 0.5
  )
  expect_equal(out$summary$post.prob, 0.08529131, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 1, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
  expect_error(
    .conj_binomial_formatter(c(1, -1))
  )
  expect_error(
    .conj_binomial_formatter(c(1, 1, 1))
  )
  expect_message(
    out <- .conj_binomial_formatter(
      list(
        c(1, 1, 4),
        c(3, 3, 10)
      )
    )
  )
  expect_equal(names(out), c("successes", "trials"))
})

test_that("conjugate single value bernoulli works", {
  set.seed(123)
  s1 <- sample(c(TRUE, FALSE), 10, replace = TRUE, prob = c(0.4, 0.6))
  s2 <- sample(c(TRUE, FALSE), 10, replace = TRUE, prob = c(0.7, 0.3))
  out <- conjugate(
    s1 = s1, s2 = s2, method = "bernoulli",
    priors = NULL,
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = 0.75
  )
  expect_equal(out$summary$post.prob, 0.3412209, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.914504, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
  expect_error(
    conjugate(c(1, 2, 3), method = "bernoulli")
  )
})

test_that("conjugate single value pareto works", {
  set.seed(123)
  s1 <- extraDistr::rpareto(10, 2, 1)
  s2 <- extraDistr::rpareto(10, 3, 1)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "pareto",
    priors = NULL,
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = 3
  )
  expect_equal(out$summary$post.prob, 0.8643824, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.01584092, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
})

test_that("conjugate single value uniform works", {
  set.seed(123)
  s1 <- runif(10, 0, 10)
  s2 <- runif(10, 0, 13)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "uniform",
    priors = NULL,
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = c(2, 3)
  )
  expect_equal(out$summary$post.prob, 0.05305783, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
})

test_that("conjugate single value von mises (1) works", {
  set.seed(123)
  s1 <- brms::rvon_mises(10, 0, 2)
  s2 <- brms::rvon_mises(10, 0, 2)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "vonmises",
    priors = list(mu = 0, kappa = 0.5, boundary = c(-pi, pi)),
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = 1
  )
  expect_equal(out$summary$post.prob, 0.4736915, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.255814, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
  out2 <- conjugate(
    s1 = s1, s2 = s2, method = "vonmises",
    priors = list(mu = 0),
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = c(1, 2)
  )
  expect_error(conjugate(
    s1 = rnorm(10, 10, 1), s2 = rnorm(10, 10, 1), method = "vonmises",
    priors = NULL,
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  ))
})

test_that("conjugate single value von mises (2) works", {
  set.seed(123)
  s1 <- rnorm(10, 50, 6)
  s2 <- rnorm(10, 60, 10)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "vonmises2",
    priors = list(mu = 0, boundary = c(0, 110)),
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = 55
  )
  expect_equal(out$summary$post.prob, 0.4529312, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.01999775, tolerance = 1e-3)
  expect_s3_class(out, "conjugate")
  expect_error(conjugate(
    s1 = s1, s2 = s2, method = "vonmises2",
    priors = NULL,
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  ))
})

test_that("conjugate single value gamma works", {
  set.seed(123)
  s1 <- rgamma(10, 2, 1)
  s2 <- rgamma(10, 1, 1)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "gamma",
    priors = NULL,
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = 2
  )
  expect_equal(out$summary$post.prob, 0.1474759, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.2627795, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
})

test_that("conjugate single value exponential works", {
  set.seed(123)
  s1 <- rexp(10, 1.2)
  s2 <- rexp(10, 1)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "exponential",
    priors = NULL,
    rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal",
    bayes_factor = c(1, 1.5)
  )
  expect_equal(out$summary$post.prob, 0.3536306, tolerance = 1e-6)
  expect_equal(out$summary$rope_prob, 0.3370408, tolerance = 1e-6)
  expect_s3_class(out, "conjugate")
})

test_that("conjugate single value lognormal vs gaussian", {
  set.seed(123)
  s1 <- rlnorm(100, log(70), log(2))
  s2 <- rnorm(100, 4, 1)
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = c("lognormal", "gaussian"), priors = list(
      list(mu = 3, sd = 5),
      list(mu = 5, n = 1, s2 = 2)
    ),
    rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal"
  )

  expect_equal(out$summary$post.prob, 0.06048735, tolerance = 1e-3)

  expect_equal(out$summary$rope_prob, 0.2325581, tolerance = 1e-3)

  expect_equal(unlist(lapply(out$posterior, function(p) {
    return(names(p))
  })), c("mu", "sd", "lognormal_sigma", "mu", "sd"))

  expect_s3_class(out, "conjugate")
})

test_that("single value bivariate conjugate uniform works", {
  set.seed(123)
  s1 <- runif(10, 1, 10)
  s2 <- runif(10, -1, 15)
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = "bivariate_uniform", priors = NULL,
    rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal"
  )
  expect_s3_class(plot(out), "ggplot")
  expect_equal(nrow(out$summary), 2)
  expect_equal(length(out$posterior), 2)
  expect_equal(names(out$posterior[[1]]), c("scale", "location_l", "location_u"))
  expect_s3_class(out, "conjugate")
  out2 <- conjugate(
    s1 = s1,
    method = "bivariate_uniform", priors = NULL,
    rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal"
  )
  expect_s3_class(out, "conjugate")
  set.seed(123)
  s1 <- runif(10, -15, -7)
  s2 <- runif(10, -10, -5)
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = "bivariate_uniform", priors = list(location_l = -10, location_u = -8, scale = 1),
    rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal"
  )
  expect_s3_class(plot(out), "ggplot")
  expect_equal(nrow(out$summary), 2)
  expect_equal(length(out$posterior), 2)
  expect_equal(names(out$posterior[[1]]), c("scale", "location_l", "location_u"))
  expect_s3_class(out, "conjugate")
})

test_that("bivariate conjugate gaussian works", {
  set.seed(123)
  s1 <- rnorm(10, 20, 5)
  s2 <- rnorm(10, 25, 5)
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = "bivariate_gaussian", priors = NULL,
    rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal"
  )
  expect_s3_class(plot(out), "ggplot")
  expect_equal(nrow(out$summary), 2)
  expect_equal(length(out$posterior), 2)
  expect_equal(names(out$posterior[[1]]), c("mu", "sd", "a", "b"))
  expect_s3_class(out, "conjugate")
})

test_that("bivariate conjugate lognormal works", {
  set.seed(123)
  s1 <- rlnorm(10, log(20), 0.25)
  s2 <- rlnorm(10, log(25), 0.4)
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = "bivariate_lognormal", priors = NULL,
    rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal"
  )
  expect_s3_class(plot(out), "ggplot")
  expect_equal(nrow(out$summary), 2)
  expect_equal(length(out$posterior), 2)
  expect_equal(names(out$posterior[[1]]), c("mu", "sd", "a", "b"))
  expect_s3_class(out, "conjugate")
  out2 <- conjugate(
    s1 = s1,
    method = "bivariate_lognormal", priors = NULL,
    rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89
  )
  p <- plot(out2)
  expect_s3_class(p, "ggplot")
})
