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
    priors = list(mu = 40, n = 1, s2 = 100),
    plot = FALSE, rope_range = c(-8, 8), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal", support = seq(20, 100, length.out = 10000)
  )

  expect_equal(out$summary$post.prob, 0.5864103, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 0.7396922, tolerance = 1e-6)

  expect_equal(names(out), c("summary", "posterior"))
})

test_that("conjugate multi value T works", {
  set.seed(123)
  mv_gauss <- mvSim(
    dists = list(
      rnorm = list(mean = 100, sd = 15),
      rnorm = list(mean = 70, sd = 10)
    ),
    n_samples = c(15, 20)
  )

  out <- conjugate(
    s1 = mv_gauss[1:15, -1], s2 = mv_gauss[16:35, -1], method = "t",
    priors = list(mu = 50, n = 1, s2 = 100),
    plot = FALSE, rope_range = c(-5, 5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal", support = NULL
  )

  expect_equal(out$summary$post.prob, 0.001036255, tolerance = 1e-6)

  expect_true(out$summary$rope_prob < 1e-5)

  expect_equal(names(out), c("summary", "posterior"))
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
    priors = list(mu = 40, n = 1, s2 = 100),
    plot = FALSE, rope_range = c(-10, 10), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal", support = NULL
  )


  expect_equal(out$summary$post.prob, 0.9198604, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 0.45927, tolerance = 1e-5)

  expect_equal(names(out), c("summary", "posterior"))
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
    priors = list(mu = 0, n = 1, s2 = 20),
    plot = FALSE, rope_range = c(-5, 5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal", support = NULL
  )

  expect_equal(out$summary$post.prob, 0.7153165, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 0.1824514, tolerance = 0.015)

  expect_equal(names(out), c("summary", "posterior"))
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
    priors = list(a = 0.5, b = 0.5),
    plot = FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )

  expect_equal(out$summary$post.prob, 0.02229246, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 0.1351534, tolerance = 1e-6)

  expect_equal(names(out), c("summary", "posterior"))
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
    priors = list(a = 0.5, b = 0.5),
    plot = FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )

  expect_equal(out$summary$post.prob, 0.1575291, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 0.4059094, tolerance = 0.0001)

  expect_equal(names(out), c("summary", "posterior"))
})

test_that("conjugate single value lognormal works", {
  set.seed(123)
  s1 <- rlnorm(100, log(130), log(1.3))
  s2 <- rlnorm(100, log(100), log(2))
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = "lognormal", priors = list(mu = 3, sd = 3),
    plot = FALSE, rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal", support = NULL
  )

  expect_equal(out$summary$post.prob, 0.3224215, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 0.8657454, tolerance = 1e-6)

  expect_equal(names(out), c("summary", "posterior"))
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
    priors = list(mu = 3, sd = 3),
    plot = FALSE, rope_range = c(-1, 1), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal", support = NULL
  )

  expect_equal(out$summary$post.prob, 0.005472403, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 1, tolerance = 0.0001)

  expect_equal(names(out), c("summary", "posterior"))
})


test_that("conjugate single value poisson works", {
  set.seed(123)
  s1 <- rpois(20, 10)
  s2 <- rpois(20, 8)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "poisson",
    priors = list(a = 0.5, b = 0.5),
    plot = FALSE, rope_range = c(-1, 1), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )

  expect_equal(out$summary$post.prob, 0.09622298, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 0.05594877, tolerance = 1e-6)

  expect_equal(names(out), c("summary", "posterior"))
})


test_that("conjugate single value negative binomial works", {
  set.seed(123)
  s1 <- rnbinom(20, 10, 0.5)
  s2 <- rnbinom(20, 10, 0.25)
  out <- conjugate(
    s1 = s1, s2 = s2, method = "negbin",
    priors = list(r = 10, a = 0.5, b = 0.5),
    plot = FALSE, rope_range = c(-0.5, 0.5), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal"
  )

  expect_equal(out$summary$post.prob, 6.569111e-09, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 1, tolerance = 1e-6)

  expect_equal(names(out), c("summary", "posterior"))
})



test_that("generic conjugate plotting works", {
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
  out <- conjugate(
    s1 = s1, s2 = s2, method = "t",
    priors = list(mu = c(0, 0), n = c(1, 1), s2 = c(20, 20)),
    plot = TRUE, rope_range = c(-8, 8), rope_ci = 0.89,
    cred.int.level = 0.89, hypothesis = "equal", support = NULL
  )

  expect_equal(names(out), c("summary", "posterior", "plot"))

  expect_s3_class(out$plot, "ggplot")
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
    plot = FALSE, rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal", support = NULL
  )

  expect_equal(out$summary$post.prob, 0.6857498, tolerance = 1e-3)

  expect_equal(out$summary$rope_prob, 0.7193574, tolerance = 1e-3)

  expect_equal(unlist(lapply(out$posterior, function(p) {
    names(p)
  })), c("mu", "sd", "lognormal_sigma", "mu", "n", "s2"))

  expect_equal(names(out), c("summary", "posterior"))
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
    cred.int.level = 0.89, hypothesis = "equal", support = NULL
  )

  expect_equal(out$summary$post.prob, 0.167973, tolerance = 1e-6)

  expect_equal(out$summary$rope_prob, 0.28, tolerance = 0.01)

  expect_equal(unlist(lapply(out$posterior, function(p) {
    names(p)
  })), c("mu", "sd", "lognormal_sigma", "mu", "n", "s2"))

  expect_equal(names(out), c("summary", "posterior"))
})

test_that("bivariate conjugate uniform works", {
  set.seed(123)
  s1 <- runif(10, -50, 10)
  s2 <- runif(10, -30, 15)
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = "bivariate_uniform", priors = list(location_l = 0, location_u = 1, scale = 1),
    plot = TRUE, rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal", support = NULL
  )
  expect_s3_class(out$plot, "ggplot")
  expect_equal(nrow(out$summary), 2)
  expect_equal(length(out$posterior), 2)
  expect_equal(names(out$posterior[[1]]), c("scale", "location_l", "location_u"))
  expect_equal(names(out), c("summary", "posterior", "plot"))
})

test_that("bivariate conjugate gaussian works", {
  set.seed(123)
  s1 <- rnorm(10, 20, 5)
  s2 <- rnorm(10, 25, 5)
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = "bivariate_gaussian", priors = list(mu = 20, sd = 10, a = 1, b = 1),
    plot = TRUE, rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal", support = NULL
  )
  expect_s3_class(out$plot, "ggplot")
  expect_equal(nrow(out$summary), 2)
  expect_equal(length(out$posterior), 2)
  expect_equal(names(out$posterior[[1]]), c("mu", "sd", "a", "b"))
  expect_equal(names(out), c("summary", "posterior", "plot"))
})

test_that("bivariate conjugate lognormal works", {
  devtools::load_all("~/pcvr")
  set.seed(123)
  s1 <- rlnorm(10, log(20), 0.25)
  s2 <- rlnorm(10, log(25), 0.4)
  out <- conjugate(
    s1 = s1, s2 = s2,
    method = "bivariate_lognormal", priors = list(mu = 20, sd = 10, a = 1, b = 1),
    plot = TRUE, rope_range = c(-1, 1), rope_ci = 0.89, cred.int.level = 0.89,
    hypothesis = "equal", support = NULL
  )
  expect_s3_class(out$plot, "ggplot")
  expect_equal(nrow(out$summary), 2)
  expect_equal(length(out$posterior), 2)
  expect_equal(names(out$posterior[[1]]), c("mu", "sd", "a", "b"))
  expect_equal(names(out), c("summary", "posterior", "plot"))
})
