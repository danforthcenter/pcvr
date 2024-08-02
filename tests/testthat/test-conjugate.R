if (!interactive()) pdf(NULL)

test_that("conjugate HDE helpers work", {
  expect_equal(.betaHDE(1, 2), 0)
  expect_equal(.betaHDE(2, 1), 1)
  expect_equal(.betaHDE(10, 10), 0.5)
  expect_equal(.gammaHDE(1, 1), 0)
  expect_equal(.gammaHDE(0, 1), 0)
  expect_equal(.gammaHDE(10, 10), 90)
})

test_that("conjugate MV Sample formatting helper works", {
  mv1 <- matrix(1, 10, 10)
  expect_warning(out1 <- .mvSampleFormatting(mv1))
  expect_equal(colnames(out1), paste0("b", 1:10))
})

test_that("conjugate raises expected errors", {
  set.seed(123)
  s1 <- rnorm(10, 10, 1)
  res1 <- conjugate(
    s1 = s1, method = "t",
    priors = NULL,
    plot = TRUE,
    cred.int.level = 0.89, hypothesis = "bad", rope_range = c(-1, 1)
  )
  expect_equal(res1$summary$HDE_1, 9.158751, tolerance = 1e-6)
  s2 <- rnorm(10, 15, 2)
  expect_error(conjugate(
    s1 = s1, s2 = s2, method = "t",
    priors = NULL,
    plot = TRUE,
    cred.int.level = 0.89, hypothesis = "bad", rope_range = c(-1, 1)
  ))
  expect_error(conjugate(
    s1 = s1, s2 = s2, method = "t",
    priors = NULL,
    plot = TRUE,
    cred.int.level = 0.89, hypothesis = "lesser", rope_range = 1
  ))
  df <- data.frame(value = c(s1, s2), group = rep(c("a", "b"), each = 10))
  expect_error(conjugate(
    value ~ group, s2, method = "t",
    priors = NULL,
    plot = TRUE,
    cred.int.level = 0.89, hypothesis = "lesser", rope_range = 1
  ))
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
