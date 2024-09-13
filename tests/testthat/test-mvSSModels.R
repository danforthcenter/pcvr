library(testthat)

test_that("Test spectral index helpers", {
  indices <- c(
    "none", "ari", "ci_rededge", "cri550", "cri700",
    "egi", "evi", "gdvi", "mari", "mcari", "mtci", "ndre",
    "ndvi", "pri", "psnd_chlorophyll_a", "psnd_chlorophyll_b",
    "psnd_caroteniods", "psri", "pssr_chlorophyll_a",
    "pssr_chlorophyll_b", "pssr_caroteniods", "rgri",
    "rvsi", "savi", "sipi", "sr", "vari", "vi_green", "wi",
    "fvfm", "fqfm"
  )
  for (index in indices) {
    fun <- get(paste0(".", index, "_mvss_hlp"))
    expect_equal(names(fun()), c("trunc", "family"))
  }
})

#* `Non-Longitudinal Multi-Value Trait Models`

set.seed(123)
mv_df <- mvSim(dists = list(rnorm = list(mean = 100, sd = 30)), wide = FALSE)
mv_df$group <- rep(c("a", "b"), times = 900)
mv_df <- mv_df[mv_df$value > 0, ]
mv_df$label <- as.numeric(gsub("sim_", "", mv_df$variable))


test_that("Test brms mv trait non-longitudinal model skew model", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("mnormt")
  skip_on_cran()
  ss1 <- mvSS(
    model = "linear", form = label | value ~ group, df = mv_df,
    start = list("A" = 5), type = "brms", spectral_index = "ci_rededge"
  )
  expect_equal(ss1$family, "skew_normal")
  mod1 <- fitGrowth(ss1, backend = "cmdstanr", iter = 1000, chains = 1, cores = 1)
  expect_s3_class(mod1, "brmsfit")
  p <- growthPlot(mod1, ss1$pcvrForm, df = ss1$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test brms mv trait non-longitudinal model", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  ss1 <- mvSS(
    model = "linear", form = label | value ~ group, df = mv_df,
    start = list("A" = 5), type = "brms", spectral_index = "none"
  )
  expect_equal(ss1$family, "student")
  mod1 <- fitGrowth(ss1, backend = "cmdstanr", iter = 1000, chains = 1, cores = 1)
  expect_s3_class(mod1, "brmsfit")
  p <- growthPlot(mod1, ss1$pcvrForm, df = ss1$df)
  expect_s3_class(p, "ggplot")
  p2 <- growthPlot(mod1, ss1$pcvrForm, df = ss1$df, groups = "a")
  expect_s3_class(p2, "ggplot")
})

test_that("Test nls mv trait non-longitudinal model", {
  skip_on_cran()
  ss1 <- mvSS(
    model = "linear", form = label | value ~ group, df = mv_df,
    start = list("A" = 5), type = "nls", spectral_index = "ci_rededge"
  )
  mod1 <- fitGrowth(ss1)
  expect_s3_class(mod1, "lm")
  p <- growthPlot(mod1, ss1$pcvrForm, df = ss1$df, groupFill = TRUE)
  expect_s3_class(p, "ggplot")
  p2 <- growthPlot(mod1, ss1$pcvrForm, df = ss1$df, groups = "groupa")
  expect_s3_class(p, "ggplot")
})

test_that("Test nlrq mv trait non-longitudinal model", {
  skip_on_cran()
  ss1 <- mvSS(
    model = "linear", form = label | value ~ group, df = mv_df, tau = 0.5,
    start = list("A" = 5), type = "nlrq", spectral_index = "none"
  )
  mod1 <- fitGrowth(ss1)
  expect_s3_class(mod1, "rq")
  p <- growthPlot(mod1, ss1$pcvrForm, df = ss1$df)
  expect_s3_class(p, "ggplot")

  ss2 <- mvSS(
    model = "linear", form = label | value ~ group, df = mv_df, tau = seq(0.3, 0.7, 0.1),
    start = list("A" = 5), type = "nlrq", spectral_index = "none"
  )
  suppressWarnings(mod2 <- fitGrowth(ss2))
  expect_s3_class(mod2[[1]], "rq")
  p2 <- growthPlot(mod2, ss2$pcvrForm, df = ss2$df, groups = "groupa")
  expect_s3_class(p2, "ggplot")
})

#* `Longitudinal Multi-Value Trait Models`

m1 <- mvSim(
  dists = list(
    rnorm = list(mean = 100, sd = 30),
    rnorm = list(mean = 110, sd = 25),
    rnorm = list(mean = 120, sd = 20),
    rnorm = list(mean = 135, sd = 15)
  ),
  wide = FALSE, n = 6
)
m1$time <- rep(1:4, times = 6 * 180)
m2 <- mvSim(
  dists = list(
    rnorm = list(mean = 85, sd = 25),
    rnorm = list(mean = 95, sd = 20),
    rnorm = list(mean = 105, sd = 15),
    rnorm = list(mean = 110, sd = 15)
  ),
  wide = FALSE, n = 6
)
m2$time <- rep(1:4, times = 6 * 180)
mv_df2 <- rbind(m1, m2)
mv_df2$group <- rep(c("a", "b"), each = 4320)
mv_df2 <- mv_df2[mv_df2$value > 0, ]
mv_df2$label <- as.numeric(gsub("sim_", "", mv_df2$variable))

test_that("Test brms mv trait longitudinal model", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  ss_mv1 <- mvSS(
    model = "linear", form = label | value ~ time | group, df = mv_df2,
    start = list("A" = 50), type = "brms", spectral_index = "none"
  )
  fit <- fitGrowth(ss_mv1, backend = "cmdstanr", iter = 600, chains = 1, cores = 1)
  expect_s3_class(fit, "brmsfit")
  p <- growthPlot(fit, ss_mv1$pcvrForm, df = ss_mv1$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test nls mv trait longitudinal model", {
  skip_on_cran()
  ss_mv1 <- mvSS(
    model = "linear", form = label | value ~ time | group, df = mv_df2,
    type = "nls", spectral_index = "ci_rededge"
  )
  fit <- fitGrowth(ss_mv1)
  expect_s3_class(fit, "nls")
  p <- growthPlot(fit, ss_mv1$pcvrForm, df = ss_mv1$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test nlrq mv trait longitudinal model", {
  skip_on_cran()
  ss_mv1 <- mvSS(
    model = "linear", form = label | value ~ time | group, df = mv_df2,
    type = "nlrq", spectral_index = "ci_rededge"
  )
  fit <- fitGrowth(ss_mv1)
  expect_s3_class(fit, "nlrq")
  p <- growthPlot(fit, ss_mv1$pcvrForm, df = ss_mv1$df)
  expect_s3_class(p, "ggplot")

  ss_mv2 <- mvSS(
    model = "linear", form = label | value ~ time | group, df = mv_df2,
    type = "nlrq", spectral_index = "ci_rededge", tau = c(0.4, 0.5, 0.6)
  )
  fit2 <- fitGrowth(ss_mv2)
  expect_s3_class(fit2[[1]], "nlrq")
  p2 <- growthPlot(fit2, ss_mv2$pcvrForm, df = ss_mv2$df)
  expect_s3_class(p2, "ggplot")
})
