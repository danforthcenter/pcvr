set.seed(123)
df <- growthSim("logistic",
  n = 20, t = 25,
  params = list(
    "A" = c(200, 190, 175, 160),
    "B" = c(13, 11, 12, 14),
    "C" = c(3, 3.25, 2.8, 3.1)
  )
)
df$group1 <- ifelse(df$group %in% c("a", "b"), "a", "b")
df$group2 <- ifelse(df$group %in% c("a", "c"), "c", "d")
df$group3 <- sample(c("e", "f"), nrow(df), replace = TRUE)

test_that("Logistic brms model with multiple groups", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  options(cmdstanr_warn_inits = FALSE)
  set.seed(123)
  ss <- growthSS(
    model = "logistic", form = y ~ time | id / group1 + group2, sigma = "spline",
    list("A" = 130, "B" = 10, "C" = 3),
    df = df, type = "brms"
  )
  expect_s3_class(ss, "pcvrss")
  fit <- fitGrowth(ss,
    backend = "cmdstanr", iter = 500, chains = 1, cores = 1,
    refresh = 0, silent = 2
  )
  expect_s3_class(fit, "brmsfit")
  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, groups = c("a", "d"))
  expect_s3_class(p2, "ggplot")
  p3 <- brmViolin(fit, ss, hypothesis = ".../A_group1a > 1.05")
  expect_s3_class(p3, "ggplot")
  p4 <- brmViolin(fit, ss, "A_group1a/A_group1b > 1.05")
  expect_s3_class(p4, "ggplot")
})

test_that("Logistic nls model with multiple groups", {
  set.seed(123)
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | id / group1 + group2,
    df = df, type = "nls"
  ))
  expect_s3_class(ss, "pcvrss")
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")
  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, groups = "a.d")
  expect_s3_class(p2, "ggplot")
})

test_that("Gam nls model with multiple groups", {
  set.seed(123)
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group1 + group2,
    df = df, type = "nls"
  ))
  expect_s3_class(ss, "pcvrss")
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "lm")
  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, groups = "a.d")
  expect_s3_class(p2, "ggplot")
})

test_that("Logistic nlrq model with multiple groups", {
  skip_if_not_installed("quantreg")
  set.seed(123)
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | id / group1 + group2, tau = 0.5,
    df = df, type = "nlrq"
  ))
  expect_s3_class(ss, "pcvrss")
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nlrq")
  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, groups = c("b", "c"))
  expect_s3_class(p2, "ggplot")
})

test_that("Gam nlrq model with multiple groups", {
  skip_if_not_installed("quantreg")
  set.seed(123)
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group1 + group2, tau = 0.5,
    df = df, type = "nlrq"
  ))
  expect_s3_class(ss, "pcvrss")
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "rq")
  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, groups = "a.d")
  expect_s3_class(p2, "ggplot")
})

test_that("Logistic nlme model with multiple groups", {
  skip_if_not_installed("nlme")
  skip_on_cran()
  set.seed(123)
  ss <- growthSS(
    model = "logistic", form = y ~ time | id / group1 + group2, sigma = "exp",
    df = df, type = "nlme"
  )
  expect_s3_class(ss, "pcvrss")
  if (file.exists("/home/josh/Desktop")) {
    fit <- suppressWarnings(fitGrowth(ss))
    expect_s3_class(fit, "nlme")
    p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    expect_s3_class(p, "ggplot")
    p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, groups = c("a", "d"))
    expect_s3_class(p2, "ggplot")
  }
})

test_that("Gam nlme model with multiple groups", {
  skip_if_not_installed("nlme")
  skip_on_cran()
  set.seed(123)
  ss <- growthSS(
    model = "gam", form = y ~ time | id / group1 + group2, sigma = "exp",
    df = df, type = "nlme"
  )
  expect_s3_class(ss, "pcvrss")
  if (file.exists("/home/josh/Desktop")) {
    fit <- suppressWarnings(fitGrowth(ss))
    expect_s3_class(fit, "lme")
    p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    expect_s3_class(p, "ggplot")
    p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, groups = "a.d")
    expect_s3_class(p2, "ggplot")
  }
})

test_that("Logistic mgcv model with multiple groups", {
  skip_if_not_installed("mgcv")
  skip_on_cran()
  set.seed(123)
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group1 + group2,
    df = df, type = "mgcv"
  ))
  expect_s3_class(ss, "pcvrss")
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "gam")
  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, groups = "a.d")
  expect_s3_class(p2, "ggplot")
})
