if (!interactive()) pdf(NULL)

#* ************************************************************
#* *************** `growthSim options` ***************
#* ************************************************************

test_that("Count data can be made by growthSim", {
  # use lowercase parameters
  df <- suppressWarnings(
    growthSim("count:gompertz",
      n = 20, t = 25,
      params = list("a" = c(100, 90), "b" = 10, "c" = 0.25)
    )
  )
  expect_equal(any(is.na(as.integer(df$y))), FALSE)
  # use unnamed parameters
  df <- growthSim("count:gompertz",
    n = 20, t = 25,
    params = list(100, 10, 0.25)
  )
  expect_equal(any(is.na(as.integer(df$y))), FALSE)
})

test_that("Fixed Changepoint data can be made by growthSim", {
  expect_error(
    growthSim(
      model = "gompertz + linear", n = 20, t = 50,
      params = list(100, 10, 0.25, 25, 3)
    )
  )
  df <- growthSim(
    model = "gompertz + linear", n = 20, t = 50,
    params = list(
      "gompertz1A" = 100, "gompertz1B" = 10, "gompertz1C" = 0.25,
      "fixedChangePoint1" = 25, "linear2A" = 3
    )
  )
  expect_equal(any(is.na(as.integer(df$y))), FALSE)
})

#* ************************************************************
#* *************** `Logistic growth modeling` ***************
#* ************************************************************

set.seed(123)
logistic_df <- growthSim("logistic",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
)

test_that("Test Logistic nls modeling", {
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | id / group,
    df = logistic_df, type = "nls"
  ))
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(
      184.201800401145, 184.201800401145, 12.0514357556166, 12.0514357556166,
      3.34892124655795, 3.34892124655795
    )
  )

  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")

  nls_p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
    ggplot2::labs(title = "nls")
  expect_s3_class(nls_p, "ggplot")

  test_res <- testGrowth(ss, fit, test = "A")$anova
  expect_s3_class(test_res, "anova")

  test_res <- testGrowth(ss, fit, test = list("A1 - A2 *1.1", "(B1+1) - B2", "C1 - (C2-0.5)"))
  expect_equal(dim(test_res), c(3, 5))
})

test_that("Test Logistic nlrq modeling", {
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | id / group,
    df = logistic_df, type = "nlrq", tau = c(0.5, 0.8)
  ))
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(
      184.201800401145, 184.201800401145, 12.0514357556166, 12.0514357556166,
      3.34892124655795, 3.34892124655795
    )
  )

  fit <- fitGrowth(ss)
  expect_s3_class(fit[[1]], "nlrq")

  nlrq_p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, groupFill = "plasma") +
    ggplot2::labs(title = "nlrq")
  expect_s3_class(nlrq_p, "ggplot")

  test_res <- suppressWarnings(testGrowth(ss, fit, test = "A")$`0.5`)
  expect_s3_class(test_res, "anova")

  test_res <- testGrowth(ss, fit = fit, test = "a|0.5|A > b|0.5|A")
  expect_equal(dim(test_res), c(2, 7))
})

test_that("Test Logistic nlme modeling", {
  ss <- growthSS(
    model = "logistic", form = y ~ time | id / group, sigma = "power",
    df = logistic_df, type = "nlme"
  )
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(
      184.201800401145, 184.201800401145, 12.0514357556166, 12.0514357556166,
      3.34892124655795, 3.34892124655795
    )
  )

  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "nlme")

  nlme_p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  nlme_p <- nlme_p +
    ggplot2::labs(title = "nlme")
  expect_s3_class(nlme_p, "ggplot")

  test_res <- suppressWarnings(testGrowth(ss, fit, test = "A")$anova)
  expect_s3_class(test_res, "anova.lme")

  test_res <- testGrowth(fit = fit, test = list(
    "A.groupa - A.groupb *1.1",
    "(B.groupa+1) - B.groupb",
    "C.groupa - (C.groupb-0.5)"
  ))
  expect_equal(dim(test_res), c(3, 5))
})

test_that("Test Logistic brms model setup", {
  ss <- growthSS(
    model = "logistic", form = y ~ time | id / group, sigma = "gompertz",
    list("A" = 130, "B" = 12, "C" = 3, "sigmaA" = 20, "sigmaB" = 15, "sigmaC" = 0.25),
    df = logistic_df, type = "brms"
  )
  expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "sigmaA", "sigmaB", "sigmaC"))

  expect_s3_class(ss$formula, "brmsformula")
})

#* ************************************************************
#* *************** `Testing pcvrFormula options` ***************
#* ************************************************************

test_that("Test Logistic nls modeling without individuals", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | group,
    df = logistic_df, type = "nls"
  ))
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(
      184.201800401145, 184.201800401145, 12.0514357556166, 12.0514357556166,
      3.34892124655795, 3.34892124655795
    )
  )

  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")

  nls_p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(nls_p2, "ggplot")
})

test_that("Test Logistic nls modeling without individuals or groups", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time,
    df = logistic_df, type = "nls"
  ))
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(184.201800401145, 12.0514357556166, 3.34892124655795)
  )

  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")

  nls_p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(nls_p2, "ggplot")
})

test_that("Test Logistic nlrq modeling without individuals", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time | group,
    df = logistic_df, type = "nlrq", tau = 0.5
  ))
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(
      184.201800401145, 184.201800401145, 12.0514357556166, 12.0514357556166,
      3.34892124655795, 3.34892124655795
    )
  )

  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nlrq")

  nlrq_p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(nlrq_p2, "ggplot")
})

test_that("Test Logistic nlrq modeling without individuals or groups", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "logistic", form = y ~ time,
    df = logistic_df, type = "nlrq", tau = 0.5
  ))
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(184.201800401145, 12.0514357556166, 3.34892124655795)
  )

  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nlrq")

  nlrq_p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(nlrq_p2, "ggplot")
})

test_that("Test Logistic nlme modeling without individuals", {
  skip_on_cran()
  ss <- growthSS(
    model = "logistic", form = y ~ time | group, sigma = "power", # failing on this so far
    df = logistic_df, type = "nlme"
  )
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(
      184.201800401145, 184.201800401145, 12.0514357556166, 12.0514357556166,
      3.34892124655795, 3.34892124655795
    )
  )

  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "nlme")

  nlme_p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(nlme_p2, "ggplot")
})

test_that("Test Logistic nlme modeling without individuals or groups", {
  skip_on_cran()
  ss <- growthSS(
    model = "logistic", form = y ~ time, sigma = "power", # failing on this so far
    df = logistic_df, type = "nlme"
  )
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(184.201800401145, 12.0514357556166, 3.34892124655795)
  )

  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "nlme")

  nlme_p2 <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(nlme_p2, "ggplot")
})


#* ************************************************************
#* *************** `Monomolecular growth modeling` ***************
#* ************************************************************

set.seed(123)
mono_df <- growthSim("monomolecular",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(0.08, 0.1))
)

test_that("Test monomolecular nls modeling", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "monomolecular", form = y ~ time | id / group,
    df = mono_df, type = "nls"
  ))
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(58.0855457770439, 58.0855457770439, 0.0457652584520121, 0.0457652584520121)
  )

  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")

  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test monomolecular nlrq modeling", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "monomolecular", form = y ~ time | id / group,
    df = mono_df, type = "nlrq"
  ))
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(58.0855457770439, 58.0855457770439, 0.0457652584520121, 0.0457652584520121)
  )

  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nlrq")

  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test monomolecular nlme modeling", {
  skip_on_cran()
  ss <- growthSS(
    model = "monomolecular", form = y ~ time | id / group, sigma = "power",
    df = mono_df, type = "nlme"
  )
  expect_equal(
    as.numeric(unlist(ss$start)),
    c(58.0855457770439, 58.0855457770439, 0.0457652584520121, 0.0457652584520121)
  )

  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "nlme")

  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test monomolecular brms model setup", {
  skip_on_cran()
  ss <- growthSS(
    model = "monomolecular", form = y ~ time | id / group, sigma = "spline",
    list("A" = 130, "B" = 0.1),
    df = mono_df, type = "brms"
  )
  expect_equal(ss$prior$nlpar, c("", "", "A", "B"))

  expect_s3_class(ss$formula, "brmsformula")
})

#* ************************************************************
#* *************** `Logarithmic growth modeling` ***************
#* ************************************************************

set.seed(123)
lgrthmc_df <- growthSim("logarithmic",
  n = 20, t = 25,
  params = list("A" = c(5, 7))
)

test_that("Test logarithmic nls modeling", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "logarithmic", form = y ~ time | id / group,
    df = lgrthmc_df, type = "nls"
  ))
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")
  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test logarithmic nlrq modeling", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "logarithmic", form = y ~ time | id / group,
    df = lgrthmc_df, type = "nlrq"
  ))
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nlrq")
  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test logarithmic nlme modeling", {
  skip_on_cran()
  ss <- growthSS(
    model = "logarithmic", form = y ~ time | id / group, sigma = "exp",
    df = lgrthmc_df, type = "nlme"
  )
  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "nlme")
  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test logarithmic brms model setup", {
  skip_on_cran()
  ss <- growthSS(
    model = "logarithmic", form = y ~ time | id / group, sigma = "spline",
    list("A" = 3),
    df = lgrthmc_df, type = "brms"
  )
  expect_equal(ss$prior$nlpar, c("", "", "A"))
  expect_s3_class(ss$formula, "brmsformula")
})

#* ************************************************************
#* *************** `general additive growth modeling` ***************
#* ************************************************************

set.seed(123)
gomp_df <- growthSim("gompertz",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(0.2, 0.25))
)


test_that("Test nls gam modeling", {
  skip_on_cran()
  ss1 <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group,
    df = gomp_df, type = "nls"
  ))
  expect_equal(as.character(ss1$formula), as.character(y ~ bs(time) * group))

  fit <- fitGrowth(ss1)
  expect_s3_class(fit, "lm")

  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")

  av <- testGrowth(ss = ss, fit)$anova
  expect_s3_class(av, "anova")
})

test_that("Test nlrq gam modeling", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group,
    df = gomp_df, type = "nlrq"
  ))
  expect_equal(as.character(ss$formula), as.character(y ~ bs(time) * group))

  fit <- fitGrowth(ss)
  expect_s3_class(fit, "rq")

  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  av <- suppressWarnings(testGrowth(ss = ss, fit)$anova)
  expect_s3_class(av, "anova.rq")
})

test_that("Test nlme gam", {
  skip_on_cran()
  ss <- growthSS(
    model = "gam", form = y ~ time | id / group, sigma = "exp",
    df = gomp_df, type = "nlme"
  )
  expect_equal(as.character(ss$formula$model), as.character(y ~ time * group))

  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "lme")

  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  av <- testGrowth(ss = ss, fit)$anova
  expect_s3_class(av, "anova.lme")
})

test_that("Test mgcv gam", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "gam", form = y ~ time | id / group,
    df = gomp_df, type = "mgcv"
  ))
  expect_equal(as.character(ss$formula), as.character(y ~ 0 + group + s(time, by = group)))

  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "gam")

  p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  p2 <- gam_diff(
    model = fit, g1 = "a", g2 = "b", plot = TRUE
  )
  expect_s3_class(p2$plot, "ggplot")

  av <- testGrowth(ss = ss, fit)$anova
  expect_s3_class(av, "anova")
})

test_that("Test gam brms model setup", {
  skip_on_cran()
  ss <- growthSS(
    model = "gam", form = y ~ time | id / group, sigma = "homo",
    df = gomp_df, type = "brms"
  )
  expect_s3_class(ss$formula, "brmsformula")
})

#* ************************************************************
#* ******************** `decay modeling` ********************
#* ************************************************************


test_that("Test logistic decay", {
  skip_on_cran()
  df <- simdf <- growthSim("logistic decay",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  ss <- growthSS(
    model = "logistic decay", form = y ~ time | id / group, sigma = "none",
    df = simdf, start = NULL, type = "nlme"
  )
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "lme")

  ss <- growthSS(
    model = "logistic decay", form = y ~ time | id / group, sigma = "none",
    df = simdf, start = NULL, type = "nls"
  )
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")

  ss <- growthSS(
    model = "logistic decay", form = y ~ time | id / group, sigma = "none",
    df = simdf, start = NULL, type = "nlrq"
  )
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nlrq")
})

#* ************************************************************
#* ******************** `time-to-event modeling` ********************
#* ************************************************************

test_that("Test survreg", {
  skip_on_cran()
  model <- "survival weibull"
  form <- y > 100 ~ time | id / group
  df <- growthSim("logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  ss <- growthSS(model = model, form = form, df = df, type = "survreg")
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "survreg")
  p <- growthPlot(fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
  test <- testGrowth(ss, fit)
  expect_s3_class(test, "survdiff")
})


#* ************************************************************
#* *************** `Models with Intercepts` ***************
#* ************************************************************

set.seed(123)
logistic_df <- growthSim("logistic",
  n = 20, t = 25,
  params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
)
logistic_df$y <- logistic_df$y + 20

test_that("Test Intercept Logistic nls modeling", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "int_logistic", form = y ~ time | id / group,
    df = logistic_df, type = "nls"
  ))
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")

  nls_p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
    ggplot2::labs(title = "nls")
  expect_s3_class(nls_p, "ggplot")

  test_res <- testGrowth(ss, fit, test = "A")$anova
  expect_s3_class(test_res, "anova")

  test_res <- testGrowth(ss, fit, test = list("A1 - A2 *1.1", "(B1+1) - B2", "C1 - (C2-0.5)"))
  expect_equal(dim(test_res), c(3, 5))
})

test_that("Test Intercept Logistic nlrq modeling", {
  skip_on_cran()
  ss <- suppressMessages(growthSS(
    model = "int_logistic", form = y ~ time | id / group,
    df = logistic_df, type = "nlrq", tau = 0.5
  ))

  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nlrq")

  nlrq_p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df) +
    ggplot2::labs(title = "nlrq")
  expect_s3_class(nlrq_p, "ggplot")
})

test_that("Test Intercept Logistic nlme modeling", {
  skip_on_cran()
  ss <- growthSS(
    model = "int_logistic", form = y ~ time | id / group, sigma = "power",
    df = logistic_df, type = "nlme"
  )

  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "nlme")

  nlme_p <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  nlme_p <- nlme_p +
    ggplot2::labs(title = "nlme")
  expect_s3_class(nlme_p, "ggplot")

  test_res <- suppressWarnings(testGrowth(ss, fit, test = "A")$anova)
  expect_s3_class(test_res, "anova.lme")

  test_res <- testGrowth(fit = fit, test = list(
    "A.groupa - A.groupb *1.1",
    "(B.groupa+1) - B.groupb",
    "C.groupa - (C.groupb-0.5)"
  ))
  expect_equal(dim(test_res), c(3, 5))
})

test_that("Test Intercept Monomolecular nls modeling", {
  skip_on_cran()
  set.seed(123)
  simdf <- growthSim(
    "monomolecular",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(0.08, 0.1))
  )
  simdf$y <- simdf$y + ifelse(simdf$group == "a", 10, 15)
  ss <- growthSS(
    model = "int_monomolecular", form = y ~ time | id / group,
    df = simdf, start = NULL, type = "nls"
  )
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")
})

test_that("Test Intercept linear nls modeling", {
  skip_on_cran()
  set.seed(123)
  simdf <- growthSim(
    "linear",
    n = 20, t = 25,
    params = list("A" = c(3, 4))
  )
  simdf$y <- simdf$y + ifelse(simdf$group == "a", 10, 15)
  ss <- growthSS(
    model = "int_linear", form = y ~ time | id / group,
    df = simdf, start = NULL, type = "nls"
  )
  fit <- fitGrowth(ss)
  coef(fit)
  expect_s3_class(fit, "nls")
})

#* ************************************************************
#* *************** `Dose-Response modeling` ***************
#* ************************************************************

test_that("Test Bragg in nls", {
  skip_on_cran()
  set.seed(123)
  simdf <- growthSim(
    "bragg",
    n = 20, t = 100,
    list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
  )
  ss <- growthSS(
    model = "bragg", form = y ~ time | id / group,
    df = simdf, start = NULL, type = "nls"
  )
  fit <- fitGrowth(ss)
  coef(fit)
  expect_s3_class(fit, "nls")
})

test_that("Test Bragg specification (not fitting) in nlme", {
  skip_on_cran()
  set.seed(123)
  simdf <- growthSim(
    "bragg",
    n = 20, t = 100,
    list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
  )
  ss <- growthSS(
    model = "bragg", form = y ~ time | id / group,
    df = simdf, start = NULL, type = "nlme"
  )
  expect_type(ss$formula, "list")
})

test_that("Test lorentz in nls", {
  skip_on_cran()
  set.seed(123)
  simdf <- growthSim(
    "lorentz",
    n = 20, t = 100,
    list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
  )
  ss <- growthSS(
    model = "lorentz", form = y ~ time | id / group,
    df = simdf, start = NULL, type = "nls"
  )
  fit <- fitGrowth(ss)
  coef(fit)
  expect_s3_class(fit, "nls")
})

test_that("Test lorentz specification (not fitting) in nlme", {
  skip_on_cran()
  set.seed(123)
  simdf <- growthSim(
    "lorentz",
    n = 20, t = 100,
    list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
  )
  ss <- growthSS(
    model = "lorentz", form = y ~ time | id / group,
    df = simdf, start = NULL, type = "nlme"
  )
  expect_type(ss$formula, "list")
})
