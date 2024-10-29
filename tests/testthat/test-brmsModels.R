#* there are lots of options and until one obviously breaks I am not going to try to test all of them.

test_that("inv_logit function", {
  expect_equal(inv_logit(100), 1)
})

test_that("Logistic brms model pipeline", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  set.seed(123)
  simdf <- growthSim(
    "logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  ss <- growthSS(
    model = "logistic", form = y ~ time | id / group, sigma = "gam",
    list("A" = 130, "B" = 10, "C" = 3),
    df = simdf, type = "brms"
  )
  expect_equal(ss$prior$nlpar, c("", "", "A", "B", "C"))

  fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
  expect_s3_class(fit, "brmsfit")

  plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(plot, "ggplot")
  plot1.5 <- growthPlot(fit = fit, form = y ~ time | group, groups = "a", df = ss$df)
  expect_s3_class(plot1.5, "ggplot")
  plot2 <- brmViolin(fit, ss, ".../A_groupa > 1.05")
  expect_s3_class(plot2, "ggplot")
  plot2.5 <- brmViolin(fit, ss, "B_groupb/B_groupa > 1.05")
  expect_s3_class(plot2.5, "ggplot")
  ss2 <- growthSS(
    model = "gompertz", form = y ~ time | id / group, sigma = "logistic",
    list("A" = 130, "B" = 10, "C" = 1, "sigmaA" = 20, "sigmaB" = 10, "sigmaC" = 2),
    df = simdf, type = "brms"
  )
  fit2 <- fitGrowth(ss2, backend = "cmdstanr", iter = 500, chains = 1, cores = 1, sample_prior = "only")
  expect_message(
    cd <- combineDraws(fit, fit2)
  )
  expect_equal(dim(cd), c(250, 21))
  cd <- combineDraws(fit, fit)
  expect_equal(dim(cd), c(250, 16))
  fit_df <- as.data.frame(fit)
  fit_df <- fit_df[, grepl("^b_", colnames(fit_df))]
  cd <- combineDraws(fit, fit_df)
  expect_equal(dim(cd), c(250, 16))
  expect_message(
    cd <- combineDraws(fit, fit_df[1:100, ])
  )
  expect_equal(dim(cd), c(250, 16))
  fit3 <- fitGrowth(ss2, backend = "cmdstanr", iter = 400, chains = 1, cores = 1, sample_prior = "only")
  expect_message(
    cd <- combineDraws(fit, fit3)
  )
  expect_equal(dim(cd), c(250, 21))
  expect_error(combineDraws(fit, list()))
  fit2 <- fit1 <- fit
  fit1$data <- fit1$data[fit1$data$time < 10, ]
  plot3 <- distributionPlot(list(fit1, fit2),
    form = ss$pcvrForm, d = ss$df,
    priors = list(
      "A" = rlnorm(500, log(130), 0.25),
      "B" = rlnorm(500, log(12), 0.25),
      "C" = rlnorm(500, log(3), 0.25)
    )
  )
  expect_s3_class(plot3, "ggplot")
  plots4 <- distributionPlot(
    list(fit1, fit2),
    form = ss$pcvrForm, d = ss$df,
    priors = NULL,
    patch = FALSE
  )
  test <- testGrowth(ss, fit, "A_groupa > A_groupb")
  expect_s3_class(test, "brmshypothesis")
  ss <- growthSS(
    model = "logistic", form = y ~ time | id / group, sigma = "logistic",
    list(
      "A" = 130, "B" = 10, "C" = 3,
      "sigmaA" = 20, "sigmaB" = 10, "sigmaC" = 3
    ),
    df = simdf, type = "brms"
  )
  pp1 <- plotPrior(ss)
  expect_s3_class(pp1, "ggplot")
  ss2 <- ss
  ss2$prior <- data.frame()
  expect_error(
    err <- plotPrior(ss2)
  )
  pp2 <- plotPrior(
    priors = list("A" = c(100, 130), "B" = c(10, 8), "C" = c(0.2, 0.1)),
    type = "logistic",
    n = 200, t = 25
  )
  expect_s3_class(pp2$simulated, "ggplot")
  pp3 <- plotPrior(
    priors = list("A" = c(100, 130), "B" = c(0.08, 0.1)),
    type = "monomolecular",
    n = 200, t = 25
  )
  expect_s3_class(pp3$simulated, "ggplot")
  pp4 <- plotPrior(
    priors = list("A" = c(101, 11), "B" = c(0.12, 0.15)),
    type = "exponential",
    n = 200, t = 25
  )
  expect_s3_class(pp4$simulated, "ggplot")
  suppressWarnings(barg_output1 <- barg(fit, ss))
  fit_2 <- fit
  fit_list <- list(fit, fit_2)
  x <- suppressWarnings(barg(fit_list, list(ss, ss)))
})

test_that("distPlot works with many models", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  load(url("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/brmsFits.rdata"))
  fits <- list(fit_3, fit_15)
  form <- y ~ time | id / group
  priors <- list(
    "phi1" = rlnorm(2000, log(130), 0.25),
    "phi2" = rlnorm(2000, log(12), 0.25),
    "phi3" = rlnorm(2000, log(3), 0.25)
  )
  from3to25 <- list(
    fit_3, fit_5, fit_7, fit_9, fit_11,
    fit_13, fit_15, fit_17, fit_19, fit_21, fit_23, fit_25
  )
  plot <- distributionPlot(
    fits = from3to25, form = y ~ time | id / group,
    params = c("A", "B", "C"), d = simdf, priors = priors
  )
  expect_s3_class(plot, "ggplot")
})

test_that("brms model warns about priors", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  set.seed(123)
  simdf <- growthSim(
    "linear",
    n = 20, t = 25,
    params = list("A" = c(1))
  )
  ss <- growthSS(
    model = "linear", form = y ~ time, sigma = "spline",
    df = simdf, type = "brms"
  )
  ss <- ss[-which(names(ss) == "prior")]
  expect_warning(
    fit <- fitGrowth(
      ss,
      backend = "cmdstanr",
      iter = 100, chains = 1, cores = 1
    )
  )
  plot <- growthPlot(fit, form = ss$pcvrForm)
  expect_s3_class(plot, "ggplot")
})

test_that("Hierarchical Model Works", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  set.seed(123)
  simdf <- growthSim(
    "logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  simdf$covar <- rnorm(nrow(simdf), 10, 1)
  ss <- growthSS(
    model = "logistic", form = y ~ time + covar | id / group, sigma = "logistic",
    list(
      "AI" = 100, "AA" = 5, "B" = 10, "C" = 3,
      "sigmaA" = 10, "sigmaB" = 10, "sigmaC" = 3
    ),
    df = simdf, type = "brms",
    hierarchy = list("A" = "int_linear")
  )
  lapply(ss, head)
  fit <- fitGrowth(ss, iter = 600, cores = 1, chains = 1, backend = "cmdstanr")
  expect_s3_class(fit, "brmsfit")
  p <- growthPlot(fit, ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Changepoint model can be specified", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  set.seed(123)
  noise <- do.call(rbind, lapply(1:30, function(i) {
    chngpt <- c(20, 21)
    rbind(
      data.frame(
        id = paste0("id_", i), time = 1:chngpt[1], group = "a",
        y = c(runif(chngpt[1] - 1, 0, 20), rnorm(1, 5, 1))
      ),
      data.frame(
        id = paste0("id_", i), time = 1:chngpt[2], group = "b",
        y = c(runif(chngpt[2] - 1, 0, 20), rnorm(1, 5, 1))
      )
    )
  }))
  noise2 <- do.call(rbind, lapply(1:30, function(i) {
    start1 <- max(noise[noise$id == paste0("id_", i) & noise$group == "a", "time"])
    start2 <- max(noise[noise$id == paste0("id_", i) & noise$group == "b", "time"])
    rbind(
      data.frame(
        id = paste0("id_", i), time = start1:40, group = "a",
        y = c(runif(length(start1:40), 15, 50))
      ),
      data.frame(
        id = paste0("id_", i), time = start2:40, group = "b",
        y = c(runif(length(start2:40), 15, 50))
      )
    )
  }))
  simdf <- rbind(noise, noise2)
  ss <- growthSS(
    model = "int_linear + decay linear", form = y ~ time | id / group, sigma = "int + gam",
    list(
      "I" = 1, "linear1A" = 10, "fixedChangePoint1" = 20, "linear2A" = 2, "sigmaint1" = 1,
      "sigmachangePoint1" = 25
    ),
    df = simdf, type = "brms"
  )
  expect_equal(ss$prior$nlpar, c("", "linear1A", "linear2A", "I", "sigmaint1", "sigmachangePoint1"))
  ss <- growthSS(
    model = "int_linear + decay linear", form = y ~ time | id / group, sigma = "int + gam",
    list(
      "I" = 1, "linear1A" = 10, "fixedChangePoint1" = 20, "linear2A" = 2, "sigmaint1" = 1,
      "sigmachangePoint1" = 25
    ),
    df = simdf[1:10, ], type = "brms"
  )
  expect_equal(ss$prior$nlpar, c("", "linear1A", "linear2A", "I", "sigmaint1", "sigmachangePoint1"))
})

test_that("weibull survival", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  set.seed(123)
  model <- "survival"
  form <- y > 100 ~ time | id / group
  df <- growthSim(
    "logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  prior <- c(0, 5)
  ss <- growthSS(model = model, form = form, df = df, start = prior)
  expect_equal(ss$prior$coef, c("groupa", "groupb"))
  fit <- fitGrowth(ss, iter = 600, cores = 1, chains = 1, backend = "cmdstanr")
  expect_s3_class(fit, "brmsfit")
  plot <- growthPlot(fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(plot, "ggplot")
  test <- testGrowth(ss, fit, "groupa > groupb")
  expect_s3_class(test, "brmshypothesis")
})

test_that("binomial survival", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_on_cran()
  set.seed(123)
  model <- "survival binomial"
  form <- y > 100 ~ time | id / group
  df <- growthSim(
    "logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  prior <- c(0, 5)
  ss <- growthSS(model = model, form = form, df = df, start = prior)
  fit <- fitGrowth(ss, iter = 600, cores = 1, chains = 1, backend = "cmdstanr")
  expect_s3_class(fit, "brmsfit")
  plot <- growthPlot(fit, form = ss$pcvrForm, df = ss$df)
  expect_s3_class(plot, "ggplot")
})

test_that(".brmSurvSS options all work", {
  set.seed(123)
  df <- growthSim("logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  surv <- .survModelParser("survival weibull")
  ss <- suppressMessages(
    .brmsSurvSS(
      model = surv$model,
      form = y > 100 ~ time | id / group,
      df = df,
      priors = c(0, 5)
    )
  )
  expect_equal(names(ss), c("df", "family", "formula", "prior", "initfun", "pcvrForm"))
  ss2 <- suppressMessages(
    .brmsSurvSS(
      model = surv$model,
      form = y > 100 ~ time | id / group,
      df = df,
      priors = NULL
    )
  )
  expect_equal(names(ss2), c("df", "family", "formula", "prior", "initfun", "pcvrForm"))
  df2 <- df
  df2$censor <- 0 # dummy data
  df2$event <- 1 # dummy data
  df2$n_eligible <- 100 # dummy data
  df2$n_events <- 5 # dummy data
  ss3 <- suppressMessages(
    .brmsSurvSS(
      model = surv$model,
      form = y > 100 ~ time | id / group,
      df = df2,
      priors = NULL
    )
  )
  expect_equal(names(ss3), c("df", "family", "formula", "prior", "initfun", "pcvrForm"))
  ss4 <- suppressMessages(
    .brmsSurvSS(
      model = surv$model,
      form = y > 100 ~ time | id / group,
      df = df2[df2$group == "a", ],
      priors = NULL
    )
  )
  expect_equal(names(ss4), c("df", "family", "formula", "prior", "initfun", "pcvrForm"))

  surv <- .survModelParser("survival binomial")
  ss <- suppressMessages(
    .brmsSurvSS(
      model = surv$model,
      form = y > 100 ~ time | id / group,
      df = df,
      priors = c(0, 5)
    )
  )
  expect_equal(names(ss), c("df", "family", "formula", "prior", "initfun", "pcvrForm"))
  ss2 <- suppressMessages(
    .brmsSurvSS(
      model = surv$model,
      form = y > 100 ~ time | id / group,
      df = df,
      priors = NULL
    )
  )
  expect_equal(names(ss2), c("df", "family", "formula", "prior", "initfun", "pcvrForm"))
  ss3 <- suppressMessages(
    .brmsSurvSS(
      model = surv$model,
      form = y > 100 ~ time | id / group,
      df = df2,
      priors = NULL
    )
  )
  expect_equal(names(ss3), c("df", "family", "formula", "prior", "initfun", "pcvrForm"))
  ss4 <- suppressMessages(
    .brmsSurvSS(
      model = surv$model,
      form = y > 100 ~ time | id / group,
      df = df2[df2$group == "a", ],
      priors = NULL
    )
  )
  expect_equal(names(ss4), c("df", "family", "formula", "prior", "initfun", "pcvrForm"))
})

#* ***********************************
#* ***** `Not Run on the remote` *****
#* ***********************************

if (file.exists("/home/josh/Desktop/")) {
  # only run locally, don't test for each R-CMD Check
  test_that("Gompertz brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("gompertz",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(0.25, 0.25))
    )

    ss <- growthSS(
      model = "gompertz", form = y ~ time | id / group, sigma = "int",
      list("A" = 130, "B" = 10, "C" = 1),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "B", "C"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/gompertz_fitGrowth.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("Monomolecular brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("monomolecular",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(0.01, 0.08))
    )

    ss <- growthSS(
      model = "monomolecular", form = y ~ time | id / group, sigma = "int",
      list("A" = 130, "B" = 1),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "B"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/monomolecular_fitGrowth.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("Exponential brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("exponential",
      n = 20, t = 25,
      params = list("A" = c(15, 12), "B" = c(0.1, 0.085))
    )

    ss <- growthSS(
      model = "exponential", form = y ~ time | id / group, sigma = "int",
      list("A" = 10, "B" = 1),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "B"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/exponential_fitGrowth.png", plot,
      width = 10,
      height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("Power law brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("power law",
      n = 20, t = 25,
      params = list("A" = c(15, 12), "B" = c(0.75, 0.8))
    )

    ss <- growthSS(
      model = "power law", form = y ~ time | id / group, sigma = "linear",
      list("A" = 10, "B" = 1),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "B"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/powerlaw_fitGrowth.png", plot,
      width = 10,
      height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("linear brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("linear",
      n = 20, t = 25,
      params = list("A" = c(15, 12))
    )

    ss <- growthSS(
      model = "linear", form = y ~ time | id / group, sigma = "int",
      list("A" = 5),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/linear_fitGrowth.png", plot,
      width = 10,
      height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("logarithmic brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("logarithmic",
      n = 20, t = 25,
      params = list("A" = c(15, 12))
    )

    ss <- growthSS(
      model = "logarithmic", form = y ~ time | id / group, sigma = "int",
      list("A" = 5),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/logarithmic_fitGrowth.png", plot,
      width = 10,
      height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("linear sub model with prior brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("linear",
      n = 20, t = 25,
      params = list("A" = c(15, 12))
    )

    model <- "linear"
    form <- y ~ time | id / group
    sigma <- "linear"
    priors <- list("A" = 5, "sigmaA" = 2)
    df <- simdf
    type <- "brms"

    ss <- growthSS(
      model = "linear", form = y ~ time | id / group, sigma = "linear",
      list("A" = 5, "sigmaA" = 2),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "sigmaA"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/linear_fitGrowth.png", plot,
      width = 10,
      height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("GAM brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("logistic",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
    )

    ss <- growthSS(
      model = "gam", form = y ~ time | id / group, sigma = "int",
      df = simdf, type = "brms"
    )

    fit <- suppressWarnings(fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1))
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/gam_fitGrowth.png", plot,
      width = 10,
      height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })



  test_that("linear+linear brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("linear + linear",
      n = 20, t = 25,
      params = list("linear1A" = c(15), "changePoint1" = c(8), "linear2A" = c(3))
    )

    ss <- growthSS(
      model = "linear + linear", form = y ~ time, sigma = "spline",
      list("linear1A" = 10, "changePoint1" = 5, "linear2A" = 2),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "", "linear1A", "changePoint1", "linear2A"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/linearPlusLinear_fitGrowth.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("linear+logistic brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("linear + logistic",
      n = 20, t = 25,
      params = list(
        "linear1A" = c(15, 12), "changePoint1" = c(8, 6),
        "logistic2A" = c(100, 150), "logistic2B" = c(10, 8), "logistic2C" = c(3, 2.5)
      )
    )

    ss <- growthSS(
      model = "linear + logistic", form = y ~ time | id / group, sigma = "spline",
      list(
        "linear1A" = 10, "changePoint1" = 5,
        "logistic2A" = 100, "logistic2B" = 10, "logistic2C" = 3
      ),
      df = simdf, type = "brms"
    )

    expect_equal(ss$prior$nlpar, c(
      "", "", "linear1A", "changePoint1",
      "logistic2A", "logistic2B", "logistic2C"
    ))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/linearPlusLogistic_fitGrowth.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })


  test_that("linear+gam brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("linear + logistic",
      n = 20, t = 25, # using logistic data, but modeling as a gam
      params = list(
        "linear1A" = c(15, 12), "changePoint1" = c(8, 6),
        "logistic2A" = c(100, 150), "logistic2B" = c(10, 8), "logistic2C" = c(3, 2.5)
      )
    )

    ss <- growthSS(
      model = "linear + gam", form = y ~ time | id / group, sigma = "homo",
      list("linear1A" = 10, "changePoint1" = 5),
      df = simdf, type = "brms"
    )

    expect_equal(ss$prior$nlpar, c("", "linear1A", "changePoint1"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/linearPlusGAM_fitGrowth.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })


  test_that("linear + linear + linear brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("linear + linear + linear",
      n = 25, t = 50,
      params = list(
        "linear1A" = c(10, 12), "changePoint1" = c(8, 6), "linear2A" = c(1, 2),
        "changePoint2" = c(25, 30), "linear3A" = c(20, 24)
      )
    )

    ss <- growthSS(
      model = "linear + linear + linear", form = y ~ time | id / group, sigma = "spline",
      list("linear1A" = 10, "changePoint1" = 5, "linear2A" = 2, "changePoint2" = 15, "linear3A" = 5),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c(
      "", "", "linear1A", "changePoint1",
      "linear2A", "changePoint2", "linear3A"
    ))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/linearPlusLinearPlusLinear_fitGrowth.png",
      plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("Logistic brms logistic sub model pipeline", {
    set.seed(123)
    simdf <- growthSim("logistic",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
    )

    ss <- growthSS(
      model = "logistic", form = y ~ time | id / group, sigma = "logistic",
      list("A" = 130, "B" = 10, "C" = 3, "sigmaA" = 20, "sigmaB" = 10, "sigmaC" = 2),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "sigmaA", "sigmaB", "sigmaC"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1) # that's fast
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/logistic_logisticSubModel.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("Logistic brms gompertz sub model pipeline", {
    set.seed(123)
    simdf <- growthSim("logistic",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
    )

    ss <- growthSS(
      model = "logistic", form = y ~ time | id / group, sigma = "gompertz",
      list("A" = 130, "B" = 10, "C" = 3, "sigmaA" = 20, "sigmaB" = 10, "sigmaC" = 2),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "sigmaA", "sigmaB", "sigmaC"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/logistic_gompSubModel.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("Logistic brms monomolecular sub model pipeline", {
    set.seed(123)
    simdf <- growthSim("logistic",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
    )

    ss <- growthSS(
      model = "logistic", form = y ~ time | id / group, sigma = "monomolecular",
      list("A" = 130, "B" = 10, "C" = 3, "sigmaA" = 5, "sigmaB" = 0.5),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "sigmaA", "sigmaB"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/logistic_monoSubModel.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })


  test_that("int+int homoskedastic model pipeline", {
    set.seed(123)

    noise <- do.call(rbind, lapply(1:30, function(i) {
      chngpt <- rnorm(2, 18, 2)
      rbind(
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[1], group = "a",
          y = c(runif(chngpt[1] - 1, 0, 20), rnorm(1, 5, 1))
        ),
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[2], group = "b",
          y = c(runif(chngpt[2] - 1, 0, 20), rnorm(1, 5, 1))
        )
      )
    }))
    noise2 <- do.call(rbind, lapply(1:30, function(i) {
      start1 <- max(noise[noise$id == paste0("id_", i) & noise$group == "a", "time"])
      start2 <- max(noise[noise$id == paste0("id_", i) & noise$group == "b", "time"])

      rbind(
        data.frame(
          id = paste0("id_", i), time = start1:40, group = "a",
          y = c(runif(length(start1:40), 15, 50))
        ),
        data.frame(
          id = paste0("id_", i), time = start2:40, group = "b",
          y = c(runif(length(start2:40), 15, 50))
        )
      )
    }))
    simdf <- rbind(noise, noise2)

    ss <- growthSS(
      model = "int + int", form = y ~ time | id / group, sigma = "int",
      list("int1" = 10, "changePoint1" = 10, "int2" = 20),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "int1", "changePoint1", "int2"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/intPlusInt_fitGrowth.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("int+int fixed changepoint homoskedastic model pipeline", {
    set.seed(123)

    noise <- do.call(rbind, lapply(1:30, function(i) {
      chngpt <- c(20, 21)
      rbind(
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[1], group = "a",
          y = c(runif(chngpt[1] - 1, 0, 20), rnorm(1, 5, 1))
        ),
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[2], group = "b",
          y = c(runif(chngpt[2] - 1, 0, 20), rnorm(1, 5, 1))
        )
      )
    }))
    noise2 <- do.call(rbind, lapply(1:30, function(i) {
      start1 <- max(noise[noise$id == paste0("id_", i) & noise$group == "a", "time"])
      start2 <- max(noise[noise$id == paste0("id_", i) & noise$group == "b", "time"])

      rbind(
        data.frame(
          id = paste0("id_", i), time = start1:40, group = "a",
          y = c(runif(length(start1:40), 15, 50))
        ),
        data.frame(
          id = paste0("id_", i), time = start2:40, group = "b",
          y = c(runif(length(start2:40), 15, 50))
        )
      )
    }))
    simdf <- rbind(noise, noise2)

    ss <- growthSS(
      model = "int + int", form = y ~ time | id / group, sigma = "int",
      list("int1" = 10, "fixedChangePoint1" = 20, "int2" = 20),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "int1", "int2"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/intPlusInt_fixedChngpt_fitGrowth.png",
      plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })



  test_that("int+int thresholded homoskedasticity model pipeline", {
    set.seed(123)

    noise <- do.call(rbind, lapply(1:30, function(i) {
      chngpt <- rnorm(2, 18, 2)
      rbind(
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[1], group = "a",
          y = c(runif(chngpt[1] - 1, 0, 20), rnorm(1, 5, 1))
        ),
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[2], group = "b",
          y = c(runif(chngpt[2] - 1, 0, 20), rnorm(1, 5, 1))
        )
      )
    }))
    noise2 <- do.call(rbind, lapply(1:30, function(i) {
      start1 <- max(noise[noise$id == paste0("id_", i) & noise$group == "a", "time"])
      start2 <- max(noise[noise$id == paste0("id_", i) & noise$group == "b", "time"])

      rbind(
        data.frame(
          id = paste0("id_", i), time = start1:40, group = "a",
          y = c(runif(length(start1:40), 15, 50))
        ),
        data.frame(
          id = paste0("id_", i), time = start2:40, group = "b",
          y = c(runif(length(start2:40), 15, 50))
        )
      )
    }))
    simdf <- rbind(noise, noise2)

    ss <- growthSS(
      model = "int + int", form = y ~ time | id / group, sigma = "int + int",
      list(
        "int1" = 10, "changePoint1" = 10, "int2" = 20, "sigmaint1" = 10,
        "sigmachangePoint1" = 10, "sigmaint2" = 10
      ),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c(
      "", "int1", "changePoint1", "int2", "sigmaint1",
      "sigmachangePoint1", "sigmaint2"
    ))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/intPlusInt_heteroskedastic_fitGrowth.png",
      plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })




  test_that("int + linear model and submodel pipeline", {
    set.seed(123)
    noise <- do.call(rbind, lapply(1:30, function(i) {
      chngpt <- rnorm(2, 18, 2)
      rbind(
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[1], group = "a",
          y = c(runif(chngpt[1] - 1, 0, 20), rnorm(1, 5, 1))
        ),
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[2], group = "b",
          y = c(runif(chngpt[2] - 1, 0, 20), rnorm(1, 5, 1))
        )
      )
    }))
    signal <- growthSim("linear",
      n = 30, t = 20,
      params = list("A" = c(3, 5))
    )
    signal <- do.call(rbind, lapply(unique(paste0(signal$id, signal$group)), function(int) {
      noisesub <- noise[paste0(noise$id, noise$group) == int, ]
      signalSub <- signal[paste0(signal$id, signal$group) == int, ]
      y_end <- noisesub[noisesub$time == max(noisesub$time), "y"]
      signalSub$time <- signalSub$time + max(noisesub$time)
      signalSub$y <- y_end + signalSub$y
      signalSub
    }))
    simdf <- rbind(noise, signal)
    ss <- growthSS(
      model = "int + linear", form = y ~ time | id / group, sigma = "int + linear",
      list(
        "int1" = 10, "changePoint1" = 10, "linear2A" = 20, "sigmaint1" = 10,
        "sigmachangePoint1" = 10, "sigmalinear2A" = 10
      ),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c(
      "", "int1", "changePoint1", "linear2A", "sigmaint1",
      "sigmachangePoint1", "sigmalinear2A"
    ))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df, timeRange = 1:40)
    ggsave("~/scripts/fahlgren_lab/labMeetings/intPlusLinear_heteroskedIntPlusLinear_fitGrowth.png",
      plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })


  test_that("int + Logistic brms int+spline sub model pipeline", {
    set.seed(123)
    noise <- do.call(rbind, lapply(1:30, function(i) {
      chngpt <- rnorm(2, 18, 2)
      rbind(
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[1], group = "a",
          y = c(runif(chngpt[1] - 1, 0, 20), rnorm(1, 5, 1))
        ),
        data.frame(
          id = paste0("id_", i), time = 1:chngpt[2], group = "b",
          y = c(runif(chngpt[2] - 1, 0, 20), rnorm(1, 5, 1))
        )
      )
    }))
    signal <- growthSim("logistic",
      n = 20, t = 30,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
    )
    signal <- do.call(rbind, lapply(unique(paste0(signal$id, signal$group)), function(int) {
      noisesub <- noise[paste0(noise$id, noise$group) == int, ]
      signalSub <- signal[paste0(signal$id, signal$group) == int, ]
      y_end <- noisesub[noisesub$time == max(noisesub$time), "y"]
      signalSub$time <- signalSub$time + max(noisesub$time)
      signalSub$y <- y_end + signalSub$y
      signalSub
    }))
    simdf <- rbind(noise, signal)
    simdf <- simdf[simdf$time < 45, ]

    ss <- growthSS(
      model = "int+logistic", form = y ~ time | id / group, sigma = "int + spline",
      list(
        "int1" = 5, "changePoint1" = 10, "logistic2A" = 130, "logistic2B" = 10, "logistic2C" = 3,
        "sigmaint1" = 5, "sigmachangePoint1" = 15
      ),
      df = simdf, type = "brms"
    )


    expect_equal(ss$prior$nlpar, c(
      "", "int1", "changePoint1", "logistic2A", "logistic2B", "logistic2C",
      "sigmaint1", "sigmachangePoint1"
    ))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/intPluslogistic_intPlusGAMSubModel.png",
      plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })


  test_that("fixed and estimated changepoints can be mixed in growth formula", {
    simdf1 <- growthSim(
      model = "logistic", n = 20, t = 20,
      params = list("A" = c(180, 160), "B" = c(9, 11), "C" = c(3, 3.5))
    )

    simdf2 <- growthSim(
      model = "linear + linear", n = 20, t = 20,
      params = list("linear1A" = c(6, 8), "changePoint1" = c(7, 9), "linear2A" = c(15, 20))
    )

    simdf2_adj <- do.call(rbind, lapply(unique(paste0(simdf2$id, simdf2$group)), function(int) {
      p1 <- simdf1[paste0(simdf1$id, simdf1$group) == int, ]
      p2 <- simdf2[paste0(simdf2$id, simdf2$group) == int, ]
      y_end_p1 <- p1[p1$time == max(p1$time), "y"]
      p2$time <- p2$time + max(p1$time)
      p2$y <- y_end_p1 + p2$y
      p2
    }))
    simdf <- rbind(simdf1, simdf2_adj)

    ss <- growthSS(
      model = "logistic+linear+linear", form = y ~ time | id / group,
      sigma = "logistic+linear", df = simdf,
      start = list(
        "logistic1A" = 130, "logistic1B" = 10, "logistic1C" = 3.5,
        "fixedChangePoint1" = 20, "linear2A" = 5, "changePoint2" = 28, "linear3A" = 20,
        "sigmalogistic1A" = 10, "sigmalogistic1B" = 12, "sigmalogistic1C" = 20,
        "sigmafixedChangePoint1" = 20, "sigmalinear2A" = 3
      ), type = "brms"
    )

    expect_equal(ss$prior$nlpar, c(
      "", "logistic1A", "logistic1B", "logistic1C", "linear2A", "changePoint2",
      "linear3A", "sigmalogistic1A", "sigmalogistic1B", "sigmalogistic1C",
      "sigmalinear2A"
    ))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/threePart_fixedAndEstimatedChangepoint.png",
      plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("logistic decay as a segment", {
    simdf <- growthSim(
      model = "logistic + logistic decay", n = 20, t = 45,
      params = list(
        "logistic1A" = c(120, 140), "logistic1B" = c(12, 10), "logistic1C" = c(3, 3.5),
        "changePoint1" = c(20, 23),
        "logistic2A" = c(90, 100), "logistic2B" = c(11, 13), "logistic2C" = c(3, 3.5)
      )
    )
    ss <- growthSS(
      model = "logistic + logistic decay", form = y ~ time | id / group, sigma = "spline",
      list(
        "logistic1A" = 100, "logistic1B" = 10, "logistic1C" = 3, "changePoint1" = 20,
        "logistic2A" = 100, "logistic2B" = 10, "logistic2C" = 3
      ),
      df = simdf, type = "brms"
    )
    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/logistic_plus_logisticDecay.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("Test flexsurv model", {
    set.seed(123)
    model <- "survival gompertz"
    form <- y > 100 ~ time | id / group
    df <- growthSim("logistic",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
    )
    ss <- growthSS(model = model, form = form, df = df, type = "flexsurv")
    library(flexsurv)
    fit <- fitGrowth(ss)
    p <- growthPlot(fit, form = ss$pcvrForm, df = ss$df)
    expect_s3_class(p, "ggplot")
  })

  test_that("Logistic Poisson Model", {
    set.seed(123)
    form <- y ~ time | id / group
    df <- growthSim("count: logistic",
      n = 20, t = 25,
      params = list("A" = c(10, 12), "B" = c(13, 11), "C" = c(3, 3.5))
    )
    ss <- growthSS(
      model = "poisson: logistic", form = y ~ time | id / group, sigma = NULL,
      df = df, start = list("A" = 8, "B" = 10, "C" = 3)
    )
    lapply(ss, head)
    fit <- fitGrowth(ss, iter = 600, cores = 1, chains = 1, backend = "cmdstanr")
    expect_s3_class(fit, "brmsfit")
    p <- growthPlot(fit, ss$pcvrForm, df = ss$df)
    expect_s3_class(p, "ggplot")
  })

  test_that("Beta DRC Model", {
    set.seed(123)
    form <- y ~ time | id / group
    df <- growthSim(
      "beta",
      n = 20, t = 50,
      params = list(
        "A" = c(10, 10),
        "B" = c(1.25, 1.3),
        "C" = c(20, 22),
        "D" = c(5, 5),
        "E" = c(30, 32)
      )
    )
    #* consider using ss with nls to get ideas for parameters
    ss <- growthSS(
      model = "beta", form = y ~ time | id / group, sigma = NULL,
      df = df, start = list("A" = 10, "B" = 1, "C" = 15, "D" = 3, "E" = 25)
    )
    expect_s3_class(ss, "pcvrss")
    tryCatch({
      fit <- fitGrowth(ss, iter = 600, cores = 1, chains = 1, backend = "cmdstanr")
      expect_s3_class(fit, "brmsfit")
      p <- growthPlot(fit, ss$pcvrForm, df = ss$df)
      expect_s3_class(p, "ggplot")
    }, error = function(err) {})
  })
}
