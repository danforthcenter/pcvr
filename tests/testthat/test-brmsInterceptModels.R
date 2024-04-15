if (file.exists("/home/josh/Desktop/") && interactive()) {
  library(testthat)
  library(brms)
  library(ggplot2)
  library(pcvr)

  test_that("Logistic brms Model with Intercept", {
    set.seed(123)
    simdf <- growthSim("logistic",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
    )
    simdf$y <- simdf$y + 15
    ss <- growthSS(
      model = "int_logistic", form = y ~ time | id / group, sigma = "int",
      list("A" = 130, "B" = 10, "C" = 3, "I" = 10),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "", "A", "B", "C", "I"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/int_logistic_fitGrowth.png",
      plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })


  test_that("Gompertz brms model pipeline", {
    set.seed(123)
    simdf <- growthSim("gompertz",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(0.25, 0.25))
    )
    simdf$y <- simdf$y + 30
    ss <- growthSS(
      model = "int_gompertz", form = y ~ time | id / group, sigma = "int",
      list("A" = 130, "B" = 10, "C" = 1, "I" = 20),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "I"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/gompertz_fitGrowth.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })

  test_that("intercept in submodel works", {
    set.seed(123)
    simdf <- growthSim("gompertz",
      n = 20, t = 25,
      params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(0.25, 0.25))
    )
    ss <- growthSS(
      model = "gompertz", form = y ~ time | id / group, sigma = "int_linear",
      list("A" = 130, "B" = 10, "C" = 1, "sigmaI" = 1, "sigmaA" = 2),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "sigmaI", "sigmaA"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/gompertz_interceptLinearSigma_fitGrowth.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })


  test_that("intercepts work in a changepoint model", {
    set.seed(123)
    simdf <- growthSim("linear + linear",
      n = 20, t = 25,
      params = list(
        "linear1A" = c(15, 12),
        "changePoint1" = c(8, 6),
        "linear2A" = c(3, 5)
      )
    )
    simdf$y <- simdf$y + 30
    ss <- growthSS(
      model = "int_linear + linear", form = y ~ time | id / group, sigma = "spline",
      list("linear1A" = 10, "changePoint1" = 5, "linear2A" = 2, "I" = 20),
      df = simdf, type = "brms"
    )
    expect_equal(ss$prior$nlpar, c("", "", "linear1A", "changePoint1", "linear2A", "I"))

    fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
    expect_s3_class(fit, "brmsfit")

    plot <- growthPlot(fit = fit, form = ss$pcvrForm, df = ss$df)
    ggsave("~/scripts/fahlgren_lab/labMeetings/linearPlusLinear_fitGrowth.png", plot,
      width = 10, height = 6, dpi = 300, bg = "#ffffff"
    )
    expect_s3_class(plot, "ggplot")
  })
}
