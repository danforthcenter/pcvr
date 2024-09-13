if (!interactive()) {
  pdf(NULL)
}
library(testthat)
library(pcvr)

test_that(".brmsMakeSurvPriors in survSS returns priors correctly", {
  set.seed(123)
  model <- "survival weibull"
  form <- y > 100 ~ time | id / group
  df <- growthSim("logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  surv <- .survModelParser(model)
  survivalBool <- surv$survival
  model <- surv$model
  prior <- brms::set_prior("normal(0, 10)", class = "b")
  ss <- .brmsSurvSS(model, form, df, prior)
  expect_s3_class(ss$prior, "brmsprior")
})

test_that("survival data is returned correctly if provided in survial format", {
  set.seed(123)
  model <- "survival weibull"
  form <- y > 100 ~ time | id / group
  df <- growthSim("logistic",
    n = 20, t = 25,
    params = list("A" = c(200), "B" = c(13), "C" = c(3))
  )
  df$event <- 1
  surv <- .survModelParser(model)
  survivalBool <- surv$survival
  model <- surv$model
  ss <- .survSS(model, form, df)
  expect_true(identical(ss$df, df))
})

test_that(".flexSurvSS returns errors and is covered", {
  set.seed(123)
  model <- "survival gengamma"
  form <- y > 100 ~ time | id / group
  df <- growthSim(
    "logistic",
    n = 20, t = 25,
    params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
  )
  surv <- .survModelParser(model)
  survivalBool <- surv$survival
  model <- surv$model
  expect_error(.flexSurvSS("survival normalorsomething", form, df))
  ss <- .flexSurvSS(model, form, df)
  df2 <- growthSim(
    "logistic",
    n = 20, t = 25,
    params = list("A" = c(200), "B" = c(13), "C" = c(3))
  )
  df2$event <- 1
  ss2 <- .flexSurvSS(model, form, df2, anc = x ~ y)
  expect_equal(names(ss2$formula), c("f1", "f2"))
  expect_true(identical(ss2$df, df2))
})
