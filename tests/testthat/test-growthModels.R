
if(!interactive()) pdf(NULL)

#* ************************************************************
#* *************** `Logistic growth modeling` *************** 
#* ************************************************************

set.seed(123)
logistic_df<-growthSim("logistic", n=20, t=25,
   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))

test_that("Test Logistic nls modeling", {
  ss<-suppressMessages(growthSS(model = "logistic", form=y~time|id/group, 
               df=logistic_df, type = "nls"))
  expect_equal(as.numeric(unlist(ss$start)),
               c(184.201800401145, 184.201800401145, 12.0514357556166, 12.0514357556166, 
                                   3.34892124655795, 3.34892124655795) )
  
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")
  
  nls_p <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)+
    ggplot2::labs(title="nls")
  expect_s3_class(nls_p, "ggplot")
})

test_that("Test Logistic nlrq modeling", {
  ss<-suppressMessages(growthSS(model = "logistic", form=y~time|id/group,
               df=logistic_df, type = "nlrq", tau = 0.5)) # seq(0.02, 0.99, 0.02)
  expect_equal(as.numeric(unlist(ss$start)),
               c(184.201800401145, 184.201800401145, 12.0514357556166, 12.0514357556166, 
                 3.34892124655795, 3.34892124655795) )
  
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nlrq")
  
  nlrq_p <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)+
    ggplot2::labs(title="nlrq")
  expect_s3_class(nlrq_p, "ggplot")
})

test_that("Test Logistic nlme modeling", {
  ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="power",
               df=logistic_df, type = "nlme")
  expect_equal(as.numeric(unlist(ss$start)),
               c(184.201800401145, 184.201800401145, 12.0514357556166, 12.0514357556166, 
                 3.34892124655795, 3.34892124655795) )
  
  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "nlme")
  
  nlme_p <- suppressWarnings(growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df, boot = 3))
  nlme_p<-nlme_p+
    ggplot2::labs(title="nlme")
  expect_s3_class(nlme_p, "ggplot")
})

test_that("Test Logistic brms model setup", {
  ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="gompertz",
               list('A' = 130, 'B' = 12, 'C' = 3, 'subA' = 20, 'subB' = 15, 'subC' = 0.25),
               df=logistic_df, type = "brms")
  expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "subA", "subB", "subC"))

  expect_s3_class(ss$formula, "brmsformula")
})

#* ************************************************************
#* *************** `Monomolecular growth modeling` *************** 
#* ************************************************************

set.seed(123)
mono_df<-growthSim("monomolecular", n=20, t=25,
                 params = list("A"=c(200,160), "B"=c(0.08, 0.1)))

test_that("Test monomolecular nls modeling", {
  ss<-suppressMessages(growthSS(model = "monomolecular", form=y~time|id/group, 
               df=mono_df, type = "nls"))
  expect_equal(as.numeric(unlist(ss$start)),
               c(58.0855457770439, 58.0855457770439, 0.0457652584520121, 0.0457652584520121) )
  
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nls")
  
  p <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test monomolecular nlrq modeling", {
  ss<-suppressMessages(growthSS(model = "monomolecular", form=y~time|id/group,
               df=mono_df, type = "nlrq"))
  expect_equal(as.numeric(unlist(ss$start)),
               c(58.0855457770439, 58.0855457770439, 0.0457652584520121, 0.0457652584520121) )
  
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "nlrq")
  
  p <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test monomolecular nlme modeling", {
  ss<-growthSS(model = "monomolecular", form=y~time|id/group, sigma="power",
               df=mono_df, type = "nlme")
  expect_equal(as.numeric(unlist(ss$start)),
               c(58.0855457770439, 58.0855457770439, 0.0457652584520121, 0.0457652584520121) )
  
  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "nlme")
  
  p <- suppressWarnings(growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df, boot = 3))
  expect_s3_class(p, "ggplot")
})

test_that("Test monomolecular brms model setup", {
  ss<-growthSS(model = "monomolecular", form=y~time|id/group, sigma="spline",
               list('A' = 130, 'B' = 0.1),
               df=mono_df, type = "brms")
  expect_equal(ss$prior$nlpar, c("", "", "A", "B"))
  
  expect_s3_class(ss$formula, "brmsformula")
})

#* ************************************************************
#* *************** `general additive growth modeling` *************** 
#* ************************************************************

set.seed(123)
gomp_df<-growthSim("gompertz", n=20, t=25,
                       params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(0.2, 0.25)))


test_that("Test nls gam modeling", {
  ss<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, 
                                df=gomp_df, type = "nls"))
  expect_equal(as.character(ss$formula), as.character(y ~ bs(time) * group) )
  
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "lm")
  
  p <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test nlrq gam modeling", {
  ss<-suppressMessages(growthSS(model = "gam", form=y~time|id/group,
                                df=gomp_df, type = "nlrq"))
  expect_equal(as.character(ss$formula), as.character(y ~ bs(time) * group) )
  
  fit <- fitGrowth(ss)
  expect_s3_class(fit, "rq")
  
  p <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})

test_that("Test nlme gam", {
  ss<-growthSS(model = "gam", form=y~time|id/group, sigma="power",
               df=gomp_df, type = "nlme")
  expect_equal(as.character(ss$formula$model), as.character(y ~ time * group))
  
  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "lme")
  
  p <- suppressWarnings(growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df, boot = 3))
  expect_s3_class(p, "ggplot")
})

test_that("Test mgcv gam", {
  ss<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=gomp_df, type = "mgcv"))
  expect_equal(as.character(ss$formula), as.character(y ~ 0+group+s(time, by = group) ) )
  
  fit <- suppressWarnings(fitGrowth(ss))
  expect_s3_class(fit, "gam")
  
  p <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
  expect_s3_class(p, "ggplot")
})


test_that("Test gam brms model setup", {
  ss<-growthSS(model = "gam", form=y~time|id/group, sigma="homo",
               df=gomp_df, type = "brms")
  
  expect_s3_class(ss$formula, "brmsformula")
})
