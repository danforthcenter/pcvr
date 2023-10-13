
if(FALSE){
  
  #* there are lots of options and until one obviously breaks I am not going to try to test all of them. 
  library(brms)
  test_that("Logistic brms model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("logistic", n=20, t=25,
                           params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
    
    ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="spline",
                 list('A' = 130, 'B' = 10, "C" = 3),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A", "B", "C"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/logistic_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  test_that("Gompertz brms model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("gompertz", n=20, t=25,
                           params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(0.25, 0.25)))
    
    ss<-growthSS(model = "gompertz", form=y~time|id/group, sigma="homo",
                 list('A' = 130, 'B' = 10, "C" = 1),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A", "B", "C"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/gompertz_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  test_that("Monomolecular brms model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("monomolecular", n=20, t=25,
                           params = list("A"=c(200,160), "B"=c(0.01, 0.08)))
    
    ss<-growthSS(model = "monomolecular", form=y~time|id/group, sigma="homo",
                 list('A' = 130, 'B' = 1),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A", "B"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/monomolecular_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  test_that("Exponential brms model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("exponential", n=20, t=25,
                     params = list("A"=c(15, 12), "B"=c(0.1, 0.085)))
    
    ss<-growthSS(model = "exponential", form=y~time|id/group, sigma="homo",
                 list('A' = 10, 'B' = 1),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A", "B"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/exponential_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("Power law brms model pipeline", {
    set.seed(123)
    simdf<-growthSim("power law", n=20, t=25,
                     params = list("A"=c(15, 12), "B"=c(0.75, 0.8)))
    
    ss<-growthSS(model = "power law", form=y~time|id/group, sigma="linear",
                 list('A' = 10, 'B' = 1),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A", "B"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/powerlaw_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("linear brms model pipeline", {
    set.seed(123)
    simdf<-growthSim("linear", n=20, t=25,
                     params = list("A"=c(15, 12)))
    
    ss<-growthSS(model = "linear", form=y~time|id/group, sigma="homo",
                 list('A' = 5),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linear_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("GAM brms model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("logistic", n=20, t=25,
                     params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
    
    ss<-growthSS(model = "gam", form=y~time|id/group, sigma="homo",
                 df=simdf, type = "brms")
    
    fit <- suppressWarnings(fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1))
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/gam_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  
  
  test_that("linear+linear brms model pipeline", {
    set.seed(123)
    simdf<-growthSim("linear + linear", n=20, t=25,
                     params = list("linear1A"=c(15, 12), "changePoint1"=c(8, 6), "linear2A"=c(3, 5) ))
    
    ss<-growthSS(model = "linear + linear", form=y~time|id/group, sigma="spline",
                 list("linear1A"=10, "changePoint1"=5, "linear2A"=2 ),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("","","linear1A","changePoint1", "linear2A"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linearPlusLinear_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("linear+logistic brms model pipeline", {
    set.seed(123)
    simdf<-growthSim("linear + logistic", n=20, t=25,
                     params = list("linear1A"=c(15, 12), "changePoint1"=c(8, 6),
                                   "logistic2A"=c(100, 150), "logistic2B"=c(10, 8), "logistic2C"=c(3, 2.5) ))
    
    ss<-growthSS(model = "linear + logistic", form=y~time|id/group, sigma="spline",
                 list("linear1A"=10, "changePoint1"=5, 
                      "logistic2A"=100, "logistic2B"=10, "logistic2C"=3),
                 df=simdf, type = "brms")
    
    expect_equal(ss$prior$nlpar, c("","","linear1A","changePoint1", 
                                   "logistic2A", "logistic2B", "logistic2C"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linearPlusLogistic_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  
  test_that("linear+gam brms model pipeline", {
    set.seed(123)
    simdf<-growthSim("linear + logistic", n=20, t=25, # using logistic data, but modeling as a gam
                     params = list("linear1A"=c(15, 12), "changePoint1"=c(8, 6),
                                   "logistic2A"=c(100, 150), "logistic2B"=c(10, 8), "logistic2C"=c(3, 2.5) ))
    
    ss<-growthSS(model = "linear + gam", form=y~time|id/group, sigma="homo",
                 list("linear1A"=10, "changePoint1"=5),
                 df=simdf, type = "brms")
    
    expect_equal(ss$prior$nlpar, c("","","linear1A","changePoint1"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linearPlusGAM_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  
  test_that("linear + linear + linear brms model pipeline", {
    set.seed(123)
    simdf<-growthSim("linear + linear + linear", n=25, t=50,
                     params = list("linear1A"=c(10, 12), "changePoint1"=c(8, 6), "linear2A"=c(1, 2), 
                                   "changePoint2"=c(25, 30), "linear3A"=c(20, 24)))
    
    ss<-growthSS(model = "linear + linear + linear", form=y~time|id/group, sigma="spline",
                 list("linear1A"=10, "changePoint1"=5, "linear2A"=2, "changePoint2"=15, "linear3A"=5 ),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("","","linear1A","changePoint1", "linear2A","changePoint2", "linear3A"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    # plot<-plot+coord_cartesian(xlim = c(0,50), ylim = c(0,500))
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linearPlusLinearPlusLinear_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  
  
  
}

