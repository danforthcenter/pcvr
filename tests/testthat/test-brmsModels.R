
if(file.exists("/home/josh/Desktop/") & interactive()){ # only run locally, don't test for each R-CMD Check
  
  #* there are lots of options and until one obviously breaks I am not going to try to test all of them. 
  library(testthat)
  library(brms)
  # devtools::load_all("/home/josh/Desktop/stargate/fahlgren_lab/pcvr")
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
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/logistic_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
    plot2 <- postPred(fit, form=ss$pcvrForm, groups=NULL, timeRange=NULL, hyp = "abs(diff) > 20", plot=FALSE)
    expect_s3_class(plot2$hypothesis, "data.frame")
    plot2 <- postPred(fit, form=ss$pcvrForm, groups=NULL, timeRange=10:30, hyp = "abs(diff) > 20", plot=TRUE)
    expect_s3_class(plot2$plot, "ggplot")
    plot2 <- postPred(fit, form=ss$pcvrForm, groups=c("b", "a"), timeRange=NULL, hyp = "abs(diff) > 20", plot=TRUE)
    expect_s3_class(plot2$plot, "ggplot")
    plot2 <- postPred(fit, form=ss$pcvrForm, groups=c("b", "a"), timeRange=10:30, hyp = "abs(diff) > 20", plot=TRUE)
    expect_s3_class(plot2$plot, "ggplot")
    
  })
  
  test_that("Gompertz brms model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("gompertz", n=20, t=25,
                           params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(0.25, 0.25)))
    
    ss<-growthSS(model = "gompertz", form=y~time|id/group, sigma="int",
                 list('A' = 130, 'B' = 10, "C" = 1),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A", "B", "C"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/gompertz_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  test_that("Monomolecular brms model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("monomolecular", n=20, t=25,
                           params = list("A"=c(200,160), "B"=c(0.01, 0.08)))
    
    ss<-growthSS(model = "monomolecular", form=y~time|id/group, sigma="int",
                 list('A' = 130, 'B' = 1),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A", "B"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/monomolecular_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  test_that("Exponential brms model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("exponential", n=20, t=25,
                     params = list("A"=c(15, 12), "B"=c(0.1, 0.085)))
    
    ss<-growthSS(model = "exponential", form=y~time|id/group, sigma="int",
                 list('A' = 10, 'B' = 1),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A", "B"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/exponential_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
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
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/powerlaw_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("linear brms model pipeline", {
    set.seed(123)
    simdf<-growthSim("linear", n=20, t=25,
                     params = list("A"=c(15, 12)))
    
    ss<-growthSS(model = "linear", form=y~time|id/group, sigma="int",
                 list('A' = 5),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "", "A"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linear_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("linear sub model with prior brms model pipeline", {
    set.seed(123)
    simdf<-growthSim("linear", n=20, t=25,
                     params = list("A"=c(15, 12)))
    
    model = "linear"; form=y~time|id/group; sigma="linear";
    priors = list('A' = 5, "subA"=2);
    df=simdf; type = "brms"
    
    ss<-growthSS(model = "linear", form=y~time|id/group, sigma="linear",
                 list('A' = 5, "subA"=2),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c( "", "A", "subA"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linear_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("GAM brms model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("logistic", n=20, t=25,
                     params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
    
    ss<-growthSS(model = "gam", form=y~time|id/group, sigma="int",
                 df=simdf, type = "brms")
    
    fit <- suppressWarnings(fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1))
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/gam_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
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
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linearPlusLinear_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
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
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linearPlusLogistic_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
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
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linearPlusGAM_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
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
    #ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/linearPlusLinearPlusLinear_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("Logistic brms logistic sub model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("logistic", n=20, t=25,
                     params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
    
    ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="logistic",
                 list('A' = 130, 'B' = 10, "C" = 3, 'subA' = 20, 'subB' = 10, "subC" = 2),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "subA", "subB", "subC"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1) # that's fast
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/logistic_logisticSubModel.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  test_that("Logistic brms gompertz sub model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("logistic", n=20, t=25,
                     params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
    
    ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="gompertz",
                 list('A' = 130, 'B' = 10, "C" = 3, 'subA' = 20, 'subB' = 10, "subC" = 2),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "subA", "subB", "subC"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1) # that's fast
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/logistic_gompSubModel.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  test_that("Logistic brms monomolecular sub model pipeline", {
    
    set.seed(123)
    simdf<-growthSim("logistic", n=20, t=25,
                     params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
    
    ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="monomolecular",
                 list('A' = 130, 'B' = 10, "C" = 3, 'subA' = 5, 'subB' = 0.5),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "A", "B", "C", "subA", "subB"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1) # that's fast
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/logistic_monoSubModel.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  
  test_that("int+int homoskedastic model pipeline", {
    
    set.seed(123)
    
    noise<-do.call(rbind, lapply(1:30, function(i){
      chngpt <- rnorm(2, 18, 2)
      rbind(data.frame(id = paste0("id_",i), time = 1:chngpt[1], group = "a", y = c(runif(chngpt[1]-1, 0, 20), rnorm(1,5,1))),
            data.frame(id = paste0("id_",i), time = 1:chngpt[2], group = "b", y = c(runif(chngpt[2]-1, 0, 20), rnorm(1,5,1))) )
    }))
    noise2<-do.call(rbind, lapply(1:30, function(i){
      start1 <- max(noise[noise$id == paste0("id_",i) & noise$group=="a", "time"])
      start2 <- max(noise[noise$id == paste0("id_",i) & noise$group=="b", "time"])
      
      rbind(data.frame(id = paste0("id_",i), time = start1:40, group = "a", y = c(runif(length(start1:40), 15, 50))),
            data.frame(id = paste0("id_",i), time = start2:40, group = "b", y = c(runif(length(start2:40), 15, 50))) )
    }))
    simdf <- rbind(noise, noise2)
    
    # ggplot(simdf, aes(x=time, y=y, color=group, group = interaction(group, id)))+
    #   geom_line()+
    #   theme_minimal()
    
    ss<-growthSS(model = "int + int", form=y~time|id/group, sigma="int",
                 list("int1"=10, "changePoint1"=10, "int2"=20),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("","","int1","changePoint1", "int2"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/intPlusInt_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("int+int fixed changepoint homoskedastic model pipeline", {
    
    set.seed(123)
    
    noise<-do.call(rbind, lapply(1:30, function(i){
      chngpt <- 20
      rbind(data.frame(id = paste0("id_",i), time = 1:chngpt[1], group = "a", y = c(runif(chngpt[1]-1, 0, 20), rnorm(1,5,1))),
            data.frame(id = paste0("id_",i), time = 1:chngpt[2], group = "b", y = c(runif(chngpt[2]-1, 0, 20), rnorm(1,5,1))) )
    }))
    noise2<-do.call(rbind, lapply(1:30, function(i){
      start1 <- max(noise[noise$id == paste0("id_",i) & noise$group=="a", "time"])
      start2 <- max(noise[noise$id == paste0("id_",i) & noise$group=="b", "time"])
      
      rbind(data.frame(id = paste0("id_",i), time = start1:40, group = "a", y = c(runif(length(start1:40), 15, 50))),
            data.frame(id = paste0("id_",i), time = start2:40, group = "b", y = c(runif(length(start2:40), 15, 50))) )
    }))
    simdf <- rbind(noise, noise2)
    
    # ggplot(simdf, aes(x=time, y=y, color=group, group = interaction(group, id)))+
    #   geom_line()+
    #   theme_minimal()
    
    ss<-growthSS(model = "int + int", form=y~time|id/group, sigma="int",
                 list("int1"=10, "fixedChangePoint1"=20, "int2"=20),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("","","int1", "int2"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/intPlusInt_fixedChngpt_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })  
  
  
  
  test_that("int+int thresholded homoskedasticity model pipeline", {
    
    set.seed(123)
    
    noise<-do.call(rbind, lapply(1:30, function(i){
      chngpt <- rnorm(2, 18, 2)
      rbind(data.frame(id = paste0("id_",i), time = 1:chngpt[1], group = "a", y = c(runif(chngpt[1]-1, 0, 20), rnorm(1,5,1))),
            data.frame(id = paste0("id_",i), time = 1:chngpt[2], group = "b", y = c(runif(chngpt[2]-1, 0, 20), rnorm(1,5,1))) )
    }))
    noise2<-do.call(rbind, lapply(1:30, function(i){
      start1 <- max(noise[noise$id == paste0("id_",i) & noise$group=="a", "time"])
      start2 <- max(noise[noise$id == paste0("id_",i) & noise$group=="b", "time"])
      
      rbind(data.frame(id = paste0("id_",i), time = start1:40, group = "a", y = c(runif(length(start1:40), 15, 50))),
            data.frame(id = paste0("id_",i), time = start2:40, group = "b", y = c(runif(length(start2:40), 15, 50))) )
    }))
    simdf <- rbind(noise, noise2)
    
    # ggplot(simdf, aes(x=time, y=y, color=group, group = interaction(group, id)))+
    #   geom_line()+
    #   theme_minimal()
    ss<-growthSS(model = "int + int", form=y~time|id/group, sigma="int + int",
                 list("int1"=10, "changePoint1"=10, "int2"=20, "subint1"=10, "subchangePoint1"=10, "subint2"=10),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "int1", "changePoint1", "int2", "subint1", "subchangePoint1", "subint2"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/intPlusInt_heteroskedastic_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  
  
  
  test_that("int + linear model and submodel pipeline", {
    
    set.seed(123)
    noise<-do.call(rbind, lapply(1:30, function(i){
      chngpt <- rnorm(2, 18, 2)
      rbind(data.frame(id = paste0("id_",i), time = 1:chngpt[1], group = "a", y = c(runif(chngpt[1]-1, 0, 20), rnorm(1,5,1))),
            data.frame(id = paste0("id_",i), time = 1:chngpt[2], group = "b", y = c(runif(chngpt[2]-1, 0, 20), rnorm(1,5,1))) )
    }))
    signal<-growthSim("linear", n=30, t=20,
                      params = list("A"=c(3, 5) ))
    signal<-do.call(rbind, lapply(unique(paste0(signal$id, signal$group)), function(int){
      noisesub<-noise[paste0(noise$id, noise$group)==int,]
      signalSub <- signal[paste0(signal$id, signal$group) == int, ]
      y_end <- noisesub[noisesub$time == max(noisesub$time), "y"]
      signalSub$time <- signalSub$time + max(noisesub$time)
      signalSub$y <- y_end + signalSub$y
      signalSub
    }))
    simdf <- rbind(noise, signal)
    
    # ggplot(simdf, aes(x=time, y=y, color=group, group = interaction(group, id)))+
    #   geom_line()+
    #   theme_minimal()
    
    ss<-growthSS(model = "int + linear", form=y~time|id/group, sigma="int + linear",
                 list("int1"=10, "changePoint1"=10, "linear2A"=20, "subint1"=10, "subchangePoint1"=10, "sublinear2A"=10),
                 df=simdf, type = "brms")
    expect_equal(ss$prior$nlpar, c("", "int1", "changePoint1", "linear2A", "subint1", "subchangePoint1", "sublinear2A"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df, timeRange = 1:40)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/intPlusLinear_heteroskedasticIntPlusLinear_fitGrowth.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  
  test_that("int + Logistic brms int+spline sub model pipeline", {
    
    set.seed(123)
    noise<-do.call(rbind, lapply(1:30, function(i){
      chngpt <- rnorm(2, 18, 2)
      rbind(data.frame(id = paste0("id_",i), time = 1:chngpt[1], group = "a", y = c(runif(chngpt[1]-1, 0, 20), rnorm(1,5,1))),
            data.frame(id = paste0("id_",i), time = 1:chngpt[2], group = "b", y = c(runif(chngpt[2]-1, 0, 20), rnorm(1,5,1))) )
    }))
    signal<-growthSim("logistic", n=20, t=30,
                     params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
    signal<-do.call(rbind, lapply(unique(paste0(signal$id, signal$group)), function(int){
      noisesub<-noise[paste0(noise$id, noise$group)==int,]
      signalSub <- signal[paste0(signal$id, signal$group) == int, ]
      y_end <- noisesub[noisesub$time == max(noisesub$time), "y"]
      signalSub$time <- signalSub$time + max(noisesub$time)
      signalSub$y <- y_end + signalSub$y
      signalSub
    }))
    simdf <- rbind(noise, signal)
    simdf<-simdf[simdf$time<45, ]
    
    # ggplot(simdf, aes(x=time, y=y, color=group, group = interaction(group, id)))+
    #   geom_line()+
    #   theme_minimal()
    
    ss<-growthSS(model = "int+logistic", form=y~time|id/group, sigma="int + spline",
                 list("int1" = 5, "changePoint1"=10 ,'logistic2A' = 130, 'logistic2B' = 10, "logistic2C" = 3,
                      'subint1' = 5, "subchangePoint1"=15),
                 df=simdf, type = "brms")
    
    
    expect_equal(ss$prior$nlpar, c("", "int1", "changePoint1", "logistic2A", "logistic2B", "logistic2C", 
                                   "subint1", "subchangePoint1"))
    
    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/intPluslogistic_intPlusGAMSubModel.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
    
  })
  
  
  test_that("fixed and estimated changepoints can be mixed in growth formula", {
    
    simdf1 <- growthSim(model = "logistic", n=20, t=20, params = list("A"=c(180, 160), "B"=c(9,11), "C"=c(3,3.5)))
    
    simdf2 <- growthSim(model = "linear + linear", n=20, t=20,
                        params = list("linear1A"= c(6,8), "changePoint1" = c(7,9), "linear2A"=c(15, 20)))
    
    simdf2_adj<-do.call(rbind, lapply(unique(paste0(simdf2$id, simdf2$group)), function(int){
      p1<-simdf1[paste0(simdf1$id, simdf1$group)==int,]
      p2 <- simdf2[paste0(simdf2$id, simdf2$group) == int, ]
      y_end_p1 <- p1[p1$time == max(p1$time), "y"]
      p2$time <- p2$time + max(p1$time)
      p2$y <- y_end_p1 + p2$y
      p2
    }))
    simdf<-rbind(simdf1, simdf2_adj)

    ss <- growthSS(model = "logistic+linear+linear", form=y~time|id/group,
                   sigma = "logistic+linear", df = simdf,
                   start = list("logistic1A" = 130, "logistic1B" = 10, "logistic1C" = 3.5,
                                "fixedChangePoint1" = 20, "linear2A" = 5, "changePoint2" = 28, "linear3A" = 20,
                                "sublogistic1A" = 10, "sublogistic1B" = 12, "sublogistic1C" = 20,
                                "subfixedChangePoint1"=20, "sublinear2A"=3), type = "brms")
    
    expect_equal(ss$prior$nlpar, c("", "logistic1A", "logistic1B", "logistic1C", "linear2A", "changePoint2", 
                                   "linear3A", "sublogistic1A", "sublogistic1B", "sublogistic1C", 
                                   "sublinear2A"))

    fit <- fitGrowth(ss, backend="cmdstanr", iter=500, chains=1, cores=1)
    
    expect_s3_class(fit, "brmsfit")
    
    plot <- growthPlot(fit=fit, form=ss$pcvrForm, df = ss$df)
    ggsave("/home/josh/Desktop/stargate/fahlgren_lab/labMeetings/threePart_fixedAndEstimatedChangepoint.png", plot, width=10, height=6, dpi=300, bg="#ffffff")
    expect_s3_class(plot, "ggplot")
  })
  
  test_that("postPred can compare multiple models", {
    set.seed(123)
    simdf1<-growthSim("logistic", n=20, t=25,
                     params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
    
    ss1<-growthSS(model = "logistic", form=y~time|id/group, sigma="spline",
                 list('A' = 130, 'B' = 10, "C" = 3),
                 df=simdf1, type = "brms")
    
    fit1 <- fitGrowth(ss1, backend="cmdstanr", iter=500, chains=1, cores=1)
   
    simdf2<-growthSim("monomolecular", n=20, t=25,
                     params = list("A"=c(200,160), "B"=c(0.01, 0.08)))
    
    ss2<-growthSS(model = "monomolecular", form=y~time|id/group, sigma="int",
                 list('A' = 130, 'B' = 1),
                 df=simdf2, type = "brms")
    
    fit2 <- fitGrowth(ss2, backend="cmdstanr", iter=550, chains=1, cores=1)
    
    
    post <- postPred(fit = list(fit1, fit2), form=ss$pcvrForm, groups=c("a", "a"),
                     timeRange=NULL, hyp = "abs(diff) > 20", plot=TRUE)
    expect_s3_class(post$plot, "ggplot")
    
  })
  
  
}

