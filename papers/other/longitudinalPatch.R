rm(list=ls())
library(ggplot2)
devtools::load_all("~/pcvr")
#* *******************************************************
#* ***** `Figure 3: Longitudinal Modeling Pipeline`
#* *******************************************************
#* Need to add EVD (weibull, gumbel, frechet)
#* Want to add survival (in brms, nls-as survreg)? Not sure if that fits into the general patchwork well.

#* ***** `Logistic Data`
set.seed(123)
logistic_df<-growthSim("logistic", n=20, t=25,
                       params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#* `nls`
ss_logistic_nls<-suppressMessages(growthSS(model = "logistic", form=y~time|id/group, 
                                           df=logistic_df, type = "nls"))
fit_logistic_nls <- fitGrowth(ss_logistic_nls)
plot_logistic_nls <- growthPlot(fit=fit_logistic_nls, form=ss_logistic_nls$pcvrForm, df = ss_logistic_nls$df) +
  ggplot2::labs(title="nls", y = "Logistic")+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20), plot.title = element_text(size=25), legend.position="none")
#* `nlrq`
ss_logistic_nlrq<-suppressMessages(growthSS(model = "logistic", form=y~time|id/group,
                                            df=logistic_df, type = "nlrq", tau = seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
fit_logistic_nlrq <- fitGrowth(ss_logistic_nlrq)
plot_logistic_nlrq <- growthPlot(fit=fit_logistic_nlrq, form=ss_logistic_nlrq$pcvrForm, df = ss_logistic_nlrq$df)+
  ggplot2::labs(title="nlrq")+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size=25))
#* `nlme`
ss_logistic_nlme<-growthSS(model = "logistic", form=y~time|id/group, sigma="power",
                           df=logistic_df, type = "nlme")
fit_logistic_nlme <- suppressWarnings(fitGrowth(ss_logistic_nlme))
plot_logistic_nlme <- suppressWarnings(growthPlot(fit=fit_logistic_nlme, form=ss_logistic_nlme$pcvrForm, df = ss_logistic_nlme$df))+
  ggplot2::labs(title="nlme")+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size=25))
#* `mgcv`
ss_logistic_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=logistic_df, type = "mgcv"))
fit_logistic_mgcv <- suppressWarnings(fitGrowth(ss_logistic_mgcv))
plot_logistic_mgcv <- growthPlot(fit=fit_logistic_mgcv, form=ss_logistic_mgcv$pcvrForm, df = ss_logistic_mgcv$df)+
  ggplot2::labs(title="mgcv")+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size=25), legend.position="none")
#* `brms`
ss_logistic_brms<-growthSS(model = "logistic", form=y~time|id/group, sigma="spline",
                           list('A' = 130, 'B' = 10, "C" = 3),
                           df=logistic_df, type = "brms")
fit_logistic_brms <- fitGrowth(ss_logistic_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_logistic_brms <- growthPlot(fit=fit_logistic_brms, form=ss_logistic_brms$pcvrForm, df = ss_logistic_brms$df)+
  ggplot2::labs(title="brms")+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size=25))


#* ***** `Gompertz Data`
set.seed(123)
gompertz_df<-growthSim("gompertz", n=20, t=25,
                       params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(0.2, 0.25)))
#* `nls`
ss_gompertz_nls<-suppressMessages(growthSS(model = "gompertz", form=y~time|id/group, 
                                           df=gompertz_df, type = "nls"))
ss_gompertz_nls$start <- list("A" = c(200, 160), "B" = c(13, 11), "C"=c(0.2, 0.25) )
fit_gompertz_nls <- fitGrowth(ss_gompertz_nls)
plot_gompertz_nls <- growthPlot(fit=fit_gompertz_nls, form=ss_gompertz_nls$pcvrForm, df = ss_gompertz_nls$df) +
  ggplot2::labs(y = "Gompertz")+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20), legend.position="none")
#* `nlrq`
ss_gompertz_nlrq<-suppressMessages(growthSS(model = "gompertz", form=y~time|id/group,
                                            df=gompertz_df, type = "nlrq", tau = seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
ss_gompertz_nlrq$start <- list("A" = c(200, 160), "B" = c(13, 11), "C"=c(0.2, 0.25) )
fit_gompertz_nlrq <- fitGrowth(ss_gompertz_nlrq)
plot_gompertz_nlrq <- growthPlot(fit=fit_gompertz_nlrq, form=ss_gompertz_nlrq$pcvrForm, df = ss_gompertz_nlrq$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `nlme`
ss_gompertz_nlme<-growthSS(model = "gompertz", form=y~time|id/group, sigma="none",
                           df=gompertz_df, type = "nlme")
ss_gompertz_nlme$start <- c(c(200, 160), c(13, 11), c(0.2, 0.25) )
fit_gompertz_nlme <- suppressWarnings(fitGrowth(ss_gompertz_nlme))
plot_gompertz_nlme <- suppressWarnings(growthPlot(fit=fit_gompertz_nlme, form=ss_gompertz_nlme$pcvrForm, df = ss_gompertz_nlme$df))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `mgcv`
ss_gompertz_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=gompertz_df, type = "mgcv"))
fit_gompertz_mgcv <- suppressWarnings(fitGrowth(ss_gompertz_mgcv))
plot_gompertz_mgcv <- growthPlot(fit=fit_gompertz_mgcv, form=ss_gompertz_mgcv$pcvrForm, df = ss_gompertz_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none")
#* `brms`
ss_gompertz_brms<-growthSS(model = "gompertz", form=y~time|id/group, sigma="spline",
                           list('A' = 130, 'B' = 10, "C" = 1),
                           df=gompertz_df, type = "brms")
fit_gompertz_brms <- fitGrowth(ss_gompertz_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_gompertz_brms <- growthPlot(fit=fit_gompertz_brms, form=ss_gompertz_brms$pcvrForm, df = ss_gompertz_brms$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())





#* ***** `Weibull Data`
set.seed(123)
weibull_df<-growthSim("weibull", params = list("A"=c(200,160), "B"=c(1, 0.75), "C"=c(2, 3)) )
#* `nls`
ss_weibull_nls<-suppressMessages(growthSS(model = "weibull", form=y~time|id/group,
                                           df=weibull_df, type = "nls"))
fit_weibull_nls <- fitGrowth(ss_weibull_nls)
plot_weibull_nls <- growthPlot(fit=fit_weibull_nls, form=ss_weibull_nls$pcvrForm, df = ss_weibull_nls$df) +
  ggplot2::labs(y = "Weibull")+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20), legend.position="none")
# #* `nlrq`
ss_weibull_nlrq<-suppressMessages(growthSS(model = "weibull", form=y~time|id/group,
                                            df=weibull_df, type = "nlrq", tau = seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
fit_weibull_nlrq <- fitGrowth(ss_weibull_nlrq)
plot_weibull_nlrq <- growthPlot(fit=fit_weibull_nlrq, form=ss_weibull_nlrq$pcvrForm, df = ss_weibull_nlrq$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
# #* `nlme`
ss_weibull_nlme<-growthSS(model = "weibull", form=y~time|id/group, sigma="none",
                           df=weibull_df, type = "nlme")
fit_weibull_nlme <- suppressWarnings(fitGrowth(ss_weibull_nlme))
plot_weibull_nlme <- suppressWarnings(growthPlot(fit=fit_weibull_nlme, form=ss_weibull_nlme$pcvrForm, df = ss_weibull_nlme$df))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
# #* `mgcv`
ss_weibull_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=weibull_df, type = "mgcv"))
fit_weibull_mgcv <- suppressWarnings(fitGrowth(ss_weibull_mgcv))
plot_weibull_mgcv <- growthPlot(fit=fit_weibull_mgcv, form=ss_weibull_mgcv$pcvrForm, df = ss_weibull_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none")
# #* `brms`
ss_weibull_brms<-growthSS(model = "weibull", form=y~time|id/group, sigma="spline",
                           list('A' = 130, 'B' = 1, "C" = 1),
                           df=weibull_df, type = "brms")
fit_weibull_brms <- fitGrowth(ss_weibull_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_weibull_brms <- growthPlot(fit=fit_weibull_brms, form=ss_weibull_brms$pcvrForm, df = ss_weibull_brms$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())

#* ***** `Frechet Data`
set.seed(123)
frechet_df<-growthSim("frechet", params = list("A"=c(200,180), "B"=c(5, 3), "C"=c(6, 4)) )
#* `nls`
ss_frechet_nls<-suppressMessages(growthSS(model = "frechet", form=y~time|id/group,
                                          df=frechet_df, type = "nls"))
fit_frechet_nls <- fitGrowth(ss_frechet_nls)
plot_frechet_nls <- growthPlot(fit=fit_frechet_nls, form=ss_frechet_nls$pcvrForm, df = ss_frechet_nls$df) +
  ggplot2::labs(y = "Frechet")+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20), legend.position="none")
# #* `nlrq`
ss_frechet_nlrq<-suppressMessages(growthSS(model = "frechet", form=y~time|id/group,
                                           df=frechet_df, type = "nlrq", tau = seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
fit_frechet_nlrq <- fitGrowth(ss_frechet_nlrq)
plot_frechet_nlrq <- growthPlot(fit=fit_frechet_nlrq, form=ss_frechet_nlrq$pcvrForm, df = ss_frechet_nlrq$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
# #* `nlme`
ss_frechet_nlme<-growthSS(model = "frechet", form=y~time|id/group, sigma="none",
                          df=frechet_df, type = "nlme")
fit_frechet_nlme <- suppressWarnings(fitGrowth(ss_frechet_nlme))
plot_frechet_nlme <- suppressWarnings(growthPlot(fit=fit_frechet_nlme, form=ss_frechet_nlme$pcvrForm, df = ss_frechet_nlme$df))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
# #* `mgcv`
ss_frechet_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=frechet_df, type = "mgcv"))
fit_frechet_mgcv <- suppressWarnings(fitGrowth(ss_frechet_mgcv))
plot_frechet_mgcv <- growthPlot(fit=fit_frechet_mgcv, form=ss_frechet_mgcv$pcvrForm, df = ss_frechet_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none")
# #* `brms`
ss_frechet_brms<-growthSS(model = "frechet", form=y~time|id/group, sigma="spline",
                          list('A' = 130, 'B' = 4, "C" = 4),
                          df=frechet_df, type = "brms")
fit_frechet_brms <- fitGrowth(ss_frechet_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_frechet_brms <- growthPlot(fit=fit_frechet_brms, form=ss_frechet_brms$pcvrForm, df = ss_frechet_brms$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())


#* ***** `Gumbel Data`
set.seed(123)
gumbel_df <- growthSim("gumbel", params = list("A"=c(120,140), "B"=c(6, 5), "C"=c(4, 3)) )
#* `nls`
ss_gumbel_nls<-suppressMessages(growthSS(model = "gumbel", form=y~time|id/group,
                                          df=gumbel_df, type = "nls"))
fit_gumbel_nls <- fitGrowth(ss_gumbel_nls)
plot_gumbel_nls <- growthPlot(fit=fit_gumbel_nls, form=ss_gumbel_nls$pcvrForm, df = ss_gumbel_nls$df) +
  ggplot2::labs(y = "EVD (1,2,3)")+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20), legend.position="none")
#* `nlrq`
ss_gumbel_nlrq<-suppressMessages(growthSS(model = "gumbel", form=y~time|id/group,
                                           df=gumbel_df, type = "nlrq", tau = seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
fit_gumbel_nlrq <- fitGrowth(ss_gumbel_nlrq)
plot_gumbel_nlrq <- growthPlot(fit=fit_gumbel_nlrq, form=ss_gumbel_nlrq$pcvrForm, df = ss_gumbel_nlrq$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `nlme`
ss_gumbel_nlme<-growthSS(model = "gumbel", form=y~time|id/group, sigma="none",
                          df=gumbel_df, type = "nlme")
fit_gumbel_nlme <- suppressWarnings(fitGrowth(ss_gumbel_nlme))
plot_gumbel_nlme <- suppressWarnings(growthPlot(fit=fit_gumbel_nlme, form=ss_gumbel_nlme$pcvrForm, df = ss_gumbel_nlme$df))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `mgcv`
ss_gumbel_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=gumbel_df, type = "mgcv"))
fit_gumbel_mgcv <- suppressWarnings(fitGrowth(ss_gumbel_mgcv))
plot_gumbel_mgcv <- growthPlot(fit=fit_gumbel_mgcv, form=ss_gumbel_mgcv$pcvrForm, df = ss_gumbel_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none")
#* `brms`
ss_gumbel_brms<-growthSS(model = "gumbel", form=y~time|id/group, sigma="spline",
                          list('A' = 130, 'B' = 4, "C" = 4),
                          df=gumbel_df, type = "brms")
fit_gumbel_brms <- fitGrowth(ss_gumbel_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_gumbel_brms <- growthPlot(fit=fit_gumbel_brms, form=ss_gumbel_brms$pcvrForm, df = ss_gumbel_brms$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())

#* ***** `Monomolecular Data`
set.seed(123)
mm_df<-growthSim("monomolecular", n=20, t=25,
                 params = list("A"=c(200,160), "B"=c(0.08, 0.1)))
#* `nls`
ss_mm_nls<-suppressMessages(growthSS(model = "monomolecular", form=y~time|id/group, 
                                     df=mm_df, type = "nls"))
fit_mm_nls <- fitGrowth(ss_mm_nls)
plot_mm_nls <- growthPlot(fit=fit_mm_nls, form=ss_mm_nls$pcvrForm, df = ss_mm_nls$df) +
  ggplot2::labs(y = "Monomolecular")+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20), legend.position="none")
#* `nlrq`
ss_mm_nlrq<-suppressMessages(growthSS(model = "monomolecular", form=y~time|id/group,
                                      df=mm_df, type = "nlrq", tau = seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
fit_mm_nlrq <- fitGrowth(ss_mm_nlrq)
plot_mm_nlrq <- growthPlot(fit=fit_mm_nlrq, form=ss_mm_nlrq$pcvrForm, df = ss_mm_nlrq$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `nlme`
ss_mm_nlme<-growthSS(model = "monomolecular", form=y~time|id/group, sigma="power",
                     df=mm_df, type = "nlme")
fit_mm_nlme <- suppressWarnings(fitGrowth(ss_mm_nlme))
plot_mm_nlme <- suppressWarnings(growthPlot(fit=fit_mm_nlme, form=ss_mm_nlme$pcvrForm, df = ss_mm_nlme$df))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `mgcv`
ss_monomolecular_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=mm_df, type = "mgcv"))
fit_monomolecular_mgcv <- suppressWarnings(fitGrowth(ss_monomolecular_mgcv))
plot_monomolecular_mgcv <- growthPlot(fit=fit_monomolecular_mgcv, form=ss_monomolecular_mgcv$pcvrForm, df = ss_monomolecular_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none")
#* `brms`
ss_mm_brms<-growthSS(model = "monomolecular", form=y~time|id/group, sigma="spline",
                     list('A' = 130, 'B' = 0.25),
                     df=mm_df, type = "brms")
fit_mm_brms <- fitGrowth(ss_mm_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_mm_brms <- growthPlot(fit=fit_mm_brms, form=ss_mm_brms$pcvrForm, df = ss_mm_brms$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())


#* ***** `Exponential Data`
set.seed(123)
exponential_df<-growthSim("exponential", n=20, t=25,
                          params = list("A"=c(15, 20), "B"=c(0.085, 0.09)))
#* `nls`
ss_exponential_nls<-suppressMessages(growthSS(model = "exponential", form=y~time|id/group, 
                                              df=exponential_df, type = "nls"))
fit_exponential_nls <- fitGrowth(ss_exponential_nls)
plot_exponential_nls <- growthPlot(fit=fit_exponential_nls, form=ss_exponential_nls$pcvrForm, df = ss_exponential_nls$df) +
  ggplot2::labs(y = "Exponential")+
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20), legend.position="none")
#* `nlrq`
ss_exponential_nlrq<-suppressMessages(growthSS(model = "exponential", form=y~time|id/group,
                                               df=exponential_df, type = "nlrq", tau = seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
fit_exponential_nlrq <- fitGrowth(ss_exponential_nlrq)
plot_exponential_nlrq <- growthPlot(fit=fit_exponential_nlrq, form=ss_exponential_nlrq$pcvrForm, df = ss_exponential_nlrq$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `nlme`
ss_exponential_nlme<-growthSS(model = "exponential", form=y~time|id/group, sigma="power",
                              df=exponential_df, type = "nlme")
fit_exponential_nlme <- suppressWarnings(fitGrowth(ss_exponential_nlme))
plot_exponential_nlme <- suppressWarnings(growthPlot(fit=fit_exponential_nlme, form=ss_exponential_nlme$pcvrForm, df = ss_exponential_nlme$df))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `mgcv`
ss_exponential_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=exponential_df, type = "mgcv"))
fit_exponential_mgcv <- suppressWarnings(fitGrowth(ss_exponential_mgcv))
plot_exponential_mgcv <- growthPlot(fit=fit_exponential_mgcv, form=ss_exponential_mgcv$pcvrForm, df = ss_exponential_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none")
#* `brms`
ss_exponential_brms<-growthSS(model = "exponential", form=y~time|id/group, sigma="spline",
                              list('A' = 10, 'B' = 0.1),
                              df=exponential_df, type = "brms")
fit_exponential_brms <- fitGrowth(ss_exponential_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_exponential_brms <- growthPlot(fit=fit_exponential_brms, form=ss_exponential_brms$pcvrForm, df = ss_exponential_brms$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())



#* ***** `Power Law Data`
set.seed(123)
powerlaw_df<-growthSim("power law", n=20, t=25,
                       params = list("A"=c(16, 11), "B"=c(0.75, 0.7)))
#* `nls`
ss_powerlaw_nls<-suppressMessages(growthSS(model = "power law", form=y~time|id/group, 
                                           df=powerlaw_df, type = "nls"))
fit_powerlaw_nls <- fitGrowth(ss_powerlaw_nls)
plot_powerlaw_nls <- growthPlot(fit=fit_powerlaw_nls, form=ss_powerlaw_nls$pcvrForm, df = ss_powerlaw_nls$df) +
  ggplot2::labs(y = "Power Law")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20), legend.position="none")
#* `nlrq`
ss_powerlaw_nlrq<-suppressMessages(growthSS(model = "power law", form=y~time|id/group,
                                            df=powerlaw_df, type = "nlrq", tau = seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
fit_powerlaw_nlrq <- fitGrowth(ss_powerlaw_nlrq)
plot_powerlaw_nlrq <- growthPlot(fit=fit_powerlaw_nlrq, form=ss_powerlaw_nlrq$pcvrForm, df = ss_powerlaw_nlrq$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `nlme`
ss_powerlaw_nlme<-growthSS(model = "power law", form=y~time|id/group, sigma="power",
                           df=powerlaw_df, type = "nlme")
fit_powerlaw_nlme <- suppressWarnings(fitGrowth(ss_powerlaw_nlme))
plot_powerlaw_nlme <- suppressWarnings(growthPlot(fit=fit_powerlaw_nlme, form=ss_powerlaw_nlme$pcvrForm, df = ss_powerlaw_nlme$df))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `mgcv`
ss_powerlaw_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=powerlaw_df, type = "mgcv"))
fit_powerlaw_mgcv <- suppressWarnings(fitGrowth(ss_powerlaw_mgcv))
plot_powerlaw_mgcv <- growthPlot(fit=fit_powerlaw_mgcv, form=ss_powerlaw_mgcv$pcvrForm, df = ss_powerlaw_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), legend.position="none")
#* `brms`
ss_powerlaw_brms<-growthSS(model = "power law", form=y~time|id/group, sigma="spline",
                           list('A' = 10, 'B' = 1),
                           df=powerlaw_df, type = "brms")
fit_powerlaw_brms <- fitGrowth(ss_powerlaw_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_powerlaw_brms <- growthPlot(fit=fit_powerlaw_brms, form=ss_powerlaw_brms$pcvrForm, df = ss_powerlaw_brms$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())

#* ***** `Linear Data`
set.seed(123)
linear_df<-growthSim("linear", n=20, t=25,
                     params = list("A"=c(7.5, 8.3)))
#* `nls`
ss_linear_nls<-suppressMessages(growthSS(model = "linear", form=y~time|id/group, 
                                         df=linear_df, type = "nls"))
fit_linear_nls <- fitGrowth(ss_linear_nls)
plot_linear_nls <- growthPlot(fit=fit_linear_nls, form=ss_linear_nls$pcvrForm, df = ss_linear_nls$df) +
  ggplot2::labs(y = "Linear", x="time")+
  theme(axis.title.y = element_text(size = 20), legend.position="none")
#* `nlrq`
ss_linear_nlrq<-suppressMessages(growthSS(model = "linear", form=y~time|id/group,
                                          df=linear_df, type = "nlrq", tau = seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
fit_linear_nlrq <- fitGrowth(ss_linear_nlrq)
plot_linear_nlrq <- growthPlot(fit=fit_linear_nlrq, form=ss_linear_nlrq$pcvrForm, df = ss_linear_nlrq$df)+
  ggplot2::labs(x="time")+
  theme(axis.title.y = element_blank())
#* `nlme`
ss_linear_nlme<-growthSS(model = "linear", form=y~time|id/group, sigma="none",
                         df=linear_df, type = "nlme")
fit_linear_nlme <- suppressWarnings(fitGrowth(ss_linear_nlme))
plot_linear_nlme <- suppressWarnings(growthPlot(fit=fit_linear_nlme, form=ss_linear_nlme$pcvrForm, df = ss_linear_nlme$df))+
  ggplot2::labs(x="time")+
  theme(axis.title.y = element_blank())
#* `mgcv`
ss_linear_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=linear_df, type = "mgcv"))
fit_linear_mgcv <- suppressWarnings(fitGrowth(ss_linear_mgcv))
plot_linear_mgcv <- growthPlot(fit=fit_linear_mgcv, form=ss_linear_mgcv$pcvrForm, df = ss_linear_mgcv$df)+
  ggplot2::labs(x="time")+
  theme(axis.title.y = element_blank(), legend.position="none")
#* `brms`
ss_linear_brms<-growthSS(model = "linear", form=y~time|id/group, sigma="spline",
                         list('A' = 10),
                         df=linear_df, type = "brms")
fit_linear_brms <- fitGrowth(ss_linear_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_linear_brms <- growthPlot(fit=fit_linear_brms, form=ss_linear_brms$pcvrForm, df = ss_linear_brms$df)+
  ggplot2::labs(x="time")+
  theme(axis.title.y = element_blank())


#* *******************************************************
#* ***** `Patchwork assembling plots`
#* *******************************************************
#* add rows for EVDs


layout <- c(
  # row 1
  area(1,1,3,3),
  area(1,4,3,6),
  area(1,7,3,9),
  area(1,10,3,12),
  area(1,13,3,15),
  # row 2
  area(4,1,6,3),
  area(4,4,6,6),
  area(4,7,6,9),
  area(4,10,6,12),
  area(4,13,6,15),
  # row 3
  area(7,1,9,3),
  area(7,4,9,6),
  area(7,7,9,9),
  area(7,10,9,12),
  area(7,13,9,15),
  # row 4
  area(10,1,12,3),
  area(10,4,12,6),
  area(10,7,12,9),
  area(10,10,12,12),
  area(10,13,12,15),
  # row 5
  area(13,1,15,3),
  area(13,4,15,6),
  area(13,7,15,9),
  area(13,10,15,12),
  area(13,13,15,15),
  # row 6
  area(16,1,18,3),
  area(16,4,18,6),
  area(16,7,18,9),
  area(16,10,18,12),
  area(16,13,18,15),
  # row 7
  area(19,1,21,3),
  area(19,4,21,6),
  area(19,7,21,9),
  area(19,10,21,12),
  area(19,13,21,15)
)
plot(layout)

patch <- plot_logistic_nls + plot_logistic_nlrq + plot_logistic_nlme + 
  plot_logistic_mgcv + plot_logistic_brms +
  plot_gompertz_nls + plot_gompertz_nlrq + plot_gompertz_nlme + 
  plot_gompertz_mgcv + plot_gompertz_brms + 
  plot_gumbel_nls + plot_gumbel_nlrq + plot_gumbel_nlme + 
  plot_gumbel_mgcv + plot_gumbel_brms + 
  plot_mm_nls + plot_mm_nlrq + plot_mm_nlme + 
  plot_monomolecular_mgcv + plot_mm_brms +
  plot_exponential_nls + plot_exponential_nlrq + plot_exponential_nlme +
  plot_exponential_mgcv + plot_exponential_brms + 
  plot_powerlaw_nls + plot_powerlaw_nlrq + plot_powerlaw_nlme +
  plot_powerlaw_mgcv + plot_powerlaw_brms +
  plot_linear_nls + plot_linear_nlrq + plot_linear_nlme + 
  plot_linear_mgcv + plot_linear_brms +
  patchwork::plot_layout(design = layout)

ggsave("~/pcvr/papers/other/pcvr_longitudinal_patch.png",
       patch, width= 22, height=17, dpi = 300, bg = "#ffffff")

save(list=ls(), file = "~/pcvr/papers/other/pcvr_longitudinal_patch.rdata")
