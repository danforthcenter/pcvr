
#* new figures for lab meeting 10/16/2023

library(pcvr)
library(ggplot2)
library(patchwork)

#* ****************************************
#* ***** `Figure 1: pcvBox Example` *****
#* ****************************************

set.seed(678)
s1 = rnorm(20, 50,10)
s2= rnorm(20, 60,8)
df <- data.frame(y = c(s1,s2), x = c(rep("A",20), rep("B",20)))
(pb<-pcvBox(df, x="x", y="y", compare="A", showPoints = T))
ggsave("~/Desktop/stargate/fahlgren_lab/labMeetings/pcvr_pcvBox_ex.png",
       pb, width= 5, height=5, dpi = 300, bg = "#ffffff")


#* ****************************************
#* ***** `Figure 2: Conjugate Example` *****
#* ****************************************

set.seed(678)
s1 = rnorm(20, 50,10)
s2= rnorm(20, 60,8)

gaussianMeans_sv_ex <- conjugate(s1=s1, s2= s2, method="t",
                                 priors = list( mu=c(50, 50),n=c(1,1),s2=c(15,15) ),
                                 plot=TRUE, rope_range = c(-5,5), rope_ci = 0.89, 
                                 cred.int.level = 0.89, hypothesis="greater", support=NULL)

gaussianMeans_sv_ex$plot

ggsave("~/Desktop/stargate/fahlgren_lab/labMeetings/pcvr_conjugate_ex.png",
       gaussianMeans_sv_ex$plot, width= 9, height=5.5, dpi = 300, bg = "#ffffff")

#* using frem from bellwether vignette
# x <- x +theme(legend.position=c(0.9, 0.1)) +scale_x_continuous(breaks =c (5, 10, 15, 20))
# ggsave("~/Desktop/stargate/fahlgren_lab/labMeetings/pcvr_variancePartitioning.png",
#        x, width= 12, height=8, dpi = 300, bg = "#ffffff")

#* *******************************************************
#* ***** `Figure 3: Longitudinal Modeling Pipeline`
#* *******************************************************

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
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20), plot.title = element_text(size=25))
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
plot_logistic_nlme <- suppressWarnings(growthPlot(fit=fit_logistic_nlme, form=ss_logistic_nlme$pcvrForm, df = ss_logistic_nlme$df, boot = 3))+
  ggplot2::labs(title="nlme")+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size=25))
#* `mgcv`
ss_logistic_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=logistic_df, type = "mgcv"))
fit_logistic_mgcv <- suppressWarnings(fitGrowth(ss_logistic_mgcv))
plot_logistic_mgcv <- growthPlot(fit=fit_logistic_mgcv, form=ss_logistic_mgcv$pcvrForm, df = ss_logistic_mgcv$df)+
  ggplot2::labs(title="mgcv")+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title = element_text(size=25))
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
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20))
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
plot_gompertz_nlme <- suppressWarnings(growthPlot(fit=fit_gompertz_nlme, form=ss_gompertz_nlme$pcvrForm, df = ss_gompertz_nlme$df, boot = 3))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `mgcv`
ss_gompertz_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=gompertz_df, type = "mgcv"))
fit_gompertz_mgcv <- suppressWarnings(fitGrowth(ss_gompertz_mgcv))
plot_gompertz_mgcv <- growthPlot(fit=fit_gompertz_mgcv, form=ss_gompertz_mgcv$pcvrForm, df = ss_gompertz_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `brms`
ss_gompertz_brms<-growthSS(model = "gompertz", form=y~time|id/group, sigma="spline",
                           list('A' = 130, 'B' = 10, "C" = 1),
                           df=gompertz_df, type = "brms")
fit_gompertz_brms <- fitGrowth(ss_gompertz_brms, backend="cmdstanr", iter=500, chains=1, cores=1)
plot_gompertz_brms <- growthPlot(fit=fit_gompertz_brms, form=ss_gompertz_brms$pcvrForm, df = ss_gompertz_brms$df)+
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
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20))
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
plot_mm_nlme <- suppressWarnings(growthPlot(fit=fit_mm_nlme, form=ss_mm_nlme$pcvrForm, df = ss_mm_nlme$df, boot = 3))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `mgcv`
ss_monomolecular_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=mm_df, type = "mgcv"))
fit_monomolecular_mgcv <- suppressWarnings(fitGrowth(ss_monomolecular_mgcv))
plot_monomolecular_mgcv <- growthPlot(fit=fit_monomolecular_mgcv, form=ss_monomolecular_mgcv$pcvrForm, df = ss_monomolecular_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
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
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 20))
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
plot_exponential_nlme <- suppressWarnings(growthPlot(fit=fit_exponential_nlme, form=ss_exponential_nlme$pcvrForm, df = ss_exponential_nlme$df, boot = 3))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `mgcv`
ss_exponential_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=exponential_df, type = "mgcv"))
fit_exponential_mgcv <- suppressWarnings(fitGrowth(ss_exponential_mgcv))
plot_exponential_mgcv <- growthPlot(fit=fit_exponential_mgcv, form=ss_exponential_mgcv$pcvrForm, df = ss_exponential_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
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
        axis.title.y = element_text(size = 20))
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
plot_powerlaw_nlme <- suppressWarnings(growthPlot(fit=fit_powerlaw_nlme, form=ss_powerlaw_nlme$pcvrForm, df = ss_powerlaw_nlme$df, boot = 3))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
#* `mgcv`
ss_powerlaw_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=powerlaw_df, type = "mgcv"))
fit_powerlaw_mgcv <- suppressWarnings(fitGrowth(ss_powerlaw_mgcv))
plot_powerlaw_mgcv <- growthPlot(fit=fit_powerlaw_mgcv, form=ss_powerlaw_mgcv$pcvrForm, df = ss_powerlaw_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
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
  theme(axis.title.y = element_text(size = 20))
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
plot_linear_nlme <- suppressWarnings(growthPlot(fit=fit_linear_nlme, form=ss_linear_nlme$pcvrForm, df = ss_linear_nlme$df, boot = 3))+
  ggplot2::labs(x="time")+
  theme(axis.title.y = element_blank())
#* `mgcv`
ss_linear_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=linear_df, type = "mgcv"))
fit_linear_mgcv <- suppressWarnings(fitGrowth(ss_linear_mgcv))
plot_linear_mgcv <- growthPlot(fit=fit_linear_mgcv, form=ss_linear_mgcv$pcvrForm, df = ss_linear_mgcv$df)+
  ggplot2::labs(x="time")+
  theme(axis.title.y = element_blank())
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
  area(16,13,18,15)
  
)
plot(layout)

patch <- plot_logistic_nls + plot_logistic_nlrq + plot_logistic_nlme + 
  plot_logistic_mgcv + plot_logistic_brms +
  plot_gompertz_nls + plot_gompertz_nlrq + plot_gompertz_nlme + 
  plot_gompertz_mgcv + plot_gompertz_brms + 
  plot_mm_nls + plot_mm_nlrq + plot_mm_nlme + 
  plot_monomolecular_mgcv + plot_mm_brms +
  plot_exponential_nls + plot_exponential_nlrq + plot_exponential_nlme +
  plot_exponential_mgcv + plot_exponential_brms + 
  plot_powerlaw_nls + plot_powerlaw_nlrq + plot_powerlaw_nlme +
  plot_powerlaw_mgcv + plot_powerlaw_brms +
  plot_linear_nls + plot_linear_nlrq + plot_linear_nlme + 
  plot_linear_mgcv + plot_linear_brms +
  patchwork::plot_layout(design = layout)

ggsave("~/Desktop/stargate/fahlgren_lab/labMeetings/pcvr_longitudinal_patch.png",
       patch, width= 22, height=15, dpi = 300, bg = "#ffffff")


#* ****************************************
#* ***** `Figure 4: MV Traits Example` *****
#* ****************************************


set.seed(123)

simFreqs<-function(vec, group){
  s1<-hist(vec, breaks=seq(1,181,1), plot=FALSE)$counts
  s1d<-as.data.frame(cbind(data.frame(group), matrix(s1,nrow=1)))
  colnames(s1d)<-c('group', paste0("sim_",1:180))
  s1d
}

sim_df<-rbind(do.call(rbind, lapply(1:10, function(i){ simFreqs(rnorm(200, 70, 10), group="normal") })),
              do.call(rbind, lapply(1:10, function(i){ simFreqs(rlnorm(200, log(30),0.35), group="lognormal") })),
              do.call(rbind, lapply(1:10, function(i){ simFreqs( c(rlnorm(125, log(15),0.25), rnorm(75, 75,5) ), group="bimodal") })),
              do.call(rbind, lapply(1:10, function(i){ simFreqs( c(rlnorm(80, log(20),0.3), rnorm(70, 70,10), rnorm(50, 120,5) ), group="trimodal") })),
              do.call(rbind, lapply(1:10, function(i){ simFreqs( runif(200,1,180), group="uniform") }))
)

sim_df_long<-as.data.frame(data.table::melt(data.table::as.data.table(sim_df), id.vars = "group"))
sim_df_long$bin<-as.numeric(sub("sim_", "", sim_df_long$variable))

sim_plot <- ggplot(sim_df_long, aes(x=bin, y=value, fill=group), alpha=0.25)+
  geom_col(position="identity", show.legend = F)+
  pcv_theme()+
  labs(x="Color Histogram Bin")+
  theme(axis.title.y = element_blank())+
  facet_wrap(~group)



sim_emd<-pcv.emd(df = sim_df, cols="sim_", reorder=c("group"),
                 mat =FALSE, plot=TRUE, parallel = 1, raiseError=TRUE)
emd_plot <- sim_emd$plot + 
  theme(legend.position = "none", axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.line.y.left = element_blank())+
  labs(x = "EMD")

layout <- c(area(1, 1, 2,3), 
            area(2,3,2,3))
# plot(layout)

patch1 <- sim_plot + emd_plot + plot_layout(design = layout)


#emd = sim_emd$data; meta = NULL; dissim=TRUE; distCol="emd"; filter = 0.5; direction="greater"
set.seed(123)
n<-pcv.net(sim_emd$data, filter = 0.5)
net1 <- net.plot(n, fill="group")+
  labs(color = "", title = "Network 1")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_discrete(limits = c("bimodal", "lognormal", "normal",
                                  "trimodal", "uniform"))
set.seed(123)
n<-pcv.net(sim_emd$data, filter = '0.5')
net2 <- net.plot(n, fill="group")+
  labs(color = "", title = "Network 2")+
  theme(axis.line.y.left = element_line(color="gray50", linewidth=1),
        plot.title = element_text(hjust = 0.5), 
        legend.position="none")

layout <- c(area(1, 1, 4,6), 
            area(3,5,4,6),
            area(5,1, 7,3),
            area(5,4,7,6)
            )
# plot(layout)

patch2 <- sim_plot + emd_plot + net1 + net2 + plot_layout(design = layout)

ggsave("~/Desktop/stargate/fahlgren_lab/labMeetings/pcvr_emd_ex.png",
       patch2, width= 12, height=14, dpi = 300, bg = "#ffffff")









#* ***********************************
#* ***** `Excess Longitudinal` *****
#* ***********************************

simdf<-growthSim("logistic", n=20, t=25, params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
l<-ggplot(simdf,aes(time, y, group=interaction(group,id)))+
  geom_line(aes(color=group))+labs(title="Logistic")+theme_minimal()+theme(legend.position="none")

simdf<-growthSim("gompertz", n=20, t=25, params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(0.2, 0.25)))
g<-ggplot(simdf,aes(time, y, group=interaction(group,id)))+
  geom_line(aes(color=group))+labs(title="Gompertz")+theme_minimal()+theme(legend.position="none")

simdf<-growthSim("monomolecular", n=20, t=25, params = list("A"=c(200,160), "B"=c(0.08, 0.1)))
m<-ggplot(simdf,aes(time, y, group=interaction(group,id)))+
  geom_line(aes(color=group))+labs(title="Monomolecular")+theme_minimal()+theme(legend.position="none")

simdf<-growthSim("exponential", n=20, t=25, params = list("A"=c(15, 20), "B"=c(0.095, 0.095)))
e<-ggplot(simdf,aes(time, y, group=interaction(group,id)))+
  geom_line(aes(color=group))+labs(title="Exponential")+theme_minimal()+theme(legend.position="none")

simdf<-growthSim("linear", n=20, t=25, params = list("A"=c(1.1, 0.95)))
ln<-ggplot(simdf,aes(time, y, group=interaction(group,id)))+
  geom_line(aes(color=group))+labs(title="Linear")+theme_minimal()+theme(legend.position="none")

simdf<-growthSim("power law", n=20, t=25, params = list("A"=c(16, 11), "B"=c(0.75, 0.7)))
pl<-ggplot(simdf,aes(time, y, group=interaction(group,id)))+ geom_line(aes(color=group))+
  labs(title="Power Law")+theme_minimal()+theme(legend.position="none")

(l+g+m)/(e+ln+pl)

simdf<-growthSim("double logistic", n=20, t=70,
                 params = list("A"=c(70,50), "B"=c(20, 15), "C"=c(3, 3),
                               "A2"=c(160,210), "B2"=c(45, 45), "C2"=c(4,4)))

dl <- ggplot(simdf,aes(time, y, group=interaction(group,id)))+
  geom_line(aes(color=group))+labs(title="Double Logistic")+
  theme_minimal()+
  theme(legend.position = "none")

simdf<-growthSim("double gompertz", n=20, t=70,
                 params = list("A"=c(70,50), "B"=c(8, 8), "C"=c(0.2, 0.2),
                               "A2"=c(160,210), "B2"=c(35, 40), "C2"=c(0.1, 0.1)))

dg <- ggplot(simdf,aes(time, y, group=interaction(group,id)))+
  geom_line(aes(color=group))+labs(title="Double Gompertz")+
  scale_y_continuous(limits = layer_scales(dl)$y$range$range)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.title.y=element_blank())


simdf<-growthSim("linear + logistic", n=20, t=25,
                 params = list("linear1A"=c(15, 6), "changePoint1"=c(8, 12),
                               "logistic2A"=c(100, 150), "logistic2B"=c(10, 8),
                               "logistic2C"=c(3, 2.5) ))

dd <- ggplot(simdf,aes(time, y, group=interaction(group,id)))+
  geom_line(aes(color=group))+labs(title="GAM/Changepoint*")+
  scale_y_continuous(limits = layer_scales(dl)$y$range$range)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.title.y=element_blank())




(l+g+m)/(e+ln+pl)/(dl + dg + dd)








#* ***** `Double Logistic`
set.seed(123)
dl_df<-growthSim("double logistic", n=20, t=70,
                 params = list("A"=c(70,50), "B"=c(20, 15), "C"=c(3, 3),
                               "A2"=c(160,210), "B2"=c(45, 45), "C2"=c(4,4)))
#* `nls`
ss_dl_nls<-suppressMessages(growthSS(model = "double logistic", form=y~time|id/group, 
                                         df=dl_df, type = "nls"))
ss_dl_nls$start <- list("A"=c(70, 50), "B"=c(20,15), "C"=c(3,3),
                        "A2"=c(160, 210), "B2"=c(45, 45), "C2"=c(4,4))
fit_dl_nls <- fitGrowth(ss_dl_nls)
plot_dl_nls <- growthPlot(fit=fit_dl_nls, form=ss_dl_nls$pcvrForm, df = ss_dl_nls$df) +
  ggplot2::labs(y = "Double Sigmoid")+
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size = 12))
plot_dl_nls
#* `nlrq`
ss_dl_nlrq<-suppressMessages(growthSS(model = "double logistic", form=y~time|id/group,
                                          df=dl_df, type = "nlrq", tau =  seq(0.02, 0.99, 0.02))) # seq(0.02, 0.99, 0.02)
ss_dl_nlrq$start <- list("A"=c(70, 50), "B"=c(20,15), "C"=c(3,3),
                        "A2"=c(160, 210), "B2"=c(45, 45), "C2"=c(4,4))
fit_dl_nlrq <- fitGrowth(ss_dl_nlrq)
plot_dl_nlrq <- growthPlot(fit=fit_dl_nlrq, form=ss_dl_nlrq$pcvrForm, df = ss_dl_nlrq$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank() )

#* `nlme`
ss_dl_nlme<-growthSS(model = "double logistic", form=y~time|id/group, sigma="none",
                         df=dl_df, type = "nlme")
ss_dl_nlme$start <- c("A"=c(70, 50), "B"=c(20,15), "C"=c(3,3),
                        "A2"=c(160, 210), "B2"=c(45, 45), "C2"=c(4,4))
fit_dl_nlme <- suppressWarnings(fitGrowth(ss_dl_nlme))
plot_linear_nlme <- suppressWarnings(growthPlot(fit=fit_dl_nlme, form=ss_dl_nlme$pcvrForm, df = ss_dl_nlme$df, boot = 3))+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank() )
#* `mgcv`
ss_dl_mgcv<-suppressMessages(growthSS(model = "gam", form=y~time|id/group, df=dl_df, type = "mgcv"))
fit_dl_mgcv <- suppressWarnings(fitGrowth(ss_dl_mgcv))
plot_dl_mgcv <- growthPlot(fit=fit_dl_mgcv, form=ss_dl_mgcv$pcvrForm, df = ss_dl_mgcv$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank() )
#* `brms`
ss_dl_brms<-growthSS(model = "double logistic", form=y~time|id/group, sigma="gompertz",
                         list('A' = 30, "B"=15, "C"=3,
                              "A2"= 150, "B2"=30, "C2"=3,
                              "subA" = 20, "subB"= 30, "subC"=2),
                         df=dl_df, type = "brms")
fit_dl_brms <- fitGrowth(ss_dl_brms, backend="cmdstanr", iter=1000, chains=1, cores=1)
plot_dl_brms <- growthPlot(fit=fit_dl_brms, form=ss_dl_brms$pcvrForm, df = ss_dl_brms$df)+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank() )




#* expanded list of growth models


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
  plot_dl_nls + plot_dl_nlrq + plot_dl_nlme + 
  plot_dl_mgcv + plot_dl_brms + 
  plot_mm_nls + plot_mm_nlrq + plot_mm_nlme + 
  plot_monomolecular_mgcv + plot_mm_brms +
  plot_exponential_nls + plot_exponential_nlrq + plot_exponential_nlme +
  plot_exponential_mgcv + plot_exponential_brms + 
  plot_powerlaw_nls + plot_powerlaw_nlrq + plot_powerlaw_nlme +
  plot_powerlaw_mgcv + plot_powerlaw_brms +
  plot_linear_nls + plot_linear_nlrq + plot_linear_nlme + 
  plot_linear_mgcv + plot_linear_brms +
  patchwork::plot_layout(design = layout)

ggsave("~/Desktop/stargate/fahlgren_lab/labMeetings/pcvr_expanded_longitudinal_patch.png",
       patch, width= 22, height=15, dpi = 300, bg = "#ffffff")
ls()

paste0(ls(pattern = "fit_|ss_|plot_|_df"), collapse=", ")

save(exponential_df, fit_exponential_brms, fit_exponential_mgcv, fit_exponential_nlme, 
     fit_exponential_nlrq, fit_exponential_nls, fit_gompertz_brms, fit_gompertz_mgcv, 
     fit_gompertz_nlme, fit_gompertz_nlrq, fit_gompertz_nls, fit_linear_brms, fit_linear_mgcv,
     fit_linear_nlme, fit_linear_nlrq, fit_linear_nls, fit_logistic_brms, fit_logistic_mgcv,
     fit_logistic_nlme, fit_logistic_nlrq, fit_logistic_nls, fit_mm_brms, fit_mm_nlme, 
     fit_mm_nlrq, fit_mm_nls, fit_monomolecular_mgcv, fit_powerlaw_brms, fit_powerlaw_mgcv, 
     fit_powerlaw_nlme, fit_powerlaw_nlrq, fit_powerlaw_nls, gompertz_df, linear_df,
     logistic_df, mm_df, plot_exponential_brms, plot_exponential_mgcv, plot_exponential_nlme, 
     plot_exponential_nlrq, plot_exponential_nls, plot_gompertz_brms, plot_gompertz_mgcv, 
     plot_gompertz_nlme, plot_gompertz_nlrq, plot_gompertz_nls, plot_linear_brms, 
     plot_linear_mgcv, plot_linear_nlme, plot_linear_nlrq, plot_linear_nls, plot_logistic_brms, 
     plot_logistic_mgcv, plot_logistic_nlme, plot_logistic_nlrq, plot_logistic_nls, 
     plot_mm_brms, plot_mm_nlme, plot_mm_nlrq, plot_mm_nls, plot_monomolecular_mgcv, 
     plot_powerlaw_brms, plot_powerlaw_mgcv, plot_powerlaw_nlme, plot_powerlaw_nlrq, 
     plot_powerlaw_nls, powerlaw_df, ss_exponential_brms, ss_exponential_mgcv, 
     ss_exponential_nlme, ss_exponential_nlrq, ss_exponential_nls, ss_gompertz_brms, 
     ss_gompertz_mgcv, ss_gompertz_nlme, ss_gompertz_nlrq, ss_gompertz_nls, ss_linear_brms, 
     ss_linear_mgcv, ss_linear_nlme, ss_linear_nlrq, ss_linear_nls, ss_logistic_brms,
     ss_logistic_mgcv, ss_logistic_nlme, ss_logistic_nlrq, ss_logistic_nls, ss_mm_brms, 
     ss_mm_nlme, ss_mm_nlrq, ss_mm_nls, ss_monomolecular_mgcv, ss_powerlaw_brms, ss_powerlaw_mgcv,
     ss_powerlaw_nlme, ss_powerlaw_nlrq, ss_powerlaw_nls,
     file = "~/Desktop/stargate/fahlgren_lab/labMeetings/labMeeting_10-19-2023.rdata")








