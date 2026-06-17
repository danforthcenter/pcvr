setwd("~/pcvr/data/")
devtools::load_all("~/pcvr")
set.seed(123)
df <- pcvr::growthSim("logistic", n = 10, t = 20, params = list(A = 100, B = 10, C = 3))
ss <- pcvr::growthSS("logistic", y ~ time | id / group, sigma = "logistic", df = df,
  start = list(A = 100, B = 10, C = 3, sigmaA = 10, sigmaB = 10, sigmaC = 3))
fit <- pcvr::fitGrowth(ss, backend = "cmdstanr", iter = 4100, warmup = 4000, chains = 2, cores = 2)
names(fit$fit)
object.size(fit) #165128 bytes
object.size(ss) # 35368 bytes
save(fit, ss, file = "fit.rda")
# might want a survival one?
set.seed(123)
df <- growthSim("exponential",
  n = 10, t = 30,
  params = list("A" = c(1, 1), "B" = c(0.15, 0.2))
)
survss <- growthSS(
  model = "survival weibull", form = y > 25 ~ time | id / group,
  df = df, start = c(0, 5)
)
surv <- fitGrowth(survss, iter = 4100, warmup = 4000, cores = 2, chains = 2, backend = "cmdstanr")
brmSurvPlot(surv, form = survss$pcvrForm, df = survss$df)
save(surv, survss, file = "surv.rda")
