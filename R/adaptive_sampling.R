#' Adaptive sampling for models specified in pcvrss objects.


if (FALSE) {
  message("loaded adaptive sampling values")
  devtools::load_all("~/pcvr")
  library(brms)
  library(ggplot2)
  simdf <- growthSim("gompertz",
                     n = 20, t = 25,
                     params = list("A" = c(100), "B" = c(9), "C" = c(0.25))
  )
  ss <- growthSS(
    model = "gompertz", form = y ~ time | id / group,
    sigma = "int", df = simdf,
    start = list("A" = 130, "B" = 12, "C" = 3), type = "brms"
  )
  full_fit <- fitGrowth(ss, cores = 4, chains = 4, iter = 2000)
  ss <- growthSS(
    model = "gompertz", form = y ~ time | id / group,
    sigma = "int", df = simdf[simdf$time < 7, ],
    start = list("A" = 130, "B" = 12, "C" = 3), type = "brms"
  )
  (proposed_sampling_times <- unique(simdf$time[!simdf$time %in% ss$df$time]))
  threshold <- 1 # threshold for variance of Posterior of Expected Value to take a new sample.
  #* eventually I'd like this to take something more like a brms hypothesis string.
  #* sampling_rule <- "sigma_e > 10"
  #* additional_quantities <- list("sigma_e" = sd(e_mu), # idk, maybe it works something like this.
  #*                               "e_mu" = brms::posterior_epred(.x))
  fit <- NULL


pdf("~/Desktop/adaptive_design_paper/thresholding_sigma.pdf", width = 8, height = 6.5)
res <- lapply(1:30, function(u) {
  fit <- NULL
#adaptive_sampling <- function(ss, proposed_sampling_times, threshold = 1, fit = NULL, ...) {
  #* note, there should be `...` arguments passed to update.brmsfit I think.
  #* might include a `simplify` argument here that would remove any distributional models not
  #*    required for the threshold rule?
  if (is.null(fit)) {
    fit <- fitGrowth(ss, cores = 4, chains = 4, iter = 2000) # ellipses are passed here too I think
  }
  #* this will need to be parsed from the pcvrForm to take into account variable names.
  group <- "group"
  individual <- "id"
  x_var <- "time"
  y_var <- "y"
  possible_sampling_times <- c(8:25)
  predict_at <- ss$df
  threshold <- 0.6
  variances <- list()
  means <- list()
  accepted_sampling_times <- list()
  for (t in possible_sampling_times) {
    data_to_new_t <- do.call(rbind, lapply(unique(predict_at[[group]]), function(grp) {
      sub <- predict_at[predict_at[[group]] == grp, ]
      #* for each ID I want to add another timepoint
      sub <- do.call(rbind, lapply(unique(sub[[individual]]), function(id) {
        id_df <- sub[sub[[individual]] == id, ]
        id_df <- rbind(id_df, id_df[nrow(id_df), ])
        id_df[nrow(id_df), x_var] <- t
        return(id_df)
      }))
      p_iter <- posterior_predict(fit, sub)
      sub[sub[[x_var]] == t, y_var] <- apply(p_iter[, which(sub[[x_var]] == t)], 2, sample, 1)
      return(sub)
    }))
    # data_to_new_t[data_to_new_t$time %in% c(10, 11), ]
    data_to_new_t #* here I need to be making the new data to fit the model to. This is where the
    #* raw MCMC was a little easier to work with I think.
    fit_iter <- update(fit, newdata = data_to_new_t, cores = 4, chains = 4, iter = 2000) #* pass `...`
    post_pred_t <- posterior_epred(fit_iter)
    E_y_epred_t <- mean(post_pred_t[, which(data_to_new_t[[x_var]] == t)[1]])
    # Var_y_pred_t <- sapply(t, pred_var, theta = post)
    #* okay this isn't working because the epreds are not different rep to rep. It's at group level.
    Var_y_epred_t <- var(post_pred_t[, which(data_to_new_t[[x_var]] == t)[1]])
    means <- append(means, E_y_epred_t)
    variances <- append(variances, Var_y_epred_t)
    if (Var_y_epred_t > threshold) {
      predict_at <- data_to_new_t
      accepted_sampling_times <- append(accepted_sampling_times, t)
      fit <- fit_iter
    }
  }
  unlist(accepted_sampling_times)
  unlist(means)
  unlist(variances)
  
#}

# data.frame(
#   mu = unlist(means),
#   var = unlist(variances),
#   t = possible_sampling_times
#   )
# 
# apply(as.data.frame(fit), 2, summary)
# apply(as.data.frame(full_fit), 2, summary)


# growthPlot(full_fit, form = ss$pcvrForm, df = simdf)

p1 <- growthPlot(fit, form = ss$pcvrForm, timeRange = 1:25) +
  lapply(accepted_sampling_times, function(t) {
    geom_vline(xintercept = t, linewidth = 0.25, linetype = 3)
  }) +
  geom_line(data = simdf,
            aes(time, y, group = interaction(group, id)),
            color = "gray60", linewidth = 0.25,
            inherit.aes = FALSE) +
  geom_point(data = predict_at[predict_at$time > max(ss$df$time), ],
             aes(time, y), color = "black",
             size = 1,
             inherit.aes = FALSE) +
  pcv_theme()

# ggsave("~/Desktop/adaptive_design_paper/thresholding_sigma.png", width = 8, height = 6,
#        dpi = 300, bg = "#ffffff")

print(p1)
return(unlist(accepted_sampling_times))
})
dev.off()

table(unlist(res)) # histogram of 30 trials

ggplot(data.frame(x = unlist(res))) +
  geom_histogram(aes(x = x), binwidth = 1,color = "white")
ggsave("~/Desktop/adaptive_design_paper/sampling_time_frequencies_30_trials.png", width = 8, height = 6,
       dpi = 300, bg = "#ffffff")


#* derivative based version
#*
#* Calculate derivative with respect to each model parameter
#* Fit model with 7 days of data
#* For each parameter of interest
#*    Draw values from posterior
#*    Use values from posterior to calculate Information over remaining X
#*    calculate Remaining_Information as 1 - Cumul.Sum(Info_x) / sum(cumul.sum(Info)) * 
#*    * ie, (amount of change left to happen)
#*       If Remaining_Information > Threshold
#*          return Remaining_Information
#*       else
#*          pass to next parameter
#* If any parameters of interest are < Threshold of Current_Information (> threshold of remaining info)
#*    sum Remaining_Information for each parameter of interest
#*    sample proposal time T from remaining X weighted by parameter of interest (could just pick mode?)
#*       (for simulation purposes, add data from time T)
#*    Re-fit model
#*    Repeat.
#* return model
#* 

devtools::load_all("~/pcvr")
library(brms)
library(ggplot2)

derivs <- function(x, form = expression(a * exp(-b * exp(-c*x))), wrt = "x", scale = TRUE, ...) {
  y <- with(list(x = x, wrt = wrt, form = form, ...), {
    eval(D(form, wrt))
  })
  if (scale) {
    y <- y / sum(y)
  }
  return(y)
}

required_information_content <- 0.8 # get 80% of the information possible
parameters_of_interest <- c("a", "b", "c")

for (i in 1) { # test N times
  simdf <- growthSim("gompertz",
                     n = 20, t = 25,
                     params = list("A" = c(100), "B" = c(9), "C" = c(0.25))
  )
  ss_full <- growthSS(
    model = "gompertz", form = y ~ time | id / group,
    sigma = "int", df = simdf,
    start = list("A" = 130, "B" = 12, "C" = 3), type = "brms"
  )
  full_fit <- fitGrowth(ss_full, cores = 4, chains = 4, iter = 2000)
  ss <- growthSS(
    model = "gompertz", form = y ~ time | id / group,
    sigma = "int", df = simdf[simdf$time < 7, ],
    start = list("A" = 130, "B" = 12, "C" = 3), type = "brms"
  )
  proposed_sampling_times <- unique(simdf$time[!simdf$time %in% ss$df$time])
  sampled_times <- unique(ss$df$time)
  required_information_content <- c(a = 0.8, b = 0.8, c = 0.8) # get 80% of the information possible
  fit <- fitGrowth(ss, cores = 4, chains = 4, iter = 2000)
  
  collected_information <- c(a = 0, b = 0, c = 0)
  ti <- 1
  t <- proposed_sampling_times[ti]
  while (any(collected_information < required_information_content) && t < max(proposed_sampling_times)) {
    fitd <- setNames(as.data.frame(fit)[,3:5], c("a", "b", "c"))
    n_draws <- 100
    draw_index <- sample(seq_len(nrow(fitd)), n_draws)
    draws <- fitd[draw_index, ]
    primes <- as.data.frame(do.call(cbind, lapply(parameters_of_interest, function(par) {
      prime <- derivs(
        x = rep(c(sampled_times, proposed_sampling_times), each = n_draws),
        a = rep(draws$a, times = length(c(sampled_times, proposed_sampling_times))),
        b = rep(draws$b, times = length(c(sampled_times, proposed_sampling_times))),
        c = rep(draws$c, times = length(c(sampled_times, proposed_sampling_times))),
        wrt = par
      )
      return(prime)
    })))
    colnames(primes) <- parameters_of_interest
    primes$time <- rep(c(sampled_times, proposed_sampling_times), each = n_draws)
    primes$i <- rep(seq_len(n_draws), times = length(c(sampled_times, proposed_sampling_times)))
    # # gut check
    ggplot(primes, aes(x = time, group = i)) +
      geom_line(aes(y = a, color = "a")) +
      geom_line(aes(y = b, color = "b")) +
      geom_line(aes(y = c, color = "c")) +
      geom_vline(xintercept = t)
    mu_time <- aggregate(cbind(a, b, c) ~ time, primes, \(x) {mean(x)})
    cs_time <- apply(mu_time[, 2:4], 2, cumsum)
    info <- as.data.frame(cbind(ag$time, apply(cs_time, 2, \(x) {x / max(x)})))
    collected_information <- info[info[[1]] == t, parameters_of_interest]
    ti <- ti + 1
  }
}

#* derivative based version
#*
#* Calculate derivative with respect to each model parameter
#* Fit model with 7 days of data
#* For each parameter of interest
#*    Draw values from posterior
#*    Use values from posterior to calculate Information over remaining X
#*    calculate Remaining_Information as 1 - Cumul.Sum(Info_x) / sum(cumul.sum(Info)) * 
#*    * ie, (amount of change left to happen)
#*       If Remaining_Information > Threshold
#*          return Remaining_Information
#*       else
#*          pass to next parameter
#* If any parameters of interest are < Threshold of Current_Information (> threshold of remaining info)
#*    sum Remaining_Information for each parameter of interest
#*    sample proposal time T from remaining X weighted by parameter of interest (could just pick mode?)
#*       (for simulation purposes, add data from time T)
#*    Re-fit model
#*    Repeat.
#* return model
#* 















}