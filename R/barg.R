#' Function to help fulfill elements of the Bayesian Analysis Reporting Guidelines.
#'
#' The Bayesian Analysis Reporting Guidelines were put forward by Kruschke
#' (https://www.nature.com/articles/s41562-021-01177-7) to aide in reproducibility and documentation
#' of Bayesian statistical analyses that are sometimes unfamiliar to reviewers or scientists.
#' The purpose of this function is to summarize goodness of fit metrics from one or more Bayesian models
#' made by \link{growthSS} and \link{fitGrowth}. See details for explanations of those metrics and
#' the output.
#'
#' @param fit The brmsfit object or a list of brmsfit objects in the case that you split models to run
#' on subsets of the data for computational simplicity.
#' @param ss The growthSS output used to specify the model. If fit is a list then this can either be one
#' growthSS list in which case the priors are assumed to be the same for each model or it can be a list
#' of the same length as fit. Note that the only parts of this which are used are the \code{call$start}
#' which is expected to be a call, \code{pcvrForm}, and \code{df} list elements,
#' so if you have a list of brmsfit objects and no ss object you can specify a stand-in list. This
#' can also be left NULL (the default) and posterior predictive plots and prior predictive plots will
#' not be made.
#'
#'
#' @details
#'
#'
#' \itemize{
#'     \item \bold{General}: This includes chain number, length, and total divergent transitions per
#'     model. Divergent transitions are a marker that the MCMC had something go wrong.
#'     Conceptually it may be helpful to think about rolling a marble over a 3D curve then having the
#'     marble suddenly jolt in an unexpected direction, something happened that suggests a
#'     problem/misunderstood surface. In practice you want extremely few (ideally no) divergences.
#'     If you do have divergences then consider specifying more control parameters
#'     (see brms::brm or examples for \link{fitGrowth}). If the problem persists then the model may need
#'     to be simplified. For more information on MCMC and divergence see the stan
#'     manual (https://mc-stan.org/docs/2_19/reference-manual/divergent-transitions).
#'
#'     \item \bold{ESS}: ESS stands for Effective Sample Size and is a goodness of fit metric that
#'     approximates the number of independent replicates that would equate to the same amount of
#'     information as the (autocorrelated) MCMC iterations. ESS of 1000+ is often considered as a pretty
#'     stable value, but more is better. Still, 100 per chain may be plenty depending on your
#'     applications and the inference you wish to do. One of the benefits to using lots of chains and/or
#'     longer chains is that you will get more complete information and that benefit will be shown by a
#'     larger ESS. This is separated into "bulk" and "tail" to represent the middle and tails of the
#'     posterior distribution, since those can sometimes have very different sampling behavior.
#'     A summary and the total values are returned, with the summary being useful if several models are
#'     included in a list for fit argument
#'
#'     \item \bold{Rhat}: Rhat is a measure of "chain mixture". It compares the between vs within chain
#'     values to assess how well the chains mixed. If chains did not mix well then Rhat will be greater
#'     than 1, with 1.05 being a broadly
#'     agreed upon cutoff to signify a problem. Running longer chains should result in lower Rhat
#'     values. The default in brms is to run 4 chains, partially to ensure that there is a good chance
#'     to check that the chains mixed well via Rhat. A summary and the total values are returned, with
#'     the summary being useful if several models are included in a list for fit argument
#'
#'     \item \bold{NEFF}: NEFF is the NEFF ratio (Effective Sample Size over Total MCMC Sample Size).
#'     Values greater than 0.5 are generally considered good, but there is a consensus that lower can be
#'     fine down to about 0.1. A summary and the total values are returned, with the summary being
#'     useful if several models are included in a list for fit argument
#'
#'     \item \bold{mcmcTrace}: A plot of each model's MCMC traces. Ideally these should be very mixed
#'     and stationary. For more options for visualizing MCMC diagnostics see
#'     \code{bayesplot::mcmc_trace}.
#'
#'     \item \bold{priorPredictive}: A plot of data simulated from the prior using \link{plotPrior}.
#'     This should generate data that is biologically plausible for your situation, but it will
#'     probably be much more variable than your data. That is the effect of the mildly informative thick
#'     tailed lognormal priors. If you specified non-default style priors then this currently will not
#'     work.
#'
#'     \item \bold{posteriorPredictive}: A plot of each model's posterior predictive interval over time.
#'     This is the same as plots returned from \link{growthPlot} and shows 1-99% intervals in purple
#'     coming to a mean yellow trend line. These should encompass the overwhelming majority of your data
#'     and ideally match the variance pattern that you see in your data. If parts of the predicted
#'     interval are biologically impossible (area below 0, percentage about 100%, etc) then your chosen
#'     model should be reconsidered.
#'     }
#'
#'
#'
#'
#' @keywords Bayesian brms prior
#' @return A named list containing Rhat, ESS, NEFF, and Trace/Prior/Posterior Predictive plots.
#' See details for interpretation.
#' @importFrom rlang is_installed
#' @seealso \link{plotPrior} for visual prior predictive checks.
#' @examples
#' \donttest{
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group, sigma = "logistic",
#'   df = simdf, start = list(
#'     "A" = 130, "B" = 12, "C" = 3,
#'     "sigmaA" = 20, "sigmaB" = 10, "sigmaC" = 2
#'   ), type = "brms"
#' )
#' fit_test <- fitGrowth(ss,
#'   iter = 600, cores = 1, chains = 1, backend = "cmdstanr",
#'   sample_prior = "only" # only sampling from prior for speed
#' )
#' barg(fit_test, ss)
#' fit_2 <- fit_test
#' fit_list <- list(fit_test, fit_2)
#' x <- barg(fit_list, list(ss, ss))
#' }
#'
#' @export

barg <- function(fit, ss = NULL) {
  out <- list()
  #* `format everything into lists`
  if (methods::is(fit, "brmsfit")) {
    fitList <- list(fit)
  } else {
    fitList <- fit
  }
  if (!is.null(names(ss)) && names(ss)[1] == "formula") {
    ssList <- lapply(seq_along(fitList), function(i) {
      return(ss)
    })
  } else {
    ssList <- ss
  }
  #* `General Info`
  general <- do.call(rbind, lapply(seq_along(fitList), function(i) {
    fitobj <- fitList[[i]]
    ms <- summary(fitobj)
    df <- data.frame(
      chains = ms$chains,
      iter = ms$iter,
      num.divergent = rstan::get_num_divergent(fitobj$fit),
      model = i
    )
    return(df)
  }))
  out[["General"]] <- general
  #* `Rhat summary`
  rhats <- as.data.frame(do.call(rbind, lapply(fitList, function(fitobj) {
    return(brms::rhat(fitobj))
  })))
  rhat_metrics <- apply(rhats, MARGIN = 2, summary)
  rhats$model <- seq_along(fitList)
  out[["Rhat"]][["summary"]] <- rhat_metrics
  out[["Rhat"]][["complete"]] <- rhats
  #* `NEFF summary`
  neff <- as.data.frame(do.call(rbind, lapply(fitList, function(fitobj) {
    return(brms::neff_ratio(fitobj))
  })))
  neff_metrics <- apply(neff, MARGIN = 2, summary)
  neff$model <- seq_along(fitList)
  out[["NEFF"]][["summary"]] <- neff_metrics
  out[["NEFF"]][["complete"]] <- neff
  #* `ESS Summary`
  ess <- do.call(rbind, lapply(seq_along(fitList), function(i) {
    fit <- fitList[[i]]
    ms <- summary(fit)
    df <- data.frame(
      "par" = c(rownames(ms$fixed), rownames(ms$spec_pars)),
      "Bulk_ESS" = c(ms$fixed$Bulk_ESS, ms$spec_pars$Bulk_ESS),
      "Tail_ESS" = c(ms$fixed$Tail_ESS, ms$spec_pars$Tail_ESS),
      "model" = i
    )
    return(df)
  }))
  ag_b_ess <- aggregate(Bulk_ESS ~ par, ess, summary)
  tag_b_ess <- t(ag_b_ess[-1])
  colnames(tag_b_ess) <- ag_b_ess[, 1]
  ag_t_ess <- aggregate(Tail_ESS ~ par, ess, summary)
  tag_t_ess <- t(ag_t_ess[-1])
  colnames(tag_t_ess) <- ag_t_ess[, 1]
  ess_metrics <- rbind(tag_b_ess, tag_t_ess)
  out[["ESS"]][["summary"]] <- ess_metrics
  out[["ESS"]][["complete"]] <- ess
  #* `MCMC Diagnostic Plot`
  tracePlots <- lapply(fitList, function(fit) {
    p <- suppressMessages(brms::mcmc_plot(fit, type = "trace"))
    return(p)
  })
  out[["mcmcTrace"]] <- tracePlots
  #* `Prior Predictive Check`
  if (methods::is(eval(ssList[[1]]$call$start), "list")) {
    pri_preds <- lapply(seq_along(ssList), function(i) {
      ss <- ssList[[i]]
      x <- trimws(gsub("[|].*|[/].*", "", as.character(ss$pcvrForm)[3]))
      t <- max(ss$df[[x]])
      pri_pred <- plotPrior(priors = eval(ss$call$start), type = ss$model, n = 200, t = t)
      return(pri_pred$simulated)
    })
    out[["priorPredictive"]] <- pri_preds
  }
  #* `Posterior Predictive Check`
  if (methods::is(eval(ssList[[1]]$pcvrForm), "formula")) {
    post_preds <- lapply(seq_along(fitList), function(i) {
      return(brmPlot(fitList[[i]], form = ssList[[i]]$pcvrForm, df = ssList[[i]]$df))
    })
    out[["posteriorPredictive"]] <- post_preds
  }
  return(out)
}
