#' Function to help fulfill elements of the Bayesian Analysis Reporting Guidelines.
#' 
#' The Bayesian Analysis Reporting Guidelines were put forward by Kruschke (https://www.nature.com/articles/s41562-021-01177-7) 
#' to aide in reproducibility and documentation of Bayesian statistical analyses that are sometimes unfamiliar to reviewers or scientists.
#' The purpose of this function is to summarize goodness of fit metrics from one or more Bayesian models made by
#' \link{growthSS} and \link{fitGrowth}
#' 
#' @param fit The brmsfit object or a list of brmsfit objects in the case that you split models to run on subsets of the data for computational simplicity.
#' @param ss The growthSS output used to specify the model. If fit is a list then this can either be one growthSS list in which case the priors are assumed
#' to be the same for each model or it can be a list of the same length as fit.
#' 
#' @details
#' 
#' More details about ESS, Rhat, and NEFF ratio here
#' 
#' 
#' @keywords Bayesian, brms, prior
#' @return A named list containing Rhat, ESS, NEFF, and Prior/Posterior Predictive plots.
#' See details for more about Rhat, ESS, and NEFF.
#' 
#' @examples 
#' 
#' ## Not run:
#'
#' if(FALSE){
#' 
#' }
#' 
#' ## End(Not run:)
#'
#' @export

barg <- function(fit, ss){
  out <- list()
  #* `format everything into lists`
  if(methods::is(fit, "brmsfit")){
    fitList <- list(fit)
  } else{
    fitList <- fit
  }
  if(names(ss)[1]=="formula"){
    singleSS =TRUE
    ssList <- lapply(1:length(fitList), function(i) {ss})
  } else {
    singleSS = FALSE
    ssList <- ss
  }
  #* `General Info`
  general <- do.call(rbind, lapply(1:length(fitList), function(i){
    fitobj <- fitList[[i]]
    ms <- summary(fitobj)
    data.frame(chains = ms$chains,
         iter = ms$iter,
         num.divergent = rstan::get_num_divergent(fitobj$fit),
         model=i)
  }))
  #* `Rhat summary`
  rhats <- do.call(rbind, lapply(fitList, function(fitobj){
    brms::rhat(fitobj)
  }))
  rhat_metrics <- apply(rhats, MARGIN = 2, summary)
  out[["Rhat"]][["summary"]] <- rhat_metrics
  out[["Rhat"]][["complete"]] <- rhats
  #* `NEFF summary`
  neff <- do.call(rbind, lapply(fitList, function(fitobj){
    brms::neff_ratio(get(fitobj))
  }))
  neff_metrics <- apply(neff, MARGIN = 2, summary)
  out[["NEFF"]][["summary"]] <- neff_metrics
  out[["NEFF"]][["complete"]] <- neff
  #* `ESS Summary`
  # files <- dir("~/Desktop/model_fit_R_objects", full.names = TRUE)
  # for(file in files){load(file)}
  # fitList <- lapply(ls(pattern="_fit"), get)
  ess <- do.call(rbind, lapply(1:length(fitList), function(i){
    fit <- fitList[[i]]
    ms <- summary(fit)
    data.frame("par" = c(rownames(ms$fixed), rownames(ms$spec_pars)), 
               "Bulk_ESS" = c(ms$fixed$Bulk_ESS, ms$spec_pars$Bulk_ESS),
               "Tail_ESS" = c(ms$fixed$Tail_ESS, ms$spec_pars$Tail_ESS),
               "model" = i)
  }))
  ag_b_ess <- aggregate(Bulk_ESS ~ par, ess, summary)
  tag_b_ess <- t(ag_b_ess[-1])
  colnames(tag_b_ess)<-ag_b_ess[,1]
  ag_t_ess <- aggregate(Tail_ESS ~ par, ess, summary)
  tag_t_ess <- t(ag_t_ess[-1])
  colnames(tag_t_ess)<-ag_t_ess[,1]
  ess_metrics <- rbind(tag_b_ess, tag_t_ess)
  out[["ESS"]][["summary"]] <- ess_metrics
  out[["ESS"]][["complete"]] <- ess
  #* `Prior Predictive Check`
  pri_preds <- lapply(1:length(ssList), function(i){
    ss <- ssList[[i]]
    x <- trimws(gsub("[|].*|[/].*","",as.character(ss$pcvrForm)[3]))
    t <- max(ss$df[[x]])
    pri_pred <- plotPrior(priors = ss$call$start, type = ss$model, n = 200, t = t)
    pri_pred$simulated
  })
  out[["priorPredictive"]] <- pri_preds
  #* `Posterior Predictive Check`
  post_preds <- lapply(1:length(fitList), function(i){
    post_pred <- brmPlot(fitList[[i]], form = ssList[[i]]$pcvrForm, df = ssList[[i]]$df)
  })
  out[["posteriorPredictive"]] <- post_preds
  
  return(out)
}
