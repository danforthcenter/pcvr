#' Function to visualize posterior predictive distributions from brms models
#' similar to those made using growthSS outputs.
#' 
#' Models fit using \link{growthSS} inputs by \link{fitGrowth} (and similar models made through other means)
#'  can be visualized easily using this function. This function tests hypotheses in the original scale
#'  (not by growth curve parameters), visualizes the difference in predictions,
#'   and will generally be less powerful than brms::hypothesis used on growth model parameters.
#'  
#' @param fit A brmsfit object, similar to those fit with \code{\link{growthSS}} outputs. Alternatively this can be a list 
#' of two brmsfit objects, in which case they are compared. If this option is used then the groups argument should be specified.
#' @param form A formula similar to that in \code{growthSS} inputs specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' If the fit argument is a list of two brmsfit objects then this can be a list or a single form will be recycled to length 2.
#' @param groups An optional set of groups to use.
#' Defaults to NULL in which case the first two groups in the model are used. If the fit argument is a list of two
#' brmsfit objects then the first group comes from the first model and the second group comes from the second model,
#' in which case a NULL value here will use the first group from each model.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the available data has not reached some point (such as asymptotic size),
#' although prediction using splines outside of the observed range is not necessarily reliable.
#' @param hyp Optionally a character string representing a hypothesis to test. To test the difference of 
#' two groups "diff" can be used as shorthand for group1-group2. If fit is a list then group names should be prefixed with 
#' "fit1_" or "fit2_", respectively.
#' @param plot Logical, should a ggplot be returned?
#' @keywords growth-curve, brms
#' @import ggplot2
#' @import viridis
#' @importFrom stats as.formula setNames
#' @importFrom methods is
#' @examples 
#' 
#' ## Not run:
#' 
#' if(FALSE){
#' df <- growthSim("logistic", n=20, t=25,
#'   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "logistic", form=y~time|id/group,
#'             sigma="logistic", df=df,
#'             start = list("A"=130, "B"=12, "C"=3,
#'             "subA"=20, "subB"=10, "subC"=1), type="brms")
#' fit <- fitGrowth(ss, iter = 1000, cores = 2, chains = 2, backend = "cmdstanr",
#'                 control = list(adapt_delta = 0.999, max_treedepth = 20))
#' 
#' x<-postPred(fit, form=ss$pcvrForm, groups=NULL, timeRange=NULL, hyp=NULL, plot=TRUE)
#' names(x)
#' x$plot
#' }
#' 
#' ## End(Not run)
#' 
#' @return Returns a named list including a dataframe of draws from the posterior predictive distribution
#' for the specified groups, a plot if plot=TRUE and a dataframe of posterior probabilities
#' from brms::hypothesis if a hypothesis was specified.
#' 
#' @export



postPred <- function(fit, form, groups=NULL, timeRange=NULL, hyp = NULL, plot=TRUE){
  
  if(!methods::is(fit, "brmsfit")){
    fit1 <- fit[[1]]
    fit2 <- fit[[2]]
    multiModel <- TRUE
  } else{
    fit1 <- fit
    multiModel <- FALSE
  }
  
  if(multiModel){
    if(is.character(form) | methods::is(form, "formula")){
      form <- list(form, form)
    }
    
    form1 <- .parseFormula(fit1, form[[1]] )
    form2 <- .parseFormula(fit2, form[[2]] )
    
    fitData1 <-form1[[1]]
    y1 <- form1[[2]]
    x1 <- form1[[3]]
    group1 <- form1[[4]]

    fitData2 <-form2[[1]]
    y2 <- form2[[2]]
    x2 <- form2[[3]]
    group2 <- form2[[4]]
    
    if(is.null(groups)){
      groups = c(fitData1[[group1]][1], fitData2[[group2]][1])
    }
    
    if(is.null(timeRange)){
      times <- union(unique(fitData1[[x1]]), unique(fitData2[[x2]]))
    } else {
      times = timeRange
    }
    ndf1test <- stats::setNames(data.frame(min(times), groups[1]), c(x1, group1))
    ndraws1 <- nrow(predict(fit1, ndf1test, summary=FALSE))
    
    ndf2test <- stats::setNames(data.frame(min(times), groups[2]), c(x2, group2))
    ndraws2 <- nrow(predict(fit2, ndf2test, summary=FALSE))
    
    nDraws <- min(c(ndraws1, ndraws2))
    
    ppred <-do.call(rbind, lapply(times, function(i){
      ndf1 <- stats::setNames(data.frame(i, groups[1]), c("time", "group"))
      ppred1 <- as.data.frame(predict(fit1, ndf1, summary=FALSE, ndraws = nDraws))
      colnames(ppred1)<-paste0("fit1_",groups[1])
      ndf2 <- stats::setNames(data.frame(i, groups[2]), c("time", "group"))
      ppred2 <- as.data.frame(predict(fit2, ndf2, summary=FALSE, ndraws = nDraws))
      colnames(ppred2)<-paste0("fit2_",groups[2])
      tpred <- cbind(ppred1, ppred2)
      tpred[["time"]]<-i
      tpred
    }))
    ppred$diff <- ppred[[1]] - ppred[[2]]
    
  } else { # all using one model
    
    form1 <- .parseFormula(fit1, form )
    
    fitData <-form1[[1]]
    y <- form1[[2]]
    x <- form1[[3]]
    group <- form1[[4]]
    
    if(is.null(groups)){
      groups = unique(fitData[[group]])[1:2]
    }
    
    if(is.null(timeRange)){
      times <- unique(fitData[[x]])
    } else {
      times = timeRange
    }
    
    ppred <-do.call(rbind, lapply(times, function(i){
      tpred <- do.call(cbind, lapply(groups, function(g){
        ndf <- stats::setNames(data.frame(i, g), c("time", "group"))
        ppred <- as.data.frame(predict(fit, ndf, summary=FALSE))
        colnames(ppred)<-g
        ppred
      }))
      tpred[["time"]]<-i
      tpred
    }))
    ppred$diff <- ppred[[groups[1] ]] - ppred[[groups[2] ]]
  }
  
  out <- list()
  
  
  out[["posterior_predictive"]] <- ppred
  
  probs <- seq(from=99, to=1, by=-2)/100
  prob_labs <- seq(from=99, to=1, by=-2)
  
  qpreds <- as.data.frame(do.call(rbind, lapply(unique(ppred$time), function(i){
    x <- ppred[ppred$time == i, "diff"]
    q <- as.numeric(quantile(x, probs))
    names(q) <- paste0("Q", prob_labs)
    q
  })))
  qpreds$time <- unique(ppred$time)
  avg_pal <- viridis::mako(n=length(probs))
  
  if(plot){
    p <- ggplot2::ggplot(qpreds, ggplot2::aes(x = .data[["time"]] ))+
      lapply(seq(1,49,2),function(i) ggplot2::geom_ribbon(ggplot2::aes(ymin=.data[[paste0("Q",i)]],
                                                                       ymax=.data[[paste0("Q",100-i)]]),
                                                          fill=avg_pal[i],alpha=0.5))+
      pcv_theme()
    out[["plot"]]<-p
  }
  
  if(!is.null(hyp)){
    hyps_df <- do.call(rbind, lapply(unique(ppred$time), function(i){
      p<-brms::hypothesis(ppred[ppred$time==i,], hyp)$hypothesis$Post.Prob
      data.frame(time = i, prob = p)
    }))
    
    hyps_df$hyp <- hyp
    out[["hypothesis"]]<-hyps_df
    
    hyps_df$y <- 1.1 * qpreds[,1]
    if(plot){
      p <- p +
        ggplot2::geom_text(data = hyps_df, ggplot2::aes(x=.data[["time"]], y=.data[["y"]],
                                                        label = round(.data[['prob']], 2), angle=90, hjust=0))+
        ggplot2::labs(title = paste0("Posterior Probability of `", hyp, "`"))+
        ggplot2::coord_cartesian(ylim = c(NA, 1.15*max(hyps_df$y)))
      out[["plot"]]<-p
    }
    
  }
  return(out)
}


#' @description
#' Internal function for parsing formulae and brmsfit objects for plotting
#' @param fit a brmsfit object
#' @param form a pcvr style formula
#' 
#' @keywords internal
#' @noRd

.parseFormula <- function(fit, form){
  fitData<-fit$data
  parsed_form <- .parsePcvrForm(form, fitData)
  y <- parsed_form$y
  x <- parsed_form$x
  group <- parsed_form$group
  return(list(fitData, y, x, group))
}
