#' Function to visualize common \link{nlme::nlme} growth models.
#' 
#' Models fit using \link{growthSS} inputs by \link{fitGrowth} (and similar models made through other means)
#'  can be visualized easily using this function.
#' 
#' 
#' @param fit A model fit returned by \code{fitGrowth} with type="nlme".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm} part of the output) specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param df An optional dataframe to use in plotting observed growth curves on top of the model.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the avaiable data has not reached some point (such as asymptotic size).
#' @keywords growth-curve, logistic, gompertz, monomolecular, linear, exponential, power-law
#' @importFrom methods is
#' @import ggplot2
#' @importFrom stats predict
#' @importFrom nlme nlme
#' @examples 
#' 
#' ## Not run:
#' 
#' if(FALSE){
#' 
#' simdf<-growthSim("logistic", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' 
#' ss<-growthSS(model = "logistic", form=y~time|id/group,
#'  df=simdf, sigma="power", start=NULL, type="nlme")
#' dim(ss$df)
#' 
#' fit<-fitGrowth(ss)
#' 
#' 
#' 
#' }
#' 
#' ## End(Not run)
#' 
#' @return Returns a ggplot showing a brms model's credible
#' intervals and optionally the individual growth lines.
#' 
#' @export
#' 

nlsPlot<-function(fit, form, groups = NULL, df = NULL, timeRange = NULL){
  #* `get needed information from formula`
  x <-as.character(form)[3]
  if(grepl("\\|", x) & grepl("\\/",x)){
    x3<-trimws(strsplit(x, "[|]|[/]")[[1]])
    x<-x3[1]
    individual = x3[2]
    group = x3[3] 
  } else if (grepl("\\|", x)){
    x2<-trimws(strsplit(x, "[|]")[[1]])
    x<-x2[1]
    group = x2[2]
    individual=NULL
  }
  #* `filter by groups if groups != NULL`
  if(!is.null(groups)){
    df <- df[df[[groups]] %in% groups, ]
  }
  #* `make new data if timerange is not NULL`
  if(!is.null(timeRange)){
    new_data = do.call(rbind, lapply(unique(df[[group]]), function(g){
      stats::setNames(data.frame(g, timeRange), c(group, x))
    }))
  } else{
    new_data <- df
  }
  #* `Generate bootstrap predictions`
  yvals <- lapply(1:Nresample, function(i) {try(pfun(update(fm1, data=sampfun(fit, ss$df,"group", "y"))),silent=TRUE )} )
  yvals_mat <- do.call(cbind, yvals[which(unlist(lapply(yvals, is.numeric)))])
  dim(yvals_mat)
  while(ncol(yvals_mat) < Nresample){
    new<-try(pfun(update(fm1, data=sampfun(fit, ss$df,"group", "y"))), silent=TRUE)
    if(is.numeric(new)){
      yvals_mat <- cbind(yvals_mat, new)
    }
  }
  #* `get CIs from boostrap preds`
  predMat<-.nlmePlotCIs(yvals_mat, seq(from=99, to=1, by=-2)/100 )
  preds <- cbind(pframe, predMat)
  #* `plot CIs`
  avg_pal <- viridis::plasma(n=ncol(predMat))
  
  ggplot2::ggplot(preds, ggplot2::aes(x=.data[[x]]))+
    facet_wrap(~group) +
    lapply(seq(1,49,2)/100,function(i) ggplot2::geom_ribbon(ggplot2::aes(ymin=.data[[paste0("Q_",i)]],
                                                                         ymax=.data[[paste0("Q_",1-i)]]),
                                                            fill=avg_pal[i*100],alpha=0.5))+
    ggplot2::labs(x=x, y=y)+
    pcv_theme()
}

#' convenience function for pulling CIs
#' @keywords internal
#' @noRd
.nlmePlotCIs <- function(y, quantiles = c(0.025, 0.975), prefix="Q_") {
  r1 <- t(apply(y,1,quantile, quantiles))
  setNames(as.data.frame(r1),paste0(prefix,quantiles))
}

#' Resampling function for bootstrapping CIs
#' @keywords internal
#' @noRd
.nlmePlotResample <- function(fitted,data,idvar="Seed", y="y") {
  pp <- predict(fitted,levels=1)
  rr <- residuals(fitted)
  dd <- data.frame(data, pred=pp, res=rr)
  ## sample groups with replacement
  iv <- levels(data[[idvar]])
  bsamp1 <- sample(iv, size=length(iv),replace=TRUE)
  res <- do.call(rbind, lapply(bsamp1,
                   function(x) {
                     ## within groups, sample *residuals* with replacement
                     ddb <- dd[dd[[idvar]]==x,]
                     ## bootstrapped response = pred + bootstrapped residual
                     ddb[[y]] <- ddb$pred +
                       sample(ddb$res,size=nrow(ddb),replace=TRUE)
                     return(ddb)
                   }))
  return(res)
}

#' Prediction function for bootstrapped values
#' @keywords internal
#' @noRd
pfun <- function(fm) {
  predict(fm,newdata=pframe,level=0)
}


if(FALSE){
  
  library(nlme)
  library(MASS)
  
  # fm1 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
  #             data = Loblolly,
  #             fixed = Asym + R0 + lrc ~ 1,
  #             random = Asym ~ 1,
  #             start = c(Asym = 103, R0 = -8.5, lrc = -3.3))
  # 
  # xvals <-  with(Loblolly,seq(min(age),max(age),length.out=100))
  fm1<-fit
  y="y"
  xvals <-  seq(min(ss$df$time),max(ss$df$time),length.out=100)
  Nresample <- 10
  ## pick new parameter values by sampling from multivariate normal distribution based on fit
  pars.picked <- mvrnorm(nresamp, mu = fixef(fm1), Sigma = vcov(fm1))
  
  ## predicted values: useful below
  # pframe <- with(Loblolly,data.frame(age=xvals))
  pframe <- with(ss$df,data.frame(time=xvals))
  pframe <- expand.grid(time = xvals, group = unique(ss$df$group))
  pframe[[y]] <- predict(fm1, newdata=pframe, level=0)
  
  ## utility function
  get_CI <- function(y, quantiles = c(0.025, 0.975), prefix="Q_") {
    r1 <- t(apply(y,1,quantile, quantiles))
    setNames(as.data.frame(r1),paste0(prefix,quantiles))
  }
  
  ## bootstrapping
  sampfun <- function(fitted,data,idvar="Seed", y="y") {
    pp <- predict(fitted,levels=1)
    rr <- residuals(fitted)
    dd <- data.frame(data, pred=pp, res=rr)
    ## sample groups with replacement
    iv <- levels(data[[idvar]])
    bsamp1 <- sample(iv, size=length(iv),replace=TRUE)
    bsamp2 <- lapply(bsamp1,
                     function(x) {
                       ## within groups, sample *residuals* with replacement
                       ddb <- dd[dd[[idvar]]==x,]
                       ## bootstrapped response = pred + bootstrapped residual
                       ddb[[y]] <- ddb$pred +
                         sample(ddb$res,size=nrow(ddb),replace=TRUE)
                       return(ddb)
                     })
    res <- do.call(rbind,bsamp2)  ## collect results
    return(res)
  }
  
  pfun <- function(fm) {
    predict(fm,newdata=pframe,level=0)
  }
  
  yvals <- lapply(1:Nresample, function(i) {try(pfun(update(fm1, data=sampfun(fit, ss$df,"group", "y"))),silent=TRUE )} )
  yvals_mat <- do.call(cbind, yvals[which(unlist(lapply(yvals, is.numeric)))])
  dim(yvals_mat)
  while(ncol(yvals_mat) < Nresample){
    new<-try(pfun(update(fm1, data=sampfun(fit, ss$df,"group", "y"))), silent=TRUE)
    if(is.numeric(new)){
      yvals_mat <- cbind(yvals_mat, new)
    }
  }
  predMat<-get_CI(yvals_mat, seq(from=99, to=1, by=-2)/100 )
  preds <- cbind(pframe, predMat)
  avg_pal <- viridis::plasma(n=ncol(predMat))
  
  ggplot2::ggplot(preds, ggplot2::aes(x=.data[[x]]))+
    facet_wrap(~group) +
    lapply(seq(1,49,2)/100,function(i) ggplot2::geom_ribbon(ggplot2::aes(ymin=.data[[paste0("Q_",i)]],
                                                                     ymax=.data[[paste0("Q_",1-i)]]),
                                                        fill=avg_pal[i*100],alpha=0.5))+
    ggplot2::labs(x=x, y=y)+
    pcv_theme()
  
  
  
}
























