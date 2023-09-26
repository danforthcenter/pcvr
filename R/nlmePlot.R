#' Function to visualize common \link{nlme::nlme} growth models.
#' 
#' Models fit using \link{growthSS} inputs by \link{fitGrowth} (and similar models made through other means)
#'  can be visualized easily using this function. This will generally be called by \code{growthPlot}.
#' 
#' @param fit A model fit returned by \code{fitGrowth} with type="nlme".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm} part of the output) specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}
#' @param df An optional dataframe to use in plotting observed growth curves on top of the model.
#' This must be supplied for nlme models.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the avaiable data has not reached some point (such as asymptotic size).
#' @keywords growth-curve, logistic, gompertz, monomolecular, linear, exponential, power-law
#' @importFrom methods is
#' @import ggplot2
#' @importFrom stats predict update residuals
#' @importFrom nlme nlme nlme.formula
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
#' nlmePlot(fit, form = ss$pcvrForm, groups = NULL, df = ss$df, timeRange = NULL, boot=10)
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

nlmePlot<-function(fit, form, df = NULL, groups = NULL, timeRange = NULL, boot=1000){
  #* `get needed information from formula`
  x <-as.character(form)[3]
  y <-as.character(form)[2]
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
  yvals <- lapply(1:boot, function(i) {try(.predfun(stats::update(fit,
                                                                  data=.nlmePlotResample(fit, df, group, y)), new_data),
                                           silent=TRUE )} )
  yvals_mat <- do.call(cbind, yvals[which(unlist(lapply(yvals, is.numeric)))])
  
  while(ncol(yvals_mat) < boot){
    new<-try(.predfun(stats::update(fit, data=.nlmePlotResample(fit, df, group, y)), new_data), silent=TRUE)
    if(is.numeric(new)){
      yvals_mat <- cbind(yvals_mat, new)
    }
  }
  #* `get CIs from boostrap preds`
  predMat<-.nlmePlotCIs(yvals_mat, seq(from=99, to=1, by=-2)/100 )
  preds <- cbind(new_data, predMat)
  #* `plot CIs`
  avg_pal <- viridis::plasma(n=ncol(predMat))
  
  plot <- ggplot2::ggplot(preds, ggplot2::aes(x=.data[[x]]))+
    ggplot2::facet_wrap(paste0("~", group)) +
    lapply(seq(1,49,2)/100,function(i) ggplot2::geom_ribbon(ggplot2::aes(ymin=.data[[paste0("Q_",i)]],
                                                                         ymax=.data[[paste0("Q_",1-i)]]),
                                                            fill=avg_pal[i*100],alpha=0.5))+
    ggplot2::labs(x=x, y=y)+
    pcv_theme()
  return(plot)
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
.nlmePlotResample <- function(fitted,data,idvar="group", y="y") {
  pp <- stats::predict(fitted,levels=1)
  rr <- stats::residuals(fitted)
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
.predfun <- function(fm, new_data) {
  stats::predict(fm,newdata=new_data,level=0)
}

