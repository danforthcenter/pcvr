#' Function to visualize common \link{quantreg::nlrq} growth models.
#' 
#' Models fit using \link{growthSS} inputs by \link{fitGrowth} (and similar models made through other means)
#'  can be visualized easily using this function. This will generally be called by \code{growthPlot}.
#' 
#' @param fit A model fit, or list of model fits, returned by \code{fitGrowth} with type="nlrq".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm} part of the output) specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}. 
#' If the individual and group are specified then the observed growth lines are plotted.
#' @param df A dataframe to use in plotting observed growth curves on top of the model. This must be supplied for nlrq models.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the avaiable data has not reached some point (such as asymptotic size).
#' @keywords growth-curve, logistic, gompertz, monomolecular, linear, exponential, power-law
#' @importFrom methods is
#' @import ggplot2
#' @importFrom stats setNames predict
#' @importFrom viridis plasma
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
#'   tau=c(0.25, 0.5, 0.75), df=simdf, start=NULL, type="nlrq")
#' dim(ss$df)
#' 
#' fits<-fitGrowth(ss)
#' 
#' nlrqPlot(fits, form=ss$pcvrForm, df = ss$pcvrForm)
#' 
#' }
#' 
#' ## End(Not run)
#' 
#' @return Returns a ggplot showing an nlrq model's quantiles
#'  and optionally the individual growth lines.
#' 
#' @export

nlrqPlot<-function(fit, form, df = NULL, groups = NULL, timeRange = NULL){
  #fit = fits; form = ss$pcvrForm; groups = NULL; df = ss$df; timeRange = NULL
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
  #* `standardize fit class`
  if(methods::is(fit, "nlrq")){
    fit <- list(fit)
    names(fit)<-fit[[1]]$m$tau()
  }
  #return(names(fit))
  #* `add predictions and record taus`
  taus <- as.numeric(unlist(lapply(fit, function(f){f$m$tau()})))
  
  preds<-do.call(cbind, lapply(fit, function(f){
    tau <- f$m$tau()
    stats::setNames(data.frame(stats::predict(f, newdata=new_data)), paste0("Q_",tau))
  }))
  predCols <-colnames(preds)
  keep <- which(!duplicated(preds))
  plotdf <- cbind(df[keep,], preds[keep,])
  colnames(plotdf) <- c(colnames(df), colnames(preds))
  #* `define colors`
  virPal_p1<-viridis::plasma(ceiling(length(predCols)/2), direction=1, end=1)
  virPal_p2<-viridis::plasma(floor(length(predCols)/2), direction=-1, end=1)
  virPal<-c(virPal_p1,virPal_p2)#; scales::show_col(virPal)
  #* `layer for individual lines if formula was complete`
  if(!is.null(individual)){
    individual_lines<-ggplot2::geom_line(data=df, ggplot2::aes(x=.data[[x]], y=.data[[y]],
                                                      group = interaction(.data[[individual]],
                                                                          .data[[group]]) ),
                                         linewidth=0.25, color="gray40")
  } else{
    individual_lines<-list()
  }
  #* `plot`
  plot<-ggplot(plotdf, ggplot2::aes(group = interaction(.data[[group]])))+
    facet_wrap(stats::as.formula(paste0("~",group)))+
    individual_lines+
    lapply(1:length(predCols), function(i){
      ggplot2::geom_line(ggplot2::aes(x=.data[[x]], y=.data[[ predCols[i] ]]), color=virPal[i], linewidth=0.7)
    })+
    labs(x=x, y = as.character(form)[2])+
    pcv_theme()
  
  return(plot)
}

