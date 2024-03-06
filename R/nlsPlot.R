#' Function to visualize common \code{stats::nls} growth models.
#' 
#' Models fit using \link{growthSS} inputs by \link{fitGrowth} (and similar models made through other means)
#'  can be visualized easily using this function. This will generally be called by \code{growthPlot}.
#' 
#' @param fit A model fit returned by \code{fitGrowth} with type="nls".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm} part of the output) specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' If the individual and group are specified then the observed growth lines are plotted.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param df A dataframe to use in plotting observed growth curves on top of the model.
#' This must be supplied for nls models.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the avaiable data has not reached some point (such as asymptotic size).
#' @param facetGroups logical, should groups be separated in facets? Defaults to TRUE. 
#' @param groupFill logical, should groups have different colors? Defaults to FALSE. If TRUE then viridis colormaps are used in the order
#' of virMaps
#' @param virMaps order of viridis maps to use. Will be recycled to necessary length. Defaults to "plasma", but will generally be informed by growthPlot's default.
#' @keywords growth-curve, logistic, gompertz, monomolecular, linear, exponential, power-law
#' @importFrom methods is
#' @import ggplot2
#' @importFrom stats predict
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
#'  df=simdf, start=NULL, type="nls")
#' dim(ss$df)
#' 
#' fit<-fitGrowth(ss)
#' 
#' nlsPlot(fit, form=ss$pcvrForm, df = ss$pcvrForm)
#' 
#' }
#' 
#' ## End(Not run)
#' 
#' @return Returns a ggplot showing an nls model's predictions.
#' 
#' @export

nlsPlot<-function(fit, form, df = NULL, groups = NULL, timeRange = NULL, facetGroups=TRUE, groupFill=FALSE, virMaps = c("plasma")){
  #fit = fit_outer; form = ss$pcvrForm; groups = NULL; df = ss$df; timeRange = NULL
  #* `get needed information from formula`
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  if(individual == "dummyIndividual"){individual=NULL}
  group <- parsed_form$group
  df <- parsed_form$data
  #* `filter by groups if groups != NULL`
  if(!is.null(groups)){
    df <- df[df[[group]] %in% groups, ]
  }
  #* `make new data if timerange is not NULL`
  if(!is.null(timeRange)){
    new_data = do.call(rbind, lapply(unique(df[[group]]), function(g){
      stats::setNames(data.frame(g, timeRange), c(group, x))
    }))
  } else{
    new_data <- df
  }
  #* `add predictions`

  preds <- data.frame(pred=stats::predict(fit, newdata=new_data))
  keep <- which(!duplicated(preds$pred))
  plotdf <- df[keep,]
  plotdf$pred <- preds[keep,"pred"]

  #* `when implemented SE can be added here, see ?predict.nls`
  #* 
  #* `layer for individual lines if formula was complete`
  if(!is.null(individual)){
    individual_lines<-ggplot2::geom_line(data=df, ggplot2::aes(x=.data[[x]], y=.data[[y]],
                                                               group = interaction(.data[[individual]],
                                                                                   .data[[group]]) ),
                                         linewidth=0.25, color="gray40")
  } else{
    individual_lines<-list()
  }
  #* `facetGroups`
  if(facetGroups){
    facet_layer <- ggplot2::facet_wrap(stats::as.formula(paste0("~",group)))
  } else{ facet_layer <- NULL }
  #* `groupFill`
  if(groupFill){
    virVals <- unlist(lapply( rep(virMaps, length.out=length(unique( df[[group]] )) ), function(pal){
      viridis::viridis(1, begin=0.5, option = pal)
    }))
    color_scale <- ggplot2::scale_color_manual(values = virVals)
  } else{
    color_scale <- ggplot2::scale_color_manual(values = rep("#CC4678FF", length(unique( df[[group]] ))))
  }
  
  #* `plot`
  plot<-ggplot(plotdf, ggplot2::aes(group = interaction(.data[[group]])))+
    facet_layer+
    individual_lines+
    ggplot2::geom_line(ggplot2::aes(x=.data[[x]], y=.data[["pred"]], color = .data[[group]]), linewidth=0.7)+ # using middle of plasma pal
    color_scale+
    labs(x=x, y = as.character(form)[2])+
    pcv_theme()
  
  return(plot)
  
}













