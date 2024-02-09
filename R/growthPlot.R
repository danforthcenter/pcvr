#' Function to visualize common growth models.
#' 
#' Models fit using \link{growthSS} inputs by \link{fitGrowth} (and similar models made through other means)
#'  can be visualized easily using this function.
#' 
#' 
#' @param fit A model fit object (or a list of \code{nlrq} models) as returned by \code{fitGrowth}.
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm} part of the output) specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param df A dataframe to use in plotting observed growth curves on top of the model and for making predictions.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the avaiable data has not reached some point (such as asymptotic size).
#' @param facetGroups logical, should groups be separated in facets? Defaults to TRUE. 
#' @param groupFill logical, should groups have different colors? Defaults to FALSE. If TRUE then viridis colormaps are used in the order
#' c('plasma', 'mako', 'viridis', 'inferno', 'cividis', 'magma', 'turbo', 'rocket'). Alternatively this can be given as a vector of
#' viridis colormap names to use in a different order than above.
#' @keywords growth-curve, logistic, gompertz, monomolecular, linear, exponential, power-law
#' @importFrom methods is
#' @examples 
#' 
#' ## Not run:
#' 
#' if(FALSE){
#' 
#' data(bw_vignette_fit)
#' brmPlot(bw_vignette_fit, y~time|id/group, df=NULL)
#' print(load(url("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/brmsFits.rdata")))
#' brmPlot(fit_25, form = y~time|id/group)
#' brmPlot(fit_9, form = y~time|id/group)
#' brmPlot(fit_15, form = y~time|id/group)
#' 
#' }
#' 
#' ## End(Not run)
#' 
#' @return Returns a ggplot showing a brms model's credible
#' intervals and optionally the individual growth lines.
#' 
#' @export

growthPlot<-function(fit, form, groups = NULL, df = NULL, timeRange = NULL, facetGroups=TRUE, groupFill=FALSE){
  
  if(is.logical(groupFill)){
    virMaps <- c('plasma', 'mako', 'viridis', 'inferno', 'cividis', 'magma', 'turbo', 'rocket')
  } else{
    virMaps <- groupFill
    groupFill <- TRUE
  }
  
  if(methods::is(fit, "brmsfit")){
    if(as.character(fit$family)[1]=="student"){
      plot <- brmPlot(fit=fit, form=form, groups = groups, df = df, timeRange = timeRange, facetGroups = facetGroups, groupFill = groupFill, virMaps)
    } else{
      plot <- brmSurvPlot(fit=fit, form=form, groups = groups, df = df, timeRange = timeRange, facetGroups = facetGroups, groupFill = groupFill, virMaps)
    }
  } else if(methods::is(fit, "nls") | methods::is(fit, "gam") | methods::is(fit, "lm") ){
    plot <- nlsPlot(fit=fit, form=form, groups = groups, df = df, timeRange = timeRange, facetGroups = facetGroups, groupFill = groupFill, virMaps)
  } else if(methods::is(fit, "lme")){
    plot <- nlmePlot(fit, form = form, groups = groups, df = df, timeRange = timeRange, facetGroups = facetGroups, groupFill = groupFill, virMaps)
  } else if(methods::is(fit, "nlme")){
    plot <- nlmePlot(fit, form = form, groups = groups, df = df, timeRange = timeRange, facetGroups = facetGroups, groupFill = groupFill, virMaps)
  } else if(methods::is(fit, "nlrq") | is.list(fit) && methods::is(fit[[1]], "nlrq") ){
    plot <- nlrqPlot(fit=fit, form=form, groups = groups, df = df, timeRange = timeRange, facetGroups = facetGroups, groupFill = groupFill, virMaps)
  } else if(methods::is(fit, "rq") | is.list(fit) && methods::is(fit[[1]], "rq") ){
    plot <- rqPlot(fit=fit, form=form, groups = groups, df = df, timeRange = timeRange, facetGroups = facetGroups, groupFill = groupFill, virMaps)
  } else if( methods::is(fit, "survreg") ){
    plot <- survRegPlot(fit=fit, form=form, groups = groups, df = df, facetGroups = facetGroups, groupFill = groupFill, virMaps)
  } else if( methods::is(fit, "flexsurvreg") ){
    plot <- flexSurvRegPlot(fit=fit, form=form, groups = groups, df = df, facetGroups = facetGroups, groupFill = groupFill, virMaps)
  } else{
    stop("Only implemented for brms, nls, nlme, mgcv, and nlrq models from fitGrowth currently.")
  }
  
  return(plot)
}










