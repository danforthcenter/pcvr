#' Function to visualize common \link{stats::nls} growth models.
#' 
#' Models fit using \link{growthSS} inputs by \link{fitGrowth} (and similar models made through other means)
#'  can be visualized easily using this function.
#' 
#' 
#' @param fit A model fit returned by \code{fitGrowth} with type="nls".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm} part of the output) specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' If the individual and group are specified then the observed growth lines are plotted.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param df An optional dataframe to use in plotting observed growth curves on top of the model.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the avaiable data has not reached some point (such as asymptotic size).
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

nlsPlot<-function(fit, form, groups = NULL, df = NULL, timeRange = NULL){
  #fit = fit_outer; form = ss$pcvrForm; groups = NULL; df = ss$df; timeRange = NULL
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
  #* `plot`
  plot<-ggplot(plotdf, ggplot2::aes(group = interaction(.data[[group]])))+
    facet_wrap(stats::as.formula(paste0("~",group)))+
    individual_lines+
    ggplot2::geom_line(ggplot2::aes(x=.data[[x]], y=.data[["pred"]]), color="#CC4678FF", linewidth=0.7)+ # using middle of plasma pal
    labs(x=x, y = as.character(form)[2])+
    pcv_theme()
  
  return(plot)
  
}













