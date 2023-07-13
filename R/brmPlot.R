#' Function to visualize brms models similar to those made using growthSS outputs.
#' 
#' @param fit A brmsfit object, similar to those fit with \code{\link{growthSS}} outputs.
#' @param form A formula similar to that in \code{growthSS} inputs specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param df An optional dataframe to use in plotting observed growth curves on top of the model.
#' @keywords growth-curve, logistic, gompertz, monomolecular, linear, exponential, power-law
#' @import ggplot2
#' @import viridis
#' @importFrom stats as.formula
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

brmPlot<-function(fit, form, groups = NULL, df=NULL){
  fitData<-fit$data
  y=as.character(form)[2]
  x<-as.character(form)[3]
  if(grepl("\\|", x) | grepl("\\/",x)){
    x3<-trimws(strsplit(x, "[|]|[/]")[[1]])
    x<-x3[1]
    individual = x3[2]
    group = x3[3]
  } else {stop("form must specify grouping for observations. See documentation and examples.")}
  probs <- seq(from=99, to=1, by=-2)/100
  avg_pal <- viridis::plasma(n=length(probs))
  
  newData<-data.frame(x=rep(unique(fitData[[x]]), times=length(unique(fitData[[group]]))), 
                      group = rep(unique(fitData[[group]]), each = length(unique(fitData[[x]]))),
                      individual = rep(paste0("new_",1:length(unique(fitData[[group]]))), each=length(unique(fitData[[x]]))) )
  colnames(newData)<-c(x, group, individual)
  predictions <- cbind(newData, predict(fit, newData, probs=probs))
  
  if(!is.null(groups)){
    predictions<-predictions[predictions$group %in% groups, ]
    if(!is.null(df)){
      df<-df[df[[group]] %in% groups, ]
    }
  }
  p<-ggplot2::ggplot(predictions, ggplot2::aes(x=.data[[x]], y=.data$Estimate))+
    ggplot2::facet_wrap(as.formula(paste0("~",group)))+
    lapply(seq(1,49,2),function(i) ggplot2::geom_ribbon(ggplot2::aes(ymin=.data[[paste0("Q",i)]],
                                                                     ymax=.data[[paste0("Q",100-i)]]),
                                                        fill=avg_pal[i],alpha=0.5))+
    ggplot2::labs(x=x, y=y)+
    pcv_theme()
  
  if(!is.null(df)){
    p<-p+ggplot2::geom_line(data=df, ggplot2::aes(.data[[x]], .data[[y]], group=.data[[individual]]),color="gray20", linewidth=0.2)
  }
  return(p)
}

