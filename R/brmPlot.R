#' Function to visualize brms models similar to those made using growthSS outputs.
#' 
#' @param fit A brmsfit object, similar to those fit with \code{\link{growthSS}} outputs.
#' @param form A formula similar to that in \code{growthSS} inputs specifying the outcome, predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' @param df An optional dataframe to use in plotting observed growth curves on top of the model. 
#' @keywords growth curve, logistic, gompertz, monomolecular, linear, exponential, power law
#' @import ggplot2
#' @import viridis
#' @examples 
#' simdf<-growthSim("logistic", n=20, t=25, params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="spline", df=simdf, priors = list("A"=130, "B"=12, "C"=3))
#' lapply(ss,class)
#' ss$initfun()
#' fit_test <- brm(ss$formula, prior = ss$prior, data = ss$df, family = ss$family, # main componenets of the model
#'               iter = 1000, cores = 2, chains = 2, init = ss$initfun, # parameters controling chain number, chain length, parallelization and starting values
#'               control = list(adapt_delta = 0.999, max_treedepth = 20), backend = "cmdstanr") # options to increase performance
#' brmPlot(fit_test, y~time|id/group, df=NULL)
#' print(load("/home/jsumner/Desktop/stargate/bayesian_growth/earlyStoppingSim/bayesian_updating/parameterUpdating_DataAndModels_3to25.rdata"))
#' brmPlot(fit_25, form = y~time|sample/treatment)
#' brmPlot(fit_9, form = y~time|sample/treatment)
#' brmPlot(fit_15, form = y~time|sample/treatment)
#' 
#' @export

brmPlot<-function(fit, form, df=NULL){
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
  
  p<-ggplot2::ggplot(predictions, ggplot2::aes(x=.data[[x]], y=Estimate))+
    ggplot2::facet_wrap(as.formula(paste0("~",group)))+
    lapply(seq(1,49,2),function(i) ggplot2::geom_ribbon(ggplot2::aes(ymin=.data[[paste0("Q",i)]],ymax=.data[[paste0("Q",100-i)]]),fill=avg_pal[i],alpha=0.5))+
    ggplot2::labs(x=x, y=y)+
    pcv_theme()
  
  if(!is.null(df)){
    p<-p+ggplot2::geom_line(data=fitData, ggplot2::aes(.data[[x]], .data[[y]], group=.data[[individual]]),color="gray20")
  }
  return(p)
}

