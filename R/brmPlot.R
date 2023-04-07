#' Function to visualize brms models similar to those made using growthSS outputs.
#' 
#' @param fit
#' @keywords growth curve, logistic, gompertz, monomolecular, linear, exponential, power law
#' @import ggplot2
#' @import brms
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
#' @export
#' 
#' 

library(ggplot2)
simdf<-growthSim("logistic", n=20, t=25, params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="spline", df=simdf, priors = list("A"=130, "B"=12, "C"=3))
lapply(ss,class)
fit_test <- brm(ss$formula, prior = ss$prior, data = ss$df, family = ss$family, # main componenets of the model
               iter = 1000, cores = 2, chains = 2, init = ss$initfun, # parameters controling chain number, chain length, parallelization and starting values
               control = list(adapt_delta = 0.999, max_treedepth = 20), backend = "cmdstanr") # options to increase performance
brmPlot(fit_test, y~time|id/group, df=NULL)

brmPlot<-function(fit, form, df){
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
  
  p<-ggplot(predictions, aes(x=.data[[x]], y=Estimate))+
    facet_wrap(as.formula(paste0("~",group)))+
    lapply(seq(1,49,2),function(i) geom_ribbon(aes(ymin=.data[[paste0("Q",i)]],ymax=.data[[paste0("Q",100-i)]]),fill=avg_pal[i],alpha=0.5))+
    labs(x=x, y=y)+
    theme_light()+
    theme(axis.ticks.length=unit(0.2,"cm"))+
    theme(strip.background=element_rect(fill="gray50",color="gray20"),
          strip.text.x=element_text(size=14,color="white"),
          strip.text.y=element_text(size=14,color="white"))+
    theme(axis.title= element_text(size = 18))+
    theme(axis.text = element_text(size = 14))+
    theme(legend.position='top')
  
  if(!is.null(df)){
    p<-p+geom_line(data=fitData, aes(.data[[x]], .data[[y]], group=.data[[individual]]),color="gray20")
  }
  return(p)
}

