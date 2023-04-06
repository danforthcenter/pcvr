#' brms helper function to plot prior or posterior distributions for a given parameter
#' 
#' 
#' 
#' 
#' 
#' Goal:
#' take brmsfit object
#' arg: parameter
#' arg: dist (prior or posterior)
#' should extract draws and plot density/plot density as a function for the prior.
#' 
#' 
#' Another function should do the ribbon plotting pretty well automatically.

setwd("/home/jsumner/Desktop/stargate/bayesian_growth/earlyStoppingSim/bayesian_updating")
dir()
print(load("parameterUpdating_DataAndModels_3to25.rdata"))

library(brms)

fit=fit_11
parameter = "phi1"
dist = "posterior"
func<-rlnorm
n=4000

brmsDist(fit, parameter, dist= "posterior", func,n)

brmsDist<-function(fit, parameter=NULL, dist=NULL, func=NULL, n=4000){
  if(dist=="prior"){
    if(nrow(prior_draws(fit, variable = parameter))>1){ # if sample.prior = T in the model
      x<-prior_draws(fit, variable = parameter)
    } else{
      x<-as.data.frame(fit$prior)
      x<-x[nchar(x$prior)>0,]
      col<-colnames(x)[grepl(parameter, x)]
      pri<-x[x[[col]]==parameter,"prior"]
      names(pri)<-paste0("prior", x[x[[col]]==parameter,"coef"])
      #pri<-pri[nchar(names(pri))>5]
      if(is.null(func)){
        warning(paste0("Supply a function to use to draw random deviations from ", unique(pri)) )
      }else{
        x<-data.frame(draw=1:n)
        pris<-lapply(pri, function(p){
          # p="lognormal(log(130), 0.25)"
          p<-sub("[a-zA-Z]+\\(", "", p)
          p<-sub("\\)$", "", p)
          p<-trimws( strsplit(p,",")[[1]] )
          p<-lapply(p, function(n) eval(parse(text = n)))
          return(p)
        })
        for(nm in names(pri)){
          x[[which(names(pri)==nm)]]<-func(n, pris[[nm]][[1]], pris[[nm]][[2]] )
        }
        colnames(x)<-names(pri)
        }
    }
  } else{
    x<-as.data.frame(fit)
    x<-x[grepl(parameter, colnames(x))]
  }
  splits<-strsplit(colnames(x), split = "")
  mx<-max(unlist(lapply(splits,length)))
  ind<-which(unlist(lapply(1:mx, function(i) {length(unique(rapply(splits, function(j) {j[i]})))!=1 })))
  if(length(ind)>0){
    colnames(x)<-unlist(lapply(colnames(x), function(c){ substr(c,min(ind), max(c(ind, nchar(c)))) }))
  }
  x<-reshape2::melt(x, measure.vars=colnames(x), variable.name = "group")
  
  x2<-as_draws(fit)
  x2<-do.call(rbind, lapply(1:length(x2), function(i){
    nms<-names(x2[[i]])[grepl(parameter, names(x2[[i]]))]
    splits<-strsplit(nms, split = "")
    mx<-max(unlist(lapply(splits,length)))
    ind<-which(unlist(lapply(1:mx, function(i) {length(unique(rapply(splits, function(j) {j[i]})))!=1 })))
    if(length(ind)>0){
      newnames<-unlist(lapply(nms, function(c){ substr(c,min(ind), max(c(ind, nchar(c)))) }))
    } else{newnames<-nms}
    newnames
    out<-setNames(data.frame(do.call(cbind, lapply(nms, function(nm) x2[[i]][[nm]]))), newnames)
    out$chain <-i
    out$chainIter <-1:nrow(out)
    out<-reshape2::melt(out, id.vars = c("chain", "chainIter"))
    return(out)
  }))
  
  ggplot(x2, aes(x=chainIter, y = value, color=as.character(chain)))+
    facet_wrap(~variable)+
    geom_line()
  
  ggplot(x, aes(x=value, fill=group))+
    geom_density(alpha=0.7)+
    pcv_theme()
  
  return(head(x))
}















