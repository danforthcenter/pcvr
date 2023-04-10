#' Function for plotting iterations of posterior distributions
#' 
#' This should take a formula and a list of fits (which should be able to be length 1)
#' Then it should plot the growth trend and the posterior distributions at different times.
#' See examples from ~/Desktop/stargate/bayesian_growth/earlyStoppingSim/bayesian_updating/distributionUpdating.R
#' and: /home/jsumner/Desktop/stargate/bayesian_growth/earlyStoppingSim/bayesian_updating/bayesianUpdatingFigures/logisticPosteriorUpdating_withData.png
#' goal is to make a function that does that whole file worth of stuff a whole lot cleaner.
#' what does this function need:
#' list of models/model [arg: fits]
#' if sample_prior was not true in the model then priors [arg: priors]
#' a vector of parameters to make distribution plots of [arg: params]
#' a dataframe used to draw the growth curves (possibly optional, could simulate data based on the model instead?)[arg: df]
#' color palates should be infered using LETTERS and viridis to maximize row to row difference
#' this will need patchwork and some flexible grobbing
#' 
#' @param fit A brmsfit object to extract draws from
#' @param update A dataframe of previous output from this function
#' @keywords Bayesian, brms
#' @import brms
#' @import ggplot2
#' @import patchwork
#' @return A named list of elements to make it easier to fit common brms models.
#' @export
#' @examples 
#' 
#' 
print(load("/home/jsumner/Desktop/stargate/bayesian_growth/earlyStoppingSim/bayesian_updating/parameterUpdating_DataAndModels_3to25.rdata"))
library(brms)
library(ggplot2)
library(patchwork)
fits<-list(fit_3, fit_15)
form <- y~time|sample/treatment
priors<-list("phi1"=rlnorm(2000, log(130), 0.25), "phi2"=rlnorm(2000, log(12), 0.25), "phi3"=rlnorm(2000, log(3), 0.25))
params = c("phi1", "phi2", "phi3")
d <- df
maxTime = NULL
patch=T


distributionPlot<-function(fits, form, priors, params, d, maxTime, patch){
  #* ***** `Check args`
  
  #* ***** `Reused helper variables`
  y=as.character(form)[2]
  x<-as.character(form)[3]
  if(grepl("\\|", x) | grepl("\\/",x)){
    x3<-trimws(strsplit(x, "[|]|[/]")[[1]])
    x<-x3[1]
    individual = x3[2]
    group = x3[3]
  } else {stop("form must specify grouping for observations. See documentation and examples.")}
  fitData<-fits[[length(fits)]]$data
  dSplit<-split(d, d[[group]])
  startTime<-min(unlist(lapply(fits, function(ft){min(ft$data[[x]], na.rm=T)})))
  if(is.null(maxTime)){
    endTime<-max(unlist(lapply(fits, function(ft){max(ft$data[[x]], na.rm=T)})))
  }
  byTime <- mean(diff(unlist(lapply(fits, function(ft){max(ft$data[[x]], na.rm=T)}))))
  
  timeRange<-seq(startTime, endTime, byTime)
  
  virOptions<-c('C', 'G', 'B', 'D', 'A', 'H', 'E', 'F')
  palettes<-lapply(1:length(unique(fitData[[group]])), function(i) viridis::viridis(length(timeRange), begin=0.1, end=0.9, option = virOptions[i], direction = 1))
  names(palettes)<-unique(fitData[[group]])
  
  #* ***** `growth trendline plots`
  
  growthTrendPlots<-lapply(1:length(dSplit), function(i){
    dt<-dSplit[[i]]
    ggplot(dt, aes(x=.data[[x]], y = .data[[y]], color = .data[[x]], group = .data[[individual]] ))+
      geom_line(show.legend=F)+
      scale_color_viridis(begin=0.1, end=0.9, option = virOptions[i], direction=1 )+
      pcv_theme()
  })
  
  #* ***** `posterior distribution extraction`
  
  posts<-do.call(rbind, lapply(fits, function(fit){
    time<-max(fit$data[[x]], na.rm=T)
    fitDraws<-do.call(cbind, lapply(params, function(par){
      draws<-as.data.frame(fit)[grepl(par, colnames(as.data.frame(fit)))]
      splits<-strsplit(colnames(draws), split = "")
      mx<-max(unlist(lapply(splits,length)))
      ind<-which(unlist(lapply(1:mx, function(i) {length(unique(rapply(splits, function(j) {j[i]})))!=1 })))
      if(length(ind)>0){
        colnames(draws)<-paste(par,unlist(lapply(colnames(draws), function(c){ substr(c,min(ind), max(c(ind, nchar(c)))) })),sep="_")
      }
      draws
    }))
    fitDraws$time<-time
    fitDraws
  } ))
  
  #* ***** `prior distribution extraction`
  #* If the model was fit with sample_prior=T then this is super easy
  #* otherwise, I need to extract and parse the prior
  #* OR I can ask for an argument describing the prior/giving an example vector for the prior.
  
  if(all(unlist(lapply(fits, function(fit) nrow(prior_draws(fit))<1 )))){ # if no models were fit with sample_prior = T
    if(!is.null(prior)){ # if prior is supplied as argument
      USEPRIOR=T
      if(class(priors[[1]])!="list"){
        priors<-lapply(1:length(unique(d[[group]])), function(i) priors)
        names(priors)<-unique(d[[group]])
      }
      prior_df<-do.call(cbind, lapply(names(priors), function(nm){
        nmp<-priors[[nm]]
        setNames(data.frame(do.call(cbind, lapply(names(nmp), function(nmpn){
          nmp[[nmpn]]
        }))), paste0(names(nmp), "_", nm))
      }))
      prior_df[[x]]<-0
    } else{USEPRIOR = F}
  } else { #* `need to fit some mdoels with sample_prior and see how this works with them`
    lapply(fits, function(fit){
      p<-prior_draws(fit)
      head(p)
    })
    
    
    prior_df[[x]]<-0
    
    USEPRIOR=T
  }
  
  #* ***** `posterior distribution plots`
  #* need to assign ordering of factors 
  #* if USEPRIOR then join data, don't make separate geom
  
  if(USEPRIOR){
    posts<-rbind(prior_df, posts)
  }
  posts[[x]]<-factor(posts[[x]], levels = sort(as.numeric(unique(posts[[x]]))), ordered=T)
  
  postPlots<-lapply(unique(fitData[[group]]), function(groupVal){
    groupPlots<-lapply(params, function(par){
      
      p<-ggplot(posts)+
        geom_density(aes(x=.data[[paste(par, groupVal,sep="_")]],
                                       fill=.data[[x]],color=.data[[x]],
                                       group=.data[[x]]), alpha=0.5)+
        labs(x=paste(par, group, groupVal))+
        pcv_theme()+
        theme(axis.text.x.bottom = element_text(angle=0),
              legend.position = "none", axis.title.y = element_blank())
      
      if(USEPRIOR){
        p<-p+scale_fill_manual(values = c("black",palettes[[groupVal]]) )+
          scale_color_manual(values = c("black",palettes[[groupVal]]) )
      } else{
        p<-p+scale_fill_manual(values = palettes[[groupVal]]) +
          scale_color_manual(values = palettes[[groupVal]])
      }
      p
    })
  })
  
  if(patch){
    ncol_patch = 1+length(params)
    nrow_patch = length(postPlots)
    
    patchPlot = growthTrendPlots[[1]]+postPlots[[1]]
    if(length(unique(d[[group]]))>1){
      for(i in 2:length(growthTrendPlots)){
        patchPlot <- patchPlot + growthTrendPlots[[i]]+postPlots[[i]]
      }
    }
    out<-patchPlot + plot_layout(ncol=ncol_patch, nrow=nrow_patch)
  } else{
    out<-list(growthTrendPlots, postPlots)
  }
  return(out)
}

