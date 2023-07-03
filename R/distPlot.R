#' Function for plotting iterations of posterior distributions
#' 
#' @param fits A list of brmsfit objects following the same data over time. Currently checkpointing is not supported.
#' @param form A formula describing the growth model similar to \code{\link{growthSS}} and \code{\link{brmPlot}} such as: outcome ~ predictor |individual/group
#' @param priors a named list of samples from the prior distributions for each parameter in \code{params}. This is only used if sample_prior=F in the brmsfit object. If left NULL then no prior is included.
#' @param params a vector of parameters to include distribution plots of. Defaults to NULL which will use all parameters from the top level model.
#' @param d data used to fit models (this is used to plot each subject's trend line)
#' @param maxTime Optional parameter to designate a max time not observed in the models so far
#' @param patch Logical, should a patchwork plot be returned or should lists of ggplots be returned?
#' @keywords Bayesian, brms
#' @import ggplot2
#' @import patchwork
#' @import viridis
#' @return A ggplot or a list of ggplots (depending on patch).
#' @export
#' @examples 
#' 
#' ## Not run:
#' 
#' print(load(url("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/brmsFits.rdata")))
#' library(brms)
#' library(ggplot2)
#' library(patchwork)
#' fits<-list(fit_3, fit_15)
#' form <- y~time|id/group
#' priors<-list("phi1"=rlnorm(2000, log(130), 0.25), "phi2"=rlnorm(2000, log(12), 0.25), "phi3"=rlnorm(2000, log(3), 0.25))
#' params = c("A", "B", "C")
#' d <- simdf
#' maxTime = NULL
#' patch=T
#' from3to25<-list(fit_3, fit_5, fit_7, fit_9, fit_11, fit_13, fit_15, fit_17, fit_19, fit_21, fit_23, fit_25)
#' distributionPlot(fits = from3to25, form = y~time|id/group, params=params, d=d, priors=priors)
#' 
#' ## End(Not run)


distributionPlot<-function(fits, form, priors=NULL, params=NULL, d, maxTime=NULL, patch=T){
  #* ***** `Check args`
  if(missing(fits)){stop("A list of fits must be supplied")}
  if(missing(form)){stop("A formula must be supplied")}
  if(missing(d)){stop("Data used to fit the final model must be supplied")}
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
  palettes<-lapply(1:length(unique(fitData[[group]])), function(i) viridis::viridis(length(timeRange),begin=0.1,end=1,option = virOptions[i],direction = 1))
  names(palettes)<-unique(fitData[[group]])
  
  #* ***** `if params is null then pull them from growth formula`
  
  if(is.null(params)){
    fit<-fits[[1]]
    growthForm<-as.character(fit$formula[[1]])[[3]]
    
    test<-gsub(x, "", growthForm) # ;test
    test2<-gsub("exp\\(", "", test) # ; test2
    test3<-gsub("\\(1", "", test2) # ;test3
    test4<-gsub("[/]|[+]|[-]|[)]|[()]", "", test3)
    params<-strsplit(test4, "\\s+")[[1]]
    
    
    test3<-gsub("[)]|[()]","",test2)
    test3
    
    
  }
  
  #* ***** `growth trendline plots`
  
  growthTrendPlots<-lapply(1:length(dSplit), function(i){
    dt<-dSplit[[i]]
    ggplot2::ggplot(dt, ggplot2::aes(x=.data[[x]], y = .data[[y]], color = .data[[x]], group = .data[[individual]] ))+
      ggplot2::geom_line(show.legend=F)+
      viridis::scale_color_viridis(begin=0.1, end=1, option = virOptions[i], direction=1 )+
      ggplot2::scale_x_continuous(limits = c(startTime, endTime))+
      pcv_theme()
  })
  
  #* ***** `posterior distribution extraction`
  
  posts<-do.call(rbind, lapply(fits, function(fit){
    time<-max(fit$data[[x]], na.rm=T)
    fitDraws<-do.call(cbind, lapply(params, function(par){
      draws<-as.data.frame(fit)[grepl(par, colnames(as.data.frame(fit)))]
      if(nrow(prior_draws(fit))>1){draws<-draws[!grepl("^prior_", colnames(draws))]}
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
  } else { #* `need to fit some models with sample_prior and see how this works with them`
    prior_df<-prior_draws(fits[[1]])
    prior_df<-prior_df[ ,grepl(paste0("b_", paste0(params,collapse="|")), colnames(prior_df))]
    colnames(prior_df)<-gsub(group, "",colnames(prior_df))
    colnames(prior_df)<-gsub("^b_","",colnames(prior_df))
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
  
  lapply(posts,summary)
  
  xlims<-lapply(params, function(par){
    diff<-as.numeric(as.matrix(posts[,grepl(paste0("^",par,"_"),colnames(posts))]))
    c(min(diff,na.rm=T), max(diff,na.rm=T))
  })
  names(xlims)<-params
  postPlots<-lapply(unique(fitData[[group]]), function(groupVal){
    groupPlots<-lapply(params, function(par){
      
      p<-ggplot2::ggplot(posts)+
        ggplot2::geom_density(ggplot2::aes(x=.data[[paste(par, groupVal,sep="_")]],
                                       fill=.data[[x]],color=.data[[x]],
                                       group=.data[[x]]), alpha=0.8)+
        ggplot2::labs(x=paste(par, group, groupVal))+
        ggplot2::coord_cartesian(xlim=xlims[[par]])+
        pcv_theme()+
        ggplot2::theme(axis.text.x.bottom = ggplot2::element_text(angle=0),
              legend.position = "none", axis.title.y = ggplot2::element_blank())
      
      if(USEPRIOR){
        p<-p+ggplot2::scale_fill_manual(values = c("black",palettes[[groupVal]]) )+
          ggplot2::scale_color_manual(values = c("black",palettes[[groupVal]]) )
      } else{
        p<-p+ggplot2::scale_fill_manual(values = palettes[[groupVal]]) +
          ggplot2::scale_color_manual(values = palettes[[groupVal]])
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
    out<-patchPlot + patchwork::plot_layout(ncol=ncol_patch, nrow=nrow_patch)
  } else{
    out<-list(growthTrendPlots, postPlots)
  }
  return(out)
}

