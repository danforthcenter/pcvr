#' Growth data simulating function
#' 
#' @description growthSim can be used to help pick reasonable parameters for common
#'  growth models to use in prior distributions or to simulate data for example models/plots.
#' 
#' @param model One of "logistic", "gompertz", "monomolecular", "exponential",
#' "linear", "power law", "double logistic", or "double gompertz". Alternatively this can be 
#' a pseudo formula to generate data from a segmented growth curve by specifying "model1 + model2",
#' see examples and \code{\link{growthSS}}. While "gam" models are supported by \code{growthSS}
#' they are not simulated by this function.
#' @param n Number of individuals to simulate over time per each group in params
#' @param t Max time (assumed to start at 1) to simulate growth to as an integer.
#' @param params A list of numeric parameters. A, B, C notation is used in the order that parameters appear in the formula (see examples).
#' Number of groups is inferred from the length of these vectors of parameters. In the case of the "double" models there are also 
#' A2, B2, and C2 terms. Changepoints should be specified as "changePointX" or "fixedChangePointX" as in \code{\link{growthSS}}.
#' @param noise Optionally this can be used to add specific amounts of
#' noise to the input parameters by specifying a list similar to params.
#' If NULL (the default) then data is simulated with 10\% random noise like: param + N(0, 0.1*param).
#' This exists for fringe cases and should generally be left NULL.
#' 
#' @keywords growth curve, logistic, gompertz, monomolecular, linear, exponential, power-law
#' @return Returns a dataframe of example growth data following the input parameters.
#' 
#' @importFrom stats rnorm setNames
#' 
#' @examples 
#' 
#' ## Not run:
#' library(ggplot2)
#' simdf<-growthSim("logistic", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+ 
#' geom_line(aes(color=group))+labs(title="Logistic")
#' 
#' simdf<-growthSim("gompertz", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(0.2, 0.25)))
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+
#'  geom_line(aes(color=group))+labs(title="Gompertz")
#'  
#' simdf<-growthSim("double logistic", n=20, t=70,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5),
#'  "A2"=c(400,300), "B2"=c(35, 40), "C2"=c(3.25, 2.75)))
#' ggplot(simdf, aes(time, y, group=interaction(group,id)))+
#'  geom_line(aes(color=group))+labs(title="Double Logistic")
#'  
#' simdf<-growthSim("double gompertz", n=20, t=100,
#' params = list("A"=c(180,140), "B"=c(13, 11), "C"=c(0.2, 0.2),
#'  "A2"=c(400,300), "B2"=c(50,50), "C2"=c(0.1, 0.1)))
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+
#'  geom_line(aes(color=group))+labs(title="Double Gompertz")
#' 
#' simdf<-growthSim("monomolecular", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(0.08, 0.1)))
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+ 
#' geom_line(aes(color=group))+
#' labs(title="Monomolecular")
#' 
#' simdf<-growthSim("exponential", n=20, t=25,
#' params = list("A"=c(15, 20), "B"=c(0.095, 0.095)))
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+ 
#' geom_line(aes(color=group))+labs(title="Exponential")
#' 
#' simdf<-growthSim("linear", n=20, t=25,
#' params = list("A"=c(1.1, 0.95)))
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+ 
#' geom_line(aes(color=group))+labs(title="Linear")
#' 
#' simdf<-growthSim("power law", n=20, t=25,
#' params = list("A"=c(16, 11), "B"=c(0.75, 0.7)))
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+ 
#' geom_line(aes(color=group))+labs(title="Power Law")
#' 
#' # simulating models from segmented growth models
#' 
#' simdf<-growthSim(model = "linear + linear", n=20, t=25,
#' params = list("linear1A"=c(16, 11), "linear2A"=c(0.75, 0.7), "changePoint1" = c(11, 14)))
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+ 
#' geom_line(aes(color=group))+labs(title="linear + linear")
#' 
#' simdf<-growthSim(model = "linear + linear + logistic", n=20, t=50,
#'   params = list("linear1A"=c(16, 11), "linear2A"=c(3, 4), # linear slopes, very intuitive
#'    "changePoint1" = c(11, 14), "changePoint2" = c(10, 12),
#'    # changepoint1 is standard, changepoint2 happens relative to changepoint 1
#'    "logistic3A" = c(200, 210), "logistic3B" = c(20, 25), "logistic3C"=c(3,3)  ))
#'    # similar to changepoint2, the asymptote and inflection point are relative to the starting
#'    # point of the logistic growth component. This is different than the model output
#'    # if you were to fit a curve to this model using `growthSS`.
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+ 
#' geom_line(aes(color=group))+ labs(title="linear + linear + logistic")
#' 
#' ## End(Not run)
#' 
#' @details 
#'     The \code{params} argument requires some understanding of how each growth model is parameterized.
#'     Examples of each are below should help, as will the examples.
#'     \itemize{
#'     \item \bold{Logistic}: `A / (1 + exp( (B-x)/C) )`
#'     Where A is the asymptote, B is the inflection point, C is the growth rate. 
#'     \item \bold{Gompertz}: `A * exp(-B * exp(-C*x))` 
#'     Where A is the asymptote, B is the inflection point, C is the growth rate. 
#'     \item \bold{Double Logistic}: `A / (1+exp((B-x)/C)) + ((A2-A) /(1+exp((B2-x)/C2)))`
#'     Where A is the asymptote, B is the inflection point, C is the growth rate,
#'     A2 is the second asymptote, B2 is the second inflection point, and C2 is the second 
#'     growth rate.
#'     \item \bold{Double Gompertz}: `A * exp(-B * exp(-C*x)) + ((A2-A) * exp(-B2 * exp(-C2*(x-B))))`
#'     Where A is the asymptote, B is the inflection point, C is the growth rate,
#'     A2 is the second asymptote, B2 is the second inflection point, and C2 is the second 
#'     growth rate.
#'     \item \bold{Monomolecular}: `A-A * exp(-B * x)`
#'     Where A is the asymptote and B is the growth rate. 
#'     \item \bold{Exponential}: `A * exp(B * x)` 
#'     Where A is the scale parameter and B is the growth rate. 
#'     \item \bold{Linear}: `A * x` 
#'     Where A is the growth rate.
#'     \item \bold{Power Law}: `A * x^(B)` 
#'     Where A is the scale parameter and B is the growth rate.
#'     }
#'     Note that for these distributions parameters do not exist in a vacuum.
#'     Changing one will make the others look different in the resulting data.
#'     The examples are a good place to start if you are unsure what parameters to use.
#' 
#' @export
#' 
#' 

growthSim<-function(model=c("logistic", "gompertz", "double logistic", "double gompertz", "monomolecular", "exponential", "linear", "power law"), n=20, t=25, params=list(), noise=NULL){
  
  if(length(model)>1){stop("Select one model to use")}
  if(is.null(noise)){noise = lapply(params, function(i) mean(i)/10) ; wasNull = TRUE} else{wasNULL=FALSE}
  if(is.null(names(noise))){names(noise)<-c(LETTERS[1:length(noise)])}
  if(any(names(noise)%in%letters)){ names(noise)<-c(LETTERS[ which(letters%in%substr(names(noise),1,1)) ]) }
  if(is.null(names(params))){names(params)<-c(LETTERS[1:length(params)])}
  if(any(names(params)%in%letters)){ names(params)<-c(LETTERS[ which(letters%in%substr(names(params),1,1)) ]) }
  if(wasNULL & any(grepl("fixedChangePoint", names(noise)))){
    noise[grepl("fixedChangePoint", names(noise))]<-0
  }
  
  #* check that params are all the same length, if not then rep until they are
  if(!all(unlist(lapply(params,length))== max(unlist(lapply(params,length))))){
    warning("params are not uniform length, values are being recycled to fit max length")
    diffLengths<-which(!unlist(lapply(params,length))== max(unlist(lapply(params,length))))
    params[diffLengths]<-lapply(diffLengths, function(i) rep(params[[i]], length.out =max(unlist(lapply(params,length)))  ))
  }
  #* decide which internal funciton to use
  if(!grepl("\\+", model)){
    out <- .singleGrowthSim(model, n, t, params, noise)
  } else{
    out <- .multiGrowthSim(model, n, t, params, noise)
  }
  return(out)
}

#' Internal helper function to simulate growth from a series of growth models
#' @keywords internal
#' @noRd

.multiGrowthSim<- function(model, n=20, t=25, params=list(), noise=NULL){
  
  component_models <- trimws(strsplit(model, "\\+")[[1]])
  
  firstModel <- component_models[1]
  firstModel <- sub("\\s", "", firstModel)
  firstParams <- params[grepl(paste0(firstModel,"1"), names(params))]
  firstChangepoints <- params[["changePoint1"]]
  firstNoise <- noise[grepl(paste0(firstModel,"1|changePoint1"), names(noise))]
  names(firstNoise)<- sub(paste0(firstModel,"1|Point."), "", names(firstNoise))
  
  if(is.null(firstChangepoints)){stop("Simulating segmented growth models requires 'changePointX' parameters as described in growthSS.")}
  
  df1<-do.call(rbind, lapply(1:n, function(i){
    
    firstChangepoints_rand <- lapply(firstChangepoints, function(fc){
      round(rnorm(1, fc, firstNoise$change))
    })
    
    n_df <- do.call(rbind, lapply(1:length(firstChangepoints_rand), function(g){
      .singleGrowthSim(firstModel, n=1, t=firstChangepoints_rand[[g]],
                       params= stats::setNames(lapply(firstParams, function(l) l[[g]]), c(sub(paste0(firstModel,"1"), "", names(firstParams)))),
                       noise = firstNoise)
    }))
    n_df$group <- rep(letters[1:length(firstChangepoints_rand)], times = unlist(firstChangepoints_rand))
    n_df$id <- paste0("id_",i)
    n_df
  }))
  dataList <- list(df1)
  
  for(u in 2:length(component_models)){
    
    iterModel <- component_models[u]
    iterModel <- sub("\\s", "", iterModel)
    iterParams <- params[grepl(paste0(iterModel,u), names(params))]
    
    nextChangepoints <- params[[paste0("changePoint",u)]]
    iterNoise <- noise[grepl(paste0(iterModel, u,"|changePoint",u), names(noise))]
    names(iterNoise)<- sub(paste0(iterModel, u,"|Point."), "", names(iterNoise))
    
    iter_data<-do.call(rbind, lapply(1:n, function(i){
      if(is.null(nextChangepoints) | u==length(component_models) ){
        iterChangepoints_rand <- rep(t, length(iterParams[[1]]))
      } else{
        iterChangepoints_rand <- lapply(nextChangepoints, function(fc){
          round(rnorm(1, fc, iterNoise$change))
        })
      }
      n_df <- do.call(rbind, lapply(1:length(iterChangepoints_rand), function(g){
        if(u==length(component_models)){ gt = t - max(df1[df1$id==paste0("id_",i), "time"])
        }else{gt = iterChangepoints_rand[[g]]}
        
        inner_df <- .singleGrowthSim(iterModel, n=1, t= gt ,
                         params= stats::setNames(lapply(iterParams, function(l) l[[g]]), c(sub(paste0(iterModel,u), "", names(iterParams)))),
                         noise=iterNoise )
        inner_df$group <- letters[g]
        inner_df
      }))
      n_df$id <- paste0("id_",i)
      n_df
    }))
    
    prev_data <- dataList[[(u-1)]]
    
    new_data <- do.call(rbind, lapply(unique(paste0(iter_data$id, iter_data$group)), function(int){
      prev_data_sub<-prev_data[paste0(prev_data$id, prev_data$group)==int,]
      iter_data_sub <- iter_data[paste0(iter_data$id, iter_data$group) == int, ]
      y_end <- prev_data_sub[prev_data_sub$time == max(prev_data_sub$time), "y"]
      iter_data_sub$time <- iter_data_sub$time + max(prev_data_sub$time)
      iter_data_sub$y <- y_end + iter_data_sub$y
      iter_data_sub
    }))
    dataList[[(u)]] <- new_data
  }
  out <- do.call(rbind, dataList)
  out <- out[out$time < t,]
  return(out)
}




#' Internal helper function to simulate growth from a single growth model
#' @keywords internal
#' @noRd

.singleGrowthSim<- function(model, n=20, t=25, params=list(), noise=NULL){
  
  models<-c("logistic", "gompertz", "double logistic", "double gompertz",
            "monomolecular", "exponential", "linear", "power law")
  
  matched_model <- match.arg(model, models)
  
  if(matched_model=="logistic"){
    gsi<-gsi_logistic
  } else if (matched_model =="gompertz"){
    gsi<-gsi_gompertz
  } else if (matched_model=="double logistic"){
    gsi<-gsi_dou_logistic
  } else if (matched_model=="double gompertz"){
    gsi<-gsi_dou_gompertz
  } else if (matched_model=="monomolecular"){
    gsi<-gsi_monomolecular
  } else if (matched_model=="exponential"){
    gsi<-gsi_exponential
  } else if (matched_model=="linear"){
    gsi<-gsi_linear
  } else if (matched_model=="power law"){
    gsi<-gsi_powerlaw
  }
  
  out<-do.call(rbind,lapply(1:length(params[[1]]), function(i) {
    pars<-lapply(params, function(p) p[i])
    as.data.frame(rbind(do.call(rbind,lapply(1:n,function(e){
      data.frame("id"=paste0("id_",e),"group"=letters[i],"time"=1:t,"y"=gsi(1:t, pars, noise),stringsAsFactors = FALSE)}
    ))))
  }))
  
  return(out)
}


#* ************************************************************
#* ***** `gsi functions to simulate individual plants` *****
#* ************************************************************

gsi_logistic <- function(x,pars, noise){
  a_r <- pars[["A"]]+rnorm(1,mean = 0,sd=noise[["A"]])
  b_r <- pars[["B"]]+rnorm(1,mean=0,sd=noise[["B"]])
  c_r <- pars[["C"]]+rnorm(1,mean=0,sd=noise[["C"]])
  return(a_r / (1 + exp( (b_r-x)/c_r) ))
}
gsi_gompertz <- function(x,pars, noise){
  a_r <- pars[["A"]]+rnorm(1,mean = 0,sd=noise[["A"]])
  b_r <- pars[["B"]]+rnorm(1,mean=0,sd=noise[["B"]])
  c_r <- pars[["C"]]+rnorm(1,mean=0,sd=noise[["C"]])
  return(a_r*exp(-b_r*exp(-c_r*x)))
}
gsi_dou_logistic <- function(x,pars, noise){
  a_r <- pars[["A"]]+rnorm(1,mean = 0,sd=noise[["A"]])
  b_r <- pars[["B"]]+rnorm(1,mean=0,sd=noise[["B"]])
  c_r <- pars[["C"]]+rnorm(1,mean=0,sd=noise[["C"]])
  a2_r <- pars[["A2"]]+rnorm(1,mean = 0,sd=noise[["A2"]])
  b2_r <- pars[["B2"]]+rnorm(1,mean=0,sd=noise[["B2"]])
  c2_r <- pars[["C2"]]+rnorm(1,mean=0,sd=noise[["C2"]])
  return( a_r/(1+exp((b_r-x)/c_r)) + ((a2_r-a_r) /(1+exp((b2_r-x)/c2_r))) )
}
gsi_dou_gompertz <- function(x,pars, noise){
  a_r <- pars[["A"]]+rnorm(1,mean = 0,sd=noise[["A"]])
  b_r <- pars[["B"]]+rnorm(1,mean=0,sd=noise[["B"]])
  c_r <- pars[["C"]]+rnorm(1,mean=0,sd=noise[["C"]])
  a2_r <- pars[["A2"]]+rnorm(1,mean = 0,sd=noise[["A2"]])
  b2_r <- pars[["B2"]]+rnorm(1,mean=0,sd=noise[["B2"]])
  c2_r <- pars[["C2"]]+rnorm(1,mean=0,sd=noise[["C2"]])
  return( (a_r * exp(-b_r * exp(-c_r*x))) + ( (a2_r-a_r) * exp(-b2_r * exp(-c2_r*(x-b_r))) ) )
}
gsi_monomolecular <- function(x,pars, noise){
  a_r <- pars[["A"]]+rnorm(1,mean = 0,sd=noise[["A"]])
  b_r <- pars[["B"]]+rnorm(1,mean=0,sd=noise[["B"]])
  return(a_r-a_r*exp(-b_r*x))
}
gsi_exponential <- function(x,pars, noise){
  a_r <- pars[["A"]]+rnorm(1,mean = 0,sd=noise[["A"]])
  b_r <- pars[["B"]]+rnorm(1,mean=0,sd=noise[["B"]])
  return(a_r*exp(b_r * x))
}
gsi_linear <- function(x,pars, noise){
  a_r <- pars[["A"]]+rnorm(1,mean = 0,sd=noise[["A"]])
  return(a_r*x)
}
gsi_powerlaw <- function(x,pars, noise){
  a_r <- pars[["A"]]+rnorm(1,mean = 0,sd=noise[["A"]])
  b_r <- pars[["B"]]+rnorm(1,mean=0,sd=noise[["B"]])
  return(a_r * x^(b_r))
}

