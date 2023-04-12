#' Growth data simulating function
#' 
#' @description growthSim can be used to help pick reasonable parameters for common growth models to use in prior distributions.
#' 
#' maybe I should make a data simulating function for each supported distribution as well? 
#' something that takes the parameters, the number of samples per group, and returns a dataframe?
#' 
#' 
#' @param df Data frame to use. Can be wide or long format from read.pcv
#' @keywords growth curve, logistic, gompertz, monomolecular, linear, exponential, power law
#' @examples 
#' simdf<-growthSim("logistic", n=20, t=25, params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ggplot(simdf,aes(time, y, group=interaction(group,id)))+ geom_line(aes(color=group))
#' 
#' 
#' @export
#' 
#' 

growthSim<-function(model=c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law"), n=20, t=25, params=list(), noise=NULL){
  if(length(model)>1){stop("Select one model to use")}
  
  if(is.null(noise)){noise = lapply(params, function(i) mean(i)/10)}
  if(is.null(names(noise))){names(noise)<-c(LETTERS[1:length(noise)])}
  if(any(names(noise)%in%letters)){ names(noise)<-c(LETTERS[ which(letters%in%substr(names(noise),1,1)) ]) }
  
  if(is.null(names(params))){names(params)<-c(LETTERS[1:length(params)])}
  if(any(names(params)%in%letters)){ names(params)<-c(LETTERS[ which(letters%in%substr(names(params),1,1)) ]) }
  
  models<-c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law")
  if(match.arg(model, models)=="logistic"){
    gsi<-gsi_logistic
  } else if (match.arg(model, models)=="gompertz"){
    gsi<-gsi_gompertz
  } else if (match.arg(model, models)=="monomolecular"){
    gsi<-gsi_monomolecular
  } else if (match.arg(model, models)=="exponential"){
    gsi<-gsi_exponential
  } else if (match.arg(model, models)=="linear"){
    gsi<-gsi_linear
  } else if (match.arg(model, models)=="power law"){
    gsi<-gsi_powerlaw
  }
  #* check that params are all the same length, if not then rep until they are
  if(!all(unlist(lapply(params,length))== max(unlist(lapply(params,length))))){
    warning("params are not uniform length, values are being recycled to fit max length")
    diffLengths<-which(!unlist(lapply(params,length))== max(unlist(lapply(params,length))))
    params[diffLengths]<-lapply(diffLengths, function(i) rep(params[[i]], length.out =max(unlist(lapply(params,length)))  ))
  }
  
  out<-do.call(rbind,lapply(1:length(params[[1]]), function(i) {
    pars<-lapply(params, function(p) p[i])
    as.data.frame(rbind(do.call(rbind,lapply(1:n,function(e) data.frame("id"=paste0("id_",e),"group"=letters[i],"time"=1:t,"y"=gsi(1:t, pars, noise),stringsAsFactors = F)))))
  }))
  return(out)
}

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

