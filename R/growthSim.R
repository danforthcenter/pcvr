#' Growth data simulating function
#' 
#' @description growthSim can be used to help pick reasonable parameters for common growth models to use in prior distributions or to simulate data for example models/plots.
#' 
#' @param model One of "logistic", "gompertz", "monomolecular", "exponential", "linear", or "power law" 
#' @param n Number of individuals to simulate over time per each group in params
#' @param t Max time (assumed to start at T1) to simulate growth to.
#' @param params A list of numeric parameters. A, B, C notation is used in the order that parameters appear in the formula (see examples). Number of groups is inferred from the length of these vectors of parameters.
#' @param noise Optionally this can be used to add specific amounts of noise to the input parameters. If NULL (the default) then data is simulated with 10\% random noise like: param + N(0, 0.1*param)
#' 
#' @keywords growth curve, logistic, gompertz, monomolecular, linear, exponential, power-law
#' @return Returns a dataframe of example growth data following the input parameters.
#' 
#' @importFrom stats rnorm
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
    as.data.frame(rbind(do.call(rbind,lapply(1:n,function(e) data.frame("id"=paste0("id_",e),"group"=letters[i],"time"=1:t,"y"=gsi(1:t, pars, noise),stringsAsFactors = FALSE)))))
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

