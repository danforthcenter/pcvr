#' Ease of use nlrq starter function for 6 growth model parameterizations
#' 
#' Internal to growthSS
#' 
#' @examples 
#' 
#' simdf<-growthSim("logistic", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' 
#' 
#' ss<-growthSS(model = "logistic", form=y~time|id/group,
#' sigma="spline", df=simdf, priors = list("A"=130, "B"=12, "C"=3), type="nlrq")
#' lapply(ss,class)
#' ss$initfun()
#' 
#' @keywords internal
#' @noRd

.nlrqSS<-function(model, form, tau, df, priors=NULL){
  
  model = "logistic"; form=y~time|id/group
  sigma="spline"; df=simdf; priors = list("A"=130, "B"=12, "C"=3)
  #* ***** `Define choices and make empty output list`
  out<-list()
  models<-c("logistic", "gompertz", "monomolecular",
            "exponential", "linear", "power law",
            "double logistic", "double gompertz")
  #* ***** `Make nlrq formula` *****
  #* `parse form argument`
  y=as.character(form)[2]
  x<-as.character(form)[3]
  group = NULL
  if(grepl("\\|", x) & grepl("\\/",x)){
    x3<-trimws(strsplit(x, "[|]|[/]")[[1]])
    x<-x3[1]
    individual = x3[2]
    group = x3[3] 
    message("Individual is not used with type = 'nlrq'.")
    if(length(unique(df[[group]]))==1){USEGROUP=FALSE}else{USEGROUP=TRUE} # if there is only one group then ignore grouping for parameter and variance formulas
  } else if (grepl("\\|", x)){
    x2<-trimws(strsplit(x, "[|]")[[1]])
    x<-x2[1]
    group = x2[2] # individual will not be used here 
    USEGROUP=FALSE # leave x as is, in this parameterization this means no term[group] syntax
    message(paste0("Individual is not used with type = 'nlrq', interpreting ", group, ", as a group."))
  } else {
    # leave x as is, in this parameterization this means no term[group] syntax
    USEGROUP=FALSE
  }
  #* `assemble growth formula`
  matched_model <- match.arg(model, models)
  if(matched_model=="double logistic"){
    form_fun<-.nlrq_form_dou_logistic
  } else if (matched_model=="double gompertz"){
    form_fun<-.nlrq_form_dou_gompertz
  } else if(matched_model=="logistic"){
    form_fun<-.nlrq_form_logistic
    #getInitial(y ~ SSlogis(time, Asym, xmid, scal), data = simdf)
  } else if (matched_model=="gompertz"){
    form_fun<-.nlrq_form_gompertz
    #getInitial(y ~ SSlogis(time, Asym, xmid, scal), data = simdf)
  } else if (matched_model=="monomolecular"){
    form_fun<-.nlrq_form_monomolecular
  } else if (matched_model=="exponential"){
    form_fun<-.nlrq_form_exponential
  } else if (matched_model=="linear"){
    form_fun<-.nlrq_form_linear
  } else if (matched_model=="power law"){
    form_fun<-.nlrq_form_powerlaw
  }
  growthForm = form_fun(x,y, USEGROUP, group)
  # still need to turn this into a self starting formula
  
  out[["formula"]] <- growthForm
  #* nlrq actually has so few arguments that you have to specify compared to brms that this is sort of odd to have.
  
  out[["taus"]] <- tau
  
  #* if I didn't like the selfStart function I could also infer starting values 
  #* some other way and return a list of starting values? That might be difficult to do well?
  #* I think asymptote is easy, then put the inflection towards the middle,
  #* then take a linear slope/some constant per model
  #* as an idea of the growth rate? The growth rate is the hard part.
  
  
  
  
}

#* example of using selfStart

initLogistic <- function(mCall, data, LHS, ...) {
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4)
    stop("too few distinct input values to fit a logistic model")
  z <- xy[["y"]]
  ## transform to proportion, i.e. in (0,1) :
  rng <- range(z); dz <- diff(rng)
  z <- (z - rng[1L] + 0.05 * dz)/(1.1 * dz)
  xy[["z"]] <- log(z/(1 - z))		# logit transformation
  aux <- coef(lm(x ~ z, xy))
  pars <- coef(nls(y ~ 1/(1 + exp((B - x)/C)),
                   data = xy,
                   start = list(B = aux[[1L]], C = aux[[2L]]),
                   algorithm = "plinear", ...))
  setNames(pars [c(".lin", "B", "C")],
           mCall[c("A", "B", "C")])
}

.nlrq_SSlogistic <- selfStart(~ A/(1 + exp((B - x)/C)),
                       initial = initLogistic,
                       parameters = c("A", "B", "C"))

getInitial(y ~ .nlrq_SSlogistic(time, A, B, C), data = simdf)







initGompertz <- function(mCall, data, LHS, ...) {
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4)
    stop("too few distinct input values to fit a gompertz model")
  z <- xy[["y"]]
  ## transform to proportion, i.e. in (0,1) :
  rng <- range(z); dz <- diff(rng)
  z <- (z - rng[1L] + 0.05 * dz)/(1.1 * dz)
  xy[["z"]] <- log(z/(1 - z))		# logit transformation
  aux <- coef(lm(x ~ z, xy))
  pars <- coef(nls(y ~ 1 * exp(-B * exp(-C * x)),
                   data = xy,
                   start = list(B = aux[[1L]], C = aux[[2L]]),
                   algorithm = "plinear", ...))
  setNames(pars [c(".lin", "B", "C")],
           mCall[c("A", "B", "C")])
}

.nlrq_SSGompertz <- selfStart(~ A * exp(-B * exp(-C * x)),
                              initial = initGompertz,
                              parameters = c("A", "B", "C"))

getInitial(y ~ .nlrq_SSGompertz(time, A, B, C), data = simdf)











.nlrq_form_logistic<-function(x, y, USEGROUP, group){
  if(USEGROUP){
    nf<-as.formula(paste0(y," ~ A[",group,"]/(1+exp((B[",group,"]-",x,")/C[",group,"]))"))
  } else{
    nf<-as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C))"))
  }
  return(nf)
}

.nlrq_form_gompertz<-function(x, y, USEGROUP, group){
  if(USEGROUP){
    nf<-as.formula(paste0(y," ~ A[",group,"]*exp(-B[",group,"]*exp(-C[",group,"]*",x,"))"))
  } else{
    nf<-as.formula(paste0(y," ~ A*exp(-B*exp(-C*",x,"))"))
  }
  return(nf)
}

.nlrq_form_dou_logistic<-function(x, y, USEGROUP, group){
  if(USEGROUP){
    nf <- as.formula(paste0(y," ~ A[",group,"]/(1+exp((B[",group,"]-",x,")/C[",group,"]))",
                            " + ((A2[",group,"]-A[",group,"]) /(1+exp((B2[",group,"]-",x,")/C2[",group,"])))"))
  } else{
    nf <- as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C)) + ((A2-A) /(1+exp((B2-",x,")/C2)))"))
  }
  return(nf)
}

.nlrq_form_dou_gompertz<-function(x, y, USEGROUP, group){
  if(USEGROUP){
    nf <- as.formula(paste0(y," ~ A[",group,"] * exp(-B[",group,"] * exp(-C[",group,"]*",x,"))",
          " + (A2[",group,"]-A[",group,"]) * exp(-B2[",group,"] * exp(-C2[",group,"]*(",x,"-B[",group,"])))"))
  } else{
    nf <- as.formula(paste0(y," ~ A * exp(-B * exp(-C*",x,")) + (A2-A) * exp(-B2 * exp(-C2*(",x,"-B)))"))
  }
  return(nf)
}

.nlrq_form_dou_monomolecular<-function(x, y, USEGROUP, group){
  if(USEGROUP){
    nf <- as.formula(paste0(y,"~A[",group,"]-A[",group,"]*exp(-B[",group,"]*",x,")"))
    } else{
    nf <- as.formula(paste0(y,"~A-A*exp(-B*",x,")"))
  }
  return(nf)
}

.nlrq_form_dou_exponential<-function(x, y, USEGROUP, group){
  if(USEGROUP){
    nf <- as.formula(paste0(y," ~ A[",group,"]*exp(B[",group,"]*",x, ")"))
  } else{
    nf <- as.formula(paste0(y," ~ A*exp(B*",x, ")"))
  }
  return(nf)
}

.nlrq_form_dou_linear<-function(x, y, USEGROUP, group){
  if(USEGROUP){
    nf <- as.formula(paste0(y," ~ A[",group,"]*",x))
  } else{
    nf <- as.formula(paste0(y," ~ A*",x))
  }
  return(nf)
}

.nlrq_form_dou_powerlaw<-function(x, y, USEGROUP, group){
  if(USEGROUP){
    nf <- as.formula(paste0(y," ~ A[",group,"]*",x,"^B[",group,"]"))
  } else{
    nf <- as.formula(paste0(y," ~ A*",x,"^B"))
  }
  return(nf)
}









