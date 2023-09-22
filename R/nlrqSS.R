#' Ease of use nlrq starter function for 6 growth model parameterizations
#' 
#' Internal to growthSS
#' 
#' @examples 
#' 
#' simdf<-growthSim("logistic", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' 
#' ss<-.nlrqSS(model = "logistic", form=y~time|id/group,
#'   tau=0.5, df=simdf, start=NULL)
#' dim(ss$df)
#' ss[c("formula", "taus", "start", "pcvrForm")]
#' 
#' @keywords internal
#' @noRd

.nlrqSS<-function(model, form, tau, df, start=NULL, type="nlrq"){
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
    message(paste0("Individual is not used with type = '",type,"'."))
    if(length(unique(df[[group]]))==1){USEGROUP=FALSE}else{USEGROUP=TRUE} # if there is only one group then ignore grouping for parameter and variance formulas
  } else if (grepl("\\|", x)){
    x2<-trimws(strsplit(x, "[|]")[[1]])
    x<-x2[1]
    group = x2[2] # individual will not be used here 
    USEGROUP=FALSE # leave x as is, in this parameterization this means no term[group] syntax
    message(paste0("Individual is not used with type = '",type,"', interpreting ", group, ", as a group."))
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
  
  if(is.null(start)){
    if(matched_model=="double logistic"){
      warning("Double Sigmoid models are not supported as self-starting models, you will need to add starting parameters")
      startingValues<-NULL
      } else if (matched_model=="double gompertz"){
      warning("Double Sigmoid models are not supported as self-starting models, you will need to add starting parameters")
      startingValues<-NULL
    } else if(matched_model=="logistic"){
      startingValues <- .initLogistic(df, x, y)
    } else if (matched_model=="gompertz"){
      startingValues <- .initGompertz(df, x, y)
    } else if (matched_model=="monomolecular"){
      startingValues <- .initMonomolecular(df, x, y)
    } else if (matched_model=="exponential"){
      startingValues <- .initExponential(df,x,y)
    } else if (matched_model=="linear"){
      startingValues <- .initLinear(df,x,y)
    } else if (matched_model=="power law"){
      startingValues <- .initPowerLaw(df,x,y)
    }
    if((!matched_model %in% c("double logistic", "double gompertz")) & USEGROUP){
      startingValuesList <- lapply(startingValues, function(i) rep(i, length(unique(df[[group]])) ))
    } else{ # non-grouped, just make it into a list
      startingValuesList <- as.list(startingValues)
    }
  } else{
    startingValuesList <- start
  }
  out[["formula"]] <- growthForm
  if(type == "nlrq"){ out[["taus"]] <- tau }
  out[["start"]] <- startingValuesList
  out[["df"]]<-df
  out[["pcvrForm"]]<-form
  return(out)
}

#* example of using selfStart
#* `logistic self starter`
# simdf<-growthSim("logistic", n=20, t=25,
#   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
# .initLogistic(simdf, "time", "y")

.initLogistic <- function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if(nrow(xy) < 4)
    stop("too few distinct input values to fit a logistic model")
  z <- xy[["y"]]
  ## transform to proportion, i.e. in (0,1) :
  rng <- range(z); dz <- diff(rng)
  z <- (z - rng[1L] + 0.05 * dz)/(1.1 * dz)
  xy[["z"]] <- log(z/(1 - z))		# logit transformation
  aux <- stats::coef(stats::lm(x ~ z, xy))
  pars <- stats::coef(stats::nls(y ~ 1/(1 + exp((B - x)/C)),
                   data = xy,
                   start = list(B = aux[[1L]], C = aux[[2L]]),
                   algorithm = "plinear", control = stats::nls.control(warnOnly=TRUE)))
  stats::setNames(pars [c(".lin", "B", "C")], c("A", "B", "C"))
}


#* `Goempertz self starter`
#* 
# .initGompertz(simdf, "time", "y")

.initGompertz<-function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if (nrow(xy) < 4) {
    stop("too few distinct input values to fit the Gompertz model")
  }
  xyL <- xy
  xyL$y <- log(abs(xyL$y))
  pars <- stats::NLSstAsymptotic(xyL)
  pars <- stats::coef(stats::nls(y ~ exp(-B * C^x), data = xy, start = c(B = pars[["b1"]], 
                                                             C = exp(-exp(pars[["lrc"]]))),
                                 algorithm = "plinear", control = stats::nls.control(warnOnly=TRUE)))
  stats::setNames(pars[c(".lin", "B", "C")], c("A", "B", "C"))
}


#* `Monomolecular self starter`
# ex<-growthSim("monomolecular", n=20, t=25,
#                  params = list("A"=c(200,160), "B"=c(0.08, 0.1)))
# .initMonomolecular(ex, "time", "y")


.initMonomolecular <- function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if(nrow(xy) < 4)
    stop("too few distinct input values to fit a monomolecular model")
  z <- xy[["y"]]
  ## transform to proportion, i.e. in (0,1) :
  rng <- range(z); dz <- diff(rng)
  z <- (z - rng[1L] + 0.05 * dz)/(1.1 * dz)
  xy[["z"]] <- z#log(z)		# log transformation
  aux <- stats::coef(stats::lm(z ~ x, xy))
  pars <- stats::coef(stats::nls(y ~ 1 * exp(B*x), data = xy, start = list( B = aux[2L] ), 
                   algorithm = "plinear", control = stats::nls.control(warnOnly=TRUE)))
  stats::setNames(pars[c(".lin", "B")], c("A", "B"))
}

#* `Linear self starter`

# ex<-growthSim("linear", n=20, t=25,
#                  params = list("A"=c(1.1, 0.95)))
#  .initLinear(ex, "time", "y")

.initLinear <- function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if(nrow(xy) < 2)
    stop("too few distinct input values to fit a linear model")
  pars <- stats::coef(stats::lm(y~x, xy))
  stats::setNames(pars[c("x")], c("A"))
}


#* `power law self starter`

# ex<-growthSim("power law", n=20, t=25,
#                  params = list("A"=c(16, 11), "B"=c(0.75, 0.7)))
# .initPowerLaw(df, "time", "y")

.initPowerLaw <- function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if(nrow(xy) < 3)
    stop("too few distinct input values to fit a power law model")
  aux <- stats::coef(stats::lm(y ~ x, xy))
  
  pars <- stats::coef(stats::nls(y ~ 1 * x^B, data = xy, start = list( B = aux[2L] ), 
                   algorithm = "plinear", control = stats::nls.control(warnOnly=TRUE)))
  stats::setNames(pars[c(".lin", "B")], c("A", "B"))
}


#* `power law self starter`
# ex<-growthSim("exponential", n=20, t=25,
#                  params = list("A"=c(15, 20), "B"=c(0.095, 0.095)))
# .initExponential2(ex, "time","y")

.initExponential <- function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if(nrow(xy) < 3)
    stop("too few distinct input values to fit a exponential model")
  aux <- stats::coef(stats::lm(y ~ x, xy))
  
  pars <- stats::coef(stats::nls(y ~ 1 * exp(B*x), data = xy, start = list( B = aux[2L] ), 
                   algorithm = "plinear", control = stats::nls.control(warnOnly=TRUE) ) )
  stats::setNames(pars[c(".lin", "B")],c("A", "B"))
}

#* `Define growth formulas`

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

