#' Ease of use nlrq starter function for 6 growth model parameterizations
#' 
#' Internal to growthSS
#' 
#' @examples 
#' 
#' simdf<-growthSim("logistic", n=20, t=25,
#'   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' 
#' ss<-.nlrqSS(model = "logistic", form=y~time|id/group,
#'   tau=0.5, df=simdf, start=NULL)
#'   
#' ss<-.nlrqSS(model = "gam", form=y~time|id/group, df=simdf, start=NULL, tau=0.5)
#'   
#' dim(ss$df)
#' ss[c("formula", "taus", "start", "pcvrForm")]
#' @importFrom splines bs
#' @keywords internal
#' @noRd

.nlrqSS<-function(model, form, tau=0.5, df, pars=NULL, start=NULL, type="nlrq"){
  #* ***** `Define choices and make empty output list`
  out<-list()
  models<-c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law",
            "double logistic", "double gompertz", "gam", "frechet", "weibull", "gumbel")
  #* ***** `Make nlrq formula` *****
  #* `parse form argument`
  y=as.character(form)[2]
  x<-as.character(form)[3]
  group = NULL
  if(grepl("\\|", x) & grepl("\\/",x)){
    x3<-trimws(strsplit(x, "[|]|[/]")[[1]])
    x<-x3[1]
    individual = x3[2]
    group = x3[length(x3)] 
    message(paste0("Individual is not used with type = '",type,"'."))
    if(length(unique(df[[group]]))==1){USEGROUP=FALSE}else{USEGROUP=TRUE} # if there is only one group then ignore grouping for parameter and variance formulas
  } else if (grepl("\\|", x)){
    x2<-trimws(strsplit(x, "[|]")[[1]])
    x<-x2[1]
    group = x2[2] # individual will not be used here 
    USEGROUP=TRUE # leave x as is, in this parameterization this means no term[group] syntax
    message(paste0("Individual is not used with type = '",type,"', interpreting ", group, ", as a group."))
  } else {
    # leave x as is, in this parameterization this means no term[group] syntax
    USEGROUP=FALSE
  }
  if(USEGROUP){
    df[[group]]<-factor(df[[group]])
    df[[paste0(group,"_numericLabel")]]<-unclass(df[[group]])
  }
  #* `assemble growth formula`
  if(grepl("decay", model)){
    decay=TRUE
    model <- trimws(gsub("decay", "", model))
  } else{
    decay=FALSE
  }
  
  matched_model <- match.arg(model, models)
  string_formFun <- paste0(".nlrq_form_", gsub(" ", "", matched_model))
  form_fun <- match.fun(string_formFun)
  res = form_fun(x,y, USEGROUP, group, pars)
  growthForm <- res[[1]]
  pars <- res[[2]]
  
  if(decay){ growthFrom <- .nlrq_Decay(growthForm) }
  
  if(matched_model =="gam"){start<-0}
  
  if(is.null(start)){
    if(matched_model=="double logistic"){
      warning("Double Sigmoid models are not supported as self-starting models, you will need to add starting parameters. Note for these models type='brms' is recommended.")
      startingValues<-NULL
      } else if (matched_model=="double gompertz"){
      warning("Double Sigmoid models are not supported as self-starting models, you will need to add starting parameters. Note for these models type='brms' is recommended.")
      startingValues<-NULL
      } else{
        string_initFun <- paste0(".init", gsub(" ", "", matched_model))
        initFunction <- match.fun(string_initFun)
        startingValues <- initFunction(df, x, y)
    }
    if((!matched_model %in% c("double logistic", "double gompertz")) & USEGROUP){
      nms <- names(startingValues)
      startingValuesList <- lapply(names(startingValues), function(nm){
        val<-startingValues[nm]
        if(nm %in% pars){
          rep(val, length(unique(df[[group]])) ) # if this is one of pars then make starting value per group
        } else{ val }# else return one starting value
        })
      names(startingValuesList)<-nms
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
# .initlogistic(simdf, "time", "y")

.initlogistic <- function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if(nrow(xy) < 4)
    stop("too few distinct input values to fit a logistic model")
  z <- abs(xy[["y"]])
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
# .initgompertz(simdf, "time", "y")

.initgompertz<-function(df, x, y) {
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
# .initmonomolecular(ex, "time", "y")


.initmonomolecular <- function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if(nrow(xy) < 4)
    stop("too few distinct input values to fit a monomolecular model")
  z <- abs(xy[["y"]])
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
#  .initlinear(ex, "time", "y")

.initlinear <- function(df, x, y) {
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

.initpowerlaw <- function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if(nrow(xy) < 3)
    stop("too few distinct input values to fit a power law model")
  aux <- stats::coef(stats::lm(y ~ x, xy))
  
  pars <- stats::coef(stats::nls(y ~ 1 * x^B, data = xy, start = list( B = aux[2L] ), 
                   algorithm = "plinear", control = stats::nls.control(warnOnly=TRUE)))
  stats::setNames(pars[c(".lin", "B")], c("A", "B"))
}


#* `exponential self starter`
# ex<-growthSim("exponential", n=20, t=25,
#                  params = list("A"=c(15, 20), "B"=c(0.095, 0.095)))
# .initExponential(ex, "time","y")

.initexponential <- function(df, x, y) {
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if(nrow(xy) < 3)
    stop("too few distinct input values to fit a exponential model")
  aux <- stats::coef(stats::lm(y ~ x, xy))
  
  pars <- stats::coef(stats::nls(y ~ 1 * exp(B*x), data = xy, start = list( B = aux[2L] ), 
                   algorithm = "plinear", control = stats::nls.control(warnOnly=TRUE) ) )
  stats::setNames(pars[c(".lin", "B")],c("A", "B"))
}

#* `Extreme Value Distribution Self Starter (weibull, frechet, gumbel)`

.initweibull <- function(df, x, y){
    xy <- sortedXyData(df[[x]], df[[y]])
    if (nrow(xy) < 5) {
        stop("too few distinct input values to fit the EVD growth model")
    }
    if (any(xy[["x"]] < 0)) {
        stop("all 'x' values must be non-negative to fit the EVD growth model")
    }
    Rasym <- stats::NLSstRtAsymptote(xy)
    Lasym <- stats::NLSstLfAsymptote(xy)
    pars <- stats::coef(stats::lm(log(-log((Rasym - y)/(Rasym - Lasym))) ~ 
        log(x), data = xy, subset = x > 0))
    # coef(nls(y ~ cbind(1, -exp(-exp(lrc) * x^pwr)), 
    #                   data = xy, start = c(lrc = pars[[1L]], pwr = pars[[2L]]), 
    #                   algorithm = "plinear"))
    stats::setNames(c(Rasym, exp(pars)+c(1,0) ), c("A", "B", "C"))
}
.initfrechet <- .initweibull
.initgumbel <- .initweibull

#* `Define growth formulas`

.nlrq_form_logistic<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B", "C")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y," ~ A[]/(1+exp((B[]-",x,")/C[]))")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf<-as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C))"))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_gompertz<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B", "C")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y," ~ A[]*exp(-B[]*exp(-C[]*",x,"))")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf<-as.formula(paste0(y," ~ A*exp(-B*exp(-C*",x,"))"))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_doublelogistic<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B", "C", "A2", "B2", "C2")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y," ~ A[]/(1+exp((B[]-",x,")/C[]))",
                     " + ((A2[]-A[]) /(1+exp((B2[]-",x,")/C2[])))")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf <- as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C)) + ((A2-A) /(1+exp((B2-",x,")/C2)))"))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_doublegompertz<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B", "C", "A2", "B2", "C2")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y," ~ A[] * exp(-B[] * exp(-C[]*",x,"))",
                     " + (A2[]-A[]) * exp(-B2[] * exp(-C2[]*(",x,"-B[])))")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf <- as.formula(paste0(y," ~ A * exp(-B * exp(-C*",x,")) + (A2-A) * exp(-B2 * exp(-C2*(",x,"-B)))"))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_monomolecular<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y,"~A[]-A[]*exp(-B[]*",x,")")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf <- as.formula(paste0(y,"~A-A*exp(-B*",x,")"))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_exponential<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y," ~ A[]*exp(B[]*",x, ")")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf <- as.formula(paste0(y," ~ A*exp(B*",x, ")"))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_linear<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y," ~ A[]*",x)
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf <- as.formula(paste0(y," ~ A*",x))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_powerlaw<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y," ~ A[]*",x,"^B[]")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf <- as.formula(paste0(y," ~ A*",x,"^B"))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_gam<-function(x, y, USEGROUP, group, pars){
  if(USEGROUP){
    nf<-as.formula(paste0(y, " ~ bs(",x,")*",group))
  } else{
    nf <- as.formula(paste0(y, " ~ bs(",x,")"))
  }
  return(list("formula" = nf, "pars" = NULL))
}

.nlrq_Decay<-function(form){
  chars <- as.character(form)
  as.formula(paste0(chars[2], chars[1], "-(", chars[3],")" ))
}

.nlrq_form_frechet<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B", "C")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y, " ~ A[] * exp(-((",x,"-0)/C[])^(-B[]))")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf<-as.formula(paste0(y," ~ A * exp(-((",x,"-0)/C)^(-B))"))
  }
  return(list("formula" = nf, "pars" = pars))
}


.nlrq_form_weibull<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B", "C")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y, " ~ A[] * (1-exp(-(",x,"/C[])^B[]))")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf<-as.formula(paste0(y," ~ A * (1-exp(-(",x,"/C)^B))"))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_gumbel<-function(x, y, USEGROUP, group, pars){
  total_pars = c("A", "B", "C")
  if(is.null(pars)){pars <- total_pars}
  if(USEGROUP){
    str_nf <- paste0(y, " ~ A[] * exp(-exp(-(",x,"-B[])/C[]))")
    for(par in total_pars){
      if(par %in% pars){
        str_nf<-gsub(paste0(par, "\\[\\]"), paste0(par, "[",group,"]"), str_nf)
      } else{
        str_nf<-gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf<-as.formula(str_nf)
  } else{
    nf<-as.formula(paste0(y," ~ A * exp(-exp(-(",x,"-B)/C))"))
  }
  return(list("formula" = nf, "pars" = pars))
}


