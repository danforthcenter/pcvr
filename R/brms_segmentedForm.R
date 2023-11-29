#' Helper function to put together formulae for brms changepoint growth models
#' 
#' @param model A multi-part model passed from brmSS passed from \code{\link{growthSS}}
#' @param x The x variable from the pcvrForm argument in \code{\link{growthSS}}
#' @param y The y variable from the pcvrForm argument in \code{\link{growthSS}}
#' @param group The grouping variable from the pcvrForm argument in \code{\link{growthSS}}
#' @param nTimes a Number of times that are present in the data, only used for making splines have a workable number of knots.
#' @param useGroup logical, should groups be used?
#' @param priors A list describing priors in the style of \code{\link{brmSS}}, \code{\link{growthSS}}, and \code{\link{growthSim}}.
#' This is only used currently to identify fixed and estimated changepoints. If a changepoint is called "changePointX" with X being its
#' position in the formula then it will be estimated as a parameter in the model, but if the changepoint is called "fixedChangePointX"
#' then it will be passed as a numeric in the growth model. 
#'
#' @examples
#' df1<-do.call(rbind, lapply(1:30, function(i){
#' chngpt <- rnorm(2, 10, 1.5)
#' A<-growthSim("linear", n=1, t=chngpt[1], params=list("A"=c(1)))
#' B<-growthSim("linear", n=1, t=chngpt[2], params=list("A"=c(0.9)))
#' B$group <- "b"
#' x<-rbind(A, B)
#' x$id <- paste0("id_",i)
#' x
#' }))
#' df2 <- growthSim("linear", n=30, t=20, params=list("A"=c(4.1, 5)))
#' df2<-do.call(rbind, lapply(unique(paste0(df2$id, df2$group)), function(int){
#'   df1sub<-df1[paste0(df1$id, df1$group)==int,]
#'   df2sub <- df2[paste0(df2$id, df2$group) == int, ]
#'   y_end <- df1sub[df1sub$time == max(df1sub$time), "y"]
#'   df2sub$time <- df2sub$time + max(df1sub$time)
#'   df2sub$y <- y_end + df2sub$y
#'   df2sub
#' }))
#' df <- rbind(df1, df2)
#' ggplot(df, aes(time, y, group=interaction(group,id)))+ 
#'   geom_line(aes(color=group))+theme_minimal()
#' 
#' .brmsChangePointHelper(model = "linear + linear", x = "time", y="y", group="group")
#' 
#' @keywords internal
#' @noRd

.brmsChangePointHelper <- function(model, x, y, group, sigma=FALSE, nTimes=25, useGroup, priors){

  component_models <- trimws(strsplit(model, "\\+")[[1]])
  models<-c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law", "gam", "spline", "int", "homo")
  if(is.null(priors)){priors <- stats::setNames(lapply(1:length(component_models), identity ),
                                                paste0("changePoint",1:length(component_models))) }
  
  formulae <- lapply(1:length(component_models), function(i){
    iter_model <- component_models[i]
    matched_iter_model <- match.arg(iter_model, models)
    if(matched_iter_model=="logistic"){
      iter <- .logisticChngptForm(x, i, sigma)
    } else if (matched_iter_model=="gompertz"){
      iter <- .gompertzChngptForm(x, i, sigma)
    } else if (matched_iter_model=="monomolecular"){
      iter <- .monomolecularChngptForm(x, i, sigma)
    } else if (matched_iter_model=="exponential"){
      iter <- .exponentialChngptForm(x, i, sigma)
    } else if (matched_iter_model=="linear"){
      iter <- .linearChngptForm(x, i, sigma)
    } else if (matched_iter_model=="power law"){
      iter <- .powerLawChngptForm(x, i, sigma)
    } else if (matched_iter_model %in% c("gam", "spline")){
      iter <- .gamChngptForm(x, i, sigma)
      if(i != length(component_models)){
        stop("gam segments are only supported as the last segment of a multi part model")}
    } else if(matched_iter_model %in% c("int", "homo")){
      iter <- .intChngptForm(x, i, sigma)
    }
    return(iter)
  })
  
  if(sigma){y = "sigma"}

  growthForm <- paste0(y, " ~ ", formulae[[1]]$form, " * ", formulae[[1]]$cp)
  for (i in 2:length(formulae)){
    nextPhase <- paste0("+ (",formulae[[(i-1)]]$cpInt," + ", formulae[[i]]$form ,") * ", formulae[[i]]$cp)
    growthForm <- paste0(growthForm, nextPhase)
  }
  growthForm <- stats::as.formula(growthForm)
  
  if(sigma){growthForm <- brms::nlf(growthForm)}
  
  params <- unique(unlist(lapply(formulae, function(f){f$params})))
  params <- params[-length(params)]
  
  splineSegments <- which(unlist(lapply(formulae, function(fml){"splineVar"%in%names(fml)})))
  
  if(length(splineSegments)>0){
    if(sigma){prefix <- "sub"} else { prefix <- NULL}
    if(useGroup){by <- paste0(", by = ", group)
    } else {by <- "," }
    if(nTimes < 11){
      k = paste0(", k = ", nTimes)
    } else{ k = NULL }
    splineVars <-c()
    for (seg in splineSegments){
      splineVars <- c(splineVars, formulae[[seg]]$splineVar)
    }
    lhs <- paste0(splineVars, collapse = "+")
    rhs <- paste0("s(",x, by, k,")")
    splineForm <- paste0(lhs, "~", rhs)
  } else {splineForm <- NULL}
  
  
  
  return(list("growthForm"= growthForm, "pars" = params, "splineHelperForm" = splineForm))
}

#* ****************************************
#* ***** `Linear Changepoint Phase` *****
#* ****************************************

#' Linear changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#'
#' @examples
#'
#' .linearChngptForm(x="time", 1)
#' .linearChngptForm(x="time", 2)
#' .linearChngptForm(x="time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is 
#' used in starting the next growth phase from the right y value.
#' @noRd

.linearChngptForm <- function(x, position=1, sigma = FALSE){ # return f, cp, and cpInt 
  if(sigma){prefix <- "sub"} else { prefix <- NULL}
  if(position==1){
    form <- paste0(prefix,"linear",position,"A * ", x)
    cp <- paste0("inv_logit((",prefix,"changePoint1 - ",x,") * 5)")
    cpInt <- paste0("(",prefix,"linear", position,"A * ",prefix,"changePoint1)") # intercept at END of this growth phase
  } else{
    form <- paste0(prefix,"linear", position, "A * (", x,"-", paste0(prefix, "changePoint", 1:(position-1), collapse = "-"), ")")
    cp <- paste0("inv_logit((", x,"-", paste0(prefix, "changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
        paste0(prefix, "linear", i, "A * (", paste0(prefix, "changePoint", i:1, collapse="-"),")")
      }), collapse=" + "))
  }
  pars <- c(paste0(prefix, "linear", position, "A"),
            paste0(prefix, "changePoint", position))
  return(list("form" = form,
              "cp" = cp,
              "cpInt" = cpInt,
              "params" = pars))
}


#* ****************************************
#* ***** `Logistic Changepoint Phase` *****
#* ****************************************

#' Logistic changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#'
#' @examples
#'
#' .logisticChngptForm(x="time", 1)
#' .logisticChngptForm(x="time", 2)
#' .logisticChngptForm(x="time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is 
#' used in starting the next growth phase from the right y value.
#' 
#' @noRd

.logisticChngptForm <- function(x, position=1, sigma = FALSE){ # return f, cp, and cpInt 
  if(sigma){prefix <- "sub"} else { prefix <- NULL}
  if(position==1){
    form <- paste0(prefix,"logistic",position,"A / (1 + exp( (",prefix,"logistic",position,"B-(",x,"))/",prefix,"logistic",position,"C) )")
    cp <- paste0("inv_logit((",prefix,"changePoint1 - ",x,") * 5)")
    cpInt <- paste0(prefix,"logistic",position,"A / (1 + exp( (",prefix,"logistic",position,"B-(",prefix,"changePoint1))/",prefix,"logistic",position,"C) )") 
  } else{
    form <- paste0(prefix, "logistic",position,"A / (1 + exp( (",prefix,"logistic",position,
                   "B-(",x,"-",paste0(prefix,"changePoint", 1:(position-1), collapse = "-"),
                   "))/",prefix,"logistic",position,"C) )")
    cp <- paste0("inv_logit((", x,"-", paste0(prefix,"changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0(prefix,"logistic",position,"A / (1 + exp( (",prefix,"logistic",position,"B-(",paste0(prefix,"changePoint", i:1, collapse="-"),"))/",prefix,"logistic",position,"C) )")
    }), collapse=" + "))
  }
  pars <- c(paste0(prefix,"logistic", position, c("A", "B", "C")),
            paste0(prefix,"changePoint", position))
  return(list("form" = form,
              "cp" = cp,
              "cpInt" = cpInt,
              "params" = pars))
}


#* ****************************************
#* ***** `Gompertz Changepoint Phase` *****
#* ****************************************

#' Gompertz changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#'
#' @examples
#'
#' .gompertzChngptForm(x="time", 1)
#' .gompertzChngptForm(x="time", 2)
#' .gompertzChngptForm(x="time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is 
#' used in starting the next growth phase from the right y value.
#' @noRd

.gompertzChngptForm <- function(x, position=1, sigma = FALSE){ # return f, cp, and cpInt 
  if(sigma){prefix <- "sub"} else { prefix <- NULL}
  if(position==1){
    form <- paste0(prefix,"gompertz", position, "A * exp(-",prefix,"gompertz", position, "B * exp(-",prefix,"gompertz", position, "C * ", x, "))" )
    cp <- paste0("inv_logit((",prefix,"changePoint1 - ",x,") * 5)")
    cpInt <- paste0(prefix,"gompertz", position, "A * exp(-",prefix,"gompertz", position, "B * exp(-",prefix,"gompertz", position, "C * ",prefix,"changePoint1))" )
  } else{
    form <- paste0(prefix,"gompertz", position, "A * exp(-",prefix,"gompertz", position,
                   "B * exp(-",prefix,"gompertz", position, "C * (", x," - ",
                   paste0(prefix,"changePoint", 1:(position-1), collapse = "-"), ")))" )
    
    cp <- paste0("inv_logit((", x,"-", paste0(prefix, "changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0(prefix, "gompertz", position, "A * exp(-",prefix,"gompertz", position, "B * exp(-",prefix,"gompertz", position,
             "C * (", paste0(prefix,"changePoint", i:1, collapse="-"), "))" )
    }), collapse=" + "))
  }
  pars <- c(paste0(prefix,"gompertz", position, c("A", "B", "C")),
            paste0(prefix, "changePoint", position))
  return(list("form" = form,
              "cp" = cp,
              "cpInt" = cpInt,
              "params" = pars))
}

#* ****************************************
#* ***** `monomolecular Changepoint Phase` *****
#* ****************************************

#' Monomolecular changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#'
#' @examples
#'
#' .monomolecularChngptForm(x="time", 1)
#' .monomolecularChngptForm(x="time", 2)
#' .monomolecularChngptForm(x="time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is 
#' used in starting the next growth phase from the right y value.
#' 
#' @noRd

.monomolecularChngptForm <- function(x, position=1, sigma = FALSE){ # return f, cp, and cpInt 
  if(sigma){prefix <- "sub"} else { prefix <- NULL}
  if(position==1){
    form <- paste0(prefix, "monomolecular",position,"A-",prefix,"monomolecular",position,"A * exp(-",prefix,"monomolecular",position,"B * ",x,")")
    cp <- paste0("inv_logit((",prefix,"changePoint1 - ",x,") * 5)")
    cpInt <-  paste0(prefix,"monomolecular",position,"A-",prefix,"monomolecular",position,"A * exp(-",prefix,"monomolecular",position,"B * ",prefix,"changePoint1)")
  } else{
    form <- paste0(prefix,"monomolecular",position,"A-",prefix,"monomolecular",position,"A * exp(-",prefix,"monomolecular",position,"B * ",x,"-",
                   paste0(prefix,"changePoint", 1:(position-1), collapse = "-"),")")
    
    cp <- paste0("inv_logit((", x,"-", paste0(prefix, "changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0(prefix, "monomolecular",position,"A-",prefix,"monomolecular",position,"A * exp(-",prefix,"monomolecular",position,"B * ",
             paste0(prefix, "changePoint", i:1, collapse="-"),")")
      }), collapse=" + "))
  }
  pars <- c(paste0(prefix, "monomolecular", position, c("A", "B")),
            paste0(prefix, "changePoint", position))
  return(list("form" = form,
              "cp" = cp,
              "cpInt" = cpInt,
              "params" = pars))
}


#* ****************************************
#* ***** `Exponential Changepoint Phase` *****
#* ****************************************

#' Exponential changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#'
#' @examples
#'
#' .exponentialChngptForm(x="time", 1)
#' .exponentialChngptForm(x="time", 2)
#' .exponentialChngptForm(x="time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is 
#' used in starting the next growth phase from the right y value.
#' @noRd

.exponentialChngptForm <- function(x, position=1, sigma = FALSE){ # return f, cp, and cpInt 
  if(sigma){prefix <- "sub"} else { prefix <- NULL}
  if(position==1){
    form <- paste0(prefix, "exponential",position,"A * exp(",prefix,"exponential", position, "B * ",x,")")
    cp <- paste0("inv_logit((",prefix,"changePoint1 - ",x,") * 5)")
    cpInt <-  paste0(prefix,"exponential",position,"A * exp(",prefix,"exponential", position, "B * ",prefix,"changePoint1)")
  } else{
    form <- paste0(prefix,"exponential",position,"A * exp(",prefix,"exponential", position, "B * (",
                   x, "-", paste0(prefix, "changePoint", 1:(position-1), collapse = "-"),"))")
    cp <- paste0("inv_logit((", x,"-", paste0(prefix, "changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0(prefix,"exponential",position,"A * exp(",prefix,"exponential", position, "B * (",
             paste0(prefix, "changePoint", i:1, collapse="-"),"))")
    }), collapse=" + "))
  }
  pars <- c(paste0(prefix, "exponential", position, c("A", "B")),
            paste0(prefix, "changePoint", position))
  return(list("form" = form,
              "cp" = cp,
              "cpInt" = cpInt,
              "params" = pars))
}

#* ****************************************
#* ***** `Power Law Changepoint Phase` *****
#* ****************************************

#' Power Law changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#'
#' @examples
#'
#' .powerLawChngptForm(x="time", 1)
#' .powerLawChngptForm(x="time", 2)
#' .powerLawChngptForm(x="time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is 
#' used in starting the next growth phase from the right y value.
#' 
#' @noRd

.powerLawChngptForm <- function(x, position=1, sigma = FALSE){ # return f, cp, and cpInt 
  if(sigma){prefix <- "sub"} else { prefix <- NULL}
  if(position==1){
    form <- paste0(prefix, "powerLaw",position,"A * ", x, "^(",prefix,"powerLaw",position,"B)")
    cp <- paste0("inv_logit((",prefix,"changePoint1 - ",x,") * 5)")
    cpInt <-  paste0(prefix, "powerLaw",position,"A * ",prefix,"changePoint1^(",prefix,"powerLaw",position,"B)")
  } else{
    form <- paste0(prefix, "powerLaw",position,"A * ",  x, "-", paste0(prefix, "changePoint", 1:(position-1), collapse = "-"), "^(",prefix,"powerLaw",position,"B)")
    cp <- paste0("inv_logit((", x,"-", paste0(prefix, "changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0(prefix, "powerLaw",position,"A * (",paste0(prefix,"changePoint", i:1, collapse="-"),")^(",prefix,"powerLaw",position,"B)")
    }), collapse=" + "))
  }
  pars <- c(paste0(prefix,"powerLaw", position, c("A", "B")),
            paste0(prefix, "changePoint", position))
  return(list("form" = form,
              "cp" = cp,
              "cpInt" = cpInt,
              "params" = pars))
}

#* ****************************************
#* ***** `Intercept Changepoint Phase` *****
#* ****************************************

#' intercept only changepoint section function
#' 
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#'
#' @examples
#'
#' .intChngptForm(x="time", 1, nTimes = 20)
#' .intChngptForm(x="time", 2, nTimes = 20)
#' .intChngptForm(x="time", 3, nTimes = 5)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is 
#' used in starting the next growth phase from the right y value, for GAMs this is 
#' undefined and GAMs should only be used at the end of a segmented model.
#' @noRd

.intChngptForm <- function(x, position=1, sigma = FALSE){ # return f, cp, and cpInt 
  if(sigma){prefix <- "sub"} else { prefix <- NULL}
  
  if(position==1){
    form <- paste0(prefix,"int",position) 
    cp <- paste0("inv_logit((",prefix,"changePoint1 - ",x,") * 5)")
    cpInt <- paste0(prefix,"int",position)
  } else{
    form <-  paste0(prefix,"int",position) 
    cp <- paste0("inv_logit((", x,"-", paste0(prefix, "changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- paste0(prefix,"int",position)
  }
  
  pars <- c(paste0(prefix,"int",position),
            paste0(prefix, "changePoint", position))
  return(list("form" = form,
              "cp" = cp,
              "cpInt" = cpInt,
              "params" = pars))
}




#* ****************************************
#* ***** `Gam Changepoint Phase` *****
#* ****************************************

#' gam changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#'
#' @examples
#'
#' .gamChngptForm(x="time", 1, nTimes = 20)
#' .gamChngptForm(x="time", 2, nTimes = 20)
#' .gamChngptForm(x="time", 3, nTimes = 5)
#'
#' @return a list with form, cp, cpInt, and splineForm elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is 
#' used in starting the next growth phase from the right y value, for GAMs this is 
#' undefined and GAMs should only be used at the end of a segmented model.
#' "splineForm" is to use in making a spline for a predictor.
#' @noRd

.gamChngptForm <- function(x, position=1, sigma = FALSE){ # return f, cp, and cpInt
  if(sigma){prefix <- "sub"} else { prefix <- NULL}
  
  if(position==1){
    stop("GAMs are only supported as the last function of a multi-part formula")
  } else{
    form <- paste0(prefix, "spline") # paste0("s(",x, by, k,")")
    cp <- paste0("inv_logit((", x,"-", paste0(prefix,"changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- NA
  }
  pars <- c( paste0(prefix, "spline") , paste0(prefix,"changePoint", position))
  return(list("form" = form,
              "cp" = cp,
              "cpInt" = cpInt,
              "params" = pars,
              "splineVar" = paste0(prefix, "spline") ))
}



