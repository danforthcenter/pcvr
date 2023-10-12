#' Helper function to put together formulae for brms changepoint growth models
#' 
#' @param model A multi-part model passed from brmSS passed from \code{\link{growthSS}}
#' @param x The x variable from the pcvrForm argument in \code{\link{growthSS}}
#' @param y The y variable from the pcvrForm argument in \code{\link{growthSS}}
#' @param group The grouping variable from the pcvrForm argument in \code{\link{growthSS}}
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

.brmsChangePointHelper <- function(model, x, y, group){

  component_models <- trimws(strsplit(model, "\\+")[[1]])
  models<-c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law", "gam")
  
  formulae <- lapply(1:length(component_models), function(i){
    iter_model <- component_models[i]
    matched_iter_model <- match.arg(iter_model, models)
    if(matched_iter_model=="logistic"){
      iter <- .logisticChngptForm(x, i)
    } else if (matched_iter_model=="gompertz"){
      iter <- .gompertzChngptForm(x, i)
    } else if (matched_iter_model=="monomolecular"){
      iter <- .monomolecularChngptForm(x, i)
    } else if (matched_iter_model=="exponential"){
      iter <- .exponentialChngptForm(x, i)
    } else if (matched_iter_model=="linear"){
      iter <- .linearChngptForm(x, i)
    } else if (matched_iter_model=="power law"){
      iter <- .powerLawChngptForm(x, i)
    } else if (matched_iter_model=="gam"){
      iter <- .gamChngptForm(x, group, i)
      if(i != length(component_models)){
        stop("gam segments are only supported as the last segment of a multi part model")}
    }
    return(iter)
  })
  
  growthForm <- paste0(y, " ~ ", formulae[[1]]$form, " * ", formulae[[1]]$cp)
  for (i in 2:length(formulae)){
    nextPhase <- paste0("+ (",formulae[[(i-1)]]$cpInt," + ", formulae[[i]]$form ,") * ", formulae[[i]]$cp)
    growthForm <- paste0(growthForm, nextPhase)
  }
  growthForm <- stats::as.formula(growthForm)
  
  params <- unique(unlist(lapply(formulae, function(f){f$params})))
  params <- params[-length(params)]
  
  return(list("growthForm"= growthForm, "pars" = params))
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

.linearChngptForm <- function(x, position=1){ # return f, cp, and cpInt 
  if(position==1){
    form <- paste0("linear",position,"A * ", x)
    cp <- paste0("inv_logit((changePoint1 - ",x,") * 5)")
    cpInt <- paste0("(linear", position,"A * changePoint1)") # intercept at END of this growth phase
  } else{
    form <- paste0("linear", position, "A * (", x,"-", paste0("changePoint", 1:(position-1), collapse = "-"), ")")
    cp <- paste0("inv_logit((", x,"-", paste0("changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
        paste0("linear", i, "A * (", paste0("changePoint", i:1, collapse="-"),")")
      }), collapse=" + "))
  }
  pars <- c(paste0("linear", position, "A"),
            paste0("changePoint", position))
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

.logisticChngptForm <- function(x, position=1){ # return f, cp, and cpInt 
  if(position==1){
    form <- paste0("logistic",position,"A / (1 + exp( (logistic",position,"B-(",x,"))/logistic",position,"C) )")
    cp <- paste0("inv_logit((changePoint1 - ",x,") * 5)")
    cpInt <- paste0("logistic",position,"A / (1 + exp( (logistic",position,"B-(changePoint1))/logistic",position,"C) )") 
  } else{
    form <- paste0("logistic",position,"A / (1 + exp( (logistic",position,
                   "B-(",x,"-",paste0("changePoint", 1:(position-1), collapse = "-"),
                   "))/logistic",position,"C) )")
    cp <- paste0("inv_logit((", x,"-", paste0("changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0("logistic",position,"A / (1 + exp( (logistic",position,"B-(",paste0("changePoint", i:1, collapse="-"),"))/logistic",position,"C) )")
    }), collapse=" + "))
  }
  pars <- c(paste0("logistic", position, c("A", "B", "C")),
            paste0("changePoint", position))
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

.gompertzChngptForm <- function(x, position=1){ # return f, cp, and cpInt 
  if(position==1){
    form <- paste0("gompertz", position, "A * exp(-gompertz", position, "B * exp(-gompertz", position, "C * ", x, "))" )
    cp <- paste0("inv_logit((changePoint1 - ",x,") * 5)")
    cpInt <- paste0("gompertz", position, "A * exp(-gompertz", position, "B * exp(-gompertz", position, "C * changePoint1))" )
  } else{
    form <- paste0("gompertz", position, "A * exp(-gompertz", position,
                   "B * exp(-gompertz", position, "C * (", x," - ",
                   paste0("changePoint", 1:(position-1), collapse = "-"), ")))" )
    
    cp <- paste0("inv_logit((", x,"-", paste0("changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0("gompertz", position, "A * exp(-gompertz", position, "B * exp(-gompertz", position,
             "C * (", paste0("changePoint", i:1, collapse="-"), "))" )
    }), collapse=" + "))
  }
  pars <- c(paste0("gompertz", position, c("A", "B", "C")),
            paste0("changePoint", position))
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

.monomolecularChngptForm <- function(x, position=1){ # return f, cp, and cpInt 
  if(position==1){
    form <- paste0("monomolecular",position,"A-monomolecular",position,"A * exp(-monomolecular",position,"B * ",x,")")
    cp <- paste0("inv_logit((changePoint1 - ",x,") * 5)")
    cpInt <-  paste0("monomolecular",position,"A-monomolecular",position,"A * exp(-monomolecular",position,"B * changePoint1)")
  } else{
    form <- paste0("monomolecular",position,"A-monomolecular",position,"A * exp(-monomolecular",position,"B * ",x,"-",
                   paste0("changePoint", 1:(position-1), collapse = "-"),")")
    
    cp <- paste0("inv_logit((", x,"-", paste0("changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0("monomolecular",position,"A-monomolecular",position,"A * exp(-monomolecular",position,"B * ",
             paste0("changePoint", i:1, collapse="-"),")")
      }), collapse=" + "))
  }
  pars <- c(paste0("monomolecular", position, c("A", "B")),
            paste0("changePoint", position))
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

.exponentialChngptForm <- function(x, position=1){ # return f, cp, and cpInt 
  if(position==1){
    form <- paste0("exponential",position,"A * exp(exponential", position, "B * ",x,")")
    cp <- paste0("inv_logit((changePoint1 - ",x,") * 5)")
    cpInt <-  paste0("exponential",position,"A * exp(exponential", position, "B * changePoint1)")
  } else{
    form <- paste0("exponential",position,"A * exp(exponential", position, "B * (",
                   x, "-", paste0("changePoint", 1:(position-1), collapse = "-"),"))")
    cp <- paste0("inv_logit((", x,"-", paste0("changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0("exponential",position,"A * exp(exponential", position, "B * (",
             paste0("changePoint", i:1, collapse="-"),"))")
    }), collapse=" + "))
  }
  pars <- c(paste0("exponential", position, c("A", "B")),
            paste0("changePoint", position))
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

.powerLawChngptForm <- function(x, position=1){ # return f, cp, and cpInt 
  if(position==1){
    form <- paste0("powerLaw",position,"A * ", x, "^(powerLaw",position,"B)")
    cp <- paste0("inv_logit((changePoint1 - ",x,") * 5)")
    cpInt <-  paste0("powerLaw",position,"A * changePoint1^(powerLaw",position,"B)")
  } else{
    form <- paste0("powerLaw",position,"A * ",  x, "-", paste0("changePoint", 1:(position-1), collapse = "-"), "^(powerLaw",position,"B)")
    cp <- paste0("inv_logit((", x,"-", paste0("changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- do.call(paste, list(lapply(1:(position), function(i){
      paste0("powerLaw",position,"A * (",paste0("changePoint", i:1, collapse="-"),")^(powerLaw",position,"B)")
    }), collapse=" + "))
  }
  pars <- c(paste0("powerLaw", position, c("A", "B")),
            paste0("changePoint", position))
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
#' .gamChngptForm(x="time", 1)
#' .gamChngptForm(x="time", 2)
#' .gamChngptForm(x="time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is 
#' used in starting the next growth phase from the right y value, for GAMs this is 
#' undefined and GAMs should only be used at the end of a segmented model.
#' @noRd

.gamChngptForm <- function(x, group, position=1){ # return f, cp, and cpInt
  if(position==1){
    stop("GAMs are not supported as the first function of a multi-part formula")
  } else{
    form <- paste0("s(",  x, ", by = ",group,")")
    cp <- paste0("inv_logit((", x,"-", paste0("changePoint", 1:(position-1), collapse = "-"), ") * 5)")
    cpInt <- NA
  }
  pars <- c(paste0("changePoint", position))
  return(list("form" = form,
              "cp" = cp,
              "cpInt" = cpInt,
              "params" = pars))
}





