#' Ease of use nlrq starter function for 6 growth model parameterizations
#' 
#' Internal to growthSS
#' 
#' @param model One of the 8 model options in growthSS
#' @param form A pcvr style form, see growthSS
#' @param sigma One of "none", "power", or "exp", which will correspond to varIdent, varPower, or varExp respectively.
#' This is also meant to take a varFunc object, but that is untested.
#' @param df a dataframe to use to make the model.
#' @param pars optional parameters to vary by group as fixed effects.
#' @param start Starting values. These are optional unless model is a double sigmoid.
#' For any other model these will be estimated from the data if left NULL.
#' 
#' @examples 
#' 
#' if(FALSE){
#' ex<-growthSim("logistic", n=20, t=25,
#'   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ex$group<-factor(ex$group) # group MUST be a factor for nlme (typical)
#' nl1<-nlme::nlme(model = y ~ A/(1 + exp((B - time)/C)),
#'           data = ex,
#'           fixed = list(
#'             A ~ 0 + group,
#'             B ~ 0 + group,
#'             C ~ 0 + group),
#'           random = A +B +C~ 1,
#'           groups = ~group,
#'           weights = nlme::varIdent(form = ~1|group),
#'           start = c(150, 150,10, 10,3, 3) )
#' # linear sigma using power
#' nl2<-nlme::nlme(model = y ~ A/(1 + exp((B - time)/C)),
#'           data = ex,
#'           fixed = list(
#'            A ~ 0 + group,
#'             B ~ 0 + group,
#'             C ~ 0 + group),
#'           random = A ~ 1,
#'           groups = ~group,
#'           weights = nlme::varPower(form = ~time|group),
#'           start = c(150, 150,10, 10,3, 3) )
#' # linear sigma using exp
#' nl3<-nlme::nlme(model = y ~ A/(1 + exp((B - time)/C)),
#'           data = ex,
#'           fixed = list(
#'             A ~ 0 + group,
#'             B ~ 0 + group,
#'             C ~ 0 + group),
#'           random = A ~ 1,
#'           groups = ~group,
#'           weights = nlme::varExp(form = ~time|group),
#'           start = c(150, 150,10, 10,3, 3) )
#' }
#' simdf<-growthSim("logistic", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' 
#' ss<-.nlmeSS(model = "logistic", form=y~time|id/group,
#'   sigma="power", df=simdf, start=NULL)
#' dim(ss$df)
#' ss[c("formula", "taus", "start", "pcvrForm")]
#' 
#' @importFrom nlme varPower varIdent varExp nlme nlme.formula corAR1
#' @importFrom stats as.formula
#' @importFrom methods is
#' 
#' @keywords internal
#' @noRd


.nlmeSS <- function(model, form, sigma, df, pars = NULL, start=NULL){
#* `general steps`
#* ***** `Define choices and make empty output list`
out<-list()
models<-c("logistic", "gompertz", "monomolecular",
          "exponential", "linear", "power law",
          "double logistic", "double gompertz")
sigmas = c("none", "power", "exp")
#* check if sigma is class "varFunc", if it is then return it as is?

#* ***** `Make nlme model formula` *****
#* `parse form argument`
y=as.character(form)[2]
x<-as.character(form)[3]
group = NULL
if(grepl("\\|", x) & grepl("\\/",x)){
  x3<-trimws(strsplit(x, "[|]|[/]")[[1]])
  x<-x3[1]
  individual = x3[2]
  group = x3[3] 
  if(length(unique(df[[group]]))==1){USEGROUP=FALSE}else{USEGROUP=TRUE} # if there is only one group then ignore grouping for parameter and variance formulas
} else if (grepl("\\|", x)){
  x2<-trimws(strsplit(x, "[|]")[[1]])
  x<-x2[1]
  individual = x2[2]
  group="group"
  df[[group]]<-"group"
  USEGROUP=FALSE
} else {stop("form must specify grouping. See documentation and examples.")}

#* `make group a factor for nlme`
df[[group]]<-as.factor(df[[group]])
#* `make an interaction variable for autocorrelation`
df[[paste0(group, individual)]]<-interaction(df[[group]], df[[individual]])
intVar <- paste0(group, individual)

#* `assemble growth formula with FE, RE, Groups, and Weights`
matched_model <- match.arg(model, models)

if(is.character(sigma)){
  matched_sigma <- match.arg(sigma, sigmas)
  } else{ matched_sigma = sigma}

if(matched_model=="double logistic"){
  if(is.null(pars)){ pars = c("A", "B", "C", "A2", "B2", "C2") }
  form_fun<-.nlme_form_dou_logistic
} else if (matched_model=="double gompertz"){
  if(is.null(pars)){ pars = c("A", "B", "C", "A2", "B2", "C2") }
  form_fun<-.nlme_form_dou_gompertz
} else if(matched_model=="logistic"){
  if(is.null(pars)){ pars = c("A", "B", "C") }
  form_fun<-.nlme_form_logistic
} else if (matched_model=="gompertz"){
  if(is.null(pars)){ pars = c("A", "B", "C") }
  form_fun<-.nlme_form_gompertz
} else if (matched_model=="monomolecular"){
  if(is.null(pars)){ pars = c("A", "B") }
  form_fun<-.nlme_form_monomolecular
} else if (matched_model=="exponential"){
  if(is.null(pars)){ pars = c("A", "B") }
  form_fun<-.nlme_form_exponential
} else if (matched_model=="linear"){
  if(is.null(pars)){ pars = c("A") }
  form_fun<-.nlme_form_linear
} else if (matched_model=="power law"){
  if(is.null(pars)){ pars = c("A", "B") }
  form_fun<-.nlme_form_powerlaw
}
growthForm_list = form_fun(x,y, group, individual, intVar, matched_sigma, pars)

#* `Make starting values`
if(is.null(start)){
  if(matched_model=="double logistic"){
    warning("Double Sigmoid models are not supported as self-starting models, you will need to add starting parameters. Note for these models type='brms' is recommended.")
    startingValues<-NULL
  } else if (matched_model=="double gompertz"){
    warning("Double Sigmoid models are not supported as self-starting models, you will need to add starting parameters. Note for these models type='brms' is recommended.")
    startingValues<-NULL
  } else if(matched_model=="logistic"){
    startingValues <- .initLogistic(df, x, y) # see nlrqSS.R helper functions
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
  startingValuesList <- unlist(lapply(names(startingValues), function(nm){
    val<-startingValues[nm]
    if(nm %in% pars){
      rep(val, length(unique(df[[group]])) ) # if this is one of pars then make starting value per group
    } else{ val }# else return one starting value
  }))
} else{
  startingValuesList <- start
}

#* `return model components`

out[["formula"]] <- growthForm_list
out[["start"]] <- startingValuesList
out[["df"]]<-df
out[["pcvrForm"]]<-form
return(out)
}







#* `Sigma matching helper function`

.nlme_sigma_form <- function(matched_sigma, x, group){
  #* `variance formula`
  if(methods::is(matched_sigma, "varFun")){
    weights_form = matched_sigma
  } else if(matched_sigma =="none"){
    weights_form = nlme::varIdent(form = stats::as.formula(paste0("~1|", group)))
  } else if(matched_sigma=="power"){
    weights_form = nlme::varPower(form = stats::as.formula(paste0("~",x,"|",group)))
  } else if(matched_sigma =="exp"){
    weights_form = nlme::varExp(form = stats::as.formula(paste0("~",x,"|",group)))
  }
  return(weights_form)
}


#* `Define growth formulas`

.nlme_form_logistic<-function(x, y, group, individual, intVar, matched_sigma, pars){
  #* `main growth formula`
  model_form <- as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C))"))
  #* `random effects formula`
  random_form <- A + B + C ~ 1
  #* `fixed effects formula`
  total_pars <- c("A", "B", "C")
  fixed_form <- lapply(total_pars, function(par){
    if(par %in% pars){
      stats::as.formula(paste0(par," ~ 0 + ", group))
    } else{
      stats::as.formula(paste0(par," ~ 1"))
    }
  })
  #* `groups formula`
  groups_form = stats::as.formula(paste0("~", intVar))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |",intVar)))
  
  formulas <- list("model"=model_form, "random"=random_form,
                   "fixed"=fixed_form, "groups"=groups_form,
                   "weights" = weights_form, "cor_form" = correlation_form)
  return(formulas)
}

.nlme_form_gompertz<-function(x, y, group, individual, intVar, matched_sigma, pars){
  model_form <- as.formula(paste0(y," ~ A*exp(-B*exp(-C*",x,"))"))
  #* `random effects formula`
  random_form <- A + B + C ~ 1
  #* `fixed effects formula`
  total_pars <- c("A", "B", "C")
  fixed_form <- lapply(total_pars, function(par){
    if(par %in% pars){
      stats::as.formula(paste0(par," ~ 0 + ", group))
    } else{
      stats::as.formula(paste0(par," ~ 1"))
    }
  })
  #* `groups formula`
  groups_form = stats::as.formula(paste0("~", intVar))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |",intVar)))
  
  formulas <- list("model"=model_form, "random"=random_form,
                   "fixed"=fixed_form, "groups"=groups_form,
                   "weights" = weights_form, "cor_form" = correlation_form)
  return(formulas)
}

.nlme_form_dou_logistic<-function(x, y, group, individual, intVar, matched_sigma, pars){
  model_form <- as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C)) + ((A2-A) /(1+exp((B2-",x,")/C2)))"))
  #* `random effects formula`
  random_form <- A + B + C + A2 + B2 + C2 ~ 1
  #* `fixed effects formula`
  total_pars <- c("A", "B", "C", "A2", "B2", "C2")
  fixed_form <- lapply(total_pars, function(par){
    if(par %in% pars){
      stats::as.formula(paste0(par," ~ 0 + ", group))
    } else{
      stats::as.formula(paste0(par," ~ 1"))
    }
  })
  #* `groups formula`
  groups_form = stats::as.formula(paste0("~", intVar))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |",intVar)))
  
  formulas <- list("model"=model_form, "random"=random_form,
                   "fixed"=fixed_form, "groups"=groups_form,
                   "weights" = weights_form, "cor_form" = correlation_form)
  return(formulas)
}

.nlme_form_dou_gompertz<-function(x, y, group, individual, intVar, matched_sigma, pars){
    model_form <- as.formula(paste0(y," ~ A * exp(-B * exp(-C*",x,")) + (A2-A) * exp(-B2 * exp(-C2*(",x,"-B)))"))
    #* `random effects formula`
    random_form <- A + B + C + A2 + B2 + C2 ~ 1
    #* `fixed effects formula`
    total_pars <- c("A", "B", "C", "A2", "B2", "C2")
    fixed_form <- lapply(total_pars, function(par){
      if(par %in% pars){
        stats::as.formula(paste0(par," ~ 0 + ", group))
      } else{
        stats::as.formula(paste0(par," ~ 1"))
      }
    })
    #* `groups formula`
    groups_form = stats::as.formula(paste0("~", intVar))
    #* `variance formula`
    weights_form <- .nlme_sigma_form(matched_sigma, x, group)
    #* `correlation formula`
    correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |",intVar)))
    
    formulas <- list("model"=model_form, "random"=random_form,
                     "fixed"=fixed_form, "groups"=groups_form,
                     "weights" = weights_form, "cor_form" = correlation_form)
    return(formulas)
}

.nlme_form_monomolecular<-function(x, y, group, individual, intVar, matched_sigma, pars){
  model_form <- as.formula(paste0(y,"~A-A*exp(-B*",x,")"))
  #* `random effects formula`
  random_form <- A + B ~ 1
  #* `fixed effects formula`
  total_pars <- c("A", "B")
  fixed_form <- lapply(total_pars, function(par){
    if(par %in% pars){
      stats::as.formula(paste0(par," ~ 0 + ", group))
    } else{
      stats::as.formula(paste0(par," ~ 1"))
    }
  })
  #* `groups formula`
  groups_form = stats::as.formula(paste0("~", intVar))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |",intVar)))
  
  formulas <- list("model"=model_form, "random"=random_form,
                   "fixed"=fixed_form, "groups"=groups_form,
                   "weights" = weights_form, "cor_form" = correlation_form)
  return(formulas)
}

.nlme_form_exponential<-function(x, y, group, individual, intVar, matched_sigma, pars){
  model_form <- as.formula(paste0(y," ~ A*exp(B*",x, ")"))
  #* `random effects formula`
  random_form <- A + B ~ 1
  #* `fixed effects formula`
  total_pars <- c("A", "B")
  fixed_form <- lapply(total_pars, function(par){
    if(par %in% pars){
      stats::as.formula(paste0(par," ~ 0 + ", group))
    } else{
      stats::as.formula(paste0(par," ~ 1"))
    }
  })
  #* `groups formula`
  groups_form = stats::as.formula(paste0("~", intVar))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |",intVar)))
  
  formulas <- list("model"=model_form, "random"=random_form,
                   "fixed"=fixed_form, "groups"=groups_form,
                   "weights" = weights_form, "cor_form" = correlation_form)
  return(formulas)
}

.nlme_form_linear<-function(x, y, group, individual, intVar, matched_sigma, pars){
  model_form <- as.formula(paste0(y," ~ A*",x))
  #* `random effects formula`
  random_form <- A ~ 1
  #* `fixed effects formula`
  total_pars <- c("A")
  fixed_form <- lapply(total_pars, function(par){
    if(par %in% pars){
      stats::as.formula(paste0(par," ~ 0 + ", group))
    } else{
      stats::as.formula(paste0(par," ~ 1"))
    }
  })
  #* `groups formula`
  groups_form = stats::as.formula(paste0("~", intVar))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |",intVar)))
  
  formulas <- list("model"=model_form, "random"=random_form,
                   "fixed"=fixed_form, "groups"=groups_form,
                   "weights" = weights_form, "cor_form" = correlation_form)
  return(formulas)
}

.nlme_form_powerlaw<-function(x, y, group, individual, intVar, matched_sigma, pars){
  model_form <- as.formula(paste0(y," ~ A*",x,"^B"))
  #* `random effects formula`
  random_form <- A + B ~ 1
  #* `fixed effects formula`
  total_pars <- c("A", "B")
  fixed_form <- lapply(total_pars, function(par){
    if(par %in% pars){
      stats::as.formula(paste0(par," ~ 0 + ", group))
    } else{
      stats::as.formula(paste0(par," ~ 1"))
    }
  })
  #* `groups formula`
  groups_form = stats::as.formula(paste0("~", intVar))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |",intVar)))
  
  formulas <- list("model"=model_form, "random"=random_form,
                   "fixed"=fixed_form, "groups"=groups_form,
                   "weights" = weights_form, "cor_form" = correlation_form)
  return(formulas)
}

