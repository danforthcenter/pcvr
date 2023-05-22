#' Ease of use brms starter function for 6 growth model parameterizations
#' 
#' @param model The name of a model as a character string. Supported options are c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law"). See \code{\link{growthSim}} for examples of each type of growth curve.
#' @param form A formula describing the model. The left hand side should only be the outcome variable (phenotype). The right hand side needs at least the x variable (typically time). Grouping is also described in this formula using roughly lme4 style syntax, with formulas like \code{y~time|individual/group} to show that predictors should vary by \code{group} and autocorrelation between \code{individual:group} interactions should be modeled.
#' @param sigma A model for heteroskedasticity from c("homo", "linear", "spline"). 
#' @param df A dataframe to use. Must contain all the variables listed in the formula.
#' @param priors A named list of means for prior distributions. Currently this function makes lognormal priors for all growth model parameters. This is done because the values are strictly positive and the lognormal distribution is easily interpreted. If this argument is not provided then priors are not returned and a different set of priors will need to be made for the model using \code{brms::set_prior}. This works similarly to the \code{params} argument in \code{growthSim}. Names should correspond to parameter names from the \code{model} argument. A numeric vector can also be used, but specifying names is best practice for clarify. See details.
#' @keywords Bayesian, brms
#' @return A named list of elements to make it easier to fit common brms models.
#' @examples 
#' 
#' simdf<-growthSim("logistic", n=20, t=25, params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="spline", df=simdf, priors = list("A"=130, "B"=12, "C"=3))
#' lapply(ss,class)
#' ss$initfun()
#' fit_test <- brm(ss$formula, prior = ss$prior, data = ss$df, family = ss$family, # main components of the model
#'               iter = 1000, cores = 2, chains = 2, init = ss$initfun, # parameters controling chain number, chain length, parallelization and starting values
#'               control = list(adapt_delta = 0.999, max_treedepth = 20), backend = "cmdstanr") # options to increase performance
#' @export

growthSS<-function(model, form, sigma=NULL, df, priors=NULL){
  out<-list()
  sigmas<-c("homo", "linear", "spline")
  models<-c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law")
  #* ***** `Make bayesian formula` *****
    #* `parse form argument`
    y=as.character(form)[2]
    x<-as.character(form)[3]
    if(grepl("\\|", x) | grepl("\\/",x)){
      x3<-trimws(strsplit(x, "[|]|[/]")[[1]])
      x<-x3[1]
      individual = x3[2]
      group = x3[3]
    } else {stop("form must specify grouping for observations. See documentation and examples.")}
    #* `Make autocorrelation formula`
    corForm<-as.formula(paste0("~arma(~",x,"|", individual,":",group,",1,1)"))
    #* `Make parameter formula`
    if(match.arg(model, models) %in% c("logistic", "gompertz")){
      pars=c("A","B", "C")
    } else if(match.arg(model, models) %in% c("monomolecular", "exponential", "power law")){
      pars=c("A","B")
    } else if(match.arg(model, models) == "linear"){
      pars="A"
    }
    parForm<-as.formula(paste0( paste(pars,collapse="+"),"~0+",group ))
    #* `Make growth formula`
    if(match.arg(model, models)=="logistic"){
      form_fun<-form_logistic
    } else if (match.arg(model, models)=="gompertz"){
      form_fun<-form_gompertz
    } else if (match.arg(model, models)=="monomolecular"){
      form_fun<-form_monomolecular
    } else if (match.arg(model, models)=="exponential"){
      form_fun<-form_exponential
    } else if (match.arg(model, models)=="linear"){
      form_fun<-form_linear
    } else if (match.arg(model, models)=="power law"){
      form_fun<-form_powerlaw
    }
    growthForm = form_fun(x,y)
    #* `Make heteroskedasticity formula`
    if(!is.null(sigma)){
      sigmaForm<-if(match.arg(sigma, choices=sigmas)=="homo"){
        as.formula(paste0("sigma ~ 0+", group))
      } else if (match.arg(sigma, choices=sigmas)=="linear"){
        as.formula(paste0("sigma ~ ", x, "+", x, ":",group))
      } else if(match.arg(sigma, choices=sigmas)=="spline"){
        if(length(unique(df[[x]]))<11){
          as.formula(paste0("sigma ~ s(",x,", by=", group, ", k=",length(unique(df[[x]])),")"))
        }else{ as.formula(paste0("sigma ~ s(",x,", by=", group, ")")) }
      } 
      #* `Combining for brms formula`
      bayesForm<-brms::bf(formula = growthForm, sigmaForm, parForm, autocor = corForm, nl=T)
    }else{
      bayesForm<-brms::bf(formula = growthForm,parForm,autocor = corForm,nl=T)
    }
    out[["formula"]]<-bayesForm
  #* ***** `Make make priors` *****
  if(!is.null(priors)){
    # priors = list("A"=130, "B"=12, "C"=3)
    # priors = list("A"=c(130, 150), "B"=c(12,11), "C"=c(3,3))
    #* `if priors is a brmsprior`
    if(any(class(priors)=="brmsprior")){out[["prior"]]<-priors } else{
    #* `if priors is a numeric vector`
    if(is.numeric(priors)){ # priors = c(130, 12, 3) 
      if(length(priors)==length(pars)){
      warning("Assuming that prior is in order: ", paste0(pars, collapse=', '))
        priors<-as.list(priors)
        names(priors)<-pars
      } else {stop(paste0("`priors` is length ", length(priors), " while the specified model requires ", length(pars), " parameters."))}
      priors<-lapply(1:length(unique(df[[group]])), function(i) priors)
      names(priors)<-unique(df[[group]])
    }
    #* `if priors is a list`
    # priors = list("A"=c(130, 150), "B"=c(12,11), "C"=c(3,3))
    # priors = list("A"=c("a"=130, "b"=150), "B"=c(12,11), "C"=c(3,3))
    # priors = list("A"=c("a"=130, "b"=150), "B"=c(12, 11), "C"=c(3))
    # priors = list("A"=c(130, 150), "B"=c(12,11), "C"=c(3))
    if(is.list(priors)){
      if(is.null(names(priors))){
        warning("Assuming that each element in priors is in order: ", paste0(pars, collapse=', '))
        names(priors)<-pars
      }
      if(any(unlist(lapply(priors,length))>1)){ # if more than one value is specified per parameter
        ml<-max(unlist(lapply(priors, length)))
        priors<-lapply(priors, function(p) rep(p, length.out=ml))
        if(any(unlist(lapply(priors, function(p) !is.null(names(p)))))){ # if any inner values are named then apply that to all priors
          wch<-which(unlist(lapply(priors, function(p) !is.null(names(p)))))
          nms<-names(priors[[wch]])
          for(i in 1:length(priors)){names(priors[[i]])<-nms}
        }
        if(any(unlist(lapply(priors, function(p) is.null(names(p)))))){ # if no inner values were named
          for(i in 1:length(priors)){names(priors[[i]])<-unique(df[[group]])}
        }
      } else{ # else is for prior of length 1 for each element, in which case they need to replicated per groups
        l<-length(unique(df[[group]]))
        priors<-lapply(priors, rep, length.out=l)
        nms<-unique(df[[group]])
        for(i in 1:length(priors)){names(priors[[i]])<-nms}
      }
    }
    }
    if(!any(names(out)=="prior")){
      priorStanStrings<-lapply(pars, function(par) paste0("lognormal(log(", priors[[par]],"), 0.25)") )
      priorStanStrings<-unlist(priorStanStrings)
      parNames<-rep(names(priors), each = length(priors[[1]]))
      groupNames<-rep(names(priors[[1]]), length.out = length(priorStanStrings))
      names(priorStanStrings)<-paste(parNames, groupNames, sep="_")
      # priorStanStrings
      #* assemble set_prior statements from these
      prior<-brms::set_prior('student_t(3,0,5)', dpar="sigma")+brms::set_prior('gamma(2,0.1)', class="nu")
      for(nm in names(priorStanStrings)){
        dist = priorStanStrings[[nm]]
        pr = strsplit(nm, "_")[[1]][1]
        gr = paste0(group, strsplit(nm, "_")[[1]][2])
        prior<-prior+brms::set_prior(dist, coef = gr, nlpar = pr)
      }
      out[["prior"]]<-prior
    }
  }
    
  #* ***** `Make initializer function` *****
    initFun<-function(pars="?", nPerChain=1){
      init<-lapply(pars, function(i) rgamma(nPerChain,1))
      names(init)<-paste0("b_",pars)
      init
    }
    formals(initFun)$pars<-pars
    formals(initFun)$nPerChain<-length(unique(df[[group]]))
    wrapper<-function(){initFun()}
    out[["initfun"]]<-wrapper
    out[["df"]]<-df
    out[["family"]]<-"student"
  return(out)
  }
  
form_logistic<-function(x, y){
  return(as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C))")))
}
form_gompertz<-function(x, y){
  return(as.formula(paste0(y," ~ A*exp(-B*exp(-C*",x,"))")))
}
form_monomolecular<-function(x, y){
  return(as.formula(paste0(y,"~A-A*exp(-B*",x,")")))
}
form_exponential<-function(x, y){
  return(as.formula(paste0(y," ~ A*exp(B*",x, ")")))
}
form_linear<-function(x, y){
  return(as.formula(paste0(y," ~ A*",x)))
}
form_powerlaw<-function(x, y){
  return(as.formula(paste0(y," ~ A*",x,"^B")))
}










