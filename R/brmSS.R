#' Ease of use brms starter functions
#' 
#' maybe I should make a data simulating function for each supported distribution as well? 
#' something that takes the parameters, the number of samples per group, and returns a dataframe?
#' then the self starter functions should make the bf() call for each supported distribution
#' the self starter function could also return a list of inits for the model. I think that makes sense
#' as a component of the self-starter function because it is a self starting component.
#' This could optionally make a list of priors for you taking similar arguments as the data simulating function
#' and just centering a lognormal distribution on those means. The default would be not to do that though. I guess.
#' It could also take an argument for how many cores/chains to use? Not sure about that yet.
#' 
#' Then the actual brm call would be something like
#' brm(ss$bf, family = student, prior = prior, data = df, iter = 2000, 
#' cores = 4, chains = 4, backend = "cmdstanr",  
#' control = list(adapt_delta = 0.999, max_treedepth = 20),
#' init =ss$init)
#' 
#' @param df Data frame to use. Can be wide or long format from read.pcv
#' @keywords Bayesian, brms
#' @import brms
#' @import ggplot2
#' @examples 
#' 
#' simdf<-growthSim("logistic", n=20, t=25, params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "logistic", form=y~time|id/group, sigma="spline", df=simdf, priors = list("A"=130, "B"=12, "C"=3))
#' lapply(ss,class)
#' ss$initfun()
#' fit_test <- brm(ss$formula, prior = ss$prior, data = ss$df, family = ss$family, # main componenets of the model
#'               iter = 1000, cores = 4, chains = 4, init = ss$initfun, # parameters controling chain number, chain length, parallelization and starting values
#'               control = list(adapt_delta = 0.999, max_treedepth = 20), backend = "cmdstanr") # options to increase performance
#' @export
#' 
#' library(brms)
#' 
#' 
growthSS<-function(model, form, sigma, df, priors=NULL, group=NULL, individual=NULL){
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
    } else if (is.null(individual) | is.null(group)){stop("Either a formula with groupings (outcome ~ predictor|uniqueID/group) or a formula (y~x) and group='variablename1' and individual='variablename2' must be provided.")}
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
    #* `Make heteroskedasticity formula`
    sigmaForm<-if(match.arg(sigma, choices=sigmas)=="homo"){
      as.formula(paste0("sigma ~ 0+", group))
    } else if (match.arg(sigma, choices=sigmas)=="linear"){
      as.formula(paste0("sigma ~ ", x, "+", x, ":",group))
    } else if(match.arg(sigma, choices=sigmas)=="spline"){
      as.formula(paste0("sigma ~ s(",x,", by=", group, ")"))
    }
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
    #* `Combining for brms formula`
    bayesForm<-bf(formula = growthForm,
                  sigmaForm,
                  parForm,
                  autocor = corForm,
                  nl=T)
    out[["formula"]]<-bayesForm
  #* ***** `Make make priors` *****
  if(!is.null(priors)){
    # priors = list("A"=130, "B"=12, "C"=3)
    # priors = list("A"=c(130, 150), "B"=c(12,11), "C"=c(3,3))
    #* `if priors is a brmsprior`
    if(any(class(prior)=="brmsprior")){out[["prior"]]<-prior } else{
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
      prior<-set_prior('student_t(3,0,5)', dpar="sigma")+set_prior('gamma(2,0.1)', class="nu")
      for(nm in names(priorStanStrings)){
        dist = priorStanStrings[[nm]]
        pr = strsplit(nm, "_")[[1]][1]
        gr = paste0(group, strsplit(nm, "_")[[1]][2])
        prior<-prior+set_prior(dist, coef = gr, nlpar = pr)
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
  return(as.formula(paste0(y,"A*exp(-B*exp(-C*",x,"))")))
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










