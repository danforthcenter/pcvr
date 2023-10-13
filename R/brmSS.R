#' Ease of use brms starter function for 6 growth model parameterizations
#' 
#' @param model The name of a model as a character string.
#' Supported options are c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law", "double logistic", "double gompertz", "gam").
#' See \code{\link{growthSim}} for examples of each type of growth curve.
#' @param form A formula describing the model. The left hand side should only be 
#' the outcome variable (phenotype). The right hand side needs at least the x variable
#'  (typically time). Grouping is also described in this formula using roughly lme4
#'  style syntax,with formulas like \code{y~time|individual/group} to show that predictors
#'  should vary by \code{group} and autocorrelation between \code{individual:group}
#'  interactions should be modeled. If group has only one level or is not included then
#'  it will be ignored in formulas for growth and variance (this may be the case if 
#'  you split data before fitting models to be able to run more smaller models each more quickly).
#' @param sigma A model for heteroskedasticity from c("homo", "linear", "spline", "logistic", "gompertz"). See details. 
#' @param df A dataframe to use. Must contain all the variables listed in the formula.
#' @param priors A named list of means for prior distributions.
#'  Currently this function makes lognormal priors for all growth model parameters.
#'   This is done because the values are strictly positive and the lognormal distribution
#'   is easily interpreted. If this argument is not provided then priors are not 
#'   returned and a different set of priors will need to be made for the model using
#'   \code{brms::set_prior}. This works similarly to the \code{params} argument
#'   in \code{growthSim}. Names should correspond to parameter names from the
#'   \code{model} argument. A numeric vector can also be used, but specifying
#'   names is best practice for clarity. Additionally, due to a limitation in
#'   \code{brms} currently lower bounds cannot be set for priors for specific groups.
#'   If priors include multiple groups (\code{priors = list(A = c(10,15), ...)}) then
#'   you will see warnings after the model is fit about not having specified a lower
#'   bound explicitly. Those warnings can safely be ignored and will be addressed if
#'   the necessary features are added to \code{brms}. For GAMs priors are not created by
#'   this function but can still be provided as a \code{brmsprior} object.
#'   See details for guidance.
#' @keywords Bayesian, brms
#' 
#' @importFrom stats as.formula rgamma
#' 
#' @details 
#' 
#' Default priors are not provided, but these can serve as starting points for each distribution. 
#' You are encouraged to use \code{growthSim} to consider what kind 
#' of trendlines result from changes to your prior and for interpretation of each parameter.
#' You should not looking back and forth at your data trying to match your
#'  observed growth exactly with a prior distribution,
#' rather this should be informed by an understanding of the plants you
#'  are using and expectations based on previous research. 
#'  For the "double" models the parameter interpretation is the same
#'  as for their non-double counterparts except that there are A and A2, etc.
#'  It is strongly recommended to familiarize yourself with the double sigmoid 
#'  distributions using growthSim before attempting to model one. Additionally,
#'  those distributions are intended for use with long delays in an experiment,
#'  think stress recovery experiments, not for minor hiccups in plant growth.
#' 
#' \itemize{
#'    \item \bold{Logistic}: \code{list('A' = 130, 'B' = 12, 'C' = 3)}
#'     \item \bold{Gompertz}: \code{list('A' = 130, 'B' = 12, 'C' = 1.25)}
#'     \item \bold{Double Logistic}: \code{list('A' = 130, 'B' = 12, 'C' = 3, 'A2' = 200, 'B2' = 25, 'C2' = 1)}
#'     \item \bold{Double Gompertz}: \code{list('A' = 130, 'B' = 12, 'C' = 0.25, 'A2' = 220, 'B2' = 30, 'C2' = 0.1)}
#'     \item \bold{Monomolecular}: \code{list('A' = 130, 'B' = 2)}
#'     \item \bold{Exponential}: \code{list('A' = 15, 'B' = 0.1)}
#'     \item \bold{Linear}: \code{list('A' = 1)}
#'     \item \bold{Power Law}: \code{list('A' = 13, 'B' = 2)}
#' }
#' 
#' 
#' 
#' The \code{sigma} argument optionally specifies a sub model to account for heteroskedasticity.
#' Currently there are four supported sub models described below.
#' 
#' \itemize{
#'    \item \bold{homo}: \code{sigma ~ 1}, fitting only a global or per group intercept to sigma.
#'     \item \bold{linear}: \code{sigma ~ time}, modeling sigma with a linear relationship to time
#'      and possibly with an interaction term per groups.
#'     \item \bold{spline}: \code{sigma ~ s(time)}, modeling sigma using a smoothing function through `mgcv::s`, possibly by group.
#'     \item \bold{gompertz}: \code{sigma ~ subA * exp(-subB * exp(-subC * x))},
#'      modeling sigma as a gompertz function of time, possibly by group. Note that you
#'      should specify priors for the parameters in this sub model by adding them into the \code{priors}
#'      argument, such as \code{list(..., subA = 20, subB = 15, subC = 0.25)}. If you do not specify priors
#'      then default (flat) priors will be used, which is liable to cause fitting problems and less
#'      accurate results. Looking at your data and making a semi-informed estimate of the total variance at the end
#'      of the experiment can help set a reasonable prior for subA, while subB and subC can generally be the 
#'      same as B and C in a gompertz growth model of the same data. These priors will have fat tails so they
#'      are pretty forgiving.
#' }
#' 
#' 
#' @return A named list of elements to make it easier to fit common brms models.
#' \code{formula}: A \code{brms::bf} formula specifying the growth model, autocorrelation, variance submodel,
#' and models for each variable in the growth model.
#' \code{prior}: A brmsprior/data.frame object.
#' \code{initfun}: A function to randomly initialize chains using a random draw from a gamma
#' distribution (confines initial values to positive and makes correct number
#' of initial values for chains and groups). For "gam" models this initializes all chains at 0.
#' \code{df} The data input, possibly with dummy variables added if needed.
#' \code{family} The model family, currently this will always be "student".
#' \code{pcvrForm} The form argument unchanged. This is returned so that
#' it can be used later on in model visualization. Often it may be a good idea
#' to save the output of this function with the fit model, so having this can
#' be useful later on.
#' @examples 
#' 
#' simdf<-growthSim("logistic", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-.brmSS(model = "logistic", form=y~time|id/group,
#' sigma="spline", df=simdf, priors = list("A"=130, "B"=12, "C"=3))
#' lapply(ss,class)
#' ss$initfun()
#' 
#' 
#' @keywords internal
#' @noRd

.brmSS<-function(model, form, sigma=NULL, df, priors=NULL){
  out<-list()
  sigmas<-c("homo", "linear", "spline", "gompertz")
  models<-c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law", "double logistic", "double gompertz", "gam")
  #* ***** `Make bayesian formula` *****
    #* `parse form argument`
    y=as.character(form)[2]
    x<-as.character(form)[3]
    
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
    } else {stop("form must specify grouping at least for individual level for observations ( y ~ x|individual ). See documentation and examples.")}
    #* `convert group to character to avoid unexpected factor stuff`
    df[[group]]<-as.character(df[[group]])
    #* `Make autocorrelation formula`
    corForm<-as.formula(paste0("~arma(~",x,"|", individual,":",group,",1,1)"))
    #* `match args`
    if(!grepl("\\+", model)){
      matched_model <- match.arg(model, models)
    }
    matched_sigma <- match.arg(sigma, choices=sigmas)
    #* `Make growth formula`
    spline_pars = NULL
    if(grepl("\\+", model)){
      
      chngptHelper_list <- .brmsChangePointHelper(model, x, y, group)
      growthForm <- chngptHelper_list$growthForm
      pars <- chngptHelper_list$pars
      if("splineDummy" %in% pars){
        spline_pars <- "splineDummy"
        pars <- pars[-which(pars=="splineDummy")]
      }
      
    } else{
        if(matched_model=="double logistic"){
          form_fun<-.brms_form_dou_logistic
          pars=c("A","B", "C", "A2","B2", "C2")
        } else if (matched_model=="double gompertz"){
          form_fun<-.brms_form_dou_gompertz
          pars=c("A","B", "C", "A2","B2", "C2")
        } else if(matched_model=="logistic"){
          form_fun<-.brms_form_logistic
          pars=c("A","B", "C")
        } else if (matched_model=="gompertz"){
          form_fun<-.brms_form_gompertz
          pars=c("A","B", "C")
        } else if (matched_model=="monomolecular"){
          form_fun<-.brms_form_monomolecular
          pars=c("A","B")
        } else if (matched_model=="exponential"){
          form_fun<-.brms_form_exponential
          pars=c("A","B")
        } else if (matched_model=="linear"){
          form_fun<-.brms_form_linear
          pars=c("A")
        } else if (matched_model=="power law"){
          form_fun<-.brms_form_powerlaw
          pars=c("A","B")
        } else if (matched_model=="gam"){
          form_fun<-.brms_form_gam
          pars <- NULL
        }
        growthForm = form_fun(x,y, group) # group only used in gam model
    }
    #* `Make parameter formula`
    #* could add a pars argument then set up parForm from those.
    #* I think that would change how priors would have to work
    #* and that seems like more trouble than it is worth right 
    #* now at least.
    if(matched_sigma == "gompertz"){ # add nl pars for sigma form
      pars <- c(pars, "subA", "subB", "subC")
    }
    if(!is.null(pars)){
      if(USEGROUP){ parForm<-as.formula(paste0( paste(pars,collapse="+"),"~0+",group ))
      } else { parForm<-as.formula(paste0( paste(pars,collapse="+"),"~1" )) }
    } else {parForm=NULL} 
    
    if(!is.null(spline_pars)){
      splineParForm <- as.formula(paste0(spline_pars, " ~ s(", x, ", by=",group, ")"))
    }
    
    #* `Make heteroskedasticity formula`
    if(!is.null(sigma)){
      sigmaForm<-if(matched_sigma=="homo"){
        if(USEGROUP){ as.formula(paste0("sigma ~ 0+", group))
          } else{ as.formula(paste0("sigma ~ 1")) }
        } else if (matched_sigma=="linear"){
          if(USEGROUP){
            as.formula(paste0("sigma ~ ", x, "+", x, ":",group))
          } else{ as.formula(paste0("sigma ~ ", x)) }
        } else if(matched_sigma=="gompertz"){
          brms::nlf(as.formula(paste0("sigma ~ subA*exp(-subB*exp(-subC*",x,"))")))
          } else if(matched_sigma=="spline"){
            
        if(length(unique(df[[x]]))<11){
          if(USEGROUP){ as.formula(paste0("sigma ~ s(",x,", by=", group, ", k=",length(unique(df[[x]])),")"))
            }else{ as.formula(paste0("sigma ~ s(",x,", k=",length(unique(df[[x]])),")")) }
        }else{ 
          if(USEGROUP){ as.formula(paste0("sigma ~ s(",x,", by=", group, ")"))
          } else{ as.formula(paste0("sigma ~ s(",x, ")")) }
        }
      } 
      #* `Combining for brms formula`
      #* might be better to do this with do.call, these are getting cumbersome.
      if(!is.null(parForm)){
        if(!is.null(spline_pars)){
          bayesForm<-brms::bf(formula = growthForm, sigmaForm, parForm, splineParForm, autocor = corForm, nl=TRUE)
        } else{bayesForm<-brms::bf(formula = growthForm, sigmaForm, parForm, autocor = corForm, nl=TRUE)}
      } else{
        bayesForm<-brms::bf(formula = growthForm, sigmaForm, autocor = corForm)
      }
    }else{
      if(!is.null(parForm)){
        if(!is.null(spline_pars)){
          bayesForm<-brms::bf(formula = growthForm, parForm, splineParForm, autocor = corForm, nl=TRUE)
        }else{bayesForm<-brms::bf(formula = growthForm, parForm, autocor = corForm, nl=TRUE)}
      } else{
        bayesForm<-brms::bf(formula = growthForm, autocor = corForm)
      }
    }
    out[["formula"]]<-bayesForm
  #* ***** `Make priors` *****
  #groupedPriors = TRUE # default this to TRUE in case of brmsprior class
  if(!is.null(priors)){
    # priors = c(150, 30, 0.125)
    # priors = list("A"=130, "B"=12, "C"=3)
    # priors = list("A"=c(130, 150), "B"=c(12,11), "C"=c(3,3))
    #* `if priors is a brmsprior`
    if( any(methods::is(priors, "brmsprior")) ){out[["prior"]]<-priors } else{
    #* `if priors is a numeric vector`
    if(is.numeric(priors)){ # priors = c(130, 12, 3) 
      if(length(priors)==length(pars)){
      warning("Assuming that prior is in order: ", paste0(pars, collapse=', '))
        priors<-as.list(priors)
        names(priors)<-pars
      } else {stop(paste0("`priors` is length ", length(priors), " while the specified model requires ", length(pars), " parameters."))}
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
      groupedPriors<-any(unlist(lapply(priors,length))>1) # if any prior has multiple means then groupedPriors is TRUE
      
      if(groupedPriors){ # if more than one value is specified per parameter
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
        # this should also handle non-grouped formulae
        l<-length(unique(df[[group]]))
        priors<-lapply(priors, rep, length.out=l)
        nms<-unique(df[[group]])
        if(USEGROUP){
         for(i in 1:length(priors)){names(priors[[i]])<-nms}
        }
      }
    }
      
    }
    if(!any(names(out)=="prior")){
      priorStanStrings<-lapply(pars, function(par) paste0("lognormal(log(", priors[[par]],"), 0.25)") )
      priorStanStrings<-unlist(priorStanStrings)
      parNames<-rep(names(priors), each = length(priors[[1]]))
      if(USEGROUP){
        groupNames<-rep(names(priors[[1]]), length.out = length(priorStanStrings))
        names(priorStanStrings)<-paste(parNames, groupNames, sep="_")
      } else {
        names(priorStanStrings)<-parNames
      }
      
      if(matched_sigma =="homo" & !USEGROUP){
        prior<-brms::set_prior('student_t(3,0,5)', dpar="sigma", class="Intercept")+
          brms::set_prior('gamma(2,0.1)', class="nu")
      } else if(matched_sigma == "gompertz"){
        prior<-brms::set_prior('gamma(2,0.1)', class="nu")
      } else{
        prior<-brms::set_prior('student_t(3,0,5)', dpar="sigma")+
          brms::set_prior('gamma(2,0.1)', class="nu") 
      }
      
      for(nm in names(priorStanStrings)){
        dist = priorStanStrings[[nm]]
        pr = strsplit(nm, "_")[[1]][1]
        if(USEGROUP & groupedPriors){ # if there are groups and they have different priors
          gr = paste0(group, strsplit(nm, "_")[[1]][2])
          prior<-prior+brms::set_prior(dist, coef = gr, nlpar = pr) # currently cannot set lb for prior with coef
          # there is a clunky workaround but it wouldn't work with expected data types
          # https://github.com/paul-buerkner/brms/issues/86
        } else{
          prior<-prior+brms::set_prior(dist, nlpar = pr, lb = 0)
        }
      }
      prior <-unique(prior)
      out[["prior"]]<-prior
      # priors is still a list with ln centers
    }
  }
    
  #* ***** `Make initializer function` *****
    if(!is.null(pars)){
    initFun<-function(pars="?", nPerChain=1){
      init<-lapply(pars, function(i) array(rgamma(nPerChain,1)))
      names(init)<-paste0("b_",pars)
      init
    }
    formals(initFun)$pars<-pars
    formals(initFun)$nPerChain<-length(unique(df[[group]]))
    wrapper<-function(){initFun()}
    } else{
      wrapper = 0
    }
    
    out[["initfun"]]<-wrapper
    out[["df"]]<-df
    out[["family"]]<-"student"
    out[["pcvrForm"]]<-form
  return(out)
  }
  
.brms_form_logistic<-function(x, y, group){
  return(stats::as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C))")))
}
.brms_form_gompertz<-function(x, y, group){
  return(stats::as.formula(paste0(y," ~ A*exp(-B*exp(-C*",x,"))")))
}
.brms_form_dou_logistic<-function(x, y, group){
  return(stats::as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C)) + ((A2-A) /(1+exp((B2-",x,")/C2)))")))
}
.brms_form_dou_gompertz<-function(x, y, group){
  return(stats::as.formula(paste0(y," ~ A * exp(-B * exp(-C*",x,")) + (A2-A) * exp(-B2 * exp(-C2*(",x,"-B)))")))
}
.brms_form_monomolecular<-function(x, y, group){
  return(stats::as.formula(paste0(y,"~A-A*exp(-B*",x,")")))
}
.brms_form_exponential<-function(x, y, group){
  return(stats::as.formula(paste0(y," ~ A*exp(B*",x, ")")))
}
.brms_form_linear<-function(x, y, group){
  return(stats::as.formula(paste0(y," ~ A*",x)))
}
.brms_form_powerlaw<-function(x, y, group){
  return(stats::as.formula(paste0(y," ~ A*",x,"^B")))
}
.brms_form_gam<-function(x, y, group){
  return(stats::as.formula(paste0(y," ~ s(",x,", by=",group, ")")))
}









