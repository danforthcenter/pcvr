#' Ease of use growth model helper function for 6 model parameterizations
#' 
#' @param model The name of a model as a character string.
#' Supported options are c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law", "double logistic", "double gompertz").
#' See \code{\link{growthSim}} for examples of each type of growth curve.
#' @param form A formula describing the model. The left hand side should only be 
#' the outcome variable (phenotype). The right hand side needs at least the x variable
#'  (typically time). Grouping is also described in this formula using roughly lme4
#'  style syntax,with formulas like \code{y~time|individual/group} to show that predictors
#'  should vary by \code{group} and autocorrelation between \code{individual:group}
#'  interactions should be modeled. If group has only one level or is not included then
#'  it will be ignored in formulas for growth and variance (this may be the case if 
#'  you split data before fitting models to be able to run more smaller models each more quickly).
#' @param sigma A model for heteroskedasticity. This argument is only used with "brms" and "nlme" models.
#' When type="brms" this is turned into a formula and the supported options are
#' c("homo", "linear", "spline", "logistic", "gompertz"), if left NULL then "spline" will be used.
#'  When type ="nlme" the options are more limited to c("none", "power", "exp"),
#'   corresponding to using \code{nlme::varIdent}, \code{nlme::varPower},
#' or \code{nlme::varExp} respectively where "power" is the default.
#' @param df A dataframe to use. Must contain all the variables listed in the formula.
#' @param pars Optionally specify which parameters should change by group. Not this is model
#' dependent and is not implemented for brms models due to their more flexible hypothesis testing.
#' @param start An optional named list of starting values OR means for prior distributions.
#' If this is not provided then starting values are picked with \code{stats::selfStart}.
#'  When type = "brms" these must be provided and are treated as the means of
#'   lognormal priors for all growth model parameters.
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
#'   the necessary features are added to \code{brms}. See details for guidance.
#' @param type Type of model to fit, options are "brms", "nlrq", and "nls".
#' @param tau A vector of quantiles to fit for nlrq models.
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
#' Currently there are four supported brms sub models described below.
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
#' There are also three supported submodel options for nlme models, but a \code{varFunc} object can also be supplied,
#' see \code{?nlme::varClasses}.
#' 
#' \itemize{
#'    \item \bold{homo}: \code{varIdent(1|group)}, which models a constant variance separately for each group. 
#'     \item \bold{power}: \code{varPower(x|group)}, which models variance as a power of x per group.
#'     \item \bold{exp}: \code{varExp(x|group)}, which models variance as an exponent of x per group.
#' }
#' 
#' 
#' @return A named list of elements to make it easier to fit non linear growth models with several R packages.
#' 
#' For \code{brms} models the output contains:
#' 
#' \code{formula}: A \code{brms::bf} formula specifying the growth model, autocorrelation, variance submodel,
#' and models for each variable in the growth model.
#' \code{prior}: A brmsprior/data.frame object.
#' \code{initfun}: A function to randomly initialize chains using a random draw from a gamma
#' distribution (confines initial values to positive and makes correct number
#' of initial values for chains and groups).
#' \code{df} The data input, possibly with dummy variables added if needed.
#' \code{family} The model family, currently this will always be "student".
#' \code{pcvrForm} The form argument unchanged. This is returned so that
#' it can be used later on in model visualization. Often it may be a good idea
#' to save the output of this function with the fit model, so having this can
#' be useful later on.
#' 
#' For \code{quantreg::nlrq} models the output contains:
#' 
#' \code{formula}: An \code{nls} style formula specifying the growth model with groups if specified.
#' \code{taus}: The quantiles to be fit
#' \code{start}: The starting values, typically these will be generated from the growth model and your data in a similar way as shown in \code{stats::selfStart} models.
#' \code{df} The input data for the model.
#' \code{pcvrForm} The form argument unchanged.
#' 
#' For \code{nls} models the output is the same as for \code{quantreg::nlrq} models but without \code{taus} returned.
#' 
#' For \code{nlme::nlme} models the output contains:
#' 
#' \code{formula}: An list of \code{nlme} style formulas specifying the model, fixed and random effects, random effect grouping, and variance model (weights).
#' \code{start}: The starting values, typically these will be generated from the growth model and your data in a similar way as shown in \code{stats::selfStart} models.
#' \code{df} The input data for the model.
#' \code{pcvrForm} The form argument unchanged.
#' 
#' @examples 
#' 
#' ## Not run:
#' 
#' simdf<-growthSim("logistic", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "logistic", form=y~time|id/group,
#'              sigma="spline", df=simdf,
#'              start = list("A"=130, "B"=12, "C"=3), type="brms")
#' lapply(ss,class)
#' ss$initfun()
#' 
#' if(FALSE){
#' fit_test <- fitGrowth(ss, iter = 1000, cores = 2, chains = 2, backend = "cmdstanr",
#'  control = list(adapt_delta = 0.999, max_treedepth = 20))
#' }
#' 
#' 
#' # formulas and priors will look different if there is only one group in the data
#' 
#' ex<-growthSim("linear", n=20, t=25, params=list("A"=2))
#' ex_ss <- growthSS(model = "linear", form = y~time|id/group, sigma = "spline",
#'                 df = ex, start = list("A" = 1), type="brms")
#'                 
#' ex_ss$prior # no coef level grouping for priors
#' ex_ss$formula # intercept only model for A
#' 
#' ex2<-growthSim("linear", n=20, t=25, params=list("A"=c(2, 2.5)))
#' ex2_ss <- growthSS(model = "linear", form = y~time|id/group, sigma = "spline",
#'                    df = ex2, start = list("A" = 1), type="brms")
#' ex2_ss$prior # has coef level grouping for priors
#' ex2_ss$formula # specifies an A intercept for each group and splines by group for sigma
#' 
#' ## End(Not run)              
#'               
#' @export

growthSS<-function(model, form, sigma=NULL, df, start=NULL, pars=NULL, type="brms", tau = 0.5){
  type_matched = match.arg(type, choices = c("brms", "nlrq", "nls", "nlme"))
  if(type_matched=="brms"){
    if(is.null(sigma)){sigma="spline"}
    res <- .brmSS(model=model, form=form, sigma=sigma, df=df, priors = start)
  } else if(type_matched %in% c("nlrq", "nls")){
    res <- .nlrqSS(model=model, form=form, tau=tau, df=df, start=start, pars=pars, type=type)
  } else if(type_matched=="nlme"){
    if(is.null(sigma)){sigma="power"}
    res <- .nlmeSS(model=model, form=form, sigma=sigma, df=df, pars=pars, start=start)
  } else{stop("Type must match one of brms, nlrq, nls, or nlme")}
  res$type = type
  res$call = match.call()
  return(res)
}


