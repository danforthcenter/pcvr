#' Ease of use growth model helper function. Output from this should be passed to \link{fitGrowth} to
#' fit the specified model.
#'
#' @param model The name of a model as a character string.
#'  Supported options are c("logistic", "gompertz", "weibull", "frechet", "gumbel", "monomolecular",
#'  "exponential", "linear", "power law",
#'  "double logistic", "double gompertz", "gam", "int"), with "int" representing an intercept only model
#'  which is only used in brms (and is expected to only be used in threshold models or to model
#'  homoskedasticity.)
#'  See \code{\link{growthSim}} for examples of each type of single parameterized growth curve
#'  ("gam" is not supported in \code{growthSim}).
#'  You can also specify decay models by including the "decay" keyword with the model name.
#'  Note that in brms models the entire formula is negated for decay models so that lognormal priors can
#'  still be used when at least some coefficients would be negative.
#'  With type="brms" you can also specify segmented models by combining model names with a plus sign
#'  such as "linear + linear". In a segmented model the names for parameters do not follow the normal
#'  "A", "B", "C" notation, instead they are named for the type of model, the position in the formula,
#'  then for the parameter of that model. There will also be parameters to represent the time when
#'  growth switches from one model to another called either "changepointX" or "fixedChangePointX".
#'  All "changePointX" terms are estimated as parameters of the model.
#'  "fixedChangePointX" parameters are not estimated and are kept as the numeric value given in the
#'  priors, this is useful if your experiment has an intervention at a set time which you expect to
#'  change the growth process acutely.
#'  For the "linear + linear" example this would yield parameters "linear1A", "changePoint1"
#'  (or "fixedChangePoint1"), and "linear2A". A "linear + gompertz" model would have
#'  "linear1A", "changePoint1", "gompertz2A", "gompertz2B", and "gompertz2C" for parameters.
#'  Note that double sigmoid models are not supported as parts of segmented models and gams
#'  can currently only be included as the last part of a segmented model. When using a changepoint model
#'  it may be worth using segments that are simpler to fit
#'  (gompertz instead of EVD options, for instance). Currently "homo" and "int" are treated the same
#'  and "spline" and "gam" are interchangeable. Time-to-event models can be specified using the
#'  "survival" keyword, see details for an explanation of the changes that entails.
#'  Similarly, using the brms backend response distributions (see \code{brms::brmsfamily})
#'  can be specified in the model as "family: model" so that a model
#'  of logistic increasing counts may be written as \code{model = "poisson: logistic"}.
#' @param form A formula describing the model. The left hand side should only be
#' the outcome variable (phenotype), and a cutoff if you are making a survival model (see details).
#' The right hand side needs at least the x variable
#'  (typically time). Grouping is also described in this formula using roughly lme4
#'  style syntax,with formulas like \code{y~time|individual/group} to show that predictors
#'  should vary by \code{group} and autocorrelation between \code{individual:group}
#'  interactions should be modeled. Note that autocorrelation is only modeled with the "brms"
#'  backend in this way. "nlme" requires random effects and correlations to use the same grouping,
#'  so autocorrelation using the "nlme" backend works at the group level, so will slightly underestimate
#'  the autocorrelation at the individual level. If group has only one level or is not included then
#'  it will be ignored in formulas for growth and variance (this may be the case if
#'  you split data before fitting models to be able to run more smaller models each more quickly).
#' @param sigma Other models for distributional parameters.
#' This argument is only used with "brms" and "nlme" models and is handled differently for each.
#' When type="brms" this can be supplied as a model or as a list of models.
#' It is turned into a formula (or list of formulas) with an entry corresponding to each distributional
#' parameter (after the mean) of the growth model family.
#' If no family was specified (\code{model="logistic"} for instance) then the student T distribution
#' is used, with additional distributional parameters sigma and nu.
#' To check the naming of distributional parameters in each response family use
#' \code{brms::brmsfamily("family")$dpars}. The supported options are the same as the model options
#' (including threshold models).
#' For distributional parameters that do not have a formula specified they will be modeled as
#' intercept only (not by group).
#' Parameter names are the same as those in the main model but with the distributional parameter name
#' as a prefix. Additionally, if a linear model is used for sigma then it can be modeled with or without
#' a prior, if a prior is specified ("sigmaA") then a non-linear formula is used and the "sigmaA"
#' parameter will be included in the output instead of the default "sigma" term.
#' In the rare case that you wish to model the mean and the 3rd distributional parameter but not
#' the 2nd then \code{sigma = list("not_estimated", "model")} would allow for that.
#'  When type ="nlme" the options are more limited to c("none", "power", "exp"), corresponding to using
#' \code{nlme::varIdent}, \code{nlme::varPower}, or \code{nlme::varExp} respectively where "power"
#' is the default.
#' @param df A dataframe to use. Must contain all the variables listed in the formula.
#' @param pars Optionally specify which parameters should change by group. Not this is model
#' dependent and is not implemented for brms models due to their more flexible hypothesis testing.
#' @param start An optional named list of starting values OR means for prior distributions.
#' If this is not provided then starting values are picked with \code{stats::selfStart}.
#'  When type = "brms" these must be provided and are treated as the means of
#'   lognormal priors for all growth model parameters and T_5(mu, 3) priors for changepoint parameters.
#'   This is done because the values are strictly positive and the lognormal distribution
#'   is easily interpreted. The changepoint priors are T distributions for symmetry, 5 DF
#'   having been chosen for heavy but not unmanageable tails.
#'   If this argument is not provided then priors are not
#'   returned and a different set of priors will need to be made for the model using
#'   \code{brms::set_prior}. When specifying starting values/prior means
#'   think of this as being similar to the \code{params} argument
#'   in \code{growthSim}. Names should correspond to parameter names from the
#'   \code{model} argument. A numeric vector can also be used, but specifying
#'   names is best practice for clarity. Additionally, due to a limitation in
#'   \code{brms} currently lower bounds cannot be set for priors for specific groups.
#'   If priors include multiple groups (\code{start = list(A = c(10,15), ...)}) then
#'   you will see warnings after the model is fit about not having specified a lower
#'   bound explicitly. Those warnings can safely be ignored and will be addressed if
#'   the necessary features are added to \code{brms}. See details for guidance.
#' @param type Type of model to fit, options are "brms", "nlrq", "nlme", "nls", and "mgcv".
#' Note that the "mgcv" option only supports "gam" models.
#' Survival models can use the "survreg" model type
#' (this will be called if any non-brms/flexsurv type is given) or the "flexsurv" model type
#' which requires the flexsurv package to be installed.
#' Note that for non-brms models variables in the model will be labeled by the factor level of the
#' group, not necessarily by the group name.
#' This is done for ease of use with different modeling functions, the levels are alphabetically sorted
#' and can be checked using:
#' \code{table(ss$df$group, ss$df$group_numericLabel)}.
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
#' The \link{plotPrior} function can be used to do prior predictive checks.
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
#'     \item \bold{Weibull}: \code{list('A' = 130, 'B' = 2, 'C' = 2)}
#'     \item \bold{Frechet}: \code{list('A' = 130, 'B' = 5, 'C' = 6)}
#'     \item \bold{Gumbel}: \code{list('A' = 130, 'B' = 6, 'C' = 4)}
#'     \item \bold{Double Logistic}: \code{list('A' = 130, 'B' = 12, 'C' = 3,
#'     'A2' = 200, 'B2' = 25, 'C2' = 1)}
#'     \item \bold{Double Gompertz}: \code{list('A' = 130, 'B' = 12, 'C' = 0.25,
#'     'A2' = 220, 'B2' = 30, 'C2' = 0.1)}
#'     \item \bold{Monomolecular}: \code{list('A' = 130, 'B' = 2)}
#'     \item \bold{Exponential}: \code{list('A' = 15, 'B' = 0.1)}
#'     \item \bold{Linear}: \code{list('A' = 1)}
#'     \item \bold{Power Law}: \code{list('A' = 13, 'B' = 2)}
#' }
#'
#' See details below about parameterization for each model option.
#' \itemize{
#'  \item \bold{Logistic}: `A / (1 + exp( (B-x)/C) )`
#'  Where A is the asymptote, B is the inflection point, C is the growth rate.
#'  \item \bold{Gompertz}: `A * exp(-B * exp(-C*x))`
#'  Where A is the asymptote, B is the inflection point, C is the growth rate.
#'  \item \bold{Weibull}: `A * (1-exp(-(x/C)^B))`
#'  Where A is the asymptote, B is the weibull shape parameter, C is the weibull scale parameter.
#'  \item \bold{Frechet}: `A * exp(-((x-0)/C)^(-B))`
#'  Where A is the asymptote, B is the frechet shape parameter, C is the frechet scale parameter.
#'  Note that the location parameter (conventionally m) is 0 in these models for simplicity but is still
#'  included in the formula.
#'  \item \bold{Gumbel}: `A * exp(-exp(-(x-B)/C))`
#'  Where A is the asymptote, B is the inflection point (location), C is the growth rate (scale).
#'  \item \bold{Double Logistic}: `A / (1+exp((B-x)/C)) + ((A2-A) /(1+exp((B2-x)/C2)))`
#'  Where A is the asymptote, B is the inflection point, C is the growth rate,
#'  A2 is the second asymptote, B2 is the second inflection point, and C2 is the second
#'  growth rate.
#'  \item \bold{Double Gompertz}: `A * exp(-B * exp(-C*x)) + ((A2-A) * exp(-B2 * exp(-C2*(x-B))))`
#'  Where A is the asymptote, B is the inflection point, C is the growth rate,
#'  A2 is the second asymptote, B2 is the second inflection point, and C2 is the second
#'  growth rate.
#'  \item \bold{Monomolecular}: `A-A * exp(-B * x)`
#'  Where A is the asymptote and B is the growth rate.
#'  \item \bold{Exponential}: `A * exp(B * x)`
#'  Where A is the scale parameter and B is the growth rate.
#'  \item \bold{Linear}: `A * x`
#'  Where A is the growth rate.
#'  \item \bold{Power Law}: `A * x^(B)`
#'  Where A is the scale parameter and B is the growth rate.
#'  }
#'  Note that for these distributions parameters do not exist in a vacuum.
#'  Changing one will make the others look different in the resulting data.
#'  The \code{growthSim} function can be helpful in familiarizing further with these distributions.
#'
#' Using the \code{brms} backend the \code{sigma} argument optionally specifies a sub model to account
#' for heteroskedasticity.
#' Any combination of models (except for decay models) can be specified in the \code{sigma} term.
#' If you need variance to raise and lower then a gam/spline is the most appropriate option.
#'
#' Using the \code{brms} backend a model with lots of parameters may be difficult to estimate if there
#' are lots of groups.
#' If you have very many levels of your "group" variable in a complex model then consider fitting models
#' to subsets of the "group" variable and using \link{combineDraws} to make a data.frame for
#' hypothesis testing.
#'
#' Limits on the Y variable can be specified in the \code{brms} backend. This should generally be
#' unnecessary and will make the model slower to fit and potentially more difficult to set priors on.
#' If you do have a limited phenotype (besides the normal positive constraint for growth models)
#' then this may be helpful, one situation may be canopy coverage percentage which is naturally bounded
#' at an upper and lower limit.
#' To specify these limits add square brackets to the Y term with upper and lower limits such as
#' \code{"y[0,100] ~ time|id/group"}.
#'
#' There are also three supported submodel options for \code{nlme} models, but a \code{varFunc} object
#' can also be supplied, see \code{?nlme::varClasses}.
#'
#' \itemize{
#'    \item \bold{none}: \code{varIdent(1|group)}, which models a constant variance separately for each
#'    group.
#'     \item \bold{power}: \code{varPower(x|group)}, which models variance as a power of x per group.
#'     \item \bold{exp}: \code{varExp(x|group)}, which models variance as an exponent of x per group.
#' }
#'
#' Survival models can be fit using the "survival" keyword in the model specification.
#' Using the "brms" backend (type argument) you can specify "weibull" (the default) or "binomial" for
#' the distribution to use in that model so that the final model string would be "survival binomial" or
#' "survival weibull" which is equivalent to "survival". Time to event data is very different than
#' standard phenotype data, so the formula argument should include a cutoff for the Y variable to count
#' as an "event". For example, if you were checking germination using area and wanted to use 50 pixels
#' as a germinated plant your formula would be \code{area > 50 ~ time|id/group}.
#' Internally the input dataframe will be converted to time-to-event data based on that formula.
#' Alternatively you can make your own time to event data and supply that to growthSS. In that case your
#' data should have columns called "n_events"
#' (number of individuals experiencing the event at this time) and "n_eligible"
#' (number of individuals who had not experienced the event at least up to this time)
#' for the binomial model family OR "event" (binary 1,0 for TRUE, FALSE) for the Weibull model family.
#' Note that since these are linear models using different model families the priors are handled
#' differently. For survival models the default priors are weak regularizing priors (Normal(0,5))
#' on all parameters. If you wish to specify your own priors you can supply them as brmsprior objects
#' or as a list such as \code{priors = list("group1" = c(0,3), "group2" = c(0,1))} where the order of
#' values is Mu, Sigma.
#' Any non-brms backend will instead use \code{survival::survreg} to fit the model unless the
#' "flexsurv" type is specified.
#' Distributions will be passed to \code{survreg} where options are "weibull", "exponential",
#' "gaussian", "logistic","lognormal" and "loglogistic" if type = "survreg" or to
#' \code{flexsurv::flexsurvreg} if type = "flexsurv" where options are "gengamma", "gengamma.orig",
#' "genf", "genf.orig", "weibull", "gamma", "exp", "llogis", "lnorm", "gompertz", "exponential",
#' and "lognormal". In \code{flexsurvreg} distributional modeling is supported and additional
#' formula can be passed as a list to the sigma argument of growthSS in the same way as to the anc
#' argument of \code{flexsurv::flexsurvreg}.
#' Further additional arguments should be supplied via \code{fitGrowth} if desired.
#'
#'
#'
#' @return A named list of elements to make it easier to fit non linear growth models with several R
#' packages.
#'
#' For \code{brms} models the output contains:
#'
#' \code{formula}: A \code{brms::bf} formula specifying the growth model, autocorrelation,
#' variance submodel, and models for each variable in the growth model.
#' \code{prior}: A brmsprior/data.frame object.
#' \code{initfun}: A function to randomly initialize chains using a random draw from a gamma
#' distribution (confines initial values to positive and makes correct number
#' of initial values for chains and groups).
#' \code{df} The data input, with dummy variables added if needed and a column to link groups to their
#' factor levels.
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
#' \code{start}: The starting values, typically these will be generated from the growth model and your
#' data in a similar way as shown in \code{stats::selfStart} models.
#' \code{df} The input data for the model.
#' \code{pcvrForm} The form argument unchanged.
#'
#' For \code{nls} models the output is the same as for \code{quantreg::nlrq} models but without
#' \code{taus} returned.
#'
#' For \code{nlme::nlme} models the output contains:
#'
#' \code{formula}: An list of \code{nlme} style formulas specifying the model, fixed and random effects,
#' random effect grouping, and variance model (weights).
#' \code{start}: The starting values, typically these will be generated from the growth model and your
#' data in a similar way as shown in \code{stats::selfStart} models.
#' \code{df} The input data for the model.
#' \code{pcvrForm} The form argument unchanged.
#'
#' For all models the type and model are also returned for simplicity downstream.
#'
#' @examples
#'
#' ## Not run:
#'
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group,
#'   sigma = "spline", df = simdf,
#'   start = list("A" = 130, "B" = 12, "C" = 3), type = "brms"
#' )
#' lapply(ss, class)
#' ss$initfun()
#'
#' if (FALSE) {
#'   fit_test <- fitGrowth(ss,
#'     iter = 1000, cores = 2, chains = 2, backend = "cmdstanr",
#'     control = list(adapt_delta = 0.999, max_treedepth = 20)
#'   )
#' }
#'
#'
#' # formulas and priors will look different if there is only one group in the data
#'
#' ex <- growthSim("linear", n = 20, t = 25, params = list("A" = 2))
#' ex_ss <- growthSS(
#'   model = "linear", form = y ~ time | id / group, sigma = "spline",
#'   df = ex, start = list("A" = 1), type = "brms"
#' )
#'
#' ex_ss$prior # no coef level grouping for priors
#' ex_ss$formula # intercept only model for A
#'
#' ex2 <- growthSim("linear", n = 20, t = 25, params = list("A" = c(2, 2.5)))
#' ex2_ss <- growthSS(
#'   model = "linear", form = y ~ time | id / group, sigma = "spline",
#'   df = ex2, start = list("A" = 1), type = "brms"
#' )
#' ex2_ss$prior # has coef level grouping for priors
#' ex2_ss$formula # specifies an A intercept for each group and splines by group for sigma
#'
#' ## End(Not run)
#'
#' @export

growthSS <- function(model, form, sigma = NULL, df, start = NULL,
                     pars = NULL, type = "brms", tau = 0.5) {
  type_matched <- match.arg(type, choices = c(
    "brms", "nlrq", "nls",
    "nlme", "mgcv", "survreg",
    "flexsurv"
  ))

  surv <- .survModelParser(model)
  survivalBool <- surv$survival
  model <- surv$model

  if (survivalBool) {
    if (type_matched == "brms") {
      res <- .brmsSurvSS(model = model, form = form, df = df, priors = start)
      res$type <- type_matched
    } else if (type_matched == "flexsurv") {
      res <- .flexSurvSS(model = model, form = form, df = df, anc = sigma)
      res$type <- "flexsurv"
    } else {
      res <- .survSS(model = model, form = form, df = df)
      res$type <- "survreg"
    }
  } else {
    if (type_matched == "brms") {
      res <- .brmSS(model = model, form = form, sigma = sigma, df = df, priors = start)
    } else if (type_matched %in% c("nlrq", "nls")) {
      res <- .nlrqSS(
        model = model, form = form, tau = tau, df = df, start = start, pars = pars,
        type = type
      )
    } else if (type_matched == "nlme") {
      if (is.null(sigma)) {
        sigma <- "power"
      }
      res <- .nlmeSS(model = model, form = form, sigma = sigma, df = df, pars = pars, start = start)
    } else if (type_matched == "mgcv") {
      res <- .mgcvSS(model = model, form = form, df = df)
    } else {
      stop("Type must match one of brms, nlrq, nls, nlme, or mgcv")
    }
    res$type <- type
  }
  res$model <- model
  res$call <- match.call()
  return(res)
}
