#' Ease of use brms starter function for 6 growth model parameterizations
#'
#' @param model The name of a model as a character string.
#' Supported options are "logistic", "gompertz", "frechet", "gumbel", "weibull",
#' "monomolecular", "exponential", "linear", "logarithmic",
#' "power law", "double logistic", "double gompertz", and "gam".
#' See \code{\link{growthSim}} for examples of each type of growth curve.
#' @param form A formula describing the model. The left hand side should only be
#' the outcome variable (phenotype). The right hand side needs at least the x variable
#'  (typically time). Grouping is also described in this formula using roughly lme4
#'  style syntax,with formulas like \code{y~time|individual/group} to show that predictors
#'  should vary by \code{group} and autocorrelation between \code{individual:group}
#'  interactions should be modeled. If group has only one level or is not included then
#'  it will be ignored in formulas for growth and variance (this may be the case if
#'  you split data before fitting models to be able to run more smaller models each more quickly).
#' @param sigma A model for heteroskedasticity from the same list of options as the model argument.
#' @param df A dataframe to use. Must contain all the variables listed in the formula.
#' @param priors A named list of means for prior distributions.
#'  Currently this function makes lognormal priors for all growth model parameters
#'  and T_5(mu, 3) priors for changepoint parameters.
#'   This is done because the values are strictly positive and the lognormal distribution
#'   is easily interpreted. The changepoint priors are T distributions for symmetry, 5 DF
#'   having been chosen for heavy but not unmanageable tails.
#'   If this argument is not provided then priors are not
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
#' @param int Logical, should an intercept term be included?
#' @keywords Bayesian, brms
#'
#' @importFrom stats as.formula rgamma
#'
#' @details
#'
#' Default informative priors are not provided,
#' but these can serve as starting points for each distribution.
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
#'
#'
#' The \code{sigma} argument optionally specifies a sub model to account for heteroskedasticity.
#' Currently there are four supported sub models described below.
#'
#' \itemize{
#'    \item \bold{homo}: \code{sigma ~ 1}, fitting only a global or per group intercept to sigma.
#'     \item \bold{linear}: \code{sigma ~ time}, modeling sigma with a linear relationship to time
#'      and possibly with an interaction term per groups.
#'     \item \bold{spline}: \code{sigma ~ s(time)}, modeling sigma using a smoothing function through
#'     `mgcv::s`, possibly by group.
#'     \item \bold{gompertz}: \code{sigma ~ subA * exp(-subB * exp(-subC * x))},
#'      modeling sigma as a gompertz function of time, possibly by group. Note that you
#'      should specify priors for the parameters in this sub model by adding them into the \code{priors}
#'      argument, such as \code{list(..., subA = 20, subB = 15, subC = 0.25)}. If you do not specify
#'      priors then default (flat) priors will be used, which is liable to cause fitting problems
#'      and less accurate results. Looking at your data and making a semi-informed estimate of the
#'      total variance at the end of the experiment can help set a reasonable prior for subA, while
#'      subB and subC can generally be the same as B and C in a gompertz growth model of the same data.
#'      These priors will have fat tails so they are pretty forgiving.
#' }
#'
#'
#' @return A named list of elements to make it easier to fit common brms models.
#' \code{formula}: A \code{brms::bf} formula specifying the growth model, autocorrelation, variance
#' submodel, and models for each variable in the growth model.
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
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- .brmSS(
#'   model = "logistic", form = y ~ time | id / group,
#'   sigma = "spline", df = simdf, priors = list("A" = 130, "B" = 12, "C" = 3)
#' )
#' lapply(ss, class)
#' ss$initfun()
#'
#' @keywords internal
#' @noRd

.brmSS <- function(model, form, sigma = NULL, df, priors = NULL, int = FALSE, hierarchy = NULL) {
  out <- list()
  models <- c(
    "int", "logistic", "gompertz", "monomolecular", "exponential", "linear", "power law", "logarithmic",
    "double logistic", "double gompertz", "gam", "spline", "homo", "frechet", "gumbel", "weibull",
    "not_estimated", "bragg", "lorentz", "beta"
  )
  #* ***** `Make bayesian formula` *****
  #* `parse form argument`
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  if (grepl("\\[", y)) {
    lims <- gsub(".*\\[|\\]", "", y)
    lims <- as.numeric(strsplit(lims, ",")[[1]])
    lower <- lims[1]
    upper <- lims[2]
    y <- trimws(gsub("\\[.*", paste0("|trunc(lb=", lower, ", ub=", upper, ")"), y))
  }
  x <- parsed_form$x
  individual <- parsed_form$individual
  group <- parsed_form$group
  USEGROUP <- parsed_form$USEG
  USEINDIVIDUAL <- parsed_form$USEID
  df <- parsed_form$data
  hierarchical_predictor <- parsed_form$hierarchical_predictor

  #* `convert group to character to avoid unexpected factor stuff`
  df[, group] <- lapply(group, function(grp) {
    as.character(df[[grp]])
  })
  #* `if there are gams involved and multiple groups then make a group interaction variable`
  if (length(group) > 1 && any(grepl("spline|gam", c(model, sigma)))) {
    df[[paste(group, collapse = ".")]] <- interaction(df[, group])
  }
  #* `Make autocorrelation formula`
  if (USEINDIVIDUAL) {
    corForm <- as.formula(paste0("~arma(~", x, "|",
                                 paste(c(individual, group), collapse = ":"),
                                 ",1,1)"))
  } else {
    corForm <- NULL
  }
  #* `get family specific elements`
  family_res <- .brmFamilyHelper(model)
  model <- family_res$rhs
  family <- family_res$family
  dpars <- family_res$dpars
  #* `Reformat sigma if it is a formula`
  sigma <- .sigmaHelper(sigma, dpars, family, models)
  #* `match args`
  matchGrowthModelRes <- .matchGrowthModel(model, models)
  matched_model <- matchGrowthModelRes[["model"]]
  decay <- matchGrowthModelRes[["decay"]]
  #* `Make growth formula`
  nTimes <- min(unlist(lapply(split(df, interaction(df[, group])), function(d) {
    length(unique(d[[x]]))
  }))) # for spline knots
  splineHelperForm <- NULL
  if (grepl("\\+", model)) {
    chngptHelperList <- .brmsChangePointHelper(model, x, y, group,
      dpar = FALSE, nTimes = nTimes,
      useGroup = USEGROUP, priors = priors, int = int
    )
    growthForm <- chngptHelperList$growthForm
    pars <- chngptHelperList$pars
    splineHelperForm <- chngptHelperList$splineHelperForm
  } else {
    matched_model <- gsub("homo", "int", matched_model)
    matched_model <- gsub("spline", "gam", matched_model)

    stringBrmsFormFun <- paste0(".brms_form_", gsub(" ", "", matched_model))
    brmsFormFun <- match.fun(stringBrmsFormFun)
    formRes <- brmsFormFun(x, y, group,
      dpar = FALSE, nTimes = nTimes,
      useGroup = USEGROUP, prior = priors, int = int
    )
    if (decay) {
      formRes <- .brms_form_decay(formRes, int)
    }
    pars <- formRes$pars
    growthForm <- formRes$form
  }

  #* `Make distributional parameter formulas`
  #* Note there is always a sigma after .sigmaHelper (it just might be homoskedastic)
  dpar_res <- lapply(seq_along(sigma), function(i) {
    dpar <- names(sigma)[[i]]
    model <- sigma[[i]]
    intModelRes <- .intModelHelper(model)
    model <- intModelRes$model
    sigmaInt <- intModelRes$int
    .brmDparHelper(dpar, model, x, group, nTimes, USEGROUP, priors, sigmaInt)
  })
  names(dpar_res) <- names(sigma)
  dparForm <- unlist(lapply(dpar_res, function(res) {
    res$dparForm
  }))
  dparSplineHelperForm <- unlist(lapply(dpar_res, function(res) {
    res$dparSplineHelperForm
  }))
  dpar_pars <- unlist(lapply(dpar_res, function(res) {
    res$dpar_pars
  }))
  names(dpar_pars) <- NULL
  pars <- append(pars, dpar_pars)
  pars <- pars[!grepl("spline", pars)]

  #* `Make hierarchical parameter model formulas`
  if (!is.null(hierarchical_predictor)) {
    if (is.null(hierarchy)) {
      warning(paste0(
        "hierarchy argument not provided, assuming linear models with intercepts ",
        "for all ", matched_model, " model parameters (",
        paste(pars, collapse = ", "),
        ")."
      ))
      hierarchy <- lapply(pars, function(p) {
        "int_linear"
      })
      names(hierarchy) <- pars
    }

    hrc_res <- lapply(names(hierarchy), function(pname) {
      hrc_model <- hierarchy[[pname]]
      intModelRes <- .intModelHelper(hrc_model)
      hrc_model <- intModelRes$model
      hrc_int <- intModelRes$int
      .brmDparHelper(
        dpar = pname, model = hrc_model, x = hierarchical_predictor,
        group, nTimes, USEGROUP, priors, int = hrc_int, force_nl = TRUE
      )
      #* here passing `pname` to the `dpar` argument of .brmDparHelper will make
      #* .brmDparHelper add that name as a prefix on all of the existing model parameters.
      #* Since all the parameter names are unique coming into this they will be unique coming
      #* out of this as well.
      #* I'm also passing hierarchical_predictor to the x argument
    })
    names(hrc_res) <- names(hierarchy)
    hrcForm <- unlist(lapply(hrc_res, function(res) {
      res$dparForm
    }))
    hrcSplineHelperForm <- unlist(lapply(hrc_res, function(res) {
      res$dparSplineHelperForm
    }))
    hrc_pars <- unlist(lapply(hrc_res, function(res) {
      res$dpar_pars
    }))
    names(hrc_pars) <- NULL
    pars <- append(pars, hrc_pars)
    pars <- pars[-which(pars %in% names(hierarchy))] #* remove pars that are now estimated by other
    #* sub models from the later steps.
  } else {
    hrcForm <- NULL
    hrcSplineHelperForm <- NULL
  }
  pars <- pars[!grepl("spline", pars)]

  #* `Make parameter grouping formulae`

  if (as.logical(length(pars))) {
    if (USEGROUP) {
      parForm <- as.formula(paste0(paste(pars, collapse = "+"), "~0+", paste(group, collapse = "*")))
    } else {
      parForm <- as.formula(paste0(paste(pars, collapse = "+"), "~1"))
    }
  } else {
    parForm <- NULL
  }

  #* `Combine formulas into brms.formula object`
  if (is.null(pars)) {
    NL <- FALSE
  } else {
    NL <- TRUE
  }
  bf_args <- list(
    "formula" = growthForm, dparForm, hrcForm, parForm,
    dparSplineHelperForm, hrcSplineHelperForm, splineHelperForm,
    "autocor" = corForm, "nl" = NL
  )
  bf_args <- bf_args[!unlist(lapply(bf_args, is.null))]
  bayesForm <- do.call(brms::bf, args = bf_args)

  out[["formula"]] <- bayesForm
  #* ***** `Make priors` *****
  out[["prior"]] <- .makePriors(priors, pars, df, group, USEGROUP, sigma, family, bayesForm)
  #* ***** `Make initializer function` *****
  if (as.logical(length(pars))) {
    initFun <- function(pars = "?", nPerChain = 1) {
      init <- lapply(pars, function(i) array(rgamma(nPerChain, 1)))
      names(init) <- paste0("b_", pars)
      init
    }
    formals(initFun)$pars <- pars
    formals(initFun)$nPerChain <- length(table(df[, group]))
    wrapper <- function() {
      initFun()
    }
  } else {
    wrapper <- 0
  }

  #* ***** `Raise Message for complex models` *****

  if (length(pars) * length(unique(interaction(df[, group]))) > 50) {
    message(paste0(
      "This model will estimate >50 parameters (excluding any smooth terms). \n\n",
      "If the MCMC is very slow then consider fitting separate models and using `combineDraws()` ",
      "to make a data.frame for hypothesis testing."
    ))
  }

  #* ***** `Return Components` *****
  out[["initfun"]] <- wrapper
  out[["df"]] <- df
  out[["family"]] <- family
  out[["pcvrForm"]] <- form
  return(out)
}
