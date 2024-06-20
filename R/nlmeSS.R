#' Ease of use nlme starter function for standard growth model parameterizations
#'
#' Internal to growthSS
#'
#' @param model One of the 8 model options in growthSS
#' @param form A pcvr style form, see growthSS
#' @param sigma One of "int", "power", or "exp", which will correspond to varIdent, varPower, or varExp
#' respectively.
#' This can also take a varFunc object.
#' @param df a dataframe to use to make the model.
#' @param pars optional parameters to vary by group as fixed effects.
#' @param start Starting values. These are optional unless model is a double sigmoid.
#' For any other model these will be estimated from the data if left NULL.
#'
#' @examples
#'
#' if (FALSE) {
#'   ex <- growthSim("logistic",
#'     n = 20, t = 25,
#'     params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#'   )
#'   ex$group <- factor(ex$group) # group MUST be a factor for nlme (typical)
#'   nl1 <- nlme::nlme(
#'     model = y ~ A / (1 + exp((B - time) / C)),
#'     data = ex,
#'     fixed = list(
#'       A ~ 0 + group,
#'       B ~ 0 + group,
#'       C ~ 0 + group
#'     ),
#'     random = A + B + C ~ 1,
#'     groups = ~group,
#'     weights = nlme::varIdent(form = ~ 1 | group),
#'     start = c(150, 150, 10, 10, 3, 3)
#'   )
#'   # linear sigma using power
#'   nl2 <- nlme::nlme(
#'     model = y ~ A / (1 + exp((B - time) / C)),
#'     data = ex,
#'     fixed = list(
#'       A ~ 0 + group,
#'       B ~ 0 + group,
#'       C ~ 0 + group
#'     ),
#'     random = A ~ 1,
#'     groups = ~group,
#'     weights = nlme::varPower(form = ~ time | group),
#'     start = c(150, 150, 10, 10, 3, 3)
#'   )
#'   # linear sigma using exp
#'   nl3 <- nlme::nlme(
#'     model = y ~ A / (1 + exp((B - time) / C)),
#'     data = ex,
#'     fixed = list(
#'       A ~ 0 + group,
#'       B ~ 0 + group,
#'       C ~ 0 + group
#'     ),
#'     random = A ~ 1,
#'     groups = ~group,
#'     weights = nlme::varExp(form = ~ time | group),
#'     start = c(150, 150, 10, 10, 3, 3)
#'   )
#' }
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#'
#' ss <- .nlmeSS(
#'   model = "logistic", form = y ~ time | id / group,
#'   sigma = "power", df = simdf, start = NULL
#' )
#'
#' ss <- .nlmeSS(
#'   model = "gam", form = y ~ time | id / group,
#'   sigma = "power", df = simdf, start = NULL
#' )
#'
#' dim(ss$df)
#' ss[c("formula", "taus", "start", "pcvrForm")]
#'
#' @import lmeSplines
#' @importFrom nlme varPower varIdent varExp nlme nlme.formula corAR1 pdIdent
#' @importFrom stats as.formula setNames
#' @importFrom methods is
#'
#' @keywords internal
#' @noRd


.nlmeSS <- function(model, form, sigma, df, pars = NULL, start = NULL, int = FALSE) {
  #* `general steps`
  #* ***** `Define choices and make empty output list`
  out <- list()
  models <- c(
    "logistic", "gompertz", "monomolecular",
    "exponential", "linear", "power law",
    "double logistic", "double gompertz", "gam",
    "frechet", "weibull", "gumbel", "logarithmic", "bragg", "lorentz", "beta"
  )
  sigmas <- c("none", "int", "power", "exp")
  #* check if sigma is class "varFunc", if it is then return it as is?

  #* ***** `Make nlme model formula` *****
  #* `parse form argument`
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  group <- parsed_form$group
  df <- parsed_form$data

  #* `make group a factor for nlme`
  df[[group]] <- as.factor(df[[group]])
  df[[paste0(group, "_numericLabel")]] <- unclass(df[[group]])
  #* `make an interaction variable for autocorrelation`
  #* Note that nlme does not allow random effects and correlations to apply at different "scales"
  #* so A,B,C can either vary by this interaction variable to have autocorrelation accurately modeled
  #* OR A,B,C can be estimated by group and autocorrelation can be by group. Currently this option is
  #* used. This note is kept here for reference.
  df[[paste0(group, individual)]] <- interaction(df[[group]], df[[individual]])

  #* `assemble growth formula with FE, RE, Groups, and Weights`
  if (grepl("decay", model)) {
    decay <- TRUE
    model <- trimws(gsub("decay", "", model))
  } else {
    decay <- FALSE
  }

  matched_model <- match.arg(model, models)

  if (is.character(sigma)) {
    matched_sigma <- match.arg(sigma, sigmas)
  } else {
    matched_sigma <- sigma
  }

  stringFormFun <- paste0(".nlme_form_", gsub(" ", "", matched_model))
  form_fun <- match.fun(stringFormFun)
  if (matched_model == "gam") { # gam takes some extra work
    df$splines <- lmeSplines::smspline(as.formula(paste0("~ ", x)), data = df) # spline setup
    start <- 0 # no starting values (no parameters)
    pars <- df # there are no pars, this is just to pass df to pdIdent for the splines
  }
  growthFormList <- form_fun(x, y, group, individual, matched_sigma, pars, int)
  pars <- growthFormList$pars
  growthFormList <- growthFormList[!grepl("pars", names(growthFormList))]
  if (decay) {
    growthFormList <- .nlmeDecay(growthFormList)
  }

  #* `Make starting values`
  if (matched_model == "gam") {
    start <- 0
  }
  if (is.null(start)) {
    if (matched_model == "double logistic") {
      warning(paste0(
        "Double Sigmoid models are not supported as self-starting models, ",
        "you will need to add starting parameters.",
        "Note for these models type='brms' is recommended."
      ))
      startingValues <- NULL
    } else if (matched_model == "double gompertz") {
      warning(paste0(
        "Double Sigmoid models are not supported as self-starting models, ",
        "you will need to add starting parameters.",
        " Note for these models type='brms' is recommended."
      ))
      startingValues <- NULL
    } else {
      stringInitFun <- paste0(".init", gsub(" ", "", matched_model))
      initFunction <- match.fun(stringInitFun)
      startingValues <- initFunction(df, x, y, int)
    }
    startingValuesList <- unlist(lapply(names(startingValues), function(nm) {
      val <- startingValues[nm]
      if (nm %in% pars) {
        rep(val, length(unique(df[[group]])))
        # if this is one of pars then make starting value per group
      } else {
        val
      } # else return one starting value
    }))
  } else {
    startingValuesList <- start
  }

  #* `return model components`

  out[["formula"]] <- growthFormList
  out[["start"]] <- startingValuesList
  out[["df"]] <- df
  out[["pcvrForm"]] <- form
  return(out)
}



#* `Sigma matching helper function`

.nlme_sigma_form <- function(matched_sigma, x, group) {
  if (group == "dummyGroup") {
    group <- NULL
  }
  #* `variance formula`
  if (methods::is(matched_sigma, "varFun")) {
    weights_form <- matched_sigma
  } else if (matched_sigma %in% c("int", "none")) {
    weights_form <- nlme::varIdent(
      form = stats::as.formula(paste(c("~1", group), collapse = "|"))
    )
  } else if (matched_sigma == "power") {
    weights_form <- nlme::varPower(
      form = stats::as.formula(paste(c(paste0("~", x), group), collapse = "|"))
    )
  } else if (matched_sigma == "exp") {
    weights_form <- nlme::varExp(
      form = stats::as.formula(paste(c(paste0("~", x), group), collapse = "|"))
    )
  }
  if (is.null(group)) {

  }
  return(weights_form)
}


#* `Define growth formulas`

.nlme_form_logistic <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ I + (A/(1+exp((B-", x, ")/C)))"))
  } else {
    total_pars <- c("A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ A/(1+exp((B-", x, ")/C))"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_gompertz <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ I + (A*exp(-B*exp(-C*", x, ")))"))
  } else {
    total_pars <- c("A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ A*exp(-B*exp(-C*", x, "))"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_doublelogistic <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C", "A2", "B2", "C2")
    model_form <- as.formula(paste0(
      y, " ~ I + (A/(1+exp((B-", x,
      ")/C)) + ((A2-A) /(1+exp((B2-", x,
      ")/C2))))"
    ))
  } else {
    total_pars <- c("A", "B", "C", "A2", "B2", "C2")
    model_form <- as.formula(paste0(
      y, " ~ A/(1+exp((B-", x,
      ")/C)) + ((A2-A) /(1+exp((B2-", x,
      ")/C2)))"
    ))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_doublegompertz <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C", "A2", "B2", "C2")
    model_form <- as.formula(paste0(
      y, " ~ I + (A * exp(-B * exp(-C*", x,
      ")) + (A2-A) * exp(-B2 * exp(-C2*(", x, "-B))))"
    ))
  } else {
    total_pars <- c("A", "B", "C", "A2", "B2", "C2")
    model_form <- as.formula(paste0(
      y, " ~ A * exp(-B * exp(-C*", x,
      ")) + (A2-A) * exp(-B2 * exp(-C2*(", x, "-B)))"
    ))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_monomolecular <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B")
    model_form <- as.formula(paste0(y, "~I + (A-A*exp(-B*", x, "))"))
  } else {
    total_pars <- c("A", "B")
    model_form <- as.formula(paste0(y, "~A-A*exp(-B*", x, ")"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_exponential <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B")
    model_form <- as.formula(paste0(y, " ~ I + (A*exp(B*", x, "))"))
  } else {
    total_pars <- c("A", "B")
    model_form <- as.formula(paste0(y, " ~ A*exp(B*", x, ")"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_linear <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A")
    model_form <- as.formula(paste0(y, " ~ I + A*", x))
  } else {
    total_pars <- c("A")
    model_form <- as.formula(paste0(y, " ~ A*", x))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_logarithmic <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A")
    model_form <- as.formula(paste0(y, " ~ I + A*log(", x, ")"))
  } else {
    total_pars <- c("A")
    model_form <- as.formula(paste0(y, " ~ A*log(", x, ")"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}


.nlme_form_powerlaw <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B")
    model_form <- as.formula(paste0(y, " ~ I + (A*", x, "^B)"))
  } else {
    total_pars <- c("A", "B")
    model_form <- as.formula(paste0(y, " ~ A*", x, "^B"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}



.nlme_form_gam <- function(x, y, group, individual, matched_sigma, pars, int) {
  model_form <- as.formula(paste0(y, " ~", x, "*", group))
  #* `random effects formula`
  random_form <- stats::setNames(list(nlme::pdIdent(~ splines - 1, data = pars)), group)
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  #* nlme insists that correlation and RE formulas use the same grouping,
  #* so i will not be able to account for individual autocorrelation in the GAM option
  correlation_form <- nlme::corAR1(0.8, form = stats::as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlmeDecay <- function(formList) {
  modelForm <- formList$model
  chars <- as.character(modelForm)
  formList$model <- as.formula(paste0(chars[2], chars[1], "-(", chars[3], ")"))
  formList
}



.nlme_form_weibull <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ I + (A * (1-exp(-(", x, "/C)^B)))"))
  } else {
    total_pars <- c("A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ A * (1-exp(-(", x, "/C)^B))"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_frechet <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ I + (A * exp(-((", x, "-0)/C)^(-B)))"))
  } else {
    total_pars <- c("A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ A * exp(-((", x, "-0)/C)^(-B))"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}


.nlme_form_gumbel <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ I + (A * exp(-exp(-(", x, "-B)/C)))"))
  } else {
    total_pars <- c("A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ A * exp(-exp(-(", x, "-B)/C))"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_bragg <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ I + A * exp(-B * (", x, " - C)^2)"))
  } else {
    total_pars <- c("A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ A * exp(-B * (", x, " - C)^2)"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_lorentz <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ I + A / (1 + B * (", x, " - C) ^ 2)"))
  } else {
    total_pars <- c("A", "B", "C")
    model_form <- as.formula(paste0(y, " ~ A / (1 + B * (", x, " - C) ^ 2)"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}

.nlme_form_beta <- function(x, y, group, individual, matched_sigma, pars, int) {
  #* `Define parameters and main growth formula`
  if (int) {
    total_pars <- c("I", "A", "B", "C", "D", "E")
    model_form <- as.formula(paste0(y, " ~ I + A * (((x - D) / (C - D)) * ((E - ", x,
                                    ") / (E - C)) ^ ((E - C) / (C - D))) ^ B"))
  } else {
    total_pars <- c("A", "B", "C", "D", "E")
    model_form <- as.formula(paste0(y, " ~ A * (((x - D) / (C - D)) * ((E - ", x,
                                    ") / (E - C)) ^ ((E - C) / (C - D))) ^ B"))
  }
  #* `random effects formula`
  random_form <- as.formula(paste0(paste0(total_pars, collapse = " + "), "~ 1"))
  #* `fixed effects formula`
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (is.null(group) || group == "dummyGroup") {
    pars <- ""
  }
  fixed_form <- lapply(total_pars, function(par) {
    if (par %in% pars) {
      stats::as.formula(paste0(par, " ~ 0 + ", group))
    } else {
      stats::as.formula(paste0(par, " ~ 1"))
    }
  })
  #* `groups formula`
  groups_form <- stats::as.formula(paste0("~", group))
  #* `variance formula`
  weights_form <- .nlme_sigma_form(matched_sigma, x, group)
  #* `correlation formula`
  correlation_form <- nlme::corAR1(0.8, form = as.formula(paste0("~ 1 |", group)))

  formulas <- list(
    "model" = model_form, "random" = random_form,
    "fixed" = fixed_form, "groups" = groups_form,
    "weights" = weights_form, "cor_form" = correlation_form, "pars" = pars
  )
  return(formulas)
}
