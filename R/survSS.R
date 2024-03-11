#' Function to parse survival model specifications in growthSS for brms type models
#'
#' @param model a survival distribution to use (currently "binomial" and "weibull" are supported)
#' @param form a formula in pcvr syntax
#' @param df a dataframe to use
#' @param priors priors specified per details in growthSS
#'
#' @return A list of elements to pass to fitGrowth
#'
#' @examples
#' set.seed(123)
#' model = "survival weibull"
#' form <- y > 100 ~ time | id / group
#' df <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' surv <- .survModelParser(model)
#' survivalBool <- surv$survival
#' model <- surv$model
#' prior <- c(0, 5)
#' ss <- .brmsSurvSS(model, form, df, prior)
#' lapply(ss, head)
#'
#' @keywords internal
#' @noRd

.brmsSurvSS <- function(model = NULL, form = NULL, df = NULL, priors = NULL) {
  out <- list()
  #* `make survival data`
  fixData <- TRUE
  if (grepl("binomial", model)) {
    if (all(c("n_events", "n_eligible") %in% colnames(df))) {
      fixData <- FALSE
    }
  } else { # weibull
    if (all(c("event", "censor") %in% colnames(df))) {
      fixData <- FALSE
    }
  }
  if (fixData) {
    makeSurvDataRet <- .makeSurvData(df, form, model)
    out_df <- makeSurvDataRet$data
    out[["df"]] <- out_df
  } else {
    makeSurvDataRet <- list()
    getGroup <- trimws(strsplit(as.character(form)[3], "[|]|[/]")[[1]])
    makeSurvDataRet$group <- getGroup[length(getGroup)]
    makeSurvDataRet$x <- getGroup[1]
    out[["df"]] <- df
  }

  #* `make bayesian formula`
  if (model == "binomial") {
    form_ret <- .brmsurv_binomial_formula(x = makeSurvDataRet$x, group = makeSurvDataRet$group,
                                          df = out_df)
  } else if (model == "weibull") {
    form_ret <- .brmsurv_weibull_formula(x = makeSurvDataRet$x, group = makeSurvDataRet$group,
                                         df = out_df)
  }
  out[["family"]] <- form_ret$family
  out[["formula"]] <- form_ret$formula
  #* `make priors if none specified`
  out[["priors"]] <- .brmsMakeSurvPriors(priors, out_df, makeSurvDataRet)
  #* `set initialization to 0 for all chains`
  out[["initfun"]] <- 0
  out[["pcvrForm"]] <- form
  return(out)
}


#' Helper function to make priors for brms survival models
#' @return A list with a formula and a model family
#' @keywords internal
#' @noRd

.brmsMakeSurvPriors <- function(priors, out_df, makeSurvDataRet, form_ret) {
  if (is.null(priors)) {
    return(brms::prior_string("normal(0, 5)", class = "b"))
  } else if (any(methods::is(priors, "brmsprior"))) {
    return(priors)
  } else if (is.numeric(priors)) {
    priors <- rep(priors, length.out = 2 * length(unique(out_df[[makeSurvDataRet$group]])))
    priors <- stats::setNames(
      lapply(seq(1, length(priors), 2), function(i) {
        c(priors[i], priors[i + 1])
      }),
      unique(out_df[[makeSurvDataRet$group]])
    )
    message(
      "Prior is numeric, replicating to ", length(unique(out_df[[makeSurvDataRet$group]])),
      " length 2 elements (mu, sd) and assuming order ",
      paste(unique(out_df[[makeSurvDataRet$group]]),
            collapse = ", ")
    )
  }
  pars <- brms::get_prior(formula = form_ret$formula, data = out_df, family = form_ret$family)$coef
  pars <- pars[which(nchar(pars) > 0)]
  if (length(priors) != length(pars)) {
    message(paste0(
      "Priors and parameters are not the same length. Output will assume that priors are for groups",
      " and are in order: ", paste(unique(out_df[[makeSurvDataRet$group]]), collapse = ", ")
    ))
    priors <- stats::setNames(
      rep(priors, length.out = length(unique(out_df[[makeSurvDataRet$group]]))),
      unique(out_df[[makeSurvDataRet$group]])
    )
  }
  prior_obj <- NULL
  for (g in unique(out_df[[makeSurvDataRet$group]])) {
    sub_pars <- pars[grepl(paste0(makeSurvDataRet$group, g), pars)]
    for (param in sub_pars) {
      if (is.null(prior_obj)) {
        prior_obj <- brms::set_prior(
          prior = paste0("normal(", priors[[g]][1], ",", priors[[g]][2], ")"),
          coef = param
        )
      } else {
        prior_obj <- prior_obj + brms::set_prior(
          prior = paste0("normal(", priors[[g]][1], ",", priors[[g]][2], ")"),
          coef = param
        )
      }
    }
  }
  return(prior_obj)
}


#' Helper function to make formulas for brms binomial survival models
#' @return A list with a formula and a model family
#' @keywords internal
#' @noRd

.brmsurv_binomial_formula <- function(y = "n_events", x = "time", total = "n_eligible",
                                      group = "groups", df = NULL) {
  #* make formula
  if (is.null(group) || length(unique(df[[group]])) == 1) {
    formula <- stats::as.formula(paste0(y, " | trials(", total, ") ~ 0 + ", x))
  } else {
    formula <- stats::as.formula(paste0(y, " | trials(", total, ") ~ 0 + ", x, ":", group))
  }
  #* specify family
  family <- "binomial" # using default links
  return(list("formula" = formula, "family" = family))
}

#' Helper function to make formulas for brms weibull survival models
#' @return A list with a formula and a model family
#' @keywords internal
#' @noRd

.brmsurv_weibull_formula <- function(y = "event", x = "time", censor = "censor",
                                     group = "groups", df = NULL) {
  #* make formula
  if (is.null(group) || length(unique(df[[group]])) == 1) {
    formula <- stats::as.formula(paste0(x, " | cens(", censor, ") ~ 1"))
  } else {
    formula <- stats::as.formula(paste0(x, " | cens(", censor, ") ~ 0 + ", group))
  }
  #* specify family
  family <- "weibull" # using default links
  return(list("formula" = formula, "family" = family))
}




#' Function to parse survival model specifications in growthSS for modeling in the survival package
#'
#' @return a list of components passed to fitGrowth
#'
#' @examples
#'
#' set.seed(123)
#' model = "survival weibull"
#' form <- y > 100 ~ time | id / group
#' df <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' surv <- .survModelParser(model)
#' survivalBool <- surv$survival
#' model <- surv$model
#' ss <- .survSS(model, form, df)
#' lapply(ss, head)
#'
#' @importFrom survival survreg Surv
#'
#' @keywords internal
#' @noRd

.survSS <- function(model = NULL, form = NULL, df = NULL) {
  out <- list()
  #* `make survival data`
  fixData <- TRUE
  if (all(c("event") %in% colnames(df))) {
    fixData <- FALSE
  }
  if (fixData) {
    makeSurvDataRet <- .makeSurvData(df, form, model = "weibull")
    out_df <- makeSurvDataRet$data
    out_df[[makeSurvDataRet$group]] <- factor(out_df[[makeSurvDataRet$group]])
    out_df[[paste0(makeSurvDataRet$group, "_numericLabel")]] <- unclass(out_df[[makeSurvDataRet$group]])
    out[["df"]] <- out_df
  } else {
    makeSurvDataRet <- list()
    getGroup <- trimws(strsplit(as.character(form)[3], "[|]|[/]")[[1]])
    makeSurvDataRet$group <- getGroup[length(getGroup)]
    makeSurvDataRet$x <- getGroup[1]
    out[["df"]] <- df
  }
  #* `make survreg formula`
  x <- makeSurvDataRet$x
  group <- makeSurvDataRet$group
  if (is.null(group) || length(unique(df[[group]])) == 1) {
    formula <- stats::as.formula(paste0("Surv(", x, ", event) ~ 1"))
  } else {
    formula <- stats::as.formula(paste0("Surv(", x, ", event) ~ 1 + group"))
  }
  out[["formula"]] <- formula
  #* `return all`
  out[["pcvrForm"]] <- form
  out[["distribution"]] <- model
  return(out)
}


#' Function to parse survival model specifications in growthSS for modeling in the flexsurv package
#'
#' @return a list of components passed to fitGrowth
#'
#' @examples
#'
#' set.seed(123)
#' model = "survival gengamma"
#' form <- y > 100 ~ time | id / group
#' df <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' surv <- .survModelParser(model)
#' survivalBool <- surv$survival
#' model <- surv$model
#' ss <- .flexSurvSS(model, form, df)
#' lapply(ss, head)
#'
#' @importFrom survival survreg Surv
#'
#' @keywords internal
#' @noRd

.flexSurvSS <- function(model = NULL, form = NULL, df = NULL, anc = NULL) {
  out <- list()
  distributions <- c(
    "gengamma", "gengamma.orig", "genf", "genf.orig",
    "weibull", "gamma", "exp", "llogis", "lnorm", "gompertz",
    "exponential", "lognormal"
  )
  if (!model %in% distributions) {
    stop(paste0(
      "Supported distributions for flexsurv models are ",
      paste(distributions, collapse = ", "),
      ".\nIf you are using a custom distribution please call flexsurvreg directly."
    ))
  }
  #* `make survival data`
  fixData <- TRUE
  if (all(c("event") %in% colnames(df))) {
    fixData <- FALSE
  }
  if (fixData) {
    makeSurvDataRet <- .makeSurvData(df, form, model = "weibull")
    out_df <- makeSurvDataRet$data
    out_df[[makeSurvDataRet$group]] <- factor(out_df[[makeSurvDataRet$group]])
    out_df[[paste0(makeSurvDataRet$group, "_numericLabel")]] <- unclass(out_df[[makeSurvDataRet$group]])
    out[["df"]] <- out_df
  } else {
    makeSurvDataRet <- list()
    getGroup <- trimws(strsplit(as.character(form)[3], "[|]|[/]")[[1]])
    makeSurvDataRet$group <- getGroup[length(getGroup)]
    makeSurvDataRet$x <- getGroup[1]
    out[["df"]] <- df
  }
  #* `make main survival formula`
  x <- makeSurvDataRet$x
  group <- makeSurvDataRet$group
  if (is.null(group) || length(unique(df[[group]])) == 1) {
    formula <- stats::as.formula(paste0("Surv(", x, ", event) ~ 1"))
  } else {
    formula <- stats::as.formula(paste0("Surv(", x, ", event) ~ 1 + group"))
  }
  out[["formula"]][["f1"]] <- formula
  #* `collect additional formulas`
  if (!is.null(anc)) {
    out[["formula"]][["f2"]] <- anc
  } else {
    out[["formula"]][["f2"]] <- NULL
  }

  #* `return all`
  out[["pcvrForm"]] <- form
  out[["distribution"]] <- model
  return(out)
}
