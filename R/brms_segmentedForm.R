#' Helper function to put together formulae for brms changepoint growth models
#'
#' @param model A multi-part model passed from brmSS passed from \code{\link{growthSS}}
#' @param x The x variable from the pcvrForm argument in \code{\link{growthSS}}
#' @param y The y variable from the pcvrForm argument in \code{\link{growthSS}}
#' @param group The grouping variable from the pcvrForm argument in \code{\link{growthSS}}
#' @param dpar Logical, is this a distributional parameter formula (TRUE) or part of the main growth
#' formula (FALSE)?
#' @param nTimes a Number of times that are present in the data, only used for making splines have a
#' workable number of knots.
#' @param useGroup logical, should groups be used?
#' @param priors A list describing priors in the style of \code{\link{brmSS}}, \code{\link{growthSS}},
#' and \code{\link{growthSim}}. This is only used currently to identify fixed and estimated
#' changepoints. If a changepoint is called "changePointX" with X being its position in the formula
#' then it will be estimated as a parameter in the model, but if the changepoint is called
#' "fixedChangePointX" then it will be passed as a numeric in the growth model.
#'
#' @examples
#' df1 <- do.call(rbind, lapply(1:30, function(i) {
#'   chngpt <- rnorm(2, 10, 1.5)
#'   A <- growthSim("linear", n = 1, t = chngpt[1], params = list("A" = c(1)))
#'   B <- growthSim("linear", n = 1, t = chngpt[2], params = list("A" = c(0.9)))
#'   B$group <- "b"
#'   x <- rbind(A, B)
#'   x$id <- paste0("id_", i)
#'   x
#' }))
#' df2 <- growthSim("linear", n = 30, t = 20, params = list("A" = c(4.1, 5)))
#' df2 <- do.call(rbind, lapply(unique(paste0(df2$id, df2$group)), function(int) {
#'   df1sub <- df1[paste0(df1$id, df1$group) == int, ]
#'   df2sub <- df2[paste0(df2$id, df2$group) == int, ]
#'   y_end <- df1sub[df1sub$time == max(df1sub$time), "y"]
#'   df2sub$time <- df2sub$time + max(df1sub$time)
#'   df2sub$y <- y_end + df2sub$y
#'   df2sub
#' }))
#' df <- rbind(df1, df2)
#' ggplot(df, aes(time, y, group = interaction(group, id))) +
#'   geom_line(aes(color = group)) +
#'   theme_minimal()
#'
#' .brmsChangePointHelper(model = "linear + linear", x = "time", y = "y", group = "group")
#'
#' @keywords internal
#' @noRd

.brmsChangePointHelper <- function(model, x, y, group, dpar = FALSE, nTimes = 25, useGroup, priors) {
  component_models <- trimws(strsplit(model, "\\+")[[1]])
  models <- c(
    "logistic", "gompertz", "monomolecular", "exponential", "linear", "power law", "gam",
    "spline", "int", "homo", "weibull", "frechet", "gumbel"
  )
  if (is.null(priors)) {
    priors <- stats::setNames(
      lapply(seq_along(component_models), identity),
      paste0("changePoint", seq_along(component_models))
    )
  }

  if (dpar) {
    prefix <- y
  } else {
    prefix <- NULL
  }

  mainGrowthModelPriorStrings <- paste(
    paste0("^", gsub(
      " ", "",
      c(models, "changePoint", "fixedChangePoint")
    )),
    collapse = "|"
  )
  if (dpar) {
    priors <- priors[grepl(prefix, names(priors))]
  } else {
    priors <- priors[grepl(mainGrowthModelPriorStrings, names(priors))]
  }
  # else should be any prior whose name starts with a model name,
  # with changePoint or with fixedChangePoint.

  formulae <- lapply(seq_along(component_models), function(i) {
    iter_model <- component_models[i]

    if (grepl("decay", iter_model)) {
      decay <- TRUE
      iter_model <- trimws(gsub("decay", "", iter_model))
    } else {
      decay <- FALSE
    }

    matched_iter_model <- match.arg(iter_model, models)
    if (matched_iter_model == "homo") {
      matched_iter_model <- "int"
    } # recoding
    if (matched_iter_model == "spline") {
      matched_iter_model <- "gam"
    }
    chngptFormFun <- match.fun(paste0(".", gsub(" ", "", matched_iter_model), "ChngptForm"))
    iter <- chngptFormFun(x, i, dpar = prefix, priors)
    if (decay) {
      iter <- .decayChngptForm(iter)
    }
    return(iter)
  })

  growthForm <- paste0(y, " ~ ", formulae[[1]]$form, " * ", formulae[[1]]$cp)
  #* Make cpInt cumulative
  for (i in 2:length(formulae)) {
    cumulativeCpInt <- do.call(paste, list(lapply(1:i, function(o) {
      formulae[[o]]$cpInt
    }), collapse = " + "))
    formulae[[i]]$cpInt <- cumulativeCpInt
  }
  # assemble segments into complete formula
  for (i in 2:length(formulae)) {
    nextPhase <- paste0(
      "+ (", formulae[[(i - 1)]]$cpInt, " + ", formulae[[i]]$form, ") * ",
      formulae[[i]]$cp
    )
    growthForm <- paste0(growthForm, nextPhase)
  }
  growthForm <- stats::as.formula(growthForm)

  if (dpar) {
    growthForm <- brms::nlf(growthForm)
  }

  params <- unique(unlist(lapply(formulae, function(f) {
    f$params
  })))
  params <- params[-length(params)]

  splineSegments <- which(unlist(lapply(formulae, function(fml) {
    "splineVar" %in% names(fml)
  })))

  splineForm <- .handleSplineSegments(splineSegments, useGroup, group, nTimes, formulae, x)

  return(list("growthForm" = growthForm, "pars" = params, "splineHelperForm" = splineForm))
}

#* ****************************************
#* ***** `Handle Spline Segments` *****
#* ****************************************

#' spline formula helper function
#' @keywords internal
#' @noRd

.handleSplineSegments <- function(splineSegments, useGroup, group, nTimes, formulae, x) {
  if (length(splineSegments) > 0) {
    if (useGroup) {
      by <- paste0(", by = ", group)
    } else {
      by <- ","
    }
    if (nTimes < 11) {
      k <- paste0(", k = ", nTimes)
    } else {
      k <- NULL
    }
    splineVars <- c()
    for (seg in splineSegments) {
      splineVars <- c(splineVars, formulae[[seg]]$splineVar)
    }
    lhs <- paste0(splineVars, collapse = "+")
    rhs <- paste0("s(", x, by, k, ")")
    splineForm <- paste0(lhs, "~", rhs)
  } else {
    splineForm <- NULL
  }
  return(splineForm)
}


#* ****************************************
#* ***** `Linear Changepoint Phase` *****
#* ****************************************

#' Linear changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .linearChngptForm(x = "time", 1)
#' .linearChngptForm(x = "time", 2)
#' .linearChngptForm(x = "time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value.
#' @noRd

.linearChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value, even in sub model
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }

    form <- paste0(prefix, "linear", position, "A * ", x)
    cp <- paste0("inv_logit((", chngptPrefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0("(", prefix, "linear", position, "A * ", chngptPrefix, changePointObj, ")")
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prevAndCurrentChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts))
                                                    %in% c(position, position - 1))]
    #* per location where "fixed" is in the prior name, replace the name with that number.
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }
    pac_fixed_index <- which(grepl("fixed", prevAndCurrentChangePoints))
    if (length(pac_fixed_index) > 0) {
      prevAndCurrentChangePoints[pac_fixed_index] <- as.numeric(
        priors[[prevAndCurrentChangePoints[pac_fixed_index]]]
      )
    }

    form <- paste0(
      prefix, "linear", position, "A * (", x, "-",
      paste0(prevChangePoints, collapse = "-"), ")"
    )
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(prefix, "linear", position, "A * (", paste0(rev(prevAndCurrentChangePoints),
                                                                collapse = "-"),
                    ")")
    #* cpInt would be wrong for the last position but it isn't used.
  }

  pars <- paste0(prefix, "linear", position, "A")
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }

  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}


#* ****************************************
#* ***** `Logistic Changepoint Phase` *****
#* ****************************************

#' Logistic changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .logisticChngptForm(x = "time", 1)
#' .logisticChngptForm(x = "time", 2)
#' .logisticChngptForm(x = "time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value.
#'
#' @noRd

.logisticChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }

    form <- paste0(
      prefix, "logistic", position, "A / (1 + exp( (", prefix, "logistic", position,
      "B-(", x, "))/", prefix, "logistic", position, "C) )"
    )
    cp <- paste0("inv_logit((", chngptPrefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0(
      prefix, "logistic", position, "A / (1 + exp( (", prefix, "logistic", position,
      "B-(", chngptPrefix, changePointObj, "))/", prefix, "logistic", position, "C) )"
    )
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prevAndCurrentChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts))
                                                    %in% c(position, position - 1))]
    #* per location where "fixed" is in the prior name, replace the name with that number.
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }
    pac_fixed_index <- which(grepl("fixed", prevAndCurrentChangePoints))
    if (length(pac_fixed_index) > 0) {
      prevAndCurrentChangePoints[pac_fixed_index] <- as.numeric(
        priors[[prevAndCurrentChangePoints[pac_fixed_index]]]
      )
    }

    form <- paste0(
      prefix, "logistic", position, "A / (1 + exp( (", prefix, "logistic", position,
      "B-(", x, "-", paste0(prevChangePoints, collapse = "-"),
      "))/", prefix, "logistic", position, "C) )"
    )
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(
      prefix, "logistic", position, "A / (1 + exp( (", prefix, "logistic", position,
      "B-(", paste0(rev(prevAndCurrentChangePoints), collapse = "-"), "))/", prefix, "logistic",
      position, "C) )"
    )
  }
  pars <- paste0(prefix, "logistic", position, c("A", "B", "C"))
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }
  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}


#* ****************************************
#* ***** `Gompertz Changepoint Phase` *****
#* ****************************************

#' Gompertz changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .gompertzChngptForm(x = "time", 1)
#' .gompertzChngptForm(x = "time", 2)
#' .gompertzChngptForm(x = "time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value.
#' @noRd

.gompertzChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  if (dpar) {
    prefix <- chngptPrefix <- "sub"
  } else {
    prefix <- chngptPrefix <- NULL
  }

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }
    form <- paste0(
      prefix, "gompertz", position, "A * exp(-", prefix, "gompertz", position,
      "B * exp(-", prefix, "gompertz", position, "C * ", x, "))"
    )
    cp <- paste0("inv_logit((", prefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0(
      prefix, "gompertz", position, "A * exp(-", prefix, "gompertz", position,
      "B * exp(-", prefix, "gompertz", position, "C * ", chngptPrefix, changePointObj,
      "))"
    )
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prevAndCurrentChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts))
                                                    %in% c(position, position - 1))]
    #* per location where "fixed" is in the prior name, replace the name with that number.
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }
    pac_fixed_index <- which(grepl("fixed", prevAndCurrentChangePoints))
    if (length(pac_fixed_index) > 0) {
      prevAndCurrentChangePoints[pac_fixed_index] <- as.numeric(
        priors[[prevAndCurrentChangePoints[pac_fixed_index]]]
      )
    }

    form <- paste0(
      prefix, "gompertz", position, "A * exp(-", prefix, "gompertz", position,
      "B * exp(-", prefix, "gompertz", position, "C * (", x, " - ",
      paste0(prevChangePoints, collapse = "-"), ")))"
    )

    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(
      prefix, "gompertz", position, "A * exp(-", prefix, "gompertz", position, "B * exp(-", prefix,
      "gompertz", position,
      "C * (", paste0(rev(prevAndCurrentChangePoints), collapse = "-"), "))"
    )
  }
  pars <- paste0(prefix, "gompertz", position, c("A", "B", "C"))
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }
  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}

#* ****************************************
#* ***** `monomolecular Changepoint Phase` *****
#* ****************************************

#' Monomolecular changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .monomolecularChngptForm(x = "time", 1)
#' .monomolecularChngptForm(x = "time", 2)
#' .monomolecularChngptForm(x = "time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value.
#'
#' @noRd

.monomolecularChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }
    form <- paste0(
      prefix, "monomolecular", position, "A-", prefix, "monomolecular", position,
      "A * exp(-", prefix, "monomolecular", position, "B * ", x, ")"
    )
    cp <- paste0("inv_logit((", chngptPrefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0(
      prefix, "monomolecular", position, "A-", prefix, "monomolecular", position,
      "A * exp(-", prefix, "monomolecular", position, "B * ", chngptPrefix,
      changePointObj, ")"
    )
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prevAndCurrentChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts))
                                                    %in% c(position, position - 1))]
    #* per location where "fixed" is in the prior name, replace the name with that number.
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }
    pac_fixed_index <- which(grepl("fixed", prevAndCurrentChangePoints))
    if (length(pac_fixed_index) > 0) {
      prevAndCurrentChangePoints[pac_fixed_index] <- as.numeric(
        priors[[prevAndCurrentChangePoints[pac_fixed_index]]]
      )
    }

    form <- paste0(
      prefix, "monomolecular", position, "A-", prefix, "monomolecular", position, "A * exp(-",
      prefix, "monomolecular", position, "B * ", x, "-",
      paste0(prevChangePoints, collapse = "-"), ")"
    )
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(
      prefix, "monomolecular", position, "A-", prefix, "monomolecular", position, "A * exp(-",
      prefix, "monomolecular", position, "B * ",
      paste0(rev(prevAndCurrentChangePoints), collapse = "-"), ")"
    )
  }
  pars <- paste0(prefix, "monomolecular", position, c("A", "B"))
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }
  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}


#* ****************************************
#* ***** `Exponential Changepoint Phase` *****
#* ****************************************

#' Exponential changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .exponentialChngptForm(x = "time", 1)
#' .exponentialChngptForm(x = "time", 2)
#' .exponentialChngptForm(x = "time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value.
#' @noRd

.exponentialChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }

    form <- paste0(
      prefix, "exponential", position, "A * exp(", prefix, "exponential", position,
      "B * ", x, ")"
    )
    cp <- paste0("inv_logit((", chngptPrefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0(
      prefix, "exponential", position, "A * exp(", prefix, "exponential", position,
      "B * ", chngptPrefix, changePointObj, ")"
    )
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prevAndCurrentChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts))
                                                    %in% c(position, position - 1))]
    #* per location where "fixed" is in the prior name, replace the name with that number.
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }
    pac_fixed_index <- which(grepl("fixed", prevAndCurrentChangePoints))
    if (length(pac_fixed_index) > 0) {
      prevAndCurrentChangePoints[pac_fixed_index] <- as.numeric(
        priors[[prevAndCurrentChangePoints[pac_fixed_index]]]
      )
    }

    form <- paste0(
      prefix, "exponential", position, "A * exp(", prefix, "exponential", position, "B * (",
      x, "-", paste0(prevChangePoints, collapse = "-"), "))"
    )
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(
      prefix, "exponential", position, "A * exp(", prefix, "exponential", position, "B * (",
      paste0(rev(prevAndCurrentChangePoints), collapse = "-"), "))"
    )
  }
  pars <- paste0(prefix, "exponential", position, c("A", "B")) # this needs to be conditional on fixed
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }
  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}

#* ****************************************
#* ***** `Power Law Changepoint Phase` *****
#* ****************************************

#' Power Law changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .powerlawChngptForm(x = "time", 1)
#' .powerlawChngptForm(x = "time", 2)
#' .powerlawChngptForm(x = "time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value.
#'
#' @noRd

.powerlawChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }
    form <- paste0(prefix, "powerLaw", position, "A * ", x, "^(", prefix, "powerLaw", position, "B)")
    cp <- paste0("inv_logit((", chngptPrefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0(
      prefix, "powerLaw", position, "A * ", chngptPrefix, changePointObj, "^(", prefix,
      "powerLaw", position, "B)"
    )
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prevAndCurrentChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts))
                                                    %in% c(position, position - 1))]
    #* per location where "fixed" is in the prior name, replace the name with that number.
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }
    pac_fixed_index <- which(grepl("fixed", prevAndCurrentChangePoints))
    if (length(pac_fixed_index) > 0) {
      prevAndCurrentChangePoints[pac_fixed_index] <- as.numeric(
        priors[[prevAndCurrentChangePoints[pac_fixed_index]]]
      )
    }

    form <- paste0(
      prefix, "powerLaw", position, "A * ", x, "-",
      paste0(prevChangePoints, collapse = "-"), "^(", prefix, "powerLaw", position, "B)"
    )
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(
      prefix, "powerLaw", position, "A * (", paste0(rev(prevAndCurrentChangePoints), collapse = "-"),
      ")^(", prefix, "powerLaw", position, "B)"
    )
  }
  pars <- paste0(prefix, "powerLaw", position, c("A", "B"))
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }
  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}

#* ****************************************
#* ***** `Intercept Changepoint Phase` *****
#* ****************************************

#' intercept only changepoint section function
#'
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .intChngptForm(x = "time", 1, nTimes = 20)
#' .intChngptForm(x = "time", 2, nTimes = 20)
#' .intChngptForm(x = "time", 3, nTimes = 5)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value, for GAMs this is
#' undefined and GAMs should only be used at the end of a segmented model.
#' @noRd

.intChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }
    form <- paste0(prefix, "int", position)
    cp <- paste0("inv_logit((", chngptPrefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0(prefix, "int", position)
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }

    form <- paste0(prefix, "int", position)
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(prefix, "int", position)
  }

  pars <- paste0(prefix, "int", position)
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }

  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}




#* ****************************************
#* ***** `Gam Changepoint Phase` *****
#* ****************************************

#' gam changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .gamChngptForm(x = "time", 1, nTimes = 20)
#' .gamChngptForm(x = "time", 2, nTimes = 20)
#' .gamChngptForm(x = "time", 3, nTimes = 5)
#'
#' @return a list with form, cp, cpInt, and splineForm elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value, for GAMs this is
#' undefined and GAMs should only be used at the end of a segmented model.
#' "splineForm" is to use in making a spline for a predictor.
#' @noRd

.gamChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    stop("GAMs are only supported as the last function of a multi-part formula")
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }

    form <- paste0(prefix, "spline")
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- NA
  }
  pars <- paste0(prefix, "spline")
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }
  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars,
    "splineVar" = paste0(prefix, "spline")
  ))
}

#* ****************************************
#* ***** `Gam Changepoint Phase` *****
#* ****************************************

#' flip any model to a decay model
#'
#' @param phaseList A list returned from some *ChngptForm function
#'
#' @return a list with form, cp, cpInt and params for a decay segment to a model
#' @noRd

.decayChngptForm <- function(phaseList) {
  phaseList$form <- paste0("-", phaseList$form)
  phaseList
}

#* ****************************************
#* ***** `Weibull Changepoint Phase` *****
#* ****************************************

#' Weibull changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .weibullChngptForm(x = "time", 1)
#' .weibullChngptForm(x = "time", 2)
#' .weibullChngptForm(x = "time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value.
#'
#' @noRd

.weibullChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }
    form <- paste0(
      prefix, "weibull", position, "A * (1-exp(-(", x, "/", prefix,
      "weibull", position, "C)^", prefix, "weibull", position, "B))"
    )
    cp <- paste0("inv_logit((", chngptPrefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0(
      prefix, "weibull", position, "A * (1-exp(-(", chngptPrefix, changePointObj,
      "/", prefix, "weibull", position, "C)^", prefix, "weibull", position, "B))"
    )
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prevAndCurrentChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts))
                                                    %in% c(position, position - 1))]
    #* per location where "fixed" is in the prior name, replace the name with that number.
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }
    pac_fixed_index <- which(grepl("fixed", prevAndCurrentChangePoints))
    if (length(pac_fixed_index) > 0) {
      prevAndCurrentChangePoints[pac_fixed_index] <- as.numeric(
        priors[[prevAndCurrentChangePoints[pac_fixed_index]]]
      )
    }
    form <- paste0(
      prefix, "weibull", position, "A * (1-exp(-(", x, "-", paste0(prevChangePoints, collapse = "-"),
      "/", prefix, "weibull", position, "C)^", prefix, "weibull", position, "B))"
    )
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(
      prefix, "weibull", position, "A * (1-exp(-(", paste0(rev(prevAndCurrentChangePoints),
        collapse = "-"
      ),
      "/", prefix, "weibull", position, "C)^", prefix, "weibull", position, "B))"
    )
  }
  pars <- paste0(prefix, "weibull", position, c("A", "B", "C"))
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }
  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}


#* ****************************************
#* ***** `Frechet Changepoint Phase` *****
#* ****************************************

#' Frechet changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .frechetChngptForm(x = "time", 1)
#' .frechetChngptForm(x = "time", 2)
#' .frechetChngptForm(x = "time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value.
#'
#' @noRd

.frechetChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt

  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }
    form <- paste0(
      prefix, "frechet", position, "A * exp(-((", x, "-0)/", prefix, "frechet",
      position, "C)^(-", prefix, "frechet", position, "B))"
    )
    cp <- paste0("inv_logit((", chngptPrefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0(
      prefix, "frechet", position, "A * exp(-((", chngptPrefix, changePointObj, "-0)/",
      prefix, "frechet", position, "C)^(-", prefix, "frechet", position, "B))"
    )
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prevAndCurrentChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts))
                                                    %in% c(position, position - 1))]
    #* per location where "fixed" is in the prior name, replace the name with that number.
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }
    pac_fixed_index <- which(grepl("fixed", prevAndCurrentChangePoints))
    if (length(pac_fixed_index) > 0) {
      prevAndCurrentChangePoints[pac_fixed_index] <- as.numeric(
        priors[[prevAndCurrentChangePoints[pac_fixed_index]]]
      )
    }

    form <- paste0(
      prefix, "frechet", position, "A * exp(-((", x, "-", paste0(prevChangePoints, collapse = "-"),
      "-0)/", prefix, "frechet", position, "C)^(-", prefix, "frechet", position, "B))"
    )
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(
      prefix, "frechet", position, "A * exp(-((", paste0(rev(prevAndCurrentChangePoints),
        collapse = "-"
      ),
      "-0)/", prefix, "frechet", position, "C)^(-", prefix, "frechet", position, "B))"
    )
  }
  pars <- paste0(prefix, "frechet", position, c("A", "B", "C"))
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }
  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}


#* ****************************************
#* ***** `Gumbel Changepoint Phase` *****
#* ****************************************

#' Gumbel changepoint section function
#'
#' @param x X variable name
#' @param position Position in growth formula ("1" + "2" + "3"... etc)
#' @param dpar string or NULL, if string should be the name of the distributional parameter
#' @param priors a list of prior distributions (used for fixed vs estimated changepoints)
#'
#' @examples
#'
#' .gumbelChngptForm(x = "time", 1)
#' .gumbelChngptForm(x = "time", 2)
#' .gumbelChngptForm(x = "time", 3)
#'
#' @return a list with form, cp, and cpInt elements. "form" is the growth formula
#' for this phase of the model. "cp" is the inv_logit function defining when this
#' phase should happen. "cpInt" is the value at the end of this growth phase and is
#' used in starting the next growth phase from the right y value.
#'
#' @noRd

.gumbelChngptForm <- function(x, position = 1, dpar = NULL, priors) { # return f, cp, and cpInt


  prefix <- chngptPrefix <- dpar

  if (any(grepl(paste0("fixedChangePoint", position), names(priors)))) {
    changePointObj <- as.numeric(priors[[paste0(prefix, "fixedChangePoint", position)]])[1]
    fixed <- TRUE
    chngptPrefix <- NULL # never a prefix if the changepoint is a fixed value
  } else {
    fixed <- FALSE
  }

  if (position == 1) {
    if (!fixed) {
      changePointObj <- "changePoint1"
    }
    form <- paste0(
      prefix, "gumbel", position, "A * exp(-exp(-(", x, "-",
      prefix, "gumbel", position, "B)/", prefix, "gumbel", position, "C))"
    )
    cp <- paste0("inv_logit((", chngptPrefix, changePointObj, " - ", x, ") * 5)")
    cpInt <- paste0(
      prefix, "gumbel", position, "A * exp(-exp(-(", chngptPrefix, changePointObj, "-",
      prefix, "gumbel", position, "B)/", prefix, "gumbel", position, "C))"
    )
  } else {
    all_chngpts <- names(priors)[grepl("fixedChangePoint*|changePoint*", names(priors))]
    prevChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts)) < position)]
    prevAndCurrentChangePoints <- all_chngpts[which(as.numeric(sub(".*hangePoint", "", all_chngpts))
                                                    %in% c(position, position - 1))]
    #* per location where "fixed" is in the prior name, replace the name with that number.
    prev_fixed_index <- which(grepl("fixed", prevChangePoints))
    if (length(prev_fixed_index) > 0) {
      prevChangePoints[prev_fixed_index] <- as.numeric(priors[[prevChangePoints[prev_fixed_index]]])
    }
    pac_fixed_index <- which(grepl("fixed", prevAndCurrentChangePoints))
    if (length(pac_fixed_index) > 0) {
      prevAndCurrentChangePoints[pac_fixed_index] <- as.numeric(
        priors[[prevAndCurrentChangePoints[pac_fixed_index]]]
      )
    }

    form <- paste0(
      prefix, "gumbel", position, "A * exp(-exp(-(", x, "-", paste0(prevChangePoints, collapse = "-"),
      "-", prefix, "gumbel", position, "B)/", prefix, "gumbel", position, "C))"
    )
    cp <- paste0("inv_logit((", x, "-", paste0(prevChangePoints, collapse = "-"), ") * 5)")
    cpInt <- paste0(
      prefix, "gumbel", position, "A * exp(-exp(-(", paste0(rev(prevAndCurrentChangePoints),
        collapse = "-"
      ),
      "-", prefix, "gumbel", position, "B)/", prefix, "gumbel", position, "C))"
    )
  }
  pars <- paste0(prefix, "gumbel", position, c("A", "B", "C"))
  if (!fixed) {
    pars <- c(pars, paste0(chngptPrefix, "changePoint", position))
  }
  return(list(
    "form" = form,
    "cp" = cp,
    "cpInt" = cpInt,
    "params" = pars
  ))
}
