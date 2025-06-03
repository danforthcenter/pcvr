#' Ease of use nlrq/nls starter function for standard growth model parameterizations
#'
#' Internal to growthSS
#'
#' @examples
#'
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#'
#' ss <- .nlrqSS(
#'   model = "logistic", form = y ~ time | id / group,
#'   tau = 0.5, df = simdf, start = NULL
#' )
#' model <- "logistic"
#' form <- y ~ time | id / group
#' tau <- 0.5
#' df <- simdf
#' start <- NULL
#' pars <- NULL
#' ss <- .nlrqSS(model = "gam", form = y ~ time | id / group, df = simdf, start = NULL, tau = 0.5)
#'
#' dim(ss$df)
#' ss[c("formula", "taus", "start", "pcvrForm")]
#' @importFrom splines bs
#' @importFrom stats sortedXyData
#' @keywords internal
#' @noRd

.nlrqSS <- function(model, form, tau = 0.5, df, pars = NULL, start = NULL, type = "nlrq", int = FALSE) {
  #* ***** `Define choices and make empty output list`
  out <- list()
  models <- .available_models()
  #* ***** `Make nlrq formula` *****
  #* `parse form argument`
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  x <- parsed_form$x
  group <- parsed_form$group
  USEGROUP <- parsed_form$USEG
  if (parsed_form$USEID) {
    message(paste0("Individual is not used with type = '", type, "'."))
  }
  df <- parsed_form$data
  if (model == "gam") {
    df[[paste(group, collapse = ".")]] <- interaction(df[, group])
    group <- paste(group, collapse = ".")
  }
  #* `if there are multiple groups combine them for nlrq/nls models`
  #* Also assign numeric labels to factor groups
  if (USEGROUP) {
    new_group <- paste0(group, collapse = ".")
    df[[new_group]] <- interaction(df[, group])
    group <- new_group
    df[[paste0(group, "_numericLabel")]] <- unclass(df[[group]])
  }
  #* `assemble growth formula`
  if (grepl("decay", model)) {
    decay <- TRUE
    model <- trimws(gsub("decay", "", model))
  } else {
    decay <- FALSE
  }

  matched_model <- match.arg(model, models)
  stringFormFun <- paste0(".nlrq_form_", gsub(" ", "", matched_model))
  form_fun <- match.fun(stringFormFun)
  res <- form_fun(x, y, USEGROUP, group, pars, int)
  growthForm <- res[[1]]
  pars <- res[[2]]

  if (decay) {
    growthForm <- .nlrqDecay(growthForm)
  }

  if (matched_model == "gam") {
    start <- 0
  }

  if (is.null(start)) {
    if (grepl("double", matched_model)) {
      warning(paste0(
        "Double Sigmoid models are not supported as self-starting models,",
        " you will need to add starting parameters.",
        " Note for these models type='brms' is recommended."
      ))
      startingValues <- NULL
    } else {
      stringInitFun <- paste0(".init", gsub(" ", "", matched_model))
      initFunction <- match.fun(stringInitFun)
      startingValues <- initFunction(df, x, y, int)
    }
    if ((!matched_model %in% c("double logistic", "double gompertz")) && USEGROUP) {
      nms <- names(startingValues)
      startingValuesList <- lapply(names(startingValues), function(nm) {
        val <- startingValues[nm]
        if (nm %in% pars) {
          return(rep(val, length(unique(df[[group]]))))
          # if this is one of pars then make starting value per group
        } else {
          return(val)
        } # else return one starting value
      })
      names(startingValuesList) <- nms
    } else { # non-grouped, just make it into a list
      startingValuesList <- as.list(startingValues)
    }
  } else {
    startingValuesList <- start
  }
  out[["formula"]] <- growthForm
  if (type == "nlrq") {
    out[["taus"]] <- tau
  }
  out[["start"]] <- startingValuesList
  out[["df"]] <- df
  out[["pcvrForm"]] <- form
  return(out)
}

#' example of using selfStart
#' `logistic self starter`
#' @examples
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' .initlogistic(simdf, "time", "y")
#'
#' @keywords internal
#' @noRd

.initlogistic <- function(df, x, y, int = FALSE) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if (nrow(xy) < 4) {
    stop("too few distinct input values to fit a logistic model")
  }
  z <- abs(xy[["y"]])
  ## transform to proportion, i.e. in (0,1) :
  rng <- range(z)
  dz <- diff(rng)
  z <- (z - rng[1L] + 0.05 * dz) / (1.1 * dz)
  xy[["z"]] <- log(z / (1 - z)) # logit transformation
  aux <- stats::coef(stats::lm(x ~ z, xy))
  pars <- stats::coef(stats::nls(y ~ 1 / (1 + exp((B - x) / C)),
    data = xy,
    start = list(B = aux[[1L]], C = aux[[2L]]),
    algorithm = "plinear", control = stats::nls.control(warnOnly = TRUE)
  ))
  start <- stats::setNames(pars[c(".lin", "B", "C")], c("A", "B", "C"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "A", "B", "C"))
  }
  return(start)
}


#' `4 parameter logistic self starter`
#' @examples
#' simdf <- growthSim("logistic4",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5), "D" = c(5, 20))
#' )
#' .initlogistic4(simdf, "time", "y")
#'
#' @keywords internal
#' @noRd

.initlogistic4 <- function(df, x, y, int = FALSE) {
  obs_min <- min(df[[y]], na.rm = TRUE)
  if (int) {
    df[[y]] <- df[[y]] - obs_min
  }
  start <- .initlogistic(df, x, y, int = int)
  start <- stats::setNames(append(start, obs_min), c(names(start), "D"))
  return(start)
}

#' `5 parameter logistic self starter`
#' @examples
#' simdf <- growthSim("logistic5",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11),
#'   "C" = c(3, 3.5), "D" = c(5, 20), "E" = c(1, 1.5))
#' )
#' .initlogistic5(simdf, "time", "y")
#'
#' @keywords internal
#' @noRd

.initlogistic5 <- function(df, x, y, int = FALSE) {
  obs_min <- min(df[[y]], na.rm = TRUE)
  if (int) {
    df[[y]] <- df[[y]] - obs_min
  }
  start <- .initlogistic(df, x, y, int = int)
  start <- stats::setNames(append(start, c(obs_min, 1)), c(names(start), "D", "E"))
  return(start)
}


#' `Goempertz self starter`
#'
#' .initgompertz(simdf, "time", "y")
#' @keywords internal
#' @noRd

.initgompertz <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if (nrow(xy) < 4) {
    stop("too few distinct input values to fit the Gompertz model")
  }
  xyL <- xy
  xyL$y <- log(abs(xyL$y))
  pars <- stats::NLSstAsymptotic(xyL)
  pars <- stats::coef(stats::nls(y ~ exp(-B * C^x),
    data = xy, start = c(
      B = pars[["b1"]],
      C = exp(-exp(pars[["lrc"]]))
    ),
    algorithm = "plinear", control = stats::nls.control(warnOnly = TRUE)
  ))
  start <- stats::setNames(pars[c(".lin", "B", "C")], c("A", "B", "C"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "A", "B", "C"))
  }
  return(start)
}


#' `Monomolecular self starter`
#' ex<-growthSim("monomolecular", n=20, t=25,
#'                  params = list("A"=c(200,160), "B"=c(0.08, 0.1)))
#' .initmonomolecular(ex, "time", "y")
#' @keywords internal
#' @noRd

.initmonomolecular <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if (nrow(xy) < 4) {
    stop("too few distinct input values to fit a monomolecular model")
  }
  z <- abs(xy[["y"]])
  ## transform to proportion, i.e. in (0,1) :
  rng <- range(z)
  dz <- diff(rng)
  z <- (z - rng[1L] + 0.05 * dz) / (1.1 * dz)
  xy[["z"]] <- z
  aux <- stats::coef(stats::lm(z ~ x, xy))
  pars <- stats::coef(stats::nls(y ~ 1 * exp(B * x),
    data = xy, start = list(B = aux[2L]),
    algorithm = "plinear", control = stats::nls.control(warnOnly = TRUE)
  ))
  start <- stats::setNames(pars[c(".lin", "B")], c("A", "B"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "A", "B"))
  }
  return(start)
}

#' `Linear self starter`
#' @examples
#' ex <- growthSim("linear",
#'   n = 20, t = 25,
#'   params = list("A" = c(1.1, 0.95))
#' )
#' .initlinear(ex, "time", "y")
#' @keywords internal
#' @noRd

.initlinear <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if (nrow(xy) < 2) {
    stop("too few distinct input values to fit a linear model")
  }
  pars <- stats::coef(stats::lm(y ~ x, xy))
  start <- stats::setNames(pars[c("x")], c("A"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "A"))
  }
  return(start)
}

#' `Logarithmic self starter`
#' @examples
#' ex <- growthSim("logarithmic",
#'   n = 20, t = 25,
#'   params = list("A" = c(1.1, 0.95))
#' )
#' .initlinear(ex, "time", "y")
#' @keywords internal
#' @noRd

.initlogarithmic <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if (nrow(xy) < 2) {
    stop("too few distinct input values to fit a logarithmic model")
  }
  pars <- stats::coef(stats::lm(y ~ log(x), xy))
  start <- stats::setNames(pars[c("log(x)")], c("A"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "A"))
  }
  return(start)
}


#' `power law self starter`
#' @examples
#' ex <- growthSim("power law",
#'   n = 20, t = 25,
#'   params = list("A" = c(16, 11), "B" = c(0.75, 0.7))
#' )
#' .initPowerLaw(df, "time", "y")
#' @keywords internal
#' @noRd

.initpowerlaw <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if (nrow(xy) < 3) {
    stop("too few distinct input values to fit a power law model")
  }
  aux <- stats::coef(stats::lm(y ~ x, xy))

  pars <- stats::coef(stats::nls(y ~ 1 * x^B,
    data = xy, start = list(B = aux[2L]),
    algorithm = "plinear", control = stats::nls.control(warnOnly = TRUE)
  ))
  start <- stats::setNames(pars[c(".lin", "B")], c("A", "B"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "A", "B"))
  }
  return(start)
}


#' `exponential self starter`
#' @examples
#' ex <- growthSim("exponential",
#'   n = 20, t = 25,
#'   params = list("A" = c(15, 20), "B" = c(0.095, 0.095))
#' )
#' .initexponential(ex, "time", "y")
#' @keywords internal
#' @noRd

.initexponential <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if (nrow(xy) < 3) {
    stop("too few distinct input values to fit a exponential model")
  }
  aux <- stats::coef(stats::lm(y ~ x, xy))

  pars <- stats::coef(stats::nls(y ~ 1 * exp(B * x),
    data = xy, start = list(B = aux[2L]),
    algorithm = "plinear", control = stats::nls.control(warnOnly = TRUE)
  ))
  start <- stats::setNames(pars[c(".lin", "B")], c("A", "B"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "A", "B"))
  }
  return(start)
}

#' `Extreme Value Distribution Self Starter (weibull, frechet, gumbel)`
#' @keywords internal
#' @noRd

.initweibull <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  if (nrow(xy) < 5) {
    stop("too few distinct input values to fit the EVD growth model")
  }
  rAsym <- stats::NLSstRtAsymptote(xy)
  pars <- stats::coef(stats::lm(
    log(-log((rAsym - y) / (rAsym - stats::NLSstLfAsymptote(xy)))) ~ log(x),
    data = xy, subset = x > 0
  ))
  start <- stats::setNames(c(rAsym, exp(pars) + c(1, 0)), c("A", "B", "C"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "A", "B", "C"))
  }
  return(start)
}
.initfrechet <- .initweibull
.initgumbel <- .initweibull

#' `bragg DRC self starter`
#' @examples
#' ex <- growthSim("bragg",
#'   n = 20, t = 25,
#'   params = list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
#' )
#' .initbragg(ex, "time", "y")
#' @keywords internal
#' @noRd

.initbragg <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  x <- xy[, "x"]
  y <- xy[, "y"]
  A <- max(y) # amplitude, conventionally d
  C <- x[which.max(y)] # position of midpoint, conventionally e
  pseudoY <- log((y + 1e-04) / A)
  pseudoX <- (x - C)^2
  dat <- data.frame(pseudoY = pseudoY, pseudoX = pseudoX)
  coefs <- coef(lm(pseudoY ~ pseudoX - 1, data = dat))
  B <- -coefs[1] # pseudo precision/slope at inflection, conventionally b
  start <- stats::setNames(c(B, A, C), c("B", "A", "C"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "B", "A", "C"))
  }
  return(start)
}

#' `lorentz DRC self starter`
#' @examples
#' ex <- growthSim("lorentz",
#'   n = 20, t = 25,
#'   params = list("A" = c(10, 15), "B" = c(0.01, 0.02), "C" = c(50, 60))
#' )
#' .initlorentz(ex, "time", "y")
#' @keywords internal
#' @noRd

.initlorentz <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  x <- xy[, "x"]
  y <- xy[, "y"]
  A <- max(y) # amplitude, conventionally d
  C <- x[which.max(y)] # position of midpoint, conventionally e
  pseudoY <- (A - y) / y
  pseudoX <- (x - C)^2
  dat <- data.frame(pseudoY = pseudoY, pseudoX = pseudoX)
  coefs <- coef(lm(pseudoY ~ pseudoX - 1, data = dat))
  B <- coefs[1] # pseudo precision/slope at inflection, conventionally b
  start <- stats::setNames(c(B, A, C), c("B", "A", "C"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "B", "A", "C"))
  }
  return(start)
}

#' `beta DRC self starter`
#' @examples
#' ex <- growthSim("beta",
#'   n = 20, t = 25,
#'   params = list("A" = 10, "B" = 1.2, "C" = 15, "D" = 8, "E" = 19)
#' )
#' .initbeta(ex, "time", "y")
#' @keywords internal
#' @noRd

.initbeta <- function(df, x, y, int) {
  if (int) {
    obs_min <- min(df[[y]], na.rm = TRUE)
    df[[y]] <- df[[y]] - obs_min
  }
  xy <- stats::sortedXyData(df[[x]], df[[y]])
  x <- xy[, "x"]
  y <- xy[, "y"]
  A <- max(y)
  C <- x[which.max(y)]
  firstidx <- min(which(y != 0))
  D <- ifelse(firstidx == 1,
    x[1],
    (x[firstidx] + x[(firstidx - 1)]) / 2
  )
  secidx <- max(which(y != 0))
  E <- ifelse(secidx == length(y),
    x[length(x)],
    (x[secidx] + x[(secidx + 1)]) / 2
  )
  start <- stats::setNames(c(A, 1, C, D, E), c("A", "B", "C", "D", "E"))
  if (int) {
    start <- stats::setNames(append(obs_min, start), c("I", "A", "B", "C", "D", "E"))
  }
  return(start)
}


#' `Define growth formulas`
#' @keywords internal
#' @noRd

.nlrq_form_logistic <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C")
  } else {
    total_pars <- c("A", "B", "C")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~I[] + (A[]/(1+exp((B[]-", x, ")/C[])))")
  } else {
    str_nf <- paste0(y, " ~ A[]/(1+exp((B[]-", x, ")/C[]))")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_logistic4 <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C", "D")
  } else {
    total_pars <- c("A", "B", "C", "D")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + D[] + (A[] - D[])/(1+exp((B[]-", x, ")/C[]))")
  } else {
    str_nf <- paste0(y, " ~ D[] + (A[] - D[])/(1+exp((B[]-", x, ")/C[]))")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_logistic5 <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C", "D", "E")
  } else {
    total_pars <- c("A", "B", "C", "D", "E")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + D[] + ((A[] - D[])/(1+exp((B[]-", x, ")/C[])) ^ E[])")
  } else {
    str_nf <- paste0(y, " ~ D[] + ((A[] - D[])/(1+exp((B[]-", x, ")/C[])) ^ E[])")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_gompertz <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C")
  } else {
    total_pars <- c("A", "B", "C")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + (A[]*exp(-B[]*exp(-C[]*", x, ")))")
  } else {
    str_nf <- paste0(y, " ~ A[]*exp(-B[]*exp(-C[]*", x, "))")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_doublelogistic <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C", "A2", "B2", "C2")
  } else {
    total_pars <- c("A", "B", "C", "A2", "B2", "C2")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(
      y, " ~ I[] + (A[]/(1+exp((B[]-", x, ")/C[]))",
      " + ((A2[]-A[]) /(1+exp((B2[]-", x, ")/C2[]))))"
    )
  } else {
    str_nf <- paste0(
      y, " ~ A[]/(1+exp((B[]-", x, ")/C[]))",
      " + ((A2[]-A[]) /(1+exp((B2[]-", x, ")/C2[])))"
    )
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_doublegompertz <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C", "A2", "B2", "C2")
  } else {
    total_pars <- c("A", "B", "C", "A2", "B2", "C2")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(
      y, " ~ I[] + (A[] * exp(-B[] * exp(-C[]*", x, "))",
      " + (A2[]-A[]) * exp(-B2[] * exp(-C2[]*(", x, "-B[]))))"
    )
  } else {
    str_nf <- paste0(
      y, " ~ A[] * exp(-B[] * exp(-C[]*", x, "))",
      " + (A2[]-A[]) * exp(-B2[] * exp(-C2[]*(", x, "-B[])))"
    )
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_monomolecular <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B")
  } else {
    total_pars <- c("A", "B")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, "~I[] + (A[]-A[]*exp(-B[]*", x, "))")
  } else {
    str_nf <- paste0(y, "~A[]-A[]*exp(-B[]*", x, ")")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_exponential <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B")
  } else {
    total_pars <- c("A", "B")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + (A[]*exp(B[]*", x, "))")
  } else {
    str_nf <- paste0(y, " ~ A[]*exp(B[]*", x, ")")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_linear <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A")
  } else {
    total_pars <- c("A")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + A[]*", x)
  } else {
    str_nf <- paste0(y, " ~ A[]*", x)
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_powerlaw <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B")
  } else {
    total_pars <- c("A", "B")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + (A[]*", x, "^B[])")
  } else {
    str_nf <- paste0(y, " ~ A[]*", x, "^B[]")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_gam <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (USEGROUP) {
    nf <- as.formula(paste0(y, " ~ bs(", x, ")*", group))
  } else {
    nf <- as.formula(paste0(y, " ~ bs(", x, ")"))
  }
  return(list("formula" = nf, "pars" = NULL))
}

.nlrqDecay <- function(form) {
  chars <- as.character(form)
  return(as.formula(paste0(chars[2], chars[1], "-(", chars[3], ")")))
}

.nlrq_form_frechet <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C")
  } else {
    total_pars <- c("A", "B", "C")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + (A[] * exp(-((", x, "-0)/C[])^(-B[])))")
  } else {
    str_nf <- paste0(y, " ~ A[] * exp(-((", x, "-0)/C[])^(-B[]))")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}


.nlrq_form_weibull <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C")
  } else {
    total_pars <- c("A", "B", "C")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + (A[] * (1-exp(-(", x, "/C[])^B[])))")
  } else {
    str_nf <- paste0(y, " ~ A[] * (1-exp(-(", x, "/C[])^B[]))")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_gumbel <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C")
  } else {
    total_pars <- c("A", "B", "C")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + (A[] * exp(-exp(-(", x, "-B[])/C[])))")
  } else {
    str_nf <- paste0(y, " ~ A[] * exp(-exp(-(", x, "-B[])/C[]))")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_logarithmic <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A")
  } else {
    total_pars <- c("A")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + A[] * log(", x, ")")
  } else {
    str_nf <- paste0(y, " ~ A[] * log(", x, ")")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_bragg <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C")
  } else {
    total_pars <- c("A", "B", "C")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + A[] * exp(-B[] * (", x, " - C[])^2)")
  } else {
    str_nf <- paste0(y, " ~ A[] * exp(-B[] * (", x, " - C[])^2)")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_lorentz <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C")
  } else {
    total_pars <- c("A", "B", "C")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(y, " ~ I[] + A[] / (1 + B[] * (", x, " - C[]) ^ 2)")
  } else {
    str_nf <- paste0(y, " ~ A[] / (1 + B[] * (", x, " - C[]) ^ 2)")
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}

.nlrq_form_beta <- function(x, y, USEGROUP, group, pars, int = FALSE) {
  if (int) {
    total_pars <- c("I", "A", "B", "C", "D", "E")
  } else {
    total_pars <- c("A", "B", "C", "D", "E")
  }
  if (is.null(pars)) {
    pars <- total_pars
  }
  if (int) {
    str_nf <- paste0(
      y, " ~ I[] + A[] * (((", x, " - D[]) / (C[] - D[])) * ((E[] - ", x,
      ") / (E[] - C[])) ^ ((E[] - C[]) / (C[] - D[]))) ^ B[]"
    )
  } else {
    str_nf <- paste0(
      y, " ~ A[] * (((", x, " - D[]) / (C[] - D[])) * ((E[] - ", x,
      ") / (E[] - C[])) ^ ((E[] - C[]) / (C[] - D[]))) ^ B[]"
    )
  }
  if (USEGROUP) {
    for (par in total_pars) {
      if (par %in% pars) {
        str_nf <- gsub(paste0(par, "\\[\\]"), paste0(par, "[", group, "]"), str_nf)
      } else {
        str_nf <- gsub(paste0(par, "\\[\\]"), par, str_nf)
      }
    }
    nf <- as.formula(str_nf)
  } else {
    nf <- as.formula(gsub("\\[|\\]", "", str_nf))
  }
  return(list("formula" = nf, "pars" = pars))
}
