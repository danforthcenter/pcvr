#' Helper to make pcvr default priors from several kinds of input
#' @keywords internal
#' @noRd

.makePriors <- function(priors, pars, df, group, USEGROUP, sigma, family, formula) {
  if (is.null(priors)) {
    prior <- .explicitDefaultPrior(formula, df, family)
    return(prior)
  }
  #* `if priors is a brmsprior`
  if (any(methods::is(priors, "brmsprior"))) {
    return(priors)
  }
  #* `if priors is a numeric vector`
  priors <- .fixNumericPriors(priors, pars)
  #* `if priors is a list`
  formatListPriorsRes <- .formatListPriors(priors, pars, df, group, USEGROUP)
  priors <- formatListPriorsRes[["priors"]]
  # might need something like "what group does this value come from"?
  groupedPriors <- formatListPriorsRes[["groupedPriors"]]
  #* `if multiple group variables then make interaction priors as Normal (0, sd)`
  group_interaction_priors <- .groupInteractionPriors(group, formula, df, family, priors = priors)
  #* `make lookup table for grouping and values of groups`
  lookup <- data.frame(
    group = unlist(lapply(group, function(x) rep(x, length(unique(df[[x]]))))),
    value = unlist(lapply(group, function(x) unique(df[[x]])))
  )
  #* `Arrange priors to match pars explicitly`
  if (length(setdiff(pars, names(priors))) > 0) {
    specified_pars <- intersect(names(priors), pars)
    unspecified_pars <- setdiff(pars, names(priors))
    priors <- priors[specified_pars]
    pars <- c(specified_pars, unspecified_pars)
  } else {
    priors <- priors[pars]
  }
  #* `replicate lookup table for pars`
  lookup <- do.call(rbind, lapply(pars, function(pr) {
    lookup$par <- pr
    lookup
  }))
  #* `Make stan strings`
  priorStanStrings <- .stanStringHelper(priors, pars, USEGROUP)
  #* `get default priors for intercept only distributional parameters`
  prior <- .initializePriorObject(sigma, family)
  #* `add priors for estimated parameters`
  for (i in seq_along(priorStanStrings)) {
    nm <- names(priorStanStrings)[i]
    dist <- priorStanStrings[[nm]]
    pr <- strsplit(nm, "_")[[1]][1]
    if (USEGROUP && groupedPriors) { # if there are groups and they have different priors
      grp <- group
      grp_level <- strsplit(nm, "_")[[1]][2]
      if (length(group) > 1) {
        grp <- lookup[i, "group"]
        grp_level <- lookup[i, "value"] # should always match previous grp_level definition
      }
      gr <- paste0(grp, grp_level)
      prior <- prior + brms::set_prior(dist, coef = gr, nlpar = pr)
      # currently cannot set lb for prior with coef
      # there is a clunky workaround but it wouldn't work with expected data types
      # https://github.com/paul-buerkner/brms/issues/86
    } else {
      lb <- ifelse(grepl("changePoint|I$", pr), NA, 0)
      prior <- prior + brms::set_prior(dist, nlpar = pr, lb = lb)
    }
  }
  prior <- prior[-1, ] # remove flat prior on b
  prior <- rbind(prior, group_interaction_priors)
  prior <- unique(prior)
  # could add intercept term prior here
  return(prior)
}

#' Helper function to make Normal priors on interaction terms between grouping variables
#'
#' @keywords internal
#' @noRd

.groupInteractionPriors <- function(group, formula, df, family, priors) {
  if (length(group) > 1) {
    default_prior <- .explicitDefaultPrior(formula, df, family)
    default_interaction_prior <- default_prior[grepl(":", default_prior$coef), ]
    tenth_of_priors <- lapply(priors, function(x) {
      mean(x) / 10
    })
    for (nlp in unique(default_interaction_prior$nlpar)) {
      sd <- ifelse(nlp %in% names(tenth_of_priors), tenth_of_priors[[nlp]], 3)
      default_interaction_prior[
        default_interaction_prior$nlpar == nlp, "prior"
      ] <- paste0("normal(0, ", sd, ")")
    }
    return(default_interaction_prior)
  } else {
    return(NULL)
  }
}

#' Helper function to fix numeric priors
#'
#' @keywords internal
#' @noRd

.fixNumericPriors <- function(priors, pars) {
  #* `if priors is a numeric vector`
  if (is.numeric(priors)) {
    if (length(priors) == length(pars)) {
      warning("Assuming that prior is in order: ", paste0(pars, collapse = ", "))
      priors <- as.list(priors)
      names(priors) <- pars
    } else {
      stop(paste0(
        "`priors` is length ", length(priors), " while the specified model requires ",
        length(pars), " parameters."
      ))
    }
  }
  return(priors)
}

#' Helper function to standardize list priors
#'
#' @keywords internal
#' @noRd

.formatListPriors <- function(priors, pars, df, group, USEGROUP) {
  if (is.list(priors)) {
    if (is.null(names(priors))) {
      warning("Assuming that each element in priors is in order: ", paste0(pars, collapse = ", "))
      names(priors) <- pars
    }
    priors <- priors[!grepl("fixedChangePoint", names(priors))]
    if (!all(pars %in% names(priors))) {
      warning(paste0(
        "Parameter names and prior names do not match. Priors include ",
        paste(setdiff(names(priors), pars), collapse = ", "),
        "... and parameters include ",
        paste(setdiff(pars, names(priors)), collapse = ", "),
        "... Please rename the misspecified priors."
      ))
    }
    groupedPriors <- any(unlist(lapply(priors, length)) > 1)
    # if any prior has multiple means then groupedPriors is TRUE

    if (groupedPriors) { # if more than one value is specified per parameter
      l <- sum(unlist(lapply(group, function(grp) {
        length(unique(df[[grp]]))
      })))
      ml <- max(c(l, unlist(lapply(priors, length))))
      priors <- lapply(priors, function(p) rep(p, length.out = ml))
      if (any(unlist(lapply(priors, function(p) !is.null(names(p)))))) {
        # if any inner values are named then apply that to all priors
        wch <- which(unlist(lapply(priors, function(p) !is.null(names(p)))))
        nms <- names(priors[[wch]])
        for (i in seq_along(priors)) {
          names(priors[[i]]) <- nms
        }
      }
      if (any(unlist(lapply(priors, function(p) is.null(names(p)))))) {
        # if no inner values were named
        for (i in seq_along(priors)) {
          names(priors[[i]]) <- unique(interaction(df[, group]))
        }
      }
    } else { # else is for prior of length 1 for each element,
      # in which case they need to replicated per groups
      # this should also handle non-grouped formulae
      l <- sum(unlist(lapply(group, function(grp) {
        length(unique(df[[grp]]))
      })))
      priors <- lapply(priors, rep, length.out = l)
      nms <- unlist(lapply(group, function(grp) {
        unique(df[[grp]])
      }))
      if (USEGROUP) {
        for (i in seq_along(priors)) {
          names(priors[[i]]) <- nms
        }
      }
    }
  }
  return(list("priors" = priors, "groupedPriors" = groupedPriors))
}


#' Helper function to write stan priors
#'
#' @keywords internal
#' @noRd

.stanStringHelper <- function(priors, pars, USEGROUP) {
  if (!is.null(pars)) {
    priorStanStrings <- lapply(pars, function(par) {
      if (!grepl("changePoint|I$", par)) {
        paste0("lognormal(log(", priors[[par]], "), 0.25)") # growth parameters are LN
      } else {
        # changepoints/intercepts are T_5(mu, mu / 5) by default
        paste0("student_t(5,", priors[[par]], ", ", abs(priors[[par]] / 5), ")")
      }
    })
    priorStanStrings <- unlist(priorStanStrings)
    parNames <- rep(names(priors), each = length(priors[[1]]))
    if (USEGROUP) {
      groupNames <- rep(names(priors[[1]]), length.out = length(priorStanStrings))
      names(priorStanStrings) <- paste(parNames, groupNames, sep = "_")
    } else {
      names(priorStanStrings) <- parNames
    }
    return(priorStanStrings)
  }
}


#' Helper function to write stan priors
#'
#' @keywords internal
#' @noRd

.initializePriorObject <- function(sigma, family) {
  int_only_dpars <- names(sigma[which(sigma == "not_estimated")])
  if (length(int_only_dpars) >= 1) {
    int_dpars_form <- as.formula(paste0(paste(int_only_dpars, collapse = "+"), "~1"))
  } else {
    int_dpars_form <- NULL
  }
  smooth_dpars <- names(sigma[which(sigma %in% c("gam", "spline"))])
  if (length(smooth_dpars) >= 1) {
    smooth_dpars_form <- as.formula(paste0(paste(smooth_dpars, collapse = "+"), "~s(x)"))
  } else {
    smooth_dpars_form <- NULL
  }
  flist <- list(int_dpars_form, smooth_dpars_form)
  flist <- flist[!unlist(lapply(flist, is.null))]
  if (length(flist) == 0) {
    flist <- NULL
  }
  gp <- brms::get_prior(brms::bf(y ~ x, flist = flist),
    data = data.frame(y = 1:100, x = 1:100), family = family
  )
  prior <- rbind(
    gp[1, ], gp[gp$dpar %in% smooth_dpars & gp$class == "Intercept", ],
    gp[gp$dpar %in% int_only_dpars, ]
  )
  return(prior)
}

#' Helper function to explicitly return default priors from get_prior
#'
#' @keywords internal
#' @noRd

.explicitDefaultPrior <- function(formula, df, family) {
  gp <- brms::get_prior(formula = formula, data = df, family = family)
  return(gp)
}

#' Helper function to reformat sigma argument in brmSS
#'
#' @keywords internal
#' @noRd

.sigmaHelper <- function(sigma, dpars, family, models) {
  if (is.null(sigma)) {
    sigma <- lapply(dpars, function(i) {
      "int"
    })
  }
  if (methods::is(sigma, "formula")) {
    sigma <- list(sigma)
  }

  if (length(sigma) > length(dpars)) {
    stop(paste0(
      "sigma contains ", length(sigma), " formulas.",
      "The specified family (", family, ") only has ", length(dpars),
      " valid additional distributional parameters (", paste0(dpars, collapse = ", "), ")."
    ))
  } else if (length(sigma) < length(dpars)) {
    n_to_add <- length(dpars) - length(sigma)
    sigma <- append(sigma, lapply(1:n_to_add, function(i) paste0("not_estimated")))
    names(sigma) <- dpars
  } else { # same length
    names(sigma) <- dpars
  }
  # here I am foregoing pattern matching so that it is simpler to check for intercepts later.
  if (!any(grepl("\\+", sigma))) { # no distributional changepoint models
    sigma <- lapply(sigma, identity)
  }

  return(sigma)
}

#' Helper function to match growth model
#'
#' @keywords internal
#' @noRd

.matchGrowthModel <- function(model, models) {
  if (!grepl("\\+", model)) {
    if (grepl("decay", model)) {
      decay <- TRUE
      model <- trimws(gsub("decay", "", model))
    } else {
      decay <- FALSE
    }
    matched_model <- match.arg(model, models)
  } else {
    matched_model <- model
    decay <- FALSE
  }
  return(list("model" = matched_model, "decay" = decay))
}



#' Helper function for logistic brms formulas
#'
#' @keywords internal
#' @noRd

.brms_form_logistic <- function(x, y, group, dpar = FALSE,
                                nTimes = NULL, useGroup = TRUE, prior = NULL, int = FALSE, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + ", y, "A/(1+exp((",
        y, "B-", x, ")/", y, "C))"
      )))
      pars <- paste0(y, LETTERS[c(1:3, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "A/(1+exp((",
        y, "B-", x, ")/", y, "C))"
      )))
      pars <- paste0(y, LETTERS[1:3])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + (A/(1+exp((B-", x, ")/C)))"))
      pars <- LETTERS[c(1:3, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A/(1+exp((B-", x, ")/C))"))
      pars <- LETTERS[1:3]
    }
  }
  return(list(form = form, pars = pars))
}
#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd

.brms_form_gompertz <- function(x, y, group, dpar = FALSE,
                                nTimes = NULL, useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + (",
        y, "A*exp(-", y, "B*exp(-", y, "C*", x, ")))"
      )))
      pars <- paste0(y, LETTERS[c(1:3, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ",
        y, "A*exp(-", y, "B*exp(-", y, "C*", x, "))"
      )))
      pars <- paste0(y, LETTERS[1:3])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + (A*exp(-B*exp(-C*", x, ")))"))
      pars <- LETTERS[c(1:3, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A*exp(-B*exp(-C*", x, "))"))
      pars <- LETTERS[1:3]
    }
  }
  return(list(form = form, pars = pars))
}
#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_doublelogistic <- function(x, y, group, dpar = FALSE,
                                      nTimes = NULL, useGroup = TRUE,
                                      prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + (", y,
        "A/(1+exp((", y, "B-", x, ")/", y, "C)) + ((",
        y, "A2-", y, "A) /(1+exp((", y, "B2-", x,
        ")/", y, "C2))))"
      )))
      pars <- paste0(y, c("I", "A", "B", "C", "A2", "B2", "C2"))
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y,
        "A/(1+exp((", y, "B-", x, ")/", y, "C)) + ((",
        y, "A2-", y, "A) /(1+exp((", y, "B2-", x,
        ")/", y, "C2)))"
      )))
      pars <- paste0(y, c("A", "B", "C", "A2", "B2", "C2"))
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(
        y, " ~ I + (A/(1+exp((B-", x, ")/C)) + ((A2-A) /(1+exp((B2-", x, ")/C2))))"
      ))
      pars <- c("I", "A", "B", "C", "A2", "B2", "C2")
    } else {
      form <- stats::as.formula(paste0(
        y, " ~ A/(1+exp((B-", x, ")/C)) + ((A2-A) /(1+exp((B2-", x, ")/C2)))"
      ))
      pars <- c("A", "B", "C", "A2", "B2", "C2")
    }
  }
  return(list(form = form, pars = pars))
}
#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_doublegompertz <- function(x, y, group, dpar = FALSE,
                                      nTimes = NULL, useGroup = TRUE,
                                      prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + (", y, "A * exp(-", y, "B * exp(-", y,
        "C*", x, ")) + (", y, "A2-", y, "A) * exp(-", y,
        "B2 * exp(-", y, "C2*(", x, "-", y, "B))))"
      )))
      pars <- paste0(y, c("I", "A", "B", "C", "A2", "B2", "C2"))
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "A * exp(-", y, "B * exp(-", y,
        "C*", x, ")) + (", y, "A2-", y, "A) * exp(-", y,
        "B2 * exp(-", y, "C2*(", x, "-", y, "B)))"
      )))
      pars <- paste0(y, c("A", "B", "C", "A2", "B2", "C2"))
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(
        y, " ~ I + (A * exp(-B * exp(-C*", x,
        ")) + (A2-A) * exp(-B2 * exp(-C2*(", x, "-B))))"
      ))
      pars <- c("I", "A", "B", "C", "A2", "B2", "C2")
    } else {
      form <- stats::as.formula(paste0(
        y, " ~ A * exp(-B * exp(-C*", x,
        ")) + (A2-A) * exp(-B2 * exp(-C2*(", x, "-B)))"
      ))
      pars <- c("A", "B", "C", "A2", "B2", "C2")
    }
  }
  return(list(form = form, pars = pars))
}
#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_monomolecular <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                                     useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + (",
        y, "A-", y, "A*exp(-", y, "B*", x, "))"
      )))
      pars <- paste0(y, LETTERS[c(1:2, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(y, " ~ ", y, "A-", y, "A*exp(-", y, "B*", x, ")")))
      pars <- paste0(y, LETTERS[1:2])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, "~I + (A-A*exp(-B*", x, "))"))
      pars <- LETTERS[c(1:2, 9)]
    } else {
      form <- stats::as.formula(paste0(y, "~A-A*exp(-B*", x, ")"))
      pars <- LETTERS[1:2]
    }
  }
  return(list(form = form, pars = pars))
}
#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_exponential <- function(x, y, group, dpar = FALSE,
                                   nTimes = NULL, useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ",
        y, "I + (", y, "A*exp(", y, "B*", x, "))"
      )))
      pars <- paste0(y, LETTERS[c(1:2, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(y, " ~ ", y, "A*exp(", y, "B*", x, ")")))
      pars <- paste0(y, LETTERS[1:2])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + (A*exp(B*", x, "))"))
      pars <- LETTERS[c(1:2, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A*exp(B*", x, ")"))
      pars <- LETTERS[1:2]
    }
  }
  return(list(form = form, pars = pars))
}
#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_linear <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                              useGroup = TRUE, prior = NULL, int, force_nl = FALSE) {
  if (dpar) {
    if (!is.null(prior) && any(grepl(paste0(y, "A"), names(prior))) || force_nl) {
      #* use non-linear parameterization with subA
      if (int) {
        form <- brms::nlf(stats::as.formula(paste0(y, " ~ ", y, "I + (", y, "A", "*", x, ")")))
        pars <- c(paste0(y, c("I", "A")))
      } else {
        form <- brms::nlf(stats::as.formula(paste0(y, " ~ ", y, "A", "*", x)))
        pars <- c(paste0(y, "A"))
      }
    } else {
      #* linear parameterization using x directly
      if (int) {
        form <- brms::nlf(as.formula(paste0(
          y, " ~ ", y, "I + (", x, "+", x, ":",
          paste(group, collapse = ":"), ")"
        )))
        pars <- paste0(y, "I")
      } else {
        form <- as.formula(paste0(y, " ~ ", x, "+", x, ":", paste(group, collapse = ":")))
        pars <- c()
      }
    }
  } else { # non-dpar option, with or without intercept
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + A*", x))
      pars <- c("I", "A")
    } else {
      form <- stats::as.formula(paste0(y, " ~ A*", x))
      pars <- c("A")
    }
  }
  return(list(form = form, pars = pars))
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_logarithmic <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                                   useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(y, " ~ ", y, "I + (", y, "A*log(", x, "))")))
      pars <- paste0(y, LETTERS[c(1, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(y, " ~ ", y, "A*log(", x, ")")))
      pars <- paste0(y, "A")
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + (A*log(", x, "))"))
      pars <- LETTERS[c(1, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A*log(", x, ")"))
      pars <- "A"
    }
  }
  return(list(form = form, pars = pars))
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_powerlaw <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                                useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(y, " ~ ", y, "I + (", y, "A*", x, "^", y, "B)")))
      pars <- paste0(y, LETTERS[c(1:2, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(y, " ~ ", y, "A*", x, "^", y, "B")))
      pars <- paste0(y, LETTERS[1:2])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + (A*", x, "^B)"))
      pars <- LETTERS[c(1:2, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A*", x, "^B"))
      pars <- LETTERS[1:2]
    }
  }
  return(list(form = form, pars = pars))
}
#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_gam <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                           useGroup = TRUE, prior = NULL, int, ...) {
  if (useGroup) {
    by <- paste0(", by = ", paste(group, collapse = ".")) # special variable that is made if there are
    # multiple groups and a gam involved.
  } else {
    by <- NULL
  }
  if (nTimes < 11) {
    k <- paste0(", k = ", nTimes)
  } else {
    k <- NULL
  }
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(y, " ~ ", y, "I + s(", x, by, k, ")")))
      pars <- paste0(y, "I")
    } else {
      form <- stats::as.formula(paste0(y, " ~ s(", x, by, k, ")"))
      pars <- NULL
    }
  } else {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(y, " ~ I + s(", x, by, k, ")")))
      pars <- "I"
    } else {
      form <- stats::as.formula(paste0(y, " ~ s(", x, by, k, ")"))
      pars <- NULL
    }
  }
  return(list(form = form, pars = pars))
}
#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_int <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                           useGroup = TRUE, prior = NULL, int, ...) {
  if (useGroup) {
    rhs <- paste0("0 + ", paste(group, collapse = "+"))
  } else {
    rhs <- paste0("1")
  }
  form <- stats::as.formula(paste0(y, " ~ ", rhs))
  pars <- c()
  return(list(form = form, pars = pars))
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd


.brms_form_not_estimated <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                                     useGroup = TRUE, prior = NULL, int, ...) {
  form <- stats::as.formula(paste0(y, " ~ 1"))
  pars <- c()
  return(list(form = form, pars = pars))
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd

.brms_form_decay <- function(formList, int = FALSE) {
  modelForm <- formList$form
  chars <- as.character(modelForm)
  if (!int) {
    formList$form <- as.formula(paste0(chars[2], chars[1], "-(", chars[3], ")"))
  } else {
    rhs <- chars[3]
    rhs <- trimws(gsub("I\\s?\\+", "", rhs))
    formList$form <- as.formula(paste0(chars[2], chars[1], "I - (", rhs, ")"))
  }
  formList
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd

.brms_form_frechet <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                               useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + (", y, "A * exp(-((", x, "-0)/", y, "C)^(-", y, "B)))"
      )))
      pars <- paste0(y, LETTERS[c(1:3, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "A * exp(-((", x, "-0)/", y, "C)^(-", y, "B))"
      )))
      pars <- paste0(y, LETTERS[1:3])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + (A * exp(-((", x, "-0)/C)^(-B)))"))
      pars <- LETTERS[c(1:3, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A * exp(-((", x, "-0)/C)^(-B))"))
      pars <- LETTERS[1:3]
    }
  }
  return(list(form = form, pars = pars))
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd

.brms_form_weibull <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                               useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + (", y, "A * (1-exp(-(", x, "/", y, "C)^", y, ")))"
      )))
      pars <- paste0(y, LETTERS[c(1:3, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "A * (1-exp(-(", x, "/", y, "C)^", y, "))"
      )))
      pars <- paste0(y, LETTERS[1:3])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + (A * (1-exp(-(", x, "/C)^B)))"))
      pars <- LETTERS[c(1:3, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A * (1-exp(-(", x, "/C)^B))"))
      pars <- LETTERS[1:3]
    }
  }
  return(list(form = form, pars = pars))
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd

.brms_form_gumbel <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                              useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + (", y, "A * exp(-exp( -(", x, "-", y, "B)/", y, "C)))"
      )))
      pars <- paste0(y, LETTERS[c(1:3, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "A * exp(-exp( -(", x, "-", y, "B)/", y, "C))"
      )))
      pars <- paste0(y, LETTERS[1:3])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + (A * exp(-exp( -(", x, "-B)/C)))"))
      pars <- LETTERS[c(1:3, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A * exp(-exp( -(", x, "-B)/C))"))
      pars <- LETTERS[1:3]
    }
  }
  return(list(form = form, pars = pars))
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd

.brms_form_bragg <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                             useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + ", y, "A * exp(-", y, "B * (", x, " - ", y, "C)^2)"
      )))
      pars <- paste0(y, LETTERS[c(1:3, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "A * exp(-", y, "B * (", x, " - ", y, "C)^2)"
      )))
      pars <- paste0(y, LETTERS[1:3])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + A * exp(-B * (", x, " - C)^2)"))
      pars <- LETTERS[c(1:3, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A * exp(-B * (", x, " - C)^2)"))
      pars <- LETTERS[1:3]
    }
  }
  return(list(form = form, pars = pars))
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd

.brms_form_lorentz <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                               useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + ", y, "A / (1 + ", y, "B * (", x, " - ", y, "C) ^ 2)"
      )))
      pars <- paste0(y, LETTERS[c(1:3, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "A / (1 + ", y, "B * (", x, " - ", y, "C) ^ 2)"
      )))
      pars <- paste0(y, LETTERS[1:3])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(y, " ~ I + A / (1 + B * (", x, " - C) ^ 2)"))
      pars <- LETTERS[c(1:3, 9)]
    } else {
      form <- stats::as.formula(paste0(y, " ~ A / (1 + B * (", x, " - C) ^ 2)"))
      pars <- LETTERS[1:3]
    }
  }
  return(list(form = form, pars = pars))
}

#' Helper function for brms formulas
#'
#' @keywords internal
#' @noRd

.brms_form_beta <- function(x, y, group, dpar = FALSE, nTimes = NULL,
                            useGroup = TRUE, prior = NULL, int, ...) {
  if (dpar) {
    if (int) {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "I + ", y, "A * (((", x, " - ", y, "D) / (", y, "C - ", y, "D)) * ((", y, "E - ", x,
        ") / (", y, "E - ", y, "C)) ^ ((", y, "E - ", y, "C) / (", y, "C - ", y, "D))) ^ ", y, "B"
      )))
      pars <- paste0(y, LETTERS[c(1:5, 9)])
    } else {
      form <- brms::nlf(stats::as.formula(paste0(
        y, " ~ ", y, "A * (((", x, " - ", y, "D) / (", y, "C - ", y, "D)) * ((", y, "E - ", x,
        ") / (", y, "E - ", y, "C)) ^ ((", y, "E - ", y, "C) / (", y, "C - ", y, "D))) ^ ", y, "B"
      )))
      pars <- paste0(y, LETTERS[1:5])
    }
  } else {
    if (int) {
      form <- stats::as.formula(paste0(
        y, " ~ I + A * (((", x, " - D) / (C - D)) * ((E - ", x,
        ") / (E - C)) ^ ((E - C) / (C - D))) ^ B"
      ))
      pars <- LETTERS[c(1:5, 9)]
    } else {
      form <- stats::as.formula(paste0(
        y, " ~ A * (((", x, " - D) / (C - D)) * ((E - ", x,
        ") / (E - C)) ^ ((E - C) / (C - D))) ^ B"
      ))
      pars <- LETTERS[1:5]
    }
  }
  return(list(form = form, pars = pars))
}
