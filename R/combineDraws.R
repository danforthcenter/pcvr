#' Combine Draws From brms Models
#'
#' @description Helper function for binding draws from several \code{brms} models to make a data.frame
#' for use with \code{brms::hypothesis()}. This will also check that the draws are comparable using
#' basic model metrics.
#'
#' @param ... Some number of brmsfit objects and/or dataframes of draws
#' (should generally be the same type of model fit to different data)
#' @param message Logical, should messages about possible problems be printed? Default is TRUE.
#' This will warn if models may not have converged, if there are different numbers of draws in
#' the objects, or if models have different formulations.
#' @keywords brms
#' @importFrom methods is
#' @examples
#' # note that this example will fit several bayesian models and may run for several minutes.
#' \donttest{
#'   simdf <- growthSim("logistic",
#'     n = 20, t = 25,
#'     params = list(
#'       "A" = c(200, 160, 220, 200, 140, 300),
#'       "B" = c(13, 11, 10, 9, 16, 12),
#'       "C" = c(3, 3.5, 3.2, 2.8, 3.3, 2.5)
#'     )
#'   )
#'   ss_ab <- growthSS(
#'     model = "logistic", form = y ~ time | id / group,
#'     sigma = "logistic", df = simdf[simdf$group %in% c("a", "b"), ],
#'     start = list("A" = 130, "B" = 12, "C" = 3,
#'                  "sigmaA" = 15, "sigmaB" = 10, "sigmaC" = 3), type = "brms"
#'   )
#'
#'   ss_cd <- growthSS(
#'     model = "logistic", form = y ~ time | id / group,
#'     sigma = "logistic", df = simdf[simdf$group %in% c("c", "d"), ],
#'     start = list("A" = 130, "B" = 12, "C" = 3,
#'                  "sigmaA" = 15, "sigmaB" = 10, "sigmaC" = 3), type = "brms"
#'   )
#'
#'   ss_ef <- growthSS(
#'     model = "logistic", form = y ~ time | id / group,
#'     sigma = "logistic", df = simdf[simdf$group %in% c("e", "f"), ],
#'     start = list("A" = 130, "B" = 12, "C" = 3,
#'                  "sigmaA" = 15, "sigmaB" = 10, "sigmaC" = 3), type = "brms"
#'   )
#'   ss_ef2 <- growthSS(
#'     model = "gompertz", form = y ~ time | id / group,
#'     sigma = "logistic", df = simdf[simdf$group %in% c("e", "f"), ],
#'     start = list("A" = 130, "B" = 12, "C" = 3,
#'                  "sigmaA" = 15, "sigmaB" = 10, "sigmaC" = 3), type = "brms"
#'   )
#'
#'
#'   fit_ab <- fitGrowth(ss_ab, chains = 1, cores = 1, iter = 1000)
#'   fit_ab2 <- fitGrowth(ss_ab, chains = 1, cores = 1, iter = 1200)
#'   fit_cd <- fitGrowth(ss_cd, chains = 1, cores = 1, iter = 1000)
#'   fit_ef <- fitGrowth(ss_ef, chains = 1, cores = 1, iter = 1000)
#'   fit_ef2 <- fitGrowth(ss_ef2, chains = 1, cores = 1, iter = 1000)
#'
#'   x <- combineDraws(fit_ab, fit_cd, fit_ef)
#'   draws_ef <- as.data.frame(fit_ef)
#'   draws_ef <- draws_ef[, grepl("^b_", colnames(draws_ef))]
#'   x2 <- combineDraws(fit_ab2, fit_cd, draws_ef)
#'   x3 <- combineDraws(fit_ab, fit_cd, fit_ef2)
#' }
#'
#' @return Returns a dataframe of posterior draws.
#' @export




combineDraws <- function(..., message = TRUE) {
  objects <- list(...)
  if (!all(unlist(lapply(objects, function(m) {
    methods::is(m, "brmsfit") | methods::is(m, "data.frame")
  })))) {
    stop("Only brmsfit objects and data frames are accepted")
  }

  obj_names <- sapply(substitute(list(...)), deparse)[-1]
  models <- objects[unlist(lapply(objects, function(m) {
    methods::is(m, "brmsfit")
  }))]
  model_names <- obj_names[unlist(lapply(objects, function(m) {
    methods::is(m, "brmsfit")
  }))]
  supplied_draw_dfs <- objects[unlist(lapply(objects, function(m) {
    methods::is(m, "data.frame")
  }))]
  df_names <- obj_names[unlist(lapply(objects, function(m) {
    methods::is(m, "data.frame")
  }))]

  max_fit_draws <- max(unlist(lapply(models, function(m) {
    nrow(as.data.frame(m))
  })))
  max_nrow_supplied <- max(c(0, unlist(lapply(supplied_draw_dfs, nrow))))
  limit_size <- max(c(max_fit_draws, max_nrow_supplied))

  #* `check that formulae are the same`
  if (message) {
    formulae <- unlist(lapply(models, function(m) {
      x <- as.character(m$formula$formula)
      form <- paste0(x[2], x[1], x[3])
      return(form)
    }))
    names(formulae) <- model_names
    if (length(unique(formulae)) > 1) {
      message("Some of these models have different growth formulas, consider if this is what you want.")
      message(paste0(paste(names(formulae), formulae, sep = ": "), collapse = ", "))
    }
  }
  #* `get and bind draws from models`
  new_draws <- do.call(cbind, lapply(seq_along(models), function(i) {
    m <- models[[i]]
    mn <- model_names[[i]]
    d <- as.data.frame(m)
    draws <- d[, grepl("^b_", colnames(d))]
    colnames(draws) <- gsub("^b", mn, colnames(draws))

    if (message) {
      rhats <- brms::rhat(m)
      rhats <- rhats[grepl("^b_", names(rhats))]
      if (any(rhats > 1.05)) {
        message(paste0(mn, " has Rhat values >1.05 for some model parameters.",
                       "See ?barg for possible improvements.\n"))
      }
    }

    if (nrow(draws) < limit_size) {
      if (message) {
        message(paste0(mn, " has fewer than ", limit_size, " draws and will be padded with ",
                       limit_size - nrow(draws), " NAs\n"))
      }
      draws[(nrow(draws) + 1):limit_size, ] <- NA
      draws
    }
    draws
  }))
  #* `bind any other dataframes of draws`
  if (length(supplied_draw_dfs) > 0) {
    supplied_draw_dfs <- lapply(seq_along(supplied_draw_dfs), function(i) {
      df <- supplied_draw_dfs[[i]]
      dn <- df_names[[i]]
      if (nrow(df) < limit_size) {
        if (message) {
          message(paste0(dn, " has fewer than ", limit_size, " draws and will be padded with ",
                         limit_size - nrow(df), " NAs\n"))
        }
        df[(nrow(df) + 1):limit_size, ] <- NA
      }
      df
    })
    new_draws <- do.call(cbind, args = list(supplied_draw_dfs, new_draws))
  }
  new_draws
}
