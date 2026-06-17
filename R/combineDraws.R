#' Combine Draws From brms Models
#'
#' @description Helper function for binding draws from several \code{brms} models to make a data.frame
#' for use with \code{brms::hypothesis()}. This will also check that the draws are comparable using
#' basic model metrics.
#'
#' @param ... Some number of brmsfit objects and/or dataframes of draws
#' (should generally be the same type of model fit to different data)
#' @param names Optional vector of names for the models/data.frames. The default (NULL) will use
#' the object names of arguments to \code{...} as a vector of names to uniquely identify columns in
#' the output data frame. See details for using this with \code{do.call} on a list of models.
#' @param message Logical, should messages about possible problems be printed? Default is TRUE.
#' This will warn if models may not have converged, if there are different numbers of draws in
#' the objects, or if models have different formulations.
#' @keywords brms
#' @importFrom methods is
#' @returns A data.frame of posterior draws, labeled to show which object they come from.
#' @details
#' If you fit models as part of a loop/apply function and end up with a list of models it may be
#' helpful to call this function on the list. In that case object names are not parsed well from
#' the list by default so passing the \code{names} argument is helpful and can be done as
#' \code{do.call(combineDraws, c(fits, list(names = names(fits))))}.
#'
#' @examples
#' # note that this example will fit several models using Stan and may run slowly.
#' \donttest{
#' data(fit)
#' fit_1 <- fit_2 <- fit
#' x <- combineDraws(fit_1, fit_2, fit)
#' draws_1 <- as.data.frame(fit_1)
#' draws_1 <- draws_1[, grepl("^b_", colnames(draws_ef))]
#' x2 <- combineDraws(fit_2, draws_1)
#' }
#'
#' @return Returns a dataframe of posterior draws.
#' @export

combineDraws <- function(..., names = NULL, message = TRUE) {
  objects <- list(...)
  if (!all(unlist(lapply(objects, function(m) {
    return(methods::is(m, "brmsfit") | methods::is(m, "data.frame"))
  })))) {
    stop("Only brmsfit objects and data frames are accepted")
  }
  if (is.null(names)) {
    names <- sapply(substitute(list(...)), deparse)[-1]
  }
  models <- objects[unlist(lapply(objects, function(m) {
    return(methods::is(m, "brmsfit"))
  }))]
  model_names <- names[unlist(lapply(objects, function(m) {
    return(methods::is(m, "brmsfit"))
  }))]
  supplied_draw_dfs <- objects[unlist(lapply(objects, function(m) {
    return(methods::is(m, "data.frame"))
  }))]
  df_names <- names[unlist(lapply(objects, function(m) {
    return(methods::is(m, "data.frame"))
  }))]

  max_fit_draws <- max(unlist(lapply(models, function(m) {
    return(nrow(as.data.frame(m)))
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
        message(paste0(
          mn, " has Rhat values >1.05 for some model parameters.",
          "See ?barg for possible improvements.\n"
        ))
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
    return(draws)
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
      return(df)
    })
    new_draws <- do.call(cbind, args = list(supplied_draw_dfs, new_draws))
  }
  return(new_draws)
}
