#' Class \code{pcvrss} for models specified in \code{pcvr}.
#'
#' Models specified by \link{growthSS} or \link{mvSS} are represented by a \code{pcvrss} object,
#' which contains the model type, formulas, starting values or priors, the data for the model
#' to use, and the model backend to use.
#'
#' @name pcvrss-class
#' @aliases pcvrss
#' @docType class
#'
#' @details
#' See \code{methods(class = "pcvrss")} for an overview of available methods.
#'
#' @slot formula The formula that will be used to fit the model.
#' @slot prior Priors if the model is a Bayesian model (ie using the brms backend).
#' @slot initfun Initialization function if the model is a Bayesian model.
#' @slot df The data that will be used to fit the model.
#' @slot family The model family, currently only used in the brms backend.
#' @slot pcvrForm The formula that was specified in \link{growthSS} and used in other pcvr functions.
#' @slot type The model backend.
#' @slot model The name of the main growth formula.
#' @slot call The call to \link{growthSS} or \link{mvSS}.
#' @slot start Starting values for frequentist models.
#' @slot taus Quantiles for nlrq/rq models.
#'
#' @seealso
#'   \code{\link{growthSS}},
#'   \code{\link{mvSS}}
#'
NULL

pcvrss <- function(x) {
  class(x) <- "pcvrss"
  x
}

#' Print a \code{pcvrss} object.
#'
#' @aliases print.pcvrss
#'
#' @param x An object of class \code{pcvrss}
#'  to method \code{summary} of \code{pcvrss}.
#' @param ... further arguments, passed to print.default.
#'
#' @seealso \code{\link{summary.pcvrss}}
#' @method print pcvrss
#' @export
print.pcvrss <- function(x, ...) {
  print(summary.pcvrss(x), ...)
}


#' Summarize a \code{pcvrss} object.
#'
#' @aliases summary.pcvrss
#'
#' @param object An object of class \code{pcvrss}
#'  to method \code{summary} of \code{pcvrss}.
#' @param ... further arguments, passed to print.default.
#'
#' @method summary pcvrss
#' @export

summary.pcvrss <- function(object, ...) {
  out <- object[which(names(object) %in% c("type", "family", "model", "formula", "df", "pcvrForm"))]
  class(out) <- "pcvrsssummary"
  return(out)
}

#' Print a \code{pcvrsssummary} object.
#'
#' @aliases print.pcvrsssummary
#'
#' @param x An object of class \code{pcvrsssummary}.
#' @param ... further arguments, which are currently ignored.
#'
#' @seealso \code{\link{print.pcvrsssummary}}
#' @method print pcvrsssummary
#' @export
print.pcvrsssummary <- function(x, ...) {
  model_type <- gsub("int_", "(Intercept)", x$model)
  cat(paste(model_type,
            x$type,
            x$family,
            "model:\n"))
  cat("\npcvr formula variables:\n")
  yxig <- .parsePcvrForm(x$pcvrForm)[1:4]
  non_null <- !unlist(lapply(yxig, is.null))
  non_dummy <- !grepl("dummyIndividual|dummyGroup", yxig)
  yxig <- yxig[non_null & non_dummy]
  yxig_key <- c("Outcome:", "X:", "Individual:", "Group:")
  yxig_key <- yxig_key[non_null & non_dummy]
  cat(paste(yxig_key, yxig, collapse = "\n"))
  cat("\n\nModel Formula:\n")
  if (x$type == "nlme") {
    print(x$formula$model)
    fixed <- apply(do.call(rbind,
                           lapply(x$formula$fixed, as.character))[, c(2, 1, 3)],
                   1, paste, collapse = " ")
    cat(paste(fixed, collapse = "\n"))
  } else {
    print(x$formula)
  }
  cat("\nData:\n")
  print(x$df[1:3, !grepl("dummyIndividual|dummyGroup", colnames(x$df))])
  cat(paste0("...\n"))
  cat(paste0("(", nrow(x$df), " rows)"))
  invisible(x)
}
