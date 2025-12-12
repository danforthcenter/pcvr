#' Show a model fit with pcvr in a ggplot layer
#'
#' @description
#' Add a model fit with \link{growthSS} and \link{fitGrowth} to a ggplot object.
#' The exact geom used depends on the model class (see details).
#'
#' @param mapping Set of aesthetic mappings created by \code{ggplot2::aes()}.
#' If specified and ‘inherit.aes = TRUE’ (the default),
#' it is combined with the default mapping at the top level of the plot.
#' If there is no mapping then it is filled in by default using the \code{pcvrss} object.
#' @param data The data to be displayed in this layer.
#' This behaves per normal ggplot2 expectations except
#' that if data is missing (ie, not inherited or specified) then the data from \code{ss} is used.
#' @param fit A model object returned from \code{fitGrowth}.
#' @param ss A \code{pcvrss} object. Only the "pcvrForm" and "df" elements are used.
#' @param inherit.aes Logical, should aesthetics be inherited from top level? Defaults to TRUE.
#' @param CI A vector of credible intervals to plot, defaults to 0.95.
#' @param ... Additional arguments passed to the ggplot layer.
#'
#' @details
#' These layers will behave largely like output from \code{\link{growthPlot}}, although \code{growthPlot}
#' has more arguments that directly control the plot since this stat only makes one layer.
#' The geometries used for each type of model are:
#' \itemize{
#'   \item{\strong{brms}: \code{geom_ribbon} for longitudinal plots, \code{geom_rect} for others.}
#'   \item{\strong{nlrq}: \code{geom_line}, replicated per each quantile.}
#'   \item{\strong{nlme}: \code{geom_smooth}, with ribbon based on the heteroskedastic term.}
#'   \item{\strong{nls}: \code{geom_line}, replicated per each quantile.}
#'   \item{\strong{nlrq}: \code{geom_line}, replicated per each quantile.}
#' }
#'
#' @examples
#' library(ggplot2)
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group,
#'   df = simdf, start = NULL, type = "nls"
#' )
#' fit <- fitGrowth(ss)
#' ggplot() +
#'   stat_growthss(fit = fit, ss = ss)
#'
#' @rdname stat_growthss
#' @seealso \link{growthPlot} for a self-contained plotting function
#' @keywords ggplot
#' @export

stat_growthss <- function(mapping = NULL, data = NULL,
                          fit = NULL, ss = NULL,
                          inherit.aes = TRUE, ...) {
  layer <- NULL
  model_class <- class(fit)[1]
  if (methods::is(fit, "list")) {
    model_class <- class(fit[[1]])
  }
  if (methods::is(fit, "brmsfit")) {
    is_gam <- as.logical(length(fit$basis$dpars[[fit$family$dpars[1]]]$sm))
    if (attr(fit$formula$formula, "nl") || is_gam) {
      layer <- stat_brms_model(
        mapping = mapping, data = data,
        fit = fit, ss = ss,
        inherit.aes = inherit.aes, ...
      ) # use brms stat
    } else { # linear models are survival models
      stop("stat_growthss is not implemented for survival models.")
    }
  } else {
    stat_model_function <- match.fun(paste0("stat_", model_class, "_model"))
    layer <- stat_model_function(
      mapping = mapping, data = data,
      fit = fit, ss = ss,
      inherit.aes = inherit.aes, ...
    )
  }
  return(layer)
}

#' @keywords internal
#' @noRd

stat_survreg_model <- function(...) {
  stop("stat_growthss is not implemented for survival models. Use growthPlot instead.")
}

#' @keywords internal
#' @noRd

stat_flexsurvreg_model <- stat_survreg_model
