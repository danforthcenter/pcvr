#' Example model
#'
#' @description A logistic brms model of simulated growth data.
#'
#' @format A brmsfit object and a growthss object
#' @examples
#' \donttest{
#' data(fit, package = "pcvr")
#' summary(fit)
#' }
#'
"fit"

#' Example growthSS object used to make example models
#'
#' @description A growthSS class object describing either a growth model or survival model.
#'
#' @examples
#' \donttest{
#' data(fit, package = "pcvr")
#' summary(ss)
#' }
#'
"ss"


#' Example survival model
#'
#' @description An exponential survival brms model of simulated data.
#'
#' @format A brmsfit object and a growthss object
#' @examples
#' \donttest{
#' data(surv, package = "pcvr")
#' summary(surv)
#' }
#'
"surv"
