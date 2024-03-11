#' Function to parse survival model specifications in growthSS
#'
#' @param model the model specified for growthSS
#'
#' @return A list. The first component is logical whether or not the formula is a survival model.
#' Second component is the distribution to use for survival modeling (this is only used if type=="brms")
#'
#' @examples
#'
#' .survModelParser("survival - binomial")
#' .survModelParser("survival")
#' .survModelParser("logistic")
#' .survModelParser("logistic+linear")
#'
#' @keywords internal
#' @noRd

.survModelParser <- function(model) {
  distributions <- c(
    "binomial", "gengamma", "gengamma.orig", "genf", "genf.orig",
    "weibull", "gamma", "exp", "llogis", "lnorm", "gompertz",
    "exponential", "lognormal"
  )
  if (grepl("survival", model)) {
    survival <- TRUE
    model <- trimws(gsub("survival", "", model))
    if (nchar(model) == 0) {
      model <- "weibull"
    }
    dist <- match.arg(model, distributions)
    return(list("survival" = survival, "model" = dist))
  } else {
    return(list("survival" = FALSE, "model" = model))
  }
}
