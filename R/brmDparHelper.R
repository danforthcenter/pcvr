#' Helper function to parse distributional formulas and return components for fitting a brms model
#'
#' @param model a model passed from growthSS
#' @param x the x variable
#' @param y the response variable
#' @param group the grouping variable
#' @param sigma Logical, is this a distributional parameter? Defaults to TRUE here.
#' @param nTimes Number of unique x values. Used for splines, if nTimes is too low then the spline knots
#' must be adjusted. Defaults to 25.
#' @param useGroup Logical, should groups be used?
#' @param priors Priors in growthSS syntax, passed to .brmsChangePointHelper for thresholded models.
#' @param int logical, should an intercept be included in the model?
#'
#' @return A list of elements to pass brmSS for fitting distributional models
#'
#'
#' @keywords internal
#' @noRd

.brmDparHelper <- function(dpar, model, x, group, nTimes, useGroup, priors, int = FALSE,
                           force_nl = FALSE) {
  splineDparHelperForm <- NULL
  if (grepl("\\+", model)) {
    chngptHelperList <- .brmsChangePointHelper(model, x,
      y = dpar, group, dpar = TRUE,
      nTimes, useGroup, priors, int = int
    )
    dparForm <- chngptHelperList$growthForm
    dpar_pars <- chngptHelperList$pars
    splineDparHelperForm <- chngptHelperList$splineHelperForm
  } else {
    if (model == "homo") {
      model <- "int" # recode alternate names
    } else if (model == "spline") {
      model <- "gam"
    }
    stringBrmsDparFormFun <- paste0(".brms_form_", gsub(" ", "", model))
    brmsDparFormFun <- match.fun(stringBrmsDparFormFun)
    formResDpar <- brmsDparFormFun(x, dpar, group,
      dpar = TRUE, nTimes = nTimes,
      useGroup = useGroup, prior = priors, int = int, force_nl = force_nl
    )
    dparForm <- formResDpar$form
    dpar_pars <- formResDpar$pars
  }
  return(list(
    "dparForm" = dparForm, "dpar_pars" = dpar_pars,
    "dparSplineHelperForm" = splineDparHelperForm
  ))
}
