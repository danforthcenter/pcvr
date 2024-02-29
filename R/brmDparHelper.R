#' Helper function to parse distributional formulas and return components for fitting a brms model
#' 
#' @param model a model passed from growthSS
#' @param x the x variable
#' @param y the response variable
#' @param group the grouping variable
#' @param sigma Logical, is this a distributional parameter? Defaults to TRUE here.
#' @param nTimes Number of unique x values. Used for splines, if nTimes is too low then the spline knots must be adjusted. Defaults to 25.
#' @param useGroup Logical, should groups be used?
#' @param priors Priors in growthSS syntax, passed to .brmsChangePointHelper for thresholded models.
#' 
#' @return A list of elements to pass brmSS for fitting distributional models
#' 
#' @examples 
#' .brmDparHelper("logistic")
#' .brmDparHelper("poisson: logistic")
#' .brmDparHelper("von_mises: logistic")
#' 
#' @keywords internal
#' @noRd

.brmDparHelper <- function(dpar, model, x, group, nTimes, useGroup, priors){
  splineDparHelperForm <- NULL
  if(grepl("\\+", model)){
    chngptHelper_list <- .brmsChangePointHelper(model, x, y = dpar, group, dpar = TRUE, nTimes, useGroup, priors)
    dparForm <- chngptHelper_list$growthForm
    dpar_pars <- chngptHelper_list$pars
    splineDparHelperForm <- chngptHelper_list$splineHelperForm
    
  } else{
    if(model=="homo"){model="int" # recode alternate names
    } else if(model =="spline"){model="gam"}
    string_brmsDparFormFun <- paste0(".brms_form_", gsub(" ", "", model))
    brmsDparFormFun <- match.fun(string_brmsDparFormFun)
    formResDpar <- brmsDparFormFun(x, dpar, group, dpar=TRUE, nTimes=nTimes, useGroup=useGroup, prior=priors)
    dparForm <- formResDpar$form
    dpar_pars <- formResDpar$pars
  }
  return(list("dparForm" = dparForm, "dpar_pars" = dpar_pars, "dparSplineHelperForm"=splineDparHelperForm))
}


