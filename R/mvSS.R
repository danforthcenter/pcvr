#' Ease of use multi-value trait model helper function.
#' 
#' This function provides a simplified interface to modeling multi-value traits using \link{growthSS}.
#' Output from this should be passed to \link{fitGrowth} to fit the specified model.
#'
#' @param form A formula similar to \code{label | value ~ time + id/group} where label is a column
#' of histogram bins, value is the counts within those bins, time is an optional time variable,
#' id identifies an individual, and group contains the treatment groups.
#' @param sigma Distributional models passed to \link{growthSS}.
#' @param df Data passed to \link{growthSS}.
#' @param pars Parameters to vary, passed to \link{growthSS}.
#' @param start Starting values or priors, passed to \link{growthSS}.
#' @param type Backend to use, passed to \link{growthSS}.
#' @param tau Quantile to model, passed to \link{growthSS}.
#' @param hierarchy Formulae describing any hierarchical models, see \link{growthSS}.
#' @param spectral_index Optionally, a spectral index
#' \href{https://plantcv.readthedocs.io/en/stable/spectral_index/}{from those calculated by PlantCV}.
#' If this is given then the appropriate truncation and model family (if applicable)
#' will be included for the index you are using without you having to write it in the formula.
#' @keywords multi-value
#' @return A named list of plots showing prior distributions that \code{growthSS} would use,
#' optionally with a plot of simulated growth curves using draws from those priors.
#'
#' @examples
#' mv_df <- mvSim(dists = list(rnorm = list(mean = 100, sd = 30)), wide = FALSE)
#' mv_df <- mv_df[mv_df$value > 0, ]
#' # show that this has identical output to a more complicated call to growthSS
#'
#' @export

mvSS <- function(model, form, sigma = NULL, df, start = NULL,
                 pars = NULL, type = "brms", tau = 0.5, hierarchy = NULL,
                 spectral_index = c("none", "ari", "ci_rededge", "cri550", "cri700",
                                    "egi", "evi", "gdvi", "mari", "mcari", "mtci", "ndre",
                                    "ndvi", "pri", "psnd_chlorophyll_a", "psnd_chlorophyll_b",
                                    "psnd_caroteniods", "psri", "pssr_chlorophyll_a",
                                    "pssr_chlorophyll_b", "pssr_caroteniods", "rgri",
                                    "rvsi", "savi", "sipi", "sr", "vari", "vi_green", "wi",
                                    "fvfm", "fqfm")) {
  #* `get spectral index helper function`
  #* spectral index function should just return truncation and family, then there will be a separate
  #* function that applies those changes to the form argument.
  
  #* `run spectral index helper function`
  #* 
  
  #* `run formula cleaner`
  #* This should return the final, usable formula
  
  #* `if time is missing then..? do stuff here?`
  #* this is the tricky part I think. but on the other hand if time is missing then it's
  #* really only one option for what happens next, so maybe it's not bad. Model is basically
  #* ignored in this case.
  
  #* `if time is not missing, call growthSS with new model formula`
  
  out <- growthSS(model = model, form = form, sigma = sigma, df = df, start = start,
                  pars = pars, type = type, tau = tau, hierarchy = hierarchy)
  return(out)
}




