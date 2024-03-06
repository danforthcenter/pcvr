#' Functions to prepare several families of distributional brms models
#'
#' @param model a model passed from growthSS
#'
#' @return A list of elements to pass brmSS for fitting distributional models
#'
#' @examples
#' .brmFamilyHelper("logistic")
#' .brmFamilyHelper("poisson: logistic")
#' .brmFamilyHelper("von_mises: logistic")
#'
#' @keywords internal
#' @noRd

.brmFamilyHelper <- function(model) {
  if (!grepl("[:]", model)) {
    model <- paste0("student:", model)
  }
  dist <- trimws(gsub(":.*", "", model))
  rhs <- trimws(gsub(".*:", "", model))
  family <- dist
  dpars <- brms::brmsfamily(dist)$dpars[-1]
  return(list(family = family, dpars = dpars, rhs = rhs))
}
