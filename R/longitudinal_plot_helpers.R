#' Helper to turn off faceting for dummy grouping
#' @keywords internal
#' @noRd

.no_dummy_labels <- function(group, facetGroups) {
  if (all(group == "dummyGroup")) {
    facetGroups <- FALSE
  }
  return(facetGroups)
}

#' Helper for plotting changepoint models
#' @keywords internal
#' @noRd

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}
