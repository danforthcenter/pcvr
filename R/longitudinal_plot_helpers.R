#' Helper to turn off faceting for dummy grouping
#' @keywords internal
#' @noRd

.no_dummy_labels <- function(group, facetGroups) {
  if (group == "dummyGroup") {
    facetGroups <- FALSE
  }
  return(facetGroups)
}
