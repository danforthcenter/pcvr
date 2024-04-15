#' helper function to help make growth models including an intercept term in growthSS
#' @examples
#' .intModelHelper("int_logistic")
#' .intModelHelper("logistic")
#'
#' @keywords internal
#' @noRd

.intModelHelper <- function(model) {
  if (grepl("^int_", model)) {
    int <- TRUE
    model <- gsub("^int_", "", model)
  } else {
    int <- FALSE
  }
  return(list("model" = model, "int" = int))
}
