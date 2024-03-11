#' @rdname pcv.emd
#' @export
#'
pcv.euc <- function(df, cols = NULL, reorder = NULL, include = reorder, mat = FALSE, plot = TRUE,
                    parallel = getOption("mc.cores", 1), trait = "trait", id = "image",
                    value = "value", raiseError = TRUE, method = "euc") {
  pcv.emd(df, cols, reorder, include, mat, plot, parallel, trait, id, value, raiseError, method)
}
