#' reader function called within read.pcv when large data is used
#'
#' @param inputFile Path to csv file of plantCV output, should be provided internally in read.pcv
#' @param filters filtering conditions, see read.pcv for details.
#'  Format as list("trait in area, perimeter", "other in value")
#' @param reader reader argument passed from read.pcv, defaults to read.csv.
#' @param awk Optional awk command to use instead.
#' @param ... further arguments passed from \code{read.pcv}.
#' @import data.table
#' @importFrom utils read.csv
#' @return Returns a dataframe after subsetting happens outside of R using the awk statement from
#' \code{awkHelper}.
#' @return A data.frame
#'
#' @keywords internal
#' @noRd

pcv.sub.read <- function(inputFile, filters, reader = "read.csv", awk = NULL, ...) {
  awkCommand <- awkHelper(inputFile, filters, awk)
  COLS <- colnames(read.csv(inputFile, nrows = 1))
  if (reader == "fread") {
    x <- as.data.frame(data.table::fread(cmd = awkCommand, col.names = COLS, ...))
  } else {
    readingFunction <- match.fun(reader)
    x <- suppressMessages(as.data.frame(readingFunction(pipe(awkCommand), ...)))
    colnames(x) <- COLS
  }
  return(x)
}
