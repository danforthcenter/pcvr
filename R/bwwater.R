#' Read in lemnatech watering data from metadata.json files
#'
#' @param file Path to a json file of lemnatech metadata.
#' @param envKey Character string representing the json key for environment data.
#'  By default this is set to "environment".
#'  Currently there are no situations where this makes sense to change.
#' @keywords watering json
#' @import jsonlite
#' @importFrom utils type.convert
#' @return A data frame containing the bellwether watering data
#' @examples
#' tryCatch(
#'   {
#'     w <- bw.water("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/metadata.json")
#'   },
#'   error = function(e) {
#'     message(e)
#'   }
#' )
#' @export

bw.water <- function(file = NULL, envKey = "environment") {
  meta <- jsonlite::fromJSON(txt = file)
  env <- as.data.frame(do.call(rbind, meta[[envKey]]))
  env$snapshot <- rownames(env)
  rownames(env) <- NULL
  env <- as.data.frame(apply(env, 2, as.character))
  env <- type.convert(env, as.is = TRUE)
  if ("timestamp" %in% colnames(env)) {
    tryCatch(
      {
        env$timestamp <- as.POSIXct(env$timestamp, tryFormats = c(
          "%Y-%m-%d %H:%M:%OS",
          "%Y-%m-%dT%H:%M:%OS",
          "%Y/%m/%d %H:%M:%OS",
          "%Y-%m-%d %H:%M",
          "%Y/%m/%d %H:%M",
          "%Y-%m-%d",
          "%Y/%m/%d"
        ), tz = "UTC")
        begin <- min(env$timestamp, na.rm = TRUE)
        message(paste0(
          "Using the first watering time, ", begin,
          ", as beginning of experiment to assign DAS"
        ))
        env$DAS <- as.numeric((env$timestamp - begin) / 24 / 60 / 60)
      },
      error = function(err) {},
      warning = function(warn) {}
    )
  }
  return(env)
}
