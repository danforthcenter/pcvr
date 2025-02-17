#' Helper function to check groups in data.
#'
#' @param df Data frame to use.
#' @param group Set of variables to use in grouping observations.
#' These taken together should identify a unique plant (or unique plant at a unique angle) across time.
#' @return If there are duplicates in the grouping then this will return a message with code to start
#' checking the duplicates in your data.
#'
#' @examples
#'
#' df <- growthSim("linear",
#'   n = 10, t = 10,
#'   params = list("A" = c(2, 1.5))
#' )
#' checkGroups(df, c("time", "id", "group"))
#' df$time[12] <- 3
#' checkGroups(df, c("time", "id", "group"))
#'
#' @export

checkGroups <- function(df, group) {
  tab <- table(interaction(df[, c(group)]))
  if (any(tab > 1)) {
    dataname <- deparse(substitute(df))
    nms <- names(tab)[which(as.numeric(tab) > 1)]
    dupString <- paste0(
      dataname,
      "[duplicated(interaction(",
      paste(paste0(dataname, "$", c(group)), collapse = ", "),
      ")),]"
    )
    firstDup <- paste0(
      dataname, "[interaction(",
      paste(paste0(
        dataname, "$",
        c(group)
      ), collapse = ", "), ")=='",
      nms[1], "',]"
    )
    eval(parse(text = dupString))
    w <- paste0(
      "There are ", length(nms), " observations that are not uniquely identified.",
      "\nThe max number of duplicates is ",
      max(tab, na.rm = TRUE), ".\nRun `", dupString, "` to see the duplicated rows,\n",
      " or ", firstDup, " to see the first duplicated instance."
    )
    message(w)
  } else {
    message("Grouping is unique")
  }
  invisible(NULL)
}
