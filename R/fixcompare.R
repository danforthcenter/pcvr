#' Helper function for making comparisons between groups in various functions.
#'
#' @param compare Groups to compare. By default this is set to FALSE, which corresponds to no testing.
#' Other values of compare are passed to fixCompare to make t.test comparisons using ggpubr.
#' In short, NULL will run all pairwise T tests, a single value of the X axis variable will compare that
#' level to all other levels of the X variable, alternatively this can be a list as used by
#' ggpubr: list(c("level1", "level2"), c("level1", "level3"))
#' @param dat Dataframe to use
#' @param col Column as a character string to find values of `compare` in.
#' @param likeToLike Logical, should only similar groups be compared
#'   (currently defined as part 1 and 2 of group, split on [.])
#' @keywords internal
#' @importFrom utils combn
#' @return A list of comparisons as used in the "comparisons" argument of
#' \code{ggpubr::stat_compare_means}.
#' @export
#'
fixCompare <- function(compare, dat, col, likeToLike = FALSE) {
  if (!is.null(compare)) {
    if (length(compare) == 2 && !is.list(compare)) {
      compare <- list(compare)
    }
    if (length(compare) == 1 && !is.list(compare)) {
      if (!is.factor(dat[[col]])) {
        dat[[col]] <- factor(dat[[col]])
      }
      compare <- lapply(levels(dat[[col]])[!levels(dat[[col]]) %in% compare], function(l) c(compare, l))
    }
  } else if (is.null(compare)) {
    compare <- combn(unique(dat[[col]]), 2, simplify = FALSE)
  }
  if (likeToLike) {
    index <- unlist(lapply(compare, function(comp) {
      c1a <- strsplit(as.character(comp[1]), "[.]")[[1]][1]
      c1b <- strsplit(as.character(comp[1]), "[.]")[[1]][2]
      c2a <- strsplit(as.character(comp[2]), "[.]")[[1]][1]
      c2b <- strsplit(as.character(comp[2]), "[.]")[[1]][2]
      return(any(c(c1a, c1b) %in% c(c2a, c2b)))
    }))
    compare <- compare[index]
  }

  return(compare)
}
