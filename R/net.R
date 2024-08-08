#' Network analysis of a distance matrix
#'
#' @description Easy igraph use with pcv.emd output
#'
#'
#' @param emd A long dataframe as returned by pcv.emd.
#' Currently this function is only made to work with dataframe output,
#' not distance matrix output.
#' @param meta Metadata to be carried from pcv.emd output into the network,
#' defaults to NULL which will use all metadata.
#' Type conversion will be attempted for these columns.
#' @param dissim Logical, should the distCol be inverted to make a dissimilarity value?
#' @param distCol The name of the column containing distances/dissimilarities.
#' Defaults to "emd" for compatability with pcv.emd
#' @param filter This can be either a numeric (0.5) in which case it is taken as
#' a filter where only edges with values greater than or equal to that number are
#' kept or a character string ("0.5") in which case the strongest X percentage of edges are kept.
#' This defaults to 0.5 which does some filtering, although that should not be considered
#' the best behavior for every setting. If this is NULL then your network will be
#' almost always be a single blob, if set too high there will be very few nodes.
#' Note that this filtering happens after converting to dissimilarity if dissim=TRUE.
#' @param direction Direction of filtering, can be either "greater" or "lesser".
#' @import ggplot2
#' @import igraph
#' @importFrom utils type.convert
#' @importFrom stats quantile
#'
#' @keywords emd multi-value-trait network
#' @examples
#'
#' ## Not run:
#'
#' library(extraDistr)
#' dists <- list(
#'   rmixnorm = list(mean = c(70, 150), sd = c(15, 5), alpha = c(0.3, 0.7)),
#'   rnorm = list(mean = 90, sd = 3)
#' )
#' x <- mvSim(dists = dists, n_samples = 5, counts = 1000,
#'            min_bin = 1, max_bin = 180, wide = TRUE)
#' emd_df <- pcv.emd(x,
#'                   cols = "sim", reorder = c("group"), mat = FALSE,
#'                   plot = FALSE, parallel = 1
#' )
#' net <- pcv.net(emd_df, meta = "group")
#' net2 <- pcv.net(emd_df, meta = "group", filter = "0.9", direction = "lesser")
#'
#' ## End(Not run)
#'
#' @return Returns a list containing three elements:
#' \code{nodes}: A dataframe of node data.
#' \code{edges}: A dataframe of edges between nodes.
#' \code{graph}: The network as an igraph object
#'
#' @export
#'

pcv.net <- function(emd = NULL, meta = NULL, dissim = TRUE, distCol = "emd", filter = 0.5,
                    direction = "greater") {
  if (is.data.frame(emd)) {
    #* convert to dissimilarity if metric is similarity
    if (dissim) {
      emd[[distCol]] <- 1 / emd[[distCol]]
      emd[[distCol]] <- ifelse(is.infinite(emd[[distCol]]) | is.na(emd[[distCol]]), 0, emd[[distCol]])
    }
    #* filter for edge strength
    if (!is.null(filter)) {
      if (is.character(filter)) {
        filter <- quantile(emd[[distCol]], probs = as.numeric(filter))
      }
      if (match.arg(direction, c("greater", "lesser")) == "greater") {
        emd <- emd[emd[[distCol]] > filter, ]
      } else {
        emd <- emd[emd[[distCol]] < filter, ]
      }
    }
    #* turn long data into a graph and extract nodes/edges
    g <- igraph::graph_from_data_frame(emd, directed = FALSE)
  } else {
    stop("emd must be a dataframe.")
  }
  if (is.null(meta)) {
    meta <- unique(sub("_i$|_j$", "", colnames(emd)[grepl("_i$|_j$", colnames(emd))]))
  }

  gg <- as.data.frame(igraph::layout_nicely(g))
  both <- igraph::as_data_frame(g, "both")
  gg$index <- both$vertices$name
  eg <- both$edges
  #* link metadata to nodes
  metaIndex <- lapply(meta, function(m) which(grepl(m, colnames(eg))))
  newCols <- (ncol(gg) + 1):(ncol(gg) + length(meta))
  gg[, newCols] <- lapply(metaIndex, function(m) {
    i <- m[[1]]
    j <- m[[2]]
    f <- eg[[i]][match(gg$index, eg$from)]
    to <- eg[[j]][match(gg$index, eg$to)]
    to[which(is.na(to))] <- f[which(is.na(to))]
    to
  }) # this can be NA if there is no 'from' edge connected to a node, so check 'to' edges as well.
  colnames(gg)[newCols] <- meta
  gg[, newCols] <- type.convert(gg[, newCols], as.is = TRUE)

  #* Calculate network metrics
  gg$betweenness <- igraph::betweenness(g)
  gg$degree <- igraph::degree(g)
  igraph::E(g)$weight <- eg[[distCol]] + 0.1
  gg$strength <- igraph::strength(g)
  gg$harmonic_centrality <- igraph::harmonic_centrality(g)
  gg$eigen_centrality <- igraph::eigen_centrality(g)[[1]]
  gg$authority_score <- igraph::authority_score(g)[[1]]
  gg$page_rank <- igraph::page_rank(g)[[1]]
  #* add coordinates for plotting edges
  eg$from.x <- gg$V1[match(eg$from, gg$index)]
  eg$from.y <- gg$V2[match(eg$from, gg$index)]
  eg$to.x <- gg$V1[match(eg$to, gg$index)]
  eg$to.y <- gg$V2[match(eg$to, gg$index)]
  return(list("nodes" = gg, "edges" = eg, "graph" = g))
}
