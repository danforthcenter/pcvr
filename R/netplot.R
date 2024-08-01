#' Visualizing igraph networks
#'
#' @description Easy igraph visualization with pcv.net output
#'
#'
#' @param net Network object similar to that returned from pcv.net, having dataframes named "edges"
#' and "nodes"
#' @param fill Variable name(s) from the nodes data to be used to color points. By default "strength"
#' is used.
#' @param shape Optional discrete variable name(s) from the nodes data to be used to change the shape
#' of points. If this variable is numeric it will be coerced to character.
#' @param size Size of points, defaults to 3.
#' @param edgeWeight Edge dataframe column to weight connections between nodes. Defaults to "emd"
#' for compatability with \code{pcv.emd}.
#' @param edgeFilter How should edges be filtered? This can be either a numeric (0.5)
#' in which case it is taken as a filter where only edges with values greater than or equal to
#' that number are kept or a character string ("0.5") in which case the strongest
#' X percentage of edges are kept. This defaults to NULL which does no filtering,
#' although that should not be considered the best standard behaviour. See details.
#' @import ggplot2
#' @importFrom stats quantile
#'
#' @keywords emd, earth mover's distance, multi-value trait, network
#' @examples
#'
#' ## Not run:
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
#' net.plot(net)
#' net.plot(net, edgeFilter = "0.25")
#' net.plot(net, edgeFilter = 0.25, fill = c("degree", "group"),
#'          shape = c("degree", "group"))
#' net.plot(net, edgeFilter = 0.25, fill = c("degree", "group"),
#'          shape = c("degree"))
#' ## End(Not run)
#'
#' @return Returns a ggplot of a network.
#'
#' @export
#'

net.plot <- function(net, fill = "strength", shape = NULL, size = 3, edgeWeight = "emd",
                     edgeFilter = NULL) {
  nodes <- net[["nodes"]]
  edges <- net[["edges"]]
  if (is.null(fill)) {
    fill <- "NOFILL"
    nodes$NOFILL <- "a"
  }
  if (length(fill) > 1) {
    nodes$FILL <- interaction(nodes[, fill])
    fill <- "FILL"
  }
  if (is.null(shape)) {
    shape <- "NOSHAPE"
    nodes$NOSHAPE <- "a"
  }
  if (length(shape) > 1) {
    nodes$SHAPE <- interaction(nodes[, shape])
    shape <- "SHAPE"
  }
  if (is.numeric(nodes[[shape]])) {
    nodes[[shape]] <- as.character(nodes[[shape]])
  }
  if (!is.null(edgeFilter)) {
    if (is.character(edgeFilter)) {
      cutoff <- quantile(edges[[edgeWeight]], probs = as.numeric(edgeFilter))
      edges <- edges[edges[[edgeWeight]] >= as.numeric(cutoff), ]
    } else if (is.numeric(edgeFilter)) {
      edges <- edges[edges[[edgeWeight]] >= edgeFilter, ]
    } else {
      stop("edgeFilter must be character or numeric, see ?net.plot for details.")
    }
    nodes <- nodes[nodes$index %in% c(edges$from, edges$to), ]
  }


  p <- ggplot2::ggplot(nodes) +
    ggplot2::geom_segment(data = edges, ggplot2::aes(
      x = .data$from.x, xend = .data$to.x, y = .data$from.y, yend = .data$to.y,
      linewidth = .data[[edgeWeight]]
    ), colour = "black", alpha = 0.1) +
    ggplot2::geom_point(data = nodes, size = size, ggplot2::aes(
      x = .data$V1, y = .data$V2,
      fill = .data[[fill]], color = .data[[fill]],
      shape = .data[[shape]]
    ), alpha = 1, show.legend = TRUE) +
    ggplot2::scale_linewidth(range = c(0.1, 1.5)) +
    #* note that scaling shape should work, but there is a documented ggplot2
    #* bug where this messes up the legend, so
    #* until that is fixed I will not specify fillable shapes.
    ggplot2::guides(linewidth = "none", shape = ggplot2::guide_legend(nrow = 1), fill = "none") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "bottom")
  if (length(fill == 1) && fill == "NOFILL") {
    p <- p + ggplot2::guides(color = "none")
  }
  if (length(shape) == 1 && shape == "NOSHAPE") {
    p <- p + ggplot2::guides(shape = "none")
  }
  return(p)
}
