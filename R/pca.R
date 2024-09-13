#' Function to run a PCA, plot and optionally return the data with PCA coordinates and pca object
#'
#' @param df Dataframe to ordinate
#' @param cols columns to reduce dimensions of. Can be specified with names or positions.
#' If this is length of 1 then it is treated as regex pattern to match
#' the column names that should be used.
#' @param color column name(s) used to color points in the pca plot.
#' @param facet Optional column or vector to facet plots on.
#' @param returnData Logical, should data be returned?
#' @param ncp Optional, number of principal components to return attached
#'  to dataframe if data is returned. Defaults to all.
#' @keywords pca
#' @details If data is returned then it will contain the coordinates from the
#'  PCA and will not contain the columns that were reduced.
#'
#' @import ggplot2
#' @import FactoMineR
#' @importFrom stats as.formula
#' @return A ggplot or list with a ggplot, a dataframe with the data and PCs, and the factominer
#' PCA object as elements.
#' @examples
#'
#' dists <- list(
#'   rlnorm = list(meanlog = log(40), sdlog = 0.5),
#'   rnorm = list(mean = 60, sd = 10)
#' )
#' mv <- mvSim(
#'   dists = dists, n_samples = 100, counts = 1000,
#'   min_bin = 1, max_bin = 180, wide = TRUE
#' )
#' mv$otherGroup <- sample(c("a", "b"), size = nrow(mv), replace = TRUE)
#' pcadf(mv, cols = "sim_", returnData = TRUE)
#' pcadf(mv, cols = 2:181, color = c("group", "otherGroup"), returnData = FALSE)
#'
#' @export

pcadf <- function(df = NULL, cols = NULL, color = NULL,
                  facet = NULL, returnData = TRUE, ncp = NULL) {
  if (is.character(cols) && length(cols) == 1) {
    cols <- which(grepl(cols, colnames(df)))
  }
  if (!is.null(color) && length(color) > 1) {
    df[[paste(color, collapse = ".")]] <- interaction(df[, color])
    color <- paste(color, collapse = ".")
  }
  if (is.null(ncp)) {
    ncp <- min(dim(df[, cols])) - 1
  }
  pca <- FactoMineR::PCA(df[, cols], ncp = ncp, graph = FALSE)
  pc1Var <- round(pca$eig[1, 2], 3)
  pc2Var <- round(pca$eig[2, 2], 3)
  coords <- as.data.frame(pca$ind)
  coords <- coords[, grepl("coord", colnames(coords))]
  colnames(coords) <- gsub("coord.Dim.", "pc", colnames(coords))
  if (!is.numeric(cols)) {
    cols <- which(colnames(df) %in% cols)
  }
  pca.df <- cbind(as.data.frame(df[, -cols]), coords)
  colnames(pca.df)[1:(length(colnames(df)) - length(cols))] <- colnames(df)[-cols]
  facetLayer <- NULL
  if (!is.null(facet)) {
    facetLayer <- ggplot2::facet_wrap(as.formula(paste0("~", paste(facet, collapse = "+"))))
  }

  if (is.null(color)) {
    pca.df$dummyVariableForColor <- 1
    color <- "dummyVariableForColor"
  }

  plots <- .pcaGeneralPlot(pca.df, color, facetLayer, pc1Var, pc2Var)

  if (returnData) {
    return(list("data" = pca.df, "pca" = pca, "plot" = plots))
  } else {
    return(plots)
  }
}

#' general pca plotting
#' @keywords internal
#' @noRd

.pcaGeneralPlot <- function(pca.df, color, facetLayer, pc1Var, pc2Var) {
  plots <- ggplot2::ggplot(pca.df, ggplot2::aes(
    x = .data$pc1,
    y = .data$pc2, color = .data[[color]]
  )) +
    ggplot2::geom_point() +
    ggplot2::labs(x = paste0("PC 1 (", pc1Var, "%)"), y = paste0("PC 2 (", pc2Var, "%)")) +
    pcv_theme()
  if (color == "dummyVariableForColor") {
    plots <- plots + ggplot2::theme(legend.position = "none")
  }
  plots <- plots + facetLayer
  return(plots)
}
