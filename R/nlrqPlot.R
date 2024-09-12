#' Function to visualize common \code{quantreg::nlrq} growth models.
#'
#' Models fit using \link{growthSS} inputs by \link{fitGrowth}
#' (and similar models made through other means)
#' can be visualized easily using this function. This will generally be called by \code{growthPlot}.
#'
#' @param fit A model fit, or list of model fits, returned by \code{fitGrowth} with type="nlrq".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm} part of the
#' output) specifying the outcome, predictor, and grouping structure of the data as
#' \code{outcome ~ predictor|individual/group}. If the individual and group are specified then the
#' observed growth lines are plotted.
#' @param df A dataframe to use in plotting observed growth curves on top of the model.
#' This must be supplied for nlrq models.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the avaiable data has not reached some point (such as asymptotic size).
#' @param facetGroups logical, should groups be separated in facets? Defaults to TRUE.
#' @param groupFill logical, should groups have different colors? Defaults to FALSE.
#' If TRUE then viridis colormaps are used in the order of virMaps
#' @param virMaps order of viridis maps to use. Will be recycled to necessary length.
#' Defaults to "plasma", but will generally be informed by growthPlot's default.
#' @keywords growth-curve
#' @importFrom methods is
#' @import ggplot2
#' @importFrom stats setNames predict
#' @importFrom viridis plasma
#' @examples
#'
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group,
#'   tau = c(0.5, 0.9), df = simdf, start = NULL, type = "nlrq"
#' )
#' fit <- fitGrowth(ss)
#' nlrqPlot(fit, form = ss$pcvrForm, df = ss$df, groups = "a", timeRange = 1:20)
#' nlrqPlot(fit, form = ss$pcvrForm, df = ss$df, groupFill = TRUE, virMaps = c("plasma", "viridis"))
#'
#' @return Returns a ggplot showing an nlrq model's quantiles
#'  and optionally the individual growth lines.
#'
#' @export

nlrqPlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL,
                     facetGroups = TRUE, groupFill = FALSE, virMaps = c("plasma")) {
  #* `get needed information from formula`
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  if (individual == "dummyIndividual") {
    individual <- NULL
  }
  group <- parsed_form$group
  df <- parsed_form$data
  #* `filter by groups if groups != NULL`
  if (!is.null(groups)) {
    df <- df[df[[group]] %in% groups, ]
  }
  #* `make new data if timerange is not NULL`
  if (!is.null(timeRange)) {
    new_data <- do.call(rbind, lapply(unique(df[[group]]), function(g) {
      stats::setNames(data.frame(g, timeRange), c(group, x))
    }))
  } else {
    new_data <- NULL
  }
  #* `standardize fit class`
  if (methods::is(fit, "nlrq")) {
    fit <- list(fit)
    names(fit) <- fit[[1]]$m$tau()
  }
  #* `add predictions`
  preds <- do.call(cbind, lapply(fit, function(f) {
    tau <- f$m$tau()
    stats::setNames(data.frame(stats::predict(f, newdata = new_data)), paste0("Q_", tau))
  }))
  predCols <- colnames(preds)
  keep <- which(!duplicated(preds))
  plotdf <- cbind(df[keep, ], preds[keep, ])
  colnames(plotdf) <- c(colnames(df), colnames(preds))
  
  #* `facetGroups`
  facet_layer <- NULL
  if (facetGroups) {
    facet_layer <- ggplot2::facet_wrap(stats::as.formula(paste0("~", group)))
  }
  #* `groupFill`
  if (groupFill) {
    virList <- lapply(rep(virMaps, length.out = length(unique(df[[group]]))), function(pal) {
      virpal_p1 <- viridis::viridis(ceiling(length(predCols) / 2), direction = 1, end = 1, option = pal)
      virpal_p2 <- viridis::viridis(ceiling(length(predCols) / 2),
                                    direction = -1, end = 1, option = pal)[-1]
      c(virpal_p1, virpal_p2)
    })
  } else {
    virpal_p1 <- viridis::plasma(ceiling(length(predCols) / 2), direction = 1, end = 1)
    virpal_p2 <- viridis::plasma(ceiling(length(predCols) / 2), direction = -1, end = 1)[-1]
    virpal <- c(virpal_p1, virpal_p2)
    virList <- lapply(seq_along(unique(df[[group]])), function(i) {
      virpal
    })
  }
  #* `layer for individual lines if formula was complete`
  if (!is.null(individual)) {
    individual_lines <- ggplot2::geom_line(
      data = df, ggplot2::aes(
        x = .data[[x]], y = .data[[y]],
        group = interaction(
          .data[[individual]],
          .data[[group]]
        )
      ),
      linewidth = 0.25, color = "gray40"
    )
  } else {
    individual_lines <- list()
  }
  #* `plot`
  plot <- ggplot(plotdf, ggplot2::aes(group = interaction(.data[[group]]))) +
    facet_layer +
    individual_lines +
    labs(x = x, y = as.character(form)[2]) +
    pcv_theme()
  
  for (g in seq_along(unique(plotdf[[group]]))) {
    iteration_group <- unique(plotdf[[group]])[g]
    sub <- plotdf[plotdf[[group]] == iteration_group, ]
    plot <- plot +
      lapply(seq_along(predCols), function(i) {
        ggplot2::geom_line(data = sub, ggplot2::aes(x = .data[[x]], y = .data[[predCols[i]]]),
                           color = virList[[g]][i], linewidth = 0.7)
      })
  }
  
  return(plot)
}
