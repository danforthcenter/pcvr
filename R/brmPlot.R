#' Function to visualize brms models similar to those made using growthSS outputs.
#'
#' Models fit using \link{growthSS} inputs by \link{fitGrowth} (and similar models made through other
#' means) can be visualized easily using this function. This will generally be called by
#' \code{growthPlot}.
#'
#' @param fit A brmsfit object, similar to those fit with \code{\link{growthSS}} outputs.
#' @param form A formula similar to that in \code{growthSS} inputs specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' @param df An optional dataframe to use in plotting observed growth curves on top of the model.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the available data has not reached some point (such as asymptotic size),
#' although prediction using splines outside of the observed range is not necessarily reliable.
#' @param facetGroups logical, should groups be separated in facets? Defaults to TRUE.
#' @param groupFill logical, should groups have different colors? Defaults to FALSE.
#' If TRUE then viridis colormaps are used in the order
#' of virMaps
#' @param virMaps order of viridis maps to use. Will be recycled to necessary length.
#' Defaults to "plasma", but will generally be informed by growthPlot's default.
#' @keywords growth-curve, logistic, gompertz, monomolecular, linear, exponential, power-law
#' @import ggplot2
#' @import viridis
#' @importFrom stats as.formula
#' @examples
#'
#' ## Not run:
#'
#' if (FALSE) {
#'   data(bw_vignette_fit)
#'   brmPlot(bw_vignette_fit, y ~ time | id / group, df = NULL)
#'   print(load(url("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/brmsFits.rdata")))
#'   brmPlot(fit_25, form = y ~ time | id / group)
#'   brmPlot(fit_9, form = y ~ time | id / group)
#'   brmPlot(fit_15, form = y ~ time | id / group)
#' }
#'
#' ## End(Not run)
#'
#' @return Returns a ggplot showing a brms model's credible
#' intervals and optionally the individual growth lines.
#'
#' @export

brmPlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL, facetGroups = TRUE,
                    groupFill = FALSE, virMaps = c("plasma")) {
  fitData <- fit$data
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  if (individual == "dummyIndividual") {
    individual <- NULL
  }
  group <- parsed_form$group
  df <- parsed_form$data

  #* if no group in fitdata then add a dummy group
  if (!group %in% colnames(fitData)) {
    fitData[[group]] <- "a"
  }

  probs <- seq(from = 99, to = 1, by = -2) / 100

  if (is.null(timeRange)) {
    timeRange <- unique(fitData[[x]])
  }

  newData <- data.frame(
    x = rep(timeRange, times = length(unique(fitData[[group]]))),
    group = rep(unique(fitData[[group]]), each = length(timeRange)),
    individual = rep(paste0("new_", seq_along(unique(fitData[[group]]))), each = length(timeRange))
  )
  colnames(newData) <- c(x, group, individual)
  predictions <- cbind(newData, predict(fit, newData, probs = probs))

  if (!is.null(groups)) {
    predictions <- predictions[predictions$group %in% groups, ]
    if (!is.null(df)) {
      df <- df[df[[group]] %in% groups, ]
    }
  }

  #* `facetGroups`
  if (facetGroups) {
    if (!all(fitData[[group]] == "a")) {
      facetLayer <- ggplot2::facet_wrap(as.formula(paste0("~", group)))
    } else {
      facetLayer <- NULL
    }
  } else {
    facetLayer <- NULL
  }
  #* `groupFill`
  if (groupFill) {
    virList <- lapply(rep(virMaps, length.out = length(unique(df[[group]]))), function(pal) {
      viridis::viridis(n = length(probs), option = pal)
    })
  } else {
    pal <- viridis::plasma(n = length(probs))
    virList <- lapply(seq_along(unique(df[[group]])), function(i) {
      pal
    })
  }


  p <- ggplot2::ggplot(predictions, ggplot2::aes(x = .data[[x]], y = .data$Estimate)) +
    facetLayer +
    ggplot2::labs(x = x, y = y) +
    pcv_theme()

  for (g in seq_along(unique(predictions[[group]]))) {
    iteration_group <- unique(predictions[[group]])[g]
    sub <- predictions[predictions[[group]] == iteration_group, ]
    p <- p +
      lapply(seq(1, 49, 2), function(i) {
        ggplot2::geom_ribbon(
          data = sub, ggplot2::aes(
            ymin = .data[[paste0("Q", i)]],
            ymax = .data[[paste0("Q", 100 - i)]]
          ),
          fill = virList[[g]][i], alpha = 0.5
        )
      })
  }

  if (!is.null(df) && !is.null(individual)) {
    p <- p + ggplot2::geom_line(
      data = df, ggplot2::aes(.data[[x]], .data[[y]],
        group = interaction(.data[[individual]], .data[[group]])
      ),
      color = "gray20", linewidth = 0.2
    )
  }
  return(p)
}
