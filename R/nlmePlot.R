#' Function to visualize common \code{nlme::nlme} growth models.
#'
#' Models fit using \link{growthSS} inputs by \link{fitGrowth}
#' (and similar models made through other means)
#'  can be visualized easily using this function. This will generally be called by \code{growthPlot}.
#'
#' @param fit A model fit returned by \code{fitGrowth} with type="nlme".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm} part of the
#' output) specifying the outcome, predictor, and grouping structure of the data as
#' \code{outcome ~ predictor|individual/group}
#' @param df A dataframe to use in plotting observed growth curves on top of the model.
#' This must be supplied for nlme models.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the avaiable data has not reached some point (such as asymptotic size).
#' @param facetGroups logical, should groups be separated in facets? Defaults to TRUE.
#' @param groupFill logical, should groups have different colors? Defaults to FALSE. If TRUE then
#' viridis colormaps are used in the order of virMaps.
#' @param virMaps order of viridis maps to use. Will be recycled to necessary length.
#' Defaults to "plasma", but will generally be informed by growthPlot's default.
#' @keywords growth-curve
#' @importFrom methods is
#' @import ggplot2
#' @importFrom stats predict update residuals
#' @importFrom nlme nlme nlme.formula
#' @examples
#'
#' simdf <- growthSim("logistic",
#'   n = 10, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#'
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group, sigma = "none",
#'   df = simdf, start = NULL, type = "nlme"
#' )
#'
#' fit <- fitGrowth(ss)
#'
#' nlmePlot(fit, form = ss$pcvrForm, groups = NULL, df = ss$df, timeRange = NULL)
#' nlmePlot(fit, form = ss$pcvrForm, groups = "a", df = ss$df, timeRange = 1:10, groupFill = TRUE)
#'
#' @return Returns a ggplot showing an nlme model's credible
#' intervals and optionally the individual growth lines.
#'
#' @export
#'

nlmePlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL, facetGroups = TRUE,
                     groupFill = FALSE, virMaps = c("plasma")) {
  #* `get needed information from formula`
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  if (individual == "dummyIndividual") {
    individual <- NULL
  }
  group <- parsed_form$group
  facetGroups <- .no_dummy_labels(group, facetGroups)
  df <- parsed_form$data
  #* `filter by groups if groups != NULL`
  if (!is.null(groups)) {
    df <- df[df[[group]] %in% groups, ]
  }
  intVar <- paste0(group, individual)
  #* `make new data if timerange is not NULL`
  if (!is.null(timeRange)) {
    new_data <- do.call(rbind, lapply(unique(df[[intVar]]), function(g) {
      stats::setNames(data.frame(g, timeRange), c(intVar, x))
    }))
    new_data[[group]] <- gsub("[.].*", "", new_data[[intVar]])
    new_data[[individual]] <- gsub(".*[.]", "", new_data[[intVar]])
    df <- df[df[[x]] >= min(timeRange) & df[[x]] <= max(timeRange), ]
  } else {
    new_data <- df
  }

  preds <- new_data
  preds$trendline <- round(predict(fit, preds), 4)
  preds <- preds[!duplicated(preds$trendline), ]
  preds <- .add_sigma_bounds(preds, fit, x, group)

  #* `plot`

  #* `facetGroups`
  facet_layer <- NULL
  if (facetGroups) {
    facet_layer <- ggplot2::facet_wrap(stats::as.formula(paste0("~", group)))
  }
  #* `groupFill`
  pal <- viridis::plasma(2, begin = 0.1, end = 0.9)
  virList <- lapply(seq_along(unique(df[[group]])), function(i) {
    pal
  })
  if (groupFill) {
    virList <- lapply(rep(virMaps, length.out = length(unique(df[[group]]))), function(pal) {
      viridis::viridis(2, begin = 0.1, end = 0.9, option = pal)
    })
  }
  #* `layer for individual lines if formula was complete`
  individual_lines <- list()
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
  }

  plot <- ggplot2::ggplot(preds, ggplot2::aes(x = .data[[x]], y = .data[["trendline"]])) +
    facet_layer +
    individual_lines +
    ggplot2::labs(x = x, y = y) +
    pcv_theme()

  for (g in seq_along(unique(preds[[group]]))) {
    iteration_group <- unique(preds[[group]])[g]
    sub <- preds[preds[[group]] == iteration_group, ]
    plot <- plot +
      ggplot2::geom_ribbon(
        data = sub, ggplot2::aes(
          ymin = .data[["sigma_ymin"]],
          ymax = .data[["sigma_ymax"]]
        ),
        fill = virList[[g]][1], alpha = 0.5
      ) +
      ggplot2::geom_line(data = sub, color = virList[[g]][2], linewidth = 0.75)
  }

  return(plot)
}

#' convenience function for calculating sigma upper and lower bounds
#' @keywords internal
#' @noRd
.add_sigma_bounds <- function(preds, fit, x, group) {
  res <- do.call(rbind, lapply(unique(preds[[group]]), function(grp) {
    varCoef <- as.numeric(fit$modelStruct$varStruct[grp])

    sub <- preds[preds[[group]] == grp, ]
    exes <- sub[[x]]

    if (methods::is(fit$modelStruct$varStruct, "varPower")) {
      out <- exes^(2 * varCoef)
    } else if (methods::is(fit$modelStruct$varStruct, "varExp")) {
      out <- exp(2 * varCoef * exes)
    } else if (methods::is(fit$modelStruct$varStruct, "varIdent")) {
      baseSigma <- fit$sigma
      varSummary <- summary(fit$modelStruct$varStruct)
      coefs <- data.frame(
        x = 1 / unique(attr(varSummary, "weight")),
        g = unique(attr(varSummary, "groups"))
      )
      out <- baseSigma * coefs[coefs$g == grp, "x"]
    }

    sub$sigma_ymax <- sub$trendline + 0.5 * out
    sub$sigma_ymin <- sub$trendline - 0.5 * out
    return(sub)
  }))

  return(res)
}

#' alias of nlmePlot for using lme models via class matching
#' @keywords internal
#' @noRd

lmePlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL, facetGroups = TRUE,
                    groupFill = FALSE, virMaps = c("plasma")) {
  nlmePlot(fit, form, df, groups, timeRange, facetGroups, groupFill, virMaps)
}
