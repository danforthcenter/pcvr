#' Function to visualize models made by \link{fitGrowth}.
#'
#' Models fit using \link{growthSS} inputs by \link{fitGrowth}
#' (and similar models made through other means) can be visualized easily using this function.
#'
#'
#' @param fit A model fit object (or a list of \code{nlrq} models) as returned by \code{fitGrowth}.
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm} part of the
#' output) specifying the outcome, predictor, and grouping structure of the data as
#' \code{outcome ~ predictor|individual/group}. Generally this is given directly from
#' the growthSS output (\code{ss$pcvrForm}). If the formula does not include both individuals
#' and groups then lines from the data will not be plotted which may be best if your data does not
#' specify unique individuals and your model does not include autocorrelation.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param df A dataframe to use in plotting observed growth curves on top of the model and for making
#' predictions.
#' @param timeRange An optional range of times to use. This can be used to view predictions for
#' future data if the avaiable data has not reached some point (such as asymptotic size).
#' @param facetGroups logical, should groups be separated in facets? Defaults to TRUE.
#' @param groupFill logical, should groups have different colors? Defaults to the opposite of
#' facetGroups. If TRUE then
#' viridis colormaps are used in the order c('plasma', 'mako', 'viridis', 'inferno', 'cividis', 'magma',
#' 'turbo', 'rocket'). Alternatively this can be given as a vector of
#' viridis colormap names to use in a different order than above.
#' Note that for brms models this is ignored except if used to specify a different viridis color map
#' to use.
#' @param hierarchy_value If a hierarchical model is being plotted, what value should the
#' hiearchical predictor be? If left NULL (the default) the mean value is used.
#' @keywords growth-curve
#' @importFrom methods is
#' @examples
#'
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group,
#'   df = simdf, type = "nls"
#' )
#' fit <- fitGrowth(ss)
#' growthPlot(fit, form = ss$pcvrForm, df = ss$df)
#'
#' @return Returns a ggplot showing a brms model's credible
#' intervals and optionally the individual growth lines.
#'
#' @export

growthPlot <- function(fit, form, groups = NULL, df = NULL, timeRange = NULL,
                       facetGroups = TRUE, groupFill = !facetGroups, hierarchy_value = NULL) {
  if (is.logical(groupFill)) {
    virMaps <- c("plasma", "mako", "viridis", "inferno", "cividis", "magma", "turbo", "rocket")
  } else {
    virMaps <- groupFill
    groupFill <- TRUE
  }
  model_class <- class(fit)[1]
  if (methods::is(fit, "list")) {
    model_class <- class(fit[[1]])
  }
  if (methods::is(fit, "brmsfit")) {
    if (attr(fit$formula$formula, "nl")) { # non linear models are growth models
      plot <- brmPlot(
        fit = fit, form = form, groups = groups, df = df, timeRange = timeRange,
        facetGroups = facetGroups, hierarchy_value = hierarchy_value, vir_option = virMaps[1]
      )
    } else { # linear models are survival models
      plot <- brmSurvPlot(
        fit = fit, form = form, groups = groups, df = df, timeRange = timeRange,
        facetGroups = facetGroups
      )
    }
  } else {
    plottingFunction <- match.fun(paste0(model_class, "Plot"))
    plot <- plottingFunction(
      fit = fit, form = form, groups = groups, df = df, timeRange = timeRange,
      facetGroups = facetGroups, groupFill = groupFill, virMaps
    )
  }
  return(plot)
}
