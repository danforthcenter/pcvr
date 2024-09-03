#' Function to visualize \code{survival::survreg} models fit by \code{fitGrowth}.
#'
#' Models fit using \link{growthSS} inputs by \link{fitGrowth}
#' (and similar models made through other means) can be visualized easily using this function.
#' This will generally be called by \code{growthPlot}.
#'
#' @param fit A model fit returned by \code{fitGrowth} with type="nls".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm}
#' part of the output) specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' If the individual and group are specified then the observed growth lines are plotted.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param df A dataframe to use in plotting observed growth curves on top of the model.
#' This must be supplied for nls models.
#' @param timeRange Ignored, included for compatibility with other plotting functions.
#' @param facetGroups logical, should groups be separated in facets? Defaults to TRUE.
#' @param groupFill logical, should groups have different colors? Defaults to FALSE.
#' If TRUE then viridis colormaps are used in the order of virMaps
#' @param virMaps order of viridis maps to use. Will be recycled to necessary length. Defaults to
#' "plasma", but will generally be informed by growthPlot's default.
#' @keywords survival
#' @importFrom methods is
#' @import ggplot2
#' @importFrom stats predict
#' @examples
#'
#'
#' df <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "survival weibull", form = y > 100 ~ time | id / group,
#'   df = df, type = "survreg"
#' )
#' fit <- fitGrowth(ss)
#' survregPlot(fit, form = ss$pcvrForm, df = ss$df)
#' survregPlot(fit, form = ss$pcvrForm, df = ss$df, groups = "a")
#' survregPlot(fit,
#'   form = ss$pcvrForm, df = ss$df, facetGroups = FALSE,
#'   groupFill = TRUE, virMaps = c("plasma", "mako")
#' )
#'
#' @return Returns a ggplot showing an survival model's survival function.
#'
#' @export

survregPlot <- function(fit, form, groups = NULL, df = NULL, timeRange = NULL, facetGroups = TRUE,
                        groupFill = FALSE, virMaps = c("plasma")) {
  #* `parse formula`
  parsed_form <- .parsePcvrForm(form, df)
  x <- parsed_form$x
  group <- parsed_form$group
  df <- parsed_form$data
  #* `filter by groups if groups != NULL`
  if (!is.null(groups)) {
    df <- df[df[[group]] %in% groups, ]
  } else {
    groups <- unique(df[[group]])
  }
  #* `generate predictions`
  pct <- seq(0.01, 0.99, 0.01)
  preds <- predict(fit,
    newdata = data.frame("group" = groups),
    type = "quantile", p = pct, se.fit = TRUE
  )
  preds <- lapply(preds, function(d) {
    matrix(d, nrow = length(groups), ncol = length(pct))
  })
  pred_df <- stats::setNames(as.data.frame(t(preds$fit)), c(paste0("est_", groups)))
  pred_df <- cbind(pred_df, stats::setNames(as.data.frame(t(preds$se.fit)), c(paste0("se_", groups))))
  pred_df$pct <- 1 - pct
  dt <- data.table::as.data.table(pred_df)
  dt_long <- data.table::melt(dt, id.vars = "pct")
  dt_long$group <- gsub(".*_", "", dt_long$variable)
  dt_long$variable <- gsub("_.*", "", dt_long$variable)
  preds <- data.frame(data.table::dcast(dt_long, pct + group ~ variable, value.var = "value"))
  #* `facetGroups`
  if (facetGroups) {
    facet_layer <- ggplot2::facet_wrap(stats::as.formula(paste0("~", group)))
  } else {
    facet_layer <- NULL
  }
  #* `groupFill`
  if (groupFill) {
    virVals <- lapply(rep(virMaps, length.out = length(unique(df[[group]]))), function(pal) {
      viridis::viridis(3, begin = 0.1, option = pal)
    })
    names(virVals) <- groups
    color_scale <- ggplot2::scale_color_manual(values = unlist(lapply(virVals, function(pal) pal[3])))
  } else {
    virVals <- lapply(rep("plasma", length.out = length(unique(df[[group]]))), function(pal) {
      viridis::viridis(3, begin = 0.1, option = pal)
    })
    names(virVals) <- groups
    color_scale <- ggplot2::scale_color_manual(values = unlist(lapply(virVals, function(pal) pal[3])))
  }
  #* `Make ggplot`
  p <- ggplot2::ggplot(preds, ggplot2::aes(
    x = .data[["est"]],
    y = .data[["pct"]], group = .data[[group]]
  )) +
    facet_layer +
    lapply(groups, function(grp) {
      ggplot2::geom_ribbon(data = preds[preds[[group]] == grp, ], ggplot2::aes(
        xmin = .data[["est"]] - (2 * .data[["se"]]),
        xmax = .data[["est"]] + (2 * .data[["se"]])
      ), fill = virVals[[grp]][1], alpha = 0.5)
    }) +
    lapply(groups, function(grp) {
      ggplot2::geom_ribbon(data = preds[preds[[group]] == grp, ], ggplot2::aes(
        xmin = .data[["est"]] - (1 * .data[["se"]]),
        xmax = .data[["est"]] + (1 * .data[["se"]])
      ), fill = virVals[[grp]][2], alpha = 0.5)
    }) +
    ggplot2::geom_line(ggplot2::aes(color = .data[[group]]), show.legend = FALSE) +
    color_scale +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = x, y = "Survival") +
    pcv_theme()
  #* `Add data as KM line`
  if (!is.null(df)) {
    km_df <- do.call(rbind, lapply(groups, function(grp) {
      sub <- df[df[[group]] == grp, ]
      do.call(rbind, lapply(seq(0, max(df[[x]]), 1), function(ti) {
        sum_events <- sum(c(sub[as.numeric(sub[[x]]) <= ti, "event"], 0))
        n_at_risk <- nrow(sub) - sum_events
        surv_pct <- n_at_risk / nrow(sub)
        iter <- data.frame(
          group = grp, time = ti, events = sum_events,
          at_risk = n_at_risk, surv_pct = surv_pct
        )
        colnames(iter)[1] <- group
        iter
      }))
    }))
    p <- p + ggplot2::geom_line(data = km_df, ggplot2::aes(
      x = .data[[x]],
      y = .data[["surv_pct"]],
      group = .data[[group]],
      linetype = .data[[group]]
    ), color = "black", show.legend = FALSE)
  }


  return(p)
}
