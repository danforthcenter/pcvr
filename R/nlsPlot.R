#' Function to visualize common \code{stats::nls} growth models.
#'
#' Models fit using \link{growthSS} inputs by \link{fitGrowth}
#' (and similar models made through other means) can be visualized easily using this function.
#' This will generally be called by \code{growthPlot}.
#'
#' @param fit A model fit returned by \code{fitGrowth} with type="nls".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm}
#' part of the output) specifying the outcome, predictor, and grouping structure of the data as
#' \code{outcome ~ predictor|individual/group}.
#' If the individual and group are specified then the observed growth lines are plotted.
#' @param groups An optional set of groups to keep in the plot.
#' Defaults to NULL in which case all groups in the model are plotted.
#' @param df A dataframe to use in plotting observed growth curves on top of the model.
#' This must be supplied for nls models.
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
#' @importFrom stats predict
#' @examples
#'
#'
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group,
#'   df = simdf, start = NULL, type = "nls"
#' )
#' fit <- fitGrowth(ss)
#' nlsPlot(fit, form = ss$pcvrForm, df = ss$df, groupFill = TRUE)
#' nlsPlot(fit, form = ss$pcvrForm, df = ss$df, groups = "a", timeRange = 1:10)
#'
#' @return Returns a ggplot showing an nls model's predictions.
#'
#' @export

nlsPlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL,
                    facetGroups = TRUE, groupFill = FALSE, virMaps = c("plasma")) {
  #* `get needed information from formula`
  parsed_form <- .parsePcvrForm(form, df)
  #* `pick longitudinal or non-longitudinal helper`
  if (!is.numeric(df[, parsed_form$x]) && !parsed_form$USEG && !parsed_form$USEID) {
    p <- .nlsStaticPlot(
      fit, form, df, groups, timeRange,
      facetGroups, groupFill, virMaps, parsed_form
    )
    return(p)
  }
  p <- .nlsLongitudinalPlot(
    fit, form, df, groups, timeRange,
    facetGroups, groupFill, virMaps, parsed_form
  )
  return(p)
}

#' @keywords internal
#' @noRd

.nlsStaticPlot <- function(fit, form, df, groups, timeRange,
                           facetGroups, groupFill, virMaps, parsed_form) {
  x <- parsed_form$x
  df <- parsed_form$data
  #* `when implemented SE can be added here, see ?predict.nls`
  summary_df <- as.data.frame(coef(summary(fit)))
  colnames(summary_df) <- c("est", "err", "t", "p")
  summary_df[[x]] <- rownames(summary_df)
  summary_df[1, x] <- paste0(x, unique(df[[x]])[1])
  summary_df[["est"]] <- cumsum(summary_df[["est"]])
  #* `filter by groups if groups != NULL`
  if (!is.null(groups)) {
    summary_df <- summary_df[summary_df[[x]] %in% groups, ]
  }
  #* `facetGroups`
  facet_layer <- NULL
  if (facetGroups) {
    facet_layer <- ggplot2::facet_wrap(stats::as.formula(paste0("~", x)))
  }
  #* `groupFill`
  virVals <- unlist(lapply(
    rep(virMaps, length.out = length(unique(summary_df[[x]]))),
    function(pal) {
      return(viridis::viridis(1, begin = 0.5, option = pal))
    }
  ))
  color_scale <- ggplot2::scale_color_manual(values = virVals)
  if (!groupFill) {
    color_scale <- ggplot2::scale_color_manual(values = rep("#CC4678FF", length(unique(df[[x]]))))
  }
  #* `plot`
  plot <- ggplot(summary_df, ggplot2::aes(group = interaction(.data[[x]]))) +
    facet_layer +
    ggplot2::geom_errorbar(ggplot2::aes(
      x = .data[[x]],
      ymin = .data[["est"]] - 2 * .data[["err"]],
      ymax = .data[["est"]] + 2 * .data[["err"]]
    ), width = 0.25) +
    ggplot2::geom_point(ggplot2::aes(x = .data[[x]], y = .data[["est"]], color = .data[[x]]),
      size = 4
    ) +
    color_scale +
    labs(x = x, y = as.character(form)[2]) +
    pcv_theme()
  return(plot)
}

#' @keywords internal
#' @noRd

.nlsLongitudinalPlot <- function(fit, form, df, groups, timeRange,
                                 facetGroups, groupFill, virMaps, parsed_form) {
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  if (individual == "dummyIndividual") {
    individual <- NULL
  }
  group <- parsed_form$group
  facetGroups <- .no_dummy_labels(group, facetGroups)
  df <- parsed_form$data
  df$group_interaction <- interaction(df[, group])
  #* `filter by groups if groups != NULL`
  if (!is.null(groups)) {
    keep_index_df <- Reduce(intersect, lapply(seq_along(groups), function(i) {
      grp <- groups[i]
      return(which(df[[group[i]]] %in% grp))
    }))
    df <- df[keep_index_df, ]
  }
  #* `make new data if timerange is not NULL`
  if (!is.null(timeRange)) {
    new_data <- do.call(
      expand.grid,
      append(
        list(timeRange),
        c(lapply(group, function(grp) {
          return(unique(df[[grp]]))
        }))
      )
    )
    colnames(new_data) <- c(x, group)
    df <- df[df[[x]] >= min(timeRange) & df[[x]] <= max(timeRange), ]
  } else {
    new_data <- NULL
  }
  #* `add predictions`
  preds <- data.frame(pred = stats::predict(fit, newdata = new_data))
  keep <- which(!duplicated(preds$pred))
  plotdf <- df[keep, ]
  plotdf$pred <- preds[keep, "pred"]
  #* `when implemented SE can be added here, see ?predict.nls`
  #*
  #* `layer for individual lines if formula was complete`
  individual_lines <- list()
  if (!is.null(individual)) {
    individual_lines <- ggplot2::geom_line(
      data = df, ggplot2::aes(
        x = .data[[x]], y = .data[[y]],
        group = interaction(
          .data[[individual]],
          .data[["group_interaction"]]
        )
      ),
      linewidth = 0.25, color = "gray40"
    )
  }
  #* `facetGroups`
  facet_layer <- NULL
  if (facetGroups) {
    facet_layer <- ggplot2::facet_wrap(stats::as.formula("~group_interaction"))
  }
  #* `groupFill`
  if (groupFill) {
    virVals <- unlist(lapply(
      rep(virMaps, length.out = length(unique(df[["group_interaction"]]))),
      function(pal) {
        return(viridis::viridis(1, begin = 0.5, option = pal))
      }
    ))
    color_scale <- ggplot2::scale_color_manual(values = virVals)
  } else {
    color_scale <- ggplot2::scale_color_manual(values = rep(
      "#CC4678FF",
      length(unique(df[["group_interaction"]]))
    ))
  }

  #* `plot`
  plot <- ggplot(plotdf, ggplot2::aes(group = interaction(.data[["group_interaction"]]))) +
    facet_layer +
    individual_lines +
    ggplot2::geom_line(
      ggplot2::aes(
        x = .data[[x]], y = .data[["pred"]],
        color = .data[["group_interaction"]]
      ),
      linewidth = 0.7,
      show.legend = groupFill
    ) + # using middle of plasma pal
    color_scale +
    labs(x = x, y = as.character(form)[2], color = group) +
    pcv_theme()

  return(plot)
}

#' @rdname nlsPlot
#' @examples
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "gam", form = y ~ time | id / group,
#'   df = simdf, start = NULL, type = "nls"
#' )
#' fit <- fitGrowth(ss)
#' gamPlot(fit, form = ss$pcvrForm, df = ss$df, groupFill = TRUE)
#' gamPlot(fit, form = ss$pcvrForm, df = ss$df, groups = "a", timeRange = 1:10)
#' ss <- growthSS(
#'   model = "gam", form = y ~ time | group,
#'   df = simdf, start = NULL, type = "nls"
#' )
#' fit <- fitGrowth(ss)
#' gamPlot(fit, form = ss$pcvrForm, df = ss$df, groupFill = TRUE)
#'
#' @export

gamPlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL, facetGroups = TRUE,
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
  df <- parsed_form$data
  df$group_interaction <- interaction(df[, group])
  #* `filter by groups if groups != NULL`
  if (!is.null(groups)) {
    keep_index_df <- Reduce(intersect, lapply(seq_along(groups), function(i) {
      grp <- groups[i]
      return(which(df[[group[i]]] %in% grp))
    }))
    df <- df[keep_index_df, ]
  }
  #* `make new data if timerange is not NULL`
  if (!is.null(timeRange)) {
    new_data <- do.call(
      expand.grid,
      append(
        list(timeRange),
        c(lapply(group, function(grp) {
          return(unique(df[[grp]]))
        }))
      )
    )
    colnames(new_data) <- c(x, group)
    df <- df[df[[x]] >= min(timeRange) & df[[x]] <= max(timeRange), ]
  } else {
    # note this is the only change between this and nlsPlot
    # this change is here because predict.nls sometimes acts strangely with the given data
    # but predict.gam does not accept a NULL input for the newdata argument.
    new_data <- df
  }
  #* `add predictions`

  preds <- data.frame(pred = stats::predict(fit, newdata = new_data))
  keep <- which(!duplicated(preds$pred))
  plotdf <- df[keep, ]
  plotdf$pred <- preds[keep, "pred"]

  #* `when implemented SE can be added here, see ?predict.nls`
  #*
  #* `layer for individual lines if formula was complete`
  individual_lines <- list()
  if (!is.null(individual)) {
    individual_lines <- ggplot2::geom_line(
      data = df, ggplot2::aes(
        x = .data[[x]], y = .data[[y]],
        group = interaction(
          .data[[individual]],
          .data[["group_interaction"]]
        )
      ),
      linewidth = 0.25, color = "gray40"
    )
  }
  #* `facetGroups`
  facet_layer <- NULL
  if (facetGroups) {
    facet_layer <- ggplot2::facet_wrap(stats::as.formula("~group_interaction"))
  }
  #* `groupFill`
  if (groupFill) {
    virVals <- unlist(lapply(
      rep(virMaps, length.out = length(unique(df[["group_interaction"]]))),
      function(pal) {
        return(viridis::viridis(1, begin = 0.5, option = pal))
      }
    ))
    color_scale <- ggplot2::scale_color_manual(values = virVals)
  } else {
    color_scale <- ggplot2::scale_color_manual(values = rep(
      "#CC4678FF",
      length(unique(df[["group_interaction"]]))
    ))
  }

  #* `plot`
  plot <- ggplot(plotdf, ggplot2::aes(group = interaction(.data[["group_interaction"]]))) +
    facet_layer +
    individual_lines +
    ggplot2::geom_line(
      ggplot2::aes(
        x = .data[[x]], y = .data[["pred"]],
        color = .data[["group_interaction"]]
      ),
      linewidth = 0.7,
      show.legend = groupFill
    ) + # using middle of plasma pal
    color_scale +
    labs(x = x, y = as.character(form)[2], color = group) +
    pcv_theme()

  return(plot)
}

#' @rdname nlsPlot
#' @examples
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "gam", form = y ~ time | id / group,
#'   df = simdf, start = NULL, type = "nls"
#' )
#' fit <- fitGrowth(ss)
#' lmPlot(fit, form = ss$pcvrForm, df = ss$df)
#' @export

lmPlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL, facetGroups = TRUE,
                   groupFill = FALSE, virMaps = c("plasma")) {
  p <- nlsPlot(fit, form, df, groups, timeRange, facetGroups, groupFill, virMaps)
  return(p)
}
