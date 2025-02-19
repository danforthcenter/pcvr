#' Function to visualize \code{quantreg::rq} general additive growth models.
#'
#' Models fit using \link{growthSS} inputs by \link{fitGrowth}
#' (and similar models made through other means) can be visualized easily using this function.
#' This will generally be called by \code{growthPlot}.
#'
#' @param fit A model fit, or list of model fits, returned by \code{fitGrowth}
#' with type="nlrq" and model="gam".
#' @param form A formula similar to that in \code{growthSS} inputs (or the \code{pcvrForm}
#' part of the output) specifying the outcome,
#' predictor, and grouping structure of the data as \code{outcome ~ predictor|individual/group}.
#' If the individual and group are specified then the observed growth lines are plotted.
#' @param df A dataframe to use in plotting observed growth curves on top of the model.
#' This must be supplied for rq models.
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
#'
#' simdf <- growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "gam", form = y ~ time | id / group,
#'   tau = c(0.25, 0.5, 0.75), df = simdf, start = NULL, type = "nlrq"
#' )
#' fits <- fitGrowth(ss)
#' rqPlot(fits, form = ss$pcvrForm, df = ss$df, groupFill = TRUE)
#' rqPlot(fits, form = ss$pcvrForm, df = ss$df, groups = "a", timeRange = 1:10)
#'
#' ss <- growthSS(
#'   model = "gam", form = y ~ time | group,
#'   tau = c(0.5), df = simdf, start = NULL, type = "nlrq"
#' )
#' fit <- fitGrowth(ss)
#' rqPlot(fit, form = ss$pcvrForm, df = ss$df, groupFill = TRUE)
#'
#' @return Returns a ggplot showing an rq general additive model's quantiles
#'  and optionally the individual growth lines.
#'
#' @export

rqPlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL, facetGroups = TRUE,
                   groupFill = FALSE, virMaps = c("plasma")) {
  #* `get needed information from formula`
  parsed_form <- .parsePcvrForm(form, df)
  #* `pick longitudinal or non-longitudinal helper`
  if (!is.numeric(df[, parsed_form$x]) && !parsed_form$USEG && !parsed_form$USEID) {
    p <- .rqStaticPlot(
      fit, form, df, groups, timeRange,
      facetGroups, groupFill, virMaps, parsed_form
    )
    return(p)
  }
  p <- .rqLongitudinalPlot(
    fit, form, df, groups, timeRange,
    facetGroups, groupFill, virMaps, parsed_form
  )
  return(p)
}

#' @keywords internal
#' @noRd

.rqStaticPlot <- function(fit, form, df, groups, timeRange,
                          facetGroups, groupFill, virMaps, parsed_form) {
  x <- parsed_form$x
  df <- parsed_form$data

  if (methods::is(fit, "rq")) {
    fit <- list(fit)
  }

  summary_df <- do.call(rbind, lapply(fit, function(model) {
    iter_df <- as.data.frame(coef(summary(model)))
    colnames(iter_df) <- c("est", "err", "t", "p")
    iter_df[[x]] <- rownames(iter_df)
    iter_df[1, x] <- paste0(x, unique(df[[x]])[1])
    iter_df[["est"]] <- cumsum(iter_df[["est"]])
    iter_df$tau <- model$tau
    return(iter_df)
  }))

  #* `filter by groups if groups != NULL`
  if (!is.null(groups)) {
    summary_df <- summary_df[summary_df[[x]] %in% groups, ]
  }
  #* `facetGroups`
  facet_layer <- NULL
  if (facetGroups) {
    facet_layer <- ggplot2::facet_wrap(stats::as.formula(paste0("~", x)),
      scales = "free_x"
    )
  }
  #* `groupFill`
  n_taus <- length(unique(summary_df$tau))

  virpal_p1 <- viridis::plasma(ceiling(n_taus / 2), direction = 1, end = 1)
  virpal_p2 <- viridis::plasma(ceiling(n_taus / 2), direction = -1, end = 1)[-1]
  virpal <- c(virpal_p1, virpal_p2)
  virList <- lapply(seq_along(unique(summary_df[[x]])), function(i) {
    return(virpal)
  })

  if (groupFill) {
    virList <- lapply(rep(virMaps, length.out = length(unique(summary_df[[x]]))), function(pal) {
      virpal_p1 <- viridis::viridis(ceiling(n_taus / 2),
        direction = 1, end = 1, option = pal
      )
      virpal_p2 <- viridis::viridis(ceiling(n_taus / 2),
        direction = -1, end = 1, option = pal
      )[-1]
      return(c(virpal_p1, virpal_p2))
    })
  }

  #* `plot`
  plot <- ggplot(summary_df, ggplot2::aes(group = interaction(.data[[x]]))) +
    facet_layer +
    labs(x = x, y = as.character(form)[2]) +
    pcv_theme()

  for (g in seq_along(unique(summary_df[[x]]))) {
    iteration_group <- unique(summary_df[[x]])[g]
    sub <- summary_df[summary_df[[x]] == iteration_group, ]
    for (i in seq_along(unique(sub$tau))) {
      inner_sub <- sub[sub$tau == unique(sub$tau)[i], ]
      plot <- plot +
        ggplot2::geom_errorbar(data = inner_sub, ggplot2::aes(
          x = .data[[x]],
          ymin = .data[["est"]] - 2 * .data[["err"]],
          ymax = .data[["est"]] + 2 * .data[["err"]]
        ), width = 0.15, color = virList[[g]][i]) +
        ggplot2::geom_point(
          data = inner_sub, ggplot2::aes(
            x = .data[[x]],
            y = .data[["est"]]
          ),
          color = virList[[g]][i], size = 4
        ) +
        ggplot2::geom_text(
          data = inner_sub, ggplot2::aes(
            x = .data[[x]],
            y = .data[["est"]],
            label = .data[["tau"]]
          ),
          size = 2, color = "white"
        )
    }
  }
  plot
  return(plot)
}

#' @keywords internal
#' @noRd

.rqLongitudinalPlot <- function(fit, form, df, groups, timeRange,
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
  df[[paste(group, collapse = ".")]] <- interaction(df[, group])
  group <- paste(group, collapse = ".")
  #* `filter by groups if groups != NULL`
  if (!is.null(groups)) {
    df <- df[df[[group]] %in% groups, ]
  }
  #* `make new data if timerange is not NULL`
  if (!is.null(timeRange)) {
    new_data <- do.call(rbind, lapply(unique(df[[group]]), function(g) {
      return(stats::setNames(data.frame(g, timeRange), c(group, x)))
    }))
    df <- df[df[[x]] >= min(timeRange) & df[[x]] <= max(timeRange), ]
  } else {
    new_data <- df
  }
  #* `standardize fit class`
  if (methods::is(fit, "rq")) {
    fit <- list(fit)
    names(fit) <- fit$tau
  }
  #* `add predictions and record taus`
  preds <- do.call(cbind, lapply(fit, function(f) {
    tau <- f$tau
    return(stats::setNames(data.frame(stats::predict(f, newdata = new_data)), paste0("Q_", tau)))
  }))
  predCols <- colnames(preds)
  keep <- which(!duplicated(preds))
  plotdf <- cbind(df[keep, ], preds[keep, ])
  colnames(plotdf) <- c(colnames(df), colnames(preds))
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
        direction = -1, end = 1, option = pal
      )[-1]
      return(c(virpal_p1, virpal_p2))
    })
  } else {
    virpal_p1 <- viridis::plasma(ceiling(length(predCols) / 2), direction = 1, end = 1)
    virpal_p2 <- viridis::plasma(ceiling(length(predCols) / 2), direction = -1, end = 1)[-1]
    virpal <- c(virpal_p1, virpal_p2)
    virList <- lapply(seq_along(unique(df[[group]])), function(i) {
      return(virpal)
    })
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
        line_layer <- ggplot2::geom_line(
          data = sub, ggplot2::aes(x = .data[[x]], y = .data[[predCols[i]]]),
          color = virList[[g]][i], linewidth = 0.7
        )
        return(line_layer)
      })
  }
  return(plot)
}
