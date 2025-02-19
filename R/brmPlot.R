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
#' @param hierarchy_value If a hierarchical model is being plotted, what value should the
#' hiearchical predictor be? If left NULL (the default) the mean value is used. If this is >1L
#' then the x axis will use the hierarchical variable from the model at the mean of the timeRange
#' (mean of x values in the model if timeRange is not specified).
#' @param vir_option Viridis color scale to use for plotting credible intervals. Defaults to "plasma".
#' @keywords growth-curve brms
#' @import ggplot2
#' @import viridis
#' @importFrom stats as.formula
#' @examples
#' \donttest{
#' simdf <- growthSim(
#'   "logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group, sigma = "spline",
#'   list("A" = 130, "B" = 10, "C" = 3),
#'   df = simdf, type = "brms"
#' )
#' fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
#' growthPlot(fit = fit, form = y ~ time | group, groups = "a", df = ss$df)
#' }
#'
#' @return Returns a ggplot showing a brms model's credible
#' intervals and optionally the individual growth lines.
#'
#' @export

brmPlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL, facetGroups = TRUE,
                    hierarchy_value = NULL, vir_option = "plasma") {
  fitData <- fit$data
  parsed_form <- .parsePcvrForm(form, df)
  if (!is.numeric(fitData[, parsed_form$x]) && !parsed_form$USEG && !parsed_form$USEID) {
    p <- .brmStaticPlot(fit, form, df, groups, vir_option, fitData, parsed_form)
    return(p)
  }
  p <- .brmLongitudinalPlot(
    fit, form, df, groups, timeRange, facetGroups,
    hierarchy_value, vir_option, fitData, parsed_form
  )
  return(p)
}

#' @keywords internal
#' @noRd
.brmStaticPlot <- function(fit, form, df = NULL, groups = NULL,
                           vir_option = "plasma", fitData, parsed_form) {
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  if (individual == "dummyIndividual") {
    individual <- NULL
  }
  df <- parsed_form$data
  probs <- seq(from = 99, to = 1, by = -2) / 100
  newData <- data.frame(
    x = unique(fitData[[x]])
  )
  colnames(newData) <- x
  predictions <- cbind(newData, predict(fit, newData, probs = probs))
  if (!is.null(groups)) {
    predictions <- predictions[predictions[[x]] %in% groups, ]
    if (!is.null(df)) {
      df <- df[df[[x]] %in% groups, ]
    }
  }
  #* `lengthen predictions`
  max_prime <- 0.99
  min_prime <- 0.01
  max_obs <- 49
  min_obs <- 1
  c1 <- (max_prime - min_prime) / (max_obs - min_obs)
  longPreds <- do.call(rbind, lapply(seq_len(nrow(predictions)), function(r) {
    sub <- predictions[r, ]
    lp <- do.call(rbind, lapply(seq(1, 49, 2), function(i) {
      min <- paste0("Q", i)
      max <- paste0("Q", 100 - i)
      iter <- sub[, c(x, "Estimate")]
      iter$q <- round(1 - (c1 * (i - max_obs) + max_prime), 2)
      iter$min <- sub[[min]]
      iter$max <- sub[[max]]
      return(iter)
    }))
    return(lp)
  }))
  #* `Make Numeric Groups`
  longPreds$numericGroup <- as.numeric(as.factor(longPreds[[x]]))
  #* `Make plot`
  p <- ggplot2::ggplot(longPreds, ggplot2::aes(x = .data[[x]], y = .data$Estimate)) +
    ggplot2::labs(x = x, y = y) +
    pcv_theme()
  p <- p +
    lapply(unique(longPreds$q), function(q) {
      ribbon_plot <- ggplot2::geom_rect(
        data = longPreds[longPreds$q == q, ],
        ggplot2::aes(
          xmin = .data[["numericGroup"]] - c(0.45 * (1 - .data[["q"]])),
          xmax = .data[["numericGroup"]] + c(0.45 * (1 - .data[["q"]])),
          ymin = .data[["min"]],
          ymax = .data[["max"]],
          group = .data[[x]],
          fill = .data[["q"]]
        ), alpha = 0.5
      )
      return(ribbon_plot)
    }) +
    viridis::scale_fill_viridis(direction = -1, option = vir_option) +
    ggplot2::labs(fill = "Credible\nInterval")
  return(p)
}

#' @keywords internal
#' @noRd

.brmLongitudinalPlotSetup <- function(fitData, timeRange, x,
                                      hierarchy_value, group, hierarchical_predictor) {
  x_plot_var <- x
  x_plot_label <- x
  allow_data_lines <- TRUE
  if (is.null(timeRange)) {
    timeRange <- unique(fitData[[x]])
  }
  if (length(hierarchy_value) > 1) {
    timeRange <- mean(timeRange, na.rm = TRUE)
    x_plot_var <- hierarchical_predictor
    x_plot_label <- paste0(hierarchical_predictor, " (", x, " = ", round(timeRange, 1), ")")
    allow_data_lines <- FALSE
  }
  if (!all(group %in% colnames(fitData))) {
    fitData[, group] <- ""
  }
  if (!is.null(hierarchical_predictor) && is.null(hierarchy_value)) {
    hierarchy_value <- mean(fitData[[hierarchical_predictor]])
  }
  if (length(hierarchy_value) == 1) {
    x_plot_label <- paste0(x, " (", hierarchical_predictor, " = ", round(hierarchy_value, 1), ")")
  }
  return(
    list(
      "timeRange" = timeRange, "x_plot_var" = x_plot_var,
      "x_plot_label" = x_plot_label,
      "allow_data_lines" = allow_data_lines, "fitData" = fitData, "hierarchy_value" = hierarchy_value
    )
  )
}

#' @keywords internal
#' @noRd

.brmLongitudinalPlot <- function(fit, form, df = NULL, groups = NULL,
                                 timeRange = NULL, facetGroups = TRUE,
                                 hierarchy_value = NULL, vir_option = "plasma",
                                 fitData, parsed_form) {
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  hierarchical_predictor <- parsed_form$hierarchical_predictor
  group <- parsed_form$group
  facetGroups <- .no_dummy_labels(group, facetGroups)
  df <- parsed_form$data
  probs <- seq(from = 99, to = 1, by = -2) / 100
  blp_setup <- .brmLongitudinalPlotSetup(
    fitData, timeRange, x,
    hierarchy_value, group, hierarchical_predictor
  )
  timeRange <- blp_setup$timeRange
  x_plot_var <- blp_setup$x_plot_var
  x_plot_label <- blp_setup$x_plot_label
  allow_data_lines <- blp_setup$allow_data_lines
  fitData <- blp_setup$fitData
  hierarchy_value <- blp_setup$hierarchy_value
  newDataArgs <- append(
    c(
      lapply(group, function(grp) {
        return(unique(fitData[[grp]]))
      }),
      list(
        "new_1"
      )
    ),
    append(
      list(timeRange),
      list(hierarchy_value)
    )
  )
  newDataArgs <- newDataArgs[!unlist(lapply(newDataArgs, is.null))]
  newData <- do.call(expand.grid, newDataArgs)
  colnames(newData) <- c(group, individual, x, hierarchical_predictor)
  if (length(group) > 1 && paste(group, collapse = ".") %in% colnames(fitData)) {
    newData[[paste(group, collapse = ".")]] <- interaction(newData[, group])
  }
  predictions <- cbind(newData, predict(fit, newData, probs = probs))

  if (!is.null(groups)) {
    keep_index <- Reduce(intersect, lapply(seq_along(groups), function(i) {
      grp <- groups[i]
      return(which(predictions[[group[i]]] %in% grp))
    }))
    predictions <- predictions[keep_index, ]
    if (!is.null(df)) {
      keep_index_df <- Reduce(intersect, lapply(seq_along(groups), function(i) {
        grp <- groups[i]
        return(which(df[[group[i]]] %in% grp))
      }))
      df <- df[keep_index_df, ]
    }
  }
  #* `facetGroups`
  facetLayer <- NULL
  if (facetGroups) {
    facetLayer <- ggplot2::facet_wrap(as.formula(paste0("~", paste(group, collapse = "+"))))
  }
  #* `lengthen predictions`
  max_prime <- 0.99
  min_prime <- 0.01
  max_obs <- 49
  min_obs <- 1
  c1 <- (max_prime - min_prime) / (max_obs - min_obs)
  longPreds <- do.call(rbind, lapply(seq_len(nrow(predictions)), function(r) {
    sub <- predictions[r, ]
    lp <- do.call(rbind, lapply(seq(1, 49, 2), function(i) {
      min <- paste0("Q", i)
      max <- paste0("Q", 100 - i)
      iter <- sub[, c(x, group, individual, "Estimate", hierarchical_predictor)]
      iter$q <- round(1 - (c1 * (i - max_obs) + max_prime), 2)
      iter$min <- sub[[min]]
      iter$max <- sub[[max]]
      return(iter)
    }))
    return(lp)
  }))
  longPreds$plot_group <- as.character(interaction(longPreds[, group]))
  #* `Make plot`
  p <- ggplot2::ggplot(longPreds, ggplot2::aes(x = .data[[x_plot_var]], y = .data$Estimate)) +
    facetLayer +
    ggplot2::labs(x = x_plot_label, y = y) +
    pcv_theme()
  p <- p +
    lapply(unique(longPreds$q), function(q) {
      ribbon_plot <- ggplot2::geom_ribbon(
        data = longPreds[longPreds$q == q, ],
        ggplot2::aes(
          ymin = .data[["min"]],
          ymax = .data[["max"]],
          group = .data[["plot_group"]],
          fill = .data[["q"]]
        ), alpha = 0.5
      )
      return(ribbon_plot)
    }) +
    viridis::scale_fill_viridis(direction = -1, option = vir_option) +
    ggplot2::labs(fill = "Credible\nInterval")

  if (!is.null(df) && individual != "dummyIndividual" && allow_data_lines) {
    df$plot_group <- as.character(interaction(df[, group]))
    p <- p + ggplot2::geom_line(
      data = df,
      ggplot2::aes(
        .data[[x_plot_var]], .data[[y]],
        group = interaction(.data[[individual]], .data[["plot_group"]])
      ),
      color = "gray20", linewidth = 0.2
    )
  }
  return(p)
}
