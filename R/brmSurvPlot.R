#' Function to visualize brms survival models specified using growthSS.
#'
#' Models fit using \link{growthSS} inputs by \link{fitGrowth}
#' (and similar models made through other means)
#'  can be visualized easily using this function. This will generally be called by \code{growthPlot}.
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
#'   set.seed(123)
#'   df <- growthSim("exponential",
#'     n = 20, t = 50,
#'     params = list("A" = c(1, 1), "B" = c(0.15, 0.1))
#'   )
#'   ss1 <- growthSS(
#'     model = "survival weibull", form = y > 100 ~ time | id / group,
#'     df = df, start = c(0, 5)
#'   )
#'   fit1 <- fitGrowth(ss1, iter = 600, cores = 2, chains = 2, backend = "cmdstanr")
#'   brmSurvPlot(fit1, form = ss1$pcvrForm, df = ss1$df)
#'
#'   # note that using the cumulative hazard to calculate survival is likely to underestimate
#'   # survival in these plots if events do not start immediately.
#'   ss2 <- growthSS(
#'     model = "survival binomial", form = y > 100 ~ time | id / group,
#'     df = df, start = c(-4, 3)
#'   )
#'   fit2 <- fitGrowth(ss2, iter = 600, cores = 2, chains = 2, backend = "cmdstanr")
#'   brmSurvPlot(fit2, form = ss2$pcvrForm, df = ss2$df)
#' }
#'
#' ## End(Not run)
#'
#' @return Returns a ggplot showing a brms model's credible
#' intervals and optionally the individual growth lines.
#'
#' @export

brmSurvPlot <- function(fit, form, df = NULL, groups = NULL, timeRange = NULL, facetGroups = TRUE,
                        groupFill = FALSE, virMaps = c("plasma")) {
  family <- as.character(fit$family)[1]

  if (family == "weibull") {
    p <- .weibullBrmSurvPlot(
      fit = fit, form = form, df = df, groups = groups,
      timeRange = timeRange, facetGroups = facetGroups,
      groupFill = groupFill, virMaps = c("plasma")
    )
  } else if (family == "binomial") {
    p <- .binomialBrmSurvPlot(
      fit = fit, form = form, df = df, groups = groups,
      timeRange = timeRange, facetGroups = facetGroups,
      groupFill = groupFill, virMaps = c("plasma")
    )
  }
  return(p)
}

#' Internal Plotting function for weibull survival times
#' brms predicts mu where: scale = exp(mu) / (gamma(1 + 1 / shape))
#' @keywords internal
#' @noRd

.binomialBrmSurvPlot <- function(fit, form, df = NULL, groups = NULL,
                                 timeRange = NULL, facetGroups = TRUE,
                                 groupFill = FALSE, virMaps = c("plasma")) {
  #* `pull model data`
  fitData <- fit$data
  #* `general pcvr formula parsing`
  parsed_form <- .parsePcvrForm(form, df)
  x <- parsed_form$x
  group <- parsed_form$group
  df <- parsed_form$data
  #* `set groups to use`
  if (is.null(groups)) {
    groups <- unique(fitData[[group]])
  }
  #* `pull draws and convert to survival`
  draws <- as.matrix(fit)
  s_hat <- brms::inv_logit_scaled(draws[, grepl("^b_", colnames(draws))])
  colnames(s_hat) <- gsub("b_", "haz_", colnames(s_hat))
  s_hat <- 1 - s_hat # take complement of hazard
  s_hat <- as.data.frame(
    do.call(rbind, lapply(groups, function(grp) {
      s_hat_grp <- s_hat[, grepl(paste0(group, grp, "$"), colnames(s_hat))]
      x <- do.call(cbind, lapply(seq_along(s_hat_grp), function(j) {
        do.call(rbind, lapply(seq_len(nrow(s_hat_grp)), function(i) {
          cumprod(s_hat_grp[i, 1:j])[j]
        }))
      }))
      colnames(x) <- gsub("haz_", "surv_", colnames(s_hat_grp))
      colnames(x) <- gsub(paste0(":", group, grp), "", colnames(x))
      x <- as.data.frame(x)
      x[[group]] <- grp
      x
    }))
  )

  #* `define probabilities and take quantiles`
  probs <- seq(from = 99, to = 1, by = -2) / 100
  quantiles <- as.data.frame(
    do.call(rbind, lapply(groups, function(grp) {
      s_hat_grp <- s_hat[s_hat[[group]] == grp, ]
      grp_quantile <- as.data.frame(do.call(rbind, lapply(1:(ncol(s_hat_grp) - 1), function(i) {
        quantile(s_hat_grp[, i], probs)
      })))
      colnames(grp_quantile) <- paste0("Q", seq(99, 1, -2))
      nms <- colnames(s_hat_grp)[grepl("surv_", colnames(s_hat_grp))]
      grp_quantile[[x]] <- as.numeric(gsub(paste0("surv_", x), "", nms))
      grp_quantile[[group]] <- grp
      grp_quantile
    }))
  )
  quantiles <- quantiles[quantiles[[group]] %in% groups, ]
  #* `Decide faceting`
  if (facetGroups) {
    if (!all(fitData[[group]] == "a")) {
      facetLayer <- ggplot2::facet_wrap(as.formula(paste0("~", group)))
    } else {
      facetLayer <- NULL
    }
  } else {
    facetLayer <- NULL
  }
  #* `Decide Fill Colors`
  if (groupFill) {
    virList <- lapply(rep(virMaps, length.out = length(groups)), function(pal) {
      viridis::viridis(n = length(probs), option = pal)
    })
  } else {
    pal <- viridis::plasma(n = length(probs))
    virList <- lapply(seq_along(unique(df[[group]])), function(i) {
      pal
    })
  }
  #* `Initialize Plot`
  p <- ggplot2::ggplot(quantiles, ggplot2::aes(x = .data[[x]])) +
    facetLayer +
    ggplot2::labs(x = x, y = "Survival") +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    pcv_theme()
  #* `Add Ribbons`
  for (g in seq_along(groups)) {
    iteration_group <- groups[g]
    sub <- quantiles[quantiles[[group]] == iteration_group, ]
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
  #* `Add KM Trend`
  if (!is.null(df)) {
    df$pct_surv <- (1 - df$pct_event)
    df[[x]] <- as.numeric(df[[x]])
    p <- p + ggplot2::geom_line(data = df, ggplot2::aes(
      x = .data[[x]],
      y = .data[["pct_surv"]],
      group = .data[[group]],
      linetype = .data[[group]]
    ), color = "black", show.legend = FALSE)
  }
  return(p)
}


#' Internal Plotting function for weibull survival times
#' brms predicts mu where: scale = exp(mu) / (gamma(1 + 1 / shape))
#' @keywords internal
#' @noRd

.weibullBrmSurvPlot <- function(fit, form, df = NULL, groups = NULL,
                                timeRange = NULL, facetGroups = TRUE,
                                groupFill = FALSE, virMaps = c("plasma")) {
  #* `Transform draws`
  fitData <- fit$data
  fdf <- as.data.frame(fit)
  mu_cols <- colnames(fdf)[grepl("b_", colnames(fdf))]
  scales <- do.call(cbind, lapply(mu_cols, function(col) {
    exp(fdf[[col]]) / (gamma(1 + (1 / fdf$shape)))
  }))
  colnames(scales) <- paste0("scale_", mu_cols)
  fdf <- cbind(fdf, scales)
  #* `general pcvr formula parsing`
  parsed_form <- .parsePcvrForm(form, df)
  x <- parsed_form$x
  group <- parsed_form$group
  df <- parsed_form$data
  #* `further survival formula steps`

  #* `Define Time Range`
  if (is.null(timeRange)) {
    timeRange <- seq(0, round(max(fitData[[x]]), -1), length.out = 100)
  }
  #* `set groups to use`
  if (is.null(groups)) {
    groups <- unique(fitData[[group]])
  }
  #* `Make Survival Quantiles`
  probs <- seq(from = 99, to = 1, by = -2) / 100
  quantiles <- do.call(rbind, lapply(groups, function(grp) {
    test <- fdf[, c(paste0("scale_b_", group, grp), "shape")]
    colnames(test) <- c("scale", "shape")
    metrics <- do.call(rbind, lapply(probs, function(i) {
      shape <- quantile(test[, "shape"], probs = i)
      scale <- quantile(test[, "scale"], probs = i)
      do.call(rbind, lapply(timeRange, function(t) {
        t_row <- data.frame(
          probs = i,
          h = shape / scale * (t / scale)^(shape - 1),
          cdf = 1 - exp(1)^-(t / scale)^shape,
          pdf = shape / scale * (t / scale)^(shape - 1) * exp(1)^-(t / scale)^shape
        )
        t_row[[x]] <- t
        t_row$St <- 1 - t_row$cdf
        t_row$c <- -log(t_row$St)
        t_row
      }))
    }))

    metrics$probs2 <- round(metrics$probs * 100, 2)
    metrics <- data.table::as.data.table(metrics)
    wide_metrics <- as.data.frame(data.table::dcast(metrics, time ~ paste0("Q", probs2),
      value.var = "St", fun.aggregate = mean
    ))
    wide_metrics[[group]] <- grp
    wide_metrics
  }))
  #* `Decide faceting`
  if (facetGroups) {
    if (!all(fitData[[group]] == "a")) {
      facetLayer <- ggplot2::facet_wrap(as.formula(paste0("~", group)))
    } else {
      facetLayer <- NULL
    }
  } else {
    facetLayer <- NULL
  }
  #* `Decide Fill Colors`
  if (groupFill) {
    virList <- lapply(rep(virMaps, length.out = length(groups)), function(pal) {
      viridis::viridis(n = length(probs), option = pal)
    })
  } else {
    pal <- viridis::plasma(n = length(probs))
    virList <- lapply(seq_along(unique(df[[group]])), function(i) {
      pal
    })
  }
  #* `Initialize Plot`
  p <- ggplot2::ggplot(quantiles, ggplot2::aes(x = .data[["time"]])) +
    facetLayer +
    ggplot2::labs(x = x, y = "Survival") +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    pcv_theme()
  #* `Add Ribbons`
  for (g in seq_along(groups)) {
    iteration_group <- groups[g]
    sub <- quantiles[quantiles[[group]] == iteration_group, ]
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
  #* `Add KM Trend`
  if (!is.null(df)) {
    km_df <- do.call(rbind, lapply(groups, function(grp) {
      sub <- df[df[[group]] == grp, ]
      do.call(rbind, lapply(timeRange, function(ti) {
        sum_events <- sum(c(sub[as.numeric(sub[[x]]) <= ti, "event"], 0))
        n_at_risk <- nrow(sub) - sum_events
        surv_pct <- n_at_risk / nrow(sub)
        data.frame(
          group = grp, time = ti, events = sum_events,
          at_risk = n_at_risk, surv_pct = surv_pct
        )
      }))
    }))
    p <- p + ggplot2::geom_line(data = km_df, ggplot2::aes(
      x = .data[["time"]],
      y = .data[["surv_pct"]],
      group = .data[["group"]],
      linetype = .data[["group"]]
    ), color = "black", show.legend = FALSE)
  }
  return(p)
}
