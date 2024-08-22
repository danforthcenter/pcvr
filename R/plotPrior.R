#' Check priors used in ease of use brms functions
#'
#' @param priors A named list of means for prior distributions.
#' This takes the same input as the prior argument of \code{\link{growthSS}}.
#' Alternatively, if given the output of growthSS this will preform a prior predictive check
#' and return a plot from \code{\link{growthPlot}} of that check ignoring all other arguments.
#' Note that all priors must be
#' proper in that case (non-flat) and the fit is likely to be strange looking due to how thick
#' tailed the default priors from \code{\link{growthSS}} are.
#' @param type Either "density", the default, or a model as would be specified in \code{growthSS}
#' or \code{growthSim} such as "logistic", "gompertz", "monomolecular", "exponential",
#' "linear", "power law", "double logistic", or "double gompertz".
#' If this is a model type then n draws from the prior will be simulated as growth
#' trendlines and densities will be plotted on margins for some distributions.
#' @param n Numeric, if type is a model then how many draws from the prior should be simulated?
#' @param t Numeric, time passed to growthSim. Defaults to 25 (the growthSim default).
#' @keywords Bayesian brms priors
#' @return A named list of plots showing prior distributions that \code{growthSS} would use,
#' optionally with a plot of simulated growth curves using draws from those priors.
#' @import ggplot2
#' @import patchwork
#' @importFrom stats rlnorm dlnorm
#' @examples
#'
#' ## Not run:
#'
#' set.seed(123)
#' priors <- list("A" = c(100, 130), "B" = c(10, 8), "C" = c(0.2, 0.1))
#' plotPrior(priors)
#'
#' plotPrior(priors, "gompertz")[[1]]
#'
#' ## End(Not run:)
#'
#' @export

plotPrior <- function(priors, type = "density", n = 200, t = 25) {
  if ("prior" %in% names(priors)) {
    p <- .brms_prior_predictive(priors)
    return(p)
  }
  densPlots <- lapply(seq_along(priors), function(i) {
    pri <- priors[[i]]
    nm <- names(priors)[i]

    pri_df <- do.call(rbind, lapply(seq_along(pri), function(o) {
      prio <- pri[o]
      max <- ceiling(max(rlnorm(1000, log(max(pri)), 0.25)) * 1.1)
      support <- seq(0, max, length.out = 10000)
      dens <- dlnorm(support, meanlog = log(prio), sdlog = 0.25)
      pdf <- dens / sum(dens)
      data.frame(
        support = support,
        dens = pdf,
        param = nm,
        item = as.character(o)
      )
    }))

    ggplot2::ggplot(pri_df, ggplot2::aes(
      x = .data$support, y = .data$dens,
      fill = .data$item, group = .data$item
    )) +
      ggplot2::geom_polygon(alpha = 0.5) +
      ggplot2::theme_minimal() +
      ggplot2::labs(y = "Density", title = nm, fill = "Prior")
  })
  names(densPlots) <- names(priors)

  if (type == "density") {
    out <- densPlots
  } else {
    simdf <- do.call(rbind, lapply(1:n, function(i) {
      iter_params <- .prior_sampler(priors)
      x <- growthSim(model = type, n = 1, t = t, params = iter_params)
      x$id <- paste0("id_", i)
      x
    }))

    model_plot <- ggplot2::ggplot(
      simdf,
      ggplot2::aes(
        x = .data$time, y = .data$y, group = interaction(.data$id, .data$group),
        color = .data$group
      )
    ) +
      ggplot2::geom_line(linewidth = 0.1) +
      ggplot2::theme_minimal() +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 5))) +
      ggplot2::labs(
        y = "Y", title = paste0(n, " curves simulated from prior draws"),
        color = "Prior"
      ) +
      ggplot2::theme(legend.position = "inside",
                     legend.position.inside = c(0.9, 0.9))

    if (type %in% c("logistic", "gompertz", "weibull", "frechet", "gumbel")) {
      x <- "B"
      y <- "A"
      z <- "C"
    } else if (type %in% c("monomolecular")) {
      y <- "A"
      x <- z <- NULL
    } else {
      x <- y <- z <- NULL
    }
    out <- .plotPriorMarginPlots(model_plot, densPlots, x, y, z)
  }

  return(out)
}

#' @description
#' Internal function for adding marginal plots to plotPrior
#' @keywords internal
#' @noRd

.plotPriorMarginPlots <- function(model_plot, densPlots, x, y, z) {
  xLims <- ggplot2::layer_scales(model_plot)$x$range$range
  yLims <- ggplot2::layer_scales(model_plot)$y$range$range
  model_plot_solo <- model_plot

  sum_non_null <- 0
  x_margin_plot <- NULL
  y_margin_plot <- NULL
  z_margin_plot <- NULL

  if (!is.null(y)) {
    y_margin_plot <- densPlots[[y]] +
      ggplot2::scale_y_reverse(position = "right") +
      ggplot2::scale_x_continuous(position = "top", limits = yLims) +
      ggplot2::labs(x = "Asymptote Prior") +
      ggplot2::coord_flip() +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
        legend.position = "none"
      )
    sum_non_null <- sum_non_null + 1
  }

  if (!is.null(x)) {
    x_margin_plot <- densPlots[[x]] +
      ggplot2::labs(x = "Inflection Point Prior") +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::coord_cartesian(xlim = xLims)
    sum_non_null <- sum_non_null + 1
  }

  if (!is.null(z)) {
    z_margin_plot <- densPlots[[z]] +
      ggplot2::labs(x = "Growth Rate Prior") +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::scale_x_continuous(n.breaks = 3)
    sum_non_null <- sum_non_null + 1
  }

  if (sum_non_null == 3) {
    design <- c(
      patchwork::area(1, 1, 6, 6), # model plot
      patchwork::area(7, 1, 7, 6), # x margin
      patchwork::area(1, 7, 6, 7), # y margin
      patchwork::area(7, 7, 7, 7)
    ) # "z" margin
    layout <- patchwork::plot_layout(design = design)
    model_plot <- model_plot + x_margin_plot + y_margin_plot + z_margin_plot + layout
  } else if (!is.null(y_margin_plot) && is.null(x_margin_plot)) {
    design <- c(
      patchwork::area(1, 1, 6, 6), # model plot
      patchwork::area(1, 7, 6, 7)
    ) # y margin
    model_plot <- model_plot + y_margin_plot + patchwork::plot_layout(design = design)
  }

  if (is(model_plot, "patchwork")) {
    densPlots[[length(densPlots) + 1]] <- model_plot_solo
  }
  return(list("simulated" = model_plot, "distributions" = densPlots))
}




#' @description
#' Internal function for drawing from priors
#' @param priors priors as a list
#' @keywords internal
#' @noRd

.prior_sampler <- function(priors) {
  lapply(priors, function(pri) { # draw sample from prior
    unlist(lapply(pri, function(mu) {
      rlnorm(1, log(mu), 0.25)
    }))
  })
}

#' @description
#' Internal function for sampling a growthSS model's priors only
#' @param priors a list returned by growthSS
#' @keywords internal
#' @noRd

.brms_prior_predictive <- function(priors = NULL) {
  dp <- brms::get_prior(priors$formula, priors$df, priors$family)
  ssp <- priors$prior
  dpi <- as.character(interaction(dp$coef, dp$dpar, dp$nlpar))
  sspi <- as.character(interaction(ssp$coef, ssp$dpar, ssp$nlpar))
  priors$prior <- rbind(ssp, dp[-which(dpi %in% sspi), ])
  tryCatch(
    {
      m <- suppressMessages(
        fitGrowth(priors,
          iter = 1000,
          chains = 1,
          cores = 1,
          sample_prior = "only",
          silent = 2
        )
      )
    },
    error = function(err) {
      message(paste0(
        "Error trying to sample from priors distributions.",
        "All priors must be proper (non-flat).\nAttempting to sample from: "
      ))
      print(priors$prior)
      message("The original Error message is:")
      stop(conditionMessage(err))
    }
  )
  p <- growthPlot(m, form = priors$pcvrForm, df = priors$df)
  return(p)
}
