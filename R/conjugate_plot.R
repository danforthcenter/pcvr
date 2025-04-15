#' Plot a \code{conjugate} object.
#'
#' @aliases plot.conjugate
#'
#' @param object An object of class \code{conjugate}.
#' @param ... further arguments, ignored.
#' @import ggplot2
#' @import patchwork
#' @examples
#' x <- conjugate(
#'   s1 = rnorm(10, 50, 10), s2 = rnorm(10, 60, 12), method = "t",
#'   priors = list(list(mu = 40, sd = 10), list(mu = 45, sd = 8)),
#'   rope_range = c(-5, 8), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal"
#' )
#' plot(x)
#'
#' @method plot conjugate
#' @export

plot.conjugate <- function(x) {
  dirSymbol <- NULL
  if ("hyp" %in% colnames(x$summary)) {
    dirSymbol <- switch(x$summary$hyp[1],
      "unequal" = "!=",
      "equal" = "=",
      "greater" = ">",
      "lesser" = "<"
    )
  }
  args <- as.list(x$call)
  plot <- .conj_plot(
    res = x,
    rope_df = x$rope_df,
    rope_range = eval(args$rope_range),
    rope_ci = eval(args$rope_ci),
    dirSymbol = dirSymbol,
    support = x$plot_parameters[[1]]$range,
    method = eval(args$method),
    bayes_factor = eval(args$bayes_factor)
  )
  return(plot)
}

#' ***********************************************************************************************
#' *************** `General Plotting Function` ***********************************
#' ***********************************************************************************************
#' Used to pick which kind of plotting function to use.
#' @keywords internal
#' @noRd

.conj_plot <- function(res, rope_df = NULL,
                       rope_range, rope_ci, dirSymbol = NULL, support, method, bayes_factor) {
  if (any(grepl("bivariate", method))) {
    p <- .conj_bivariate_plot(res, rope_df, rope_range, rope_ci, dirSymbol)
  } else {
    p <- .conj_general_plot(res, rope_df, rope_range, rope_ci, dirSymbol, support, bayes_factor)
  }
  return(p)
}

#' ***********************************************************************************************
#' *************** `General Plotting Function` ***********************************
#' ***********************************************************************************************


#' @keywords internal
#' @noRd
.conj_general_plot <- function(res, rope_df = NULL,
                               rope_range, rope_ci, dirSymbol = NULL, support, bayes_factor) {
  s1_plot_df <- data.frame(range = res$plot_parameters[[1]]$range)
  s1_bf_title <- NULL
  s2_bf_title <- NULL
  if (!is.null(bayes_factor)) {
    s1_bf_title <- paste0(
      ", BF: (",
      paste(bayes_factor, collapse = ", "),
      ") ", round(res$summary$bf_1, 2)
    )
  }

  p <- ggplot2::ggplot(s1_plot_df, ggplot2::aes(x = .data$range)) +
    ggplot2::stat_function(
      geom = "polygon",
      fun = eval(parse(text = res$plot_parameters[[1]]$ddist_fun)),
      args = res$plot_parameters[[1]]$parameters,
      ggplot2::aes(fill = "s1"), alpha = 0.5
    ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = res$summary$HDI_1_low),
      color = "red",
      linewidth = 1.1
    ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = res$summary$HDE_1),
      color = "red", linetype = "dashed", linewidth = 1.1
    ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = res$summary$HDI_1_high),
      color = "red",
      linewidth = 1.1
    ) +
    ggplot2::scale_fill_manual(values = "red") +
    ggplot2::labs(
      x = "Posterior Distribution of Random Variable", y = "Density", title = "Distribution of Samples",
      subtitle = paste0(
        "HDE: ", round(res$summary$HDE_1, 2),
        "\nHDI: [", round(res$summary$HDI_1_low, 2), ", ",
        round(res$summary$HDI_1_high, 2), "]", s1_bf_title
      )
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 0.5))) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.9, 0.9)
    ) +
    pcv_theme()

  if (length(res$data) == 2) {
    if (!is.null(bayes_factor)) {
      s2_bf_title <- paste0(
        ", BF: (",
        paste(bayes_factor, collapse = ", "),
        ") ", round(res$summary$bf_2, 2)
      )
    }

    if (res$summary$post.prob < 1e-5) {
      post.prob.text <- "<1e-5"
    } else {
      post.prob.text <- round(res$summary$post.prob, 5)
    }

    fill_scale <- which(sapply(p$scales$scales, function(x) {
      return("fill" %in% x$aesthetics) # avoid "replacing scale" non-messages
    }))

    p$scales$scales[[fill_scale]] <- NULL
    p <- p +
      ggplot2::stat_function(
        geom = "polygon",
        fun = eval(parse(text = res$plot_parameters[[2]]$ddist_fun)),
        args = res$plot_parameters[[2]]$parameters,
        ggplot2::aes(fill = "s2"), alpha = 0.5
      ) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = res$summary$HDI_2_low),
        color = "blue",
        linewidth = 1.1
      ) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = res$summary$HDE_2),
        color = "blue",
        linetype = "dashed", linewidth = 1.1
      ) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = res$summary$HDI_2_high),
        color = "blue",
        linewidth = 1.1
      ) +
      ggplot2::scale_fill_manual(values = c("red", "blue"), breaks = c("s1", "s2")) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 0.5))) +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)
      ) +
      ggplot2::labs(subtitle = paste0(
        "Sample 1:  ", round(res$summary$HDE_1, 2), " [", round(res$summary$HDI_1_low, 2), ", ",
        round(res$summary$HDI_1_high, 2), "]", s1_bf_title, "\n",
        "Sample 2:  ", round(res$summary$HDE_2, 2), " [", round(res$summary$HDI_2_low, 2), ", ",
        round(res$summary$HDI_2_high, 2), "]", s2_bf_title, "\n",
        "P[p1", dirSymbol, "p2] = ", post.prob.text
      ))
  }

  if (!is.null(rope_df)) {
    p <- p + ggplot2::ggplot(rope_df, ggplot2::aes(x = .data$X)) +
      ggplot2::geom_histogram(bins = 100, fill = "purple", color = "purple", alpha = 0.7) +
      ggplot2::geom_histogram(
        data = data.frame(
          "X" = rope_df[
            rope_df$X > rope_range[1] &
              rope_df$X < rope_range[2] &
              rope_df$X > res$summary$HDI_rope_low &
              rope_df$X < res$summary$HDI_rope_high,
          ]
        ),
        bins = 100, fill = "gray30", color = "gray30"
      ) +
      ggplot2::annotate("segment",
        x = rope_range[1], xend = rope_range[2], y = 0, yend = 0,
        linewidth = 2, color = "gray70"
      ) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = res$summary$HDI_rope_low), linewidth = 0.7) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = res$summary$HDE_rope),
        linetype = "dashed",
        linewidth = 0.7
      ) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = res$summary$HDI_rope_high), linewidth = 0.7) +
      ggplot2::labs(
        x = "Posterior of Difference", y = "Frequency", title = "Distribution of Difference",
        subtitle = paste0(
          "Median Difference of ", round(res$summary$HDE_rope, 2), "\n",
          100 * rope_ci, "% CI [", round(res$summary$HDI_rope_low, 2), ", ",
          round(res$summary$HDI_rope_high, 2), "]\n",
          100 * rope_ci, "% HDI in [", rope_range[1], ", ", rope_range[2], "]: ",
          round(res$summary$rope_prob, 2)
        )
      ) +
      pcv_theme() +
      ggplot2::theme(
        axis.title.y.left = ggplot2::element_blank(),
        axis.text.y.left = ggplot2::element_blank()
      ) +
      patchwork::plot_layout(widths = c(2, 1))
  }
  return(p)
}


#' ***********************************************************************************************
#' *************** `Bivariate Plotting Function` ***********************************
#' ***********************************************************************************************
#'
#' @keywords internal
#' @noRd

.conj_bivariate_plot <- function(res, rope_df = NULL, rope_range, rope_ci, dirSymbol = NULL, support) {
  TITLE <- "Joint Posterior Distribution"
  if (length(res$data) == 1) {
    margin_plot_df <- res$plot_df[[1]]
    params <- unique(margin_plot_df$param)
    #* `Make joint distribution plot`
    joint_dist_s1 <- res$posteriorDraws[[1]]
    limits <- lapply(params, function(p) {
      sub <- margin_plot_df[margin_plot_df$param == p, ]
      return(range(sub$range))
    })
    names(limits) <- params
    x_lim <- limits[[1]]
    y_lim <- limits[[2]]
    v_lines <- res$summary[res$summary$param == params[1], ]
    h_lines <- res$summary[res$summary$param == params[2], ]

    joint_p <- ggplot2::ggplot(
      joint_dist_s1,
      ggplot2::aes(
        x = .data[[params[1]]], y = .data[[params[2]]],
        group = .data[["sample"]]
      )
    ) +
      ggplot2::geom_hline(
        yintercept = h_lines$HDI_1_low, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = h_lines$HDI_1_high, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_vline(
        xintercept = v_lines$HDI_1_low, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_vline(
        xintercept = v_lines$HDI_1_high, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      geom_density_2d_filled(breaks = ~ pretty(., n = 51)[-1], alpha = 0.9) +
      ggplot2::scale_fill_viridis_d(option = "plasma") +
      ggplot2::xlim(x_lim) +
      ggplot2::ylim(y_lim) +
      pcv_theme() +
      ggplot2::theme(legend.position = "none")
    #* `Make marginal distribution plot of each parameter (x, y)`
    margin_plots <- lapply(params, function(par) {
      par_plot <- ggplot2::ggplot(
        margin_plot_df[margin_plot_df$param == par, ],
        ggplot2::aes(x = .data$range, y = .data$prob)
      ) +
        ggplot2::geom_area(
          data = margin_plot_df[margin_plot_df$param == par, ],
          alpha = 0.5, ggplot2::aes(fill = "s1")
        ) +
        ggplot2::geom_vline(
          data = data.frame(),
          ggplot2::aes(
            xintercept = res$summary[res$summary$param == par, "HDI_1_low"]
          ),
          color = "red",
          linewidth = 0.5
        ) +
        ggplot2::geom_vline(
          data = data.frame(),
          ggplot2::aes(
            xintercept = res$summary[res$summary$param == par, "HDE_1"]
          ),
          color = "red", linetype = "dashed", linewidth = 0.5
        ) +
        ggplot2::geom_vline(
          data = data.frame(),
          ggplot2::aes(
            xintercept = res$summary[res$summary$param == par, "HDI_1_high"]
          ),
          color = "red", linetype = "dashed", linewidth = 0.5
        ) +
        ggplot2::scale_fill_manual(values = "red") +
        ggplot2::xlim(limits[[par]]) +
        ggplot2::theme_void() +
        ggplot2::theme(legend.title = ggplot2::element_blank())
      return(par_plot)
    })
    #* `Write title if there is only 1 sample`
    SUBTITLE <- NULL
  } else if (length(res$data) == 2) {
    #* `Make plots for sample 2 if it exists`
    margin_plot_df <- do.call(rbind, lapply(1:2, function(i) {
      md <- res$plot_df[[i]]
      md$sample <- paste0("Sample ", i)
      return(md)
    }))
    params <- unique(margin_plot_df$param)
    joint_dist <- do.call(rbind, lapply(1:2, function(i) {
      pd <- res$posterior_draws[[i]]
      pd$sample <- paste0("Sample ", i)
      return(pd)
    }))
    #* `Define Limits`
    limits <- lapply(params, function(p) {
      sub <- margin_plot_df[margin_plot_df$param == p, ]
      return(range(sub$range))
    })
    names(limits) <- params
    x_lim <- limits[[1]]
    y_lim <- limits[[2]]
    #* `Make joint distribution plot`
    v_lines <- res$summary[res$summary$param == params[1], ]
    h_lines <- res$summary[res$summary$param == params[2], ]

    joint_p <- ggplot2::ggplot(
      joint_dist,
      ggplot2::aes(
        x = .data[[params[1]]], y = .data[[params[2]]],
        group = .data[["sample"]]
      )
    ) +
      ggplot2::geom_hline(
        yintercept = h_lines$HDI_1_low, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = h_lines$HDI_1_high, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = h_lines$HDI_2_low, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_hline(
        yintercept = h_lines$HDI_2_high, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_vline(
        xintercept = v_lines$HDI_1_low, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_vline(
        xintercept = v_lines$HDI_1_high, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_vline(
        xintercept = v_lines$HDI_2_low, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      ggplot2::geom_vline(
        xintercept = v_lines$HDI_2_high, color = "gray80", linetype = 5,
        linewidth = 0.25
      ) +
      geom_density_2d_filled(breaks = ~ pretty(., n = 51)[-1], alpha = 0.9) +
      ggplot2::scale_fill_viridis_d(option = "plasma") +
      ggplot2::xlim(x_lim) +
      ggplot2::ylim(y_lim) +
      pcv_theme() +
      ggplot2::theme(legend.position = "none")

    #* `Make marginal distribution plot of each parameter (x, y)`
    margin_plots <- lapply(params, function(par) {
      hdf <- res$summary
      hdf <- hdf[hdf$param == par, ]
      par_plot <- ggplot2::ggplot() +
        ggplot2::geom_area(
          data = margin_plot_df[margin_plot_df$param == par, ],
          alpha = 0.5, ggplot2::aes(
            x = .data$range, y = .data$prob,
            fill = .data$sample
          ),
          position = "identity"
        ) +
        ggplot2::geom_vline(
          data = data.frame(),
          ggplot2::aes(
            xintercept = hdf[, "HDI_1_low"]
          ),
          color = "red", linetype = "dashed",
          linewidth = 0.5
        ) +
        ggplot2::geom_vline(
          data = data.frame(),
          ggplot2::aes(
            xintercept = hdf[, "HDE_1"]
          ),
          color = "red", linewidth = 0.5
        ) +
        ggplot2::geom_vline(
          data = data.frame(),
          ggplot2::aes(
            xintercept = hdf[, "HDI_1_high"]
          ),
          color = "red", linetype = "dashed", linewidth = 0.5
        ) +
        ggplot2::geom_vline(
          data = data.frame(),
          ggplot2::aes(
            xintercept = hdf[, "HDI_2_low"]
          ),
          color = "blue", linetype = "dashed",
          linewidth = 0.5
        ) +
        ggplot2::geom_vline(
          data = data.frame(),
          ggplot2::aes(
            xintercept = hdf[, "HDE_2"]
          ),
          color = "blue", linewidth = 0.5
        ) +
        ggplot2::geom_vline(
          data = data.frame(),
          ggplot2::aes(
            xintercept = hdf[, "HDI_2_high"]
          ),
          color = "blue", linetype = "dashed", linewidth = 0.5
        ) +
        ggplot2::scale_fill_manual(values = c("red", "blue")) +
        ggplot2::xlim(limits[[par]]) +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "inside",
          legend.position.inside = c(0.1, 0.5),
          legend.title = ggplot2::element_blank()
        )
      return(par_plot)
    })

    post.probs <- lapply(params, function(par) {
      hdf <- res$summary
      hdf <- hdf[hdf$param == par, ]
      if (hdf$post.prob < 1e-5) {
        post.prob.text <- "<1e-5"
      } else {
        post.prob.text <- round(hdf$post.prob, 5)
      }
      return(post.prob.text)
    })
    names(post.probs) <- params

    #* `Write title if there are 2 samples`
    SUBTITLE <- paste(lapply(params, function(par) {
      par_string <- paste0(par, ": P[s1", dirSymbol[[1]], "s2] = ", post.probs[[par]])
      return(par_string)
    }), collapse = "\n")
  }
  #* `Assemble Patchwork`
  layout <- c(
    patchwork::area(2, 1, 3, 2),
    patchwork::area(1, 1, 1, 2),
    patchwork::area(2, 3, 3, 3)
  )
  margin_plots[[2]] <- margin_plots[[2]] +
    ggplot2::coord_flip() +
    ggplot2::theme(legend.position = "none")

  p <- joint_p + margin_plots[[1]] + margin_plots[[2]] +
    patchwork::plot_layout(design = layout) &
    patchwork::plot_annotation(title = TITLE, subtitle = SUBTITLE)
  return(p)
}
