#' Prior sensitivity function for conjugate outputs
#'
#' Prior sensitivity analysis can be an important part of understanding Bayesian results.
#' Priors should not be adapted after the fact to make your results "less sensitive", this is
#' only meant to help you understand how much of your result comes from your data vs the choice
#' of a prior. This function performs prior sensitivity analysis on results from the \link{conjugate}
#' function.
#'
#' @param x Results from \link{conjugate}.
#' @param priors A list of priors similar to how they are specified in conjugate but named for the
#' distribution you plan to use, see details and examples.
#' @param n The number of priors to use, defaults to 100.
#'
#' @details
#' Runs prior sensitivity analysis for a conjugate result.
#'
#' Priors here are specified using a named list. For instance, for normal priors with means between
#' 5 and 20 and standard deviations between 5 and 10 the prior argument would be
#' \code{list("rnorm" = list("mean" = c(5, 20), "sd" = c(5, 10))))}.
#' The priors that are used in sensitivity analysis are drawn randomly from within the ranges specified
#' by the provided list. If you are unsure what random-generation function to use then check the
#' \link{conjugate} docs where the distributions are listed for each method in the details section.
#'
#' @keywords internal
#' @return A patchwork of ggplots showing hypothesis test results, interpretation changes, and
#' the new priors used.
#' @examples
#' set.seed(123)
#' x <- conjugate(
#'   s1 = rnorm(10, 10, 1), s2 = rnorm(10, 13, 1.5), method = "t",
#'   priors = list(list(mu = 175, sd = 2), # purposefully poor strong priors
#'                 list(mu = 50, sd = 2)),
#'   plot = FALSE, rope_range = c(-8, 8), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "unequal",
#'   bayes_factor = c(50, 55)
#' )
#' x$summary$post.prob # result of terrible priors
#' .prior_sens_conj(x, priors = list("rnorm" = list("mean" = c(5, 20), "sd" = c(5, 10))), n = 100)
#' @importFrom stats setNames
#' @import patchwork
#' @import ggplot2
#' @noRd


.prior_sens_conj <- function(x, priors, n = 100) {
  args <- as.list(x$call)
  pri_df <- do.call(cbind, lapply(seq_along(priors[[1]]), function(i) {
    x <- seq(priors[[1]][[i]][1], priors[[1]][[i]][2], length.out = 1000)
    x[sample(1000, size = n)]
  }))
  pri_df <- as.data.frame(pri_df)
  colnames(pri_df) <- names(priors[[1]])
  probs <- unlist(lapply(seq_len(n), function(i) {
    post <- conjugate(
      s1 = x$data[[1]], # will need to be recovered from conjugate somehow
      s2 = x$data[[2]], # will need to be recovered from conjugate somehow
      priors = stats::setNames(pri_df[i, ], names(x$prior[[1]])),
      hypothesis = args$hypothesis
    )$summary$post.prob
    return(post)
  }))
  df <- data.frame(x = probs)
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["x"]])) +
    ggplot2::geom_histogram(bins = 100, color = "white") +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::percent_format())+
    ggplot2::labs(title = "Hypothesis Tests",
                  x = "Posterior Probability",
                  y = "Count") +
    pcv_theme()
  conj_data <- x$data
  conj_data <- conj_data[!is.null(conj_data)]
  methods <- rep(list(args$method), length(conj_data))
  support_priors <- rep(eval(args$priors), length(conj_data))
  
  boundaries <- do.call(cbind, lapply(seq_len(n), function(i) {
    qnts <- do.call(gsub("^r", "q", names(priors)[1]),
                    args = append(
                      list("p" = c(0.001, 0.999)),
                      as.list(pri_df[i, ])
                      )
                    )
    return(qnts)
  }))
  lower <- min(boundaries[1, ])
  upper <- max(boundaries[2, ])
  
  df2 <- data.frame(x = seq(lower, upper, length.out = 10))
  p2 <- ggplot2::ggplot(df2, aes(x = x)) +
    lapply(seq_len(n), function(i) {
      layer <- ggplot2::stat_function(geom = "polygon",
                                      fun = eval(parse(text = gsub("^r", "d", names(priors)[1]))),
                                      args = pri_df[i, ], alpha = 0.1, fill = "gray60",
                                      linewidth = 0.25, color = "black")
      return(layer)
    }) +
    ggplot2::labs(title = "Random Priors",
                  x = "Prior Distribution") +
    pcv_theme() +
    ggplot2::theme(axis.text.y.left = ggplot2::element_blank(),
                   axis.title.y.left = ggplot2::element_blank(),
                   axis.line.y.left = ggplot2::element_blank())
  
  df$cat <- cut(df$x,
                breaks = c(-0.01, 0.01, 0.05, 0.95, 0.99, 1),
                labels = c("Very Low\n(0, 0.01]", "Low\n(0.01, 0.05]", "Middle\n(0.05, 0.95]", 
                           "High\n(0.95, 0.99]", "Very High\n(0.99, 1]"))
  original <- cut(x$summary$post.prob,
                  breaks = c(-0.01, 0.01, 0.05, 0.95, 0.99, 1),
                  labels = c("Very Low\n(0, 0.01]", "Low\n(0.01, 0.05]", "Middle\n(0.05, 0.95]", 
                             "High\n(0.95, 0.99]", "Very High\n(0.99, 1]"))
  n_changes <- sum(df$cat != original)
  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = cat)) +
    ggplot2::geom_bar() +
    ggplot2::scale_x_discrete(limits = rev(c("Very Low\n(0, 0.01]", "Low\n(0.01, 0.05]",
                                             "Middle\n(0.05, 0.95]", 
                                             "High\n(0.95, 0.99]", "Very High\n(0.99, 1]"))) +
    ggplot2::annotate("point", x = original, y = 0, size = 3, color = "red") +
    ggplot2::labs(title = paste0("Interpretation Changes: ", n_changes)) +
    pcv_theme() +
    ggplot2::theme(axis.text.y.left = ggplot2::element_blank(),
                   axis.title.y.left = ggplot2::element_blank(),
                   axis.line.y.left = ggplot2::element_blank(),
                   axis.title.x.bottom = ggplot2::element_blank(),
                   axis.text.x.bottom = ggplot2::element_text(hjust = 0.5))
  design <- c(
    patchwork::area(1, 1, 2, 1),
    patchwork::area(2, 2, 2, 2),
    patchwork::area(1, 2, 1, 2)
  )
  patch <- p1 + p2 + p3 + patchwork::plot_layout(design = design)
  return(patch)
}
