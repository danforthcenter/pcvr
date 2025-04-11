#' Posterior Predictive Distribution plotting or sampling for a conjugate object
#'
#' Description, maybe exported actually
#'
#' @param x Results from \link{conjugate}.
#' @param n Optionally return n draws from each posterior predictive distribution. If this is specified
#' then draws from the posterior predictive distribution are returned instead of a plot.
#'
#' @details
#'
#' @keywords internal
#' @return A patchwork of ggplots showing hypothesis test results
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
#' conj_post_pred(x)
#' conj_post_pred(x, 10)
#'
#' @noRd

conj_post_pred <- function(x, n = 0) {
  args <- as.list(x$call)
  #* find variable space RNG function
  method <- eval(args$method)[1]
  pp_dist <- switch(method,
                    t = list("rfun" = "stats::rnorm",
                             "modelled_param" = "mean",
                             "given" = "sd"), # not sure how I want this
                    # to work yet. Any other parameters will need to be given to make the posterior
                    # predictive, but they aren't as obviously retrieved in the current framework.
                    # maybe something like "look for the known_ prefix on priors/posteriors" and here
                    # I could link the name of the known parameter to the name of the argument?
                    gaussian = list("rfun" = "stats::rnorm",
                                    "modelled_param" = "mean",
                                    "given" = "sd"),
                    beta = qbeta, # this will be a trickier one? modelled_param is a list?
                    binomial = list("rfun" = "stats::rbinom",
                                    "modelled_param" = "prob",
                                    "given" = "size"), # here size isn't "known" but I do need it..?
                    # i guess it's from the trials number.
                    lognormal = list("rfun" = "stats::rlnorm",
                                     "modelled_param" = "meanlog",
                                     "given" = "sdlog"),
                    lognormal2 = list("rfun" = "stats::rlnorm",
                                      "modelled_param" = "sdlog",
                                      "given" = "meanlog"),
                    poisson = list("rfun" = "stats::rpois",
                                   "modelled_param" = "lambda"),
                    negbin = list("rfun" = "stats::rnbinom",
                                  "modelled_param" = "prob",
                                  "given" = "size"),
                    vonmises = list("rfun" = "brms::rvon_mises",
                                    "modelled_param" = "mu",
                                    "given" = "kappa"),
                    vonmises2 = circular::qvonmises, # this will all be a trickier one? modelled_param is a list?
                    uniform = list("rfun" = "stats::runif",
                                   "modelled_param" = "max",
                                   "given" = "min"),
                    pareto = list("rfun" = "extraDistr::rpareto",
                                  "modelled_param" = "a",
                                  "given" = "b"),
                    gamma = list("rfun" = "stats::rgamma",
                                 "modelled_param" = "rate",
                                 "given" = "shape"),
                    bernoulli = list("rfun" = "stats::rbinom",
                                     "modelled_param" = "prob",
                                     "given" = "size"), # need to give size of 1 here
                    exponential = list("rfun" = "stats::rexp",
                                       "modelled_param" = "rate"),
                    # these will all be a trickier one? modelled_param is a list?
                    bivariate_uniform = qunif,
                    bivariate_gaussian = qnorm,
                    bivariate_lognormal = qlnorm
  )
  #* actually simpler way to do this now that I have a print method anyway is to grab the
  #* plot_list (called plot_parameters when output) and take from there.'
  #* NEVERMIND, this is parameter space. If I want the variable space then I need to do the
  #* "switch method" above. Still this should be helpful?
  #* Raises issue about the variable space, I have a distribution for a parameter value,
  #* I could draw randomly from that distribution (or by quantiles), get the density of
  #* those parameter's resultant curve in the variable space, then get a posterior predictive
  #* by summarizing those?
  #* It should be possible to do this analytically?
  parameter_space <- lapply(x$plot_parameters, function(l) {
    dfun <- l$ddist_fun
    dfun_split <- strsplit(dfun, "::")[[1]]
    dfun_split[2] <- gsub("^d", "q", dfun_split[2])
    qfun <- paste0(dfun_split, collapse = "::")
    dfun_split[2] <- gsub("^q", "r", dfun_split[2])
    rfun <- paste0(dfun_split, collapse = "::")
    return(list("qfun"= qfun, "dfun" = dfun, "rfun" = rfun,
                "parameters" = l$parameters, "range" = l$range))
  })
  #* draw from posterior predictive and calculate quantiles
  draw <- if (as.logical(n)) {
    n
  } else {
    10000
  }
  draw_list <- lapply(seq_along(parameter_space), function(i) {
    pp_draws <- unlist(
      lapply(seq_len(draw), function(o) {
        param <- do.call(
          eval(parse(text = parameter_space[[i]]$rfun)),
          append(
            list(n = 1),
            parameter_space[[i]]$parameters
          )
        )
        pp_draw <- do.call(
          eval(parse(text = pp_dist$rfun)),
          args = Reduce(
          append, list(
            list("n" = 1),
            stats::setNames(list(param), pp_dist$modelled_param),
            stats::setNames(x$plot_parameters[[i]]$given, pp_dist$given)
          )
        )
        )
        return(pp_draw)
      })
    )
    if (as.logical(n)) {
      return(pp_draws)
    }
    probs <- seq(0.01, 0.99, 0.02)
    qs <- stats::quantile(pp_draws, probs)
    df <- data.frame(p = probs, q = qs, sample = paste0("Sample ", i), numericSample = i)
    df$ci <- c(seq(1, 49, 2), seq(49, 1, -2))
    df$limit <- rep(c("min", "max"), each = 25)
    df <- reshape2::dcast(df, sample + numericSample + ci ~ limit, value.var = "q")
    return(df)
  })
  if (as.logical(n)) {
    return(draw_list)
  }
  df <- do.call(rbind, draw_list)
  #* plot quantiles
  lowers <- seq(1, 49, 2)
  p1 <- ggplot2::ggplot() +
    lapply(lowers, function(i) {
      ribbon_layer <- ggplot2::geom_rect(
        data = df[df$ci == i, ],
        ggplot2::aes(
          xmin = .data[["numericSample"]] - (0.85 * (i / 100)),
          xmax = .data[["numericSample"]] + (0.85 * (i / 100)),
          ymin = .data[["min"]],
          ymax = .data[["max"]],
          group = .data[["sample"]],
          fill = .data[["ci"]]
        ), alpha = 0.5
      )
      return(ribbon_layer)
    }) +
    ggplot2::scale_x_continuous(labels = unique(df$sample),
                                breaks = unique(df$numericSample)) +
    viridis::scale_fill_viridis(option = "plasma",
                                breaks = c(1, 25, 49),
                                labels = c("1-99%", "25-75%", "50%")) +
    ggplot2::labs(fill = "Quantile",
                  y = "Posterior Predictive Intervals",
                  x = "Sample") +
    pcv_theme() +
    ggplot2::theme(axis.text.x.bottom = ggplot2::element_text(hjust = 0.5))
  return(p1)
}
