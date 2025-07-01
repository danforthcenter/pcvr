#' Bayesian testing using conjugate priors and method of moments for single or multi value traits.
#'
#' @description
#' Function to perform bayesian tests and ROPE comparisons using single or multi value traits with
#' several distributions.
#'
#' @param s1 A data.frame or matrix of multi value traits or a vector of single value traits.
#' If a multi value trait is used then column names should include a number representing the "bin".
#' Alternatively for distributions other than "binomial" (which requires list data
#' with "successes" and "trials" as numeric vectors in the list, see examples)
#' this can be a formula specifying \code{outcome ~ group} where group has exactly 2
#' levels. If using wide MV trait data then the formula should specify column positions ~ grouping
#' such as \code{1:180 ~ group}.
#' This sample is shown in red if plotted.
#' @param s2 An optional second sample, or if s1 is a formula then this should be a dataframe.
#' This sample is shown in blue if plotted.
#' @param method The distribution/method to use.
#' Currently "t", "gaussian", "beta", "binomial", "lognormal", "lognormal2", "poisson",
#' "negbin" (negative binomial), "uniform", "pareto", "gamma", "bernoulli", "exponential",
#' "vonmises", and "vonmises2" are supported.
#' The count (binomial, poisson and negative binomial), bernoulli, exponential,
#' and pareto distributions are only implemented for single value traits due to their updating
#' and/or the nature of the input data.
#' The "t" and "gaussian" methods both use a T distribution with "t" testing for a difference
#' of means and "gaussian" testing for a difference in the distributions (similar to a Z test).
#' Both Von Mises options are for use with circular data (for instance hue values when the circular
#' quality of the data is relevant). Note that non-circular distributions can be compared to each other.
#' This should only be done with caution and may not be supported in all downstream functions.
#' There are also 3 bivariate conjugate priors that are supported for use with single value data.
#' Those are "bivariate_uniform", "bivariate_gaussian" and "bivariate_lognormal".
#' @param priors Prior distributions described as a list of lists.
#' If this is a single list then it will be duplicated for the second sample,
#' which is generally a good idea if both
#' samples use the same distribution (method).
#' Elements in the inner lists should be named for the parameter they represent (see examples).
#' These names vary by method (see details).
#'  By default this is NULL and weak priors (generally jeffrey's priors) are used.
#'  The \code{posterior} part of output can also be recycled as a new prior if Bayesian
#'  updating is appropriate for your use.
#' @param plot deprecated, use \code{plot} method instead.
#' @param rope_range Optional vector specifying a region of practical equivalence.
#' This interval is considered practically equivalent to no effect.
#' Kruschke (2018) suggests c(-0.1, 0.1) as a broadly reasonable ROPE for standardized parameters.
#' That range could also be rescaled by a standard deviation/magnitude for
#' non-standardized parameters, but ultimately this should be informed by your
#' setting and scientific question.
#' See Kruschke (2018) for details on ROPE and other Bayesian methods to aide
#' decision-making \doi{10.1177/2515245918771304}
#' and \doi{10.1037/a0029146}.
#' @param rope_ci The credible interval probability to use for ROPE. Defaults to 0.89.
#' @param cred.int.level The credible interval probability to use
#' in computing HDI for samples, defaults to 0.89.
#' @param hypothesis Direction of a hypothesis if two samples are provided.
#'  Options are "unequal", "equal", "greater", and "lesser",
#'   read as "sample1 greater than sample2".
#' @param bayes_factor Optional point or interval to evaluate bayes factors on. Note that this
#' generally only makes sense to use if you have informative priors where the change in odds between
#' prior and posterior is meaningful about the data. If this is non-NULL then columns of bayes factors
#' are added to the summary output. Note these are only implemented for univariate distributions.
#' @param support Deprecated
#'
#' @import bayestestR
#' @import extraDistr
#' @importFrom stats var median rnorm dnorm qbeta rbeta dbeta qbeta rlnorm
#' dlnorm qlnorm qt rgamma qgamma dgamma rpois runif
#'
#' @details
#'
#' Prior distributions default to be weakly informative and in some cases you may wish to change them.
#' \itemize{
#'    \item{\strong{"t", "gaussian", and "lognormal":} \code{priors = list(mu = 0, sd = 10)},
#'     where mu is the mean and sd is standard deviation. For the lognormal method these describe
#'     the normal distribution of the mean parameter for lognormal data and are on the log scale.}
#'    \item{\strong{"beta", "bernoulli", and "binomial":}
#'    \code{priors = list(a = 0.5, b = 0.5)},
#'     where a and b are shape parameters of the beta distribution. Note that for the binomial
#'     distribution this is used as the prior for success probability P,
#'     which is assumed to be beta distributed as in a beta-binomial distribution.}
#'    \item{\strong{"lognormal2": } \code{priors = list(a = 1, b = 1) },
#'    where a and b are the shape and scale parameters of the gamma distribution of lognormal data's
#'    precision parameter (using the alternative mu, precision parameterization).
#'    }
#'    \item{\strong{"gamma": } \code{priors = list(shape = 0.5, scale = 0.5, known_shape = 1)},
#'     where shape and scale are the respective parameters of the gamma distributed rate
#'     (inverse of scale) parameter of gamma distributed data.}
#'    \item{\strong{"poisson" and "exponential": } \code{priors = list(a = 0.5,b = 0.5)},
#'     where a and b are shape parameters of the gamma distribution.}
#'    \item{\strong{"negbin": } \code{priors = list(r = 10, a = 0.5, b = 0.5)},
#'     where r is the r parameter of the negative binomial distribution
#'     (representing the number of successes required)
#'      and where a and b are shape parameters of the beta distribution.
#'      Note that the r value is not updated.
#'       The conjugate beta prior is only valid when r is fixed and known,
#'       which is a limitation for this method.}
#'     \item{\strong{"uniform": } \code{list(scale = 0.5, location = 0.5)}, where scale is the
#'     scale parameter of the pareto distributed upper boundary and location is the location parameter
#'     of the pareto distributed upper boundary. Note that different sources will use different
#'     terminology for these parameters. These names were chosen for consistency with the
#'     \code{extraDistr} implementation of the pareto distribution. On Wikipedia the parameters are
#'     called shape and scale, corresponding to extraDistr's scale and location respectively, which
#'     can be confusing. Note that the lower boundary of the uniform is assumed to be 0.
#'     }
#'     \item{\strong{"pareto": } \code{list(a = 1, b = 1, known_location = min(data))}, where
#'     a and b are the shape and scale parameters of the gamma distribution of the pareto distribution's
#'     scale parameter. In this case location is assumed to be constant and known, which is less of
#'     a limitation than knowing r for the negative binomial method since location will generally be
#'     right around/just under the minimum of the sample data. Note that the pareto method is only
#'     implemented currently for single value traits since one of the statistics needed to update
#'     the gamma distribution here is the product of the data and we do not currently have a method
#'     to calculate a similar sufficient statistic from multi value traits.
#'     }
#'     \item{\strong{"vonmises": } \code{list(mu = 0, kappa = 0.5, boundary = c(-pi, pi),
#'     known_kappa = 1, n = 1)}, where mu is the direction of the circular distribution (the mean),
#'     kappa is the precision of the mean, boundary is a vector including the two values that are the
#'     where the circular data "wraps" around the circle, known_kappa is the fixed value of precision
#'     for the total distribution, and n is the number of prior observations. This Von Mises option
#'     updates the conjugate prior for the mean direction, which is itself Von-Mises distributed. This
#'     in some ways is analogous to the T method, but assuming a fixed variance when the mean is
#'     updated. Note that due to how the rescaling works larger circular boundaries can be slow to
#'     plot.
#'     }
#'     \item{\strong{"vonmises2": } \code{priors = list(mu = 0, kappa = 0.5,
#'     boundary = c(-pi, pi), n = 1)}, where mu and kappa are mean direction and precision of the
#'     von mises distribution, boundary is a vector including the two values that are the
#'     where the circular data "wraps" around the circle, and n is the number of prior observations.
#'     This Von-Mises implementation does not assume constant variance and instead uses MLE to estimate
#'     kappa from the data and updates the kappa prior as a weighted average of the data and the prior.
#'     The mu parameter is then updated per Von-Mises conjugacy.
#'     }
#'     \item{\strong{"bivariate_uniform": }
#'     \code{list(location_l = 1, location_u = 2, scale = 1)}, where scale is the
#'     shared scale parameter of the pareto distributed upper and lower boundaries and location l and u
#'     are the location parameters for the Lower (l) and Upper (u) boundaries of the uniform
#'     distribution. Note this uses the same terminology for the pareto distribution's parameters
#'     as the "uniform" method.
#'     }
#'     \item{\strong{"bivariate_gaussian" and "bivariate_lognormal": }
#'     \code{list(mu = 0, sd = 10, a = 1, b = 1)}, where mu and sd
#'     are the mean and standard deviation of the Normal distribution of the data's mean and a and b
#'     are the shape and scale of the gamma distribution on precision. Note that internally this uses
#'     the Mu and Precision parameterization of the normal distribution and those are the parameters
#'     shown in the plot and tested, but priors use Mu and SD for the normal distribution of the mean.
#'     }
#' }
#'
#' @examples
#' mv_ln <- mvSim(
#'   dists = list(
#'     rlnorm = list(meanlog = log(130), sdlog = log(1.2)),
#'     rlnorm = list(meanlog = log(100), sdlog = log(1.3))
#'   ),
#'   n_samples = 30
#' )
#'
#' # lognormal mv
#' ln_mv_ex <- conjugate(
#'   s1 = mv_ln[1:30, -1], s2 = mv_ln[31:60, -1], method = "lognormal",
#'   priors = list(mu = 5, sd = 2),
#'   rope_range = c(-40, 40), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal", support = NULL
#' )
#'
#' # lognormal sv
#' ln_sv_ex <- conjugate(
#'   s1 = rlnorm(100, log(130), log(1.3)), s2 = rlnorm(100, log(100), log(1.6)),
#'   method = "lognormal",
#'   priors = list(mu = 5, sd = 2),
#'   rope_range = NULL, rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal", support = NULL
#' )
#'
#' # Z test mv example
#'
#' mv_gauss <- mvSim(
#'   dists = list(
#'     rnorm = list(mean = 50, sd = 10),
#'     rnorm = list(mean = 60, sd = 12)
#'   ),
#'   n_samples = 30
#' )
#'
#' gauss_mv_ex <- conjugate(
#'   s1 = mv_gauss[1:30, -1], s2 = mv_gauss[31:60, -1], method = "gaussian",
#'   priors = list(mu = 30, sd = 10),
#'   rope_range = c(-25, 25), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal", support = NULL
#' )
#'
#' # T test sv example with two different priors
#'
#' gaussianMeans_sv_ex <- conjugate(
#'   s1 = rnorm(10, 50, 10), s2 = rnorm(10, 60, 12), method = "t",
#'   priors = list(list(mu = 40, sd = 10), list(mu = 45, sd = 8)),
#'   rope_range = c(-5, 8), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal", support = NULL
#' )
#'
#' # beta mv example
#'
#' set.seed(123)
#' mv_beta <- mvSim(
#'   dists = list(
#'     rbeta = list(shape1 = 5, shape2 = 8),
#'     rbeta = list(shape1 = 10, shape2 = 10)
#'   ),
#'   n_samples = c(30, 20)
#' )
#'
#' beta_mv_ex <- conjugate(
#'   s1 = mv_beta[1:30, -1], s2 = mv_beta[31:50, -1], method = "beta",
#'   priors = list(a = 0.5, b = 0.5),
#'   rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal",
#'   bayes_factor = 0.5 # note this may not be reasonable with these priors
#' )
#'
#' # beta sv example
#'
#' beta_sv_ex <- conjugate(
#'   s1 = rbeta(20, 5, 5), s2 = rbeta(20, 8, 5), method = "beta",
#'   priors = list(a = 0.5, b = 0.5),
#'   rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal",
#'   bayes_factor = c(0.5, 0.75) # note this may not be reasonable with these priors
#' )
#'
#' # binomial sv example
#' # note that specifying trials = 20 would also work
#' # and the number of trials will be recycled to the length of successes
#'
#' binomial_sv_ex <- conjugate(
#'   s1 = list(successes = c(15, 14, 16, 11), trials = c(20, 20, 20, 20)),
#'   s2 = list(successes = c(7, 8, 10, 5), trials = c(20, 20, 20, 20)), method = "binomial",
#'   priors = list(a = 0.5, b = 0.5),
#'   rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal"
#' )
#'
#' # poisson sv example
#'
#' poisson_sv_ex <- conjugate(
#'   s1 = rpois(20, 10), s2 = rpois(20, 8), method = "poisson",
#'   priors = list(a = 0.5, b = 0.5),
#'   rope_range = c(-1, 1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal"
#' )
#'
#' # negative binomial sv example
#' # knowing r (required number of successes) is an important caveat for this method.
#' # in the current implementation we suggest using the poisson method for data such as leaf counts
#'
#' negbin_sv_ex <- conjugate(
#'   s1 = rnbinom(20, 10, 0.5), s2 = rnbinom(20, 10, 0.25), method = "negbin",
#'   priors = list(r = 10, a = 0.5, b = 0.5),
#'   rope_range = c(-1, 1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal"
#' )
#'
#' # von mises mv example
#'
#' mv_gauss <- mvSim(
#'   dists = list(
#'     rnorm = list(mean = 50, sd = 10),
#'     rnorm = list(mean = 60, sd = 12)
#'   ),
#'   n_samples = c(30, 40)
#' )
#' vm1_ex <- conjugate(
#'   s1 = mv_gauss[1:30, -1],
#'   s2 = mv_gauss[31:70, -1],
#'   method = "vonmises",
#'   priors = list(mu = 45, kappa = 1, boundary = c(0, 180), known_kappa = 1, n = 1),
#'   rope_range = c(-1, 1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal"
#' )
#'
#' # von mises 2 sv example
#' vm2_ex <- conjugate(
#'   s1 = brms::rvon_mises(10, 2, 2),
#'   s2 = brms::rvon_mises(15, 3, 3),
#'   method = "vonmises2",
#'   priors = list(mu = 0, kappa = 0.5, boundary = c(-pi, pi), n = 1),
#'   cred.int.level = 0.95
#' )
#'
#' @return
#'
#' A \link{conjugate-class} object with slots including:
#' \itemize{
#'    \item{\strong{summary}: A data frame containing HDI/HDE values for each sample and
#'    the ROPE as well as posterior probability of the hypothesis and ROPE test (if specified).
#'    The HDE is the "Highest Density
#'    Estimate" of the posterior, that is the tallest part of the probability density function. The
#'    HDI is the Highest Density Interval, which is an interval that contains X\% of the posterior
#'    distribution, so \code{cred.int.level = 0.8} corresponds to an HDI that includes 80 percent
#'    of the posterior probability. Bayes factors are calculated as posterior/prior for each sample.}
#'    \item{\strong{posterior}: A list of updated parameters in the same format as the prior
#'     for the given method. If desired this does allow for Bayesian updating.}
#'    \item{\strong{prior}: The prior in a list with the same format as the posterior.}
#'    \item{\strong{plot}: A ggplot showing the distribution of samples and optionally the
#'    distribution of differences/ROPE}.
#'    \item{\strong{plot_parameters}}: Parameters used in making a plot of the data.
#'    Contains support range and posterior recoded to use a density function.
#'    \item{\strong{data}}: Data from s1 and s2 arguments.
#'    \item{\strong{call}}: The function call.
#' }
#' @seealso \link{barg} for additional reporting
#' @keywords bayesian conjugate priors ROPE
#' @export

conjugate <- function(s1 = NULL, s2 = NULL,
                      method = c(
                        "t", "gaussian", "beta", "binomial",
                        "lognormal", "lognormal2", "poisson", "negbin", "vonmises", "vonmises2",
                        "uniform", "pareto", "gamma", "bernoulli", "exponential", "bivariate_uniform",
                        "bivariate_gaussian", "bivariate_lognormal"
                      ),
                      priors = NULL, plot = NULL, rope_range = NULL,
                      rope_ci = 0.89, cred.int.level = 0.89, hypothesis = "equal",
                      bayes_factor = NULL, support = NULL) {
  #* `Handle formula option in s1`
  samples <- .formatSamples(s1, s2)
  s1 <- samples$s1
  s2 <- samples$s2
  #* `check length of method, replicate if there is a second sample`
  if (length(method) == 1 && !is.null(s2)) {
    method <- rep(method, 2)
  }
  if (!is.null(priors) && !methods::is(priors[[1]], "list")) {
    priors <- list(priors, priors)
  }
  samplesList <- list(s1)
  if (!is.null(s2)) {
    samplesList[[2]] <- s2
  }

  if (!missing("support")) {
    warning("support argument is deprecated.")
  }
  if (!missing("plot")) {
    warning("plot argument is deprecated, use plot.conjugate instead.")
  }
  support <- .getSupport(samplesList, method, priors) # calculate shared support

  sample_results <- lapply(seq_along(samplesList), function(i) {
    sample <- samplesList[[i]]
    prior <- priors[[i]]
    #* `Check sample class`
    if (is.matrix(sample) | is.data.frame(sample)) {
      vec_suffix <- "mv"
      sample <- .mvSampleFormatting(sample)
    } else if (is.vector(sample)) {
      vec_suffix <- "sv"
    }
    matched_arg <- match.arg(method[i], choices = c(
      "t", "gaussian", "beta", "binomial",
      "lognormal", "lognormal2", "poisson", "negbin",
      "vonmises", "vonmises2",
      "uniform", "pareto", "gamma", "bernoulli", "exponential",
      "bivariate_uniform", "bivariate_gaussian", "bivariate_lognormal"
    ))
    matched_fun <- get(paste0(".conj_", matched_arg, "_", vec_suffix))
    res <- matched_fun(sample, prior, support, cred.int.level)
    return(res)
  })
  #* `combine results into an object to return`
  out <- list()
  out$summary <- do.call(cbind, lapply(seq_along(sample_results), function(i) {
    s <- sample_results[[i]]$summary
    if (!is.null(bayes_factor)) { #* `Calculate Bayes Factors`
      bf <- .conj_bayes_factor(bayes_factor, sample_results[[i]])
      s$bf_1 <- bf
    }
    if (i == 2) {
      s <- s[, !grepl("param", colnames(s))]
      colnames(s) <- gsub("1", "2", colnames(s))
    }
    return(s)
  }))
  if (!is.null(s2)) {
    postProbRes <- .pdf.handling(sample_results[[1]]$pdf, sample_results[[2]]$pdf, hypothesis)
    out$summary <- cbind(
      out$summary,
      data.frame("hyp" = hypothesis, "post.prob" = as.numeric(postProbRes$post.prob))
    )
  }
  #* `parse output and do ROPE`
  if (!is.null(rope_range)) {
    rope_res <- .conj_rope(sample_results, rope_range, rope_ci, method)
    out$summary <- cbind(out$summary, rope_res$summary)
    out$rope_df <- rope_res$rope_df
  }
  out$posterior <- lapply(sample_results, function(s) {return(s$posterior)})
  out$prior <- lapply(sample_results, function(s) {return(s$prior)})
  out$plot_parameters <- lapply(sample_results, function(s) {return(s$plot_list)})
  #* `Special bivariate distribution things`
  if (any(grepl("bivariate", method))) {
    out$plot_df <- lapply(sample_results, function(s) s$plot_df)
    out$posterior_draws <- lapply(sample_results, function(s) s$posteriorDraws)
  }
  #* `Save Data and Call`
  out$data <- list(s1, s2)[!unlist(lapply(list(s1, s2), is.null))]
  out$call <- match.call()
  #* `Change class to conjugate`
  out <- as.conjugate(out)
  return(out)
}

#' ***********************************************************************************************
#' *************** `Formula Handling Helper function` ***********************************
#' ***********************************************************************************************

#' @keywords internal
#' @noRd

.formatSamples <- function(s1 = NULL, s2 = NULL) {
  if (methods::is(s1, "formula")) {
    if (!is.data.frame(s2)) {
      stop("If s1 is a formula then s2 must be a data.frame")
    }
    rhs <- as.character(s1)[3]
    lhs <- as.character(s1)[2]
    if (lhs %in% colnames(s2)) { # handle SV traits
      samples <- split(s2[[lhs]], s2[[rhs]])
    } else { # handle MV traits
      samples <- lapply(split(s2, s2[[rhs]]), function(d) {
        return(d[, eval(str2lang(lhs))])
      })
    }
    names(samples) <- c("s1", "s2")
  } else {
    samples <- list("s1" = s1, "s2" = s2)
  }
  samples <- .check_formatted_samples(samples)
  return(samples)
}

#' @keywords internal
#' @noRd

.check_formatted_samples <- function(samples) {
  any_one_l <- any(
    unlist(
      lapply(
        samples[!unlist(lapply(samples, is.null))],
        function(s) {
          return(length(t(s)) == 1)
        }
      )
    )
  )
  if (any_one_l) {
    stop("Conjugate samples must be >1L", call. = FALSE)
  }
  return(samples)
}

#' ***********************************************************************************************
#' *************** `MV Sample Formatting Helper function` ***********************************
#' ***********************************************************************************************

#' @keywords internal
#' @noRd

.mvSampleFormatting <- function(sample) {
  #* `Standardize sample class and names`
  if (is.matrix(sample)) {
    original_names <- colnames(sample)
    sample <- as.data.frame(sample)
    colnames(sample) <- original_names
  }
  if (is.null(colnames(sample))) {
    bins <- (seq_along(sample))
    colnames(sample) <- paste0("b", bins)
    warning(paste0("Assuming unnamed columns represent bins from ", min(bins), " to ", max(bins)))
  }
  return(sample)
}

#' ***********************************************************************************************
#' *************** `Support Calculating function` ***********************************
#' ***********************************************************************************************
#'
#'
#' @keywords internal
#' @noRd


.getSupport <- function(samplesList, method, priors) {
  #* `Check for bivarate`
  if (any(grepl("bivariate", method[1]))) {
    biv <- TRUE
  } else {
    biv <- FALSE
  }
  support_quantiles <- lapply(seq_along(samplesList), function(i) {
    sample <- samplesList[[i]]
    prior <- priors[[i]]
    #* `Check sample class`
    if (is.matrix(sample) | is.data.frame(sample)) {
      vec <- FALSE
    } else if (is.vector(sample)) {
      vec <- TRUE
    }
    vec_suffix <- if (vec) {
      "sv"
    } else {
      "mv"
    }
    matched_fun <- get(paste0(".conj_", method[i], "_", vec_suffix))
    qnts <- matched_fun(s1 = sample, priors = prior, calculatingSupport = TRUE)
    return(qnts)
  })
  if (biv) {
    pars <- names(support_quantiles[[1]])
    support <- lapply(pars, function(param) {
      qnts <- range(unlist(lapply(support_quantiles, function(sq) {
        return(sq[[param]])
      })), na.rm = TRUE)
      return(seq(qnts[1], qnts[2], length.out = 10000))
    })
    names(support) <- pars
  } else {
    qnts <- range(unlist(support_quantiles), na.rm = TRUE)
    support <- seq(qnts[1], qnts[2], length.out = 10000)
  }
  return(support)
}



#' ***********************************************************************************************
#' *************** `ROPE testing on two conjugateHelper outputs` ***********************************
#' ***********************************************************************************************
#'
#' this should take outputs from conjHelpers and compare the $posteriorDraws.
#' @keywords internal
#' @noRd
.conj_rope <- function(sample_results, rope_range = c(-0.1, 0.1),
                       rope_ci = 0.89, method) {
  #* `if bivariate then call the bivariate option`
  #* note this will return to .conj_rope but with a non-bivariate method
  if (any(grepl("bivariate", method))) {
    rope_res <- .conj_bivariate_rope(
      sample_results, rope_range,
      rope_ci, method
    )
    return(rope_res)
  }
  #* `ROPE Comparison`
  rope_res <- list()
  if (!is.null(rope_range)) {
    if (length(rope_range) == 2) {
      post1 <- sample_results[[1]]$posteriorDraws
      if (length(sample_results) == 2) {
        post2 <- sample_results[[2]]$posteriorDraws
        if (any(grepl("vonmises", method))) {
          boundary <- sample_results[[1]]$posterior$boundary
          posterior <- .conj_rope_circular_diff(post1, post2, boundary = boundary)
        } else {
          posterior <- post1 - post2
        }
      } else {
        posterior <- post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior,
        range = rope_range,
        ci_method = "HDI", ci = rope_ci
      ))
      rope_test <- data.frame(
        HDE_rope = hde_diff, HDI_rope_low = hdi_diff[1],
        HDI_rope_high = hdi_diff[2], rope_prob = rope_prob
      )
      rope_res$summary <- rope_test
      rope_res$rope_df <- data.frame("X" = posterior)
    } else {
      stop("rope must be a vector of length 2")
    }
  }
  return(rope_res)
}

#' ***********************************************************************************************
#' *************** `ROPE testing on two conjugate_bivariate_Helper outputs` ********************
#' ***********************************************************************************************
#'
#' this should take outputs from conjHelpers and compare the $posteriorDraws.
#' @keywords internal
#' @noRd
.conj_bivariate_rope <- function(sample_results, rope_range = c(-0.1, 0.1),
                                 rope_ci = 0.89, method) {
  #* `Format rope_range`
  rope_res <- list()
  if (!is.list(rope_range)) {
    rope_range <- lapply(sample_results[[1]]$summary$param, function(par) {
      return(rope_range)
    })
    names(rope_range) <- sample_results[[1]]$summary$param
  }
  #* `ROPE Comparison`
  rope_res <- lapply(names(rope_range), function(nm) {
    iter_rope_range <- rope_range[[nm]]
    sample_results_param <- lapply(sample_results, function(s_res) {
      s_res$posteriorDraws <- s_res$posteriorDraws[[nm]]
      s_res$pdf <- s_res$pdf[[nm]]
      s_res$summary <- s_res$summary[s_res$summary$param == nm, ]
      s_res$plot_df <- s_res$plot_df[s_res$plot_df$param == nm, ]
      return(s_res)
    })
    nm_res <- .conj_rope(sample_results_param,
      rope_range = iter_rope_range,
      rope_ci = rope_ci, method = "NONE"
    )
    return(nm_res)
  })
  rope_res$summary <- do.call(rbind, lapply(rope_res, function(r) {
    return(r$summary)
  }))
  return(rope_res)
}

#' ***********************************************************************************************
#' *************** `Handle PDFs for testing` ***********************************
#' ***********************************************************************************************
#' @keywords internal
#' @noRd

.pdf.handling <- function(pdf1, pdf2, hypothesis) {
  if (is.list(pdf1) && is.list(pdf2)) {
    pdf.handling.output <- as.data.frame(do.call(rbind, lapply(seq_along(pdf1), function(i) {
      pdf <- .post.prob.from.pdfs(pdf1[[i]], pdf2[[i]], hypothesis)
      return(pdf)
    })))
  } else {
    pdf.handling.output <- as.data.frame(.post.prob.from.pdfs(pdf1, pdf2, hypothesis))
  }
  return(pdf.handling.output)
}

#' ***********************************************************************************************
#' *************** `Calculate Posterior Probability given PDFs` ***********************************
#' ***********************************************************************************************

#' @keywords internal
#' @noRd
.post.prob.from.pdfs <- function(pdf1, pdf2, hypothesis) {
  if (hypothesis == "unequal") {
    post.prob <- 1 - sum(apply(cbind(pdf1, pdf2), MARGIN = 1, function(i) min(i)), na.rm = TRUE)
  } else if (hypothesis == "equal") {
    post.prob <- sum(apply(cbind(pdf1, pdf2), MARGIN = 1, function(i) min(i)), na.rm = TRUE)
  } else if (hypothesis == "lesser") { # note one sided testing is less tested generally
    direction <- pdf1 <= pdf2
    post.prob <- sum(pdf1 * direction, na.rm = TRUE)
  } else if (hypothesis == "greater") {
    direction <- pdf1 >= pdf2
    post.prob <- sum(pdf1 * direction, na.rm = TRUE)
  } else {
    stop("hypothesis must be either unequal, equal, lesser, or greater")
  }
  return(list("post.prob" = post.prob))
}
