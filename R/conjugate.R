#' Bayesian testing using conjugate priors and method of moments for single or multi value traits.
#'
#' @description
#' Function to perform bayesian tests and ROPE comparisons using single or multi value traits with
#' several distributions.
#'
#' @param s1 A data.frame or matrix of multi value traits or a vector of single value traits.
#' If a multi value trait is used then column names should include a number representing the "bin".
#' Alternatively for distributions other than "binomial" (which requires list data, see examples)
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
#' This should only be done with caution. Make sure you understand the interpretation of any
#' comparison you are doing if you specify two methods (c("gaussian", "lognormal") as an arbitrary
#' example).
#' @param priors Prior distributions described as a list of lists.
#' If this is a single list then it will be duplicated for the second sample,
#' which is generally a good idea if both
#' samples use the same distribution (method).
#' Elements in the inner lists should be named for the parameter they represent (see examples).
#' These names vary by method (see details).
#'  By default this is NULL and weak priors (generally jeffrey's priors) are used.
#'  The \code{posterior} part of output can also be recycled as a new prior if Bayesian
#'  updating is appropriate for your use.
#' @param plot Logical, should a ggplot be made and returned.
#' @param rope_range Optional vector specifying a region of practical equivalence.
#' This interval is considered practically equivalent to no effect.
#' Kruschke (2018) suggests c(-0.1, 0.1) as a broadly reasonable ROPE for standardized parameters.
#' That range could also be rescaled by a standard deviation/magnitude for
#' non-standardized parameters, but ultimately this should be informed by your
#' setting and scientific question.
#' See Kruschke (2018) for details on ROPE and other Bayesian methods to aide
#' decision-making \url{https://doi.org/10.3758/s13423-016-1221-4}
#' and \url{https://doi.org/10.1177/251524591877130}.
#' @param rope_ci The credible interval probability to use for ROPE. Defaults to 0.89.
#' @param cred.int.level The credible interval probability to use
#' in computing HDI for samples, defaults to 0.89.
#' @param hypothesis Direction of a hypothesis if two samples are provided.
#'  Options are "unequal", "equal", "greater", and "lesser",
#'   read as "sample1 greater than sample2".
#' @param support Optional support vector to include all possible values the random variable
#'  (samples) might take. This defaults to NULL in which case each method will use default
#'  behavior to attempt to calculate a dense support, but it is a good idea to supply this
#'  with some suitable vector. For example, the Beta method uses \code{seq(0.0001, 0.9999, 0.0001)}
#'  for support.
#'
#' @import bayestestR
#' @import ggplot2
#' @import patchwork
#' @import extraDistr
#' @importFrom stats var median rnorm dnorm qbeta rbeta dbeta qbeta rlnorm
#' dlnorm qlnorm qt rgamma qgamma dgamma rpois
#'
#' @details
#'
#' Prior distributions default to be weakly informative and in some cases you may wish to change them.
#' \itemize{
#'    \item{\strong{"t" and "gaussian":} \code{priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ) },
#'     where mu is the mean, n is the number of prior observations, and s2 is variance}
#'    \item{\strong{"beta", "bernoulli", and "binomial":}
#'    \code{priors = list( a=c(0.5, 0.5), b=c(0.5, 0.5) )},
#'     where a and b are shape parameters of the beta distribution. Note that for the binomial
#'     distribution this is used as the prior for success probability P,
#'     which is assumed to be beta distributed as in a beta-binomial distribution.}
#'    \item{\strong{"lognormal": } \code{priors = list(mu = 0, sd = 5) },
#'    where mu and sd describe the normal distribution of the mean parameter for lognormal data.
#'    Note that these values are on the log scale.}
#'    \item{\strong{"lognormal2": } \code{priors = list(a = 1, b = 1) },
#'    where a and b are the shape and scale parameters of the gamma distribution of lognormal data's
#'    precision parameter (using the alternative mu, precision paramterization).
#'    }
#'    \item{\strong{"gamma": } \code{priors = list(shape = 0.5, scale = 0.5, known_shape = 1)},
#'     where shape and scale are the respective parameters of the gamma distributed rate
#'     (inverse of scale) parameter of gamma distributed data.}
#'    \item{\strong{"poisson" and "exponential": } \code{priors = list(a=c(0.5,0.5),b=c(0.5,0.5))},
#'     where a and b are shape parameters of the gamma distribution.}
#'    \item{\strong{"negbin": } \code{priors = list(r=c(10,10), a=c(0.5,0.5),b=c(0.5,0.5))},
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
#'     called shape and scale, corresponding to extraDistr's scale and location respecitvely, which
#'     can be confusing. Note that the lower boundary of the uniform is assumed to be constant at 0.
#'     There is a more complicated conjugate distribution to estimate both boundaries but that has
#'     not seemed relevant so far in the author's plant phenotyping experience.
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
#' }
#'
#' See examples for plots of these prior distributions.
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
#'   plot = FALSE, rope_range = c(-40, 40), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal", support = NULL
#' )
#'
#' # lognormal sv
#' ln_sv_ex <- conjugate(
#'   s1 = rlnorm(100, log(130), log(1.3)), s2 = rlnorm(100, log(100), log(1.6)),
#'   method = "lognormal",
#'   priors = list(mu = 5, sd = 2),
#'   plot = FALSE, rope_range = NULL, rope_ci = 0.89,
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
#'   priors = list(mu = 30, n = 1, s2 = 100),
#'   plot = FALSE, rope_range = c(-25, 25), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal", support = NULL
#' )
#'
#' # T test sv example
#'
#' gaussianMeans_sv_ex <- conjugate(
#'   s1 = rnorm(10, 50, 10), s2 = rnorm(10, 60, 12), method = "t",
#'   priors = list(mu = 30, n = 1, s2 = 100),
#'   plot = FALSE, rope_range = c(-5, 8), rope_ci = 0.89,
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
#'   plot = FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal"
#' )
#'
#' # beta sv example
#'
#' beta_sv_ex <- conjugate(
#'   s1 = rbeta(20, 5, 5), s2 = rbeta(20, 8, 5), method = "beta",
#'   priors = list(a = 0.5, b = 0.5),
#'   plot = FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal"
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
#'   plot = FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal"
#' )
#'
#' # poisson sv example
#'
#' poisson_sv_ex <- conjugate(
#'   s1 = rpois(20, 10), s2 = rpois(20, 8), method = "poisson",
#'   priors = list(a = 0.5, b = 0.5),
#'   plot = FALSE, rope_range = c(-1, 1), rope_ci = 0.89,
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
#'   plot = FALSE, rope_range = c(-1, 1), rope_ci = 0.89,
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
#'   plot = FALSE, rope_range = c(-1, 1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal"
#' )
#'
#' # von mises 2 sv example
#' vm2_ex <- conjugate(
#'   s1 = brms::rvon_mises(10, 2, 2),
#'   s2 = brms::rvon_mises(15, 3, 3),
#'   method = "vonmises2",
#'   priors = list(mu = 0, kappa = 0.5, boundary = c(-pi, pi), n = 1),
#'   cred.int.level = 0.95,
#'   plot = FALSE
#' )
#'
#'
#' # Example usage with plantCV data
#' ## Not run:
#' if (FALSE) {
#'   library(data.table)
#'   wide <- read.pcv(
#'     paste0(
#'       "https://media.githubusercontent.com/media/joshqsumner/",
#'       "pcvrTestData/main/pcv4-multi-value-traits.csv"
#'     ),
#'     reader = "fread", mode = "wide"
#'   )
#'
#'   wide$genotype <- substr(wide$barcode, 3, 5)
#'   wide$genotype <- ifelse(wide$genotype == "002", "B73",
#'     ifelse(wide$genotype == "003", "W605S",
#'       ifelse(wide$genotype == "004", "MM", "Mo17")
#'     )
#'   )
#'   wide$fertilizer <- substr(wide$barcode, 8, 8)
#'   wide$fertilizer <- ifelse(wide$fertilizer == "A", "100",
#'     ifelse(wide$fertilizer == "B", "50", "0")
#'   )
#'   wide <- bw.time(wide, timeCol = "timestamp", group = "barcode")
#'
#'   mo17_sample <- wide[
#'     wide$genotype == "Mo17" &
#'       wide$DAS > 18 &
#'       wide$fertilizer == 100,
#'     grepl("hue_freq", colnames(wide))
#'   ]
#'   B73_sample <- wide[
#'     wide$genotype == "B73" &
#'       wide$DAS > 18 &
#'       wide$fertilizer == 100,
#'     grepl("hue_freq", colnames(wide))
#'   ]
#'
#'   hue_res_t <- conjugate(
#'     s1 = mo17_sample, s2 = B73_sample, method = "t",
#'     priors = list(mu = 50, n = 1, s2 = 100),
#'     plot = TRUE, rope_range = c(-10, 10), rope_ci = 0.89,
#'     cred.int.level = 0.89, hypothesis = "equal"
#'   )
#'
#'   hue_res_ln <- conjugate(
#'     s1 = mo17_sample, s2 = B73_sample, method = "lognormal",
#'     plot = TRUE, rope_range = c(-10, 10), rope_ci = 0.89,
#'     cred.int.level = 0.89, hypothesis = "equal"
#'   )
#'   # Picking the right distribution makes a difference,
#'   # checking plots of your data (see ?pcv.joyplot for mv traits)
#'   # will be useful.
#'
#'   sv <- read.pcv(
#'     paste0(
#'       "https://raw.githubusercontent.com/joshqsumner/",
#'       "pcvrTestData/main/pcv4-single-value-traits.csv"
#'     ),
#'     reader = "fread", mode = "wide"
#'   )
#'
#'   sv$genotype <- substr(sv$barcode, 3, 5)
#'   sv$genotype <- ifelse(sv$genotype == "002", "B73",
#'     ifelse(sv$genotype == "003", "W605S",
#'       ifelse(sv$genotype == "004", "MM", "Mo17")
#'     )
#'   )
#'   sv$fertilizer <- substr(sv$barcode, 8, 8)
#'   sv$fertilizer <- ifelse(sv$fertilizer == "A", "100",
#'     ifelse(sv$fertilizer == "B", "50", "0")
#'   )
#'
#'   sv <- bw.time(sv, timeCol = "timestamp", group = "barcode")
#'
#'   pixels_per_cmsq <- 42.5^2 # pixel per cm^2
#'   sv$area_cm2 <- sv$area_pixels / pixels_per_cmsq
#'
#'   mo17_area <- sv[wide$genotype == "Mo17" & wide$DAS > 18 & wide$fertilizer == 50, "area_cm2"]
#'   B73_area <- sv[wide$genotype == "B73" & wide$DAS > 18 & wide$fertilizer == 50, "area_cm2"]
#'
#'   area_res_t <- conjugate(
#'     s1 = mo17_area, s2 = B73_area, method = "t",
#'     priors = list(mu = 30, n = 1, s2 = 100),
#'     plot = TRUE, rope_range = c(-5, 5), rope_ci = 0.89,
#'     cred.int.level = 0.89, hypothesis = "equal"
#'   )
#' }
#' ## End(Not run)
#'
#' # Plots of prior distributions
#' set.seed(123)
#' plot(seq(0, 1, 0.0001), dbeta(seq(0, 1, 0.0001), 0.5, 0.5),
#'   ylab = "Density",
#'   xlab = "Support", main = "Beta", type = "l"
#' )
#'
#' plot(seq(0, 10, 0.01), dgamma(seq(0, 10, 0.01), 0.5, 0.5),
#'   ylab = "Density", xlab = "Support",
#'   main = "Poisson (gamma prior on Lambda parameter)", type = "l"
#' )
#'
#' plot(seq(-20, 20.001), dnorm(seq(-20, 20.001), 0, sqrt(100)),
#'   ylab = "Density",
#'   xlab = "Support", main = "T and Gaussian", type = "l"
#' )
#'
#' plot(seq(0, 70, 0.001), dlnorm(seq(0, 70, 0.001), log(10), log(3)),
#'   ylab = "Density",
#'   xlab = "Support", main = "Lognormal", type = "l"
#' )
#'
#' @return
#'
#' A list with named elements:
#' \itemize{
#'    \item{\strong{summary}: A data frame containing HDI/HDE values for each sample and
#'    the ROPE as well as posterior probability of the hypothesis.}
#'    \item{\strong{posterior}: A list of updated parameters in the same format as the prior
#'     for the given method. If desired this does allow for Bayesian updating.}
#'    \item{\strong{plot_df}: A data frame of probabilities along the support for each sample.
#'     This is used for making the ggplot.}
#'    \item{\strong{rope_df}: A data frame of draws from the ROPE posterior.}
#'    \item{\strong{plot}: A ggplot showing the distribution of samples and optionally the
#'    distribution of differences/ROPE}
#' }
#'
#' @keywords bayesian, conjugate, ROPE
#' @export

conjugate <- function(s1 = NULL, s2 = NULL,
                      method = c(
                        "t", "gaussian", "beta", "binomial",
                        "lognormal", "lognormal2", "poisson", "negbin", "vonmises", "vonmises2",
                        "uniform", "pareto", "gamma", "bernoulli", "exponential"
                      ),
                      priors = NULL, plot = FALSE, rope_range = NULL,
                      rope_ci = 0.89, cred.int.level = 0.89, hypothesis = "equal",
                      support = NULL) {
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

  if (is.null(support)) {
    support <- .getSupport(samplesList, method, priors) # calculate shared support
  }

  sample_results <- lapply(seq_along(samplesList), function(i) {
    sample <- samplesList[[i]]
    prior <- priors[[i]]
    #* `Check sample class`
    if (is.matrix(sample) | is.data.frame(sample)) {
      vec <- FALSE
    } else if (is.vector(sample)) {
      vec <- TRUE
    } else {
      stop("samples must be a vector, data.frame, or matrix.")
    }
    matched_arg <- match.arg(method[i], choices = c(
      "t", "gaussian", "beta", "binomial",
      "lognormal", "lognormal2", "poisson", "negbin",
      "vonmises", "vonmises2",
      "uniform", "pareto", "gamma", "bernoulli", "exponential"
    ))
    # turning off dirichlet until I decide on a new implementation that I like better
    # and a use case that isn't so ripe for abuse.
    vec_suffix <- if (vec) {
      "sv"
    } else {
      "mv"
    }
    matched_fun <- get(paste0(".conj_", matched_arg, "_", vec_suffix))
    res <- matched_fun(sample, prior, plot, support, cred.int.level)
    return(res)
  })

  #* `combine results into an object to return`
  out <- list()
  out$summary <- do.call(cbind, lapply(sample_results, function(s) {
    s$summary
  }))
  if (!is.null(s2)) {
    colnames(out$summary) <- c("HDE_1", "HDI_1_low", "HDI_1_high", "HDE_2", "HDI_2_low", "HDI_2_high")
    postProbRes <- .post.prob.from.pdfs(sample_results[[1]]$pdf, sample_results[[2]]$pdf, hypothesis)
    out$summary <- cbind(
      out$summary,
      data.frame("hyp" = hypothesis, "post.prob" = postProbRes$post.prob)
    )
    dirSymbol <- postProbRes$direction
  } else {
    dirSymbol <- NULL
  }
  #* `parse output and do ROPE`
  if (!is.null(rope_range)) {
    rope_res <- .conj_rope(sample_results, rope_range, rope_ci, plot, method)
    out$summary <- cbind(out$summary, rope_res$summary)
  } else {
    rope_res <- NULL
  }
  out$posterior <- lapply(sample_results, function(s) s$posterior)

  #* `Make plot`
  if (plot) {
    out$plot <- .conj_general_plot(sample_results, rope_res,
      res = out,
      rope_range, rope_ci, dirSymbol, support
    )
  }
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
        d[, eval(str2lang(lhs))]
      })
    }
    names(samples) <- c("s1", "s2")
    return(samples)
  } else {
    return(list("s1" = s1, "s2" = s2))
  }
}


#' ***********************************************************************************************
#' *************** `Support Calculating function` ***********************************
#' ***********************************************************************************************
#'
#'
#' @keywords internal
#' @noRd


.getSupport <- function(samplesList, method, priors) {
  support_quantiles <- lapply(seq_along(samplesList), function(i) {
    sample <- samplesList[[i]]
    prior <- priors[[i]]
    #* `Check sample class`
    if (is.matrix(sample) | is.data.frame(sample)) {
      vec <- FALSE
    } else if (is.vector(sample)) {
      vec <- TRUE
    } else {
      stop("samples must be a vector, data.frame, or matrix.")
    }

    matched_arg <- match.arg(method[i], choices = c(
      "t", "gaussian", "beta", "binomial",
      "lognormal", "lognormal2", "poisson", "negbin",
      "vonmises", "vonmises2",
      "uniform", "pareto", "gamma", "bernoulli", "exponential"
    ))
    vec_suffix <- if (vec) {
      "sv"
    } else {
      "mv"
    }
    matched_fun <- get(paste0(".conj_", matched_arg, "_", vec_suffix))
    qnts <- matched_fun(s1 = sample, priors = prior, calculatingSupport = TRUE)
    return(qnts)
  })
  qnts <- range(unlist(support_quantiles))
  support <- seq(qnts[1], qnts[2], length.out = 10000)
  return(support)
}



#' ***********************************************************************************************
#' *************** `ROPE testing on two conjugateHelper outputs` ***********************************
#' ***********************************************************************************************
#'
#' this should take outputs from conjHelpers and compare the $posteriorDraws.
#' This requires me to decide how I want things to work between the loop over methods to here.
#' @keywords internal
#' @noRd
.conj_rope <- function(sample_results, rope_range = c(-0.1, 0.1),
                       rope_ci = 0.89, plot, method) {
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
      if (plot) {
        rope_res$rope_df <- data.frame("X" = posterior)
      }
    } else {
      stop("rope must be a vector of length 2")
    }
  }
  return(rope_res)
}



#' ***********************************************************************************************
#' *************** `Calculate Posterior Probability given PDFs` ***********************************
#' ***********************************************************************************************


#' @keywords internal
#' @noRd
.post.prob.from.pdfs <- function(pdf1, pdf2, hypothesis) {
  if (hypothesis == "unequal") {
    post.prob <- 1 - sum(apply(cbind(pdf1, pdf2), MARGIN = 1, function(i) min(i)), na.rm = TRUE)
    dirSymbol <- "!="
  } else if (hypothesis == "equal") {
    post.prob <- sum(apply(cbind(pdf1, pdf2), MARGIN = 1, function(i) min(i)), na.rm = TRUE)
    dirSymbol <- "="
  } else if (hypothesis == "lesser") { # note one sided testing is less tested generally
    direction <- pdf1 <= pdf2
    post.prob <- sum(pdf1 * direction, na.rm = TRUE)
    dirSymbol <- "<"
  } else if (hypothesis == "greater") {
    direction <- pdf1 >= pdf2
    post.prob <- sum(pdf1 * direction, na.rm = TRUE)
    dirSymbol <- ">"
  } else {
    stop("hypothesis must be either unequal, equal, lesser, or greater")
  }
  return(list("post.prob" = post.prob, "direction" = dirSymbol))
}

#' ***********************************************************************************************
#' *************** `General Plotting Function` ***********************************
#' ***********************************************************************************************


#' @keywords internal
#' @noRd
.conj_general_plot <- function(sample_results, rope_res = NULL, res,
                               rope_range, rope_ci, dirSymbol = NULL, support) {
  s1_plot_df <- sample_results[[1]]$plot_df

  p <- ggplot2::ggplot(s1_plot_df, ggplot2::aes(x = .data$range, y = .data$prob)) +
    ggplot2::geom_area(data = s1_plot_df, alpha = 0.5, ggplot2::aes(fill = "s1")) +
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
        round(res$summary$HDI_1_high, 2), "]"
      )
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 0.5))) +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "inside",
                   legend.position.inside = c(0.9, 0.9)) +
    pcv_theme()

  if (length(sample_results) == 2) {
    s2_plot_df <- sample_results[[2]]$plot_df

    if (res$summary$post.prob < 1e-5) {
      post.prob.text <- "<1e-5"
    } else {
      post.prob.text <- round(res$summary$post.prob, 5)
    }

    fill_scale <- which(sapply(p$scales$scales, function(x) {
      "fill" %in% x$aesthetics # avoid "replacing scale" non-messages that suppress doesn't catch
    }))

    p$scales$scales[[fill_scale]] <- NULL
    p <- p +
      ggplot2::geom_area(data = s2_plot_df, ggplot2::aes(fill = "s2"), alpha = 0.5) +
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
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = c(0.9, 0.9)) +
      ggplot2::labs(subtitle = paste0(
        "Sample 1:  ", round(res$summary$HDE_1, 2), " [", round(res$summary$HDI_1_low, 2), ", ",
        round(res$summary$HDI_1_high, 2), "]\n",
        "Sample 2:  ", round(res$summary$HDE_2, 2), " [", round(res$summary$HDI_2_low, 2), ", ",
        round(res$summary$HDI_2_high, 2), "]\n",
        "P[p1", dirSymbol, "p2] = ", post.prob.text
      ))
  }
  #* `make x axis range if using default support`
  if (is.null(support)) {
    x_lower <- min(unlist(lapply(sample_results, function(s) {
      min(s$plot_df$range)
    }))) / 1.1
    x_upper <- max(unlist(lapply(sample_results, function(s) {
      max(s$plot_df$range)
    }))) * 1.1
    p <- p + ggplot2::coord_cartesian(xlim = c(x_lower, x_upper))
  }

  if (!is.null(rope_res)) {
    rdf <- rope_res$rope_df
    p <- p + ggplot2::ggplot(rdf, ggplot2::aes(x = .data$X)) +
      ggplot2::geom_histogram(bins = 100, fill = "purple", color = "purple", alpha = 0.7) +
      ggplot2::geom_histogram(
        data = data.frame(
          "X" = rdf[rdf$X > rope_range[1] & rdf$X < rope_range[2] &
                      rdf$X > res$summary$HDI_rope_low &
                      rdf$X < res$summary$HDI_rope_high, ]
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
          rope_ci, "% HDI in [", rope_range[1], ", ", rope_range[2], "]: ",
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
