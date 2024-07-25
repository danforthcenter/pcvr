#' @description
#' Internal function for calculating \mu and \kappa of a distribution represented by single value
#' traits.
#' @param s1 A vector of numerics generated from a circular process
#' @examples
#' .conj_vonmises2_sv(
#'   s1 = brms::rvon_mises(100, 2, 2), priors = list(mu = 0.5, kappa = 0.5),
#'   cred.int.level = 0.95,
#'   plot = FALSE
#' )
#' .conj_vonmises2_sv(
#'   s1 = rnorm(20, 90, 20),
#'   priors = list(mu = 75, kappa = 0.5, boundary = c(0, 180)),
#'   cred.int.level = 0.95,
#'   plot = TRUE
#' )
#' @keywords internal
#' @noRd

.conj_vonmises2_sv <- function(s1 = NULL, priors = NULL,
                               plot = FALSE, support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  #* `set support to NULL to avoid default length of 10000`
  support <- NULL
  #* `make default prior if none provided`
  default_prior <- list(
    mu = 0, kappa = 0.5,
    boundary = c(-pi, pi),
    n = 1
  )
  if (is.null(priors)) {
    priors <- default_prior
  }
  #* `if any elements are missing from prior then use defaults`
  priors <- stats::setNames(lapply(names(default_prior), function(nm) {
    if (nm %in% names(priors)) {
      return(priors[[nm]])
    } else {
      return(default_prior[[nm]])
    }
  }), names(default_prior))
  #* `rescale data to [-pi, pi] according to boundary`
  s1 <- .boundary.to.radians(x = s1, boundary = priors$boundary)
  #* `rescale prior on mu to [-pi, pi] according to boundary`
  mu_radians <- .boundary.to.radians(x = priors$mu, boundary = priors$boundary)
  #* `Raise error if the boundary is wrong and data is not on [-pi, pi]`
  if (any(abs(s1) > pi)) {
    stop(paste0(
      "Values must be on [-pi, pi] after rescaling. ",
      "Does the boundary element in your prior include all your data?"
    ))
  }
  #* `Define dense Support`
  if (is.null(support)) {
    if (calculatingSupport) {
      return(priors$boundary) #* this would be [-pi, pi] if using radians, but plotting will be on
      #* the original scale so we can just return the boundary and use [-pi, pi] as support here
    }
    support_boundary <- seq(min(priors$boundary), max(priors$boundary), by = 0.0005)
    support <- seq(-pi, pi, length.out = length(support_boundary))
  }
  out <- list()
  #* ***** `Updating Kappa`
  n1 <- length(s1)
  obs_kappa <- .unbiased.kappa(s1, n1)
  kappa_prime <- ((obs_kappa * n1) + (priors$kappa * priors$n)) / (n1 + priors$n)
  #* ***** `Updating vMF for mu using kappa prime`
  #* `Get weighted mean of data and prior for half tangent adjustment`
  cm <- .circular.mean(c(s1, mu_radians), w = c(rep(1, length(s1)), priors$n))
  unitCircleAdj <- ifelse(abs(cm) <= pi / 2, 0, pi)
  unitCircleAdj <- ifelse(cm > 0, 1, -1) * unitCircleAdj
  #* `Update prior parameters`
  a <- priors$kappa # kappa parameter can easily overwhelm mean in updating with low sample size
  # I do not love this currently. It seems like a should be kappa prime, but that very easily
  # overwhelms small samples in the follow formula to update mu, so instead of the sequential
  # updating I am using separate updating.
  # the formula below basically weighs the prior by kappa, so if I update kappa then weigh the prior
  # by that then it is much more biased towards the prior, I think that is worth avoiding with this
  # workaround
  b <- mu_radians
  mu_prime_atan_scale <- atan(((a * sin(b)) + sum(sin(s1))) / ((a * cos(b)) + sum(cos(s1))))
  mu_prime <- unitCircleAdj + mu_prime_atan_scale
  #* `calculate density over support`
  dens1 <- brms::dvon_mises(support, mu_prime, kappa_prime)
  pdf1 <- dens1 / sum(dens1)
  #* `calculate highest density interval`
  #* note there is no qvon_mises function, so I am using bayestestR::hdi on
  #* posterior draws and rescaled posterior draws
  draws <- brms::rvon_mises(10000, mu_prime, kappa_prime)
  hdi_v1 <- as.numeric(bayestestR::hdi(draws, ci = cred.int.level))[2:3]
  draws2 <- draws
  draws2[draws2 < 0] <- draws2[draws2 < 0] + 2 * pi
  hdi_v2 <- as.numeric(bayestestR::hdi(draws2, ci = cred.int.level))[2:3]
  hdis <- list(hdi_v1, hdi_v2)
  hdi <- hdis[[which.min(c(diff(hdi_v1), diff(hdi_v2)))]]
  hdi[hdi > pi] <- hdi[hdi > pi] - (2 * pi) # if the second hdi was narrower then fix the part beyond pi
  #* `store highest density estimate`
  hde <- mu_prime
  #* `Rescale HDI, HDE, and draws, from radians to boundary units`
  hdi_boundary <- .radians.to.boundary(hdi, target = priors$boundary)
  hde_boundary <- .radians.to.boundary(hde, target = priors$boundary)
  draws_boundary <- .radians.to.boundary(draws, target = priors$boundary)
  #* `save summary and parameters`
  out$summary <- data.frame(
    HDE_1 = hde_boundary,
    HDI_1_low = hdi_boundary[1],
    HDI_1_high = hdi_boundary[2]
  )
  out$posterior$mu <- hde_boundary # rescaled mu_prime
  out$posterior$kappa <- kappa_prime
  out$posterior$n <- priors$n + length(s1)
  out$posterior$boundary <- priors$boundary
  #* `Store Posterior Draws`
  out$posteriorDraws <- draws_boundary
  out$pdf <- pdf1
  #* `keep data for plotting`
  if (plot) {
    out$plot_df <- data.frame(
      "range" = support_boundary, "prob" = pdf1,
      "sample" = rep("Sample 1", length(support_boundary))
    )
  } # tests on this seem to work fine
  return(out)
}

#' @description
#' Internal function for calculating \mu and \kappa of a distribution represented by multi value
#' traits.
#' @param s1 A vector of numerics generated from a circular process
#' @examples
#' mv_gauss <- mvSim(
#'   dists = list(
#'     rnorm = list(mean = 50, sd = 10)
#'   ),
#'   n_samples = 30
#' )
#' .conj_vonmises2_mv(
#'   s1 = mv_gauss[, -1], priors = list(mu = 30, kappa = 1, boundary = c(0, 180)),
#'   cred.int.level = 0.95,
#'   plot = TRUE
#' )
#' @keywords internal
#' @noRd

.conj_vonmises2_mv <- function(s1 = NULL, priors = NULL,
                               plot = FALSE, support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  #* `set support to NULL to avoid default length of 10000`
  support <- NULL
  #* `Reorder columns if they are not in the numeric order`
  histColsBin <- as.numeric(sub("[a-zA-Z_.]+", "", colnames(s1)))
  bins_order <- sort(histColsBin, index.return = TRUE)$ix
  s1 <- s1[, bins_order]
  #* `Turn s1 matrix into a vector`
  X1 <- rep(histColsBin[bins_order], as.numeric(round(colSums(s1))))
  #* `make default prior if none provided`
  default_prior <- list(
    mu = 0, kappa = 0.5,
    boundary = c(-pi, pi),
    n = 1
  )
  if (is.null(priors)) {
    priors <- default_prior
  }
  #* `if any elements are missing from prior then use defaults`
  priors <- stats::setNames(lapply(names(default_prior), function(nm) {
    if (nm %in% names(priors)) {
      return(priors[[nm]])
    } else {
      return(default_prior[[nm]])
    }
  }), names(default_prior))
  #* `rescale data to [-pi, pi] according to boundary`
  X1 <- .boundary.to.radians(x = X1, boundary = priors$boundary)
  #* `rescale prior on mu to [-pi, pi] according to boundary`
  mu_radians <- .boundary.to.radians(x = priors$mu, boundary = priors$boundary)
  #* `Raise error if the boundary is wrong and data is not on [-pi, pi]`
  if (any(abs(X1) > pi)) {
    stop(paste0(
      "Values must be on [-pi, pi] after rescaling. ",
      "Does the boundary element in your prior include all your data?"
    ))
  }
  #* `Define dense Support`
  if (is.null(support)) {
    if (calculatingSupport) {
      return(priors$boundary) #* this would be [-pi, pi] if using radians, but plotting will be on
      #* the original scale so we can just return the boundary and use [-pi, pi] as support here
    }
    support_boundary <- seq(min(priors$boundary), max(priors$boundary), by = 0.0005)
    support <- seq(-pi, pi, length.out = length(support_boundary))
  }
  out <- list()
  #* ***** `Updating Kappa`
  n1 <- nrow(s1)
  obs_kappa <- .unbiased.kappa(X1, length(X1))
  kappa_prime <- ((obs_kappa * n1) + (priors$kappa * priors$n)) / (n1 + priors$n)
  #* ***** `Updating vMF for mu using kappa prime`
  #* `Get weighted mean of data and prior for half tangent adjustment`
  cm <- .circular.mean(c(X1, mu_radians), w = c(rep(nrow(s1) / length(X1), length(X1)), priors$n))
  unitCircleAdj <- ifelse(abs(cm) <= pi / 2, 0, pi)
  unitCircleAdj <- ifelse(cm > 0, 1, -1) * unitCircleAdj
  #* `Update prior parameters`
  a <- priors$kappa
  b <- mu_radians
  mu_prime_atan_scale <- atan(((a * sin(b)) + sum(sin(X1))) / ((a * cos(b)) + sum(cos(X1))))
  mu_prime <- unitCircleAdj + mu_prime_atan_scale
  #* `calculate density over support`
  dens1 <- brms::dvon_mises(support, mu_prime, kappa_prime)
  pdf1 <- dens1 / sum(dens1)
  #* `calculate highest density interval`
  #* note there is no qvon_mises function, so I am using bayestestR::hdi on
  #* posterior draws and rescaled posterior draws
  draws <- brms::rvon_mises(10000, mu_prime, kappa_prime)
  hdi_v1 <- as.numeric(bayestestR::hdi(draws, ci = cred.int.level))[2:3]
  draws2 <- draws
  draws2[draws2 < 0] <- draws2[draws2 < 0] + 2 * pi
  hdi_v2 <- as.numeric(bayestestR::hdi(draws2, ci = cred.int.level))[2:3]
  hdis <- list(hdi_v1, hdi_v2)
  hdi <- hdis[[which.min(c(diff(hdi_v1), diff(hdi_v2)))]]
  hdi[hdi > pi] <- hdi[hdi > pi] - (2 * pi) # if the second hdi was narrower then fix the part beyond pi
  #* `store highest density estimate`
  hde <- mu_prime
  #* `Rescale HDI, HDE, draws, and support from radians to boundary units`
  hdi_boundary <- .radians.to.boundary(hdi, target = priors$boundary)
  hde_boundary <- .radians.to.boundary(hde, target = priors$boundary)
  draws_boundary <- .radians.to.boundary(draws, target = priors$boundary)
  #* `save summary and parameters`
  out$summary <- data.frame(
    HDE_1 = hde_boundary,
    HDI_1_low = hdi_boundary[1],
    HDI_1_high = hdi_boundary[2]
  )
  out$posterior$mu <- hde_boundary # rescaled mu_prime
  out$posterior$kappa <- kappa_prime
  out$posterior$n <- priors$n + nrow(s1)
  out$posterior$boundary <- priors$boundary
  #* `Store Posterior Draws`
  out$posteriorDraws <- draws_boundary
  out$pdf <- pdf1
  #* `keep data for plotting`
  if (plot) {
    out$plot_df <- data.frame(
      "range" = support_boundary, "prob" = pdf1,
      "sample" = rep("Sample 1", length(support_boundary))
    )
  }
  return(out)
}




#' @description
#' Function calculating inverse of the ratio of first and zeroth order bessel functions.
#' This is a biased estimate of kappa, with small samples sizes the bias is strong.
#' @param x A vector of numerics from a von mises distribution
#' @examples
#' x <- brms::rvon_mises(100, 3, 2)
#' .bessel.inv(mean(cos(x - .circular.mean(x))))
#' @keywords internal
#' @noRd

.bessel.inv <- function(x) {
  ifelse(0 <= x & x < 0.53,
    2 * x + x^3 + (5 * x^5) / 6,
    ifelse(x < 0.85,
      -0.4 + 1.39 * x + 0.43 / (1 - x),
      1 / (x^3 - 4 * x^2 + 3 * x)
    )
  )
}

#' @description
#' helper function to estimate kappa based on data generated from a von mises distribution.
#' This is used from the CircStats package.
#' @param x sample of numeric draws from a von-mises distribution
#' @param n the number of samples used. If NULL, the default, then length(x) is used. If this is <16
#' then the bias adjustment from Best, D. and Fisher N. (1981) is used.
#' @examples
#' .unbiased.kappa(brms::rvon_mises(15, 3, 2))
#' @keywords internal
#' @noRd

.unbiased.kappa <- function(x, n = NULL) {
  if (is.null(n)) {
    n <- length(x)
  }
  mean.dir <- .circular.mean(x)
  kappa <- .bessel.inv(mean(cos(x - mean.dir)))
  if (n < 16) {
    kappa.biased <- kappa
    if (kappa.biased < 2) {
      kappa <- max(kappa.biased - 2 * (n * kappa.biased)^-1, 0)
    } else {
      kappa <- ((n - 1)^3 * kappa.biased) / (n^3 + n)
    }
  }
  kappa
}
