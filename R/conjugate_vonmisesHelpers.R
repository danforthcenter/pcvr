#' @description
#' Internal function for calculating \mu and \kappa of a distribution represented by multi value
#' traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number
#' between 0.0001 and 0.9999 representing the "bin".
#' @examples
#' mv_gauss <- mvSim(
#'   dists = list(
#'     rnorm = list(mean = 50, sd = 10)
#'   ),
#'   n_samples = 30
#' )
#' .conj_vonmises_mv(
#'   s1 = mv_gauss[, -1], priors = list(mu = 30, kappa = 1, boundary = c(0, 180)),
#'   cred.int.level = 0.95,
#'   plot = TRUE
#' )
#' @keywords internal
#' @noRd

.conj_vonmises_mv <- function(s1 = NULL, priors = NULL,
                              plot = FALSE, support = NULL, cred.int.level = NULL,
                              calculatingSupport = FALSE) {
  #* `Turn off support for consistent rescaling between boundaries and to avoid default length of 10000`
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
    known_kappa = 1, n = 1
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
  #* `Get weighted mean of data and prior for half tangent adjustment`
  cm <- .circular.mean(c(X1, mu_radians), w = c(rep(nrow(s1) / length(X1), length(X1)), priors$n))
  unitCircleAdj <- ifelse(abs(cm) <= pi / 2, 0, pi)
  unitCircleAdj <- ifelse(cm > 0, 1, -1) * unitCircleAdj
  #* `Update prior parameters`
  a <- priors$kappa
  b <- mu_radians
  kappa_known <- priors$known_kappa
  kappa_prime <- kappa_known * .unbiased.kappa(X1)
  #* workaround for samples where kappa becomes negative if using the updating from the compendium
  #* where kappa prime is kappa_known x (A x sin B) + sum of sin data
  #*
  #* compendium and other sources I have read so far do not address this situation.
  #* seems like this would come up a lot, only difference I have seen is using [-pi, pi] vs
  #* the compendiums [0, 2pi], but I don't think that should make a difference.

  mu_prime_atan_scale <- atan(((a * sin(b)) + sum(sin(X1))) /
                                ((a * cos(b)) + sum(cos(X1))))
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
  out$posterior$known_kappa <- priors$known_kappa
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
#' Internal function for calculating \mu and \kappa of a distribution represented by single value
#' traits.
#' @param s1 A vector of numerics drawn from a beta distribution.
#' @examples
#' .conj_vonmises_sv(
#'   s1 = brms::rvon_mises(100, 2, 2), priors = list(mu = 0.5, kappa = 0.5),
#'   cred.int.level = 0.95,
#'   plot = FALSE
#' )
#' .conj_vonmises_sv(
#'   s1 = rnorm(20, 90, 20),
#'   priors = list(mu = 75, kappa = 0.5, boundary = c(0, 180), known_kappa = 2),
#'   cred.int.level = 0.95,
#'   plot = TRUE
#' )
#' @keywords internal
#' @noRd

.conj_vonmises_sv <- function(s1 = NULL, priors = NULL,
                              plot = FALSE, support = NULL, cred.int.level = NULL,
                              calculatingSupport = FALSE) {
  #* `to avoid default support length of 10000 which may not span boundary well`
  support <- NULL
  #* `make default prior if none provided`
  default_prior <- list(
    mu = 0, kappa = 0.5,
    boundary = c(-pi, pi),
    known_kappa = 1, n = 1
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
  #* `Get weighted mean of data and prior for half tangent adjustment`
  cm <- .circular.mean(c(s1, mu_radians), w = c(rep(1, length(s1)), priors$n))
  unitCircleAdj <- ifelse(abs(cm) <= pi / 2, 0, pi)
  unitCircleAdj <- ifelse(cm > 0, 1, -1) * unitCircleAdj
  #* `Update prior parameters`
  a <- priors$kappa
  b <- mu_radians
  kappa_known <- priors$known_kappa
  kappa_known <- priors$known_kappa
  kappa_prime <- kappa_known * .unbiased.kappa(s1)
  #* workaround for samples where kappa becomes negative if using the updating from the compendium
  #* where kappa prime is kappa_known x (A x sin B) + sum of sin data
  #* compendium and other sources I have read so far do not address this problem.
  #* seems like this would come up a lot, only difference I have seen is using [-pi, pi] vs
  #* the compendiums [0, 2pi], but I don't think that should make a difference.
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
  #* `Rescale HDI, HDE, draws, and support from radians to boundary units`
  hdi_boundary <- .radians.to.boundary(hdi, target = priors$boundary)
  hde_boundary <- .radians.to.boundary(hde, target = priors$boundary)
  draws_boundary <- .radians.to.boundary(draws, target = priors$boundary)
  support_boundary <- .radians.to.boundary(support, target = priors$boundary)
  #* `save summary and parameters`
  out$summary <- data.frame(
    HDE_1 = hde_boundary,
    HDI_1_low = hdi_boundary[1],
    HDI_1_high = hdi_boundary[2]
  )
  out$posterior$mu <- hde_boundary # rescaled mu_prime
  out$posterior$kappa <- kappa_prime
  out$posterior$known_kappa <- priors$known_kappa
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
  }
  return(out)
}

#' @description
#' Weighted Circular Mean function for use in von mises distribution conjugate function
#' @param x A vector of numerics drawn from a beta distribution.
#' @param w optional weights vector
#' @examples
#' if (FALSE) {
#'   .circular.mean(brms::rvon_mises(20, -3.1, 4))
#' }
#' @keywords internal
#' @noRd

.circular.mean <- function(x, w = NULL) {
  if (is.null(w)) {
    w <- rep(1, length(x))
  }
  sinr <- sum(sin(x) * w)
  cosr <- sum(cos(x) * w)
  circmean <- atan2(sinr, cosr)
  return(circmean)
}

#' @description
#' Function for rescaling any circular distribution to radians based on the edges.
#' @param x A vector of numerics from a circular process
#' @param boundary The edges of the circular process as a 2L numeric
#' @param target the target circular space, should not be changed from radians (-pi, pi).
#' @examples
#' if (FALSE) {
#'   x <- seq(-6, 6, length.out = 20)
#'   .boundary.to.radians(x, c(6.3, 6.3))
#'
#'   x <- seq(0, 100, length.out = 20)
#'   .boundary.to.radians(x, c(0, 100))
#' }
#' @keywords internal
#' @noRd

.boundary.to.radians <- function(x, boundary, target = c(-pi, pi)) {
  x1 <- (target[2] - target[1]) / (boundary[2] - boundary[1]) * (x - boundary[2]) + target[2]
  return(x1)
}

#' @description
#' Convenience function for easier reading. Same as .boundary.to.radians() with different defaults.
#' @param x A vector of numerics from a von mises distribution
#' @param boundary The edges of the circular process,
#' generally these should not be changed from radians.
#' @param target the target circular space, should be from priors$boundary.
#' @examples
#' if (FALSE) {
#'   x <- brms::rvon_mises(20, 2, 3)
#'   .radians.to.boundary(x, target = c(6.3, 6.3))
#'
#'   x <- brms::rvon_mises(20, 3.1, 2)
#'   .radians.to.boundary(x, target = c(0, 100))
#' }
#' @keywords internal
#' @noRd

.radians.to.boundary <- function(x, boundary = c(-pi, pi), target = c(-100, 100)) {
  x1 <- (target[2] - target[1]) / (boundary[2] - boundary[1]) * (x - boundary[2]) + target[2]
  return(x1)
}

#' @description
#' Calculate difference in draws from posterior of a von-mises distribution that has been
#' rescaled to exist on some circular space defined by the prior boundaries.
#' @param draws1 draws from a distribution
#' @param draws2 draws from another distribution
#' @param boundary a boundary vector describing the circular space's edges. Should be from priors.
#' @examples
#' if (FALSE) {
#'   draws1 <- brms::rvon_mises(10000, 3.1, 4)
#'   draws2 <- brms::rvon_mises(10000, -3, 2)
#'   x <- .conj_rope_circular_diff(draws1, draws2)
#' }
#' @keywords internal
#' @noRd

.conj_rope_circular_diff <- function(draws1, draws2, boundary = c(-pi, pi)) {
  draws1_radians <- .boundary.to.radians(draws1, boundary = boundary, target = c(-pi, pi))
  draws2_radians <- .boundary.to.radians(draws2, boundary = boundary, target = c(-pi, pi))
  span <- 2 * pi
  x <- draws1_radians + pi
  y <- draws2_radians + pi
  diff <- (x - y) %% span
  diff <- ifelse(diff <= (span / 2), diff, diff - span)
  diff <- .radians.to.boundary(diff, target = boundary) - mean(boundary)
  return(diff)
}
