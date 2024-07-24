#' @description
#' Internal function for calculating the gamma distribution of the rate parameter in gamma distributed
#' data represented by single value traits.
#' @param s1 A vector of numerics drawn from a uniform distribution.
#' @examples
#' out <- .conj_gamma_sv(
#'   s1 = rgamma(10, 1, 2), cred.int.level = 0.95,
#'   plot = FALSE
#' )
#' lapply(out, head)
#' @keywords internal
#' @noRd
.conj_gamma_sv <- function(s1 = NULL, priors = NULL,
                           plot = FALSE, support = NULL, cred.int.level = NULL,
                           calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(shape = 0.5, scale = 0.5, known_shape = 1)
  }
  #* `Update gamma prior with sufficient statistics`
  n <- length(s1)
  S <- sum(s1)
  shape_prime <- (priors$known_shape * n) + priors$shape
  scale_prime <- priors$scale / (1 + (priors$scale * S))
  #* `Define support if it is missing`
  if (is.null(support)) {
    quantiles <- qgamma(c(0.0001, 0.9999), shape = shape_prime, scale = scale_prime)
    if (calculatingSupport) {
      return(quantiles)
    }
    support <- seq(quantiles[1], quantiles[2], length.out = 10000)
  }
  #* `Make Posterior Draws`
  out$posteriorDraws <- rgamma(10000, shape = shape_prime, scale = scale_prime)
  #* `posterior`
  dens1 <- dgamma(support, shape = shape_prime, scale = scale_prime)
  pdf1 <- dens1 / sum(dens1)
  out$pdf <- pdf1
  hde1 <- .gammaHDE(shape_prime, scale_prime)
  hdi1 <- qgamma(
    c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    shape = shape_prime, scale = scale_prime
  )
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior <- list("shape" = shape_prime, "scale" = scale_prime,
                        "known_shape" = priors$known_shape)
  #* `save s1 data for plotting`
  if (plot) {
    out$plot_df <- data.frame(
      "range" = support,
      "prob" = pdf1,
      "sample" = rep("Sample 1", length(support))
    )
  }
  return(out)
}

#' @description
#' Internal function for calculating the HDE of a gamma distribution
#' @param shape shape parameter
#' @param scale scale parameter
#' @examples
#' .gammaHDE(1, 2)
#' .gammaHDE(0, 1)
#' .gammaHDE(10, 10)
#' @keywords internal
#' @noRd

.gammaHDE <- function(shape, scale) {
  if (shape >= 1) {
    hde <- (shape - 1) * scale
  } else {
    hde <- 0
  }
  return(hde)
}
