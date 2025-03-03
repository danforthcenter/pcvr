#' @description
#' Internal function for calculating posterior distribution of \Lambda parameter for poisson data.
#' The conjugate prior on Lambda is a Gamma(A,B)
#' A=B=0.5 is a reasonable weak default prior
#'
#' So if leaf count is Poisson distributed:
#' count ~ Pois(\labmda)
#' \labmda ~ gamma(A, B)
#' A = A_[prior] + sum(x)
#' B = B_[prior] / (1+n)
#'
#' via MoM \hat(\labmda) =  = 1/n +sum^1_n(x)
#' @param s1 A vector of numerics drawn from a beta distribution.
#' @examples
#'
#' .conj_poisson_sv(
#'   s1 = rpois(20, 10), priors = list(a = c(0.5, 0.5), b = c(0.5, 0.5)),
#'   plot = FALSE
#' )
#'
#' @keywords internal
#' @noRd

.conj_poisson_sv <- function(s1 = NULL, priors = NULL,
                             plot = FALSE, support = NULL, cred.int.level = NULL,
                             calculatingSupport = FALSE) {
  #* `Check samples`
  if (any(abs(s1 - round(s1)) > .Machine$double.eps^0.5) || any(s1 < 0)) {
    stop("Only positive integers can be used in the Poisson distribution")
  }
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 0.5, b = 0.5) # gamma prior on lambda
  }

  out <- list()

  #* `Use conjugate gamma prior on lambda`
  a1_prime <- priors$a[1] + sum(s1)
  b1_prime <- priors$b[1] + length(s1)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qgamma(c(0.0001, 0.9999), a1_prime, b1_prime)
    return(quantiles)
  }
  #* `calculate density over support``
  dens1 <- dgamma(support, a1_prime, b1_prime)
  pdf1 <- dens1 / sum(dens1)

  #* `calculate highest density interval`
  hdi1 <- qgamma(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), a1_prime, b1_prime)

  #* `calculate highest density estimate``
  hde1 <- .gammaHDE(shape = a1_prime, scale = 1 / b1_prime)

  #* `save summary and parameters`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a1_prime
  out$posterior$b <- b1_prime
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- rgamma(10000, a1_prime, b1_prime)
  out$pdf <- pdf1
  #* `keep data for plotting`
  if (plot) {
    out$plot_list <- list(
      "range" = support,
      "ddist_fun" = "stats::dgamma",
      "parameters" = list("shape" = a1_prime,
                          "rate" = b1_prime)
    )
  }
  return(out)
}
