#' @description
#' Internal function for calculating the gamma distributed rate parameter of an exponential distribution
#' represented by single value traits.
#' @param s1 A vector of numerics drawn from a pareto distribution.
#' @examples
#' out <- .conj_exponential_sv(
#'   s1 = rexp(10, 3), cred.int.level = 0.95,
#'   plot = FALSE
#' )
#' lapply(out, head)
#' @keywords internal
#' @noRd
.conj_exponential_sv <- function(s1 = NULL, priors = NULL,
                                 support = NULL, cred.int.level = NULL,
                                 calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 0.5, b = 0.5)
  }
  #* `Update gamma prior with sufficient statistics`
  n <- length(s1)
  S <- sum(s1)
  a_prime <- priors$a[1] + n
  b_prime <- priors$b[1] + S
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qgamma(c(0.0001, 0.9999), a_prime, b_prime)
    return(quantiles)
  }
  #* `Make Posterior Draws`
  out$posteriorDraws <- rgamma(10000, a_prime, b_prime)
  #* `posterior`
  dens1 <- dgamma(support, a_prime, b_prime)
  pdf1 <- dens1 / sum(dens1)
  out$pdf <- pdf1
  hde1 <- .gammaHDE(shape = a_prime, scale = 1 / b_prime)
  hdi1 <- qgamma(
    c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    a_prime, b_prime
  )
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior <- list(
    "a" = a_prime,
    "b" = b_prime
  )
  out$prior <- priors
  #* `save s1 data for plotting`
  out$plot_list <- list(
    "range" = range(support),
    "ddist_fun" = "stats::dgamma",
    "priors" = list("shape" = priors$a[1],  "rate" = priors$b[1]),
    "parameters" = list("shape" = a_prime,
                        "rate" = b_prime)
  )
  return(out)
}
