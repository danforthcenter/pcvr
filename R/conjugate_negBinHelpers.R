#' @description
#' Negative binomial via conjugate beta prior GIVEN A KNOWN R
#' conjugacy:
#' if counts ~ nbinom(n, r)
#' where n is the number of successful trials
#'   and r is the prob of success per trial which is fixed and known
#'
#' if r is known then p ~ beta(A, B)
#' note that the compendium of conjugate priors seems to have a typo for this relationship,
#' but the wikipedia conjugate prior article is correct
#' \code{A' = A + r*n}
#' \code{B' = B + sum(X)}
#'
#' Using MoM:
#'
#' \bar{x} = k(1-p)/p
#' s^2 = \bar{x}/p
#' r = \bar{x}^2 / (s^2 - \bar{x})
#' p = \bar{x}/s^2
#'
#' @param s1 A vector of numerics drawn from a negative binomial distribution.
#' @examples
#' .conj_negbin_sv(
#'   s1 = rnbinom(10, 10, 0.5),
#'   priors = NULL,
#'   plot = FALSE,
#'   cred.int.level = 0.89
#' )
#' @keywords internal
#' @noRd

.conj_negbin_sv <- function(s1 = NULL, priors = NULL,
                            support = NULL, cred.int.level = NULL,
                            calculatingSupport = FALSE) {
  #* `Check samples`
  if (any(abs(s1 - round(s1)) > .Machine$double.eps^0.5) || any(s1 < 0)) {
    stop("Only positive whole numbers can be used in the Negative Binomial distribution")
  }
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(r = 10, a = 0.5, b = 0.5) # beta prior on P
    warning(paste0(
      "True value of r for negative binomial distribution has defaulted to 10,",
      " you should add a prior including r parameter."
    ))
  }

  out <- list()

  #* `Use conjugate beta prior on probability`
  #* Note that this is very sensitive to the R value being appropriate
  a1_prime <- priors$a[1] + priors$r[1] * length(s1)
  b1_prime <- priors$b[1] + sum(s1)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    return(c(0.0001, 0.9999))
  }
  #* `calculate density over support``
  dens1 <- dbeta(support, a1_prime, b1_prime)
  pdf1 <- dens1 / sum(dens1)

  #* `calculate highest density interval`
  hdi1 <- qbeta(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), a1_prime, b1_prime)

  #* `calculate highest density estimate``
  hde1 <- .betaHDE(a1_prime, b1_prime)

  #* `save summary and parameters`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$r <- priors$r[1]
  out$posterior$a <- a1_prime
  out$posterior$b <- b1_prime
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- rbeta(10000, a1_prime, b1_prime)
  out$pdf <- pdf1
  #* `keep data for plotting`
  out$plot_list <- list(
    "range" = range(support),
    "ddist_fun" = "stats::dbeta",
    "priors" = list("shape1" = priors$a[1],  "shape2" = priors$b[1]),
    "parameters" = list("shape1" = a1_prime,
                        "shape2" = b1_prime),
    "given" = list("size" = priors$r[1])
  )
  return(out)
}
