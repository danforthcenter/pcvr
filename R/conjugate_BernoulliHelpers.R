#' @description
#' Internal function for calculating the beta distribution of the success rate of a bernoulli
#' distribution represented by single value traits.
#' @param s1 A vector of numerics drawn from a uniform distribution.
#' @examples
#' out <- .conj_bernoulli_sv(
#'   s1 = sample(c(TRUE, FALSE), 10, prob = c(0.3, 0.7), replace = TRUE),
#'   cred.int.level = 0.95,
#'   plot = TRUE
#' )
#' lapply(out, head)
#' @keywords internal
#' @noRd
.conj_bernoulli_sv <- function(s1 = NULL, priors = NULL,
                               plot = FALSE, support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 0.5, b = 0.5)
  }
  if (!is.logical(s1)) {
    stop("Bernoulli data must be supplied as a logical vector")
  }
  #* `Update beta prior with sufficient statistics`
  a1_prime <- priors$a[1] + sum(s1)
  b1_prime <- priors$b[1] + sum(!s1)
  #* `Define support if it is missing`
  if (is.null(support)) {
    if (calculatingSupport) {
      return(c(0.0001, 0.9999))
    }
    support <- seq(0.0001, 0.9999, 0.0001)
  }
  out <- list()
  #* `Make Posterior Draws`
  out$posteriorDraws <- rbeta(10000, a1_prime, b1_prime)
  #* `posterior`
  dens1 <- dbeta(support, a1_prime, b1_prime)
  pdf1 <- dens1 / sum(dens1)
  out$pdf <- pdf1
  #* `calculate highest density interval`
  hdi1 <- qbeta(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), a1_prime, b1_prime)
  #* `calculate highest density estimate``
  hde1 <- .betaHDE(a1_prime, b1_prime)
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior <- list("a" = a1_prime, "b" = b1_prime)
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
