#' @description
#' Internal function for calculating \alpha and \beta of the gamma distributed variance of lognormal
#' data given an estimate of the lognormal \mu obtained via the method of moments using multi value
#' traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number
#' representing the "bin".
#'
#' @examples
#' mv_ln <- mvSim(
#'   dists = list(
#'     rlnorm = list(meanlog = log(130), sdlog = log(1.2))
#'   ),
#'   n_samples = 30
#' )
#' .conj_lognormal2_mv(
#'   s1 = mv_ln[, -1],
#'   priors = NULL,
#'   plot = FALSE,
#'   cred.int.level = 0.89,
#' )
#' @keywords internal
#' @noRd

.conj_lognormal2_mv <- function(s1 = NULL, priors = NULL,
                                plot = FALSE, support = NULL, cred.int.level = NULL,
                                calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 1, b = 1) # prior on shape, scale of precision
  }
  #* `Reorder columns if they are not in the numeric order`
  histColsBin <- as.numeric(sub("[a-zA-Z_.]+", "", colnames(s1)))
  bins_order <- sort(histColsBin, index.return = TRUE)$ix
  s1 <- s1[, bins_order]
  #* `Loop over reps, get moments for each histogram`
  rep_distributions <- lapply(seq_len(nrow(s1)), function(i) {
    X1 <- rep(histColsBin[bins_order], as.numeric(s1[i, ]))
    #* `Get mean of x1`
    x_bar <- mean(X1)
    mu_s1 <- log(x_bar / (sqrt(var(X1) / x_bar^2) + 1))
    #* `Update Gamma Distribution of precision`
    #* sufficient stats: n, ss
    ss <- nrow(s1) * mean((log(X1) - mu_s1)^2) # mean * nrow instead of sum for MV traits
    n1 <- nrow(s1)
    a_prime <- priors$a[1] + (n1 / 2)
    b_prime <- priors$b[1] + (ss / 2)
    return(list("a_prime" = a_prime, "b_prime" = b_prime, "ln_mu" = mu_s1))
  })
  #* `Unlist parameters`
  a_prime <- mean(unlist(lapply(rep_distributions, function(i) {
    return(i$a_prime)
  })))
  b_prime <- mean(unlist(lapply(rep_distributions, function(i) {
    return(i$b_prime)
  })))
  ln_mu_prime <- mean(unlist(lapply(rep_distributions, function(i) {
    return(i$ln_mu)
  })))
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qgamma(c(0.0001, 0.9999), shape = a_prime, scale = b_prime)
    return(quantiles)
  }
  #* `posterior`
  dens1 <- dgamma(support, shape = a_prime, scale = b_prime)
  pdf1 <- dens1 / sum(dens1)
  hde1 <- .gammaHDE(shape = a_prime, scale = b_prime)
  hdi1 <- qgamma(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    shape = a_prime, scale = b_prime
  )
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a_prime
  out$posterior$b <- b_prime
  out$posterior$lognormal_mu <- ln_mu_prime # returning this as a number, not a distribution
  #* `Make Posterior Draws`
  out$posteriorDraws <- rgamma(10000, shape = a_prime, scale = b_prime)
  out$pdf <- pdf1
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
#' Internal function for calculating \alpha and \beta of the gamma distributed precision of lognormal
#' data given an estimate of the lognormal \mu obtained via the method of moments using single value
#' traits.
#'
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @examples
#' .conj_lognormal2_sv(
#'   s1 = rlnorm(100, log(130), log(1.3)),
#'   priors = NULL,
#'   plot = FALSE,
#'   cred.int.level = 0.89
#' )
#' @keywords internal
#' @noRd

.conj_lognormal2_sv <- function(s1 = NULL, priors = NULL,
                                plot = FALSE, support = NULL, cred.int.level = NULL,
                                calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 1, b = 1) # prior on shape, scale of precision
  }
  #* `Get mean of s1`
  x_bar <- mean(s1)
  mu_s1 <- log(x_bar / (sqrt(var(s1) / x_bar^2) + 1))
  #* `Update Gamma Distribution of precision`
  #* sufficient stats: n, ss
  ss <- sum((log(s1) - mu_s1)^2)
  n1 <- length(s1)
  a_prime <- priors$a[1] + (n1 / 2)
  b_prime <- priors$b[1] + (ss / 2)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qgamma(c(0.0001, 0.9999), shape = a_prime, scale = b_prime)
    return(quantiles)
  }
  #* `posterior`
  dens1 <- dgamma(support, shape = a_prime, scale = b_prime)
  pdf1 <- dens1 / sum(dens1)
  hde1 <- .gammaHDE(shape = a_prime, scale = b_prime)
  hdi1 <- qgamma(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    shape = a_prime, scale = b_prime
  )
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a_prime
  out$posterior$b <- b_prime
  out$posterior$lognormal_mu <- mu_s1 # returning this as a number, not a distribution
  #* `Make Posterior Draws`
  out$posteriorDraws <- rgamma(10000, shape = a_prime, scale = b_prime)
  out$pdf <- pdf1
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
