#' @description
#' Internal function for Bayesian T Tests of gaussian data represented by single value traits.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @examples
#' \donttest{
#' .conj_t_sv(
#'   s1 = rnorm(100, 50, 10), s2 = rnorm(100, 60, 12),
#'   priors = list(mu = c(0, 0), s2 = c(10, 10)),
#'   plot = FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'   cred.int.level = 0.89, hypothesis = "equal", support = NULL
#' )
#' }
#' @keywords internal
#' @noRd
.conj_t_sv <- function(s1 = NULL, priors = NULL,
                       plot = FALSE, support = NULL, cred.int.level = NULL,
                       calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  priors <- .convert_gaussian_priors(priors)
  if (is.null(priors)) {
    priors <- list(mu = 0, sd = 10)
  }
  #* `Calculate Sufficient Statistics`
  #* Using Mu and Precision Updating, see Equation 60 of compendium
  xbar <- mean(s1)
  n <- length(s1)
  #* `Get Prior Information`
  prec_prior <- 1 / (priors$sd[1] ^ 2)
  mu_prior <- priors$mu[1]
  #* `Update N(mu', prec') of Mu in N(Mu, Sd)`
  pseudo_known_precision <- 1 / ((var(s1) * n + priors$sd[1] ^ 2) / (n + 1))
  prec_prime <- prec_prior + pseudo_known_precision * n
  mu_prime <- ((mu_prior * prec_prior) + n * pseudo_known_precision * xbar) / prec_prime
  sd_prime <- sqrt(1 / prec_prime)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qnorm(c(0.0001, 0.9999), mu_prime, sd_prime)
    return(quantiles)
  }
  dens <- dnorm(support, mu_prime, sd_prime)
  pdf1 <- dens / sum(dens)
  hde1_mean <- mu_prime
  hdi1_mean <- qnorm(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
                     mu_prime, sd_prime)
  #* `Make Summary and Posterior for output`
  out$summary <- data.frame(HDE_1 = hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- mu_prime
  out$posterior$sd <- sd_prime
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- rnorm(10000, mu_prime, sd_prime)
  out$pdf <- pdf1
  #* `Save data for plotting`
  out$plot_list <- list(
    "range" = range(support),
    "ddist_fun" = "stats::dnorm",
    "priors" = list("mean" = priors$mu[1],
                    "sd" = priors$sd[1]),
    "parameters" = list("mean" = mu_prime,
                        "sd" = sd_prime)
  )
  return(out)
}





#' @description
#' Internal function for calculating \alpha and \beta of a distribution represented by multi value
#' traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number
#' representing the "bin".
#' @examples
#' \donttest{
#' mv_gauss <- mvSim(
#'   dists = list(
#'     rnorm = list(mean = 50, sd = 10)
#'   ),
#'   n_samples = 30
#' )
#' .conj_t_mv(
#'   s1 = mv_gauss[1:30, -1],
#'   priors = NULL,
#'   plot = TRUE,
#'   cred.int.level = 0.89
#' )
#' }
#' @keywords internal
#' @noRd

.conj_t_mv <- function(s1 = NULL, priors = NULL,
                       plot = FALSE, support = NULL, cred.int.level = NULL,
                       calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  priors <- .convert_gaussian_priors(priors)
  if (is.null(priors)) {
    priors <- list(mu = 0, sd = 10)
  }
  #* `Reorder columns if they are not in the numeric order`
  histColsBin <- as.numeric(sub("[a-zA-Z_.]+", "", colnames(s1)))
  bins_order <- sort(histColsBin, index.return = TRUE)$ix
  s1 <- s1[, bins_order]
  #* `Turn s1 matrix into a vector`
  X1 <- rep(histColsBin[bins_order], as.numeric(round(colSums(s1))))
  #* `Calculate Sufficient Statistics`
  #* See notes in SV version above
  xbar <- mean(X1)
  n <- nrow(s1)
  #* `Get Prior Information`
  prec_prior <- 1 / (priors$sd[1] ^ 2)
  mu_prior <- priors$mu[1]
  #* `Update N(mu', prec') of Mu in N(Mu, Sd)`
  pseudo_known_precision <- 1 / ((var(X1) * n + priors$sd[1] ^ 2) / (n + 1))
  prec_prime <- prec_prior + pseudo_known_precision * n
  mu_prime <- ((mu_prior * prec_prior) + n * pseudo_known_precision * xbar) / prec_prime
  sd_prime <- sqrt(1 / prec_prime)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qnorm(c(0.0001, 0.9999), mu_prime, sd_prime)
    return(quantiles)
  }
  dens <- dnorm(support, mu_prime, sd_prime)
  pdf1 <- dens / sum(dens)
  hde1_mean <- mu_prime
  hdi1_mean <- qnorm(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
                     mu_prime, sd_prime)
  #* `Make Summary and Posterior for output`
  out$summary <- data.frame(HDE_1 = hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- mu_prime
  out$posterior$sd <- sd_prime
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- rnorm(10000, mu_prime, sd_prime)
  out$pdf <- pdf1
  #* `Save data for plotting`
  out$plot_list <- list(
    "range" = support,
    "ddist_fun" = "stats::dnorm",
    "priors" = list("mean" = priors$mu[1],
                    "sd" = priors$sd[1]),
    "parameters" = list("mean" = mu_prime,
                        "sd" = sd_prime)
  )
  return(out)
}
