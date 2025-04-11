#' @description
#' Internal function for Bayesian comparisons of gaussian data represented by single value traits.
#' This version uses the entire posterior distribution instead of the sampling distribution of the mean.
#' In frequentist terms this is analogous to a Z test as opposed to a T test. Generally the T test is
#' desired, but this is provided for completeness.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @examples
#' .conj_gaussian_sv(
#'   s1 = rnorm(100, 50, 10),
#'   priors = list(mu = c(0, 0), sd = c(10, 10)),
#'   plot = FALSE, support = NULL
#' )
#' @keywords internal
#' @noRd

.conj_gaussian_sv <- function(s1 = NULL, priors = NULL,
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
#' Internal function for Bayesian comparisons of gaussian data represented by multi value traits.
#' This version uses the entire posterior distribution instead of the sampling distribution of the mean.
#' In frequentist terms this is analogous to a Z test as opposed to a T test. Generally the T test is
#' desired, but this is provided for completeness.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @keywords internal
#' @noRd

.conj_gaussian_mv <- function(s1 = NULL, priors = NULL,
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
#' Helper function to catch gaussian (and T/Lognormal) priors specified with other
#' parameterizations besides standard deviation.
#' @param priors Priors specified with mu and sd, var, or prec. Var or Prec are converted to sd.
#' @keywords internal
#' @noRd

.convert_gaussian_priors <- function(priors) {
  if (any(grepl("prec", names(priors)))) {
    precs <- priors[[which(grepl("prec", names(priors)))]]
    priors[[which(grepl("prec", names(priors)))]] <- sqrt(1 / precs)
    names(priors)[which(grepl("prec", names(priors)))] <- "sd"
  } else if (any(grepl("s2|var", names(priors)))) {
    vars <- priors[[which(grepl("s2|var", names(priors)))]]
    priors[[which(grepl("s2|var", names(priors)))]] <- sqrt(vars)
    names(priors)[which(grepl("s2|var", names(priors)))] <- "sd"
  }
  return(priors)
}
