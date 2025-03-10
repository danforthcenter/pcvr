#' @description
#' Internal function for Bayesian comparisosns of gaussian data represented by single value traits.
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
  if (is.null(priors)) {
    priors <- list(mu = 0, sd = 10)
  }
  #* `Calculate Sufficient Statistics`
  #* Using Mu and Variance Updating, see Equation 12 of
  #* https://people.eecs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf
  #* alternatively equation 60 of compendium of conjugate priors shows the mu/precision
  #* formula, which is simpler
  var_data <- var(s1)
  xbar <- mean(s1)
  n <- length(s1)
  #* `Get Prior Information`
  var_null <- priors$sd[1] ^ 2
  mu_null <- priors$mu[1]
  #* `Update Normal Distribution`
  mu_prime <- (var_null / ((var_data / n) + var_null)) * xbar +
    (var_data / ((var_data / n) + var_null)) * mu_null
  #* Currently the difference between the T and Z methods
  #* is whether this is 1 over var_data or N over var_data,
  #* with N making for variance of the mean, 1 making for variance of the data
  var_prime <- ((1 / var_null) + (1 / var_data)) ^ -1
  sd_prime <- sqrt(var_prime)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qnorm(c(0.0001, 0.9999), mu_prime, sd_prime)
    return(quantiles)
  }
  dens <- qnorm(support, mu_prime, sd_prime)
  pdf1 <- dens / sum(dens)
  hde1_mean <- mu_prime
  hdi1_mean <- qnorm(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
                     mu_prime, sd_prime)
  #* `Make Summary and Posterior for output`
  out$summary <- data.frame(HDE_1 = hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- m1_n
  out$posterior$sd <- sd_prime
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- rnorm(10000, mu_prime, sd_prime)
  out$pdf <- pdf1
  #* `Save data for plotting`
  out$plot_list <- list(
    "range" = support,
    "ddist_fun" = "stats::dnorm",
    "priors" = list("mu" = priors$mu[1],
                    "sd" = priors$sd[1]),
    "parameters" = list("mu" = mu_prime,
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
  n <- nrow(s1)
  xbar <- mean(X1)
  var_data <- var(X1)
  #* `Get Prior Information`
  var_null <- priors$sd[1] ^ 2
  mu_null <- priors$mu[1]
  #* `Update Normal Distribution`
  mu_prime <- (var_null / ((var_data / n) + var_null)) * xbar +
    (var_data / ((var_data / n) + var_null)) * mu_null
  #* See comment above in SV version about T vs Z
  var_prime <- ((1 / var_null) + (1 / var_data)) ^ -1
  sd_prime <- sqrt(var_prime)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qnorm(c(0.0001, 0.9999), mu_prime, sd_prime)
    return(quantiles)
  }
  dens <- qnorm(support, mu_prime, sd_prime)
  pdf1 <- dens / sum(dens)
  hde1_mean <- mu_prime
  hdi1_mean <- qnorm(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
                     mu_prime, sd_prime)
  #* `Make Summary and Posterior for output`
  out$summary <- data.frame(HDE_1 = hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- m1_n
  out$posterior$sd <- sd_prime
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- rnorm(10000, mu_prime, sd_prime)
  out$pdf <- pdf1
  #* `Save data for plotting`
  out$plot_list <- list(
    "range" = support,
    "ddist_fun" = "stats::dnorm",
    "priors" = list("mu" = priors$mu[1],
                    "sd" = priors$sd[1]),
    "parameters" = list("mu" = mu_prime,
                        "sd" = sd_prime)
  )
  return(out)
}
