#' @description
#' Internal function for Bayesian comparisosns of gaussian data represented by single value traits.
#' This version uses the entire posterior distribution instead of the sampling distribution of the mean.
#' In frequentist terms this is analogous to a Z test as opposed to a T test. Generally the T test is
#' desired, but this is provided for completeness.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @examples
#' .conj_gaussian_sv(
#'   s1 = rnorm(100, 50, 10),
#'   priors = list(mu = c(0, 0), n = c(1, 1), s2 = c(20, 20)),
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
    priors <- list(mu = 0, n = 1, s2 = 100)
  }
  #* `Get Mean, Variance, SE, and DF from s1`

  n1 <- length(s1) # n samples
  m1 <- mean(s1) # xbar
  s2_1 <- var(s1) # variance
  v1 <- priors$n[1] - 1 # prior DF
  n1_n <- priors$n[1] + n1 # total N including prior
  m1_n <- (n1 * m1 + priors$n[1] * priors$mu[1]) / n1_n # weighted mean of prior and data
  v1_n <- v1 + n1 # degrees of freedom including data
  s2_1_n <- ((n1 - 1) * s2_1 + v1 * priors$s2[1] + priors$n[1] * n1 * (priors$mu[1] - m1)^2 / n1_n) /
    v1_n # pooled variance
  sigma_1 <- sqrt(s2_1_n)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qlst(c(0.0001, 0.9999), v1_n, m1_n, sigma_1)
    return(quantiles)
  }

  dens1 <- extraDistr::dlst(support, v1_n, m1_n, sigma_1)
  pdf1 <- dens1 / sum(dens1)
  hde1_mean <- m1_n
  hdi1_mean <- m1_n + qt(c(
    (1 - cred.int.level) / 2,
    (1 - ((1 - cred.int.level) / 2))
  ), v1_n) * sigma_1

  out$summary <- data.frame(HDE_1 = hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- m1_n
  out$posterior$n <- n1_n
  out$posterior$s2 <- s2_1_n # return variance
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- extraDistr::rlst(10000, v1_n, m1_n, sigma_1)
  out$pdf <- pdf1
  #* `Save data for plotting`
  if (plot) {
    out$plot_list <- list(
      "range" = support,
      "ddist_fun" = "stats::dnorm",
      "parameters" = list("mean" = m1_n,
                          "sd" = sigma_1)
    )
  }
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
    priors <- list(mu = 0, n = 1, s2 = 100)
  }
  #* `Reorder columns if they are not in the numeric order`
  histColsBin <- as.numeric(sub("[a-zA-Z_.]+", "", colnames(s1)))
  bins_order <- sort(histColsBin, index.return = TRUE)$ix
  s1 <- s1[, bins_order]


  #* `Turn s1 matrix into a vector`
  X1 <- rep(histColsBin[bins_order], as.numeric(round(colSums(s1))))

  #* `Get Mean, Variance, SE, and DF from s2`
  n1 <- nrow(s1) # n samples
  m1 <- mean(X1) # xbar
  s2_1 <- var(X1) # var

  v1 <- priors$n[1] - 1 # prior DF
  n1_n <- priors$n[1] + n1 # total N including prior
  m1_n <- (n1 * m1 + priors$n[1] * priors$mu[1]) / n1_n # weighted mean of prior and data
  v1_n <- v1 + n1 # degrees of freedom including data
  s2_1_n <- ((n1 - 1) * s2_1 + v1 * priors$s2[1] + priors$n[1] * n1 * (priors$mu[1] - m1)^2 / n1_n) /
    v1_n # pooled variance
  sigma_1 <- sqrt(s2_1_n) # standard deviation
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qlst(c(0.0001, 0.9999), v1_n, m1_n, sigma_1)
    return(quantiles)
  }
  dens <- extraDistr::dlst(support, v1_n, m1_n, sigma_1)
  pdf1 <- dens / sum(dens)
  hde1_mean <- m1_n
  hdi1_mean <- m1_n + qt(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), v1_n) * sigma_1

  out$summary <- data.frame(HDE_1 = hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- m1_n
  out$posterior$n <- n1_n
  out$posterior$s2 <- s2_1_n
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- extraDistr::rlst(10000, v1_n, m1_n, sigma_1)
  out$pdf <- pdf1
  #* `Save data for plotting`
  if (plot) {
    out$plot_list <- list(
      "range" = support,
      "ddist_fun" = "stats::dnorm",
      "parameters" = list("mean" = m1_n,
                          "sd" = sigma_1)
    )
  }
  return(out)
}
