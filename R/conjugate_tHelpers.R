#' @description
#' Internal function for Bayesian T Tests of gaussian data represented by single value traits.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @examples
#' \donttest{
#' .conj_t_sv(
#'   s1 = rnorm(100, 50, 10), s2 = rnorm(100, 60, 12),
#'   priors = list(mu = c(0, 0), n = c(1, 1), s2 = c(20, 20)),
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
  #* `Update N(mu', var') of Mu in N(Mu, Var)`
  mu_prime <- (var_null / ((var_data / n) + var_null)) * xbar +
    (var_data / ((var_data / n) + var_null)) * mu_null
  var_prime <- ((1 / var_null) + (n / var_data)) ^ -1
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

  out$summary <- data.frame(HDE_1 = hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- m1_n
  out$posterior$sd <- sd_prime # return variance
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- rnorm(10000, mu_prime, sd_prime)
  out$pdf <- pdf1
  #* `Save data for plotting`
  out$plot_list <- list(
    "range" = support,
    "ddist_fun" = "stats::rnorm",
    "priors" = list("mu" = priors$mu[1],
                    "sd" = priors$sd[1]),
    "parameters" = list("mu" = mu_prime,
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
  se1 <- sqrt(s2_1_n / n1_n) # standard error of the mean
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qlst(c(0.0001, 0.9999), v1_n, m1_n, se1)
    return(quantiles)
  }
  dens <- extraDistr::dlst(support, v1_n, m1_n, se1)
  pdf1 <- dens / sum(dens)
  hde1_mean <- m1_n
  hdi1_mean <- m1_n + qt(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), v1_n) * se1

  out$summary <- data.frame(HDE_1 = hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- m1_n
  out$posterior$n <- n1_n
  out$posterior$s2 <- s2_1_n
  out$prior <- priors
  #* `Make Posterior Draws`
  out$posteriorDraws <- extraDistr::rlst(10000, v1_n, m1_n, se1)
  out$pdf <- pdf1
  #* `Save data for plotting`
  out$plot_list <- list(
    "range" = support,
    "ddist_fun" = "extraDistr::dlst",
    "priors" = list("df" = max(c(2, priors$n[1])),
                    "mu" = priors$mu[1],
                    "sigma" = sqrt(priors$s2[1] / priors$n[1])),
    "parameters" = list("df" = v1_n,
                        "mu" = m1_n,
                        "sigma" = se1)
  )

  return(out)
}
