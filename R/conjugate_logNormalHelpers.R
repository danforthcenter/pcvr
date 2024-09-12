#' @description
#' Internal function for calculating \mu and \sigma of the normally distributed mean of lognormal data
#' given an estimate of the lognormal \sigma obtained via the method of moments using multi value
#' traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number
#' representing the "bin".
#'
#' @examples
#'
#' mv_ln <- mvSim(
#'   dists = list(
#'     rlnorm = list(meanlog = log(130), sdlog = log(1.2))
#'   ),
#'   n_samples = 30
#' )
#' .conj_lognormal_mv(
#'   s1 = mv_ln[1:30, -1],
#'   priors = list(mu_log = c(log(10)), n = c(1), sigma_log = c(log(3))),
#'   plot = FALSE, cred.int.level = 0.9
#' )
#'
#' @importFrom stats qnorm
#' @keywords internal
#' @noRd

.conj_lognormal_mv <- function(s1 = NULL, priors = NULL,
                               plot = FALSE, support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(mu = 0, sd = 5)
  }
  #* `Reorder columns if they are not in the numeric order`
  histColsBin <- as.numeric(sub("[a-zA-Z_.]+", "", colnames(s1)))
  bins_order <- sort(histColsBin, index.return = TRUE)$ix
  s1 <- s1[, bins_order]

  #* `Loop over reps, get moments for each histogram`

  rep_distributions <- lapply(seq_len(nrow(s1)), function(i) {
    X1 <- rep(histColsBin[bins_order], as.numeric(s1[i, ]))
    #* `Get mean of X1`
    x_bar <- mean(X1)
    mu_x1 <- log(x_bar / (sqrt(var(X1) / x_bar^2) + 1))
    #* `Get sigma of X1`
    sigma_x1 <- sqrt(log((var(X1)) / (x_bar^2) + 1))
    #* `Update Normal Distribution of Mu`
    #* sufficient stats: n, mean of log data | precision
    n <- length(X1)
    m <- priors$mu[1]
    p <- 1 / (priors$sd[1]^2) # precision
    mu_prime <- ((m * p) + (n * p * mu_x1)) / (p + (n * p))
    precision_prime <- (p + (n * p))
    var_prime <- 1 / precision_prime
    sd_prime <- sqrt(var_prime)
    return(list("mu" = mu_prime, "sd" = sd_prime, "ln_sd" = sigma_x1))
  })
  #* `Unlist parameters`
  mu_ls_prime <- mean(unlist(lapply(rep_distributions, function(i) {
    i$mu
  })))
  sigma_ls_prime <- mean(unlist(lapply(rep_distributions, function(i) {
    i$sd
  })))
  ln_sigma_prime <- mean(unlist(lapply(rep_distributions, function(i) {
    i$ln_sd
  })))

  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- stats::qnorm(c(0.0001, 0.9999), mu_ls_prime, sigma_ls_prime)
    return(quantiles)
  }
  #* `posterior`
  dens1 <- stats::dnorm(support, mu_ls_prime, sigma_ls_prime)
  pdf1 <- dens1 / sum(dens1)
  hde1 <- mu_ls_prime
  hdi1 <- stats::qnorm(
    c(
      (1 - cred.int.level) / 2,
      (1 - ((1 - cred.int.level) / 2))
    ),
    mu_ls_prime, sigma_ls_prime
  )
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$mu <- mu_ls_prime
  out$posterior$sd <- sigma_ls_prime
  out$posterior$lognormal_sigma <- ln_sigma_prime
  #* `Make Posterior Draws`
  out$posteriorDraws <- stats::rnorm(10000, mu_ls_prime, sigma_ls_prime)
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
#' Internal function for calculating \mu and \sigma of the normally distributed mean of lognormal data
#' given an estimate of the lognormal \sigma obtained via the method of moments using single value
#' traits.
#'
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @examples
#' .conj_lognormal_sv(
#'   s1 = rlnorm(100, log(130), log(1.3)),
#'   priors = list(mu = 5, sd = 5),
#'   plot = FALSE,
#'   cred.int.level = 0.89
#' )
#' @keywords internal
#' @noRd

.conj_lognormal_sv <- function(s1 = NULL, priors = NULL,
                               plot = FALSE, support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(mu = 0, sd = 5)
  }
  #* `Get mean of s1`
  x_bar <- mean(s1)
  mu_s1 <- log(x_bar / (sqrt(var(s1) / x_bar^2) + 1))
  #* `Get sigma of s1`
  sigma_s1 <- sqrt(log((var(s1)) / (x_bar^2) + 1))
  #* `Update Normal Distribution of Mu`
  #* sufficient stats: n, mean of log data | precision
  n <- length(s1)
  m <- priors$mu[1]
  p <- 1 / (priors$sd[1]^2) # precision
  mu_prime <- ((m * p) + (n * p * mu_s1)) / (p + (n * p))
  precision_prime <- (p + (n * p))
  var_prime <- 1 / precision_prime
  sd_prime <- sqrt(var_prime)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- stats::qnorm(c(0.0001, 0.9999), mu_prime, sd_prime)
    return(quantiles)
  }
  #* `posterior`
  dens1 <- dnorm(support, mu_prime, sd_prime)
  pdf1 <- dens1 / sum(dens1)
  hde1 <- mu_prime
  hdi1 <- qnorm(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), mu_prime, sd_prime)
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$mu <- mu_prime
  out$posterior$sd <- sd_prime
  out$posterior$lognormal_sigma <- sigma_s1 # returning this as a number, not a distribution
  #* `Make Posterior Draws`
  out$posteriorDraws <- rnorm(10000, mu_prime, sd_prime)
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
