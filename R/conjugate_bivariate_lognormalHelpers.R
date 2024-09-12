#' @description
#' Internal function for calculating the pareto distribution of the upper boundary of a uniform
#' distribution represented by single value traits.
#' @param s1 A vector of numerics drawn from a uniform distribution.
#' @examples
#' out <- .conj_bivariate_lognormal_sv(
#'   s1 = rlnorm(10, log(20), 1), cred.int.level = 0.95,
#'   plot = FALSE
#' )
#' lapply(out, head)
#'
#' @keywords internal
#' @noRd

.conj_bivariate_lognormal_sv <- function(s1 = NULL, priors = NULL,
                                         plot = FALSE, support = NULL, cred.int.level = NULL,
                                         calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  #* conjugate prior needs alpha, beta, mu, prec (or var or sd)
  #* precision is Gamma(alpha, beta)
  #* mu is T_[2*alpha](mu, precision)
  if (is.null(priors)) {
    priors <- list(mu = 0, sd = 10, a = 1, b = 1)
  }
  #* `Extract prior components`
  alpha <- priors$a[1]
  beta <- priors$b[1]
  mu <- priors$mu[1]
  prec <- 1 / (priors$sd^2)
  #* `Calculate sufficient statistics`
  n <- length(s1)
  x_bar <- mean(log(s1))
  ss <- sum((s1 - x_bar)^2)
  #* `Update priors with sufficient statistics`
  alpha_prime <- alpha + (n / 2)
  beta_prime <- 1 / ((1 / beta) + (ss / 2) + ((prec * n * ((x_bar - mu)^2)) / (2 * (prec + n))))
  mu_prime <- ((prec * mu) + (n * x_bar)) / (prec + n)
  prec_prime <- prec + n
  df_prime <- 2 * alpha_prime
  prec_prime_t <- alpha_prime * prec_prime * beta_prime
  sigma_prime <- sqrt(1 / prec_prime_t)
  #* `Define bivariate support if it is missing`
  if (is.null(support)) {
    quantiles_mu <- extraDistr::qlst(c(0.0001, 0.9999), df_prime, mu_prime, sigma_prime)
    quantiles_prec <- stats::qgamma(c(0.0001, 0.9999), shape = alpha_prime, scale = beta_prime)
    support_mu <- seq(quantiles_mu[1], quantiles_mu[2], length.out = 10000)
    support_prec <- seq(quantiles_prec[1], quantiles_prec[2], length.out = 10000)
    if (calculatingSupport) {
      return(list("Mu" = quantiles_mu, "Prec" = quantiles_prec))
    }
  } else {
    support_mu <- support$Mu
    support_prec <- support$Prec
  }
  #* `Make Posterior Draws`
  out$posteriorDraws <- .conj_biv_rough_sampling(
    10000, alpha_prime, beta_prime,
    mu_prime, sigma_prime, df_prime
  )
  #* `posterior`
  dens_mu <- extraDistr::dlst(support_mu, df_prime, mu_prime, sigma_prime)
  dens_prec <- stats::dgamma(support_prec, shape = alpha_prime, scale = beta_prime)

  pdf_mu <- dens_mu / sum(dens_mu)
  pdf_prec <- dens_prec / sum(dens_prec)
  out$pdf <- list("Mu" = pdf_mu, "Prec" = pdf_prec)

  hde_mu <- mu_prime
  hde_prec <- .gammaHDE(shape = alpha_prime, scale = beta_prime)
  hdi_mu <- -1 * rev(extraDistr::qlst(
    c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    df_prime, mu_prime, sigma_prime
  ))
  hdi_prec <- -1 * rev(stats::qgamma(
    c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    shape = alpha_prime, scale = beta_prime
  ))

  #* `Store summary`
  out$summary <- data.frame(
    HDE_1 = c(hde_mu, hde_prec),
    HDI_1_low = c(hdi_mu[1], hdi_prec[1]),
    HDI_1_high = c(hdi_mu[2], hdi_prec[2]),
    param = c("Mu", "Prec")
  )
  out$posterior <- list(
    "mu" = mu_prime, "sd" = sigma_prime,
    "a" = alpha_prime, "b" = beta_prime
  )
  #* `save s1 data for plotting`
  if (plot) {
    out$plot_df <- data.frame(
      "range" = c(support_mu, support_prec),
      "prob" = c(pdf_mu, pdf_prec),
      "param" = rep(c("Mu", "Prec"), each = length(support_mu)),
      "sample" = rep("Sample 1", 2 * length(support_mu))
    )
  }
  return(out)
}
