#' @description
#' Internal function for calculating the gamma distributed scale parameter of a pareto distribution
#' represented by multi value traits.
#' @param s1 A data.frame or matrix of multi value traits.
#' @examples
#' library(extraDistr)
#' s1 <- mvSim(
#'   dists = list(rpareto = list(a = 1, b = 1)),
#'   n_samples = 10,
#'   counts = 50,
#'   min_bin = 1,
#'   max_bin = 180,
#'   wide = TRUE
#' )[, -1]
#' out <- .conj_pareto_mv(s1, cred.int.level = 0.95)
#' lapply(out, head)
#'
#' @keywords internal
#' @noRd
.conj_pareto_mv <- function(s1 = NULL, priors = NULL,
                            plot = FALSE, support = NULL, cred.int.level = NULL,
                            calculatingSupport = FALSE) {
  out <- list()
  #* `N observations`
  n_obs <- nrow(s1)
  #* `Reorder columns if they are not in the numeric order`
  histColsBin <- as.numeric(sub("[a-zA-Z_.]+", "", colnames(s1)))
  bins_order <- sort(histColsBin, index.return = TRUE)$ix
  s1 <- s1[, bins_order]
  #* `Min non-zero bin`
  obs_min <- min(unlist(lapply(seq_len(n_obs), function(i) {
    col <- colnames(s1)[which(s1[i, ] > 0)][1]
    as.numeric(gsub("[a-zA-Z]_*", "", col))
  })), na.rm = TRUE)
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 0.5, b = 0.5, known_location = obs_min)
  }
  #* `Note, the product of the data is not an obvious quantity with MV trait data.`
  #* I talk about this some in my notes on 5/30/2024 but I think it might be a good chance
  #* to try the alternative mv trait method that I'd been considering.
  #* `Calculate Sufficient Statistics`
  #* This is abnormal because one of the sufficient statistics is the product of the data.
  #* That quantity does not translate well to the MV trait setting.
  #* `MLE Estimates of Pareto Parameters`
  #* Note this is being done per row of the MV data
  row_scales <- unlist(lapply(seq_len(n_obs), function(i) {
    d <- s1[i, ]
    #* `Turn s1 matrix into a vector`
    X1 <- rep(histColsBin[bins_order], as.numeric(round(colSums(d))))
    #* `Estimate parameters of pareto distribution of this row`
    scale_mle <- sum(d) / sum(log(X1))
    return(scale_mle)
  }))
  scale_estimate <- mean(row_scales)
  sv_draws <- extraDistr::rpareto(n_obs, scale_estimate, priors$known_location)
  #* `Update gamma prior with sufficient statistics`
  n <- length(sv_draws)
  m <- prod(sv_draws)
  a_prime <- priors$a + n
  b_prime <- 1 / (1 / priors$b + log(m) - n * log(priors$known_location))
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qgamma(c(0.0001, 0.9999), a_prime, b_prime)
    return(quantiles)
  }
  #* `Make Posterior Draws`
  out$posteriorDraws <- rgamma(10000, a_prime, b_prime)
  #* `posterior`
  dens1 <- dgamma(support, a_prime, b_prime)
  pdf1 <- dens1 / sum(dens1)
  out$pdf <- pdf1
  hde1 <- .gammaHDE(shape = a_prime, scale = 1 / b_prime)
  hdi1 <- qgamma(
    c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    a_prime, b_prime
  )
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior <- list(
    "a" = a_prime,
    "b" = b_prime,
    "known_location" = priors$known_location
  )
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
#' Internal function for calculating the gamma distributed scale parameter of a pareto distribution
#' represented by single value traits.
#' @param s1 A vector of numerics drawn from a pareto distribution.
#' @examples
#' out <- .conj_pareto_sv(
#'   s1 = runif(10, 1, 1), cred.int.level = 0.95,
#'   plot = FALSE
#' )
#' lapply(out, head)
#' @keywords internal
#' @noRd
.conj_pareto_sv <- function(s1 = NULL, priors = NULL,
                            plot = FALSE, support = NULL, cred.int.level = NULL,
                            calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 0.5, b = 0.5, known_location = floor(min(s1)))
  }
  #* `Update gamma prior with sufficient statistics`
  n <- length(s1)
  m <- prod(s1)
  a_prime <- priors$a + n
  b_prime <- 1 / (1 / priors$b + log(m) - n * log(priors$known_location))
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qgamma(c(0.0001, 0.9999), a_prime, b_prime)
    return(quantiles)
  }
  #* `Make Posterior Draws`
  out$posteriorDraws <- rgamma(10000, a_prime, b_prime)
  #* `posterior`
  dens1 <- dgamma(support, a_prime, b_prime)
  pdf1 <- dens1 / sum(dens1)
  out$pdf <- pdf1
  hde1 <- .gammaHDE(shape = a_prime, scale = 1 / b_prime)
  hdi1 <- qgamma(
    c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    a_prime, b_prime
  )
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior <- list(
    "a" = a_prime,
    "b" = b_prime,
    "known_location" = priors$known_location
  )
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
