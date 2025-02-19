#' @description
#' Internal function for calculating the pareto distribution of the upper boundary of a uniform
#' distribution represented by multi value traits.
#' @param s1 A data.frame or matrix of multi value traits.
#' @examples
#' s1 <- mvSim(
#'   dists = list(runif = list(min = 0, max = 150)),
#'   n_samples = 10,
#'   counts = 1000,
#'   min_bin = 1,
#'   max_bin = 180,
#'   wide = TRUE
#' )
#' out <- .conj_uniform_mv(s1, cred.int.level = 0.95)
#' lapply(out, head)
#' @importFrom utils tail
#'
#' @keywords internal
#' @noRd
.conj_uniform_mv <- function(s1 = NULL, priors = NULL,
                             plot = FALSE, support = NULL, cred.int.level = NULL,
                             calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(scale = 0.5, location = 0.5)
  }
  #* `Reorder columns if they are not in the numeric order`
  histColsBin <- as.numeric(sub("[a-zA-Z_.]+", "", colnames(s1)))
  bins_order <- sort(histColsBin, index.return = TRUE)$ix
  s1 <- s1[, bins_order]
  #* `Calculate Sufficient Statistics`
  #* `N observations`
  n_obs <- nrow(s1)
  #* `Max non-zero bin`
  max_obs <- max(unlist(lapply(seq_len(n_obs), function(i) {
    col <- utils::tail(colnames(s1)[which(s1[i, ] > 0)], 1)
    return(as.numeric(gsub("[a-zA-Z]_*", "", col)))
  })), na.rm = TRUE)
  #* `Update pareto prior with sufficient statistics`
  scale_prime <- priors$scale + n_obs
  location_prime <- max(c(max_obs, priors$location), na.rm = TRUE)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qpareto(c(0.0001, 0.9999), scale_prime, location_prime)
    return(quantiles)
  }
  #* `Make Posterior Draws`
  out$posteriorDraws <- extraDistr::rpareto(10000, scale_prime, location_prime)
  #* `posterior`
  dens1 <- extraDistr::dpareto(support, scale_prime, location_prime)
  pdf1 <- dens1 / sum(dens1)
  out$pdf <- pdf1
  # hde of location is calculated off of the posterior draws
  hde1 <- mean(as.numeric(hdi(out$posteriorDraws, ci = 0.01)[c("CI_low", "CI_high")]))
  # mode is location
  # mean is defined if scale > 1: scale_prime * location_prime / scale_prime - 1
  # median is location x root scale of 2
  hdi1 <- extraDistr::qpareto(
    c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    scale_prime, location_prime
  )
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior <- list("scale" = scale_prime, "location" = location_prime)
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
#' Internal function for calculating the pareto distribution of the upper boundary of a uniform
#' distribution represented by single value traits.
#' @param s1 A vector of numerics drawn from a uniform distribution.
#' @examples
#' out <- .conj_uniform_sv(
#'   s1 = runif(10, 0, 100), cred.int.level = 0.95,
#'   plot = FALSE
#' )
#' lapply(out, head)
#' @keywords internal
#' @noRd
.conj_uniform_sv <- function(s1 = NULL, priors = NULL,
                             plot = FALSE, support = NULL, cred.int.level = NULL,
                             calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(scale = 0.5, location = 0.5)
  }
  #* `Update pareto prior with sufficient statistics`
  scale_prime <- priors$scale + length(s1)
  location_prime <- max(c(s1, priors$location), na.rm = TRUE)
  #* `Define support if it is missing`
  if (is.null(support) && calculatingSupport) {
    quantiles <- qpareto(c(0.0001, 0.9999), scale_prime, location_prime)
    return(quantiles)
  }
  #* `Make Posterior Draws`
  out$posteriorDraws <- extraDistr::rpareto(10000, scale_prime, location_prime)
  #* `posterior`
  dens1 <- extraDistr::dpareto(support, scale_prime, location_prime)
  pdf1 <- dens1 / sum(dens1)
  out$pdf <- pdf1
  # hde of location is calculated off of the posterior draws
  hde1 <- mean(as.numeric(hdi(out$posteriorDraws, ci = 0.01)[c("CI_low", "CI_high")]))
  # mode is location
  # mean is defined if scale > 1: scale_prime * location_prime / scale_prime - 1
  # median is location x root scale of 2
  hdi1 <- extraDistr::qpareto(
    c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
    scale_prime, location_prime
  )
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior <- list("scale" = scale_prime, "location" = location_prime)
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
