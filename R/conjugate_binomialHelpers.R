#' @description
#' Internal function for calculating \alpha and \beta of a beta distributed p (probability) parameter
#' in a beta-binomial conjugate distribution.
#'
#' @param s1 A named list of integer data containing number of successes and number of trials.
#' @examples
#' .conj_binomial_sv(
#'   s1 = list(successes = c(15, 14, 16, 11), trials = 20),
#'   priors = list(a = c(0.5, 0.5), b = c(0.5, 0.5)),
#'   plot = FALSE, cred.int.level = 0.95
#' )
#' @keywords internal
#' @noRd

.conj_binomial_sv <- function(s1 = NULL, priors = NULL,
                              plot = FALSE, support = NULL, cred.int.level = NULL,
                              calculatingSupport = FALSE) {
  #* `check stopping conditions`
  s1 <- .conj_binomial_formatter(s1)
  #* `separate data into counts and trials`
  s1_successes <- s1$successes
  s1_trials <- s1$trials
  #* `Replicate trials numbers if too short`
  if (length(s1_trials) < length(s1_successes)) {
    s1_trials <- rep(s1_trials, length(s1_successes))
  }

  #* `make default prior if none provided`
  #* `p parameter is beta distributed`
  if (is.null(priors)) {
    priors <- list(a = 0.5, b = 0.5)
  }
  #* `Define dense Support`
  #* `p parameter is beta distributed`
  if (is.null(support)) {
    if (calculatingSupport) {
      return(c(0.0001, 0.9999))
    }
    support <- seq(0.0001, 0.9999, 0.0001)
  }
  out <- list()
  #* `Update priors with observed counts`
  a1_prime <- priors$a[1] + sum(s1_successes)
  b1_prime <- priors$b[1] + sum(s1_trials) - sum(s1_successes)

  #* `calculate density over support``
  dens1 <- dbeta(support, a1_prime, b1_prime)
  pdf1 <- dens1 / sum(dens1)

  #* `calculate highest density interval`
  hdi1 <- qbeta(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), a1_prime, b1_prime)

  #* `calculate highest density estimate``
  if (a1_prime <= 1 && b1_prime > 1) {
    hde1 <- 0
  } else if (a1_prime > 1 && b1_prime <= 1) {
    hde1 <- 1
  } else {
    hde1 <- (a1_prime - 1) / (a1_prime + b1_prime - 2)
  }
  #* `save summary and parameters`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a1_prime
  out$posterior$b <- b1_prime
  #* `Make Posterior Draws`
  out$posteriorDraws <- rbeta(10000, a1_prime, b1_prime)
  out$pdf <- pdf1
  #* `keep data for plotting`
  if (plot) {
    out$plot_df <- data.frame(
      "range" = support, "prob" = pdf1,
      "sample" = rep("Sample 1", length(support))
    )
  }
  return(out)
}

#' Error condition handling helper function
#' @param s1 The sample given to conjugate and conjugate_binomial_sv
#' @keywords internal
#' @noRd

.conj_binomial_formatter <- function(s1) {
  #* `Check data for stopping conditions`
  if (any(unlist(s1) < 0)) {
    stop(paste0(
      "Binomial method requires successes and trials in sample data",
      "as a list of two integer vectors. See examples."
    ))
  }
  if (length(s1) != 2) {
    stop(paste0(
      "Binomial method requires successes and trials in sample data",
      "as a list of two integer vectors. See examples."
    ))
  }
  if (any(is.null(names(s1)), names(s1) != c("successes", "trials"))) {
    message("Assuming sample data is in order 'successes', 'trials'.")
    names(s1) <- c("successes", "trials")
  }
  return(s1)
}
