#' @description
#' Internal function for calculating \alpha and \beta of a distribution represented by multi value
#' traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number
#' between 0.0001 and 0.9999 representing the "bin".
#' @examples
#'
#' mv_beta <- mvSim(
#'   dists = list(
#'     rbeta = list(shape1 = 5, shape2 = 8),
#'   ),
#'   n_samples = c(30)
#' )
#' .conj_beta_mv(
#'   s1 = mv_beta[1:30, -1], priors = list(a = c(0.5), b = c(0.5)),
#'   cred.int.level = 0.9,
#'   plot = TRUE
#' )
#'
#' @keywords internal
#' @noRd
.conj_beta_mv <- function(s1 = NULL, priors = NULL,
                          plot = FALSE, support = NULL, cred.int.level = NULL,
                          calculatingSupport = FALSE) {
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 0.5, b = 0.5)
  }
  #* `Define dense Support`
  if (is.null(support) && calculatingSupport) {
    return(c(0.0001, 0.9999))
  }
  out <- list()
  #* `Reorder columns if they are not in the numeric order`
  histColsBin <- as.numeric(sub("[a-zA-Z_.]+", "", colnames(s1)))
  if (any(histColsBin > 1 || histColsBin < 0)) {
    stop("Beta Distribution is only defined on [0,1]")
  }
  bins_order <- sort(histColsBin, index.return = TRUE)$ix
  s1 <- s1[, bins_order]

  #* `Turn matrix into a vector`
  X1 <- rep(histColsBin[bins_order], as.numeric(round(colSums(s1))))

  #* `get parameters for s1 using method of moments``
  #* y ~ Beta(\alpha, \beta)
  #* \alpha ~ \bar{y}( ( (\bar{y} * (1-\bar{y}))/\bar(var) )-1 )
  #* \beta ~ (1-\bar{y})( ( (\bar{y} * (1-\bar{y}))/\bar(var) )-1 )
  mu1 <- mean(X1) #' \bar{y}
  nu1 <- var(X1) / (nrow(s1) - 1) #' \bar{var} the unbiased sample variance
  alpha1 <- mu1 * ((mu1 * (1 - mu1)) / (nu1) - 1)
  beta1 <- (1 - mu1) * ((mu1 * (1 - mu1)) / (nu1) - 1)

  #* `Add priors`
  a1_prime <- alpha1 + priors$a[1]
  b1_prime <- beta1 + priors$b[1]

  #* `calculate density`
  dens1 <- dbeta(support, a1_prime, b1_prime)
  pdf1 <- dens1 / sum(dens1)

  #* `calculate highest density interval`
  hdi1 <- qbeta(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), a1_prime, b1_prime)

  #* `calculate highest density estimate`
  hde1 <- .betaHDE(a1_prime, b1_prime)

  #* `save summary and parameters`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a1_prime
  out$posterior$b <- b1_prime
  #* `Make Posterior Draws`
  out$posteriorDraws <- rbeta(10000, a1_prime, b1_prime)
  out$pdf <- pdf1
  #* `keep data for plotting`
  if (plot) {
    out$plot_df <- data.frame("range" = support, "prob" = pdf1,
                              "sample" = rep("Sample 1", length(support)))
  }

  return(out)
}



#' @description
#' Internal function for calculating \alpha and \beta of a distribution represented by multi value
#' traits.
#' @param s1 A vector of numerics drawn from a beta distribution.
#' @examples
#' .conj_beta_sv(
#'   s1 = rbeta(100, 5, 10),
#'   priors = list(a = c(0.5, 0.5), b = c(0.5, 0.5)),
#'   cred.int.level = 0.9,
#'   plot = FALSE
#' )
#' @keywords internal
#' @noRd
.conj_beta_sv <- function(s1 = NULL, priors = NULL,
                          plot = FALSE, support = NULL, cred.int.level = NULL,
                          calculatingSupport = FALSE) {
  if (any(c(s1) > 1 || c(s1) < 0)) {
    stop("Beta Distribution is only defined on [0,1]")
  }
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 0.5, b = 0.5)
  }
  #* `Define dense Support`
  if (is.null(support) && calculatingSupport) {
    return(c(0.0001, 0.9999))
  }
  out <- list()

  #* `get parameters for s1 using method of moments``
  #* y ~ Beta(\alpha, \beta)
  #* \alpha ~ \bar{y}( ( (\bar{y} * (1-\bar{y}))/\bar(var) )-1 )
  #* \beta ~ (1-\bar{y})( ( (\bar{y} * (1-\bar{y}))/\bar(var) )-1 )
  mu1 <- mean(s1) #' \bar{y}
  nu1 <- var(s1) / (length(s1) - 1) #' \bar{var} the unbiased sample variance
  alpha1 <- mu1 * ((mu1 * (1 - mu1)) / (nu1) - 1)
  beta1 <- (1 - mu1) * ((mu1 * (1 - mu1)) / (nu1) - 1)

  #* `Add priors in`
  a1_prime <- priors$a[1] + alpha1
  b1_prime <- priors$b[1] + beta1

  #* `calculate density over support``
  dens1 <- dbeta(support, a1_prime, b1_prime)
  pdf1 <- dens1 / sum(dens1)

  #* `calculate highest density interval`
  hdi1 <- qbeta(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), a1_prime, b1_prime)

  #* `calculate highest density estimate`
  hde1 <- .betaHDE(a1_prime, b1_prime)

  #* `save summary and parameters`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a1_prime
  out$posterior$b <- b1_prime
  #* `Make Posterior Draws`
  out$posteriorDraws <- rbeta(10000, a1_prime, b1_prime)
  out$pdf <- pdf1
  #* `keep data for plotting`
  if (plot) {
    out$plot_df <- data.frame("range" = support, "prob" = pdf1,
                              "sample" = rep("Sample 1", length(support)))
  }
  return(out)
}

#' @description
#' Internal function for calculating the HDE of a beta distribution
#' @param alpha alpha parameter
#' @param beta beta parameter
#' @examples
#' .betaHDE(1, 2)
#' .betaHDE(2, 1)
#' .betaHDE(10, 10)
#' @keywords internal
#' @noRd

.betaHDE <- function(alpha, beta) {
  if (alpha <= 1 && beta > 1) {
    hde <- 0
  } else if (alpha > 1 && beta <= 1) {
    hde <- 1
  } else {
    hde <- (alpha - 1) / (alpha + beta - 2)
  }
  return(hde)
}
