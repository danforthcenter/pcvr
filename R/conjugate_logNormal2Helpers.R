#' @description
#' Internal function for calculating \alpha and \beta of the gamma distributed variance of lognormal
#' data given an estimate of the lognormal \mu obtained via the method of moments using multi value
#' traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number
#' representing the "bin".
#'
#' @examples
#' if (FALSE) {
#'   makeMvLn <- function(bins = 500, mu_log, sigma_log) {
#'     setNames(data.frame(matrix(hist(rlnorm(2000, mu_log, sigma_log),
#'       breaks = seq(1, bins, 5),
#'       plot = FALSE
#'     )$counts, nrow = 1)), paste0("b", seq(1, bins, 5))[-1])
#'   }
#'   set.seed(123)
#'   mv_ln <- rbind(
#'     do.call(rbind, lapply(1:30, function(i) {
#'       makeMvLn(mu_log = log(130), sigma_log = log(1.3))
#'     })),
#'     do.call(rbind, lapply(1:30, function(i) {
#'       makeMvLn(mu_log = log(100), sigma_log = log(1.2))
#'     }))
#'   )
#'
#'   .conj_lognormal_mv(
#'     s1 = mv_ln[1:30, ], s2 = mv_ln[31:60, ],
#'     priors = list(mu_log = c(log(10), log(10)), n = c(1, 1), sigma_log = c(log(3), log(3))),
#'     plot = FALSE, rope_range = NULL, rope_ci = 0.89,
#'     cred.int.level = 0.89, hypothesis = "equal", support = NULL
#'   )
#' }
#' @keywords internal
#' @noRd

.conj_lognormal2_mv <- function(s1 = NULL, priors = NULL,
                               plot = FALSE, support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 1, b = 1) # prior on shape, scale of precision
  }
  #* `Standardize sample 1 class and names`
  if (is.null(colnames(s1))) {
    bins <- seq_along(s1)
    colnames(s1) <- paste0("b", bins)
    warning(paste0("Assuming unnamed columns represent bins from ", min(bins), " to ", max(bins)))
  }
  if (is.matrix(s1)) {
    s1 <- as.data.frame(s1)
  }
  #* `Reorder columns if they are not in the numeric order`
  histColsBin <- as.numeric(sub("[a-zA-Z_.]+", "", colnames(s1)))
  bins_order <- sort(histColsBin, index.return = TRUE)$ix
  s1 <- s1[, bins_order]
  #* `Loop over reps, get moments for each histogram`
  rep_distributions <- lapply(seq_len(nrow(s1)), function(i) {
    X1 <- rep(histColsBin[bins_order], as.numeric(s1[i, ]))
    #* `Get mean of x1`
    x_bar <- mean(X1)
    mu_s1 <- log(x_bar / (sqrt(var(X1) / x_bar^2) +1) )
    #* `Get sigma of s1`
    sigma_s1 <- sqrt(log((var(X1)) / (x_bar ^ 2) + 1))
    prec_s1 <- 1 / (sigma_s1 ^ 2)
    #* `Update Gamma Distribution of precision`
    #* sufficient stats: n, ss
    ss <- nrow(s1) * mean((log(X1) - mu_s1)^2) # mean * nrow instead of sum for MV traits
    n1 <- nrow(s1)
    a_prime <- priors$a + (n1 / 2)
    b_prime <- priors$b + (ss / 2)
    return(list("a_prime" = a_prime, "b_prime" = b_prime, "ln_mu" = mu_s1))
  })
  #* `Unlist parameters`
  n1 <- nrow(s1)
  a_prime <- mean(unlist(lapply(rep_distributions, function(i) {
    i$a_prime
  })))
  b_prime <- mean(unlist(lapply(rep_distributions, function(i) {
    i$b_prime
  })))
  ln_mu_prime <- mean(unlist(lapply(rep_distributions, function(i) {
    i$ln_mu
  })))
  #* `Define support if it is missing`
  if (is.null(support)) {
    quantiles <- qgamma(c(0.0001, 0.9999), shape = a_prime, scale = b_prime)
    if (calculatingSupport) {
      return(quantiles)
    }
    support <- seq(quantiles[1], quantiles[2], length.out = 10000)
  }
  #* `posterior`
  dens1 <- dgamma(support, shape = a_prime, scale = b_prime)
  pdf1 <- dens1 / sum(dens1)
  if (b_prime <= 1 && a_prime > 1) {
    hde1 <- qgamma(0.5, shape = a_prime, scale = b_prime)
  } else if (b_prime == 0) {
    hde1 <- qgamma(0.5, shape = a_prime, scale = b_prime)
  } else {
    hde1 <- (b_prime - 1) * a_prime # note, using shape instead of rate (inverse) HDE
  }
  hdi1 <- qgamma(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
                 shape = a_prime, scale = b_prime)
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a_prime
  out$posterior$b <- b_prime
  out$posterior$lognormal_mu <- ln_mu_prime # returning this as a number, not a distribution
  #* `Make Posterior Draws`
  out$posteriorDraws <- rgamma(10000, shape = a_prime, scale = b_prime)
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
#' Internal function for calculating \alpha and \beta of the gamma distributed precision of lognormal
#' data given an estimate of the lognormal \mu obtained via the method of moments using single value
#' traits.
#'
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @examples
#' if (FALSE) {
#'   .conj_lognormal_sv(
#'     s1 = rlnorm(100, log(130), log(1.3)), s2 = rlnorm(100, log(100), log(1.6)),
#'     priors = list(a = 1, b = 1),
#'     plot = FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'     cred.int.level = 0.89, hypothesis = "equal", support = NULL
#'   )
#' }
#' @keywords internal
#' @noRd

.conj_lognormal2_sv <- function(s1 = NULL, priors = NULL,
                               plot = FALSE, support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(a = 1, b = 1) # prior on shape, scale of precision
  }
  if (length(s1) > 1) {
    #* `Get mean of s1`
    x_bar <- mean(s1)
    mu_s1 <- log(x_bar / (sqrt(var(s1) / x_bar^2) +1) )
    #* `Get sigma of s1`
    sigma_s1 <- sqrt(log((var(s1)) / (x_bar ^ 2) + 1))
    prec_s1 <- 1 / (sigma_s1 ^ 2)
    #* `Update Gamma Distribution of precision`
    #* sufficient stats: n, ss
    ss <- sum((log(s1) - mu_s1)^2)
    n1 <- length(s1)
    a_prime <- priors$a + (n1 / 2)
    b_prime <- priors$b + (ss / 2)
    #* `Define support if it is missing`
    if (is.null(support)) {
      quantiles <- qgamma(c(0.0001, 0.9999), shape = a_prime, scale = b_prime)
      if (calculatingSupport) {
        return(quantiles)
      }
      support <- seq(quantiles[1], quantiles[2], length.out = 10000)
    }
    #* `posterior`
    dens1 <- dgamma(support, shape = a_prime, scale = b_prime)
    pdf1 <- dens1 / sum(dens1)
    if (b_prime <= 1 && a_prime > 1) {
      hde1 <- qgamma(0.5, shape = a_prime, scale = b_prime)
    } else if (b_prime == 0) {
      hde1 <- qgamma(0.5, shape = a_prime, scale = b_prime)
    } else {
      hde1 <- (b_prime - 1) * a_prime # note, using shape instead of rate (inverse) HDE
    }
    hdi1 <- qgamma(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
                  shape = a_prime, scale = b_prime)
    #* `Store summary`
    out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
    out$posterior$a <- a_prime
    out$posterior$b <- b_prime
    out$posterior$lognormal_mu <- mu_s1 # returning this as a number, not a distribution
    #* `Make Posterior Draws`
    out$posteriorDraws <- rgamma(10000, shape = a_prime, scale = b_prime)
    out$pdf <- pdf1
    #* `save s1 data for plotting`
    if (plot) {
      out$plot_df <- data.frame(
        "range" = support,
        "prob" = pdf1,
        "sample" = rep("Sample 1", length(support))
      )
    }
  } else {
    stop("s1 must be a numeric of length 2 or greater")
  }
  return(out)
}
