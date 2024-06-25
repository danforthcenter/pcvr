#' @description
#' Internal function for Bayesian comparison of log-normal data represented by multi value traits.
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

.conj_lognormal_mv <- function(s1 = NULL, priors = NULL,
                               plot = FALSE, support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(mu_log = log(10), n = 1, sigma_log = log(3))
  }
  #* `Standardize sample 1 class and names`
  if (is.null(colnames(s1))) {
    bins <- (seq_along(s1))
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
    xbar_1 <- mean(X1)
    s2_1 <- var(X1) / (nrow(s1) - 1)
    #* `Add prior distribution and lognormal method of moments`
    mu_ls <- log(xbar_1 / sqrt((s2_1 / xbar_1^2) + 1))
    sigma_ls <- sqrt(log((s2_1 / xbar_1^2) + 1))
    return(list("mu_log" = mu_ls, "sigma_log" = sigma_ls))
  })
  n1 <- nrow(s1)
  mu_ls <- unlist(lapply(rep_distributions, function(i) {i$mu_log}))
  sigma_ls <- unlist(lapply(rep_distributions, function(i) {i$sigma_log}))

  n_prime <- n1 + priors$n[1]
  mu_ls_prime <- (sum(mu_ls) + (priors$mu_log[1] * priors$n[1])) / n_prime
  sigma_ls_prime <- (sum(sigma_ls) + (priors$sigma_log[1] * priors$n[1])) / n_prime

  # #* `Turn s1 matrix into a vector`
  # X1 <- rep(histColsBin[bins_order], as.numeric(round(colSums(s1))))
  # 
  # #* `Get mean and variance from s1`
  # xbar_1 <- mean(X1)
  # s2_1 <- var(X1) / (nrow(s1) - 1)
  # n1 <- nrow(s1)
  # #* `Add prior distribution and lognormal method of moments`
  # mu_ls <- log(xbar_1 / sqrt((s2_1 / xbar_1^2) + 1))
  # sigma_ls <- sqrt(log((s2_1 / xbar_1^2) + 1))
  # 
  # n_prime <- n1 + priors$n[1]
  # mu_ls_prime <- ((mu_ls * n1) + (priors$mu_log[1] * priors$n[1])) / n1_n
  # sigma_ls_prime <- ((sigma_ls * n1) + (priors$sigma_log[1] * priors$n[1])) / n1_n
  #* `Define support if it is missing`
  if (is.null(support)) {
    quantiles <- qlnorm(c(0.0001, 0.9999), mu_ls_prime, sigma_ls_prime)
    if (calculatingSupport) {
      return(quantiles)
    }
    support <- seq(quantiles[1], quantiles[2], length.out = 10000)
  }
  #* `posterior`
  dens1 <- dlnorm(support, mu_ls_prime, sigma_ls_prime)
  pdf1 <- dens1 / sum(dens1)
  hde1 <- exp(mu_ls_prime)
  hdi1 <- qlnorm(c((1 - cred.int.level) / 2,
                   (1 - ((1 - cred.int.level) / 2))),
                 mu_ls_prime, sigma_ls_prime)
  #* `Store summary`
  out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$mu_log <- mu_ls_prime
  out$posterior$n <- n1_prime
  out$posterior$sigma_log <- sigma_ls_prime
  #* `Make Posterior Draws`
  out$posteriorDraws <- rlnorm(10000, mu_ls_prime, sigma_ls_prime)
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
#' Internal function for Bayesian comparison of log-normal data represented by single value traits.
#' Lognormal method of moments
#'
#' \bar{x} \sim e^{\mu + \sigma^2 / 2}
#' \s^2 \sim (e^{\alpha^2} -1) \cdot e^{2 \cdot \mu + \sigma^2}
#'
#' Calculate Sigma:
#'
#' \bar{x}^2 \sim e^{2 \cdot \mu + \sigma^2}
#' s^2 \sim (e^\alpha^2 -1) \cdot \bar{x}^2
#' (s^2 / \bar{x}^2) +1 \sim e^\alpha^2
#' \alpha^2 \sim ln((s^2 / \bar{x}^2) +1 )
#' \alpha \sim \sqrt{ln((s^2 / \bar{x}^2) +1 )}
#'
#' Calculate Mu:
#'
#' \bar{x}^2 \sim e^{2 \cdot \mu + \sigma^2}
#' 2 \cdot ln(\bar{x}) \sim 2 \cdot \mu + \sigma^2
#' \mu \sim ln(\bar{x} / \sqrt{(s^2 / \bar{x}^2) +1 })
#'
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @examples
#' if (FALSE) {
#'   .conj_lognormal_sv(
#'     s1 = rlnorm(100, log(130), log(1.3)), s2 = rlnorm(100, log(100), log(1.6)),
#'     priors = list(mu_log = c(log(10), log(10)), n = c(1, 1), sigma_log = c(log(3), log(3))),
#'     plot = FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'     cred.int.level = 0.89, hypothesis = "equal", support = NULL
#'   )
#' }
#' @keywords internal
#' @noRd

.conj_lognormal_sv <- function(s1 = NULL, priors = NULL,
                               plot = FALSE, support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list(mu_log = log(10), n = 1, sigma_log = log(3))
  }
  #* `Get mean and variance from s1`
  if (length(s1) > 1) {
    xbar_1 <- mean(s1)
    s2_1 <- var(s1) / (length(s1) - 1)
    n1 <- length(s1)
    #* `Add prior distribution and lognormal method of moments`
    mu_ls <- log(xbar_1 / sqrt((s2_1 / xbar_1^2) + 1))
    sigma_ls <- sqrt(log((s2_1 / xbar_1^2) + 1))
    n1_n <- n1 + priors$n[1]
    mu_ls_n <- ((mu_ls * n1) + (priors$mu_log[1] * priors$n[1])) / n1_n
    sigma_ls_n <- ((sigma_ls * n1) + (priors$sigma_log[1] * priors$n[1])) / n1_n
    #* `Define support if it is missing`
    if (is.null(support)) {
      quantiles <- qlnorm(c(0.0001, 0.9999), mu_ls_n, sigma_ls_n)
      if (calculatingSupport) {
        return(quantiles)
      }
      support <- seq(quantiles[1], quantiles[2], length.out = 10000)
    }
    #* `posterior`
    dens1 <- dlnorm(support, mu_ls_n, sigma_ls_n)
    pdf1 <- dens1 / sum(dens1)
    hde1 <- exp(mu_ls_n)
    hdi1 <- qlnorm(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))), mu_ls_n, sigma_ls_n)
    #* `Store summary`
    out$summary <- data.frame(HDE_1 = hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
    out$posterior$mu_log <- mu_ls_n
    out$posterior$n <- n1_n
    out$posterior$sigma_log <- sigma_ls_n
    #* `Make Posterior Draws`
    out$posteriorDraws <- rlnorm(10000, mu_ls_n, sigma_ls_n)
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
