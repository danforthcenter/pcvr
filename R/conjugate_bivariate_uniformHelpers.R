#' @description
#' Internal function for calculating the pareto distribution of the upper boundary of a uniform
#' distribution represented by single value traits.
#' @param s1 A vector of numerics drawn from a uniform distribution.
#' @examples
#' if (FALSE) {
#'   out <- .conj_bivariate_uniform_sv(
#'     s1 = runif(10, -20, 100), cred.int.level = 0.95,
#'     plot = FALSE
#'   )
#'   lapply(out, head)
#' }
#' @keywords internal
#' @noRd

.conj_bivariate_uniform_sv <- function(s1 = NULL, priors = NULL,
                                       plot = FALSE, support = NULL, cred.int.level = NULL,
                                       calculatingSupport = FALSE) {
  out <- list()
  #* `make default prior if none provided`
  #* conjugate prior needs r1, r2, and alpha
  #* which are locations and a shared shape (scale) paramter
  if (is.null(priors)) {
    priors <- list(location_l = 0, location_u = 1, scale = 1)
  }
  #* `Update bivariate pareto prior with sufficient statistics`
  location_u_prime <- max(c(s1, priors$location_u[1]))
  location_l_prime <- min(c(s1, priors$location_l[1]))
  scale_prime <- priors$scale[1] + length(s1)
  #* `Define bivariate support if it is missing`
  #* `this will require some thought since there are two directions.`
  #* first problem is that qpareto requires the parameters to be > 0,
  #* but the lower boundary can be negative, or could be something like 1
  #* with a tail that becomes negative. I guess changing the center is the
  #* way to account for those potential problems.
  if (is.null(support)) {
    (quantiles_u <- qpareto(c(0.0001, 0.9999), scale_prime, abs(location_u_prime)))
    (quantiles_l <- qpareto(c(0.0001, 0.9999), scale_prime, abs(location_l_prime)))
    support_l <- seq(quantiles_l[1], quantiles_l[2], length.out = 10000)
    support_u <- seq(quantiles_u[1], quantiles_u[2], length.out = 10000)

    if (location_l_prime < 0) {
      quantiles_l <- -1 * rev(quantiles_l)
      support_l <- seq(quantiles_l[1], quantiles_l[2], length.out = 10000)
    }
    if (location_u_prime < 0) {
      quantiles_u <- -1 * rev(quantiles_u)
      support_u <- seq(quantiles_u[1], quantiles_u[2], length.out = 10000)
    }

    if (calculatingSupport) {
      return(list("A" = quantiles_l, "B" = quantiles_u))
    }
  } else {
    support_l <- support$A
    support_u <- support$B
  }
  #* `Make Posterior Draws`
  out$posteriorDraws <- .conj_cond_inv_rpareto(
    10000, location_l_prime, location_u_prime,
    scale_prime
  )
  #* `posterior`
  #* this also needs to handle the possibility of negative locations
  if (location_l_prime < 0) {
    dens_l <- extraDistr::dpareto(abs(support_l), scale_prime, abs(location_l_prime))
  } else {
    dens_l <- rev(extraDistr::dpareto(support_l, scale_prime, location_l_prime))
  }

  if (location_u_prime < 0) {
    dens_u <- extraDistr::dpareto(abs(support_u), scale_prime, abs(location_u_prime))
  } else {
    dens_u <- extraDistr::dpareto(support_u, scale_prime, location_u_prime)
  }

  pdf_l <- dens_l / sum(dens_l)
  pdf_u <- dens_u / sum(dens_u)
  out$pdf <- list("A" = pdf_l, "B" = pdf_u)

  hde_l <- location_l_prime
  hde_u <- location_u_prime

  if (location_l_prime < 0) {
    hdi_l <- -1 * rev(extraDistr::qpareto(
      c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
      scale_prime, abs(location_l_prime)
    ))
  } else {
    hdi_l <- extraDistr::qpareto(
      c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
      scale_prime, location_l_prime
    )
  }
  if (location_u_prime < 0) {
    hdi_u <- -1 * rev(extraDistr::qpareto(
      c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
      scale_prime, abs(location_u_prime)
    ))
  } else {
    hdi_u <- extraDistr::qpareto(
      c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
      scale_prime, location_u_prime
    )
  }

  #* `Store summary`
  out$summary <- data.frame(
    HDE_1 = c(hde_l, hde_u),
    HDI_1_low = c(hdi_l[1], hdi_u[1]),
    HDI_1_high = c(hdi_l[2], hdi_u[2]),
    param = c("A", "B")
  )
  out$posterior <- list(
    "scale" = scale_prime, "location_l" = location_l_prime,
    "location_u" = location_u_prime
  )
  #* `save s1 data for plotting`
  if (plot) {
    out$plot_df <- data.frame(
      "range" = c(support_l, support_u),
      "prob" = c(pdf_l, pdf_u),
      "param" = rep(c("A", "B"), each = length(support_u)),
      "sample" = rep("Sample 1", 2 * length(support_u))
    )
  }
  return(out)
}

.conj_cond_inv_rpareto <- function(n, r1, r2, scale) {
  u <- runif(n, min = 0, max = 1)
  # pareto quantile function
  x2 <- r2 / (u^(1 / scale))
  u <- runif(n, min = 0, max = 1)
  # this is a displaced origin pareto
  # also using quantile function of the marginal x1 | x2
  x1 <- r1 + (r1 / r2) * x2 * (1 / (u^(1 / (scale + 1))) - 1)
  return(cbind.data.frame("A" = x1, "B" = x2))
}
