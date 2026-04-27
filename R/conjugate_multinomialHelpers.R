#' @description
#' Internal function for calculating the dirichlet distribution underlying
#' multinomially distributed single value traits (counts).
#'
#' Note that this function returns posterior draws and densities from the marginal beta
#' distributions, NOT from the dirichlet distribution.
#'
#' @param s1 A named vector/list of numerics drawn from a multinomial distribution.
#' A data.frame will trigger the mv option,
#' which will take column sums of the data then pass to the sv method.
#' @examples
#' out <- .conj_dirichlet_sv(
#'   s1 = list("A" = 10, "B" = 10, "C" = 5),
#'   cred.int.level = 0.95,
#'   plot = TRUE
#' )
#' lapply(out, head)
#'
#' s1 = list("A" = 10, "B" = 10, "C" = 5)
#' priors = NULL
#' support = seq(0.0001, 0.9999, length.out = 10000)
#' calculatingSupport = FALSE
#' cred.int.level = 0.95
#' hypothesis = "A > B"
#' sample_number = 1
#'
#'
#' @details
#' See Examples 1.4, 1.6, and 1.7 for thoughts on default dirichlet prior here
#' https://arxiv.org/pdf/1504.02689 , updating rule defined in
#' The Compendium of Conjugate Priors (https://www.johndcook.com/CompendiumOfConjugatePriors.pdf)
#' Section 5.1
#' @keywords internal
#' @noRd

.conj_multinomial_sv <- function(s1 = NULL, priors = NULL,
                                 support = NULL, cred.int.level = NULL,
                                 calculatingSupport = FALSE, ...) {
  nms <- names(s1)
  s1 <- unlist(s1)

  #* `Define support if it is missing`
  #* support is the same as beta binomial
  if (is.null(support) && calculatingSupport) {
    return(c(0.0001, 0.9999))
  }
  #* `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list("alpha" = rep(1 / length(s1), length(s1)))
  }
  names(priors$alpha) <- nms

  #* `update dirichlet prior with sufficient statistics (counts)`
  alpha_prime <- priors$alpha + s1
  #* `Make Posterior draws`
  #* Posterior draws used for ROPE comparisons mainly.
  #* need to consider if it's better to draw from rdirichlet or rbeta
  #* I think I want betas
  marginal_post <- do.call(cbind, lapply(seq_along(alpha_prime), function(i) {
    return(rbeta(10000, alpha_prime[i], sum(alpha_prime[-i])))
  }))

  #* `Calculate density over support`
  #* marginal beta densities, probably don't need all of them, just the ones
  #* specified in the hypothesis
  marginal_dens <- do.call(cbind, lapply(seq_along(alpha_prime), function(i) {
    d <- dbeta(support, alpha_prime[i], sum(alpha_prime[-i]))
    return(d / sum(d))
  }))

  #* `calculate HDI`
  #* note that this might change some since we could have HDE/HDI for every
  #* marginal beta
  hdis <- do.call(cbind, lapply(seq_along(alpha_prime), function(i) {
    hdi <- qbeta(
      c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
      alpha_prime[i], sum(alpha_prime[-i])
    )
    return(hdi)
  }))

  #* `calculate HDE`
  #* note that this might change some since we could have HDE/HDI for every
  #* marginal beta
  hdes <- do.call(cbind, lapply(seq_along(alpha_prime), function(i) {
    hde <- .betaHDE(alpha_prime[i], sum(alpha_prime[-i]))
    return(hde)
  }))
  colnames(hdes) <- colnames(hdis) <- colnames(marginal_dens) <- colnames(marginal_post) <- nms

  #* `Store summary`
  #* note that this might change some since we could have HDE/HDI for every
  #* marginal beta
  out <- list()
  out$summary <- data.frame(
    HDE_1 = as.numeric(hdes), HDI_1_low = hdis[1, ], HDI_1_high = hdis[2, ]
  )
  rownames(out$summary) <- nms
  out$posterior$alpha <- alpha_prime
  out$prior <- priors
  out$posteriorDraws <- marginal_post
  out$pdf <- marginal_dens

  #* `Save s1 data for plotting`
  #* note that this might change significantly from other distributions
  #* because we want to test/plot the marginal betas to keep the dimensions
  #* low enough to understand (and make the math tractable?)
  out$plot_list <- list(
    "range" = range(support),
    "ddist_fun" = "stats::dbeta",
    "priors" = list("shape1" = priors$alpha,
      "shape2" = setNames(unlist(lapply(seq_along(priors$alpha), function(i) {
        return(sum(priors$alpha[-i]))
      })), nms)
    ),
    "parameters" = list("shape1" = alpha_prime,
      "shape2" = setNames(unlist(lapply(seq_along(alpha_prime), function(i) {
        return(sum(alpha_prime[-i]))
      })), nms)
    ),
    "given" = NULL
  )
  return(out)
}

#' Multinomial "multi-value data", aka a dataframe of counts
#' @param s1 A data frame of counts, columns should represent categories
#' @keywords internal
#' @noRd

.conj_multinomial_mv <- function(s1 = NULL, priors = NULL,
                                 support = NULL, cred.int.level = NULL,
                                 calculatingSupport = FALSE, ...) {
  nms <- colnames(s1)
  s1 <- colSums(s1)
  names(s1) <- nms
  out <- .conj_multinomial_sv(s1, priors, support, cred.int.level, calculatingSupport)
  return(out)
}


#' @keywords internal
#' @noRd

.multinomial.parse.hypothesis <- function(hypothesis) {
  g1 <- trimws(gsub("(.*?)([>|<|=|!]{1,2})(.*)", "\\1", hypothesis))
  g2 <- trimws(gsub("(.*?)([>|<|=|!]{1,2})(.*)", "\\3", hypothesis))
  hyp <- trimws(gsub("(.*?)([>|<|=|!]{1,2})(.*)", "\\2", hypothesis))
  hyp <- switch({hyp},
    ">" = "greater", "<" = "lesser", "==" = "equal", "=" = "equal", "unequal"
  )
  return(list(g1, g2, hyp))
}


#' @keywords internal
#' @noRd

.multinomial.pdf.handling <- function(sample_results, hypothesis) {
  parsed <- .multinomial.parse.hypothesis(hypothesis)
  g1 <- parsed[[1]]
  g2 <- parsed[[2]]
  hyp <- parsed[[3]]
  if (length(sample_results) == 2) {
    pdf.handling.output <- .post.prob.from.pdfs(
      sample_results[[1]]$pdf[, g1],
      sample_results[[2]]$pdf[, g2],
      hyp
    )
  } else {
    pdf.handling.output <- .post.prob.from.pdfs(
      sample_results[[1]]$pdf[, g1],
      sample_results[[1]]$pdf[, g2],
      hyp
    )
  }
  return(
    list(
      "pdf.handling.output" = pdf.handling.output,
      "hyp" = hyp
    )
  )
}


#' @keywords internal
#' @noRd
.multinomial.conj.plot.format <- function(res) {
  hypl <- .multinomial.parse.hypothesis(res$call$hypothesis)
  g1 <- hypl[[1]]
  g2 <- hypl[[2]]
  new_res <- res
  if (length(res$plot_parameters) == 2) {
    # keep only relevant marginal beta parameters
    new_res$plot_parameters[[1]]$parameters$shape1 <- res$plot_parameters[[1]]$parameters$shape1[g1]
    new_res$plot_parameters[[2]]$parameters$shape1 <- res$plot_parameters[[2]]$parameters$shape1[g2]
    new_res$plot_parameters[[1]]$parameters$shape2 <- res$plot_parameters[[1]]$parameters$shape2[g1]
    new_res$plot_parameters[[2]]$parameters$shape2 <- res$plot_parameters[[2]]$parameters$shape2[g2]
    # keep only relevant marginal priors
    new_res$plot_parameters[[1]]$priors$shape1 <- res$plot_parameters[[1]]$priors$shape1[g1]
    new_res$plot_parameters[[2]]$priors$shape1 <- res$plot_parameters[[2]]$priors$shape1[g2]
    new_res$plot_parameters[[1]]$priors$shape2 <- res$plot_parameters[[1]]$priors$shape2[g1]
    new_res$plot_parameters[[2]]$priors$shape2 <- res$plot_parameters[[2]]$priors$shape2[g2]
    # subset summary to draw HDE/HDI correctly
    new_res$summary <- cbind(
      res$summary[g1, grepl("_1", colnames(res$summary))],
      res$summary[g2, grepl("_2", colnames(res$summary))],
      res$summary[1, !grepl("[1|2]", colnames(res$summary))]
    )
  } else {
    # duplicate plot parameters to 2L
    new_res$plot_parameters <- list(new_res$plot_parameters[[1]], new_res$plot_parameters[[1]])
    new_res$data <- list(new_res$data[[1]], new_res$data[[1]])
    # keep only relevant marginal beta parameters
    new_res$plot_parameters[[1]]$parameters$shape1 <- new_res$plot_parameters[[1]]$parameters$shape1[g1]
    new_res$plot_parameters[[2]]$parameters$shape1 <- new_res$plot_parameters[[2]]$parameters$shape1[g2]
    new_res$plot_parameters[[1]]$parameters$shape2 <- new_res$plot_parameters[[1]]$parameters$shape2[g1]
    new_res$plot_parameters[[2]]$parameters$shape2 <- new_res$plot_parameters[[2]]$parameters$shape2[g2]
    # keep only relevant marginal priors
    new_res$plot_parameters[[1]]$priors$shape1 <- new_res$plot_parameters[[1]]$priors$shape1[g1]
    new_res$plot_parameters[[2]]$priors$shape1 <- new_res$plot_parameters[[2]]$priors$shape1[g2]
    new_res$plot_parameters[[1]]$priors$shape2 <- new_res$plot_parameters[[1]]$priors$shape2[g1]
    new_res$plot_parameters[[2]]$priors$shape2 <- new_res$plot_parameters[[2]]$priors$shape2[g2]
    # subset summary to draw HDE/HDI correctly
    g1_cols <- res$summary[g1, grepl("_1", colnames(res$summary))]
    g2_cols <- res$summary[g2, grepl("_1", colnames(res$summary))]
    colnames(g2_cols) <- gsub("_1", "_2", colnames(g2_cols))
    new_res$summary <- cbind(
      g1_cols,
      g2_cols,
      res$summary[1, !grepl("[1|2]", colnames(res$summary))]
    )
  }
  return(new_res)
}


#' @keywords internal
#' @noRd

.multinomial.rope.format <- function(sample_results, hypothesis) {
  hypl <- .multinomial.parse.hypothesis(hypothesis)
  g1 <- hypl[[1]]
  g2 <- hypl[[2]]
  post1 <- sample_results[[1]]$posteriorDraws[, g1]
  #* for ROPE comparison everything expects 2L results, so force that.
  if (length(sample_results) == 1) {
    post2 <- sample_results[[1]]$posteriorDraws[, g2]
    new_sample_results <- list(sample_results, sample_results)
  } else {
    post2 <- sample_results[[2]]$posteriorDraws[, g2]
    new_sample_results <- sample_results
  }
  new_sample_results[[1]]$posteriorDraws <- post1
  new_sample_results[[2]]$posteriorDraws <- post2
  return(new_sample_results)
}
