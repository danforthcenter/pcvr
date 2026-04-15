#' @description
#' Internal function for calculating the dirichlet distribution underlying
#' multinomially distributed single value traits (counts).
#'
#' Note that this function returns posterior draws and densities from the marginal beta
#' distributions, NOT from the dirichlet distribution.
#'
#' @param s1 A named vector/list of numerics drawn from a multinomial distribution. A data.frame will trigger the mv option, which will take column sums of the data then pass to the sv method.
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

  #' `Define support if it is missing`
  #' support is the same as beta binomial
  if (is.null(support) && calculatingSupport) {
    return(c(0.0001, 0.9999))
  }
  #' `Check stopping conditions (similar to binomial)`
  
  #' `make default prior if none provided`
  if (is.null(priors)) {
    priors <- list("alpha" = rep(1 / length(s1), length(s1)))
  }
  #' I think the default prior should be a flat dirichlet
  #' lowest entropy would mean very low precision
  #' interesting comments on jeffrey's prior here: https://arxiv.org/pdf/1504.02689
  #` NOTE see examples 1.4 and 1.6 in linked paper for criticism of the jeffrey prior
  #' in this context given large K and small N (which is probably an unlikely scenario
  #' with the data that I'm designing this for, so we could just mention in the docs)
  #' They actually come around to the explicit dirichlet example in ex 1.7 (although
  #' on first reading I thought they were talking about dirichlet's the entire time?)
  #' where `Alpha_i = 1/K for all i in K`.
  #' I think this would be a fine prior for our settings, it's going to be very low
  #' precision and flat. Besides, the use case we have in mind is:
  #' Prior_0 is updated to Posterior_Before by Before-intervention data
  #' Posterior_Before is updated to Posterior_After by After-Intervention data
  #' Marginal Betas from Posterior_After and Posterior_Before are compared.
  #'
  #' NOTE as a sample input thing I don't think we are building for that scenario,
  #' the way I imagine using that would be:
  #' S1 = before_data, S2 = before_data + after_data
  #' so that the S2 group is effectively updated for the before data, then again for
  #' the after data.

  #' `update dirichlet prior with sufficient statistics (counts)`
  alpha_prime <- priors$alpha + s1
  
  #' `Parse hypothesis to get marginal distributions`
  #' We don't need to parse the hypothesis here since we'll return the entire set
  #' of marginal posterior draws and densities.
  

  #' `Make Posterior draws`
  #' Posterior draws used for ROPE comparisons mainly.
  #' need to consider if it's better to draw from rdirichlet or rbeta
  #' I think I want betas
  marginal_post <- do.call(cbind, lapply(seq_along(alpha_prime), function(i) {
    return(rbeta(10000, alpha_prime[i], sum(alpha_prime[-i])))
  }))

  #' `Calculate density over support`
  #' marginal beta densities, probably don't need all of them, just the ones
  #' specified in the hypothesis
  marginal_dens <- do.call(cbind, lapply(seq_along(alpha_prime), function(i) {
    d <- dbeta(support, alpha_prime[i], sum(alpha_prime[-i]))
    return(d / sum(d))
  }))

  #' `calculate HDI`
  #' note that this might change some since we could have HDE/HDI for every
  #' marginal beta
  hdis <- do.call(cbind, lapply(seq_along(alpha_prime), function(i) {
    hdi <- qbeta(c((1 - cred.int.level) / 2, (1 - ((1 - cred.int.level) / 2))),
      alpha_prime[i], sum(alpha_prime[-i]))
    return(hdi)
  }))

  #' `calculate HDE`
  #' note that this might change some since we could have HDE/HDI for every
  #' marginal beta
  hdes <- do.call(cbind, lapply(seq_along(alpha_prime), function(i) {
    hde <- .betaHDE(alpha_prime[i], sum(alpha_prime[-i]))
    return(hde)
  }))
  colnames(hdes) <- colnames(hdis) <- colnames(marginal_dens) <- colnames(marginal_post) <- nms

  #' `Store summary`
  #' note that this might change some since we could have HDE/HDI for every
  #' marginal beta
  out <- list()
  out$summary <- data.frame(
    HDE_1 = as.numeric(hdes), HDI_1_low = hdis[1, ], HDI_1_high = hdis[2, ]
  )
  rownames(out$summary) <- nms
  out$posterior$alpha <- alpha_prime
  out$prior <- priors
  out$posteriorDraws <- marginal_post
  out$pdf <- marginal_dens

  #' `Save s1 data for plotting`
  #' note that this might change significantly from other distributions
  #' because we want to test/plot the marginal betas to keep the dimensions
  #' low enough to understand (and make the math tractable?)
  out$plot_list <- list(
    "range" = range(support),
    "ddist_fun" = "stats::dbeta",
    "priors" = list("shape1" = priors$alpha,
      "shape2" = unlist(lapply(seq_along(priors$alpha), function(i) {
        return(sum(priors$alpha[-i]))
      }))
    ),
    "parameters" = list("shape1" = alpha_prime,
      "shape2" = unlist(lapply(seq_along(alpha_prime), function(i) {
        return(sum(alpha_prime[-i]))
      }))
    ),
    "given" = NULL
  )
  return(out)
}

#' Might define the plot helper functions here as well, I don't think I have anything else
#' that is quite in this situation for how it's using marginals, that's why I have the
#' bivariate options for some distributions but a shared bivariate plotting function.
#' Overall example would be drawn from `conjugate_plot::.conj_general_plot` and might just
#' be pre-processing.
#' Still need to decide how those things are done, thinking that it would be specified in
#' the hypothesis, so instead of `(s1) ">" (s2)` I'd have the input be `(s1)K_1 > (s2)K_3`?
#' I think that if there are 2 samples then it would assume your hypothesis has to be between
#' the samples (previous example) but if there is only 1 sample that you're comparing categories
#' in that sample.
#' 
#' This will take some real review of how hypotheses are handled to implement.
#' Currently the hypothesis is passed to `.post.prob.from.pdfs` and is just used for logic in picking
#' which comparison to make between PDFs, so I could have the dirichlet helper parse the hypothesis
#' string and only return PDFs of the marginal betas that are requested? I think then I could return
#' two for each sample based on the hypothesis, then just always take the first from s1 and the
#' second from s2? That would mean I'd need a helper higher upstream to parse that hypothesis since
#' by the time the `.conj_distHelper` function sees data it only has 1 sample.
#' Not 100% sure that works right now but that feels reasonable.

#' @keywords internal
#' @noRd
.multinomial.pdf.handling <- function(sample_results, hypothesis) {
  g1 <- trimws(gsub("(.*?)([>|<|=|!]{1,2})(.*)", "\\1", hypothesis))
  g2 <- trimws(gsub("(.*?)([>|<|=|!]{1,2})(.*)", "\\3", hypothesis))
  hyp <- trimws(gsub("(.*?)([>|<|=|!]{1,2})(.*)", "\\2", hypothesis))
  hyp <- switch({hyp},
    ">" = "greater", "<" = "lesser", "==" = "equal", "=" = "equal", "unequal"
  )
  if (length(sample_results) == 2) {
    pdf.handling.output <- .post.prob.from.pdfs(
      sample_results[[1]]$pdf[[g1]],
      sample_results[[2]]$pdf[[g2]],
      hyp
    )
  } else {
    pdf.handling.output <- .post.prob.from.pdfs(
      sample_results[[1]]$pdf[[g1]],
      sample_results[[1]]$pdf[[g2]],
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

