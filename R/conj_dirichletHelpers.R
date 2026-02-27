#' @description
#' Internal function for calculating the dirichlet distribution underlying
#' multinomially distributed single value traits (counts).
#' @param s1 A named vector/list/data.frame of numerics drawn from a multinomial distribution.
#' @examples
#' out <- .conj_dirichlet_sv(
#'   s1 = list("A" = 10, "B" = 10, "C" = 5),
#'   cred.int.level = 0.95,
#'   plot = TRUE
#' )
#' lapply(out, head)
#' @details
#' See Examples 1.4, 1.6, and 1.7 for thoughts on default dirichlet prior here
#' https://arxiv.org/pdf/1504.02689
#' @keywords internal
#' @noRd
.conj_multinomial_sv <- function(s1 = NULL, priors = NULL,
                               support = NULL, cred.int.level = NULL,
                               calculatingSupport = FALSE) {
  #' `make default prior if none provided`
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

  #' `Define support if it is missing`
  #' support is the same as beta binomial

  #' `Make Posterior draws`
  #' Posterior draws used for ROPE comparisons mainly.
  #' need to consider if it's better to draw from rdirichlet or rbeta

  #' `calculate HDI`
  #' note that this might change some since we'd have HDE/HDI for every
  #' marginal beta

  
  #' `calculate HDE`
  #' note that this might change some since we'd have HDE/HDI for every
  #' marginal beta

  
  #' `Store summary`
  #' note that this might change some since we'd have HDE/HDI for every
  #' marginal beta

  #' `Save s1 data for plotting`
  #' note that this might change significantly from other distributions
  #' because we want to test/plot the marginal betas to keep the dimensions
  #' low enough to understand (and make the math tractable?)
}
