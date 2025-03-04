#' Calculating Bayes Factors in conjugate
#'
#' @description
#' Function to calcualte Bayes Factors using single or multi value traits with
#' several distributions in the conjugate function.
#' @param bayes_factor bayes factor range/point hypothesis passed from conjugate
#' @param s_res results from conjugate function thus far, currently the prior
#' and plot_list (for the distribution function name) elements are used. Internally this object is
#' called `sample_results` in conjugate and only has one sample at a time passed to this function.
#' @examples
#' sample_results <- list(
#' "prior" = list("a" = 1, "b" = 1),
#' "plot_list" = list(
#'   "ddist_fun" = "stats::dbeta",
#'   "parameters" = list("shape1" = 144, "shape2" = 96)
#' )
#' )
#' bayes_factor <- NULL
#' .conj_bayes_factor(bayes_factor, sample_results)
#' bayes_factor <- 0.5
#' .conj_bayes_factor(bayes_factor, sample_results)
#' bayes_factor <- c(0.4, 0.6)
#' .conj_bayes_factor(bayes_factor, sample_results)
#'
#' @keywords internal
#' @noRd

.conj_bayes_factor <- function(bayes_factor, s_res) {
  if (is.null(bayes_factor)) {
    return(NULL)
  } else if (length(bayes_factor) == 1) {
    # point hypothesis
    post_args <- append(bayes_factor, s_res$plot_list$parameters)
    names(post_args)[1] <- "x"
    prior_args <- append(bayes_factor, s_res$prior)
    names(prior_args) <- names(post_args)
    ddist <- s_res$plot_list$ddist_fun
    fn_split <- strsplit(ddist, "::")[[1]]
    fn <- get(fn_split[[2]], envir = asNamespace(fn_split[[1]]), mode = "function")
    post_dens <- do.call(fn, post_args)
    prior_dens <- do.call(fn, prior_args)
    bayes_factors <- post_dens / prior_dens
  } else if (length(bayes_factor) == 2) {
    # calculate prior and posterior tail regions
    tail_probs <- lapply(seq_along(bayes_factor), function(i) {
      q <- bayes_factor[i]
      lower_tail <- ifelse(i == 1, TRUE, FALSE)
      post_args <- append(c(q, lower_tail), s_res$plot_list$parameters)
      names(post_args)[1:2] <- c("q", "lower.tail")
      prior_args <- append(c(q, lower_tail), s_res$prior)
      names(prior_args) <- names(post_args)
      # identify functions
      ddist <- s_res$plot_list$ddist_fun
      pdist <- gsub("::d", "::p", ddist)
      fn_split <- strsplit(pdist, "::")[[1]]
      fn <- get(fn_split[[2]], envir = asNamespace(fn_split[[1]]), mode = "function")
      # divide posterior probability into regions
      post_region <- do.call(fn, post_args)
      prior_region <- do.call(fn, prior_args)
      return(
        list("post" = post_region, "prior" = prior_region)
      )
    })
    # calculate prior and posterior odds
    prior_tail_probs <- unlist(lapply(tail_probs, function(l) {return(l$prior)}))
    post_tail_probs <- unlist(lapply(tail_probs, function(l) {return(l$post)}))
    prior_null_prob <- 1 - sum(prior_tail_probs)
    post_null_prob <- 1 - sum(post_tail_probs)
    prior_odds <- sum(prior_tail_probs) / prior_null_prob
    post_odds <- sum(post_tail_probs) / post_null_prob
    # calculate bayes factor
    bayes_factors <- post_odds / prior_odds
  } else {
    stop("bayes_factor must be a point hypothesis or an interval")
  }
  return(bayes_factors)
}
