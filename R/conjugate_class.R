#' Class \code{conjugate} for output from the \code{pcvr::conjugate} function.
#'
#' Comparisons made by the \code{conjugate} function return objects of this class
#' containing parameters of the prior and posterior distributions, hypothesis tests,
#' ROPE tests, Bayes Factors, and plots of the posterior.
#'
#' @name conjugate-class
#' @docType class
#'
#' @details
#' See \code{methods(class = "conjugate")} for an overview of available methods.
#'
#' @slot summary Summary data frame of results
#' @slot posterior Posterior distribution as a list of named lists
#' @slot prior Prior distribution as a list of named lists
#' @slot plot Optionally a plot of the distributions and their differences
#' @slot data The data from s1 and s2 arguments to \link{conjugate}.
#' @slot call Matched call to \link{conjugate}.
#'
#' @seealso
#'   \code{\link{conjugate}}
#'
NULL

as.conjugate <- function(x) {
  class(x) <- "conjugate"
  return(x)
}

#' Print a \code{conjugate} object.
#'
#' @aliases print.conjugate
#'
#' @param x An object of class \code{conjugate}.
#' @param ... further arguments, passed to print.default.
#'
#' @seealso \code{\link{summary.conjugate}}
#' @method print conjugate
#' @export
print.conjugate <- function(x, ...) {
  return(print(summary.conjugate(x), ...))
}


#' Summarize a \code{conjugate} object.
#'
#' @aliases summary.conjugate
#'
#' @param object An object of class \code{conjugate}.
#' @param ... further arguments, passed to print.default.
#'
#' @method summary conjugate
#' @export

summary.conjugate <- function(object, ...) {
  out <- object[which(names(object) %in% c("summary", "posterior", "prior", "call"))]
  class(out) <- "conjugatesummary"
  return(out)
}

#' Print a \code{conjugatesummary} object.
#'
#' @aliases print.conjugatesummary
#'
#' @param x An object of class \code{conjugatesummary}.
#' @param ... further arguments, which are currently ignored.
#'
#' @seealso \code{\link{print.conjugatesummary}}
#' @method print conjugatesummary
#' @export

print.conjugatesummary <- function(x, ...) {
  call_list <- as.list(x$call)
  method <- eval(call_list$method)
  if (length(method) == 1) {
    #* `Method`
    method_list <- switch(method[1],
      t = list("Normal", "Mu", "T"),
      gaussian = list("Normal", "Mu", "Normal"),
      beta = list("Beta", "Mean", "Beta"),
      binomial = list("Beta", "Rate", "Binomial"),
      lognormal = list("Normal", "Mu", "Lognormal"),
      lognormal2 = list("Gamma", "Precision", "Lognormal"),
      poisson = list("Gamma", "Lambda", "Poisson"),
      negbin = list("Beta", "Rate", "Negative Binomial"),
      vonmises = list("Von Mises", "Direction", "Von Mises"),
      vonmises2 = list("Von Mises", "Direction", "Von Mises"),
      uniform = list("Pareto", "Upper Boundary", "Uniform"),
      pareto = list("Gamma", "Scale", "Pareto"),
      gamma = list("Gamma", "Rate", "Gamma"),
      bernoulli = list("Beta", "Rate", "Bernoulli"),
      exponential = list("Gamma", "Rate", "Exponential"),
      bivariate_uniform = list("Bivariate Pareto", "Boundaries", "Uniform"),
      bivariate_gaussian = list("Normal/Gamma", "Mu/Sd", "Normal"),
      bivariate_lognormal = list("Normal/Gamma", "Mu/Sd", "Lognormal")
    )
    method_statement <- paste0(
      method_list[[1]], # conjugate distribution
      " distributed ",
      method_list[[2]], # conjugate parameter
      " parameter of ",
      method_list[[3]], # data distribution
      " distributed data."
    )
    cat(method_statement)
    cat("\n\n")
    #* `Parameters`
    lapply(seq_along(x$posterior), function(i) {
      prior <- x$prior[[i]]
      post <- x$posterior[[i]]
      prior_statement <- paste0(
        "Sample ", i, " Prior ", method_list[[1]], "(",
        paste(
          paste(names(prior), round(unlist(prior), 3), sep = " = "),
          collapse = ", "
        ),
        ")\n"
      )
      posterior_statement <- paste0(
        "\tPosterior ", method_list[[1]], "(",
        paste(
          paste(names(post), round(unlist(post), 3), sep = " = "),
          collapse = ", "
        ),
        ")\n"
      )
      cat(prior_statement)
      cat(posterior_statement)
      return(invisible(list(prior_statement, posterior_statement)))
    })
  }
  #* `Hypothesis`
  if ("post.prob" %in% colnames(x$summary)) {
    cat("\n")
    hyp_statement <- paste0(
      "Posterior probability that S1",
      switch(x$summary$hyp[1],
        unequal = " is not equal to ",
        equal = " is equal to ",
        lesser = " is less than ",
        greater = " is greater than "
      ),
      "S2 = ",
      100 * round(x$summary$post.prob, 5), "%"
    )
    cat(hyp_statement)
    cat("\n")
  }
  #* `ROPE`
  if ("rope_prob" %in% colnames(x$summary)) {
    cat("\n")
    rope_ci <- eval(call_list$rope_ci)
    rope_range <- eval(call_list$rope_range)
    rope_message <- paste0(
      "Probability of the difference between ",
      method_list[[2]],
      " parameters being within [",
      paste(rope_range, collapse = ":"),
      "] using a ",
      100 * rope_ci, "% Credible Interval is ",
      100 * round(x$summary$rope_prob[1], 5),
      "% with an average difference of ",
      round(x$summary$HDE_rope, 3)
    )
    cat(rope_message)
    cat("\n\n")
  }
  #* `Bayes Factors`
  if ("bf_1" %in% colnames(x$summary)) {
    bayes_factor <- eval(call_list$bayes_factor)
    lapply(seq_along(x$posterior), function(i) {
      bf_statement <- paste0(
        "Sample ", i, " Bayes Factor ",
        ifelse(length(bayes_factor) == 1, "at ", "in "),
        paste(bayes_factor, collapse = " to "),
        " = ", round(x$summary[[paste0("bf_", i)]][1], 3), "\n"
      )
      cat(bf_statement)
      return(invisible(bf_statement))
    })
  }
  cat("\n")
  print(x$summary)
  return(invisible(x))
}
