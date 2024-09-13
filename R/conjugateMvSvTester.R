#' @description
#' Internal function for testing conjugate methods for MV and SV traits, useful for testing only.
#' @param method a conjugate method with mv and sv options
#' @param prior A prior for that conjugate method
#' @param generating a list similar to those for mvSim that will generate data. See example.
#' @examples
#' method = "vonmises"
#' prior = list(mu = 0, kappa = 1, known_kappa = 1, boundary = c(0, 180), n = 1)
#' generating = list(
#'   s1 = list(f = "rnorm", n = 20, mean = 20, sd = 5),
#'   s2 = list(f = "rnorm", n = 20, mean = 170, sd = 5)
#' )
#' .conjugate.mv.sv.testing(method, prior, generating)
#'
#' @keywords internal
#' @noRd

.conjugate.mv.sv.testing <- function(method, prior, generating) {
  s1 <- do.call(generating[[1]][[1]], generating[[1]][-1])
  s2 <- do.call(generating[[2]][[1]], generating[[2]][-1])

  dists1 <- stats::setNames(list(generating[[1]][-1]), generating[[1]][[1]])
  m1 <- mvSim(dists1, n_samples = generating[[1]]$n)[, -1]
  dists2 <- stats::setNames(list(generating[[2]][-1]), generating[[2]][[1]])
  m2 <- mvSim(dists2, n_samples = generating[[2]]$n)[, -1]

  ps <- conjugate(s1, s2, method = method, priors = prior, plot = FALSE, cred.int.level = 0.95)
  pm <- conjugate(m1, m2, method = method, priors = prior, plot = FALSE, cred.int.level = 0.95)
  ps$posterior[[1]]$datatype <- "sv"
  ps$posterior[[2]]$datatype <- "sv"
  pm$posterior[[1]]$datatype <- "mv"
  pm$posterior[[2]]$datatype <- "mv"
  ps$summary$datatype <- "sv"
  pm$summary$datatype <- "mv"
  posteriors <- do.call(rbind, c(ps$posterior, pm$posterior))
  summaries <- rbind(ps$summary, pm$summary)
  data_used <- list(sv_s1 = s1, sv_s2 = s2, mv_s1 = m1, mv_s2 = m2)
  return(list("posteriors" = posteriors, "summaries" = summaries, "data_used" = data_used))
}
