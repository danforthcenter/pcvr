#' Multi Value Trait simulating function
#'
#' @description mvSim can be used to simulate data for example models/plots.
#'
#' @param n_samples Number of samples per distribution to generate. Defaults to 10, can be >1L.
#' @param counts Number of counts per histogram, defaults to 1000.
#' @param max_bin The number of bins to return. Note that this is also the max value that will be
#' accepted in the distribution functions, with higher numbers being shrunk to this value.
#' Defaults to 180.
#' @param min_bin The minumum bin number. This can be thought of as the minimum value that will
#' be accepted in the distribution functions, with lower numbers being raised to this value.
#' Note that bin arguments are both ignored in the case of "rbeta" and treated as 0,1.
#' @param dists A list of lists, with names corresponding to random deviate generating functions
#' and arguments to the function in the list values (see examples). Note that the n argument
#' does not need to be provided.
#' @param wide Boolean, should data be returned in wide format (the default)?
#' If FALSE then long data is returned.
#' @param binwidth How wide should bins be? Defaults to 1.
#' @param t Number of timepoints to simulate. Defaults to NULL in which case data is simulated as
#' non-longitudinal. Note that currently the first non \code{n} argument of the data simulating
#' function is assumed to be the parameter changing over time (ie, mean in rnorm, meanlog in rlnorm).
#' @param model A type of growth model, passed to \link{growthSim}. This is only used if t is
#' specified.
#' @param params Parameters for the growth model, passed to \link{growthSim}. This is also only used
#' if t is specified. Note growth will start from the values specified in dists. See examples.
#' @keywords multi-value
#' @return Returns a dataframe of example multi-value trait data simulated from specified distributions.
#'
#' @importFrom graphics hist
#' @importFrom data.table melt as.data.table
#'
#' @examples
#'
#' library(extraDistr) # for rmixnorm
#' library(ggplot2)
#' dists <- list(
#'   rmixnorm = list(mean = c(70, 150), sd = c(15, 5), alpha = c(0.3, 0.7)),
#'   rnorm = list(mean = 90, sd = 3)
#' )
#' x <- mvSim(dists = dists, wide = FALSE)
#' dim(x)
#' x2 <- mvSim(dists = dists)
#' dim(x2)
#'
#' ggplot(x, aes(
#'   x = as.numeric(sub("sim_", "", variable)),
#'   y = value, group = interaction(group, id), fill = group
#' )) +
#'   geom_col(position = "identity", alpha = 0.25) +
#'   pcv_theme() +
#'   labs(x = "bin")
#' dists = list(rnorm = list(mean = 30, sd = 15), rnorm = list(mean = 25, sd = 10))
#' x3 <- mvSim(
#'   dists = dists, wide = FALSE, # here we make longitudinal data
#'   t = 10, model = "linear", params = list("A" = c(10, 5))
#' )
#' ggplot(x3, aes(
#'   x = as.numeric(sub("sim_", "", variable)),
#'   y = value, group = interaction(group, id), fill = group
#' )) +
#'   facet_wrap(~times) +
#'   geom_col(position = "identity", alpha = 0.25) +
#'   pcv_theme() +
#'   labs(x = "bin")
#'
#' @export

mvSim <- function(dists = list(rnorm = list(mean = 100, sd = 15)),
                  n_samples = 10, counts = 1000, min_bin = 1, max_bin = 180, wide = TRUE,
                  binwidth = 1, t = NULL, model = "linear", params = list("A" = 10)) {
  if (length(n_samples) == 1) {
    n_samples <- rep(n_samples, length(dists))
  }
  if (!is.null(t)) {
    n_samples <- rep(n_samples, t)
  }
  names(dists) <- paste(names(dists), seq_along(dists), sep = "_")
  res <- .longitudinalDists(dists, t, model, params)
  dists2 <- res$dist
  times <- res$times
  vecs <- .makeVecs(dists2, counts, n_samples)
  out <- .simFreqs(vecs, max_bin, min_bin, binwidth, times)
  if (!wide) {
    out <- cbind(
      data.frame(id = seq_len(nrow(out))),
      out
    )
    out <- as.data.frame(
      data.table::melt(
        data.table::as.data.table(out),
        id.vars = seq_len(
          min(
            which(
              grepl("sim_", colnames(out))
            )
          ) - 1
        )
      )
    )
  }
  return(out)
}

#' internal vector making helper function
#' @keywords internal
#' @noRd

.makeVecs <- function(dists, counts, n_samples) {
  funs <- sub("_[0-9]+$", "", names(dists))
  out <- lapply(seq_along(funs), function(i) {
    f <- funs[i]
    fun <- match.fun(f)
    dists[[i]]$n <- counts
    lapply(seq_len(n_samples[i]), function(e) {
      do.call(fun, args = dists[[i]])
    })
  })
  names(out) <- names(dists)
  return(out)
}

#' internal histogram making helper function
#' @keywords internal
#' @noRd

.simFreqs <- function(vecs, max_bin, min_bin, binwidth, times) {
  out <- do.call(rbind, lapply(seq_along(vecs), function(i) {
    vecName <- names(vecs)[i]
    vec <- vecs[[i]]
    do.call(rbind, lapply(vec, function(v) {
      if (grepl("rbeta", vecName)) {
        v[v > 1] <- 1
        v[v < 0] <- 0
        s1 <- hist(v, breaks = seq(0, 1, length.out = 100), plot = FALSE)$counts
        s1d <- as.data.frame(cbind(data.frame(vecName), matrix(s1, nrow = 1)))
        colnames(s1d) <- c("group", paste0("sim_", seq(0.01, 0.99, 0.01)))
      } else {
        v[v > max_bin] <- max_bin
        v[v < min_bin] <- min_bin
        s1 <- hist(v, breaks = seq(min_bin, (max_bin + binwidth), binwidth), plot = FALSE)$counts
        s1d <- as.data.frame(cbind(data.frame(vecName), matrix(s1, nrow = 1)))
        colnames(s1d) <- c("group", paste0("sim_", seq(min_bin, max_bin, binwidth)))
      }
      s1d
    }))
  }))
  if (any(as.logical(times))) {
    out <- cbind(
      data.frame(times = rep(times, each = unique(unlist(lapply(vecs, length)))[1])),
      out
    )
  }
  return(out)
}

#' interal function to replicate dists for longitudinal data
#' @keywords internal
#' @noRd

.longitudinalDists <- function(dists, t, model, params) {
  times <- 0
  if (!is.null(t)) {
    pred <- growthSim(model = model, n = 1, t = t, params = params)[, c("group", "y")]
    pred$group <- sapply(pred$group, function(i) {
      which(letters == i)
    })
    times <- rep(seq_len(t), times = length(dists))
    nms <- sub("_[0-9]+$", "", names(dists))
    dists <- Reduce(append, lapply(seq_along(nms), function(i) {
      nm <- nms[i]
      args <- names(formals(nm))
      args <- args[-which(args == "n")]
      vary <- args[1]
      temporal_means <- dists[[i]][[vary]] + pred[pred$group == i, "y"]
      tdists <- lapply(temporal_means, function(tm) {
        diter <- dists[[i]]
        diter[[vary]] <- tm
        diter
      })
      names(tdists) <- rep(names(dists)[i], t)
      tdists
    }))
  }
  return(list("dist" = dists, "times" = times))
}
