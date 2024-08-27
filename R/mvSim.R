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
#' n_samples = 10
#' counts = 1000
#' min_bin = 0
#' max_bin = 180
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
#'
#' \dontrun{
#'   ggplot(data = data.frame(x = 1:180), aes(x = x)) +
#'     facet_grid(dist ~ group) +
#'     lapply(c("a", "b"), function(grp) {
#'       lapply(1:length(dists), function(di) {
#'         rfun <- names(dists)[di]
#'         rargs <- append(dists[[di]], 1000)
#'         names(rargs)[length(rargs)] <- "n"
#'         v1 <- do.call(rfun, rargs)
#'         v1[v1 > max_bin] <- max_bin
#'         v1[v1 < min_bin] <- min_bin
#'         d <- data.frame(group = grp, x = v1, dist = di)
#'         geom_histogram(data = d, aes(x = x, fill = di),
#'         alpha = 1, binwidth = 1, show.legend = FALSE)
#'       })
#'     }) +
#'     pcv_theme()
#' }
#'
#'
#'
#' ## End(Not run)
#' @export

mvSim <- function(dists = list(rnorm = list(mean = 100, sd = 15)),
                  n_samples = 10, counts = 1000, min_bin = 1, max_bin = 180, wide = TRUE,
                  binwidth = 1) {
  if (length(n_samples) == 1) {
    n_samples <- rep(n_samples, length(dists))
  }
  vecs <- .makeVecs(dists, counts, n_samples)
  out <- .simFreqs(vecs, max_bin, min_bin, binwidth)
  if (!wide) {
    out$id <- seq_len(nrow(out))
    out <- as.data.frame(data.table::melt(data.table::as.data.table(out), id.vars = c("group", "id")))
  }
  return(out)
}

#' internal vector making helper function
#' @keywords internal
#' @noRd

.makeVecs <- function(dists, counts, n_samples) {
  funs <- names(dists)
  out <- lapply(seq_along(funs), function(i) {
    f <- funs[i]
    fun <- match.fun(f)
    dists[[i]]$n <- counts
    lapply(1:n_samples[i], function(e) {
      do.call(fun, args = dists[[i]])
    })
  })
  names(out) <- names(dists)
  return(out)
}

#' internal histogram making helper function
#' @keywords internal
#' @noRd

.simFreqs <- function(vecs, max_bin, min_bin, binwidth) {
  do.call(rbind, lapply(seq_along(vecs), function(i) {
    vecName <- names(vecs)[i]
    vec <- vecs[[i]]
    do.call(rbind, lapply(vec, function(v) {
      if (vecName == "rbeta") {
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
}
