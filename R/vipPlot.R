#' Plot Variable Influence on Projection
#'
#' @description This function is used to visualize variable influence on projection (vip) from a plsr
#' model.
#'
#' @param plsrObject Output from pcv.plsr
#' @param i An index from the plsrObject to use if the plsrObject contains models for several outcomes.
#'    Can be a name or a position. Defaults to 1.
#' @param mean Logical, should the mean be plotted (TRUE)
#'      or should the components be shown individually (FALSE, the default).
#' @param removePattern A pattern to remove to make the wavelength column into a numeric.
#'
#' @import ggplot2
#' @return A ggplot showing variable influence on projection
#' @keywords PLSR
#' @examples
#'
#' ## Not run:
#' if (rlang::is_installed("pls"))
#' dists <- list(
#'   rlnorm = list(meanlog = log(40), sdlog = 0.5),
#'   rlnorm = list(meanlog = log(60), sdlog = 0.35)
#' )
#' mv <- mvSim(dists = dists, n_samples = 100, counts = 1000,
#'             min_bin = 1, max_bin = 180, wide = TRUE)
#' sv <- growthSim("logistic",
#'                 n = 5, t = 20,
#'                 params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' d <- cbind(sv, mv[, -1])
#' x <- pcv.plsr(df = d, resps = "y", spectra = grepl("^sim_", colnames(d)))
#' plotVIP(x)
#' }
#'
#' ## End(Not run)
#'
#' @export

plotVIP <- function(plsrObject, i = 1, mean = FALSE, removePattern = ".*_") {
  d <- plsrObject[[i]]$vip_df
  d$spectra <- as.numeric(sub(removePattern, "", d$wavelength))
  if (mean) {
    d$meanVIP <- rowMeans(as.data.frame(d[, colnames(d)[grepl("VIP_[0-9]+?$", colnames(d))]]))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$spectra, y = .data$meanVIP)) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 1, linetype = 5) +
      ggplot2::labs(title = paste0(plsrObject[[i]]$model_performance$outcome), y = "Mean VIP",
                    x = "Spectra")
  } else {
    cols <- c(which(grepl("VIP", colnames(d))), which(colnames(d) == "spectra"))
    d2 <- as.data.frame(data.table::melt(data.table::as.data.table(d[, cols]), id.vars = "spectra",
                                         value.name = "VIP"))
    d2$component <- factor(as.numeric(sub("VIP_", "", d2$variable)))
    p <- ggplot2::ggplot(d2, ggplot2::aes(x = .data$spectra, y = .data$VIP, group = .data$component,
                                          color = .data$component)) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 1, linetype = 5) +
      ggplot2::guides(alpha = "none",
                      color = ggplot2::guide_legend(override.aes = list(linewidth = 3))) +
      ggplot2::labs(title = paste0(plsrObject[[i]]$model_performance$outcome),
                    y = "VIP", x = "Spectra", color = "Component") +
      pcv_theme() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 16))
  }
  return(p)
}
