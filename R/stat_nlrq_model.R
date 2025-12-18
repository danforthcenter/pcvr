#' @import ggplot2
#' @rdname stat_growthss
#' @keywords ggplot
#' @export

stat_nlrq_model <- function(mapping = NULL, data = NULL,
                            fit = NULL, ss = NULL,
                            inherit.aes = TRUE, ...) {
  # These would normally be arguments to a stat layer but they should not be changed
  geom <- "line"
  position <- "identity"
  na.rm <- FALSE
  show.legend <- c("color" = TRUE)
  parsed_form <- .parsePcvrForm(ss$pcvrForm, ss$df)
  # Multiple quantiles make fit a list of nlrq fits, standardize format
  if (methods::is(fit, "nlrq")) {
    fit <- list(fit)
    names(fit) <- fit[[1]]$m$tau()
  } # after here `fit` will always be a list of nlrq models.
  # make layer for each of the taus
  # those layers can just use statNlsMod.
  lyrs <- lapply(fit, function(one_fit) {
    tau <- one_fit$m$tau()
    # get elements to replace NULL defaults in case they are missing
    # if color is not otherwise specified then use tau from the fit
    if (is.null(data) || is.null(mapping)) {
      data <- data %||% parsed_form$data
      mapping <- mapping %||% ggplot2::aes(x = .data[[parsed_form$x]], color = tau)
    }
    lyr <- ggplot2::layer(
      stat = statNlsMod, data = data, mapping = mapping, geom = geom,
      position = position, show.legend = show.legend, inherit.aes = inherit.aes,
      params = list(na.rm = na.rm, fit = one_fit, parsed_form = parsed_form, ...)
    )
    return(lyr)
  })
  return(lyrs)
}
