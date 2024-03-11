#' Default theme for ggplots made by pcvr functions.
#'
#' @import ggplot2
#' @importFrom ggplot2 %+replace%
#' @export
#'
pcv_theme <- function() {
  ggplot2::theme_minimal() %+replace%
    ggplot2::theme(
      axis.text.x.bottom = ggplot2::element_text(hjust = 1),
      axis.line.y.left = ggplot2::element_line(),
      axis.line.x.bottom = ggplot2::element_line(),
      strip.background = ggplot2::element_rect(fill = "gray50", color = "gray20"),
      strip.text.x = ggplot2::element_text(size = 14, color = "white"),
      strip.text.y = ggplot2::element_text(size = 14, color = "white")
    )
}
