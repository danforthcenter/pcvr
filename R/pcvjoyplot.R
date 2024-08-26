#' Make Joyplots for multi value trait plantCV data
#'
#' @param df Data frame to use. Long or wide format is accepted.
#' @param index If the data is long then this is a multi value trait as a
#' character string that must be present in `trait`.
#' If the data is wide then this is a string used to find column names to use from the wide data.
#'  In the wide case this should include the entire
#'   trait name (ie, "hue_frequencies" instead of "hue_freq").
#' @param group A length 1 or 2 character vector.
#' This is used for faceting the joyplot and identifying groups for testing.
#' If this is length 1 then no faceting is done.
#' @param y Optionally a variable to use on the y axis. This is useful when you
#' have three variables to display. This argument will change faceting behavior to
#' add an additional layer of faceting (single length group will be faceted,
#' length 2 group will be faceted group1 ~ group2).
#' @param id Optionally a variable to show the outline of different replicates.
#' Note that ggridges::geom_density_ridges_gradient does not support transparency,
#' so if fillx is TRUE then only the outer line will show individual IDs.
#' @param bin Column containing histogram (multi value trait) bins. Defaults to "label".
#' @param freq Column containing histogram counts. Defaults to "value"
#' @param trait Column containing phenotype names. Defaults to "trait".
#' @param fillx Logical, whether or not to use \code{ggridges::geom_density_ridges_gradient}.
#'  Default is T, if F then \code{ggridges::geom_density_ridges} is used instead,
#'   with arbitrary fill. Note that \code{ggridges::geom_density_ridges_gradient}
#'   may issue a message about deprecated ggplot2 features.
#' @keywords multi-value-trait
#' @import ggplot2
#' @import ggridges
#' @import data.table
#' @importFrom stats setNames density aggregate as.formula ks.test
#'
#'
#' @return Returns a ggplot.
#'
#' @examples
#'
#' ## Not run:
#'
#' library(extraDistr)
#' dists <- list(
#'   rmixnorm = list(mean = c(70, 150), sd = c(15, 5), alpha = c(0.3, 0.7)),
#'   rnorm = list(mean = 90, sd = 20),
#'   rlnorm = list(meanlog = log(40), sdlog = 0.5)
#' )
#' x_wide <- mvSim(dists = dists, n_samples = 5, counts = 1000,
#'                 min_bin = 1, max_bin = 180, wide = TRUE)
#' pcv.joyplot(x_wide, index = "sim", group = "group")
#' x_long <- mvSim(dists = dists, n_samples = 5, counts = 1000,
#'                 min_bin = 1, max_bin = 180, wide = FALSE)
#' x_long$trait <- "x"
#' p <- pcv.joyplot(x_long, bin = "variable", group = "group")
#' # we might want to display hues as their hue
#' p + ggplot2::scale_fill_gradientn(colors = scales::hue_pal(l = 65)(360))
#' x_long$group2 <- "example"
#' pcv.joyplot(x_long, bin = "variable", y = "group", fillx = FALSE)
#' ## End(Not run)
#'
#' @export

pcv.joyplot <- function(df = NULL, index = NULL, group = NULL, y = NULL, id = NULL,
                        bin = "label", freq = "value", trait = "trait", fillx = TRUE) {
  #* ***** `general calculated values`

  if (!is.null(trait) && trait %in% colnames(df)) {
    mode <- "long" # if there is a trait column then use long options,
  } else {
    mode <- "wide"
  } # else use wide options

  sub <- .joyPlotFormatData(mode, df, index, trait, bin, freq, group, y, id)

  if (is.null(group)) {
    group <- "dummy"
    df$dummy <- "dummy"
    sub$dummy <- "dummy"
  }

  joyPlotFacetHelperResult <- .joyPlotFacetHelper(y, group, sub)
  facet_layer <- joyPlotFacetHelperResult[["facet"]]
  sub <- joyPlotFacetHelperResult[["data"]]
  sub$grouping <- interaction(sub[, c(y, group)], drop = TRUE)

  #* `if ID is null then aggregate, else draw with ID`
  if (is.null(id)) {
    sub <- stats::aggregate(freq ~ ., data = sub, FUN = mean, na.rm = TRUE)
    gg <- ggplot2::ggplot(sub)
  } else {
    sub$id <- sub[[id]]
    gg <- ggplot2::ggplot(sub, ggplot2::aes(alpha = 0.5, group = interaction(id, y, grouping)))
  }

  ggridgeLayer <- if (fillx) {
    x <- NULL # to make R CMD check happy with stat(x)
    list(
      suppressMessages(ggridges::geom_density_ridges_gradient(
        ggplot2::aes(
          x = .data$bin, y = .data$y,
          height = .data$freq, fill = ggplot2::after_stat(x)
        ),
        show.legend = FALSE, stat = "identity", rel_min_height = 0.001
      )),
      ggplot2::scale_fill_viridis_c(
        option = "plasma"
      )
    )
  } else {
    list(
      suppressMessages(ggridges::geom_density_ridges2(
        ggplot2::aes(
          x = .data$bin, y = .data$y,
          height = .data$freq, fill = .data[[group]], color = .data[[group]]
        ),
        show.legend = FALSE, stat = "identity"
      )),
      ggplot2::scale_color_viridis_d(option = "viridis"),
      ggplot2::scale_fill_viridis_d(option = "viridis")
    )
  }
  p <- gg +
    facet_layer +
    ggridgeLayer +
    ggplot2::scale_x_continuous(n.breaks = 5, labels = ~ round(., 1)) +
    ggplot2::labs(x = index, y = c(y, group)[1]) +
    pcv_theme() +
    ggplot2::theme(legend.position = "none")
  return(p)
}


#' ***********************************************************************************************
#' *************** `format data` ****************************************
#' ***********************************************************************************************
#'
#' @description
#' Internal function for formatting MV trait data
#'
#' @keywords internal
#' @noRd

.joyPlotFacetHelper <- function(y, group, sub) {
  if (!is.null(y)) {
    if (length(group) == 1) {
      sub$y <- sub[[y]]
      facet_layer <- ggplot2::facet_grid(as.formula(paste0("~", group[1])))
    }
    if (length(group) == 2) {
      sub$y <- sub[[y]]
      facet_layer <- ggplot2::facet_grid(as.formula(paste0(group[1], "~", group[2])))
    }
    sub$y <- as.character(sub$y)
  } else { # if y is not provided then one less layer of faceting
    if (length(group) == 1) {
      sub$y <- sub[[group]]
      facet_layer <- list()
    }
    if (length(group) == 2) {
      sub$y <- sub[[group[1]]]
      facet_layer <- ggplot2::facet_grid(as.formula(paste0("~", group[2])))
    }
  }
  return(list("data" = sub, "facet" = facet_layer))
}


#' ***********************************************************************************************
#' *************** `format data` ****************************************
#' ***********************************************************************************************
#'
#' @description
#' Internal function for formatting MV trait data
#'
#' @keywords internal
#' @noRd

.joyPlotFormatData <- function(mode, df, index, trait, bin, freq, group, y, id) {
  #* if long data then subset rows where trait is correct
  if (mode == "long") {
    if (is.null(index)) {
      sub <- df
    } else {
      sub <- df[df[[trait]] == index, ]
    }
    if (length(unique(sub[[trait]])) > 1) {
      warning("More than one trait found, consider an `index` argument")
    }
    sub$bin <- as.numeric(sub[[bin]])
    sub$freq <- as.numeric(sub[[freq]])
    sub <- sub[, c(group, y, id, "bin", "freq", trait)]
  } else if (mode == "wide") { # if wide then get column names that contain index string
    #* subset data to only have index columns
    #* turn the data longer
    sub_wide <- data.table::as.data.table(
      df[, which(colnames(df) %in% c(group, y, id) | grepl(index, colnames(df)))]
    )
    sub <- as.data.frame(data.table::melt(sub_wide,
      id.vars = c(group, y, id),
      variable.name = trait, value.name = freq
    ))
    sub[[bin]] <- sub(index, "", sub[[trait]])
    sub$bin <- as.numeric(regmatches(sub[[bin]], regexpr("[0-9].*", sub[[bin]])))
    sub[[trait]] <- index
    sub$freq <- as.numeric(sub[[freq]])
    sub <- sub[, c(group, y, id, "bin", "freq", trait)]
  }
  return(sub)
}
