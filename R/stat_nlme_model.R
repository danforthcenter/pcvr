#' @import ggplot2
#' @rdname stat_growthss
#' @keywords ggplot
#' @export

stat_nlme_model <- function(mapping = NULL, data = NULL,
                            fit = NULL, ss = NULL,
                            inherit.aes = TRUE, ...) {
  # These would normally be arguments to a stat layer but they should not be changed
  geom <- "smooth"
  position <- "identity"
  na.rm <- FALSE
  show.legend <- c("color" = TRUE)
  # get elements to replace NULL defaults in case they are missing
  if (is.null(data) || is.null(mapping)) {
    parsed_form <- .parsePcvrForm(ss$pcvrForm, ss$df)
    data <- data %||% parsed_form$data
    mapping <- mapping %||% ggplot2::aes(x = .data[[parsed_form$x]])
    show.legend <- c("color" = parsed_form$USEG)
  }
  # make layer for each of the intervals
  lyr <- ggplot2::layer(
    stat = statNlmeMod, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, fit = fit, parsed_form = parsed_form, se = TRUE, ...)
  )
  return(lyr)
}

#' @export
#' @rdname stat_growthss
stat_lme_model <- function(mapping = NULL, data = NULL,
                           fit = NULL, ss = NULL,
                           inherit.aes = TRUE, ...) {
  lyr <- stat_nlme_model(mapping, data, fit, ss, inherit.aes, ...)
  return(lyr)
}

statNlmeMod <- ggplot2::ggproto("StatNlme", Stat,
  # `specify that there will be extra params`
  extra_params = c("na.rm", "fit", "parsed_form"),
  # `data setup function`
  setup_data = function(data, params) {
    if ("data" %in% names(params$parsed_form)) {
      parsed_form <- params$parsed_form
      mod_data <- parsed_form$data
      mod_data <- mod_data[, unlist(parsed_form[c("x", "group")])]

      message(paste(
        paste(colnames(mod_data), collapse = ", "), "\n",
        nrow(mod_data), "\n\n", paste(
          do.call(
            cbind,
            lapply(mod_data, \(x) paste(unique(x), collapse = ", "))
          ),
          collapse = " -- "
        )
      ))

      colnames(mod_data) <- c("x", "MOD_GROUP")
      mod_data <- mod_data[!duplicated(mod_data), ]
      data <- plyr::join(mod_data, data, type = "left", match = "all", by = "x")
      if (length(unique(data$PANEL)) > 1) {
        data <- data[data$PANEL == as.numeric(as.factor(data$MOD_GROUP)), ]
      }
      if (!parsed_form$USEG) {
        data$MOD_GROUP <- parsed_form$group
      }
      # add predictions
      fit <- params$fit
      xvar <- parsed_form$x
      group <- parsed_form$group

      nd <- data[, c("x", "MOD_GROUP", "PANEL")]
      nd <- nd[!duplicated(nd), ]
      colnames(nd) <- c(xvar, group, "PANEL")
      preds <- data.frame(pred = round(as.numeric(stats::predict(fit, newdata = nd)), 4))
      message("made it past predict")
      keep <- which(!duplicated(preds$pred))
      data_with_preds <- nd[keep, ]
      data_with_preds$trendline <- preds[keep, "pred"]
      message("adding sigma boundaries")
      data_with_preds <- .add_sigma_bounds(data_with_preds, fit, xvar, group)
      message("made sigma bounds")
      data <- data_with_preds[, c(xvar, group, "PANEL", "trendline", "sigma_ymin", "sigma_ymax")]
      colnames(data) <- c("x", "MOD_GROUP", "PANEL", "pred", "ymin", "ymax")
    }
    return(data)
  },
  compute_panel = function(self, data, scales, ...) {
    if (ggplot2:::empty(data)) return(ggplot2:::data_frame0())
    groups <- split(data, data[["MOD_GROUP"]])
    stats <- lapply(groups, function(groupdf) {
      d <- self$compute_group(
        data = groupdf, scales = scales, ...
      )
      return(d)
    })
    non_constant_columns <- character(0)
    stats <- mapply(function(new, old) {
      if (ggplot2:::empty(new)) {
        return(ggplot2:::data_frame0())
      }
      old <- old[, !(names(old) %in% names(new)), drop = FALSE]
      non_constant <- vapply(old, vctrs::vec_unique_count, integer(1)) > 1L
      non_constant_columns <<- c(non_constant_columns, names(old)[non_constant])
      vc <- vctrs:::vec_cbind(
        new,
        old[rep(1, nrow(new)), , drop = FALSE]
      )
      return(vc)
    }, stats, groups, SIMPLIFY = FALSE)

    non_constant_columns <- ggplot2:::unique0(non_constant_columns)
    dropped <- non_constant_columns[!non_constant_columns %in% self$dropped_aes]

    if (length(dropped) > 0) {
      warning(paste0(
        "The ", paste(dropped, collapse = ", "), " aesthetics were dropped,\n",
        " did you forget to specify a group aesthetic or convert a numerical variable into a factor?"
      ))
    }
    data_new <- ggplot2:::vec_rbind0(!!!stats)
    return(
      data_new[, !names(data_new) %in% non_constant_columns, drop = FALSE]
    )
  },
  compute_group = function(data, scales, ...) {
    return(data)
  },
  # set defaults for several aesthetics, all have to be after stat is calculated
  default_aes = aes(
    y = after_stat(pred),
    ymin = after_stat(ymin),
    ymax = after_stat(ymax),
    color = after_stat(MOD_GROUP),
    group = after_stat(MOD_GROUP),
    x = after_stat(x)
  )
)
