#' Show nls/lm model in ggplot layer
#'
#' @description
#' Add predicted mean trendline from a model fit with \link{growthSS} and \link{fitGrowth}
#' to a ggplot object.
#'
#' @param mapping Set of aesthetic mappings created by \code{ggplot2::aes()}.
#' If specified and ‘inherit.aes = TRUE’ (the default),
#' it is combined with the default mapping at the top level of the plot.
#' If there is no mapping then it is filled in by default using the \code{ss} object.
#' @param data The data to be displayed in this layer. This behaves per normal ggplot2 expectations
#' except that if data is missing (ie, not inherited or specified) then the data from \code{ss} is used.
#' @param fit A brmsfit object, typically returned from \code{fitGrowth}.
#' @param ss A \code{pcvrss} object. Only the "pcvrForm" and "df" elements are used.
#' @param inherit.aes Logical, should aesthetics be inherited from top level? Defaults to TRUE.
#' @param ... Additional arguments passed to ggplot2::layer.
#'
#' @import ggplot2
#' @rdname stat_growthss
#' @keywords ggplot
#' @export

stat_nls_model <- function(mapping = NULL, data = NULL,
                           fit = NULL, ss = NULL,
                           inherit.aes = TRUE, ...) {
  # These would normally be arguments to a stat layer but they should not be changed
  geom <- "line"
  position <- "identity"
  na.rm <- FALSE
  show.legend <- c("color" = TRUE)
  parsed_form <- .parsePcvrForm(ss$pcvrForm, ss$df)
  # get elements to replace NULL defaults in case they are missing
  if (is.null(data) || is.null(mapping)) {
    data <- data %||% parsed_form$data
    mapping <- mapping %||% ggplot2::aes(x = .data[[parsed_form$x]])
  }
  lyr <- ggplot2::layer(
    stat = statNlsMod, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, fit = fit, parsed_form = parsed_form, ...)
  )
  return(lyr)
}

statNlsMod <- ggplot2::ggproto("StatNls", Stat,
  # `specify that there will be extra params`
  extra_params = c("na.rm", "fit", "parsed_form"),
  # `data setup function`
  setup_data = function(data, params) {
    #' leaving the same as statBrmsMod for now.
    if ("data" %in% names(params$parsed_form)) {
      parsed_form <- params$parsed_form
      mod_data <- parsed_form$data
      mod_data <- mod_data[, unlist(parsed_form[c("x", "group")])]
      colnames(mod_data) <- c("x", "MOD_GROUP")
      mod_data <- mod_data[!duplicated(mod_data), ]
      data <- plyr::join(mod_data, data, type = "left", match = "all", by = "x")
      if (length(unique(data$PANEL)) > 1 && parsed_form$USEG) {
        data <- data[data$PANEL == as.numeric(as.factor(data$MOD_GROUP)), ]
      }
      if (!parsed_form$USEG) {
        data$MOD_GROUP <- ""
      }
    }
    return(data)
  },
  #' NOTE I think this can stay the same from brms version
  compute_panel = function(self, data, scales, fit, parsed_form, probs, ...) {
    if (ggplot2:::empty(data)) return(ggplot2:::data_frame0())
    groups <- split(data, data[["MOD_GROUP"]], drop = TRUE)
    stats <- lapply(seq_along(groups), function(i) {
      d <- self$compute_group(
        data = groups[[i]], scales = scales, fit = fit, parsed_form = parsed_form, ...
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
    dropped <- non_constant_columns[!non_constant_columns %in% c(
      self$dropped_aes, unlist(parsed_form[c("individual", "x")])
    )]

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
  #' `make data out of model per a given aes-group, should only be 1 per panel`
  #' this is the heavily customized component which makes data for ribbons from
  #' the model and ss objects.
  compute_group = function(data, scales, fit = NULL, parsed_form = NULL, probs = NULL, ...) {
    yvar <- parsed_form$y
    xvar <- parsed_form$x
    group <- parsed_form$group
    # make data to use drawing posterior predictions
    nd <- data[, c("x", "MOD_GROUP", "PANEL")]
    nd <- nd[!duplicated(nd), ]
    colnames(nd) <- c(xvar, group, "PANEL")
    #* `add predictions`
    preds <- data.frame(pred = stats::predict(fit, newdata = nd))
    keep <- which(!duplicated(preds$pred))
    plotdf <- nd[keep, ]
    plotdf$pred <- preds[keep, "pred"]
    grpdf <- plotdf[, c(xvar, group, "PANEL", "pred")]
    colnames(grpdf) <- c("x", "model_group", "PANEL", "pred")
    return(grpdf)
  },
  # set defaults for several aesthetics, all have to be after stat is calculated
  default_aes = aes(
    y = after_stat(pred),
    color = after_stat(model_group),
    group = after_stat(model_group),
    x = after_stat(x)
  )
)


# the only difference in gam vs nls in growthPlot is about gam's not accepting
# null newdata, that is handled already with the stat here so it is fine to use
# the same function.
stat_gam_model <- stat_nls_model
