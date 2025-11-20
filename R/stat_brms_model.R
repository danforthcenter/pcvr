#' Show brms model in ggplot layer
#'
#' @description
#' Add posterior predictive distribution from a brms model fit with \link{growthSS} and \{fitGrowth}
#' to a ggplot object.
#'
#' @param mapping Set of aesthetic mappings created by \code{ggplot2::aes()}. If specified
#' and ‘inherit.aes = TRUE’ (the default), it is combined with the default mapping at the top level of the plot.
#' If there is no mapping then it is filled in by default using the \code{ss} object.
#' @param data The data to be displayed in this layer. This behaves per normal ggplot2 expectations except
#' that if data is missing (ie, not inherited or specified) then the data from \code{ss} is used.
#' @param fit A brmsfit object, typically returned from \code{fitGrowth}.
#' @param ss A \code{pcvrss} object. Only the "pcvrForm" and "df" elements are used.
#' @param probs A vector of Credible Intervals to plot, defaults to 0.95.
#' @param inherit.aes Logical, should aesthetics be inherited from top level? Defaults to TRUE.
#' @param ... Additional arguments passed to ggplot2::layer.
#'
#' @import ggplot2
#' @import vctrs
#' @importFrom plyr join
#' @importFrom cli cli_warn
#' @importFrom rlang try_fetch, inject
#' @details
#' @examples
#'
#' @rdname stat_growthSS
#' @seealso \link{growthPlot} for a self-contained plotting function
#' @keywords ggplot
#' @export

stat_brms_model <- function(mapping = NULL, data = NULL,
                          fit = NULL, ss = NULL, CI = 0.95,
                          inherit.aes = TRUE, ...) {
  # These would normally be arguments to a stat layer but they should not be changed
  geom = "ribbon"
  position = "identity"
  na.rm = FALSE
  show.legend = c("fill" = TRUE, "alpha" = FALSE) 
  # get elements to replace NULL defaults in case they are missing
  if (is.null(data) || is.null(mapping)) {
    parsed_form <- .parsePcvrForm(ss$pcvrForm, ss$df)
    data <- data %||% parsed_form$data
    mapping <- mapping %||% ggplot2::aes(x = .data[[parsed_form$x]])
  }
  # format credible intervals into a list of c(min, max) probs for predictions
  formatted_prob_list <- lapply(rev(sort(cis)), function(i) {
    c(((1 - i) / 2), (i + (1 - i) / 2))
  })
  # make layer for each of the intervals
  lapply(formatted_prob_list, function(prob_pair) {
    ggplot2::layer(
      stat = statBrmsMod, data = data, mapping = mapping, geom = geom,
      position = position, show.legend = show.legend, inherit.aes = inherit.aes,
      params = list(na.rm = na.rm, fit = fit, parsed_form = parsed_form, probs = prob_pair, ...)
    )
  })
}


"%||%" <- function(a, b) {
  if (!is.null(a)) {a} else {b}
}

#' ~I don't know the export convention but it seems like this is "soft hidden"~
statBrmsMod <- ggplot2::ggproto("StatBrm", Stat,
  # `specify that there will be extra params`
  extra_params = c("na.rm", "fit", "parsed_form", "probs"),
  # `data setup function`
  setup_data = function(data, params){
    #' possible that ss is not a pcvrss object for compatibility with other brms models
    #' if "df" is part of it then work with that otherwise use general data.
    if ("data" %in% names(params$parsed_form)) {
      parsed_form <- params$parsed_form
      mod_data <- parsed_form$data
      mod_data <- mod_data[, unlist(parsed_form[c("x", "group")])]
      colnames(mod_data) <- c("x", "MOD_GROUP")
      mod_data <- mod_data[!duplicated(mod_data), ]
      data <- plyr::join(mod_data, data, type = "left", match = "all", by = "x")
      if (length(unique(data$PANEL)) > 1) {
        data <- data[data$PANEL == as.numeric(as.factor(data$MOD_GROUP)), ]
      }
    }
    return(data)
  },
  #' NOTE ggplot2:::Stat$compute_layer can use the default from ggproto
  #' `make plot within a given panel of the ggplot (a facet)`
  #' this is mostly the same as the default ggproto compute_panel function,
  #' but it takes more named args and passes them to compute_group and
  #' avoids warning about individual/time columns.
  compute_panel = function(self, data, scales, fit, parsed_form, probs, ...){
    if (ggplot2:::empty(data)) return(ggplot2:::data_frame0())
    groups <- split(data, data[["MOD_GROUP"]])
    stats <- lapply(groups, function(groupdf) {
      self$compute_group(data = groupdf, scales = scales,
        fit = fit, parsed_form = parsed_form, probs = probs,
        ...)
    })
    non_constant_columns <- character(0)
    stats <- mapply(function(new, old) {
      if (ggplot2:::empty(new)) {return(ggplot2:::data_frame0())}
      old <- old[, !(names(old) %in% names(new)), drop = FALSE]
      non_constant <- vapply(old, vctrs::vec_unique_count, integer(1)) > 1L
      non_constant_columns <<- c(non_constant_columns, names(old)[non_constant])
      vctrs:::vec_cbind(
        new,
        old[rep(1, nrow(new)), , drop = FALSE]
      )
    }, stats, groups, SIMPLIFY = FALSE)

    non_constant_columns <- ggplot2:::unique0(non_constant_columns)
    dropped <- non_constant_columns[!non_constant_columns %in% c(
      self$dropped_aes, unlist(parsed_form[c("individual", "x")])
    )]

    if (length(dropped) > 0) {
      cli::cli_warn(c(
        "The following aesthetics were dropped during statistical transformation: {.field {dropped}}.",
        "i" = "This can happen when ggplot fails to infer the correct grouping structure in the data.",
        "i" = "Did you forget to specify a {.code group} aesthetic or to convert a numerical variable into a factor?"
      ))
    }
    data_new <- ggplot2:::vec_rbind0(!!!stats)
    data_new[, !names(data_new) %in% non_constant_columns, drop = FALSE]
  },
  #' `make data out of model per a given aes-group, should only be 1 per panel`
  #' this is the heavily customized component which makes data for ribbons from
  #' the model and ss objects.
  compute_group = function(data, scales,
                           fit = NULL, parsed_form = NULL, probs = NULL,
                           ...){
    yvar <- parsed_form$y
    xvar <- parsed_form$x
    group <- parsed_form$group
    # make data to use drawing posterior predictions
    nd <- data[, c("x", "MOD_GROUP", "PANEL")]
    nd <- nd[!duplicated(nd), ]
    colnames(nd) <- c(xvar, group, "PANEL")
    # make predictions
    mod_data <- cbind(nd, predict(fit, newdata = nd, probs = probs))
    # lengthen predictions as in brmPlot
    longPreds <- do.call(rbind, lapply(seq_len(nrow(mod_data)), function(r) {
      sub <- mod_data[r, ]
      lp <- do.call(rbind, lapply(
        head(probs, floor(length(probs) / 2)) * 100,
        function(i) {
        min <- paste0("Q", i)
        max <- paste0("Q", 100 - i)
        iter <- sub[, c(xvar, group, "Estimate", "PANEL")]
        iter$q <- abs(((100 - i) - i))
        iter$min <- sub[[min]]
        iter$max <- sub[[max]]
        return(iter)
      }))
      return(lp)
    }))
    # select columns and rename
    grpdf <- longPreds[, c(xvar, group, "Estimate", "PANEL", "q", "min", "max")]
    colnames(grpdf) <- c("x", "MOD_GROUP", "y", "PANEL", "Cred.Int", "ymin", "ymax")
    return(grpdf)
  },
  # set defaults for several aesthetics, all have to be after stat is calculated
  default_aes = aes(
    ymin = after_stat(ymin),
    ymax = after_stat(ymax),
    fill = after_stat(Cred.Int),
    alpha = after_stat(0.5),
    group = after_stat(MOD_GROUP),
    x = after_stat(x)
  )
)
