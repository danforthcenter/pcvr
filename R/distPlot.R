#' Function for plotting iterations of posterior distributions
#'
#' @param fits A list of brmsfit objects following the same data over time.
#'    Currently checkpointing is not supported.
#' @param form A formula describing the growth model similar to \code{\link{growthSS}}
#'    and \code{\link{brmPlot}} such as: outcome ~ predictor |individual/group
#' @param df data used to fit models (this is used to plot each subject's trend line).
#' @param priors a named list of samples from the prior distributions for each parameter in
#'     \code{params}. This is only used if sample_prior=FALSE in the brmsfit object.
#'     If left NULL then no prior is included.
#' @param params a vector of parameters to include distribution plots of.
#'     Defaults to NULL which will use all parameters from the top level model. Note that these
#'     parameters have to be estimated per each group in the model, if you have interecept only terms
#'     (estimated once across all groups) then manually specify params to not include those.
#' @param maxTime Optional parameter to designate a max time not observed in the models so far
#' @param patch Logical, should a patchwork plot be returned or should lists of ggplots be returned?
#' @param virOptions A vector of names or letters for which viridis maps to use for each group.
#' @keywords Bayesian brms
#' @import ggplot2
#' @import patchwork
#' @importFrom methods is
#' @importFrom stats setNames
#' @import viridis
#' @return A ggplot or a list of ggplots (depending on patch).
#' @export
#' @examples
#' \donttest{
#' f <- "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/brmsFits.rdata"
#' tryCatch(
#'   {
#'     print(load(url(f)))
#'     library(brms)
#'     library(ggplot2)
#'     library(patchwork)
#'     fits <- list(fit_3, fit_15)
#'     form <- y~time | id / group
#'     priors <- list(
#'       "phi1" = rlnorm(2000, log(130), 0.25),
#'       "phi2" = rlnorm(2000, log(12), 0.25),
#'       "phi3" = rlnorm(2000, log(3), 0.25)
#'     )
#'     params <- c("A", "B", "C")
#'     d <- simdf
#'     maxTime <- NULL
#'     patch <- TRUE
#'     from3to25 <- list(
#'       fit_3, fit_5, fit_7, fit_9, fit_11,
#'       fit_13, fit_15, fit_17, fit_19, fit_21, fit_23, fit_25
#'     )
#'     distributionPlot(
#'       fits = from3to25, form = y ~ time | id / group,
#'       params = params, d = d, priors = priors, patch = FALSE
#'     )
#'     distributionPlot(
#'       fits = from3to25, form = y ~ time | id / group,
#'       params = params, d = d, patch = FALSE
#'     )
#'   },
#'   error = function(e) {
#'     message(e)
#'   }
#' )
#' }
#'

distributionPlot <- function(fits, form, df, priors = NULL,
                             params = NULL, maxTime = NULL, patch = TRUE,
                             virOptions = c("plasma", "mako", "viridis", "cividis",
                                            "magma", "turbo", "inferno", "rocket")) {
  #* ***** `Reused helper variables`
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  group <- parsed_form$group
  d <- parsed_form$data
  fitData <- fits[[length(fits)]]$data
  dSplit <- split(d, d[[group]])
  startTime <- min(unlist(lapply(fits, function(ft) {
    return(min(ft$data[[x]], na.rm = TRUE))
  })))
  if (is.null(maxTime)) {
    maxTime <- max(unlist(lapply(fits, function(ft) {
      return(max(ft$data[[x]], na.rm = TRUE))
    })))
  }
  palettes <- lapply(
    seq_along(unique(fitData[[group]])),
    function(i) {
      group_pal <- viridis::viridis(length(fits),
        begin = 0.1,
        end = 1, option = virOptions[i], direction = 1
      )
      return(group_pal)
    }
  )
  names(palettes) <- unique(fitData[[group]])

  #* ***** `if params is null then pull them from growth formula`

  if (is.null(params)) {
    fit <- fits[[1]]
    growthForm <- as.character(fit$formula[[1]])[[3]]
    # note, simplify to proper regex
    test <- gsub(x, "", growthForm)
    test2 <- gsub("exp\\(", "", test)
    test3 <- gsub("\\(1", "", test2)
    test4 <- gsub("[/]|[+]|[-]|[)]|[()]|[*]|\\^", "", test3)
    params <- strsplit(test4, "\\s+")[[1]]
  }

  #* ***** `growth trendline plots`

  growthTrendPlots <- lapply(seq_along(dSplit), function(i) {
    dt <- dSplit[[i]]
    p <- ggplot2::ggplot(dt, ggplot2::aes(
      x = .data[[x]], y = .data[[y]], color = .data[[x]],
      group = .data[[individual]]
    )) +
      ggplot2::geom_line(show.legend = FALSE) +
      viridis::scale_color_viridis(begin = 0.1, end = 1, option = virOptions[i], direction = 1) +
      ggplot2::scale_x_continuous(limits = c(startTime, maxTime)) +
      pcv_theme()
    return(p)
  })
  growth_lims <- do.call(rbind, lapply(growthTrendPlots, function(p) {
    x_lims <- ggplot2::layer_scales(p)$x$range$range
    y_lims <- ggplot2::layer_scales(p)$y$range$range
    return(matrix(c(x_lims, y_lims), nrow = 1))
  }))
  x_min <- min(growth_lims[, 1])
  x_max <- max(growth_lims[, 2])
  y_min <- min(growth_lims[, 3])
  y_max <- max(growth_lims[, 4])
  growthTrendPlots <- lapply(growthTrendPlots, function(p) {
    p_update <- p +
      ggplot2::coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max))
    return(p_update)
  })

  #* ***** `posterior distribution extraction`

  posts <- do.call(rbind, lapply(fits, function(fit) {
    time <- max(fit$data[[x]], na.rm = TRUE)
    fitDraws <- do.call(cbind, lapply(params, function(par) {
      draws <- as.data.frame(fit)[grepl(paste0("b_", par), colnames(as.data.frame(fit)))]
      if (nrow(brms::prior_draws(fit)) > 1) {
        draws <- draws[!grepl("^prior_", colnames(draws))]
      }
      colnames(draws) <- gsub("^b_", "", colnames(draws))
      splits <- strsplit(colnames(draws), split = "")
      mx <- max(unlist(lapply(splits, length)))
      ind <- which(unlist(lapply(seq_len(mx), function(i) {
        l_over_1 <- length(unique(rapply(splits, function(j) {
          return(j[i])
        }))) != 1
        return(l_over_1)
      })))
      if (length(ind) > 0) {
        colnames(draws) <- paste(par, unlist(lapply(colnames(draws), function(c) {
          return(substr(c, min(ind), max(c(ind, nchar(c)))))
        })), sep = "_")
      }
      return(draws)
    }))
    fitDraws[[x]] <- time
    return(fitDraws)
  }))

  #* ***** `prior distribution extraction`
  distPlotPriorExtractionRes <- .distPlotPriorExtraction(fits, priors, d, group, params, x)
  prior_df <- distPlotPriorExtractionRes[["prior_df"]]
  USEPRIOR <- distPlotPriorExtractionRes[["UP"]]

  #* ***** `posterior distribution plots`

  if (USEPRIOR) {
    posts <- rbind(prior_df, posts)
  }
  posts[[x]] <- factor(posts[[x]], levels = sort(as.numeric(unique(posts[[x]]))), ordered = TRUE)

  xlims <- lapply(params, function(par) {
    diff <- as.numeric(as.matrix(posts[, grepl(paste0("^", par, "_"), colnames(posts))]))
    rng <- c(min(diff, na.rm = TRUE), max(diff, na.rm = TRUE))
    return(rng)
  })
  names(xlims) <- params
  postPlots <- lapply(unique(fitData[[group]]), function(groupVal) {
    groupPlots <- lapply(params, function(par) {
      p <- ggplot2::ggplot(posts) +
        ggplot2::geom_density(ggplot2::aes(
          x = .data[[paste(par, groupVal, sep = "_")]],
          fill = .data[[x]], color = .data[[x]],
          group = .data[[x]]
        ), alpha = 0.8) +
        ggplot2::labs(x = paste(par, group, groupVal)) +
        ggplot2::coord_cartesian(xlim = xlims[[par]]) +
        pcv_theme() +
        ggplot2::theme(
          axis.text.x.bottom = ggplot2::element_text(angle = 0),
          legend.position = "none", axis.title.y = ggplot2::element_blank()
        )

      if (USEPRIOR) {
        p <- p + ggplot2::scale_fill_manual(values = c("black", palettes[[groupVal]])) +
          ggplot2::scale_color_manual(values = c("black", palettes[[groupVal]]))
      } else {
        p <- p + ggplot2::scale_fill_manual(values = palettes[[groupVal]]) +
          ggplot2::scale_color_manual(values = palettes[[groupVal]])
      }
      return(p)
    })
    return(groupPlots)
  })

  if (patch) {
    ncol_patch <- 1 + length(params)
    nrow_patch <- length(postPlots)

    patchPlot <- growthTrendPlots[[1]] + postPlots[[1]]
    if (length(unique(d[[group]])) > 1) {
      for (i in 2:length(growthTrendPlots)) {
        patchPlot <- patchPlot + growthTrendPlots[[i]] + postPlots[[i]]
      }
    }
    out <- patchPlot + patchwork::plot_layout(ncol = ncol_patch, nrow = nrow_patch)
  } else {
    out <- list(growthTrendPlots, postPlots)
  }
  return(out)
}

#' Prior extraction in distPlot
#' @keywords internal
#' @noRd

.distPlotPriorExtraction <- function(fits, priors, d, group, params, x) {
  if (is.null(priors)) {
    return(list("prior_df" = NULL, "UP" = FALSE))
  }
  if (all(unlist(lapply(fits, function(fit) nrow(brms::prior_draws(fit)) < 1)))) {
    # if no models were fit with sample_prior
    USEPRIOR <- FALSE
    if (!is.null(priors)) { # if prior is supplied as argument
      USEPRIOR <- TRUE
      if (!methods::is(priors[[1]], "list")) {
        priors <- lapply(seq_along(unique(d[[group]])), function(i) priors)
        names(priors) <- unique(d[[group]])
      }
      prior_df <- do.call(cbind, lapply(names(priors), function(nm) {
        nmp <- priors[[nm]]
        nm_res <- setNames(data.frame(do.call(cbind, lapply(names(nmp), function(nmpn) {
          return(nmp[[nmpn]])
        }))), paste0(names(nmp), "_", nm))
        return(nm_res)
      }))
      prior_df[[x]] <- 0
    }
  } else { #* `need to fit some models with sample_prior and see how this works with them`
    prior_df <- brms::prior_draws(fits[[1]])
    prior_df <- prior_df[, grepl(paste0("b_", paste0(params, collapse = "|")), colnames(prior_df))]
    colnames(prior_df) <- gsub(group, "", colnames(prior_df))
    colnames(prior_df) <- gsub("^b_", "", colnames(prior_df))
    prior_df[[x]] <- 0
    USEPRIOR <- TRUE
  }
  return(list("prior_df" = prior_df, "UP" = USEPRIOR))
}
