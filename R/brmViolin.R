#' Function to visualize hypotheses tested on brms models similar to those made using growthSS outputs.
#'
#'
#'
#' @param fit A brmsfit object or a list of brmsfit objects
#' @param params A list of parameters to use from the fit.
#' Defaults to NULL in which case all growth model parameters are used.
#' @param hyp A character string defining the hypothesis to be tested.
#' Defaults to "num/denom > 1.05". The "num" and "denom" names should be kept,
#' but the direction and magnitude can be changed.
#' @param compareX Which groups in the model should be compared as numerator options?
#' Defaults to NULL in which case a plot will not be made but the data will be returned.
#' @param againstY Which group in the model should be used as the
#' denominator (typically a control group to compare against)?
#' Defaults to NULL in which case a plot will not be made but the data will be returned.
#' @param group_sep A regex pattern to match the separation of grouping
#' terms in the models group term (see the formula argument of \code{\link{growthSS}}).
#' The default uses "[.]" to break on a single period.
#' @param groups_into A vector of column names to make after groups are split by group_sep.
#' If this or groups_sep are NULL then no groups are assumed.
#' @param x The variable to be plotted on the x axis (should be from groups_into).
#' @param facet The variable to be used to facet the
#' ggplot (should be another option from groups_into).
#' If left NULL then the plot will only be faceted by params.
#' Note that with the nature of againstY this faceting is often redundant
#' but it does add labels which are helpful for keeping results organized..
#' @param cores Optional number of cores to run hypotheses in parallel.
#' Defaults to 1 unless the "mc.cores" option is set.
#' @param returnData Logical, should data be returned?
#' This is treated as TRUE if a plot will not be generated but
#' otherwise defaults to FALSE.
#'
#' @keywords brms, ggplot2
#'
#' @import ggplot2
#' @import viridis
#' @import parallel
#' @importFrom utils combn
#' @importFrom stats as.formula setNames
#'
#' @examples
#'
#' ## Not run:
#'
#' if (FALSE) {
#' set.seed(123)
#' simdf <- growthSim(
#'   "logistic", n = 20, t = 25,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group, sigma = "spline",
#'   list("A" = 130, "B" = 10, "C" = 3),
#'   df = simdf, type = "brms"
#' )
#'
#' fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
#' brmViolin(fit, hyp = "num/denom>1.05",
#'           compareX = "a",
#'           againstY = "b", returnData = TRUE)
#' }
#' ## End(Not run)
#'
#' @return Returns a ggplot showing a brms model's posterior distributions
#' as violins and filled by posterior probability of some hypothesis.
#'
#' @export

brmViolin <- function(fit, params = NULL, hyp = "num/denom>1.05", compareX = "a", againstY = "b",
                      group_sep = "[.]", groups_into = c("group2"), x = NULL, facet = NULL,
                      cores = getOption("mc.cores", 1), returnData = FALSE) {
  #* `parse arguments`
  model <- fit
  brmViolinParseArgumentsRes <- .brmViolinParseArguments(
    compareX, againstY,
    model, params, group_sep,
    groups_into
  )
  compareFew <- brmViolinParseArgumentsRes[["compareFew"]]
  compareX <- brmViolinParseArgumentsRes[["compareX"]]
  model <- brmViolinParseArgumentsRes[["model"]]
  params <- brmViolinParseArgumentsRes[["params"]]
  useGroups <- brmViolinParseArgumentsRes[["useGroups"]]
  if (is.null(x)) {
    x <- groups_into[1]
  }
  #* `get draws`
  draws <- do.call(cbind, lapply(model, function(mod) { # extract draws for relevant parameters
    mdf <- as.data.frame(mod[["fit"]])
    colPattern <- paste0("^b_[", paste0(params, collapse = "|"), "]")
    mdf <- mdf[, grepl(colPattern, colnames(mdf))]
    colnames(mdf) <- sub("^b_", "", colnames(mdf))
    mdf
  }))
  #* take "grouping" string from formula (pform) and using that with the parameter names I can find
  #* groups for any models?
  #* For models fit with only one group of many from data this would fail, I'll need another option
  #* for that.
  #* probably fine to build this out first then think about how that edge case works.
  group_string <- trimws(
    strsplit(
      as.character(model[[1]]$formula$pforms[[params[1]]])[3], "[+]"
    )[[1]]
  )[2]

  groupings <- unique(sub(paste0(".*", group_string), "", colnames(draws)))
  p1 <- combn(groupings, 2, simplify = FALSE)
  p2 <- lapply(p1, rev)
  p3 <- lapply(unique(groupings), function(g) c(g, g))
  comparisons <- c(p1, p2, p3)
  if (compareFew) {
    comparisons <- comparisons[unlist(lapply(comparisons, function(comp) {
      comp[1] %in% compareX & comp[2] == againstY
    }))]
  }

  colnames(draws) <- sub(group_string, "", colnames(draws))

  hyps_df <- do.call(rbind, lapply(params, function(param) {
    param_df <- do.call(rbind, parallel::mclapply(comparisons, function(comp) {
      num <- paste(param, comp[1], sep = "_")
      denom <- paste(param, comp[2], sep = "_")
      temp <- draws
      temp$num <- temp[[num]]
      temp$denom <- temp[[denom]]
      x <- as.data.frame(brms::hypothesis(temp, paste0(hyp))$h)
      x$param <- param
      x$num <- sub("grouping", "", comp[1])
      x$denom <- sub("grouping", "", comp[2])
      x[, c("Post.Prob", "param", "num", "denom")]
    }, mc.cores = cores))
    param_df
  }))

  hyps_df$discrete_post_prob <- factor(
    ifelse(hyps_df$Post.Prob >= 0.99, "A",
      ifelse(hyps_df$Post.Prob >= 0.95, "B",
        ifelse(hyps_df$Post.Prob >= 0.85, "C",
          ifelse(hyps_df$Post.Prob >= 0.75, "D", "E")
        )
      )
    ),
    levels = c("A", "B", "C", "D", "E"), ordered = TRUE
  )

  longdraw <- as.data.frame(data.table::melt(data.table::as.data.table(draws),
    measure.vars = colnames(draws), value.name = "draw"
  ))
  longdraw$param <- substr(longdraw$variable, 1, 1)
  longdraw$group <- sub(paste0("[", paste0(params, collapse = "|"), "]_"), "", longdraw$variable)
  if (useGroups) {
    group_meta <- do.call(rbind, parallel::mclapply(longdraw$group,
      function(g) {
        setNames(
          data.frame(t(strsplit(g, group_sep)[[1]])),
          groups_into
        )
      },
      mc.cores = cores
    ))
    for (col in groups_into) {
      if (suppressWarnings(any(is.na(as.numeric(group_meta[[col]]))))) {
        group_meta[[col]] <- factor(group_meta[[col]])
      } else {
        group_meta[[col]] <- factor(group_meta[[col]],
          levels = sort(as.numeric(unique(group_meta[[col]]))),
          ordered = TRUE
        )
      }
    }

    longdraw <- cbind(longdraw, group_meta)
  }
  ldj <- merge(longdraw, hyps_df, by.x = c("group", "param"), by.y = c("num", "param"))

  if (!compareFew) {
    return(ldj)
  }

  virPal <- unlist(lapply(c(1, 0.9, 0.75, 0.5, 0.25), function(i) viridis::plasma(1, 1, i)))

  # still need to separate "grouping" into x and y variables.

  if (is.null(facet)) {
    facet_layer <- ggplot2::facet_wrap(~param, scales = "free_y")
  } else {
    facet_layer <- ggplot2::facet_grid(as.formula(paste0("param~", facet)), scales = "free_y")
  }

  violinPlot <- ggplot2::ggplot(ldj, ggplot2::aes(
    x = .data[[x]], y = .data[["draw"]],
    fill = .data[["discrete_post_prob"]]
  )) +
    facet_layer +
    ggplot2::geom_violin() +
    lapply(unique(ldj$param), function(p) {
      x <- data.frame(param = p, mean = mean(ldj[ldj$param == p & ldj$group == ldj$denom, "draw"]))
      ggplot2::geom_hline(data = x, ggplot2::aes(yintercept = mean), linetype = 5, linewidth = 0.5)
    }) +
    ggplot2::scale_fill_manual(
      values = virPal, breaks = c("A", "B", "C", "D", "E"),
      labels = c(">99%", ">95%", ">85%", ">75%", "<75%"), drop = TRUE
    ) +
    ggplot2::labs(y = "Posterior Distribution", x = x, fill = "Discrete Posterior Probability") +
    pcv_theme() +
    ggplot2::theme(
      legend.position = "bottom", axis.text.x.bottom = ggplot2::element_text(angle = 0, hjust = 0.5),
      panel.border = ggplot2::element_rect(fill = NA)
    )
  if (returnData) {
    return(list("plot" = violinPlot, "data" = ldj))
  } else {
    return(violinPlot)
  }
}


.brmViolinParseArguments <- function(compareX, againstY, model, params, group_sep, groups_into) {
  compareFew <- FALSE
  if (!is.null(compareX) && !is.null(againstY)) {
    compareFew <- TRUE
    if (!againstY %in% compareX) {
      compareX <- c(compareX, againstY)
    } # make sure all violins will be filled if both are provided
  }

  if (any(unlist(lapply(model, class)) != "brmsfit")) {
    # if only one brmsfit is given then make a list of itself
    model <- list(model)
  }

  if (is.null(params)) { # if params aren't given then grab all
    nlPars <- names(model[[1]]$formula$pforms)
    params <- nlPars[-which(grepl("sigma|nu", nlPars))]
  }
  useGroups <- TRUE
  if (is.null(group_sep) || is.null(groups_into)) {
    useGroups <- FALSE
  }
  return(list(
    "compareFew" = compareFew,
    "compareX" = compareX,
    "model" = model,
    "params" = params,
    "useGroups" = useGroups
  ))
}
