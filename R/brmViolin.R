#' Function to visualize hypotheses tested on brms models similar to those made using growthSS outputs.
#'
#' @param fit A brmsfit object or a dataframe of draws. If you need to combine multiple models then
#' use \link{combineDraws} to merge their draws into a single dataframe for testing.
#' @param ss A \code{pcvrss} object. The only component that is currently used is the pcvrForm.
#' @param hypothesis A hypothesis expressed as a character string in the style of that used by
#' \code{brms::hypothesis} and \link{testGrowth}. Note that currently this only supports hypotheses
#' using two parameters from the model at a time (ie, "groupA / groupB > 1.05" works but
#' "(groupA / groupB) - (groupC / groupD) > 1" does not). In the hypothesis "..." can be used to mean
#' "all groups for this parameter" so that the hypothesis "... / A_group1 > 1.05" would include all
#' the "A" coefficients for groups 1:N along the x axis, see examples.
#'
#' @keywords brms Bayesian pcvrss
#'
#' @import ggplot2
#' @importFrom data.table melt
#' @importFrom viridis plasma
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' simdf <- growthSim(
#'   "logistic",
#'   n = 20, t = 25,
#'   params = list("A" = c(200, 180, 190, 160), "B" = c(13, 11, 10, 10), "C" = c(3, 3, 3.25, 3.5))
#' )
#' ss <- growthSS(
#'   model = "logistic", form = y ~ time | id / group, sigma = "int",
#'   list("A" = 130, "B" = 10, "C" = 3),
#'   df = simdf, type = "brms"
#' )
#'
#' fit <- fitGrowth(ss, backend = "cmdstanr", iter = 500, chains = 1, cores = 1)
#' brmViolin(fit, ss, ".../A_groupd > 1.05") # all groups used
#' brmViolin(fit, ss, "A_groupa/A_groupd > 1.05") # only these two groups
#' }
#'
#' @return Returns a ggplot showing a brms model's posterior distributions
#' as violins and filled by posterior probability of some hypothesis.
#'
#' @export

brmViolin <- function(fit, ss, hypothesis) {
  fitdf <- as.data.frame(fit)
  if (grepl("\\.\\.\\.", hypothesis)) {
    ldj <- .brmViolinEllipsis(fitdf, ss, hypothesis)
  } else {
    ldj <- .brmViolinFixedHypothesis(fitdf, ss, hypothesis)
  }
  virPal <- unlist(lapply(c(1, 0.9, 0.75, 0.5, 0.25), function(i) viridis::plasma(1, 1, i)))
  violinPlot <- ggplot2::ggplot(ldj, ggplot2::aes(
    x = .data[["variable"]], y = .data[["draw"]],
    fill = .data[["discrete_post_prob"]]
  )) +
    ggplot2::geom_violin(show.legend = TRUE) +
    ggplot2::geom_hline(
      yintercept = mean(ldj[ldj$ref, "draw"]),
      linetype = 5, linewidth = 0.5
    ) +
    ggplot2::scale_fill_manual(
      values = virPal,
      labels = c(">99%", ">95%", ">85%", ">75%", "<75%"), drop = FALSE
    ) +
    ggplot2::labs(
      y = "Posterior Distribution",
      x = hypothesis,
      fill = "Discrete Posterior Probability"
    ) +
    pcv_theme() +
    ggplot2::theme(
      legend.position = "bottom", axis.text.x.bottom = ggplot2::element_text(angle = 0, hjust = 0.5),
      panel.border = ggplot2::element_rect(fill = NA)
    )
  return(violinPlot)
}

#' @noRd
#' @keywords internal

.brmViolinEllipsis <- function(fitdf, ss, hypothesis) {
  math_or_underscore <- "\\+|\\/|\\-|\\*|\\^|>|<|=|_"
  hsplit <- trimws(strsplit(hypothesis, math_or_underscore)[[1]])
  hsplit <- hsplit[as.logical(nchar(hsplit))]
  hsplit <- hsplit[!grepl("\\.\\.\\.", hsplit)]
  hsplit <- hsplit[suppressWarnings(is.na(as.numeric(hsplit)))]
  grouping <- .parsePcvrForm(ss$pcvrForm)$group
  which_grouping <- which(unlist(lapply(grouping, function(grp) {
    any(grepl(grp, hsplit))
  })))
  grouping <- grouping[which_grouping]
  par <- hsplit[which(!grepl(grouping, hsplit))]
  ref_group <- hsplit[which(grepl(grouping, hsplit))]
  colnames(fitdf) <- sub("^b_", "", colnames(fitdf))
  draws <- fitdf[, grepl(paste0(par, ".*", grouping), colnames(fitdf))]
  if (!grepl("[:]", grouping)) {
    draws <- draws[, !grepl("[:]", colnames(draws))] # drop interaction terms
  }
  if (length(par) > 0) {
    colnames(draws) <- sub(paste0(par, "_"), "", colnames(draws))
    colnames(draws) <- paste0(par, "_", colnames(draws))
  }
  hyps_df <- do.call(rbind, lapply(colnames(draws), function(grp) {
    iter_hyp <- gsub("\\.{3,}", grp, hypothesis)
    x <- as.data.frame(brms::hypothesis(draws, iter_hyp)$h)
    x$ref_group <- ifelse(as.logical(length(par)), paste0(par, "_", ref_group), ref_group)
    x$comparison <- grp
    x[, c("Post.Prob", "ref_group", "comparison")]
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
  longdraw$group <- sub(paste0(par, "_"), "", longdraw$variable)
  longdraw$ref <- grepl(ref_group, longdraw$variable)
  ldj <- merge(longdraw, hyps_df, by.x = c("variable"), by.y = c("comparison"))
  return(ldj)
}


#' @noRd
#' @keywords internal

.brmViolinFixedHypothesis <- function(fitdf, ss, hypothesis) {
  math <- "\\+|\\/|\\-|\\*|\\^|>|<|=" # foreseeable math operators
  hsplit <- trimws(strsplit(hypothesis, math)[[1]])
  hsplit <- hsplit[suppressWarnings(is.na(as.numeric(hsplit)))]
  draws <- fitdf
  colnames(draws) <- sub("^b_", "", colnames(draws))
  draws <- draws[, hsplit]
  hyps_df <- as.data.frame(brms::hypothesis(draws, hypothesis = hypothesis)$h)
  hyps_df$ref_group <- hsplit[2]
  hyps_df$comparison <- hsplit[1]
  hyps_df <- rbind(
    hyps_df[, c("Post.Prob", "ref_group", "comparison")],
    data.frame(Post.Prob = 0.5, ref_group = hsplit[2], comparison = hsplit[2])
  )
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
  longdraw$ref <- grepl(hsplit[2], longdraw$variable)
  ldj <- merge(longdraw, hyps_df, by.x = c("variable"), by.y = c("comparison"))
  return(ldj)
}
