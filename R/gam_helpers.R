#' Helper function for visualizing differences in GAMs fit with \code{mgcv::gam}
#'
#' Note that using GAMs will be less useful than fitting parameterized models as supported by
#' \code{growthSS} and \code{fitGrowth} for common applications in plant phenotyping.
#'
#' @param model A model fit with smooth terms by \code{mgcv::gam}
#' @param newdata A data.frame of new data to use to make predictions. If this is left NULL
#' (the default) then
#' an attempt is made to make newdata using the first smooth term in the formula.
#' See examples for guidance on making appropriate newdata
#' @param g1 A character string for the level of byVar to use as the first group to compare,
#' if plot=TRUE then this will be shown in blue.
#' @param g2 The second group to compare (comparison will be g1 - g2). If plot=TRUE then this will be
#' shown in red.
#' @param byVar Categorical variable name used to separate splines as a string.
#' @param smoothVar The variable that splines were used on. This will often be a time variable.
#' @param cis Confidence interval levels, can be multiple. For example, 0.95 would return Q_0.025 and
#' Q_0.975 columns, and c(0.9, 0.95) would return Q_0.025, Q_0.05, Q_0.95, and Q_0.975 columns.
#' Defaults to \code{seq(0.05, 0.95, 0.05)}
#' @param unconditional Logical, should unconditional variance-covariance be used in calculating
#' standard errors. Defaults to TRUE.
#' @param plot Logical, should a plot of the difference be returned? Defaults to TRUE.
#'
#' @keywords gam
#' @importFrom mgcv gam s
#' @importFrom stats vcov predict df.residual qt
#' @importFrom viridis viridis
#' @import ggplot2
#' @import patchwork
#' @return A dataframe or a list containing a ggplot and a dataframe
#' @examples
#'
#' ex <- pcvr::growthSim("logistic",
#'   n = 20, t = 25,
#'   params = list(
#'     "A" = c(200, 160),
#'     "B" = c(13, 11),
#'     "C" = c(3, 3.5)
#'   )
#' )
#'
#' m <- mgcv::gam(y ~ group + s(time, by = factor(group)), data = ex)
#'
#' support <- expand.grid(
#'   time = seq(min(ex$time), max(ex$time), length = 400),
#'   group = factor(unique(ex$group))
#' )
#'
#' out <- gam_diff(
#'   model = m, newdata = support, g1 = "a", g2 = "b",
#'   byVar = "group", smoothVar = "time", plot = TRUE
#' )
#' dim(out$data)
#' out$plot
#' out2 <- gam_diff(
#'   model = m, g1 = "a", g2 = "b", byVar = NULL, smoothVar = NULL, plot = TRUE
#' )
#' @export

gam_diff <- function(model, newdata = NULL, g1, g2, byVar = NULL, smoothVar = NULL,
                     cis = seq(0.05, 0.95, 0.05), unconditional = TRUE, plot = TRUE) {
  form <- model$formula
  rhs <- as.character(form)[3]
  rg <- regexpr("s\\([a-zA-Z0-9.]+", rhs)
  xTerm <- sub("s\\(", "", regmatches(rhs, rg)[1])
  rg2 <- regexpr("by\\s?=\\s?.*[a-zA-Z0-9.]+", rhs)
  byTerm <- sub(".*\\(", "", sub("by\\s?=\\s?", "", regmatches(rhs, rg2)[1]))
  mdf <- model$model

  if (is.null(newdata)) {
    newdata <- expand.grid(
      x = seq(min(mdf[[xTerm]]), max(mdf[[xTerm]]),
        length.out = round(diff(range(mdf[[xTerm]])) * 500)
      ),
      g = factor(unique(mdf[[byTerm]]))
    )
    colnames(newdata) <- c(xTerm, byTerm)
  }
  if (is.null(byVar)) {
    byVar <- byTerm
  }
  if (is.null(smoothVar)) {
    smoothVar <- xTerm
  }

  xp <- stats::predict(model, newdata = newdata, type = "lpmatrix")
  c1 <- grepl(g1, colnames(xp))
  c2 <- grepl(g2, colnames(xp))
  r1 <- newdata[[byVar]] == g1
  r2 <- newdata[[byVar]] == g2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, !(c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl("^s\\(", colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% stats::vcov(model, unconditional = unconditional)) * X))
  df.resid <- stats::df.residual(model)

  cis_diff_df <- do.call(cbind, lapply(cis, function(ci) {
    crit <- stats::qt(1 - ((1 - ci) / 2), df.resid, lower.tail = TRUE)
    upr <- dif + (crit * se)
    lwr <- dif - (crit * se)
    setNames(data.frame(upr, lwr), paste0(
      c("Q_diff_", "Q_diff_"),
      round(c(1 - ((1 - ci) / 2), 1 - (1 - ((1 - ci) / 2))), 3)
    ))
  }))
  cis_diff_df <- cis_diff_df[, order(colnames(cis_diff_df))]

  out_df <- cbind(data.frame(g1 = g1, g2 = g2, mu = dif, se = se), cis_diff_df)
  if (plot) {
    #* `g1 model CIs`
    x_g1 <- xp[r1, ]
    ## keep only relevant splines and intercept term
    c1_g1 <- grepl(g1, colnames(xp)) | grepl("\\(Intercept\\)", colnames(xp))
    x_g1[, !(c1_g1)] <- 0
    ## zero out the parametric cols
    x_g1[, !grepl("\\(Intercept\\)|^s\\(", colnames(xp))] <- 0
    g1_vals <- x_g1 %*% coef(model)
    se_g1 <- sqrt(rowSums((x_g1 %*% stats::vcov(model, unconditional = unconditional)) * x_g1))
    cis_g1_df <- do.call(cbind, lapply(cis, function(ci) {
      crit <- stats::qt(1 - ((1 - ci) / 2), df.resid, lower.tail = TRUE)
      upr <- g1_vals + (crit * se_g1)
      lwr <- g1_vals - (crit * se_g1)
      setNames(data.frame(upr, lwr), paste0(
        c("Q_g1_", "Q_g1_"),
        round(c(1 - ((1 - ci) / 2), 1 - (1 - ((1 - ci) / 2))), 3)
      ))
    }))
    cis_g1_df <- cis_g1_df[, order(colnames(cis_g1_df))]
    #* `g2 model CIs`
    x_g2 <- xp[r2, ]
    ## keep only relevant splines and intercept term
    c2_g2 <- grepl(g2, colnames(xp)) | grepl("\\(Intercept\\)", colnames(xp))
    x_g2[, !(c2_g2)] <- 0
    ## zero out the parametric cols
    x_g2[, !grepl("\\(Intercept\\)|^s\\(", colnames(xp))] <- 0
    g2_vals <- x_g2 %*% coef(model)
    se_g2 <- sqrt(rowSums((x_g2 %*% stats::vcov(model, unconditional = unconditional)) * x_g2))
    cis_g2_df <- do.call(cbind, lapply(cis, function(ci) {
      crit <- stats::qt(1 - ((1 - ci) / 2), df.resid, lower.tail = TRUE)
      upr <- g2_vals + (crit * se_g2)
      lwr <- g2_vals - (crit * se_g2)
      setNames(data.frame(upr, lwr), paste0(
        c("Q_g2_", "Q_g2_"),
        round(c(1 - ((1 - ci) / 2), 1 - (1 - ((1 - ci) / 2))), 3)
      ))
    }))
    cis_g2_df <- cis_g2_df[, order(colnames(cis_g2_df))]

    out_df <- cbind(out_df, cis_g1_df, cis_g2_df)
  }

  out_df$df.resid <- df.resid

  smoothVarRange <- range(newdata[[smoothVar]], na.rm = TRUE)
  smoothVarOut <- seq(min(smoothVarRange), max(smoothVarRange), length.out = length(dif))
  out_df[[smoothVar]] <- smoothVarOut
  out <- out_df
  if (plot) {
    p_diff <- .plot_gam_diff(out_df, name_pattern = "Q_diff_") +
      ggplot2::theme(
        axis.title.x.bottom = ggplot2::element_blank(),
        axis.text.x.bottom = ggplot2::element_blank()
      )
    p_model <- .plot_gam_diff(out_df,
      name_pattern = "Q_g1_",
      name_pattern2 = "Q_g2_"
    )
    layout_obj <- patchwork::plot_layout(
      design = c(
        patchwork::area(1, 1, 4, 6),
        patchwork::area(5, 1, 6, 6)
      )
    )
    patchPlot <- p_model / p_diff + layout_obj
    out <- list("data" = out_df, "plot" = patchPlot)
  }
  return(out)
}

#' ***********************************************************************************************
#' *************** `Plot difference in smooths` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for plotting spline_diff output
#'
#' @keywords internal
#' @noRd

.plot_gam_diff <- function(df, name_pattern = "Q_diff_", name_pattern2 = NULL) {
  x <- colnames(df)[ncol(df)]
  nms <- colnames(df)
  nms <- as.numeric(sub(name_pattern, "", nms[grepl(paste0("^", name_pattern), nms)]))
  cis <- numeric()
  i <- 1
  while (length(nms) > 1) {
    cis[i] <- max(nms) - min(nms)
    nms <- nms[-which.max(nms)]
    nms <- nms[-which.min(nms)]
    i <- i + 1
  }
  cis <- rev(cis)
  virPal <- viridis::viridis(n = length(cis), option = "mako", direction = -1, begin = 0.1)

  if (is.null(name_pattern2)) {
    lineLayer <- list(
      ggplot2::geom_line(ggplot2::aes(y = .data[["mu"]])),
      ggplot2::geom_hline(yintercept = 0, linetype = 5),
      ggplot2::labs(y = paste0(df[1, "g1"], " - ", df[1, "g2"]))
    )
  } else {
    lineLayer <- ggplot2::labs(y = "Model Prediction")
  }

  splinePlot <- ggplot2::ggplot(df, aes(x = .data[[x]])) +
    lapply(rev(seq_along(cis)), function(i) {
      ci <- cis[i]
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[[paste0(name_pattern, 1 - (1 - ((1 - ci) / 2)))]],
          ymax = .data[[paste0(name_pattern, 1 - ((1 - ci) / 2))]]
        ),
        fill = virPal[i], alpha = 0.5
      )
    }) +
    lineLayer +
    pcv_theme()

  if (!is.null(name_pattern2)) {
    virPal2 <- viridis::viridis(n = length(cis), option = "inferno", direction = -1, begin = 0.1)
    splinePlot <- splinePlot +
      lapply(rev(seq_along(cis)), function(i) {
        ci <- cis[i]
        ggplot2::geom_ribbon(
          ggplot2::aes(
            ymin = .data[[paste0(name_pattern2, 1 - (1 - ((1 - ci) / 2)))]],
            ymax = .data[[paste0(name_pattern2, 1 - ((1 - ci) / 2))]]
          ),
          fill = virPal2[i], alpha = 0.5
        )
      })
  }

  return(splinePlot)
}
