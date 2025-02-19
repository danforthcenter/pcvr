#' Variance partitioning using Full Random Effects Models
#'
#' Variance partitioning for phenotypes (over time) using fully random effects models
#'
#' @param df Dataframe containing phenotypes and design variables, optionally over time.
#' @param des Design variables to partition variance for as a character vector.
#' @param phenotypes Phenotype column names (data is assumed to be in wide format)
#' as a character vector.
#' @param timeCol A column of the data that denotes time for longitudinal experiments.
#'   If left NULL (the default) then all data is assumed to be from one timepoint.
#' @param cor Logical, should a correlation plot be made? Defaults to TRUE.
#' @param returnData Logical, should the used to make plots be returned? Defaults to FALSE.
#' @param combine Logical, should plots be combined with patchwork?
#'  Defaults to TRUE, which works well when there is a single timepoint being used.
#' @param markSingular Logical, should singular fits be marked in the variance explained plot?
#'  This is FALSE by default but it is good practice to check with TRUE in some situations.
#'   If TRUE this will add white markings to the plot where models had singular fits,
#'    which is the most common problem with this type of model.
#' @param time If the data contains multiple timepoints then which should be used?
#'  This can be left NULL which will use the maximum time if \code{timeCol} is specified.
#'   If a single number is provided then that time value will be used.
#'    Multiple numbers will include those timepoints.
#'     The string "all" will include all timepoints.
#' @param time_format Format for non-integer time, passed to \code{strptime},
#' defaults to "\%Y-\%m-\%d".
#' @param ... Additional arguments passed to \code{lme4::lmer}.
#'
#' @import lme4
#' @import ggplot2
#' @importFrom scales percent_format
#' @importFrom stats as.formula na.omit
#' @import patchwork
#'
#' @return Returns either a plot (if returnData=FALSE) or a list with a plot and
#' data/a list of dataframes (depending on returnData and cor).
#'
#' @examples
#'
#'
#' library(data.table)
#' set.seed(456)
#' df <- data.frame(
#'   genotype = rep(c("g1", "g2"), each = 10),
#'   treatment = rep(c("C", "T"), times = 10),
#'   time = rep(c(1:5), times = 2),
#'   date_time = rep(paste0("2024-08-", 21:25), times = 2),
#'   pheno1 = rnorm(20, 10, 1),
#'   pheno2 = sort(rnorm(20, 5, 1)),
#'   pheno3 = sort(runif(20))
#' )
#' out <- frem(df, des = "genotype", phenotypes = c("pheno1", "pheno2", "pheno3"), returnData = TRUE)
#' lapply(out, class)
#' frem(df,
#'   des = c("genotype", "treatment"), phenotypes = c("pheno1", "pheno2", "pheno3"),
#'   cor = FALSE
#' )
#' frem(df,
#'   des = "genotype", phenotypes = c("pheno1", "pheno2", "pheno3"),
#'   combine = FALSE, timeCol = "time", time = "all"
#' )
#' frem(df,
#'   des = "genotype", phenotypes = c("pheno1", "pheno2", "pheno3"),
#'   combine = TRUE, timeCol = "time", time = 1
#' )
#' frem(df,
#'   des = "genotype", phenotypes = c("pheno1", "pheno2", "pheno3"),
#'   cor = FALSE, timeCol = "time", time = 3:5, markSingular = TRUE
#' )
#' df[df$time == 3, "genotype"] <- "g1"
#' frem(df,
#'   des = "genotype", phenotypes = c("pheno1", "pheno2", "pheno3"),
#'   cor = FALSE, timeCol = "date_time", time = "all", markSingular = TRUE
#' )
#'
#' @export

frem <- function(df, des, phenotypes, timeCol = NULL, cor = TRUE, returnData = FALSE, combine = TRUE,
                 markSingular = FALSE, time = NULL, time_format = "%Y-%m-%d", ...) {
  dummyX <- FALSE
  if (is.null(timeCol)) {
    timeCol <- "dummy_x_axis"
    df[[timeCol]] <- 1
    dummyX <- TRUE
  }
  #* `Format a time column if non-integer`
  formatted <- .formatNonIntegerTime(df, timeCol, format = time_format, index = NULL, digits = 0)
  df <- formatted$data
  timeCol <- formatted$timeCol
  #* `Make formulas`
  ext <- FALSE
  if (length(des) == 2) {
    ind_fmla <- paste0("(1|", des[1], ")+(1|", des[2], ")+(1|", des[1], ":", des[2], ")")
  } else {
    ind_fmla <- paste(paste0("(1|", des, ")"), collapse = "+")
    ext <- TRUE
  }
  #* `Find time and subset data`
  if (is.null(time)) {
    dat <- na.omit(df[df[[timeCol]] == max(df[[timeCol]]), c(des, phenotypes, timeCol)])
    LONGITUDINAL <- FALSE
  } else if (is.numeric(time) && length(time) == 1) {
    dat <- na.omit(df[df[[timeCol]] == time, c(des, phenotypes, timeCol)])
    LONGITUDINAL <- FALSE
  } else if (is.numeric(time) && length(time) > 1) {
    LONGITUDINAL <- TRUE
    dat <- na.omit(df[df[[timeCol]] %in% time, c(des, phenotypes, timeCol)])
  } else if (time == "all") {
    LONGITUDINAL <- TRUE
    dat <- na.omit(df[, c(des, phenotypes, timeCol)])
  }

  #* `Partition Variance`

  H2 <- .partitionVarianceFrem(dat, timeCol, phenotypes, ind_fmla, ext, des, ...)

  if (!ext) {
    colnames(H2) <- c(des[1], des[2], "Interaction", "Unexplained", timeCol, "singular", "Phenotypes")
  } else {
    colnames(H2) <- c(des, "Unexplained", timeCol, "singular", "Phenotypes")
  }
  ordering <- H2[H2[[timeCol]] == max(H2[[timeCol]]), ]
  H2$Phenotypes <- ordered(H2$Phenotypes, levels = ordering$Phenotypes[order(ordering$Unexplained)])
  h2_melt <- data.frame(data.table::melt(as.data.table(H2), id = c("Phenotypes", "singular", timeCol)))

  if (!ext) {
    h2_melt$variable <- ordered(h2_melt$variable,
      levels = c("Unexplained", des[1], des[2], "Interaction")
    )
  } else {
    h2_melt$variable <- ordered(h2_melt$variable,
      levels = c("Unexplained", des)
    )
  }
  anova_dat <- h2_melt

  #* `Plot Variance`

  plotHelperOutputs <- .fremPlotHelper(
    LONGITUDINAL, anova_dat, markSingular, dummyX, timeCol, dat, phenotypes,
    cor, combine
  )
  plot <- plotHelperOutputs[["plot"]]
  x <- plotHelperOutputs[["x"]]
  cor <- plotHelperOutputs[["cor"]]
  out <- .fremCollectOutputs(returnData, cor, H2, x, plot)
  return(out)
}

#' helper function to plot FREM results
#' @keywords internal
#' @noRd

.fremCollectOutputs <- function(returnData, cor, H2, x, plot) {
  if (returnData) {
    if (cor) {
      out_data <- list("variance" = H2, "cor" = x)
    } else {
      out_data <- H2
    }
    out <- list("plot" = plot, "data" = out_data)
  } else {
    out <- plot
  }
  return(out)
}


#' helper function to plot FREM results
#' @keywords internal
#' @noRd

.fremPlotHelper <- function(LONGITUDINAL, anova_dat, markSingular, dummyX, timeCol, dat, phenotypes,
                            cor, combine) {
  x <- NULL
  if (!LONGITUDINAL) {
    p <- ggplot2::ggplot(data = anova_dat) +
      ggplot2::geom_col(ggplot2::aes(y = .data$Phenotypes, x = .data$value, fill = .data$variable)) +
      ggplot2::xlab("Variance Explained") +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "", reverse = TRUE)) +
      ggplot2::theme_minimal() +
      ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), labels = scales::percent_format()) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 14),
        axis.title.y = ggplot2::element_blank()
      ) +
      ggplot2::theme(axis.ticks.length = ggplot2::unit(0.2, "cm")) +
      ggplot2::theme(legend.position = "bottom")

    if (markSingular) {
      p <- p + ggplot2::geom_point(
        data = anova_dat[as.logical(anova_dat$singular) & anova_dat$variable == "Unexplained", ],
        aes(x = 0.99, y = .data$Phenotypes), color = "white", shape = 0
      )
    }
    if (dummyX) {
      p <- p + ggplot2::theme(axis.title.x.bottom = ggplot2::element_blank())
    }
  } else {
    p <- ggplot(data = anova_dat, aes(x = .data[[timeCol]], y = .data$value, fill = .data$variable)) +
      ggplot2::geom_area() +
      ggplot2::facet_wrap(~ .data$Phenotypes) +
      ggplot2::ylab("Variance Explained") +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "", reverse = TRUE)) +
      ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), labels = scales::percent_format()) +
      ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), labels = ~ round(., 1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_blank(),
        legend.position = "bottom"
      )
    if (markSingular) {
      p <- p + ggplot2::geom_vline(
        data = anova_dat[anova_dat$variable == "Unexplained" & as.logical(anova_dat$singular), ],
        ggplot2::aes(xintercept = .data[[timeCol]]), color = "white", linetype = 5, linewidth = 0.1
      )
    }
    if (dummyX) {
      p <- p + ggplot2::theme(axis.title.x.bottom = ggplot2::element_blank())
    }
  }

  #* `Get Correlations`
  if (length(phenotypes) == 1) {
    cor <- FALSE
  }
  if (cor) {
    corr <- cor(apply(dat[, (colnames(dat) %in% phenotypes)], 2, as.numeric),
      use = "complete.obs", method = "spearman"
    )
    unexp <- anova_dat[anova_dat$variable == "Unexplained" &
                         anova_dat[[timeCol]] == max(anova_dat[[timeCol]]), ]
    corr <- corr[
      as.character(unexp$Phenotypes[order(unexp$value)]),
      as.character(unexp$Phenotypes[order(unexp$value)])
    ]
    x <- na.omit(
      as.data.frame(
        suppressWarnings(data.table::melt(as.data.table(corr), variable.name = "Var2"))
      )
    )
    x$Var1 <- rep(rownames(corr), length.out = nrow(x))

    x$Var1 <- ordered(x$Var1, levels = unexp$Phenotypes[order(unexp$value)])
    x$Var2 <- ordered(x$Var2, levels = unexp$Phenotypes[order(unexp$value)])

    p2 <- ggplot2::ggplot(x, ggplot2::aes(.data$Var1, .data$Var2)) +
      ggplot2::geom_point(ggplot2::aes(color = .data$value), size = 4) +
      ggplot2::scale_color_gradient2(
        limits = c(-1, 1),
        midpoint = 0
      ) +
      ggplot2::theme_minimal() +
      ggplot2::guides(color = ggplot2::guide_colourbar(barwidth = 15, title = NULL)) +
      ggplot2::labs(x = "Correlations", y = "") +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x.bottom = ggplot2::element_text(angle = 90, hjust = 1)
      )
    if (combine) {
      p2 <- p2 + ggplot2::theme(axis.text.y.left = ggplot2::element_blank())
    }
  }

  if (cor) {
    if (combine) {
      plot <- p + p2
    } else {
      plot <- list(p, p2)
    }
  } else {
    plot <- p
  }
  return(list("plot" = plot, "x" = x, "cor" = cor))
}

#' helper function to run models in frem
#' @keywords internal
#' @noRd

.partitionVarianceFrem <- function(dat, timeCol, phenotypes, ind_fmla, ext, des, ...) {
  H2 <- data.frame(do.call(rbind, lapply(sort(unique(dat[[timeCol]])), function(tm) {
    sub <- dat[dat[[timeCol]] == tm, ]
    des_bools <- unlist(lapply(des, function(var) {
      des_numbers <- as.numeric(table(sub[[var]]))
      return(sum(des_numbers != 0) > 1)
    }))
    if (any(!des_bools)) {
      message(paste("Skipping", timeCol, tm, "as grouping contains a variable that is singular"))
      return(NULL)
    }
    pheno_df <- do.call(rbind, lapply(phenotypes, function(e) {
      fmla <- as.formula(paste0("as.numeric(", e, ") ~ ", ind_fmla))
      model <- suppressMessages(lme4::lmer(fmla, data = sub, ...))
      if (length(model@optinfo$conv$lme4) >= 1) {
        singular <- any(grepl("isSingular", model@optinfo$conv$lme4$messages))
      } else {
        singular <- FALSE
      }
      re <- lme4::VarCorr(model)
      res <- attr(lme4::VarCorr(model), "sc")^2

      if (!ext) {
        interaction.var <- as.numeric(attr(re[[which(grepl(":", names(re)))]], "stddev"))^2
        des1.var <- as.numeric(attr(re[[des[1]]], "stddev"))^2
        des2.var <- as.numeric(attr(re[[des[2]]], "stddev"))^2

        tot.var <- sum(as.numeric(re), res)
        unexp <- 1 - sum(as.numeric(re)) / sum(as.numeric(re), res)

        h2 <- c(
          (des1.var / tot.var),
          (des2.var / tot.var),
          (interaction.var / tot.var),
          unexp,
          tm,
          singular
        )
      } else {
        var <- lapply(des, function(i) {
          return(as.numeric(attr(re[[i]], "stddev"))^2)
        })

        tot.var <- sum(as.numeric(re), res)
        unexp <- 1 - sum(as.numeric(re)) / sum(as.numeric(re), res)

        h2 <- c(unlist(var) / tot.var, unexp, tm, singular)
      }
      return(h2)
    }))
    return(pheno_df)
  })))
  H2$Phenotypes <- rep(phenotypes, length.out = nrow(H2))
  return(H2)
}
