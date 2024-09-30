#' Remove outliers from bellwether data using cook's distance
#'
#' @param df Data frame to use. Can be in long or wide format.
#' @param phenotype Column to use to classify outliers. If this is length > 1 then
#' it is taken as the multi-value traits to use. See examples.
#' @param naTo0 Logical, should NA values to changed to 0.
#' @param group  Grouping variables to find outliers as a character vector.
#' This is typically time  and design variables (DAS, genotype, treatment, etc).
#' These are used as predictors for `phenotype` in a generalized linear model.
#' @param cutoff Cutoff for something being an "outlier" expressed as a multiplier
#'  on the mean of Cooks Distance for this data. This defaults to 5, with higher values
#'  being more conservative (removing less of the data).
#' @param outlierMethod Method to be used in detecting outliers.
#'  Currently "cooks" and "mahalanobis" distances are supported, with "mahalanobis" only
#'  being supported for multi-value traits.
#' @param plotgroup Grouping variables for drawing plots if plot=TRUE.
#' Typically this is an identifier for images of a plant
#' over time and defaults to c('barcode',"rotation").
#' @param plot Logical, if TRUE then a list is returned with a ggplot and a dataframe.
#' @param x Optional specification for x axis variable if plot is true.
#' If left NULL (the default) then the first element of `group` is used.
#' @param traitCol Column with phenotype names, defaults to "trait".
#' This should generally not need to be changed from the default.
#'    If this and valueCol are present in colnames(df) then the data
#'    is assumed to be in long format.
#' @param valueCol Column with phenotype values, defaults to "value".
#' This should generally not need to be changed from the default.
#' @param labelCol Column with phenotype labels for long data, defaults to "label".
#' This should generally not need to be changed from the default.
#' @param idCol Column(s) that identify individuals over time.
#' Defaults to plotGroup.
#' @param ncp Optionally specify the number of principle components to be used for MV data outlier
#' detection with cooks distance. If left NULL (the default) then 3 will be used.
#' @param separate Optionally separate the data by some variable to speed up the modeling step.
#' If you have a design variable with
#' very many levels then it may be helpful to separate by that variable.
#' Note this will subset the data for each model so it will change
#' the outlier removal (generally to be more conservative).
#' @keywords ggplot outliers
#' @import ggplot2
#' @import data.table
#' @importFrom stats complete.cases cooks.distance glm as.formula lm mahalanobis cov
#' @examples
#'
#'
#' sv <- growthSim("logistic",
#'   n = 5, t = 20,
#'   params = list("A" = c(200, 160), "B" = c(13, 11), "C" = c(3, 3.5))
#' )
#' sv[130, ]$y <- 500
#' sv_res <- bw.outliers(
#'   df = sv, phenotype = "y", naTo0 = FALSE, cutoff = 15,
#'   group = c("time", "group"), outlierMethod = "cooks",
#'   plotgroup = "id", plot = TRUE
#' )
#' sv_res$plot
#' \donttest{
#' tryCatch(
#'   { # in case offline
#'     library(data.table)
#'     mvw <- read.pcv(paste0(
#'       "https://media.githubusercontent.com/media/joshqsumner/",
#'       "pcvrTestData/main/pcv4-multi-value-traits.csv"
#'     ), mode = "wide", reader = "fread")
#'     mvw$genotype <- substr(mvw$barcode, 3, 5)
#'     mvw$genotype <- ifelse(mvw$genotype == "002", "B73",
#'       ifelse(mvw$genotype == "003", "W605S",
#'         ifelse(mvw$genotype == "004", "MM", "Mo17")
#'       )
#'     )
#'     mvw$fertilizer <- substr(mvw$barcode, 8, 8)
#'     mvw$fertilizer <- ifelse(mvw$fertilizer == "A", "100",
#'       ifelse(mvw$fertilizer == "B", "50", "0")
#'     )
#'     mvw <- bw.time(mvw, timeCol = "timestamp", group = "barcode", plot = FALSE)
#'
#'     phenotypes <- which(grepl("hue_freq", colnames(mvw)))
#'
#'     mvw2 <- bw.outliers(
#'       df = mvw, phenotype = phenotypes, naTo0 = FALSE, outlierMethod = "cooks",
#'       group = c("DAS", "genotype", "fertilizer"), cutoff = 3, plotgroup = c("barcode", "rotation")
#'     )
#'
#'
#'     mvl <- read.pcv(paste0(
#'       "https://media.githubusercontent.com/media/joshqsumner/",
#'       "pcvrTestData/main/pcv4-multi-value-traits.csv"
#'     ), mode = "long")
#'     mvl$genotype <- substr(mvl$barcode, 3, 5)
#'     mvl$genotype <- ifelse(mvl$genotype == "002", "B73",
#'       ifelse(mvl$genotype == "003", "W605S",
#'         ifelse(mvl$genotype == "004", "MM", "Mo17")
#'       )
#'     )
#'     mvl$fertilizer <- substr(mvl$barcode, 8, 8)
#'     mvl$fertilizer <- ifelse(mvl$fertilizer == "A", "100",
#'       ifelse(mvl$fertilizer == "B", "50", "0")
#'     )
#'     mvl <- bw.time(mvl, timeCol = "timestamp", group = "barcode", plot = FALSE)
#'
#'     mvl2 <- bw.outliers(
#'       df = mvl, phenotype = "hue_frequencies", naTo0 = FALSE, outlierMethod = "cooks",
#'       group = c("DAS", "genotype", "fertilizer"), cutoff = 3, plotgroup = c("barcode", "rotation")
#'     )
#'   },
#'   error = function(e) {
#'     message(e)
#'   }
#' )
#' }
#'
#' @return The input dataframe with outliers removed and optionally a plot
#' (if a plot is returned then output is a list).
#' @export


bw.outliers <- function(df = NULL,
                        phenotype,
                        naTo0 = FALSE,
                        group = c(),
                        cutoff = 3,
                        outlierMethod = "cooks",
                        plotgroup = c("barcode", "rotation"),
                        plot = TRUE, x = NULL, traitCol = "trait", valueCol = "value",
                        labelCol = "label", idCol = NULL, ncp = NULL, separate = NULL) {
  wide <- .detectWideData(df, traitCol, valueCol)
  if (is.null(phenotype)) {
    stop("A phenotype must be provided")
  }

  if ((wide && length(phenotype) > 1) ||
    (!wide && length(unique(
      interaction(df[df[[traitCol]] == phenotype, colnames(df) %in% c(traitCol, labelCol)])
    )) > 1)) {
    mv <- TRUE
  } else {
    mv <- FALSE
  }

  if (is.null(idCol)) {
    idCol <- plotgroup
  }

  if (!is.null(separate)) {
    dfList <- split(df, df[[separate]])
    group <- group[!grepl(separate, group)]
  } else {
    dfList <- list(df)
  }

  resList <- .applyOutlierMethod(
    dfList, wide, mv, outlierMethod,
    naTo0, phenotype, group, cutoff,
    ncp, traitCol, valueCol, labelCol, idCol
  )

  df <- do.call(rbind, lapply(resList, function(res) {
    res[["data"]]
  }))
  pctRm <- do.call(rbind, lapply(seq_along(resList), function(i) {
    data.frame(i = i, pctRm = resList[[i]][["pctRm"]])
  }))

  out <- df[which(!df$outlier), -which(grepl("outlier", colnames(df)))]

  removedInteractions <- do.call(rbind, lapply(
    resList,
    function(d) {
      tab <- as.data.frame(table(d$data[, c(separate, group)]))
      colnames(tab) <- c(separate, group, "Freq")
      return(tab[tab$Freq == 0, ])
    }
  ))
  if (nrow(removedInteractions) > 0) {
    warning(paste0(nrow(removedInteractions), " groupings had all observations removed"))
  }

  if (plot) {
    p <- .outlierPlottingHelper(wide, mv, df, plotgroup, group, x, phenotype, traitCol, valueCol, pctRm)
    out <- list("data" = out, "plot" = p)
  }

  return(out)
}


#' ***********************************************************************************************
#' *************** `Apply outlier methods` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function to apply outlier methods
#'
#' @keywords internal
#' @noRd

.applyOutlierMethod <- function(dfList, wide, mv, outlierMethod,
                                naTo0, phenotype, group, cutoff,
                                ncp, traitCol, valueCol, labelCol, idCol) {
  resList <- lapply(dfList, function(df) {
    mv_label <- if (mv) {
      "mv"
    } else {
      "sv"
    }
    wide_label <- if (wide) {
      ".wide"
    } else {
      ".long"
    }
    matched_fun <- get(paste(wide_label, mv_label, outlierMethod, "outliers", sep = "_"))
    res <- matched_fun(df, naTo0, phenotype, group, cutoff, ncp, traitCol, valueCol, labelCol, idCol)
    return(res)
  })
  return(resList)
}


#' ***********************************************************************************************
#' *************** `wide data detection` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in wide SV data.
#'
#' @keywords internal
#' @noRd

.detectWideData <- function(df, traitCol, valueCol) {
  if (all(c(traitCol, valueCol) %in% colnames(df))) {
    wide <- FALSE
  } else {
    wide <- TRUE
  }
  return(wide)
}


#' ***********************************************************************************************
#' *************** `plotting helper` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in wide SV data.
#'
#' @keywords internal
#' @noRd

.outlierPlottingHelper <- function(wide, mv, df, plotgroup, group,
                                   x, phenotype, traitCol, valueCol, pctRm) {
  if (is.null(x)) {
    x <- group[1]
  }
  if (wide && !mv) {
    df$grouping <- interaction(df[, plotgroup])
    outPlotData <- df[!df$outlier, ]
    rmdfPlotData <- df[df$outlier, ]
    p <- ggplot2::ggplot() +
      ggplot2::facet_wrap(stats::as.formula(paste0("~", paste(group[-1], collapse = "+")))) +
      ggplot2::geom_line(data = df, ggplot2::aes(
        x = .data[[x]], y = .data[[phenotype]],
        group = .data[["grouping"]]
      ), linewidth = 0.25) +
      ggplot2::labs(title = paste0("~", round(mean(pctRm$pctRm), 3), "% Removed")) +
      pcv_theme()

    yLims <- ggplot2::layer_scales(p)$y$range$range

    p <- p +
      ggplot2::geom_point(
        data = rmdfPlotData, ggplot2::aes(
          x = .data[[x]],
          y = .data[[phenotype]]
        ),
        color = "red", size = 0.5
      ) +
      ggplot2::coord_cartesian(ylim = yLims)
  } else if (!wide && !mv) {
    plotdf <- df[df[[traitCol]] == phenotype, ]
    plotdf$grouping <- interaction(plotdf[, plotgroup])
    outPlotData <- plotdf[!plotdf$outlier, ]
    rmdfPlotData <- plotdf[plotdf$outlier, ]
    p <- ggplot2::ggplot() +
      ggplot2::facet_wrap(stats::as.formula(paste0("~", paste(group[-1], collapse = "+")))) +
      ggplot2::geom_line(data = plotdf, ggplot2::aes(
        x = .data[[x]], y = .data[[valueCol]],
        group = .data[["grouping"]]
      ), linewidth = 0.25) +
      ggplot2::labs(title = paste0("~", round(mean(pctRm$pctRm), 3), "% Removed")) +
      pcv_theme()

    yLims <- ggplot2::layer_scales(p)$y$range$range

    p <- p +
      ggplot2::geom_point(
        data = rmdfPlotData, ggplot2::aes(
          x = .data[[x]],
          y = .data[[valueCol]]
        ),
        color = "red", size = 0.5
      ) +
      ggplot2::coord_cartesian(ylim = yLims)
  } else if (wide && mv) {
    plotdf <- suppressWarnings(
      as.data.frame(
        data.table::melt(data.table::as.data.table(df),
          measure.vars = phenotype,
          variable.name = traitCol,
          value.name = valueCol
        )
      )
    )

    plotdf$bin <- as.numeric(regmatches(plotdf$trait, regexpr("[0-9]+", plotdf$trait)))

    plotdf$grouping <- interaction(plotdf[, plotgroup])
    outPlotData <- plotdf[!plotdf$outlier, ]
    rmdfPlotData <- plotdf[plotdf$outlier, ]

    p <- ggplot2::ggplot() +
      ggplot2::facet_wrap(stats::as.formula(paste0("~", paste(group[-1], collapse = "+")))) +
      ggplot2::geom_col(
        data = rmdfPlotData, ggplot2::aes(x = .data[["bin"]], y = .data[[valueCol]]),
        position = "identity",
        fill = "red", alpha = 0.25
      ) +
      ggplot2::geom_col(
        data = outPlotData, ggplot2::aes(x = .data[["bin"]], y = .data[[valueCol]]),
        position = "identity",
        alpha = 0.25
      ) +
      ggplot2::labs(title = paste0("~", round(mean(pctRm$pctRm), 3), "% Removed")) +
      pcv_theme()
  } else if (!wide && mv) {
    plotdf <- df
    plotdf$grouping <- interaction(plotdf[, plotgroup])
    outPlotData <- plotdf[!plotdf$outlier, ]
    rmdfPlotData <- plotdf[plotdf$outlier, ]

    p <- ggplot2::ggplot() +
      ggplot2::facet_wrap(stats::as.formula(paste0("~", paste(group[-1], collapse = "+")))) +
      ggplot2::geom_col(
        data = rmdfPlotData, ggplot2::aes(x = .data[[traitCol]], y = .data[[valueCol]]),
        position = "identity",
        fill = "red", alpha = 0.25
      ) +
      ggplot2::geom_col(
        data = outPlotData, ggplot2::aes(x = .data[[traitCol]], y = .data[[valueCol]]),
        position = "identity",
        alpha = 0.25
      ) +
      ggplot2::labs(title = paste0("~", round(mean(pctRm$pctRm), 3), "% Removed")) +
      pcv_theme()
  }
  return(p)
}


#' ***********************************************************************************************
#' *************** `wide SV cooks distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in wide SV data.
#'
#' @keywords internal
#' @noRd

.wide_sv_cooks_outliers <- function(df, naTo0, phenotype, group, cutoff, ncp,
                                    traitCol, valueCol, labelCol, idCol) {
  if (naTo0) {
    df[[phenotype]][is.na(df[[phenotype]])] <- 0
  }
  df <- df[complete.cases(df[, c(phenotype, group)]), ]
  outlierForm <- paste("as.numeric(", phenotype, ")~", paste(paste0("as.factor(", group, ")"),
    collapse = ":"
  ))
  cooksd <- cooks.distance(glm(data = df, as.formula(outlierForm)))
  outlierCutoff <- cutoff * mean(cooksd, na.rm = TRUE)
  cooksd[is.na(cooksd)] <- outlierCutoff - 0.1 # keeping NAs
  cooksd_df <- data.frame("outlier" = cooksd)
  df <- cbind(df, cooksd_df)
  df$outlier <- df$outlier > outlierCutoff
  pctRm <- 100 * (round(nrow(df[df$outlier, ]) / nrow(df), 5))

  return(list("data" = df, "pctRm" = pctRm))
}


#' ***********************************************************************************************
#' *************** `long SV cooks distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in wide SV data.
#'
#' @keywords internal
#' @noRd

.long_sv_cooks_outliers <- function(df, naTo0, phenotype, group, cutoff, ncp,
                                    traitCol, valueCol, labelCol, idCol) {
  if (naTo0) {
    df[df[[traitCol]] == phenotype, valueCol][is.na(df[df[[traitCol]] == phenotype, valueCol])] <- 0
  }
  subdf <- df[complete.cases(df[
    df[[traitCol]] == phenotype,
    c(valueCol, traitCol, group)
  ]) & df[[traitCol]] == phenotype, ]
  outlierForm <- paste("as.numeric(", valueCol, ")~", paste(paste0("as.factor(", group, ")"),
    collapse = ":"
  ))
  cooksd <- cooks.distance(glm(data = subdf, as.formula(outlierForm)))
  outlierCutoff <- cutoff * mean(cooksd, na.rm = TRUE)
  cooksd[is.na(cooksd)] <- outlierCutoff - 0.1 # keeping NAs by assigning a value below cutoff.
  cooksd_df <- data.frame("outlier" = cooksd)
  subdf <- cbind(subdf, cooksd_df)
  subdf <- subdf[, c(group, idCol, "outlier")]
  subdf <- subdf[!duplicated(subdf[, c(group, idCol)]), ]
  subdf$outlier <- subdf$outlier > outlierCutoff
  pctRm <- 100 * (round(nrow(subdf[subdf$outlier, ]) / nrow(subdf), 5))
  #* take IDs using plotgroup and label all phenotype rows
  df <- merge(df, subdf, all.x = TRUE)
  return(list("data" = df, "pctRm" = pctRm))
}

#' ***********************************************************************************************
#' *************** `wide MV cooks distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in wide MV data.
#'
#' @keywords internal
#' @noRd

.wide_mv_cooks_outliers <- function(df, naTo0, phenotype, group, cutoff, ncp,
                                    traitCol, valueCol, labelCol, idCol) {
  if (naTo0) {
    df[, phenotype][is.na(df[, phenotype])] <- 0
  }

  phenos_df <- df[, phenotype]
  if (is.null(ncp)) {
    ncp <- min(min(dim(phenos_df)) - 1, 3)
  }
  pca <- FactoMineR::PCA(phenos_df, ncp = ncp, graph = FALSE)
  coords <- as.data.frame(pca$ind)
  coords <- coords[, grepl("coord", colnames(coords))]
  colnames(coords) <- gsub("coord.Dim.", "pc", colnames(coords))
  pca_cols <- colnames(coords)
  df <- cbind(df, coords)

  df <- df[complete.cases(df[, c(pca_cols, group)]), ]

  outlierForm <- paste(
    "cbind(", paste0("pc", 1:ncp, collapse = ","), ")~",
    paste(paste0("as.factor(", group, ")"), collapse = ":")
  )
  cooksd <- cooks.distance(lm(data = df, as.formula(outlierForm)))

  df <- df[, -which(colnames(df) %in% c(paste0("pc", 1:ncp)))]

  if (length(cutoff) == 1) {
    cutoff <- rep(cutoff, ncp)
  }
  outlierCutoffs <- cutoff * colMeans(cooksd, na.rm = TRUE)

  outlierMatrix <- do.call(cbind, lapply(seq_len(ncol(cooksd)), function(i) {
    cooks_vec <- cooksd[, i] # this is causing a problem
    cooks_vec[is.na(cooks_vec)] <- outlierCutoffs[i] - 0.1
    setNames(data.frame(cooks_vec > outlierCutoffs[i]), paste0("outlier_", i))
  }))
  outlierMatrix$outlier <- unlist(lapply(seq_len(nrow(outlierMatrix)), function(i) {
    any(outlierMatrix[i, ]) # could be a more nuanced rule
  }))

  df <- cbind(df, outlierMatrix)
  pctRm <- 100 * (round(nrow(df[df$outlier, ]) / nrow(df), 5))
  return(list("data" = df, "pctRm" = pctRm))
}


#' ***********************************************************************************************
#' *************** `long MV cooks distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in long MV data.
#'
#' @keywords internal
#' @noRd

.long_mv_cooks_outliers <- function(df, naTo0, phenotype, group, cutoff, ncp,
                                    traitCol, valueCol, labelCol, idCol) {
  #* widen data
  dcast_form <- as.formula(paste0("... ~ ", traitCol, "+", labelCol))
  dfw <- as.data.frame(data.table::dcast(data.table::as.data.table(df[df[[traitCol]] == phenotype, ]),
    dcast_form,
    value.var = valueCol, sep = "."
  ))
  phenotypew <- which(grepl(phenotype, colnames(dfw)))
  #* call .wide method on dfw
  wide_res <- .wide_mv_cooks_outliers(dfw, naTo0, phenotypew, group, cutoff, ncp)
  pctRm <- wide_res[["pctRm"]]
  sub_df <- wide_res[["data"]]
  #* label long data based on .wide output
  kept <- unique(as.character(interaction(sub_df[!sub_df$outlier, c(group, idCol)])))
  df_ids <- as.character(interaction(df[, c(group, idCol)]))
  df$outlier <- !df_ids %in% kept

  return(list("data" = df, "pctRm" = pctRm))
}

#' ***********************************************************************************************
#' *************** `wide MV mahalanobis distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in wide MV data.
#'
#' @keywords internal
#' @noRd


.wide_mv_mahalanobis_outliers <- function(df, naTo0, phenotype, group, cutoff, ncp,
                                          traitCol, valueCol, labelCol, idCol) {
  if (naTo0) {
    df[, phenotype][is.na(df[, phenotype])] <- 0
  }

  phenos_df <- df[, phenotype]
  phenos_df <- phenos_df[, colSums(phenos_df) > 1]
  mahala_center <- colMeans(phenos_df, na.rm = TRUE)
  mahala_cov <- stats::cov(phenos_df)
  m <- stats::mahalanobis(phenos_df, mahala_center, mahala_cov)

  df$mahal <- m
  group_inter <- unique(as.character(interaction(df[, group])))

  df_out <- do.call(rbind, lapply(group_inter, function(grp) {
    subMeta <- df[interaction(df[, group]) == grp, ]
    subMeta$outlier <- ifelse(subMeta$mahal > cutoff * mean(subMeta$mahal, na.rm = TRUE), TRUE, FALSE)
    subMeta
  }))

  pctRm <- 100 * (round(nrow(df_out[df_out$outlier, ]) / nrow(df_out), 5))

  return(list("data" = df_out, "pctRm" = pctRm))
}

#' ***********************************************************************************************
#' *************** `long MV mahalanobis distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in long MV data.
#'
#' @keywords internal
#' @noRd

.long_mv_mahalanobis_outliers <- function(df, naTo0, phenotype, group, cutoff,
                                          ncp, traitCol, valueCol, labelCol, idCol) {
  #* widen data
  dcast_form <- as.formula(paste0("... ~ ", traitCol, "+", labelCol))
  dfw <- as.data.frame(data.table::dcast(data.table::as.data.table(df[df[[traitCol]] == phenotype, ]),
    dcast_form,
    value.var = valueCol, sep = "."
  ))
  phenotypew <- which(grepl(phenotype, colnames(dfw)))
  #* call .wide method on dfw
  wide_res <- .wide_mv_mahalanobis_outliers(dfw, naTo0, phenotypew, group, cutoff)
  pctRm <- wide_res[["pctRm"]]
  sub_df <- wide_res[["data"]]
  #* label long data based on .wide output
  kept <- unique(as.character(interaction(sub_df[!sub_df$outlier, c(group, idCol)])))
  df_ids <- as.character(interaction(df[, c(group, idCol)]))
  df$outlier <- !df_ids %in% kept

  return(list("data" = df, "pctRm" = pctRm))
}
