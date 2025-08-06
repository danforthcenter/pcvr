#' Calculate pseudo water use efficiency from phenotype and watering data
#'
#' @description Rate based water use efficiency (WUE) is the change in biomass per unit of water
#' metabolized. Using image based phenotypes and watering data we can calculate pseudo-WUE (pwue) over
#' time. Here area_pixels is used as a proxy for biomass and transpiration is approximated using
#' watering data. The equation is then
#' \eqn{
#' \frac{P_{t} - P_{t-1}}{W_{t_{end-1}}-W_{t_{start}} }}{P_[t] - P_[t-1] / W_[t_(end-1)]-W_[t_start]
#' },
#' where P is the phenotype and W is the weight before watering.
#'
#' Absolute value based WUE is the amount of water used to sustain a plants biomass over a given period.
#' The equation is then
#' \eqn{\frac{P_{t}}{W_{t_{end-1}}-W_{t_{start}} }}{P_[t] / W_[t_(end-1)]-W_[t_start]}
#'
#' @param df Dataframe containing wide single-value phenotype data.
#'     This should already be aggregated to one row per plant per day (angles/rotations combined).
#' @param w Watering data as returned from pcv.water.
#' @param pheno Phenotype column name, defaults to "area_pixels"
#' @param time Variable(s) that identify a plant on a given day. If variable names are different
#' between df and w then this can be a vector of two names.
#' @param id Variable(s) that identify a plant over time. Defaults to \code{"barcode"}.
#' If variable names are different between df and w then this can be a vector of two names.
#' @param offset Optionally you can specify how long before imaging a watering should not be taken into
#' account. This defaults to 0, meaning that if a plant were watered directly before being imaged then
#' that water would be counted towards WUE between the current image and the prior one.
#' This argument is taken to be in seconds.
#' @param pre_watering Column containing weight before watering in \code{w},
#' defaults to "weight_before".
#' @param post_watering Column containing weight after watering in \code{w},
#' defaults to "weight_after".
#' @param method Which method to use, options are "rate", "abs", and "ndt". The "rate" method considers
#' WUE as the change in a phenotype divided by the amount of water added. The "abs" method considers WUE
#' as the amount of water used by a plant given its absolute size.
#' The "ndt" method calculates normalized daily transpiration,
#' which is the reciprocal of the "abs" method. The "rate" method is for questions more
#' related to efficiency in using water to grow while "abs"/"ndt" are more suited to questions about
#' how efficient a plant is at maintaining size given some amount of water or how much water it uses
#' at a given size.
#' @keywords WUE
#' @import data.table
#' @return A data frame containing the watering data and
#'     to phenotype data with new columns for change in the phenotype (\code{pheno_diff}),
#'     amount of water used (\code{total_water}) over the interval between phenotype measurements
#'     (water added post to pre phenotype measurement), \code{start} and \code{end} times for the
#'     interval as well as their difference (\code{timeLengthSeconds}), and pseudo water use efficiency
#'     (\code{pWUE}).
#' @examples
#'
#' set.seed(123)
#' weight_before <- sort(round(rnorm(20, 100, 10), 0))
#' weight_after <- sort(round(rnorm(20, 120, 10), 0))
#' df <- data.frame(
#'   time = seq_along(weight_before),
#'   area_pixels = round(130 / (1 + exp( (12 - seq_along(weight_before))/3) ), 0),
#'   weight_before,
#'   weight_after,
#'   barcode = 1,
#'   other = "x"
#' )
#' ex <- pwue(df, time = "time", method = "rate", id = c("barcode", "other"))
#' w <- df[, c("time", "weight_before", "weight_after", "barcode")]
#' ex2 <- pwue(df, w, time = c("time", "time"), method = "abs")
#' ex3 <- pwue(df, w, time = c("time", "time"), method = "ndt")
#'
#' @export

pwue <- function(df, w = NULL, pheno = "area_pixels", time = "timestamp", id = "barcode",
                 offset = 0, pre_watering = "weight_before", post_watering = "weight_after",
                 method = "rate") {
  if (length(time) == 2) {
    time1 <- time[1]
    time2 <- time[2]
  } else {
    time1 <- time2 <- time
  }
  if (is.null(w)) {
    w <- df[, c(time2, id, pre_watering, post_watering)]
  }
  if (length(id) > 1) {
    df$temporary_pwue_id_column <- interaction(df[, c(id)])
    w$temporary_pwue_id_column <- interaction(w[, c(id)])
    id <- "temporary_pwue_id_column"
  }
  #* order data
  w <- data.table::setorderv(data.table::as.data.table(w), cols = c(id, time2))
  df <- data.table::setorderv(data.table::as.data.table(df), cols = c(id, time1))
  data.table::setkeyv(w, id)
  data.table::setkeyv(df, id)
  ids <- intersect(unique(w[[id]]), unique(df[[id]]))
  matched_method <- match.arg(method, choices = c("rate", "abs", "ndt"))
  #* apply method
  matched_fun <- get(paste0(".", matched_method, "WUE"))
  out <- matched_fun(ids, w, df, offset, time1, time2, pheno, id, pre_watering, post_watering)
  out <- as.data.frame(out)
  out <- out[, which(colnames(out) != "temporary_pwue_id_column")]
  # return data
  return(out)
}

#' Function to calculate rate based WUE
#' @keywords internal
#' @noRd

.rateWUE <- function(ids, w, df, offset, time1, time2, pheno, id, pre_watering, post_watering) {
  out <- do.call(rbind, lapply(ids, function(iter_id) { # per id...
    w_i <- w[list(iter_id)]
    df_i <- df[list(iter_id)]
    #* reorder watering and pheno data
    w_i <- data.table::setorderv(w_i, cols = c(time2))
    df_i <- data.table::setorderv(df_i, cols = c(time1))
    #* get post_watering_weight_[t-1] and pre_watering_weight_[t] difference as column
    w_i$deltaw <- data.table::shift(w_i[[post_watering]],
                                    n = 1, type = "lag") - w_i[[pre_watering]]
    #* get unique imaging times
    imaging_times <- unique(df_i[[time1]])
    #* set water from before first image to zero so that
    #* offset does not grab water from before imaging starts.
    w_i[["deltaw"]] <- as.numeric(ifelse(w_i[[time2]] < imaging_times[1], 0, w_i[["deltaw"]]))
    #* per imaging time
    wue_i <- do.call(rbind, lapply(seq_along(imaging_times), function(t_i) {
      start <- if (t_i == 1) {
        NA
      } else {
        imaging_times[(t_i - 1)] - offset
      }
      startNonOffset <- if (t_i == 1) {
        NA
      } else {
        imaging_times[(t_i - 1)]
      }
      end <- imaging_times[t_i] - offset
      endNonOffset <- imaging_times[t_i]
      if (!is.na(start)) {
        w_i_t <- w_i[w_i[[time2]] > start & w_i[[time2]] <= end, ]
        total_water_i <- sum(w_i_t[["deltaw"]])
        pheno_diff <- max(c(
          as.numeric(df_i[df_i[[time1]] == endNonOffset, get(pheno)] - df_i[
            df_i[[time1]] == startNonOffset, get(pheno)
          ]), 0
        ))
      } else {
        total_water_i <- NA
        pheno_diff <- NA
      }
      # make data to return
      row <- data.frame(
        total_water = total_water_i,
        pheno_diff = pheno_diff,
        start = startNonOffset,
        end = endNonOffset,
        timeLengthSeconds = as.numeric(end) - as.numeric(start),
        offset = offset
      )
      row$pWUE <- row$pheno_diff / row$total_water
      return(row)
    }))
    iter_out <- cbind(df_i, wue_i)
    return(iter_out)
  }))
  return(out)
}

#' Function to calculate absolute value based WUE
#' @keywords internal
#' @noRd


.absWUE <- function(ids, w, df, offset, time1, time2, pheno, id, pre_watering, post_watering) {
  out <- do.call(rbind, lapply(ids, function(iter_id) { # per id...
    w_i <- w[list(iter_id)]
    df_i <- df[list(iter_id)]
    #* reorder watering and pheno data
    w_i <- data.table::setorderv(w_i, cols = c(time2))
    df_i <- data.table::setorderv(df_i, cols = c(time1))
    #* get post_watering_weight_[t-1] and pre_watering_weight_[t] difference as column
    w_i$deltaw <- data.table::shift(w_i[[post_watering]],
                                    n = 1, type = "lag") - w_i[[pre_watering]]
    #* get unique imaging times
    imaging_times <- unique(df_i[[time1]])
    #* set water from before first image to zero so that offset
    #* does not grab water from before imaging starts.
    w_i[["deltaw"]] <- ifelse(w_i[[time2]] < imaging_times[1], 0, w_i[["deltaw"]])
    #* per imaging time
    wue_i <- do.call(rbind, lapply(seq_along(imaging_times), function(t_i) {
      start <- if (t_i == 1) {
        NA
      } else {
        imaging_times[(t_i - 1)] - offset
      }
      startNonOffset <- if (t_i == 1) {
        NA
      } else {
        imaging_times[(t_i - 1)]
      }
      end <- imaging_times[t_i] - offset
      endNonOffset <- imaging_times[t_i]
      # here we use the phenotype at this time, not the change in phenotype
      if (!is.na(start)) {
        w_i_t <- w_i[w_i[[time2]] > start & w_i[[time2]] <= end, ]
        total_water_i <- sum(w_i_t[["deltaw"]])
        pheno_iter <- max(df_i[df_i[[time1]] == endNonOffset, get(pheno)], na.rm = TRUE)
      } else {
        total_water_i <- NA
        pheno_iter <- NA
      }
      # make data to return
      row <- data.frame(
        total_water = total_water_i,
        pheno_iter = pheno_iter,
        start = startNonOffset,
        end = endNonOffset,
        timeLengthSeconds = as.numeric(end) - as.numeric(start),
        offset = offset
      )
      row$pWUE <- row$pheno_iter / row$total_water
      return(row)
    }))
    iter_out <- cbind(df_i, wue_i)
    return(iter_out)
  }))
  return(out)
}

#' Function to calculate normalized daily transpiration
#' @keywords internal
#' @noRd


.ndtWUE <- function(ids, w, df, offset, time1, time2, pheno, id, pre_watering, post_watering) {
  out <- .absWUE(ids, w, df, offset, time1, time2, pheno, id, pre_watering, post_watering)
  out$normalized_daily_transpiration <- 1 / out$pWUE
  return(out)
}
