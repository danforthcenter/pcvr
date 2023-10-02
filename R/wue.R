#' Calculate pseudo water use efficiency from phenotype and watering data
#' 
#' @description Rate based water use efficiency (WUE) is the change in biomass per unit of water metabolized.
#' Using image based phenotypes and watering data we can calculate pseudo-WUE (pwue) over time.
#' Here area_pixels is used as a proxy for biomass and transpiration is approximated using watering data.
#' The equation is then \eqn{\frac{P_{t} - P_{t-1}}{W_{t_{end-1}}-W_{t_{start}} }}{P_[t] - P_[t-1] / W_[t_(end-1)]-W_[t_start]},
#' where P is the phenotype and W is the weight before watering.
#' 
#' Absolute value based WUE is the amount of water used to sustain a plants biomass over a given period.
#' The equation is then \eqn{\frac{P_{t}}{W_{t_{end-1}}-W_{t_{start}} }}{P_[t] / W_[t_(end-1)]-W_[t_start]}
#' 
#' @param df Dataframe containing wide single-value phenotype data.
#'     This should already be aggregated to one row per plant per day (angles/rotations combined).
#' @param w Watering data as returned from bw.water.
#' @param pheno Phenotype column name, defaults to "area_pixels"
#' @param time Variable(s) that identify a plant on a given day.
#'     Defaults to \code{c("barcode", "DAS")}.
#' @param id Variable(s) that identify a plant over time. Defaults to \code{"barcode"}.
#' @param offset Optionally you can specify how long before imaging a watering should not be taken into account. 
#' This defaults to 0, meaning that if a plant were watered directly before being imaged then that water would
#' be counted towards WUE between the current image and the prior one. This argument is taken to be in seconds.
#' @param waterCol Column containing watering amounts in \code{w}. This defaults to "watering_amount".
#' @param method Which method to use, options are "rate" and "abs". The "rate" method considers WUE as the
#'  change in a phenotype divided by the amount of water added. The "abs" method considers WUE as the amount of
#'  water used by a plant given its absolute size. The former is for questions more related to efficiency in using 
#' water to grow while the latter is more suited to questions about how efficient a plant is at maintaining size given
#' some amount of water.
#' @keywords WUE 
#' @import data.table
#' @return A data frame containing the bellwether watering data joined
#'     to phenotype data with new columns for change in the phenotype,
#'     change in the pre-watering weight, and pseudo-water use efficiency (pWUE).
#' @examples
#' 
#' library(data.table)
#' sv<-read.pcv(
#'    "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv",
#'    reader="fread")
#' sv$genotype = substr(sv$barcode, 3,5)
#' sv$genotype = ifelse(sv$genotype == "002", "B73",
#'      ifelse(sv$genotype == "003", "W605S",
#'      ifelse(sv$genotype == "004", "MM", "Mo17")))
#' sv$fertilizer = substr(sv$barcode, 8, 8)
#' sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
#'            ifelse(sv$fertilizer == "B", "50", "0"))
#' sv<-bw.time(sv, plantingDelay = 0,
#'       phenotype="area_pixels", cutoff=10, timeCol="timestamp",
#'       group=c("barcode", "rotation"), plot=FALSE)
#' phenotypes <- colnames(sv)[19:35]
#' phenoForm<-paste0("cbind(", paste0(phenotypes, collapse=", "), ")")
#' groupForm<-"DAS+DAP+barcode+genotype+fertilizer+timestamp"
#' sv<-aggregate(as.formula(paste0(phenoForm, "~", groupForm)), data=sv, mean, na.rm=TRUE)
#' sv<-bw.outliers(sv, phenotype="area_pixels",
#'       group = c("DAS", "genotype", "fertilizer"), plotgroup = c("barcode"), plot=FALSE)
#' water<-bw.water("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/metadata.json")
#' 
#' x<-pwue(df = sv, w = water, pheno = "area_pixels", 
#'   time="timestamp", id = "barcode", offset = 0,
#'   waterCol = "water_amount", method = "rate")
#'   
#' 
#' @export

pwue<-function(df, w, pheno="area_pixels", time="timestamp", id="barcode",
               offset = 0, waterCol="water_amount", method="rate"){
  if(length(time)==2){
    time1 = time[1]
    time2 = time[2]
  } else{
    time1 <- time2 <- time
  }
  if( !time1 %in% colnames(df) | !time2 %in% colnames(w) ){
    stop(paste0(paste0(time, collapse=", "), " must be in colnames of df and w"))
  }
  # if(!methods::is(w[[time2]], "character")){
  #   w[[time2]] <- as.Date(w[[time2]])
  # }
  
  w<-data.table::setorderv(data.table::as.data.table(w), cols = c(id,time2))
  w <- w[w[[waterCol]]>0, ]
  df<-data.table::setorderv(data.table::as.data.table(df), cols = c(id,time1))
  ids <- intersect(unique(w[,get(id)]), unique(df[,get(id)]))
  matched_method = match.arg(method, choices = c("rate", "abs"))
  
  if(matched_method == "abs"){
    out <- .absWUE(ids, w, df, offset, time1, time2, pheno, id, waterCol)
  } else if (matched_method=="rate"){
    out <- .rateWUE(ids, w, df, offset, time1, time2, pheno, id, waterCol)
  }
  return(as.data.frame(out))
}

#' Function to calculate rate based WUE
#' @keywords internal
#' @noRd

.rateWUE <- function(ids, w, df, offset, time1, time2, pheno, id, waterCol){
  out <- do.call(rbind, lapply(ids, function(iter_id){ # per id...
    print(iter_id)
    w_i <- w[w[,get(id)]==iter_id, ]
    df_i <- df[df[,get(id)]==iter_id, ]
    #* reorder watering and pheno data
    w_i <- data.table::setorderv(w_i, cols = c(time2))
    df_i <- data.table::setorderv(df_i, cols = c(time1))
    #* get unique imaging times
    imaging_times <- unique(df_i[[time1]])
    #* set water from before first image to zero so that offset does not grab water from before imaging starts.
    w_i[[waterCol]] <- as.numeric(ifelse(w_i[[time2]] < imaging_times[1], 0, w_i[[waterCol]]))
    #* per imaging time
    wue_i <- do.call(rbind, lapply(1:length(imaging_times), function(t_i){
      start <- if(t_i==1){NA} else { imaging_times[(t_i-1)] - offset }
      start_nonOffset <- if(t_i==1){NA} else { imaging_times[(t_i-1)] }
      end <- imaging_times[t_i] - offset 
      end_nonOffset <- imaging_times[t_i]
      print(t_i)
      if(!is.na(start)){
        w_i_t <- w_i[w_i[[time2]]>start & w_i[[time2]]<end, ]
        total_water_i <- max(c(sum(w_i_t[[waterCol]]), 1))
        pheno_diff <- max(c(as.numeric( df_i[df_i[[time1]] == end_nonOffset, get(pheno)] -
                                          df_i[df_i[[time1]]==start_nonOffset, get(pheno)] ), 0))
      }else{ total_water_i = NA ; pheno_diff = NA}
      
      row <- data.frame(total_water = total_water_i,
                        pheno_diff= pheno_diff,
                        start = start_nonOffset,
                        end = end_nonOffset,
                        timeLengthSeconds = as.numeric(end) - as.numeric(start),
                        offset = offset)
      row$pWUE <- row$pheno_diff / row$total_water
      return(row)
      
    }))
    iter_out <- cbind(df_i, wue_i)
    return(iter_out)
  }))
}

#' Function to calculate absolute value based WUE
#' @keywords internal
#' @noRd


.absWUE <- function(ids, w, df, offset, time1, time2, pheno, id, waterCol){
  out <- do.call(rbind, lapply(ids, function(iter_id){ # per id...
    w_i <- w[w[[id]]==iter_id, ]
    df_i <- df[df[[id]]==iter_id, ]
    #* reorder watering and pheno data
    w_i <- data.table::setorderv(w_i, cols = c(time2))
    df_i <- data.table::setorderv(df_i, cols = c(time1))
    #* get unique imaging times
    imaging_times <- unique(df_i[[time1]])
    #* set water from before first image to zero so that offset does not grab water from before imaging starts.
    w_i[[waterCol]] <- ifelse(w_i[[time2]] < imaging_times[1], 0, w_i[[time2]])
    #* per imaging time
    wue_i <- do.call(rbind, lapply(1:length(imaging_times), function(t_i){
      start <- if(t_i==1){NA} else {imaging_times[(t_i-1)] - offset }
      start_nonOffset <- if(t_i==1){NA} else {imaging_times[(t_i-1)] }
      end <- imaging_times[t_i] - offset 
      end_nonOffset <- imaging_times[t_i]
      
      if(!is.na(start)){
        w_i_t <- w_i[w_i[[time2]]>start & w_i[[time2]]<end, ]
        total_water_i <- max(c(sum(w_i_t[[waterCol]]), 1))
        pheno_iter <- max(df_i[df_i[[time1]] == end_nonOffset, get(pheno)], na.rm=TRUE)
        
      }else{ total_water_i = NA ; pheno_iter = NA}
      
      row <- data.frame(total_water = total_water_i,
                        pheno_iter= pheno_iter,
                        start = start_nonOffset,
                        end = end_nonOffset,
                        timeLengthSeconds = as.numeric(end) - as.numeric(start),
                        offset = offset)
      row$pWUE <- row$pheno_iter / row$total_water
      return(row)
    }))
    iter_out <- cbind(df_i, wue_i)
    return(iter_out)
  }))
}

