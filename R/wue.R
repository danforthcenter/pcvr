#' Calculate pseudo water use efficiency from phenotype and watering data
#' 
#' @description Water use efficiency (WUE) is the change in biomass per unit of water metabolized.
#' Using image based phenotypes and watering data we can calculate pseudo-WUE (pwue) over time.
#' Here area_pixels is used as a proxy for biomass and transpiration is approximated using watering data.
#' The equation is then \eqn{\frac{P_{t} - P_{t-1}}{W_{t_{end-1}}-W_{t_{start}} }}{P_[t] - P_[t-1] / W_[t_(end-1)]-W_[t_start]},
#' where P is the phenotype and W is the weight before watering.
#' 
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
#' @keywords read.csv, pcv, wide, long
#' @import jsonlite
#' @return A data frame containing the bellwether watering data joined
#'     to phenotype data with new columns for change in the phenotype,
#'     change in the pre-watering weight, and pseudo-water use efficiency (pWUE).
#' @examples
#' 
#' ## Not run:
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
#'       group=c("barcode", "rotation"), plot=TRUE)
#' phenotypes <- colnames(sv)[19:35]
#' phenoForm<-paste0("cbind(", paste0(phenotypes, collapse=", "), ")")
#' groupForm<-"DAS+DAP+barcode+genotype+fertilizer+timestamp"
#' sv<-aggregate(as.formula(paste0(phenoForm, "~", groupForm)), data=sv, mean, na.rm=TRUE)
#' sv<-bw.outliers(sv, phenotype="area_pixels",
#'       group = c("DAS", "genotype", "fertilizer"), plotgroup = c("barcode"))
#' water<-bw.water("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/metadata.json")
#' 
#' df<-sv
#' w<-water
#' pheno="area_pixels"
#' time="timestamp"
#' id =c("barcode")
#' offset = 0
#' waterCol="water_amount"
#' 
#' x<-pwue(df, w, pheno, time, id, offset, waterCol)
#' library(ggplot2)
#' ggplot(x, aes(x=timestamp, y=pWUE, group = barcode))+
#'   geom_line()+
#'   theme_minimal()
#'  
#' ## End(Not run)
#' 
#' 
#' @export

pwue<-function(df, w, pheno="area_pixels", time="timestamp", id="barcode",
               offset = 0, waterCol="water_amount"){
  
  if(length(time)==2){
    time1 = time[1]
    time2 = time[2]
  } else{
    time1 <- time2 <- time
  }
  if( !time1 %in% colnames(df) | !time2 %in% colnames(w) ){
    stop(paste0(paste0(time, collapse=", "), " must be in colnames of df and w"))
  }
  w<-data.table::setorderv(data.table::as.data.table(w), cols = c(id,time2))
  w <- w[w[[waterCol]]>0, ]
  df<-data.table::setorderv(data.table::as.data.table(df), cols = c(id,time1))
  
  ids <- intersect(unique(w[[id]]), unique(df[[id]]))
  
  out <- do.call(rbind, lapply(ids, function(iter_id){
    
    w_i <- w[w[[id]]==iter_id, ]
    df_i <- df[df[[id]]==iter_id, ]
    
    w_i <- data.table::setorderv(w_i, cols = c(time2))
    df_i <- data.table::setorderv(df_i, cols = c(time1))
    
    imaging_times <- unique(df_i[[time1]])
    
    wue_i <- do.call(rbind, lapply(1:length(imaging_times), function(t_i){
      
      start <- if(t_i==1){NA} else {imaging_times[(t_i-1)] }
      end <- imaging_times[t_i] - offset 
      
      if(!is.na(start)){
      w_i_t <- w_i[w_i[[time]]>start & w_i[[time]]<end, ]
      total_water_i <- max(c(sum(w_i_t[[waterCol]]), 1))
      pheno_diff <- max(c(as.numeric( df_i[df_i[[time1]] == imaging_times[t_i], get(pheno)] -
                                 df_i[df_i[[time1]]==start, get(pheno)] ), 0))
      }else{
        total_water_i = NA ; pheno_diff = NA
        start <- imaging_times[t_i] - offset 
      }
      
      row <- data.frame(total_water = total_water_i,
                        pheno_diff= pheno_diff,
                        start = start,
                        end = imaging_times[t_i])
      row$pWUE <- row$pheno_diff / row$total_water
      if(offset != 0){row$end_offset = end}
      return(row)
    }))
    
    iter_out <- cbind(df_i, wue_i)
    return(iter_out)
    
  }))
  return(out)
}



# pwue<-function(df, w, pheno="area_pixels", time="DAS", id="barcode"){
#   #* `format watering data`
#   w<-w[, !grepl("^timestamp$|^local_time$",colnames(w))]
#   w$snapshot_sorter <- as.numeric(sub("snapshot", "", w$snapshot))
#   #w<-data.table::setorderv(w, cols = c(id,time, 'snapshot_sorter')) # arrange plants by time
#   w<-w[!w$weight_before<0,] # remove errors
#   #* `format phenotype data`
#   df<-df[, !grepl("^timestamp$|^local_time$",colnames(df))]
#   #* `join datasets`
#   x<-data.table::as.data.table(merge(df, w, by=intersect(colnames(df), colnames(w)), all.x=TRUE))
#   x<-data.table::setorderv(x, cols = c(id,time, "snapshot_sorter")) # x can have duplicate pheno rows if w has >1 watering per day.
#   #* `calculate delta values and pwue`
#   x_deltas<-data.table::rbindlist(lapply(split(x, by=id), function(d){
#     d<-data.table::setorderv(d, cols = c(time))
#       # could lead the weight before or lag the weight after
#     #* `water transpired/lost`
#     d$water_used_between_waterings <- d$weight_after - data.table::shift(d$weight_before, n=1, type="lead")
#     #* `sum if there are multiple waterings per day`
#     if( any(duplicated(d[[time]])) ){
#       
#       d<-data.table::rbindlist(lapply(unique(d[[time]]), function(tm){
#         sub<-d[d[[time]]==tm, ]
#         sub<-data.table::setorderv(sub, cols = c("snapshot_sorter"))
#         out<-sub[1,]
#         if(all(is.na(sub$water_used_between_waterings))){
#           out$water_used_between_waterings <- NA
#         } else{
#           out$water_used_between_waterings <- sum(sub$water_used_between_waterings, na.rm=TRUE) 
#         }
#         out[[pheno]] <- mean(sub[[pheno]], na.rm=TRUE)
#         out
#       }))
#       
#     }
#     #* `calculate delta phenotype`
#     d[[paste0("delta_",pheno)]] <- pmax(d[[pheno]] - data.table::shift(d[[pheno]], n=1, type="lag"), 0, na.rm=FALSE)
#     #* `calculate pseudo-WUE`
#     d$pWUE <- d[[paste0("delta_",pheno)]] / d$water_used_between_waterings
#     d
#   }))
#   data.table::setDF(x_deltas)
#   return(x_deltas)
# }

