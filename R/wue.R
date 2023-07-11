#' Calculate pseudo water use efficiency from phenotype and watering data
#' 
#' @description Water use efficiency (WUE) is the change in biomass per unit of water metabolized.
#' Using image based phenotypes and watering data we can calculate pseudo-WUE (pwue) over time.
#' Here area is used as a proxy for biomass and transpiration is approximated using watering data.
#' The equation is then \eqn{\frac{P_{t} - P_{t-1}}{W_{t_{end-1}}-W_{t_{start}} }}{P_[t] - P_[t-1] / W_[t_(end-1)]-W_[t_start]},
#' where P is the phenotype and W is the weight before watering.
#' 
#' 
#' @param df Dataframe containing wide single-value phenotype data.
#'     This should already be aggregated to one row per plant per day (angles/rotations combined).
#' @param w Watering data as returned from bw.water.
#' @param pheno Phenotype column name, defaults to "area.pixels"
#' @param time Variable(s) that identify a plant on a given day.
#'     Defaults to \code{c("barcode", "DAS")}.
#' @param id Variable(s) that identify a plant over time. Defaults to \code{"barcode"}.
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
#'    "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv",
#'    mode="wide", singleValueOnly = T, reader="fread")
#' sv$genotype = substr(sv$barcode, 3,5)
#' sv$genotype = ifelse(sv$genotype == "002", "B73",
#'      ifelse(sv$genotype == "003", "W605S",
#'      ifelse(sv$genotype == "004", "MM", "Mo17")))
#' sv$fertilizer = substr(sv$barcode, 8, 8)
#' sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
#'            ifelse(sv$fertilizer == "B", "50", "0"))
#' sv<-bw.time(sv, plantingDelay = 0,
#'       phenotype="area.pixels", cutoff=10, timeCol="timestamp",
#'       group=c("barcode", "rotation"), plot=T)
#' phenotypes <- c('area.pixels', 'convex_hull_area.pixels',
#'             'convex_hull_vertices', 'ellipse_angle.degrees',
#'             'ellipse_eccentricity', 'ellipse_major_axis.pixels',
#'             'ellipse_minor_axis.pixels', 'height.pixels',
#'             'hue_circular_mean.degrees', 'hue_circular_std.degrees',
#'             'hue_median.degrees', 'longest_path.pixels', 'perimeter.pixels',
#'             'solidity', 'width.pixels')
#' phenoForm<-paste0("cbind(", paste0(phenotypes, collapse=", "), ")")
#' groupForm<-"DAS+DAP+barcode+genotype+fertilizer"
#' sv<-aggregate(as.formula(paste0(phenoForm, "~", groupForm)), data=sv, mean, na.rm=T)
#' sv<-bw.outliers(sv, phenotype="area.pixels",
#'       group = c("DAS", "genotype", "fertilizer"), plotgroup = c("barcode"))
#' water<-bw.water("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/metadata.json")
#' 
#' df<-sv
#' w<-water
#' pheno="area.pixels"
#' id =c("barcode")
#' time="DAS"
#' 
#' x<-pwue(df, w, pheno, time, id)
#' library(ggplot2)
#' ggplot(x, aes(x=DAS, y=pWUE, group = barcode))+
#'   geom_line()+
#'   theme_minimal()
#'  
#' head(x[which(x$pWUE< -100),])
#' head(x[x$pWUE < -100,]) # note the counter-intuitive handling of NAs in the index column.
#' # Either use data.table or wrap conditions in which().
#' 
#' ## End(Not run)
#' 
#' 
#' @export

pwue<-function(df, w, pheno="area.pixels", time="DAS", id="barcode"){
  #* `format watering data`
  w<-w[, !grepl("^timestamp$|^local_time$",colnames(w))]
  w$snapshot_sorter <- as.numeric(sub("snapshot", "", w$snapshot))
  #w<-data.table::setorderv(w, cols = c(id,time, 'snapshot_sorter')) # arrange plants by time
  w<-w[!w$weight_before<0,] # remove errors
  #* `format phenotype data`
  df<-df[, !grepl("^timestamp$|^local_time$",colnames(df))]
  #* `join datasets`
  x<-data.table::as.data.table(plyr::join(df, w, by=intersect(colnames(df), colnames(w)), type="left", match="all"))
  x<-data.table::setorderv(x, cols = c(id,time, "snapshot_sorter")) # x can have duplicate pheno rows if w has >1 watering per day.
  #* `calculate delta values and pwue`
  x_deltas<-data.table::rbindlist(lapply(split(x, by=id), function(d){
    d<-data.table::setorderv(d, cols = c(time))
      # could lead the weight before or lag the weight after
    #* `water transpired/lost`
    d$water_used_between_waterings <- d$weight_after - data.table::shift(d$weight_before, n=1, type="lead")
    #* `sum if there are multiple waterings per day`
    if( any(duplicated(d[[time]])) ){
      d<-data.table::rbindlist(lapply(unique(d[[time]]), function(tm){
        sub<-d[d[[time]]==tm, ]
        sub<-data.table::setorderv(sub, cols = c("snapshot_sorter"))
        out<-sub[1,]
        if(all(is.na(sub$water_used_between_waterings))){
          out$water_used_between_waterings <- NA
        } else{
          out$water_used_between_waterings <- sum(sub$water_used_between_waterings, na.rm=T) 
        }
        out[[pheno]] <- mean(sub[[pheno]], na.rm=T)
        out
      }))
    }
    #* `calculate delta phenotype`
    d[[paste0("delta_",pheno)]] <- pmax(d[[pheno]] - data.table::shift(d[[pheno]], n=1, type="lag"), 0, na.rm=F)
    #* `calculate pseudo-WUE`
    d$pWUE <- d[[paste0("delta_",pheno)]] / d$water_used_between_waterings
    d
  }))
  data.table::setDF(x_deltas)
  return(x_deltas)
}

