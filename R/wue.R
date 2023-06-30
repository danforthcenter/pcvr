#' Calculate pseudo water use efficiency from phenotype and watering data
#' 
#' @description Water use efficiency (WUE) is the change in biomass per unit of water metabolized.
#' Using image based phenotypes and watering data we can calculate pseudo-WUE (pwue) over time.
#' Here area is used as a proxy for biomass and transpiration is approximated using watering data.
#' The equation is then \eqn{\frac{P_[t] - P_[t-1]}{W_[t]-W_[t-1] }}{P_[t] - P_[t-1] / W_[t]-W_[t-1]},
#' where P is the phenotype and W is the weight before watering.
#' 
#' @param df Dataframe containing wide single-value phenotype data. This should already be aggregated to one row per plant per day (angles/rotations combined).
#' @param w Watering data as returned from bw.water.
#' @param pheno Phenotype column name, defaults to "area.pixels"
#' @param timeGroup Variable(s) that identify a plant on a given day. Defaults to \code{c("barcode", "DAS")}.
#' @param id Variable(s) that identify a plant over time. Defaults to \code{"barcode"}.
#' @keywords read.csv, pcv, wide, long
#' @import jsonlite
#' @return A data frame containing the bellwether watering data joined to phenotype data with new columns for change in the phenotype, change in the pre-watering weight, and pseudo-water use efficiency (pWUE).
#' @examples
#' 
#' ## Not run:
#' 
#' library(data.table)
#' sv<-read.pcv("https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv", mode="wide", singleValueOnly = T, reader="fread")
#' sv$genotype = substr(sv$barcode, 3,5)
#' sv$genotype = ifelse(sv$genotype == "002", "B73", ifelse(sv$genotype == "003", "W605S", ifelse(sv$genotype == "004", "MM", "Mo17")))
#' sv$fertilizer = substr(sv$barcode, 8, 8)
#' sv$fertilizer = ifelse(sv$fertilizer == "A", "100", ifelse(sv$fertilizer == "B", "50", "0"))
#' sv<-bw.time(sv, plantingDelay = 0, phenotype="area.pixels", cutoff=10, timeCol="timestamp", group=c("barcode", "rotation"), plot=T)
#' phenotypes <- c('area.pixels', 'convex_hull_area.pixels', 'convex_hull_vertices', 'ellipse_angle.degrees', 'ellipse_eccentricity', 'ellipse_major_axis.pixels', 'ellipse_minor_axis.pixels', 'height.pixels', 'hue_circular_mean.degrees', 'hue_circular_std.degrees', 'hue_median.degrees', 'longest_path.pixels', 'perimeter.pixels', 'solidity', 'width.pixels')
#' phenoForm<-paste0("cbind(", paste0(phenotypes, collapse=", "), ")")
#' groupForm<-"DAS+DAP+barcode+genotype+fertilizer"
#' sv<-aggregate(as.formula(paste0(phenoForm, "~", groupForm)), data=sv, mean, na.rm=T)
#' sv<-bw.outliers(sv, phenotype="area.pixels", group = c("DAS", "genotype", "fertilizer"), plotgroup = c("barcode", "rotation"))
#' water<-bw.water("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/metadata.json")
#' 
#' df<-sv
#' w<-water
#' pheno="area.pixels"
#' timeGroup = c("barcode", "DAS")
#' id =c("barcode")
#' 
#' x<-pwue(sv, water, "area.pixels", c("barcode", "DAS"), c("barcode"))
#' library(ggplot2)
#' ggplot(x, aes(x=DAS, y=pWUE, group = barcode))+
#'   geom_line()+
#'   theme_minimal()
#' 
#' ## End(Not run)
#' 
#' 
#' @export

pwue<-function(df, w, pheno="area.pixels", timeGroup=c("barcode", "DAS"), id="barcode"){
  #* `format watering data`
  w<-data.table::as.data.table(w[, !grepl("^timestamp$|^local_time$",colnames(w))])
  w<-data.table::setorderv(w, cols = timeGroup) # arrange plants by time
  w<-w[!w$weight_before<0,] # remove errors
  
  #* `group by barcode+DAS to pseudo aggregate waterings`
  if(max(as.numeric(table(w[,..timeGroup]))) > 1){ # if multiple waterings then aggregate
    ags<-lapply(split(w, by=timeGroup), function(wd){
      const_cols<-names(which(unlist(lapply(wd, function(i) {length(unique(i))} ))==1)) # take columns that are constant in subset
      out<-wd[1, ..const_cols]
      out$water_amount <- sum(wd$water_amount, na.rm=T)
      out$weight_before <- min(wd$weight_before, na.rm=T)
      out$weight_after <- max(wd$weight_after, na.rm=T)
      out
    })
    const_cols<-Reduce(intersect, lapply(ags,colnames)) # column names could be wrong if n=1 for anything
    w <- do.call(rbind, lapply(ags, function(wd){ wd[, ..const_cols]}))
  }
  #* `format phenotype data`
  df<-data.table::as.data.table(df[, !grepl("^timestamp$|^local_time$",colnames(w))])
  
  #* `join data`
  x<-data.table::as.data.table(plyr::join(df, w, by=intersect(colnames(df), colnames(w))))
  x<-data.table::setorderv(x, cols = timeGroup)
  
  #* how to handle 0s is not obvious to me yet. They'll make Inf once I do the division, so maybe just make them
  #* a very small number?
  #* x$water_amount <- ifelse(x$water_amount <=0, 0.1, x$water_amount)
  
  x<-do.call(rbind, lapply(split(x, by=id), function(d){
    d$deltaWeight_before <- d$weight_before - data.table::shift(d$weight_before, n=1, type="lag")
    d[[paste0("delta_",pheno)]] <- d[[pheno]] - data.table::shift(d[[pheno]], n=1, type="lag")
    #* by water amount this is delta_pheno / water_previous
    d$pWUE <- d[[paste0("delta_",pheno)]] / d$deltaWeight_before
    d
  }))
  
  return(x)
}

