#' Time conversion and plotting for bellwether data
#' 
#' @param df Data frame to use, this can be in wide or long format as specified by the wide argument.
#' @param mode One of "DAS", "DAP" or "DAE" (Days After Planting and Days After Emergence).
#' Defaults to NULL in which case all columns are added.
#' Note that if timeCol is not an integer then DAS is always returned.
#' @param plantingDelay If `mode` includes "DAP" then `plantingDelay` is used to adjust "DAS"
#' @param phenotype If `mode` includes "DAE" then this is the phenotype used to classify emergence. 
#' @param cutoff If `mode` inlcludes "DAE" then this value is used to classify emergence.
#' Defaults to 1, meaning an image with a value of 1 or more for `phenotype` has "emerged".
#' @param timeCol Column of input time values, defaults to "timestamp".
#' If this is not an integer then it is assumed to be
#' a timestamp in the format of the format argument.
#' @param group  Grouping variables to specify unique plants as a character vector.
#' This defaults to "Barcodes". These taken together should identify a unique plant across time,
#' although often "angle" or "rotation" should be added.
#' @param plot Logical, should plots of the new time variables be printed?
#' @param wide Logical, is data in wide format? Defaults to TRUE.
#' @param format An R POSIXct format, defaults to lemnatech standard format.
#' This is only used if timeCol is not an integer.
#' @param traitCol Column with phenotype names, defaults to "trait".
#' This should generally not need to be changed from the default.
#' @param valueCol Column with phenotype values, defaults to "value".
#' This should generally not need to be changed from the default.
#' @keywords Bellwether, ggplot
#' @import ggplot2
#' @return The input dataframe with new integer columns for different ways
#' of describing time in the experiment.
#' @export
#' @examples 
#' 
#' ## Not run:
#' 
#' sv<-read.pcv(
#' "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv",
#'  mode="wide", singleValueOnly =TRUE, reader="fread")
#' sv$genotype = substr(sv$barcode, 3,5)
#' sv$genotype = ifelse(sv$genotype == "002", "B73",
#'                      ifelse(sv$genotype == "003", "W605S",
#'                            ifelse(sv$genotype == "004", "MM", "Mo17")))
#' sv$fertilizer = substr(sv$barcode, 8, 8)
#' sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
#'                    ifelse(sv$fertilizer == "B", "50", "0"))
#' sv<-bw.time(sv, plantingDelay = 0, phenotype="area.pixels", cutoff=10,
#'  timeCol="timestamp", group=c("barcode", "rotation"), plot=FALSE)
#' 
#' 
#' svl<-read.pcv(
#' "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv",
#'  mode="long", singleValueOnly =TRUE, reader="fread")
#' svl$genotype = substr(svl$barcode, 3,5)
#' svl$genotype = ifelse(svl$genotype == "002", "B73",
#'                      ifelse(svl$genotype == "003", "W605S",
#'                            ifelse(svl$genotype == "004", "MM", "Mo17")))
#' svl$fertilizer = substr(svl$barcode, 8, 8)
#' svl$fertilizer = ifelse(svl$fertilizer == "A", "100",
#'                    ifelse(svl$fertilizer == "B", "50", "0"))
#' svl<-bw.time(svl, plantingDelay = 0, phenotype="area", cutoff=10, timeCol="timestamp",
#'  group=c("barcode", "rotation"), plot=FALSE,wide=FALSE)
#' 
#' ## End(Not run)
#'
bw.time<-function(df = NULL, mode=NULL, plantingDelay = NULL,
                  phenotype=NULL, cutoff=1, timeCol="timestamp",
                  group="Barcodes", plot=TRUE, wide=TRUE, format="%Y-%m-%d %H:%M:%S", traitCol="trait", valueCol="value"){
  if(is.null(mode) || !mode %in% c("DAS", "DAP", "DAE")){ mode=c("DAP", "DAE") }
  if(is.null(plantingDelay) & "DAP" %in% mode){mode<-mode[-which(mode=="DAP")]}
  if(is.null(phenotype) & "DAE" %in% mode){mode<-mode[-which(mode=="DAE")]}
  
  if(!is.integer(df[[timeCol]])){
    df[[timeCol]]<-as.POSIXct(strptime(df[[timeCol]],format = format))
    beg <- min(df[[timeCol]], na.rm=TRUE)
    df$DAS <- floor(as.numeric((df[[timeCol]] - beg)/60/60/24))
    timeCol="DAS"
  }
  if("DAP" %in% mode){
    df$DAP = df[[timeCol]] + plantingDelay 
  }
  if("DAE" %in% mode & wide){
    DAE_SPLIT <- interaction(df[, group])
    df<-do.call(rbind, lapply( split(df, DAE_SPLIT), function(d){
      subd<-d[d[[phenotype]] >= cutoff & !is.na(d[[phenotype]]) , ]
      if(nrow(subd)==0){subd<-data.frame(DAS=max(df[[timeCol]])+1) ; colnames(subd)<-timeCol} # if all NA area then remove all rows
      d$DAE<-d[[timeCol]] - min(subd[[timeCol]], na.rm=TRUE)
      d
    }))
  } else if ("DAE" %in% mode & !wide){
    DAE_SPLIT <- interaction(df[, group])
    df<-do.call(rbind, lapply( split(df, DAE_SPLIT), function(d){
      subd<-d[ d[[traitCol]]==phenotype & d[[valueCol]]>=cutoff & !is.na(d[[valueCol]])  , ]
      if(nrow(subd)==0){subd<-data.frame(DAS=max(df[[timeCol]])+1) ; colnames(subd)<-timeCol}
      d$DAE<-d[[timeCol]] - min(subd[[timeCol]], na.rm=TRUE)
      d
    }))
  }
  rownames(df)<-NULL
  if(plot & !is.null(phenotype) & wide){
    plotDat<-df
    plotDat$plotGroup<-interaction(plotDat[,c(group)])
    for(m in mode){
      p<-ggplot2::ggplot(plotDat, ggplot2::aes(x=.data[[m]],y=.data[[phenotype]], group=.data$plotGroup))+
        ggplot2::geom_line()+
        ggplot2::labs(x=m, y=phenotype, title = m)+
        pcv_theme()
      print(p)
    }
  } else if(plot & !is.null(phenotype) & !wide){
    #* plot long data
    plotDat<-df[df[[traitCol]]==phenotype, ]
    plotDat$plotGroup<-interaction(plotDat[,c(group)])
    for(m in mode){
      p<-ggplot2::ggplot(plotDat, ggplot2::aes(x=.data[[m]],y=.data[[valueCol]], group=.data$plotGroup))+
        ggplot2::geom_line()+
        ggplot2::labs(x=m, y=phenotype, title = m)+
        pcv_theme()
      print(p)
    }
  }
  return(df)
}
