#' Time conversion and plotting for bellwether data
#' 
#' @param df Data frame to use. Should be in the format of output from read.pcv.bw
#' @param mode One of "DAP" or "DAE" (Days After Planting and Days After Emergence). Defaults to NULL in which case both columns are added. 
#' @param plantingDelay If `mode` includes "DAP" then `plantingDelay` is used to adjust "DAS"
#' @param phenotype If `mode` includes "DAE" then this is the phenotype used to classify emergence. 
#' @param cutoff If `mode` inlcludes "DAE" then this value is used to classify emergence. Defaults to 1, meaning an image with a value of 1 or more for `phenotype` has "emerged".
#' @param timeCol Column of input time values, defaults to "DAS"
#' @param group  Grouping variables to specify unique plants as a character vector. This defaults to "Barcodes". These taken together should identify a unique plant across time.
#' @param plot Logical, should plots of the new time variables be printed?
#' @keywords Bellwether, ggplot
#' @import ggplot2
#' @export
#' @examples 
#' 
#'
bw.time<-function(df = NULL, mode=NULL, plantingDelay = 4,
                  phenotype=NULL, cutoff=1, timeCol="DAS",
                  group="Barcodes", plot=T ){
  if(is.null(mode) || !mode %in% c("DAP", "DAE")){ mode=c("DAP", "DAE") }
  if("DAP" %in% mode){
    df$DAP = df[[timeCol]] + plantingDelay 
  }
  if("DAE" %in% mode){
    df<-do.call(rbind, lapply( split(df, as.formula(paste0("~",paste( group,collapse="+" )))), function(d){
      subd<-d[d[[phenotype]] >= 1 & !is.na(d[[phenotype]]) , ]
      if(nrow(subd)==0){subd<-data.frame(DAS=max(d[[timeCol]])+1) ; colnames(subd)<-timeCol} # if all NA area then remove all rows
      d$DAE<-d[[timeCol]] - min(subd[[timeCol]], na.rm=T)
      d#[d$DAE>=0, ]
    }))
  }
  rownames(df)<-NULL
  if(plot & !is.null(phenotype)){
    plotGroup<-c(group, "angle")
    plotDat<-df
    plotDat$plotGroup<-interaction(plotDat[,c(group,"angle")])
    p<-ggplot2::ggplot(plotDat, ggplot2::aes(x=.data[[timeCol]], y=.data[[phenotype]], group=.data$plotGroup))+
      ggplot2::geom_line()+
      pcv_theme()
    print(p)
    for(m in mode){
      p<-ggplot2::ggplot(plotDat, ggplot2::aes(x=.data[[m]],y=.data[[phenotype]], group=.data$plotGroup))+
        ggplot2::geom_line()+
        ggplot2::labs(x=m, y=phenotype, title = m)+
        pcv_theme()
      print(p)
    }
  }
  return(df)
}
