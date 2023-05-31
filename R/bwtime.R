#' Time conversion and plotting for bellwether data
#' 
#' @param df Data frame to use. Should be in the format of output from read.pcv.bw
#' @param mode One of "DAS", "DAP" or "DAE" (Days After Planting and Days After Emergence). Defaults to NULL in which case all columns are added. Note that if timeCol is not an integer then DAS is always returned.
#' @param plantingDelay If `mode` includes "DAP" then `plantingDelay` is used to adjust "DAS"
#' @param phenotype If `mode` includes "DAE" then this is the phenotype used to classify emergence. 
#' @param cutoff If `mode` inlcludes "DAE" then this value is used to classify emergence. Defaults to 1, meaning an image with a value of 1 or more for `phenotype` has "emerged".
#' @param timeCol Column of input time values, defaults to "timestamp". If this is not an integer then it is assumed to be a timestamp in the format of the format argument.
#' @param group  Grouping variables to specify unique plants as a character vector. This defaults to "Barcodes". These taken together should identify a unique plant across time, although often "angle" or "rotation" should be added.
#' @param plot Logical, should plots of the new time variables be printed?
#' @param format An R POSIXct format, defaults to lemnatech standard format. This is only used if timeCol is not an integer.
#' @keywords Bellwether, ggplot
#' @import ggplot2
#' @return The input dataframe with new integer columns for different ways of describing time in the experiment.
#' @export
#' @examples 
#' bw<-read.pcv.bw( file="https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/bwTestPhenos.csv", snapshotFile="https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/bwTestSnapshot.csv", designFile="https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/bwTestDesign.csv",metaCol="meta",metaForm="vis_view_angle_zoom_horizontal_gain_exposure_v_new_n_rep",joinSnapshot="id",conversions = list(area=13.2*3.7/46856) )
#' bw<-bw.time(df = bw, mode=NULL, plantingDelay = 4, phenotype="area_adj", timeCol="DAS", group="Barcodes", plot=T)
#'
bw.time<-function(df = NULL, mode=NULL, plantingDelay = NULL,
                  phenotype=NULL, cutoff=1, timeCol="timestamp",
                  group="Barcodes", plot=T, format="%Y-%m-%d %H:%M:%S" ){
  if(is.null(mode) || !mode %in% c("DAS", "DAP", "DAE")){ mode=c("DAP", "DAE") }
  if(is.null(plantingDelay) & "DAP" %in% mode){mode<-mode[-which(mode=="DAP")]}
  if(is.null(phenotype) & "DAE" %in% mode){mode<-mode[-which(mode=="DAE")]}
  
  if(!is.integer(df[[timeCol]])){
    df[[timeCol]]<-as.POSIXct(strptime(df[[timeCol]],format = format))
    beg <- min(df[[timeCol]], na.rm=T)
    df$DAS <- floor(as.numeric((df[[timeCol]] - beg)/60/60/24))
    timeCol="DAS"
  }
  if("DAP" %in% mode){
    df$DAP = df[[timeCol]] + plantingDelay 
  }
  if("DAE" %in% mode){
    df<-do.call(rbind, lapply( split(df, as.formula(paste0("~",paste( group,collapse="+" )))), function(d){
      subd<-d[d[[phenotype]] >= cutoff & !is.na(d[[phenotype]]) , ]
      if(nrow(subd)==0){subd<-data.frame(DAS=max(df[[timeCol]])+1) ; colnames(subd)<-timeCol} # if all NA area then remove all rows
      d$DAE<-d[[timeCol]] - min(subd[[timeCol]], na.rm=T)
      d#[d$DAE>=0, ]
    }))
  }
  rownames(df)<-NULL
  if(plot & !is.null(phenotype)){
    plotDat<-df
    plotDat$plotGroup<-interaction(plotDat[,c(group)])
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
