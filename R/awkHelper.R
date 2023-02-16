#' subset helper function for use reading in large data, called in pcv.sub.read
#' 
#' @param inputFile Path to csv file of plantCV output, should be provided internally in read.pcv
#' @param filters filtering conditions, see read.pcv for details. Format as list("trait in area, perimeter", "other in value")
#' @param awk Optional awk command to use instead.
#' @keywords read.csv, pcv, wide, long
#' 
#' 
awkHelper<-function(inputFile, filters, awk=NULL){
  if(is.null(awk)){
    if(!is.list(filters)){filters<-list(filters)}
    awkStart<-"awk -F "
    awkDelim <-"',' "
    awkFiltStart<-"'{ if ("
    awkFiltEnd<-") { print } }' "
    COLS = colnames(read.csv(inputFile, nrows=1))
    awkFilts<-lapply(filters, function(filt){
      filtCol = strsplit(filt," ")[[1]][1]
      filt<-gsub("( = )|( is )", " in ", filt)
      values = trimws(strsplit(trimws(strsplit(filt,"in")[[1]][-1]),",")[[1]])
      paste(paste0("($", which(COLS == filtCol), '=="', values,'")' ), collapse=" || ")
    })
    awkFilt = paste(paste("(",awkFilts,")"), collapse = " && ")
    awkCommand<-capture.output(cat(awkStart, awkDelim, awkFiltStart, awkFilt, awkFiltEnd, inputFile))
  } else {awkCommand = awk}
  return(awkCommand)
}