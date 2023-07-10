#' subset helper function for use reading in large data, called in pcv.sub.read
#' 
#' 
#' 
#' @param inputFile Path to csv file of plantCV output, should be provided internally in read.pcv
#' @param filters filtering conditions, see read.pcv for details. Format as list("trait in area, perimeter", "other contains stringToMatch")
#' @param awk Optional awk command to use instead.
#' @keywords read.csv, pcv, wide, long
#' @details awkHelper attempts to make awk commands from human readable input. Currently when filters are supplied the input file has quotes removed by `sed` then is piped into awk, so an equivalent command line statement may be: sed 's/\"//g' pcvrTest2.csv | awk -F ','  '{ if (NR==1 || $18=="area") { print } }'
#' @return Returns a character string representing a unix style awk statement which is typically passed to \code{pipe} or used as a connection in \code{data.table::fread}.
#' @examples 
#' ## Not run: 
#' inputFile = "localCSVfile.csv"
#' filters = "trait contains area, perimeter"
#' cat(awkHelper(inputFile, filters)) # `sed 's/"//g' /home/jsumner/Desktop/stargate/fahlgren_lab/pcvrTestData/pcvrTest1.csv |  awk -F  ','  '{ if ( ( ($18 ~ /area|perimeter/) ) ) { print } }'`
#' filters = list("trait in area, width, height")
#' cat(awkHelper(inputFile, filters)) # `sed 's/"//g' /home/jsumner/Desktop/stargate/fahlgren_lab/pcvrTestData/pcvrTest1.csv |  awk -F  ','  '{ if ( ( ($18=="area") || ($18=="width") || ($18=="height") ) ) { print } }'`
#' ## End(Not run)
#' @export
awkHelper<-function(inputFile, filters, awk=NULL){
  if(is.null(awk)){
    if(!is.list(filters)){filters<-list(filters)}
    sed<- paste0("sed 's/\"//g' ", inputFile, " | ")
    awkStart<-"awk -F "
    awkDelim <-"',' "
    awkFiltStart<-"'{ if ("
    awkFiltEnd<-") { print } }'"
    COLS = colnames(read.csv(inputFile, nrows=1))
    awkFilts<-lapply(filters, function(filt){
      filtCol = strsplit(filt," ")[[1]][1]
      affector<-strsplit(filt," ")[[1]][2]
      values<-trimws(gsub(",", " ", strsplit(filt," ")[[1]][-c(1:2)]))
      if(affector %in% c("in", "is", "=")){
        paste(paste0("($", which(COLS == filtCol), '=="', values,'")' ), collapse=" || ")
      } else if(affector=="contains"){
        valReg<-paste0(values, collapse="|")
        paste0("($", which(COLS==filtCol), " ~ /", valReg, "/)")
      }
    })
    awkFilt = paste(paste("(",awkFilts,")"), collapse = " && ")
    awkCommand<-capture.output(cat(sed, awkStart, awkDelim, awkFiltStart, awkFilt, awkFiltEnd))
  } else {awkCommand = awk}
  return(awkCommand)
}
