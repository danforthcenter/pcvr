#' reader function called within read.pcv when large data is used
#' 
#' @param inputFile Path to csv file of plantCV output, should be provided internally in read.pcv
#' @param filters filtering conditions, see read.pcv for details. Format as list("trait in area, perimeter", "other in value")
#' @param awk Optional awk command to use instead.
#' @import data.table
#' @import vroom
#' @keywords read.csv, pcv, wide, long
#' @export
#' 
pcv.sub.read<-function(inputFile, filters, reader = "read.csv", awk=NULL, ...){
  awkCommand<-awkHelper(inputFile, filters, awk)
  COLS = colnames(read.csv(inputFile, nrows=1))
  if(reader=="vroom"){
    x<-as.data.frame(vroom::vroom(pipe(awkCommand), show_col_types = FALSE, delim = ",", col_names = COLS, ...))
  } else if (reader== "fread"){
    x<-as.data.frame(data.table::fread(cmd=awkCommand, col.names = COLS, ...))
  } else {
    readingFunction <- match.fun(reader)
    x<-suppressMessages(as.data.frame(readingFunction( pipe(awkCommand), ...)))
    colnames(x)<-COLS
  }
  return(x)
}
