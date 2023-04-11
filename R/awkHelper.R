#' subset helper function for use reading in large data, called in pcv.sub.read
#' 
#' Idea with awkHelperRegex branch is to add a new option for the in/is/= syntax, probably "contains"
#' which would be used like "trait contains frequency" ~ grepl(".*?frequency.*?")
#' Since a part of awkHelper only involves reading in the first line of the file it would
#' not work to say "find all the things with this string", then parse them into || statements for awk.
#' Instead I should use the regex options in awk through `\string\`
#' 
#' 
#' 
#' @param inputFile Path to csv file of plantCV output, should be provided internally in read.pcv
#' @param filters filtering conditions, see read.pcv for details. Format as list("trait in area, perimeter", "other in value")
#' @param awk Optional awk command to use instead.
#' @keywords read.csv, pcv, wide, long
#' @details awkHelper attempts to make awk commands from human readable input. Currently when filters are supplied the input file has quotes removed by `sed` then is piped into awk, so an equivalent command line statement may be: sed 's/\"//g' pcvrTest2.csv | awk -F ','  '{ if (NR==1 || $18=="area") { print } }'
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
      filt<-gsub("( = )|( is )", " in ", filt)
      values = trimws(strsplit(trimws(strsplit(filt," in ")[[1]][-1]),",")[[1]])
      paste(paste0("($", which(COLS == filtCol), '=="', values,'")' ), collapse=" || ")
    })
    awkFilt = paste(paste("(",awkFilts,")"), collapse = " && ")
    awkCommand<-capture.output(cat(sed, awkStart, awkDelim, awkFiltStart, awkFilt, awkFiltEnd))
  } else {awkCommand = awk}
  return(awkCommand)
}
# sed 's/"//g' /home/jsumner/Desktop/stargate/fahlgren_lab/pcvrTestData/pcvrTest1.csv |  awk -F  ','  '{ if ($18 ~ /area|perimeter/)  { print } }'
# sed 's/"//g' /home/jsumner/Desktop/stargate/fahlgren_lab/pcvrTestData/pcvrTest1.csv |  awk -F  ','  '{ if ($18 ~ /^area$|^perimeter$/)  { print } }'|head