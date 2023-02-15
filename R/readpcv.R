#' Read in plantCV csv output in wide or long format
#' 
#' @param filepath Path to csv file of plantCV output.
#' @param mode One of "wide" or "long", partial string matching is supported. Priority here is given to reducing dependencies, so long to wide transformation is done using lapply. This can cause problems if there are missing values in `traitCol`
#' @param singleValueOnly Logical, should only single value traits be returned? If `mode` = "wide" then this should almost always be set to T.
#' @param traitCol Column with phenotype names, defaults to "trait". This should generally not need to be changed from the default.
#' @param labelCol Column with phenotype labels (units), defaults to "label". This should generally not need to be changed from the default. This is used with `traitCol` when `mode`="wide" to identify unique traits since some may be ambiguous (ellipseCenter.x vs ellipseCenter.y, etc)
#' @param valueCol Column with phenotype values, defaults to "value". This should generally not need to be changed from the default.
#' @param multiValPattern If `singleValueOnly`=TRUE then this is used to identify multi value traits. By default this is "hist|frequencies". If this argument has length of 1 then it is taken as either a single phenotype or a regex pattern to find values of `trait` that are multi-value phenotypes. Alternatively this can be a vector of phenotype names to remove (see examples).
#' @param reader The function to use to read in data, defaults to "read.csv". Other useful options are "vroom" and "fread", from the vroom and data.table packages, respectively. With files that are still very large after subsetting "fread" or "vroom" should be used.
#' @param filters If a very large pcv output file is read then it may be desireable to subset it before reading it into R, either for ease of use or because of RAM limitations. The filter argument works with "COLUMN in VALUES" syntax. This can either be a character vector or a list of character vectors. In these vectors there needs to be a column name, one of " in ", " is ", or " = ", then a set of comma delimited values to filter that column for (see examples). Note that this and `awk` both use awk through pipe(). This functionality will not work on a windows system. 
#' @param awk As an alternative to `filters` a direct call to awk can be supplied here, in which case that call will be used through pipe().
#' @param ... Other arguments passed to the reader function. In the case of 'vroom' and 'fread' there are several defaults provided already which can be overwritten with these extra arguments.
#' @keywords read.csv, pcv, wide, long
#' @import data.table
#' @import vroom

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

#' @export
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

#' @export
read.pcv<-function(filepath, mode="wide", singleValueOnly=T,
                   traitCol="trait", labelCol="label", valueCol="value",
                   multiValPattern = "hist|frequencies", reader="read.csv", filters=NULL, awk=NULL, ...){
  
  if(is.null(filters) & is.null(awk)){
    readingFunction<-match.fun(reader)
    df1<-readingFunction(filepath, ...)
  } else{ 
    df1<-pcv.sub.read(inputFile=filepath, filters=filters, reader = reader, awk=awk, ...)  
    }

  if(singleValueOnly){
    if(length(multiValPattern)==1){ df1<-df1[!grepl(multiValPattern, df1[[traitCol]]), ]
    } else { df1<-df1[!df1[[traitCol]] %in% multiValPattern, ] }
    }
  if(match.arg(mode, c("wide","long"))=="wide" ){
    long<-df1
    long[[traitCol]]<-ifelse(long[[labelCol]]=="none", long[[traitCol]],
                             paste(long[[traitCol]], long[[labelCol]], sep='.'))
    wide<-long[long[[traitCol]]==unique(long[[traitCol]])[1],]
    nc<-ncol(wide)
    wide[,seq(nc+1, nc+length(unique(long[[traitCol]])),1) ]<-lapply(unique(long[[traitCol]]), function(i) long[long[[traitCol]]==i, valueCol])
    colnames(wide)<-c(colnames(wide)[1:nc], unique(long[[traitCol]]))
    out<-wide
  } else{out<-df1
  if(!is.null(traitCol)){
    if(traitCol %in% colnames(out)){
      out[[traitCol]]<-gsub("/", ".over.", out[[traitCol]])
      out[[traitCol]]<-gsub("\\'", "", out[[traitCol]]) 
    }
  }
  }
  colnames(out)<-gsub("/", ".over.", colnames(out))
  colnames(out)<-gsub("\\'", "", colnames(out))
  return(out)
}

#' @examples 
#' 
#' file = "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest1.csv"
#' df1<-read.pcv(file, "wide", T, multiValPattern = "hist|frequencies")
#' df1b<-read.pcv(file, "wide", T, multiValPattern = c("index_frequencies_index_ari", "index_frequencies_index_ci_rededge", "index_frequencies_index_ndvi", "npq_hist_NPQ", "yii_hist_Fq'/Fm'", "yii_hist_Fv/Fm"))
#' identical(df1, df1b)
#' df2<-read.pcv(file, "long", T)
#' dim(df2)
#' # Note only data stored on a Unix style system can be subset before reading in.
#' # For DDPSC employees there are larger datasets on stargate that better show the benefit of subsetting before reading data in.
#' fileBig = "/shares/mgehan_share/kmurphy/maize_2022/bellwether/results_vis_SV/07252022_VIS_SV_MG001_results.csv"
#' start<-Sys.time()
#' x3a<-pcv.sub.read(inputFile=fileBig, reader = "vroom",
#'                   filters = list("trait in area, perimeter", "barcode is Ea008AA114352"))
#' Sys.time()-start
#' start<-Sys.time()
#' x3b<-pcv.sub.read(inputFile=fileBig, reader = "fread",
#'                   filters = list("trait in area, perimeter", "barcode is Ea008AA114352"))
#' Sys.time()-start
#' start<-Sys.time()
#' x3c<-pcv.sub.read(inputFile=fileBig, reader = "read.csv",
#'                   filters = list("trait in area, perimeter", "barcode is Ea008AA114352"))
#' Sys.time()-start
#' dim(x3a)
#' dim(x3b)
#' dim(x3c)
