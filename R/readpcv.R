#' Read in plantCV csv output in wide or long format
#' 
#' @param filepath Path to csv file of plantCV output.
#' @param mode NULL (the default) or one of "wide" or "long", partial string matching is supported.
#'    This controls whether data is \strong{returned} in long or wide format. If left NULL then
#'    the output format will be the same as the input format.
#' @param traitCol Column with phenotype names, defaults to "trait".
#'   This should generally not need to be changed from the default. This, 
#'   labelCol, and valueCol are used to determine if data are in long format in their
#'   raw state (the csv file itself).
#' @param labelCol Column with phenotype labels (units), defaults to "label".
#'   This should generally not need to be changed from the default.
#'   This is used with `traitCol` when `mode`="wide" to identify
#'   unique traits since some may be ambiguous
#' (ellipseCenter.x vs ellipseCenter.y, bins of histograms, etc)
#' @param valueCol Column with phenotype values, defaults to "value".
#'   This should generally not need to be changed from the default.
#' @param reader The function to use to read in data,
#'   defaults to NULL in which case `data.table::fread` is used if filters are in place
#'   and `read.csv` is used otherwise.
#'   Other useful options are "vroom" and "fread", from the vroom and data.table packages, respectively.
#'   With files that are still very large after subsetting "fread" or "vroom" should be used.
#'   Note that if you use `read.csv` with filters in place then you will need to specify `header=FALSE`
#'   so that the piped output from awk is read correctly.
#' @param filters If a very large pcv output file is read then it may be desireable
#'   to subset it before reading it into R, either for ease of use or because of RAM limitations.
#'   The filter argument works with "COLUMN in VALUES" syntax. This can either be a character vector
#'   or a list of character vectors. In these vectors there needs to be a column name,
#'   one of " in ", " is ", or " = " to match the string exactly, or "contains"
#'   to match with awk style regex, then a set of comma delimited values to filter
#'   that column for (see examples). Note that this and `awk` both use awk through pipe().
#'   This functionality will not work on a windows system. 
#' @param awk As an alternative to `filters` a direct call to awk can be supplied here,
#'   in which case that call will be used through pipe().
#' @param ... Other arguments passed to the reader function.
#'   In the cases of 'vroom' and 'fread' there are several defaults provided already
#'   which can be overwritten with these extra arguments.
#'   
#' @details
#' In plantCV version 4 the single value traits are returned in wide format from \code{json2csv}
#' and the multi value traits are returned in long format. When data is read in using read.pcv 
#' the traitCol, valueCol, and labelCol arguments are checked to determine if the data is in long 
#' format. This is done to keep compatibility with interim versions of plantcv output where all outputs
#' were in a single long format file. 
#' 
#' With the current implementation and plantcv output you can read wide or long format files into
#' wide or long format in R. Keep in mind that the 'mode' argument controls the format that will be returned in R,
#' not the format that the data saved as in your csv file.
#'   
#' @keywords read.csv, pcv, wide, long
#' @return Returns a data.frame in wide or long format.
#' @importFrom stats as.formula
#' @import data.table
#' @examples 
#' 
#' ## Not run: 
#' 
#' if(FALSE){
#' mv = "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/pcv4-multi-value-traits.csv"
#' sv = "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv"
#' 
#' w2w <- read.pcv(sv, mode = "wide", reader="fread")
#' dim(w2w)
#' 
#' w2l <- read.pcv(sv, mode = "long", reader="fread")
#' dim(w2l)
#' 
#' l2w <- read.pcv(mv, mode = "wide", reader="fread")
#' dim(l2w)
#' 
#' l2l <- read.pcv(mv, mode = "long", reader="fread")
#' dim(l2l)
#' 
#' 
#' # Note only data stored on a Unix style system can be subset before reading in.
#' # For DDPSC employees there are larger datasets on stargate that
#' # better show the benefit of subsetting before reading data in.
#' 
#' fileBig="/shares/mgehan_share/llima/Maize_Project_2022/nir_maize_first_exp_results.csv"
#' # library(vroom)
#' start<-Sys.time()
#' x3a<-pcv.sub.read(inputFile=fileBig, reader = "vroom",
#'   filters = list("trait in area, perimeter", "barcode is Ea008AA114352"))
#' Sys.time()-start
#' # library(data.table)
#' start<-Sys.time()
#' x3b<-pcv.sub.read(inputFile=fileBig, reader = "fread",
#'   filters = list("trait in area, perimeter", "barcode is Ea008AA114352"))
#' Sys.time()-start
#' start<-Sys.time()
#' x3c<-pcv.sub.read(inputFile=fileBig, reader = "read.csv",
#'   filters = list("trait in area, perimeter", "barcode is Ea008AA114352"))
#' Sys.time()-start
#' dim(x3a)
#' dim(x3b)
#' dim(x3c)
#' 
#' # There may be situations where you want to use wide mv traits which can read in easily:
#' x4<-read.pcv(fileBig, reader="fread",
#'   filters = list("trait in blue_frequencies"),
#'   mode="wide")
#' }
#' ## End(Not run) 
#' 
#' @export

read.pcv<-function(filepath, mode=NULL, 
                   traitCol="trait", labelCol="label", valueCol="value",
                   reader=NULL, filters=NULL, awk=NULL, ...){
  if(is.null(filters) & is.null(awk)){
    if(is.null(reader)){reader="read.csv"}
    if(reader!="fread"){
      readingFunction<-match.fun(reader)
    }else{readingFunction<-data.table::fread}
    df1<-as.data.frame(readingFunction(filepath, ...))
  } else{
    if(is.null(reader)){reader="fread"}
    df1<-pcv.sub.read(inputFile=filepath, filters=filters, reader = reader, awk=awk, ...)  
    if(nrow(df1)<1){ stop(paste0("0 Rows returned using awk statement:\n", awkHelper(filepath, filters),
                                 "\nMost common issues are misspellings or not including a column name and affector." )) }
    }
  #* `check original data format`
  
  if(all(c(traitCol, valueCol, labelCol) %in% colnames(df1))){
    startsLong = TRUE
  } else if(!any(c(traitCol, valueCol, labelCol) %in% colnames(df1))){
    startsLong = FALSE } else{
      found <- c('traitCol', 'valueCol', 'labelCol')[which(c(traitCol, valueCol, labelCol) %in% colnames(df1))]
      warning(paste0( paste(found, collapse = ", "), " found in column names of data but either all or none of traitCol, valueCol, and labelCol are expected."))
    }
  
  if(is.null(mode)){
    if(startsLong){outputMode = "long"
    } else {outputMode = "wide"}
  } else {
    outputMode <- match.arg(mode, c("wide","long"))
  }
  #* `if data is long and mode is wide`
  
  if(outputMode=="wide" & startsLong ){
    long<-df1
    if(substr(colnames(long)[1],1,1)=="X" & length(unique(long[[1]]))==nrow(long)){long<-long[,-1]}
    long<-long[!is.na(long[[valueCol]]),]
    long[[labelCol]]<-ifelse(is.na(long[[labelCol]]), "none", long[[labelCol]])
    wide<-as.data.frame(data.table::dcast(data.table::as.data.table(long), as.formula(paste0("... ~ ", traitCol, "+", labelCol)), value.var = valueCol, sep="."))
    colnames(wide)<-sub(".none$","",colnames(wide))
    if(any(grepl("hist|frequencies", colnames(wide)))){ # reorder the MV traits by their bins
      #* get a list of the unique non-numeric parts
      histCols <- colnames(wide)[grepl("hist|frequencies", colnames(wide))]
      unique_mvTraits<-unique(gsub("[.]+$", "",gsub("[0-9]+","",histCols)))
      #* for each unique non-numeric part, sort the names
      mvCols_reordered<-unlist(lapply(unique_mvTraits, function(umt){
        iterCols = histCols[grepl(umt, histCols)]
        iterCols_numeric = as.numeric( gsub(paste0(umt, "."), "", iterCols) )
        bins_order<-sort(iterCols_numeric, index.return=TRUE)$ix
        iterCols[bins_order]
      }))
      #* combine the histCols and the other columns, in the new order.
      sv_and_meta_cols<-colnames(wide)[!grepl("hist|frequencies", colnames(wide))]
      wide<-wide[,c(sv_and_meta_cols, mvCols_reordered)]
    }
    out<-wide
    
    #* `if data is long and mode is long`
  } else if(outputMode=="long" & startsLong ){
    out<-df1
    if(!is.null(traitCol)){
      if(traitCol %in% colnames(out)){
        out[[traitCol]]<-gsub("/", ".over.", out[[traitCol]])
        out[[traitCol]]<-gsub("\\'", "", out[[traitCol]]) 
      }
    }
    #* `if data is wide and mode is wide (single value traits only)`
  } else if(outputMode=="wide" & !startsLong){
    out<-df1
    #* `if data is wide and mode is long (single value traits only)`
  } else if(outputMode=="long" & !startsLong){
    #* ***** `find phenotype columns as section of numerics at end of data`
    sequence <- seq(ncol(df1), 1, -1) 
    numeric_cols <- as.numeric(which(unlist(lapply(df1, is.numeric))))
    pheno_position_in_seq<-which(unlist(lapply(1:length(numeric_cols), function(i){
      sequence[i] == rev(numeric_cols)[i]
    })))
    pheno_cols <- rev(sequence[pheno_position_in_seq])
    #* ***** `melt data`
    #* note this will warn about numeric vs integer so I am suppressing that since it should always be fine to do that.
    out <- suppressWarnings(as.data.frame(data.table::melt(data.table::as.data.table(df1), measure.vars = pheno_cols,
                          variable.name = traitCol, value.name = valueCol)))
  }
  
  colnames(out)<-gsub("/", ".over.", colnames(out))
  colnames(out)<-gsub("\\'", "", colnames(out))
  return(out)
}
