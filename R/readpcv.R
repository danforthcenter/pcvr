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
#'   This is used with traitCol when \code{mode="wide"} to identify
#'   unique traits since some may be ambiguous
#' (ellipseCenter.x vs ellipseCenter.y, bins of histograms, etc)
#' @param valueCol Column with phenotype values, defaults to "value".
#'   This should generally not need to be changed from the default.
#' @param reader The function to use to read in data,
#'   defaults to NULL in which case \code{data.table::fread} is used if filters are in place
#'   and \code{read.csv} is used otherwise.
#'   Note that if you use \code{read.csv} with filters in place then you will need to specify
#'   \code{header=FALSE} so that the piped output from awk is read correctly.
#'   If fread is too slow for your needs then \code{vroom::vroom()} may be useful.
#' @param filters If a very large pcv output file is read then it may be desireable
#'   to subset it before reading it into R, either for ease of use or because of RAM limitations.
#'   The filter argument works with "COLUMN in VALUES" syntax. This can either be a character vector
#'   or a list of character vectors. In these vectors there needs to be a column name,
#'   one of " in ", " is ", or " = " to match the string exactly, or "contains"
#'   to match with awk style regex, then a set of comma delimited values to filter
#'   that column for (see examples). Note that this and awk both use awk through \code{pipe()}.
#'   This functionality will not work on a windows system.
#' @param awk As an alternative to filters a direct call to awk can be supplied here,
#'   in which case that call will be used through \code{pipe()}.
#' @param ... Other arguments passed to the reader function.
#'   In the case of 'fread' there are several defaults provided already
#'   which can be overwritten with these extra arguments.
#'
#' @details
#' In plantCV version 4 the single value traits are returned in wide format from \code{json2csv}
#' and the multi value traits are returned in long format. Briefly plantCV data was returned as one
#' long table which sparked the emphasis in this function on reading data quickly and parsing it
#' outside of R. With the current plantCV output these options are largely unnecessary.
#' When data is read in using read.pcv the traitCol, valueCol, and labelCol arguments are checked
#' to determine if the data is in long format. This is done to keep compatibility with interim
#' versions of plantcv output where all outputs were in a single long format file.
#'
#' With the current implementation and plantcv output you can read wide or long format files into
#' wide or long format in R. Keep in mind that the 'mode' argument controls the format that will be
#' returned in R, not the format that the data saved as in your csv file.
#'
#' @keywords read.csv pcv4
#' @return Returns a data.frame in wide or long format.
#' @importFrom stats as.formula
#' @import data.table
#' @examples
#' \donttest{
#' tryCatch({
#' mv <- paste0(
#'   "https://media.githubusercontent.com/media/joshqsumner/",
#'   "pcvrTestData/main/pcv4-multi-value-traits.csv"
#' )
#' sv <- paste0(
#'   "https://raw.githubusercontent.com/joshqsumner/",
#'   "pcvrTestData/main/pcv4-single-value-traits.csv"
#' )
#'
#' w2w <- read.pcv(sv, mode = "wide", reader = "fread")
#' dim(w2w)
#'
#' w2l <- read.pcv(sv, mode = "long", reader = "fread")
#' dim(w2l)
#'
#' l2w <- read.pcv(mv, mode = "wide", reader = "fread")
#' dim(l2w)
#'
#' l2l <- read.pcv(mv, mode = "long", reader = "fread")
#' dim(l2l)
#' }, error = function(e) {message(e)})
#' }
#'
#' @export

read.pcv <- function(filepath, mode = NULL,
                     traitCol = "trait", labelCol = "label", valueCol = "value",
                     reader = NULL, filters = NULL, awk = NULL, ...) {
  if (is.null(filters) && is.null(awk)) {
    if (is.null(reader)) {
      reader <- "read.csv"
    }
    if (reader != "fread") {
      readingFunction <- match.fun(reader)
    } else {
      readingFunction <- data.table::fread
    }
    df1 <- as.data.frame(readingFunction(filepath, ...))
  } else {
    if (is.null(reader)) {
      reader <- "fread"
    }
    df1 <- pcv.sub.read(inputFile = filepath, filters = filters, reader = reader, awk = awk, ...)
    if (nrow(df1) < 1) {
      stop(paste0(
        "0 Rows returned using awk statement:\n", awkHelper(filepath, filters),
        "\nMost common issues are misspellings or not including a column name and affector."
      ))
    }
  }
  #* `check original data format`
  checkDataStateRes <- .readpcvCheckDataState(df1, traitCol, valueCol, labelCol, mode)
  startsLong <- checkDataStateRes[["startsLong"]]
  outputMode <- checkDataStateRes[["outputMode"]]
  #* `if data is long and mode is wide`
  out <- .readpcvReshapeData(df1, outputMode, startsLong, traitCol, valueCol, labelCol)
  colnames(out) <- gsub("/", ".over.", colnames(out))
  colnames(out) <- gsub("\\'", "", colnames(out))
  return(out)
}


#' @description
#' Internal function for checking input data to
#' @param priors priors as a list
#' @keywords internal
#' @noRd

.readpcvCheckDataState <- function(df1, traitCol, valueCol, labelCol, mode) {
  if (all(c(traitCol, valueCol, labelCol) %in% colnames(df1))) {
    startsLong <- TRUE
  } else if (!any(c(traitCol, valueCol, labelCol) %in% colnames(df1))) {
    startsLong <- FALSE
  } else {
    found <- c(
      "traitCol",
      "valueCol",
      "labelCol"
    )[which(c(traitCol, valueCol, labelCol) %in% colnames(df1))]
    warning(paste0(
      paste(found, collapse = ", "),
      " found in column names of data but either all or none of traitCol,",
      " valueCol, and labelCol are expected. Data will be returned as is."
    ))
    startsLong <- FALSE
  }

  if (is.null(mode)) {
    if (startsLong) {
      outputMode <- "long"
    } else {
      outputMode <- "wide"
    }
  } else {
    outputMode <- match.arg(mode, c("wide", "long"))
  }
  return(list(
    "startsLong" = startsLong,
    "outputMode" = outputMode
  ))
}

#' @description
#' Internal function for reshaping data in read.pcv
#' @param priors priors as a list
#' @keywords internal
#' @noRd

.readpcvReshapeData <- function(df1, outputMode, startsLong, traitCol, valueCol, labelCol) {
  if (startsLong) {
    #* `if data is long and mode is wide`
    if (outputMode == "wide") {
      long <- df1
      if (substr(colnames(long)[1], 1, 1) == "X" && length(unique(long[[1]])) == nrow(long)) {
        long <- long[, -1]
      }
      long <- long[!is.na(long[[valueCol]]), ]
      long[[labelCol]] <- ifelse(is.na(long[[labelCol]]), "none", long[[labelCol]])
      wide <- as.data.frame(data.table::dcast(data.table::as.data.table(long),
        as.formula(paste0("... ~ ", traitCol, "+", labelCol)),
        value.var = valueCol, sep = "."
      ))
      colnames(wide) <- sub(".none$", "", colnames(wide))
      if (any(grepl("hist|frequencies", colnames(wide)))) { # reorder the MV traits by their bins
        #* get a list of the unique non-numeric parts
        histCols <- colnames(wide)[grepl("hist|frequencies", colnames(wide))]
        uniqueMvTraits <- unique(gsub("[.]+$", "", gsub("[0-9]+", "", histCols)))
        #* for each unique non-numeric part, sort the names
        mvColsReordered <- unlist(lapply(uniqueMvTraits, function(umt) {
          iterCols <- histCols[grepl(umt, histCols)]
          iterColsNumeric <- as.numeric(gsub(paste0(umt, "."), "", iterCols))
          bins_order <- sort(iterColsNumeric, index.return = TRUE)$ix
          iterCols[bins_order]
        }))
        #* combine the histCols and the other columns, in the new order.
        sv_and_meta_cols <- colnames(wide)[!grepl("hist|frequencies", colnames(wide))]
        wide <- wide[, c(sv_and_meta_cols, mvColsReordered)]
      }
      out <- wide

      #* `if data is long and mode is long`
    } else if (outputMode == "long") {
      out <- df1
      out[[traitCol]] <- gsub("/", ".over.", out[[traitCol]])
      out[[traitCol]] <- gsub("\\'", "", out[[traitCol]])
      #* `if data is wide and mode is wide (single value traits only)`
    }
  } else {
    if (outputMode == "wide") {
      out <- df1
      #* `if data is wide and mode is long (single value traits only)`
    } else if (outputMode == "long") {
      #* ***** `find phenotype columns as section of numerics at end of data`
      sequence <- seq(ncol(df1), 1, -1)
      numeric_cols <- as.numeric(which(unlist(lapply(df1, is.numeric))))
      pheno_position_in_seq <- which(unlist(lapply(seq_along(numeric_cols), function(i) {
        sequence[i] == rev(numeric_cols)[i]
      })))
      pheno_cols <- rev(sequence[pheno_position_in_seq])
      #* ***** `melt data`
      #* note this will warn about numeric vs integer
      #* so I am suppressing that since it should always be fine to do that.
      out <- suppressWarnings(
        as.data.frame(
          data.table::melt(
            data.table::as.data.table(df1),
            measure.vars = pheno_cols, variable.name = traitCol,
            value.name = valueCol
          )
        )
      )
    }
  }
  return(out)
}
