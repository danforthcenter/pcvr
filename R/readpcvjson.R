#' Read in plantCV json keeping either single value, multi value, or all traits.
#'
#'
#' @param file Path to the plantCV output json.
#' @param dtSuffix Suffix for columns denoting data type.
#'  This will be one of the innermost keys in the json. Default is "datatype"
#' @param valSuffix Suffix for columns denoting phenotype values.
#' This will be one of the innermost keys in the json. Default is "value"
#' @param labSuffix Suffix for columns denoting phenotype labels.
#' This will be one of the innermost keys in the json. Default is "label"
#' @param output What phenotypes should be returned?
#' This can be a vector of names present in "label" or one of the built in options.
#' Built in options are "sv" for single value traits only,
#' "mv" for multi value traits only, and "all" for all traits.
#' Note that by default tuple value data is output with single value data in a wide format.
#' @import jsonlite
#' @keywords read.csv, pcv, bellwether
#' @return Returns a dataframe.
#' @export
#' @examples
#'
#' ## Not run:
#'
#' if (FALSE) {
#'   file <- "path/to/your/json/file.json"
#'   x <- read.pcv.jsn(file, output = "sv")
#'   y <- read.pcv.jsn(file, output = c("area", "perimeter"))
#'   z <- read.pcv.jsn(file, output = "mv")
#' }
#'
#' ## End(Not run)
read.pcv.jsn <- function(file, dtSuffix = "datatype", valSuffix = "value",
                         labSuffix = "label", output = "sv") {
  jsn <- as.data.frame(jsonlite::fromJSON(file, flatten = TRUE))
  fieldTypes <- unique(unlist(lapply(colnames(jsn), function(c) {
    x <- strsplit(c, "[.]")[[1]]
    x[length(x)]
  })))

  if (!all(c(dtSuffix, valSuffix, labSuffix) %in% fieldTypes)) {
    names <- c(dtSuffix, valSuffix, labSuffix)
    stop(paste0(
      "'", names[which(!names %in% fieldTypes)],
      "' not in column suffixes. Suffixes in data are: ",
      paste(fieldTypes, collapse = ", ")
    ))
  }
  dt <- paste0(".", dtSuffix, "$")
  val <- paste0(".", valSuffix, "$")
  lab <- paste0(".", labSuffix, "$")

  if (!all(output %in% c("sv", "mv", "all"))) {
    jsn <- jsn[, grepl(paste0(c(paste0("[.]", output, "[.]"), "metadata"),
                              collapse = "|"), colnames(jsn))]
    output <- "all"
  }

  output <- match.arg(output, choices = c("sv", "mv", "all"))

  #* `go over all datatype columns and make lists of each kind of column`

  dtcols <- colnames(jsn)[grepl(dt, colnames(jsn))]
  datatypes <- unlist(lapply(dtcols, function(c) {
    strsplit(unique(jsn[!is.na(jsn[, c]), c]), "[']")[[1]][2]
  }))
  dtl <- split(gsub(dt, "", dtcols), datatypes)

  #* `Keeping only the json key values that I need`

  sub <- jsn[, grepl(paste0(c(lab, val), collapse = "|"), colnames(jsn))]

  #* `separate single/multi value traits and metadata`

  listCols <- if (is.null(dtl[["list"]])) {
    "NOLISTCOLUMNS"
  } else {
    dtl[["list"]]
  }
  tupleCols <- if (is.null(dtl[["tuple"]])) {
    "NOTUPLECOLUMNS"
  } else {
    dtl[["tuple"]]
  }
  sv <- sub[, !grepl(paste0(listCols, collapse = "|"), colnames(sub))]
  sv <- sv[, !grepl("metadata", colnames(sv))]
  if (length(tupleCols) > 0) {
    tv <- sv[, grepl(paste0(tupleCols, collapse = "|"), colnames(sv))]
  }
  sv <- sv[, !grepl(paste0(tupleCols, collapse = "|"), colnames(sv))]
  mv <- sub[, grepl(paste0(listCols, collapse = "|"), colnames(sub))]
  mv <- mv[, !grepl("metadata", colnames(mv))]
  meta <- sub[, grepl("metadata", colnames(sub))]
  meta <- meta[, !grepl(lab, colnames(meta))]
  colnames(meta) <- gsub(paste0(".*metadata.|", val), "", colnames(meta))

  #* `format single value column names`
  if (output %in% c("sv", "all")) {
    colnames(sv) <- unlist(lapply(colnames(sv), function(c) {
      x <- strsplit(c, "[.]")[[1]]
      l <- length(x)
      return(paste0(x[l - 1], ".", x[l]))
    }))

    svl <- sv[, grepl(lab, colnames(sv))]
    labs <- unlist(lapply(svl, function(c) unique(c[!is.na(c)])))
    sv <- sv[, grepl(val, colnames(sv))]
    colnames(sv) <- paste(gsub(val, "", colnames(sv)), labs, sep = ".")
    colnames(sv) <- gsub(".none$", "", colnames(sv))

    #* `format tuple data`
    if (length(tupleCols) > 0 && tupleCols != "NOTUPLECOLUMNS") {
      colnames(tv) <- unlist(lapply(colnames(tv), function(c) {
        x <- strsplit(c, "[.]")[[1]]
        l <- length(x)
        return(paste0(x[l - 1], ".", x[l]))
      }))
      tupL <- lapply(tv, function(tup) {
        tup <- tup[!is.null(tup) & !is.na(tup)]
        length(unlist(tup[1]))
      })

      tv <- do.call(cbind, lapply(colnames(tv)[grepl(val, colnames(tv))], function(c) {
        l <- tupL[[c]]
        labs <- unique(unlist(tv[[gsub(val, ".label", c)]]))
        tupleExpanded <- data.frame(do.call(rbind, lapply(tv[[c]], function(c) {
          if (any(is.null(c)) || any(is.na(c))) {
            list(rep(NULL, l))
          } else {
            c
          }
        })))
        colnames(tupleExpanded) <- paste0(gsub(val, "", c), ".", labs)
        tupleExpanded
      }))
      if (nrow(tv) == nrow(sv)) {
        sv <- cbind(sv, tv)
      }
    }
  }
  if (output %in% c("mv", "all")) {
    #* `format list data`
    colnames(mv) <- unlist(lapply(colnames(mv), function(c) {
      x <- strsplit(c, "[.]")[[1]]
      l <- length(x)
      return(paste0(x[l - 1], ".", x[l]))
    }))

    mvL <- do.call(rbind, lapply(colnames(mv)[grepl(val, colnames(mv))], function(c) {
      lng <- do.call(rbind, lapply(seq_len(nrow(mv)), function(i) {
        data.frame(
          row = i, value = unlist(mv[i, c]),
          label = unlist(mv[i, gsub(val, paste0(".", labSuffix), c)]),
          trait = gsub(val, "", c)
        )
      }))
      lng
    }))
  }

  #* `Add metadata back to values`
  ret <- .jsonAddMetadata(output, mvL, meta, sv)
  return(ret)
}

#' @description
#' Internal function for adding metadata back to json
#' @param priors priors as a list
#' @keywords internal
#' @noRd

.jsonAddMetadata <- function(output, mvL, meta, sv) {
  if (output %in% c("mv", "all")) {
    mv <- cbind(mvL, meta[mvL$row, ])
  }
  if (output %in% c("sv", "all")) {
    sv <- cbind(meta, sv)
  }
  if (output == "sv") {
    ret <- sv
  } else if (output == "mv") {
    ret <- mv
  } else if (output == "all") {
    ret <- list()
    if (all(dim(sv) > 0)) {
      ret[["SV"]] <- sv
    }
    if (all(dim(mv) > 0)) {
      ret[["MV"]] <- mv
    }
    if (length(ret) == 1) {
      ret <- ret[[1]]
    }
  }
  return(ret)
}
