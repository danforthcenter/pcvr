#' Read in plantCV csv from bellwether phenotyper style experiments analyzed with plantCV versions <4.
#'
#' @param file Path to the version 3 plantCV output containing phenotypes.
#' @param snapshotFile path to the snapshot info metadata file, typically called SnapshotInfo.csv.
#' This needs to have a column name corresponding to `joinSnapshot` (defaults to "id")
#' which can be used to join the snapshot data to the phenotype data.
#' Generally this joining will happen through a parsed section of the file path
#' to each image present in the phenotype data.
#' This means that including a duplicate name in `metaForm` will be overwritten
#' by parsing image paths, so `metaForm` and `joinSnapshot` should not have duplicated names.
#'  If there is a timestamp column in the snapshot data then it will
#'  be converted to datetime (assuming a "Y-m-d H:M:S" format)
#'  and used to calculate days after starting (DAS) and hours.
#' @param designFile path to a csv file which contains experimental design information
#' (treatments, genotypes, etc) and which will be joined to phenotype
#' and snapshot data through all shared columns.
#' @param metaCol a column name from the phenotype data read in with the `file` argument.
#'  Generally for bellwether experiments this will correspond to an image path.
#'  The name is split on "/" characters with the last segment being taken and parsed into
#'  some number of sections based on `metaForm`.
#' @param metaForm A character string or character vector of column names to parse `metaCol` into.
#'  The number of names needs to match with length of `metaCol` when parsed.
#'   If a character string is provided then it is assumed to be underscore delimited,
#'   so do if you need underscores in a column name then use `c("column_one", "column_two",...)`
#'   instead of `column_one_column_two_...`.
#' @param joinSnapshot Column name create in phenotype data to use in joining snapshot data.
#' By default this will attempt to make an "id" column, which is parsed from a snapshot
#' folder in `metaCol` ("/shares/sinc/data/Phenotyper/SINC1/ImagesNew/**snapshot1403**/").
#'  An error will be raised if this column is not present in the snapshot data.
#' @param conversions A named list of phenotypes that should be rescaled by the value in the list.
#' For instance, at zoom 1  `list(area = 13.2 * 3.7/46856)`
#' will convert from pixels to square cm in the 5MP bellwether camera.
#' @param mode The mode to read data in with through read.pcv.
#' The default is "long" because this function is built for pcv3 output,
#' which was generally a wider format to start with than pcv4 output.
#' @param ... Other arguments passed to \code{read.pcv}.
#' @keywords read.csv pcv3
#' @importFrom utils read.csv
#' @return Returns a dataframe potentially with several files merged into it.
#'
#' @examples
#' \donttest{
#' base_url <- "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/"
#' bw <- read.pcv.3(
#'   file = paste0(base_url, "pcv3Phenos.csv"),
#'   metaCol = NULL,
#'   reader = "fread"
#' )
#' bw <- read.pcv.3(
#'   file = paste0(base_url, "pcv3Phenos.csv"),
#'   metaCol = "meta", metaForm = "vis_view_angle_zoom_horizontal_gain_exposure_v_new_n_rep",
#'   joinSnapshot = "id",
#'   reader = "fread"
#' )
#' bw <- read.pcv.3(
#'   file = paste0(base_url, "pcv3Phenos.csv"),
#'   snapshotFile = paste0(base_url, "pcv3Snapshot.csv"),
#'   designFile = paste0(base_url, "pcv3Design.csv"),
#'   metaCol = "meta", metaForm = "vis_view_angle_zoom_horizontal_gain_exposure_v_new_n_rep",
#'   joinSnapshot = "id", conversions = list(area = 13.2 * 3.7 / 46856),
#'   reader = "fread"
#' )
#' }
#'
#' @export

read.pcv.3 <- function(
    file = NULL,
    snapshotFile = NULL,
    designFile = NULL,
    metaCol = "meta",
    metaForm = "vis_view_angle_zoom_horizontal_gain_exposure_v_new_n_rep",
    joinSnapshot = "id",
    conversions = NULL,
    mode = "long",
    ...) {
  phenos <- read.pcv(filepath = file, mode = mode, ...)

  #* `parse metadata`
  if (!is.null(metaCol)) {
    metaToParse <- unlist(lapply(phenos[[metaCol]], function(meta) {
      x <- strsplit(meta, "/")[[1]]
      sub("[.]png|[.]jpg|[.]jpeg|[.]INF", "", x[length(x)])
    }))
    phenoMeta <- do.call(rbind, lapply(metaToParse, function(meta) {
      strsplit(meta, "[_]|[.]|[-]")[[1]]
    }))
    if (!is.null(metaForm)) {
      if (length(metaForm == 1)) {
        metaColNames <- strsplit(metaForm, "_")[[1]]
      } else {
        metaColNames <- metaForm
      }

      colnames(phenoMeta) <- metaColNames
    }
    phenos <- cbind(phenos, phenoMeta)
    if (!is.null(joinSnapshot)) {
      phenos[[joinSnapshot]] <- sapply(
        phenos[[metaCol]],
        function(i) {
          strsplit(strsplit(i, "/snapshot")[[1]][2], "/")[[1]][1]
        }
      )
    }
  }
  if (is.list(conversions)) {
    for (pheno in names(conversions)) {
      phenos[[paste0(pheno, "_adj")]] <- phenos[[pheno]] * conversions[[pheno]]
    }
  }
  #* `Add snapshot data`
  if (!is.null(snapshotFile)) {
    snp <- read.csv(snapshotFile)
    if (is.null(joinSnapshot)) {
      if (!joinSnapshot %in% colnames(snp)) {
        stop(paste0("joinSnapshot (", joinSnapshot, ") not in snapshot data column names"))
      }
    }
    if (any(colnames(snp) == "tiles")) {
      snp <- snp[snp$tiles != "", ]
    }
    if (any(grepl("barcode", colnames(snp)))) {
      colnames(snp)[which(grepl("barcode", colnames(snp), ignore.case = TRUE))] <- "Barcodes"
    }
    phenos <- merge(phenos, snp, by = joinSnapshot)
    #* `parse time data`
    if (any(grepl("time", colnames(snp)))) {
      timeCol <- match.arg("time", colnames(snp))
      tryCatch(expr = {
        phenos[[timeCol]] <- as.POSIXct(strptime(phenos[[timeCol]], format = "%Y-%m-%d %H:%M:%S"))
        beg <- min(phenos[[timeCol]])
        phenos$DAS <- floor(as.numeric((phenos[[timeCol]] - beg) / 60 / 60 / 24))
        phenos$hour <- as.numeric(format(phenos[[timeCol]], "%H"))
      }, error = function(err) {
        warning("Error raised while parsing time, skipping time parsing.")
      })
    }
  }
  #* `Add design data`
  if (!is.null(designFile)) {
    des <- read.csv(designFile)
    bycol <- colnames(phenos)[colnames(phenos) %in% colnames(des)]
    phenos <- merge(phenos, des, by = bycol, all.x = TRUE)
  }
  return(phenos)
}
