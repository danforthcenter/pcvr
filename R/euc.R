#' Euclidean Distance between spectral histograms
#' 
#' @description pcv.euc can be used to calculate Euclidean Distance between pairwise histograms 
#' in a wide dataframe of multi value traits. The is expected to be used with output from \code{mv_ag}.
#' See also \link{pcv.emd}.
#' 
#' @param df Data frame to use with multi value traits in wide format or long format
#' @param cols Columns to use. Defaults to NULL in which case all columns are used.
#' Single strings will be used to regex a pattern in column names (see examples).
#'  A vector of names, positions, or booleans will also work.
#'  For long data this is taken as a regex pattern (or full name)
#'  to use in filtering the trait column.
#' @param reorder Should data be reordered to put similar rows together in the resulting plot?
#'  This takes a vector of column names of length 1 or more (see examples).
#' @param include if a long dataframe is returned then these columns will be added to the dataframe,
#'  labelled for i and j (the row positions for compared histograms).
#'  If a matrix is returned then this information is stored in the row names.
#'  This defaults to \link{rescale}.
#' @param mat Logical, should data be returned as an nrow x nrow matrix or as a long dataframe?
#'  By Default this is FALSE and a long dataframe is returned.
#'  Both options are comparable in terms of speed,
#'  although for large datasets the matrix version may be slightly faster.
#' @param plot Logical, should a plot be returned? For a matrix this is made with heatmap(),
#'  for a dataframe this uses ggplot.
#' @param parallel Number of cores to use. Defaults to 1 unless the "mc.cores" option is set.
#' @param trait Column name for long data to identify traits. This defaults to "trait". If this and value are 
#' in the column names of the data then it is assumed to be in long format, otherwise it is assumed to be in 
#' wide format.
#' @param id A vector of column names that uniquely identifies observations if the
#' data is in long format. Defaults to "image".
#' @param value A column name for the values to be drawn from in long data.
#' Defaults to "value".
#' @param raiseError Logical, should warnings/errors be raised for potentially large output?
#'  It is easy to ask for very many comparisons with this function so the goal of this argument 
#'  is to catch a few of those and give estimates of how much time something may take.
#'   If the function is expected to take very long then a warning or an error is raised.
#'    If this is set to FALSE then no time estimates are made.
#' @import ggplot2
#' @import parallel
#' @return A dataframe/matrix (if plot=FALSE) or a list with a dataframe/matrix and a ggplot (if plot=TRUE).
#'  The returned data contains pairwise EMD values.
#' 
#' @keywords euclidean distance, multi-value trait, histogram
#' @examples 
#' 
#' ## Not run:
#' 
#' makeHist<-function(mu, sd){hist(rnorm(10000,mu,sd), breaks=seq(1,100,1), plot=FALSE)$counts}
#' test<-as.data.frame(do.call(rbind, lapply(seq(30,54,6), function(d) {
#'     x<-as.data.frame(do.call(rbind, lapply(1:5, function(i) makeHist(mu=d, sd=5))))
#'     x$Mu = round(d,-1)
#'     x})))
#' test<-test[sample(rownames(test), nrow(test), replace=FALSE),]
#' test$meta1<-rep(LETTERS[1:3], length.out = nrow(test))
#' test$meta2<-rep(LETTERS[4:5], length.out = nrow(test))
#' 
#' x<-pcv.euc(df=test, cols="V", reorder="Mu",
#'    include = c("meta1", "meta2"), mat =FALSE,
#'    plot=FALSE, parallel = 1)
#' head(x)
#' 
#' if(FALSE){
#' file = paste0("https://media.githubusercontent.com/media/joshqsumner/",
#'               "pcvrTestData/main/pcv4-multi-value-traits.csv")
#' df1<-read.pcv(file, "wide")
#' 
#' df1$genotype = substr(df1$barcode, 3,5)
#' df1$genotype = ifelse(df1$genotype == "002", "B73",
#'                      ifelse(df1$genotype == "003", "W605S",
#'                             ifelse(df1$genotype == "004", "MM", "Mo17")))
#' df1$fertilizer = substr(df1$barcode, 8, 8)
#' df1$fertilizer = ifelse(df1$fertilizer == "A", "100",
#'                        ifelse(df1$fertilizer == "B", "50", "0"))
#' 
#' w<-pcv.euc(df1, cols="hue_frequencies", reorder=c("fertilizer", "genotype"), 
#'   mat =FALSE, plot=TRUE, parallel = 1)
#' }
#' ## End(Not run)
#' 
#' @export
#' 
pcv.euc <- function(df, cols = NULL, reorder = NULL, include = reorder, mat = FALSE, plot = TRUE,
                    parallel = getOption("mc.cores", 1), trait = "trait", id = "image",
                    value = "value", raiseError = TRUE, method = "euc") {
  pcv.emd(df, cols, reorder, include, mat, plot, parallel, trait, id, value, raiseError, method)
}

