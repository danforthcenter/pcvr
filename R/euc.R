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
pcv.euc<-function(df, cols=NULL, reorder=NULL, include=reorder, mat=FALSE, plot =TRUE,
                  parallel = getOption("mc.cores",1), trait="trait", id="image", 
                  value="value", raiseError=TRUE){
  # df = df1; cols="ndvi_"; reorder=c("treatment", "genotype"); mat =FALSE; plot=TRUE; parallel = 1; include=reorder; trait=FALSE; id="image";value="value"
  # df_long<-read.pcv(file, "long", FALSE)
  # df = df_long; cols="index_frequencies_index_ndvi"; reorder=c("treatment", "genotype"); mat =FALSE; plot=TRUE; parallel = 1; include=reorder;
  # longTrait="trait"; id="image"; value="value"
  
  if(all(c(trait, value) %in% colnames(df))){
    long = TRUE
    traitCol = trait
  } else{
    long=FALSE
  }
  
  if(!is.null(reorder)){ df<-df[order(interaction(df[,reorder])),] }
  if(long){
    df<-df[grepl(cols, df[[traitCol]]), ]
    #* if nrow df is too high then do calculation for time per core based on estimates and report time with error.
    if(raiseError){
      eT_sec = 0.008*((nrow(df)/parallel)^2)
      eT_min = eT_sec/60
      eT_hour = eT_min/60
      if(eT_sec <= 300){message(paste0("Estimated time of calculation is roughly ", round(eT_sec,1), " seconds using ", parallel, " cores in parallel."))
      } else if(eT_min < 60){warning(paste0("Estimated time of calculation is roughly ", round(eT_min,2), " minutes using ", parallel, " cores in parallel."))
      } else if(eT_min > 60){stop(paste0("Stopping, estimated time of calculation is roughly ", round(eT_hour,2), " hours using ", parallel, " cores in parallel.",
                                         "\nIf you wish to proceed then rerun this command with raiseError=FALSE"))}
    }
    df$INNER_ID_EUC<-interaction(df[,id], drop=TRUE)
    if(mat){ # make dist matrix
      mat_obj <- matrix(0, nrow=length(unique(df$INNER_ID_EUC)), ncol = length(unique(df$INNER_ID_EUC)))
      values<- unlist(lapply(unique(df$INNER_ID_EUC), function(i){parallel::mclapply(unique(df$INNER_ID_EUC), function(j){
          if(i==j){0}else{
            euc1d(as.numeric(df[df$INNER_ID_EUC==as.character(i), value]),
                  as.numeric(df[df$INNER_ID_EUC==as.character(j), value])) }
        }, mc.cores=parallel)}))
      mat_obj[lower.tri(mat_obj)] <- values
      tmat_obj <- t(mat_obj)
      mat_obj[upper.tri(mat_obj)] <- tmat_obj[upper.tri(tmat_obj)]
      rownames(mat_obj) <- colnames(mat_obj) <- unique(df$INNER_ID_EUC)
      out_data <- mat_obj
      
    }else{ # make long data
      out_data<-do.call(rbind, lapply(1:length(unique(df$INNER_ID_EUC)), function(i_n){
        do.call(rbind, parallel::mclapply(1:length(unique(df$INNER_ID_EUC)), function(j_n){
          i <- unique(df$INNER_ID_EUC)[i_n]
          j <- unique(df$INNER_ID_EUC)[j_n]
          eucOut<- NULL
          if(i_n==j_n){eucOut=0}else if(i_n<j_n){eucOut=euc1d(as.numeric(df[df$INNER_ID_EUC==as.character(i), value]),
                                              as.numeric(df[df$INNER_ID_EUC==as.character(j), value]) )}
          if(!is.null(eucOut)){
            if(!is.null(include)){
              x<-rbind(data.frame(i=i, j=j, euc = eucOut,
                                                df[df$INNER_ID_EUC==as.character(i), include][1,],
                                                df[df$INNER_ID_EUC==as.character(j), include][1,]),
                       data.frame(i=j, j=i, euc = eucOut,
                                  df[df$INNER_ID_EUC==as.character(j), include][1,],
                                  df[df$INNER_ID_EUC==as.character(i), include][1,]) )
              
            colnames(x)<-c("i", "j", "euc", paste0(include,"_i"), paste0(include, "_j"))
            } else {x<-data.frame(i=c(i,j), j=c(j,i), euc = eucOut)}
            x
          }
        }, mc.cores=parallel))
      }))
    }
  } else{
    if(is.null(cols)){cols<-colnames(df)
    }else if(is.character(cols) && length(cols)==1){cols<-grepl(cols, colnames(df))}
    if(raiseError){
      eT_sec = 0.0025*((nrow(df)/parallel)^2)
      eT_min = eT_sec/60
      eT_hour = eT_min/60
      if(eT_sec <= 300){message(paste0("Estimated time of calculation is roughly ", round(eT_sec,1), " seconds using ", parallel, " cores in parallel."))
      } else if(eT_min < 60){warning(paste0("Estimated time of calculation is roughly ", round(eT_min,2), " minutes using ", parallel, " cores in parallel."))
      } else if(eT_min > 60){stop(paste0("Stopping, estimated time of calculation is roughly ", round(eT_hour,2), " hours using ", parallel, " cores in parallel.",
                                         "\nIf you wish to proceed then rerun this command with raiseError=FALSE"))}
    }
    if(mat){# make dist matrix
      mat_obj <- matrix(0, nrow=nrow(df), ncol = nrow(df) )
      values<-unlist(lapply(1:nrow(df), function(i){parallel::mclapply(1:nrow(df), function(j){
          if(i==j){0}else{euc1d(as.numeric(df[i, cols]), as.numeric(df[j, cols]))}
        }, mc.cores = parallel)}))
      mat_obj[lower.tri(mat_obj)] <- values
      tmat_obj <- t(mat_obj)
      mat_obj[upper.tri(mat_obj)] <- tmat_obj[upper.tri(tmat_obj)]
      rownames(mat_obj) <- colnames(mat_obj) <- 1:nrow(df)
      out_data <- mat_obj
      if(!is.null(include)){
        rownames(out_data)<-interaction(df[,include])
      }
    }else{# make long dataframe
      out_data<-do.call(rbind, lapply(1:nrow(df), function(i){
        do.call(rbind, parallel::mclapply(1:nrow(df), function(j){
          eucOut <- NULL
          if(i==j){eucOut=0}else if(i<j){eucOut=euc1d(as.numeric(df[i, cols]), as.numeric(df[j, cols]))}
          if(!is.null(eucOut)){
            if(!is.null(include)){
              x<-rbind(data.frame(i=i, j=j, euc = eucOut, df[i, include], df[j,include]),
                       data.frame(i=j, j=i, euc = eucOut, df[j, include], df[i,include]) )
              
              colnames(x)<-c("i", "j", "euc", paste0(include,"_i"), paste0(include, "_j"))
            } else {x<-data.frame(i=c(i,j), j=c(j,i), euc = eucOut)}
            x
          }
        }, mc.cores=parallel))
      }))
    }
  }
  if(plot){
    if(mat){
      p<-stats::heatmap(out_data)
    }else{
      p<-ggplot2::ggplot(out_data, ggplot2::aes(x=.data$i, y=.data$j, fill=.data$euc))+
        ggplot2::geom_tile(color=NA)+
        ggplot2::labs(fill = "Euclidean Distance")+
        ggplot2::theme_minimal()+
        ggplot2::theme(axis.line.x.bottom = ggplot2::element_line(), axis.line.y.left = ggplot2::element_line(),
                       legend.position="bottom")
    }
  }
  if(plot){outList<-list("data" = out_data, "plot"=p)} else{outList<-out_data}
  return(outList)
}



#' Euclidean Distance between spectral histograms
#' 
#' @description euc1d computes euclidean distance between two samples.
#' 
#' @param s1 Histogram as a numeric vector of counts per position.
#' @param s2 Histogram as a numeric vector of counts per position. Must be the same length as s1.
#' 
#' @importFrom stats dist
#' @keywords internal
#' @return Returns euclidean distance as a numeric
#' 
#' @noRd

euc1d<-function(s1, s2){
  if(length(s1)!=length(s2)){stop("Samples must be from the same histogram and be of the same length")}
  s1<-s1/sum(s1)
  s2<-s2/sum(s2)
  mat <- matrix(c(s1,s2), nrow=2, byrow=TRUE)
  euc <- as.numeric(stats::dist(mat, method="euclidean"))
  return(euc)
}

