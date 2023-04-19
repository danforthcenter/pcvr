#' Earth Mover's Distance between spectral histograms
#' 
#' @description pcv.emd can be used to calculate Earth Mover's Distance between pairwise histograms in a wide dataframe of multi value traits.
#' 
#' @param df Data frame to use with multi value traits in wide format
#' @param cols Columns to use. Defaults to NULL in which case all columns are used. Single strings will be used to regex a pattern in column names (see examples). A vector of names, positions, or booleans will also work.
#' @param reorder Should data be reordered to put similar rows together in the resulting plot? This takes a vector of column names of length 1 or more (see examples).
#' @param include if a long dataframe is returned then these columns will be added to the dataframe, labelled for i and j (the row positions for compared histograms).
#' @param mat Logical, should data be returned as an nrow x nrow matrix or as a long dataframe? By Default this is FALSE and a long dataframe is returned. Both options are comparable in terms of speed, although for large datasets the matrix version may be slightly faster.
#' @param plot Logical, should a plot be returned? For a matrix this is made with image(), for a dataframe this uses ggplot.
#' @param parallel Number of cores to use. If this is above 1 then \code{parallel::mclapply} is used with this number of cores.
#' @import ggplot2
#' 
#' @keywords emd, earth mover's distance, multi-value trait, histogram
#' @examples 
#' makeHist<-function(mu, sd){hist(rnorm(10000,mu,sd), breaks=seq(1,100,1), plot=F)$counts}
#' test<-as.data.frame(do.call(rbind, lapply(seq(30,54,3), function(d) {
#'     x<-as.data.frame(do.call(rbind, lapply(1:10, function(i) makeHist(mu=d, sd=5))))
#'     x$Mu = round(d,-1)
#'     x})))
#' test<-test[sample(rownames(test), nrow(test), replace=F),] #* reorder randomly for similarity to real data
#' test$meta1<-rep(LETTERS[1:3], length.out = nrow(df))
#' test$meta2<-rep(LETTERS[4:5], length.out = nrow(df))
#' pcv.emd(test, cols="V", reorder="Mu", mat =F, plot=F, parallel = 1)
#' x<-pcv.emd(df=test, cols="V", reorder="Mu", include = c("meta1", "meta2"), mat =F, plot=F, parallel = 1)
#' head(x)
#' 
#' file = "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest1.csv"
#' df1<-read.pcv(file, "wide", T, multiValPattern = c("index_frequencies_index_ari", "index_frequencies_index_ci_rededge", "npq_hist_NPQ", "yii_hist_Fq'/Fm'", "yii_hist_Fv/Fm"))
#' colnames(df1)<-sub("index_frequencies_index_ndvi.", "ndvi_", colnames(df1))
#' pcv.emd(df1, cols="ndvi_", reorder=c("treatment", "genotype"), mat =F, plot=T, parallel = 1)
#' 
#' @export
#' 
pcv.emd<-function(df, cols=NULL, reorder=NULL, include=NULL, mat=F, plot = T, parallel = 1){
  if(!is.null(reorder)){ df<-df[order(interaction(df[,reorder])),] }
  if(is.null(cols)){cols<-colnames(df)
    }else if(is.character(cols) && length(cols)==1){cols<-grepl(cols, colnames(df))}
  if(parallel > 1){innerLapply <- function(...){parallel::mclapply(..., mc.cores=parallel)}}else{innerLapply<-lapply}
  if(mat){
    out_data<-matrix(0, nrow=nrow(df), ncol = nrow(df))
    for(i in 1:nrow(df)){
      for(j in 1:nrow(df)){
        if(i==j){emd<-0}else{emd<-emd1d(as.numeric(df[i, cols]), as.numeric(df[j, cols]))}
        out_data[i,j]<-emd
      }
    }
    if(plot){
      p<-image(out_data)
    }
  }else{
    out_data<-do.call(rbind, lapply(1:nrow(df), function(i){
      do.call(rbind, innerLapply(1:nrow(df), function(j){
        if(i==j){emdOut=0}else{emdOut=emd1d(as.numeric(df[i, cols]), as.numeric(df[j, cols]))}
        if(!is.null(include)){x<-data.frame(i=i, j=j, emd = emdOut, df[i, include], df[j,include])
          colnames(x)<-c("i", "j", "emd", paste0(include,"_i"), paste0(include, "_j"))
          } else {x<-data.frame(i=i, j=j, emd = emdOut)}
        x
      }))
    }))
    if(plot){
      p<-ggplot2::ggplot(out_data, ggplot2::aes(x=i, y=j, fill=emd))+
        ggplot2::geom_tile(color=NA)+
        ggplot2::labs(fill = "Earth Mover's Distance")+
        ggplot2::theme_minimal()+
        ggplot2::theme(axis.line.x.bottom = ggplot2::element_line(), axis.line.y.left = ggplot2::element_line(),
                       legend.position="bottom")
    }
  }
  if(plot){outList<-list("data" = out_data, "plot"=p)} else{outList<-out_data}
  return(outList)
}
