#' Earth Mover's Distance between spectral histograms
#' 
#' @description pcv.emd can be used to calculate Earth Mover's Distance between pairwise histograms in a wide dataframe of multi value traits.
#' 
#' @param df Data frame to use with multi value traits in wide format or long format
#' @param cols Columns to use. Defaults to NULL in which case all columns are used. Single strings will be used to regex a pattern in column names (see examples). A vector of names, positions, or booleans will also work. For long data this is taken as a regex pattern (or full name) to use in filtering the longTrait column.
#' @param reorder Should data be reordered to put similar rows together in the resulting plot? This takes a vector of column names of length 1 or more (see examples).
#' @param include if a long dataframe is returned then these columns will be added to the dataframe, labelled for i and j (the row positions for compared histograms). If a matrix is returned then this information is stored in the row names. This defaults to \link{rescale}.
#' @param mat Logical, should data be returned as an nrow x nrow matrix or as a long dataframe? By Default this is FALSE and a long dataframe is returned. Both options are comparable in terms of speed, although for large datasets the matrix version may be slightly faster.
#' @param plot Logical, should a plot be returned? For a matrix this is made with image(), for a dataframe this uses ggplot.
#' @param parallel Number of cores to use. Defaults to 1 unless the "mc.cores" option is set.
#' @param longTrait Defaults to NULL, in which case the data is assumed to be in long format. If this is a character string then it is taken as a column name of long data and the other arguments will assume data is long.
#' @param id A vector of column names that uniquely identifies observations if the data is in long format. Defaults to "image".
#' @param value A column name for the values to be drawn from in long data. Defaults to "value".
#' @param raiseError Logical, should warnings/errors be raised for potentially large output? It is easy to ask for very many comparisons with this function so the goal of this argument is to catch a few of those and give estimates of how much time something may take. If the function is expected to take very long then a warning or an error is raised. If this is set to FALSE then no time estimates are made.
#' @import ggplot2
#' @return A dataframe/matrix (if plot=F) or a list with a dataframe/matrix and a ggplot (if plot=T). The returned data contains pairwise EMD values.
#' 
#' @keywords emd, earth mover's distance, multi-value trait, histogram
#' @examples 
#' makeHist<-function(mu, sd){hist(rnorm(10000,mu,sd), breaks=seq(1,100,1), plot=F)$counts}
#' test<-as.data.frame(do.call(rbind, lapply(seq(30,54,3), function(d) {
#'     x<-as.data.frame(do.call(rbind, lapply(1:10, function(i) makeHist(mu=d, sd=5))))
#'     x$Mu = round(d,-1)
#'     x})))
#' test<-test[sample(rownames(test), nrow(test), replace=F),] #* reorder randomly for similarity to real data
#' test$meta1<-rep(LETTERS[1:3], length.out = nrow(test))
#' test$meta2<-rep(LETTERS[4:5], length.out = nrow(test))
#' pcv.emd(test, cols="V", reorder="Mu", mat =F, plot=F, parallel = 1)
#' x<-pcv.emd(df=test, cols="V", reorder="Mu", include = c("meta1", "meta2"), mat =F, plot=F, parallel = 1)
#' head(x)
#' 
#' file = "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest1.csv"
#' df1<-read.pcv(file, "wide", T, multiValPattern = c("index_frequencies_index_ari", "index_frequencies_index_ci_rededge", "npq_hist_NPQ", "yii_hist_Fq'/Fm'", "yii_hist_Fv/Fm"))
#' colnames(df1)<-sub("index_frequencies_index_ndvi.", "ndvi_", colnames(df1))
#' w<-pcv.emd(df1, cols="ndvi_", reorder=c("treatment", "genotype"), mat =F, plot=T, parallel = 1)
#' df_long<-read.pcv(file, "long", F)
#' l<-pcv.emd(df = df_long, cols="index_frequencies_index_ndvi", reorder=c("treatment", "genotype"), mat =F, plot=T, longTrait="trait", id="image", value="value")
#' l$plot + theme(axis.text = element_blank())
#' 
#' #* Note on computational complexity
#' #* This scales as O^2, see the plot below for some idea of the time for different input data sizes.
#' emdTime<-function(x, n=1){
#' x^2 / n * 0.0023
#' }
#' plot(x=c(18,36,54,72, 108, 135), y = c(0.74, 2.89, 6.86, 10.99, 26.25, 42.44), xlab="N Input Images", ylab="time (seconds)") # benchmarked test data
#' lines(x=1:150, y = emdTime(1:150)) # exponential function
#' 
#' plot(x=1:1000, y=emdTime(1:1000), type="l", xlab="N Input Images", ylab="time (seconds)")
#' 
#' @export
#' 
pcv.emd<-function(df, cols=NULL, reorder=NULL, include=reorder, mat=F, plot = T, parallel = getOption("mc.cores",1), longTrait=NULL, id="image", value="value", raiseError=T){
  # df = df1; cols="ndvi_"; reorder=c("treatment", "genotype"); mat =F; plot=T; parallel = 1; include=reorder; longTrait=F; id="image";value="value"
  # df_long<-read.pcv(file, "long", F)
  # df = df_long; cols="index_frequencies_index_ndvi"; reorder=c("treatment", "genotype"); mat =F; plot=T; parallel = 1; include=reorder;
  # longTrait="trait"; id="image"; value="value"
  if(!is.null(longTrait)){traitCol = longTrait ; long=T}else{long=F}
  if(!is.null(reorder)){ df<-df[order(interaction(df[,reorder])),] }
  if(long){
    df<-df[grepl(cols, df[[traitCol]]), ]
    #* if nrow df is too high then do calculation for time per core based on estimates and report time with error.
    if(raiseError){
      eT_sec = 0.0025*((nrow(df)/parallel)^2)
      eT_min = eT_sec/60
      eT_hour = eT_min/60
      if(eT_sec <= 300){message(paste0("Estimated time of calculation is roughly ", round(eT_sec,1), " seconds using ", parallel, " cores in parallel."))
      } else if(eT_min < 60){warning(paste0("Estimated time of calculation is roughly ", round(eT_min,2), " minutes using ", parallel, " cores in parallel."))
          } else if(eT_min > 60){stop(paste0("Stopping, estimated time of calculation is roughly ", round(eT_hour,2), " hours using ", parallel, " cores in parallel.",
                                    "\nIf you wish to proceed then rerun this command with raiseError=F"))}
    }
    df$INNER_ID_EMD<-interaction(df[,id], drop=T)
    if(mat){ # make dist matrix
      out_data<-matrix(
        unlist(lapply(unique(df$INNER_ID_EMD), function(i){innerLapply(unique(df$INNER_ID_EMD), function(j){
          if(i==j){0}else{emd1d(as.numeric(df[df$INNER_ID_EMD==as.character(i), value]), as.numeric(df[df$INNER_ID_EMD==as.character(j), value]))}
        })})), nrow=length(unique(df$INNER_ID_EMD)), ncol = length(unique(df$INNER_ID_EMD)) )
    }else{ # make long data
      out_data<-do.call(rbind, lapply(unique(df$INNER_ID_EMD), function(i){
        do.call(rbind, parallel::mclapply(unique(df$INNER_ID_EMD), function(j){
          if(i==j){emdOut=0}else{emdOut=emd1d(as.numeric(df[df$INNER_ID_EMD==as.character(i), value]), as.numeric(df[df$INNER_ID_EMD==as.character(j), value]) )}
          if(!is.null(include)){x<-data.frame(i=i, j=j, emd = emdOut, df[df$INNER_ID_EMD==as.character(i), include], df[df$INNER_ID_EMD==as.character(j), include])
          colnames(x)<-c("i", "j", "emd", paste0(include,"_i"), paste0(include, "_j"))
          } else {x<-data.frame(i=i, j=j, emd = emdOut)}
          x
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
                                         "\nIf you wish to proceed then rerun this command with raiseError=F"))}
    }
      if(mat){# make dist matrix
        out_data<-matrix(
          unlist(lapply(1:nrow(df), function(i){innerLapply(1:nrow(df), function(j){
            if(i==j){0}else{emd1d(as.numeric(df[i, cols]), as.numeric(df[j, cols]))}
          })})), nrow=nrow(df), ncol = nrow(df) )
        if(!is.null(include)){
          rownames(out_data)<-interaction(df[,include])
        }
      }else{# make long dataframe
        out_data<-do.call(rbind, lapply(1:nrow(df), function(i){
          do.call(rbind, parallel::mclapply(1:nrow(df), function(j){
            if(i==j){emdOut=0}else{emdOut=emd1d(as.numeric(df[i, cols]), as.numeric(df[j, cols]))}
            if(!is.null(include)){x<-data.frame(i=i, j=j, emd = emdOut, df[i, include], df[j,include])
              colnames(x)<-c("i", "j", "emd", paste0(include,"_i"), paste0(include, "_j"))
              } else {x<-data.frame(i=i, j=j, emd = emdOut)}
            x
          }, mc.cores=parallel))
        }))
      }
  }
  if(plot){
    if(mat){
      p<-image(out_data)
    }else{
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
