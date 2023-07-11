#' Function to run a PCA, plot and optionally return the data with PCA coordinates and pca object
#' 
#' @param df Dataframe to ordinate
#' @param cols columns to reduce dimensions of. Can be specified with names or positions.
#' If this is length of 1 then it is treated as regex pattern to match
#' the column names that should be used.
#' @param color column name(s) used to color points in the pca plot.
#' @param trace Optional column to use to show changes over by way of geom_path,
#'  generally this would be a time variable.
#' @param facet Optional column or vector to facet/split plots on. 
#' If this is a character then it is taken as a column name. 
#' Numeric vectors are taken as values of trace to split on, 
#' if trace=NULL and this is numeric then plots will be split by the split argument.
#'  If this is a vector then a list of plots is returned.
#' @param returnData Logical, should data be returned?
#' @param ncp Optional, number of principal components to return attached
#'  to dataframe if data is returned. Defaults to all.
#' @param split Time variable to split data on if a facet on time is used.
#'  See bellwether vignette for examples.
#' @keywords pca
#' @details If data is returned then it will contain the coordinates from the
#'  PCA and will not contain the columns that were reduced.
#' 
#' @import ggplot2
#' @import FactoMineR
#' @importFrom stats as.formula
#' @examples
#' 
#' ## Not run:
#' 
#' hue_wide<-read.pcv(
#'   "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv",
#'    mode="wide", singleValueOnly = T, multiValPattern = "hist", reader="fread")
#' hue_wide$genotype = substr(hue_wide$barcode, 3,5)
#' hue_wide$genotype = ifelse(hue_wide$genotype == "002", "B73",
#'                            ifelse(hue_wide$genotype == "003", "W605S",
#'                                   ifelse(hue_wide$genotype == "004", "MM", "Mo17")))
#' hue_wide$fertilizer = substr(hue_wide$barcode, 8, 8)
#' hue_wide$fertilizer = ifelse(hue_wide$fertilizer == "A", "100",
#'                              ifelse(hue_wide$fertilizer == "B", "50", "0"))
#' hue_wide<-bw.time(hue_wide,timeCol="timestamp", group="barcode")
#' 
#' pcadf(hue_wide, cols = "hue_frequencies", color = "time_geno_fert", returnData=FALSE)
#' 
#' ## End(Not run)
#' 
#' @export

pcadf<-function(df=NULL, cols=NULL, color=NULL, trace=NULL,facet=NULL, returnData=TRUE, ncp=NULL, split="DAS"){
  if(is.character(cols) & length(cols)==1){
    cols<-which(grepl(cols,colnames(df)))
  }
  if(!is.null(color) && length(color)>1){
    df[[paste(color,collapse=".")]] = interaction(df[,color])
    color=paste(color,collapse=".")
  }
  if(is.null(ncp)){ncp=min(dim(df[,cols]))-1}
  pca<-FactoMineR::PCA(df[,cols], ncp=ncp, graph=FALSE)
  pc1Var<-round(pca$eig[1,2], 3)
  pc2Var<-round(pca$eig[2,2], 3)
  coords<-as.data.frame(pca$ind)
  coords<-coords[,grepl("coord", colnames(coords))]
  colnames(coords)<-gsub("coord.Dim.", "pc", colnames(coords))
  if(!is.numeric(cols)){ cols<-which(colnames(df) %in% cols) }
  pca.df<-cbind(df[, -cols], coords)
  if(!is.null(trace)){
    pca.df<-pca.df[ sort(pca.df[[trace]], index.return=TRUE)$ix ,]
  }
  traceSplit=FALSE
  FACET=FALSE
  additive=TRUE
  if(!is.null(facet)){
    if(is.character(facet)){
      FACET=TRUE
      facetLayer = ggplot2::facet_wrap(as.formula(paste0("~",paste(facet, collapse="+"))))
    }else if(is.numeric(facet) ){ # & !is.null(trace)
      traceSplit=TRUE
      traceSplits = facet
      traceDraw=TRUE
      if(is.null(trace)){
        trace=split
        traceDraw=FALSE
        additive=FALSE
      }
    }
  }
  
  if(is.null(color)){pca.df$dummyVariableForColor = 1; color="dummyVariableForColor"}
  
  if(traceSplit){
  
    traceLimit_x<-c( min(pca.df$pc1)-(min(pca.df$pc1)*.05), max(pca.df$pc1)+(max(pca.df$pc1)*.05))
    traceLimit_y<-c( min(pca.df$pc2)-(min(pca.df$pc2)*.05), max(pca.df$pc2)+(max(pca.df$pc2)*.05))
    
  plots<-lapply(1:length(traceSplits), function(i){
    maxTrace=traceSplits[i]
    if(i==1){prevTrace = min(pca.df[[trace]])}else{prevTrace = traceSplits[i-1]}
    if(additive){
      if(i==length(traceSplits)){through="+"
      }else{through = paste0(" through ",traceSplits[i+1])}
      TITLE = paste0(trace, ": ", maxTrace, through)
    } else{TITLE = paste(trace, maxTrace)}
    
    if(additive){ pca.df.sub<-pca.df[pca.df[[trace]]<=maxTrace, ] 
      } else{ pca.df.sub<-pca.df[pca.df[[trace]]==maxTrace, ] }
    
    p<-ggplot2::ggplot(pca.df.sub,
                           ggplot2::aes(x=.data$pc1, y=.data$pc2, color = .data[[color]]))+
      ggplot2::geom_point()+
      ggplot2::labs(x=paste0("PC 1 (",pc1Var,"%)"),y=paste0("PC 2 (",pc2Var,"%)"),
                    title = TITLE)+
      pcv_theme()+
      coord_cartesian(xlim = traceLimit_x, ylim=traceLimit_y)
    if(color=="dummyVariableForColor"){
      p<-p+ggplot2::theme(legend.position="none")
    }
    if(!is.null(trace) & traceDraw){
      p<-p+ggplot2::geom_path(data = pca.df[pca.df[[trace]]>=prevTrace & pca.df[[trace]]<=maxTrace, ], ggplot2::aes(group=.data[[color]]),linewidth=0.15, show.legend = F)
    }
    if(FACET){
      p<-p+facetLayer
    }
    return(p)
  })
  
  } else{
    plots<-ggplot2::ggplot(pca.df, ggplot2::aes(x=.data$pc1, y=.data$pc2, color = .data[[color]]))+
      ggplot2::geom_point()+
      ggplot2::labs(x=paste0("PC 1 (",pc1Var,"%)"),y=paste0("PC 2 (",pc2Var,"%)"))+
      pcv_theme()
    if(color=="dummyVariableForColor"){
      plots<-plots+ggplot2::theme(legend.position="none")
    }
    if(!is.null(trace)){
      plots<-plots+ggplot2::geom_path(ggplot2::aes(group=.data[[color]]),linewidth=0.15, show.legend = F)
    }
    if(FACET){
      plots<-plots+facetLayer
    }
  }
  
  if(returnData){return(list("data"=pca.df, "pca"=pca, "plot"=plots))}else{return(plots)}
}
