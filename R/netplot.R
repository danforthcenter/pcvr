#' Visualizing igraph networks
#' 
#' @description Easy igraph visualization with pcv.net output
#' 
#' 
#' @param net Network object similar to that returned from pcv.net, having dataframes named "edges" and "nodes" 
#' @param fill Variable name(s) from the nodes data to be used to color points. By default "strength" is used.
#' @param shape Optional discrete variable name(s) from the nodes data to be used to change the shape of points.
#' @param size Size of points, defaults to 3.
#' @param edgeWeight Edge dataframe column to weight connections between nodes. Defaults to "emd" for compatability with \code{pcv.emd}.
#' @param edgeFilter How should edges be filtered? This can be either a numeric (0.5) in which case it is taken as a filter where only edges with values greater than or equal to that number are kept or a character string ("0.5") in which case the strongest X percentage of edges are kept. This defaults to NULL which does no filtering, although that should not be considered the best standard behaviour. See details.
#' @import ggplot2
#' 
#' @keywords emd, earth mover's distance, multi-value trait, network
#' @examples 
#' 
#' file = "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest1.csv"
#' df1<-read.pcv(file, "wide", T, multiValPattern = c("index_frequencies_index_ari", "index_frequencies_index_ci_rededge", "npq_hist_NPQ", "yii_hist_Fq'/Fm'", "yii_hist_Fv/Fm"))
#' colnames(df1)<-sub("index_frequencies_index_ndvi.", "ndvi_", colnames(df1))
#' emd_df<-pcv.emd(df1, cols="ndvi_", reorder=c("treatment", "genotype", "X"), mat =F, plot=F, parallel = 1)
#' network<-pcv.net(emd_df, meta = c("treatment", "genotype"))
#' net=network; fill = "strength"; shape = "genotype"; size=5; edgeWeight="emd"
#' net.plot(network, fill = "strength", shape = "genotype", size=5, edgeFilter=0.5)
#' net.plot(network, fill = "strength", shape = "genotype", size=5, edgeFilter="0.5")
#' 
#' @return Returns a ggplot of a network.
#' 
#' @export
#' 

net.plot<-function(net, fill="strength", shape=NULL, size = 3, edgeWeight="emd", edgeFilter = NULL){
  nodes<-net[["nodes"]]
  edges<-net[["edges"]]
  if(is.null(fill)){fill="NOFILL"; edges$NOFILL="a"}
  if(length(fill)>1){
    edges$FILL = interaction(edges[,fill])
    fill="FILL"
  }
  if(is.null(shape)){shape = "NOSHAPE"; nodes$NOSHAPE="a"}
  if(length(shape)>1){
    nodes$SHAPE = interaction(nodes[,shape])
    shape="SHAPE"
  }
  if(!is.null(edgeFilter)){
    if(is.character(edgeFilter)){
      cutoff<-quantile(edges[[edgeWeight]], probs = as.numeric(edgeFilter))
      edges<-edges[edges[[edgeWeight]] >= as.numeric(cutoff), ]
    } else if(is.numeric(edgeFilter)){
      edges<-edges[edges[[edgeWeight]] >= edgeFilter, ]
    } else{stop("edgeFilter must be character or numeric, see ?net.plot for details.")}
    nodes<-nodes[nodes$index %in% c(edges$from, edges$to),]
    }
  
  
  p<-ggplot2::ggplot(nodes)+
    ggplot2::geom_segment(data=edges,ggplot2::aes(x=from.x, xend = to.x, y=from.y, yend = to.y, linewidth=.data[[edgeWeight]]),colour="black",alpha=0.1) +
    ggplot2::geom_point(data=nodes, size=size, ggplot2::aes(x=V1,y=V2, fill = .data[[fill]], color=.data[[fill]], shape=.data[[shape]]), alpha=1, show.legend=T)+
    ggplot2::scale_linewidth(range=c(0.1,1.5))+
    #* note that scaling shape should work, but there is a documented ggplot2 bug where this messes up the legend, so 
    #* until that is fixed I will not specify fillable shapes.
    #scale_shape_discrete(breaks=c(21:(21-1+length(unique(nodes[[shape]]))) ), guide="legend", position="bottom")+
    ggplot2::guides(linewidth="none", shape=ggplot2::guide_legend(nrow=1), fill="none")+
    ggplot2::theme_void()+
    ggplot2::theme(legend.position="bottom")
  if(length(fill==1) && fill=="NOFILL"){ p<-p+ggplot2::guides(color="none") }
  if(length(shape)==1 && shape=="NOSHAPE"){ p<-p+ggplot2::guides(shape="none") }
  return(p)
}




