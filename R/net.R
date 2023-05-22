#' Network analysis of a distance matrix
#' 
#' @description Easy igraph use with pcv.emd output
#' 
#' 
#' @param emd A long dataframe as returned by pcv.emd. Currently this function is only made to work with dataframe output, not distance matrix output.
#' @param meta Metadata to be carried from pcv.emd output into the network, defaults to NULL which will use all metadata.
#' @param dissim Logical, should the distCol be inverted to make a dissimilarity value?
#' @param distCol The name of the column containing distances/dissimilarities. Defaults to "emd" for compatability with pcv.emd
#' @param ... Further arguments passed to igraph::graph_from_dataframe or graph_from_adjacency_matrix
#' @import ggplot2
#' @import igraph
#' 
#' @keywords emd, earth mover's distance, multi-value trait, network
#' @examples 
#' 
#' file = "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest1.csv"
#' df1<-read.pcv(file, "wide", T, multiValPattern = c("index_frequencies_index_ari", "index_frequencies_index_ci_rededge", "npq_hist_NPQ", "yii_hist_Fq'/Fm'", "yii_hist_Fv/Fm"))
#' colnames(df1)<-sub("index_frequencies_index_ndvi.", "ndvi_", colnames(df1))
#' emd_df<-pcv.emd(df1, cols="ndvi_", reorder=c("treatment", "genotype", "X"), mat =F, plot=F, parallel = 1)
#' net<-pcv.net(emd_df, meta = c("treatment", "genotype"))
#' 
#' @export
#' 

pcv.net<-function(emd = NULL, meta = NULL, dissim=T, distCol="emd"){
  #* emd = emd_df; meta = c("treatment", "genotype"); dissim=T; distCol="emd"
  #* format distCol into a similarity col
  if(is.data.frame(emd)){
    if(dissim){
      emd[[distCol]]<-1/emd[[distCol]]
      emd[[distCol]]<-ifelse(is.infinite(emd[[distCol]])|is.na(emd[[distCol]]), 0, emd[[distCol]])
    }
    #* turn long data into a graph and extract nodes/edges
    g<-igraph::graph_from_data_frame(emd, directed=F)
  } 
  gg<-as.data.frame(igraph::layout.auto(g))
  eg<-igraph::get.data.frame(g)
  #* link metadata to nodes
  gg$index <- 1:nrow(gg)
  metaIndex<-lapply(meta, function(m) which(grepl(m, colnames(eg)))[1])
  newCols<-(ncol(gg)+1):(ncol(gg)+length(meta))
  gg[,newCols] <- lapply(metaIndex, function(i) eg[[i]][match(gg$index, eg$from)])
  colnames(gg)[newCols]<-meta
  #* Calculate network metrics
  gg$betweenness<-igraph::betweenness(g)
  gg$degree<-igraph::degree(g)
  igraph::E(g)$weight<-eg[[distCol]]+0.1
  gg$strength<-igraph::strength(g)
  gg$harmonic_centrality<-igraph::harmonic_centrality(g)
  gg$eigen_centrality<-igraph::eigen_centrality(g)[[1]]
  gg$authority_score<-igraph::authority_score(g)[[1]]
  gg$page_rank<-igraph::page_rank(g)[[1]]
  #* add coordinates for plotting edges
  eg$from.x <- gg$V1[match(eg$from, gg$index)]  
  eg$from.y <- gg$V2[match(eg$from, gg$index)]
  eg$to.x <- gg$V1[match(eg$to, gg$index)] 
  eg$to.y <- gg$V2[match(eg$to, gg$index)]
  return(list("nodes" = gg, "edges" = eg, "graph" = g))
}

