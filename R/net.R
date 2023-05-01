#' Network analysis of a distance matrix
#' 
#' @description Easy igraph use with pcv.emd output
#' 
#' This should take a long dataframe or a distance matrix from pcv.emd, get similarity from dissimilarity,
#' then make a network and plot it, returning a ggplot and a list of dataframes representing nodes and edges.
#' 
#' 
#' @param emd
#' @param meta Metadata to be carried from pcv.emd output into the network, defaults to NULL which will use all metadata.
#' @param dissim Logical, should the distCol be inverted to make a dissimilarity value?
#' @param distCol The name of the column containing distances/dissimilarities. Defaults to "emd" for compatability with pcv.emd
#' @param ... Further arguments passed to igraph::graph_from_dataframe or graph_from_adjacency_matrix
#' @import ggplot2
#' @import igraph
#' 
#' @keywords emd, earth mover's distance, multi-value trait, histogram
#' @examples 
#' 
#' file = "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest1.csv"
#' df1<-read.pcv(file, "wide", T, multiValPattern = c("index_frequencies_index_ari", "index_frequencies_index_ci_rededge", "npq_hist_NPQ", "yii_hist_Fq'/Fm'", "yii_hist_Fv/Fm"))
#' colnames(df1)<-sub("index_frequencies_index_ndvi.", "ndvi_", colnames(df1))
#' emd<-pcv.emd(df1, cols="ndvi_", reorder=c("treatment", "genotype", "X"), mat =F, plot=F, parallel = 1)
#' emdMat<-pcv.emd(df1, cols="ndvi_", reorder=c("treatment", "genotype", "X"), mat =T, plot=F, parallel = 1)
#' 
#' @export
#' 


pcv.net<-function(emd = NULL, meta = NULL, dissim=F, distCol="emd"){}


emd$dissim <- 1/emd$emd
emd$dissim<-ifelse(is.infinite(emd$dissim), NA, emd$dissim)
g<-igraph::graph_from_data_frame(emd, directed=F)
gg<-as.data.frame(layout.auto(g))

eg<-igraph::get.data.frame(g)
gg$index <- 1:nrow(gg)
gg$trt <- eg$Treatment_i[match(gg$index, eg$from)]  
gg$geno <- eg$Genotype_i[match(gg$index, eg$from)]  

eg$from.x <- gg$V1[match(eg$from, gg$index)]  
eg$from.y <- gg$V2[match(eg$from, gg$index)]
eg$to.x <- gg$V1[match(eg$to, gg$index)] 
eg$to.y <- gg$V2[match(eg$to, gg$index)]


ggplot(gg)+
  geom_segment(data=eg,aes(x=from.x, xend = to.x, y=from.y, yend = to.y, linewidth=dissim),colour="black",alpha=0.1) +
  geom_point(data=gg,size=3,shape=21,color="black",aes(x=V1,y=V2, fill = geno),alpha=1)+
  theme_void()






