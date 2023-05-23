#' Earth Mover's Distance between spectral histograms
#' 
#' @description emd1d computes 1 dimension EMD for two samples.
#' 
#' @param s1 Histogram as a numeric vector of counts per position.
#' @param s2 Histogram as a numeric vector of counts per position. Must be the same length as s1.
#' 
#' @import ggplot2
#' 
#' @keywords emd, earth mover's distance, multi-value trait, histogram
#' @return Returns EMD between two samples as a numeric.
#' @examples 
#' 
#' set.seed(123) 
#' s1<-hist(rnorm(10000,50,10), breaks=seq(1,100,1))$counts
#' s2<-hist(rlnorm(9000,log(30),0.25), breaks=seq(1,100,1))$counts
#' plot(s2,type="l"); lines(s1)
#' emd1d(s1,s2)
#' 
#' 
#' @export
#' 
#' 
emd1d<-function(s1, s2){
  if(length(s1)!=length(s2)){stop("Samples must be from the same histogram and be of the same length")}
  s1<-s1/sum(s1)
  s2<-s2/sum(s2)
  emd_iter<-numeric(length(s1)+1)
  for(i in 1:length(s1)){ emd_iter[i+1]<-s1[i]-s2[i]+emd_iter[i] }
  return(sum(abs(emd_iter)))
}
