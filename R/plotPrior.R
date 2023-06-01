#' Check priors used in ease of use brms functions
#' 
#' @param priors A named list of means for prior distributions. This does not support vectors for each element of the list or brmsprior objects, but otherwise should take the same input as \code{\link{growthSS}}.
#' @keywords Bayesian, brms
#' @return A named list of plots showing prior distributions that \code{growthSS} would make with these inputs.
#' @import ggplot2
#' @examples 
#' 
#' priors = list("A" = 130, "B" = 10, "C" = 1.2)
#' priorCheck(priors)
#' 
#' @export

plotPrior<-function(priors){
  plots<-lapply(1:length(priors), function(i){
    pri=priors[[i]]
    nm=names(priors)[i]
    max=round(max(rlnorm(1000, log(pri), 0.25))*1.1, 0)
    ggplot2::ggplot(data = data.frame(x=seq(0,max,length.out=10000)), ggplot2::aes(x=x))+
      ggplot2::stat_function(fun=dlnorm, geom="polygon", fill="gray40",color="black",
                             args=list(meanlog = log(pri), sdlog=0.25, log=F))+
      ggplot2::theme_minimal()+ggplot2::labs(y="Density", title=nm)
  })
  names(plots)<-names(priors)
  return(plots)
}
