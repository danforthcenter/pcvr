#' @description
#' Internal function for calculating posterior distribution of \Lambda parameter for poisson data.
#' The conjugate prior on Lambda is a Gamma(A,B)
#' A=B=0.5 is a reasonable weak default prior
#' 
#' So if leaf count is Poisson distributed:
#' count ~ Pois(\labmda)
#' \labmda ~ gamma(A, B)
#' A = A_[prior] + sum(x)
#' B = B_[prior] / (1+n)
#' 
#' via MoM \hat(\labmda) =  = 1/n +sum^1_n(x)
#' @param s1 A vector of numerics drawn from a beta distribution.
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' if(FALSE){
#' .conj_poisson_sv(s1 = rpois(20, 10), priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=FALSE)
#' }
#' @keywords internal
#' @noRd

.conj_poisson_sv<-function(s1 = NULL, priors = NULL,
                           plot=FALSE, support=NULL, cred.int.level = NULL){
  #* `Check samples`
  if(any(abs(c(s1,s2)-round(c(s1,s2)))>.Machine$double.eps^0.5) | any(c(s1,s2)<0) ){stop("Only positive whole numbers can be used in the Poisson distribution")}
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list(a=c(0.5,0.5),b=c(0.5,0.5)) # gamma prior on lambda
  }
  #* `Define dense Support if missing`
  if(is.null(support)){
    upper <- max(c(s1,s2))*2
    support <-seq(0, upper, 0.01)
  }
  out <- list()
  
  #* `Use conjugate gamma prior on lambda`
  a1_prime <- priors$a[1] + sum(s1)
  b1_prime <- priors$b[1] + length(s1)
  
  #* `calculate density over support``
  dens1 <- dgamma(support, a1_prime, b1_prime)
  pdf1 <- dens1/sum(dens1)
  
  #* `calculate highest density interval`
  hdi1 <- qgamma(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),a1_prime,b1_prime)
  
  #* `calculate highest density estimate``
  if(a1_prime <= 1 & b1_prime > 1){
    hde1 <- 0
  }else if(a1_prime > 1 & b1_prime <= 1){
    hde1 <- Inf
  }else{
    hde1 <- (a1_prime-1)/b1_prime 
  }
  
  #* `save summary and parameters`
  out$summary <- data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a1_prime
  out$posterior$b <- b1_prime
  #* `Make Posterior Draws`
  out$posteriorDraws =rgamma(10000,a1_prime,b1_prime)
  out$pdf <- pdf1
  #* `keep data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=support, "prob"=pdf1, "sample"=rep("Sample 1",length(support) ))
  }
  return(out)
}



