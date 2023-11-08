#' @description
#' Negative binomial via conjugate beta prior GIVEN A KNOWN R
#' conjugacy: 
#' if counts ~ nbinom(n, r)
#' where n is the number of successful trials
#'   and r is the prob of success per trial which is fixed and known
#' 
#' if r is known then p ~ beta(A, B)
#' note that the compendium of conjugate priors seems to have a typo for this relationship,
#' but the wikipedia conjugate prior article is correct
#' A' = A + r*n
#' B' = B + sum(x)
#'
#' Using MoM:
#' 
#' \bar{x} = k(1-p)/p
#' s^2 = \bar{x}/p
#' r = \bar{x}^2 / (s^2 - \bar{x})
#' p = \bar{x}/s^2
#' 
#' @param s1 A vector of numerics drawn from a negative binomial distribution.
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' if(FALSE){
#' set.seed(123)
#' true_p1 = 0.7
#' true_p2 = 0.25
#' true_r1 = 10
#' true_r2 = 10
#' s1<-rnbinom(10, true_r1, true_p1)
#' s2<-rnbinom(10, true_r2, true_p2)
#' priors = list(a=c(0.5, 0.5), b=c(0.5, 0.5))
#' s1
#' s2
#' xbar1<-(true_r1*(1-true_p1))/true_p1 # MoM Mean1: 10
#' xbar2<-(true_r2*(1-true_p2))/true_p2 # MoM Mean2: 30
#' 
#' var1 <- (true_r1*(1-true_p1))/true_p1/true_p1 # MoM var1: 20
#' var2 <- (true_r2*(1-true_p2))/true_p2/true_p2 # MoM var2: 120
#' 
#' (obs_xbar1 = mean(s1))
#' (obs_var1 = var(s1))
#' (f1 = xbar1^2 / (var1 - xbar1)) # number of failures in getting r trials
#' (p1 = xbar1/var1)
#' (a1_prime <- 0.5 + true_r1 * length(s1))
#' (b1_prime_wiki <- 0.5 + sum(s1) )
#' # (b1_prime_comp <- 0.5 + sum(s1) - (true_r1 * length(s1))) # -0.5
#' 
#' plot(density(rbeta(1000,a1_prime, b1_prime_wiki)), xlim=c(0,1))
#' # lines(density(rbeta(1000,a1_prime, b1_prime_comp)), xlim=c(0,1), col="blue")
#' abline(v=p1, col="red")
#' 
#' .conj_negbin_sv(s1 = rnbinom(10, 10, 0.5), s2= rnbinom(100, 10, 0.25), priors = list(r = c(10,10), a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'               cred.int.level = 0.89, hypothesis="equal")
#' }
#' @keywords internal
#' @noRd

.conj_negbin_sv<-function(s1 = NULL, priors = NULL,
                          plot=FALSE, support=NULL, cred.int.level = NULL){
  
  #* `Check samples`
  if(any(abs(c(s1,s2)-round(c(s1,s2)))>.Machine$double.eps^0.5) | any(c(s1,s2)<0) ){stop("Only positive whole numbers can be used in the Negative Binomial distribution")}
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list(r=c(10, 10), a=c(0.5,0.5),b=c(0.5,0.5)) # beta prior on P
    warning("True value of r for negative binomial distribution has defaulted to 10, you should add a prior including r parameter.")
  }
  
  out <- list()
  
  #* `Use conjugate beta prior on probability`
  #* Note that this is very sensitive to the R value being appropriate
  a1_prime <- priors$a[1] + priors$r[1] * length(s1)
  b1_prime <- priors$b[1] + sum(s1) 
  #* `Define support if it is missing`
  if(is.null(support)){
    support<- seq(0.0001,0.9999, 0.0001)
  }
  #* `calculate density over support``
  dens1 <- dbeta(support, a1_prime, b1_prime)
  pdf1 <- dens1/sum(dens1)
  
  #* `calculate highest density interval`
  hdi1 <- qbeta(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),a1_prime,b1_prime)
  
  #* `calculate highest density estimate``
  if(a1_prime <= 1 & b1_prime > 1){
    hde1 <- 0
  }else if(a1_prime > 1 & b1_prime <= 1){
    hde1 <- 1
  }else{
    hde1 <- (a1_prime-1)/(a1_prime+b1_prime-2) # ~ p1 ~ xbar1/var1 
  }
  
  #* `save summary and parameters`
  out$summary <- data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$r <- priors$r[1]
  out$posterior$a <- a1_prime
  out$posterior$b <- b1_prime
  #* `Make Posterior Draws`
  out$posteriorDraws = rbeta(10000,a1_prime,b1_prime)
  out$pdf <- pdf1
  #* `keep data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=support, "prob"=pdf1, "sample"=rep("Sample 1",length(support) ))
  }
  return(out)
}
