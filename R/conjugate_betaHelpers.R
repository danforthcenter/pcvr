#' @description
#' Internal function for calculating \alpha and \beta of a distribution represented by multi value traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number between 0.0001 and 0.9999 representing the "bin".
#' @examples
#' if(FALSE){
#' makeMvBeta<-function(n=100,a,b){
#'   setNames(data.frame(matrix(hist(rbeta(2000,a,b),
#'   breaks=seq(0,1,length.out=n), plot=FALSE)$counts, nrow=1)),paste0("b0.",1:(n-1)))
#' }
#' 
#' mv_beta<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvBeta(n=100, a=5, b=8 )})),
#'                do.call(rbind, lapply(1:30, function(i){makeMvBeta(n=100, a=10, b=3 )})))
#' 
#' .conj_beta_mv(s1 = mv_beta[1:30,], priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=FALSE)
#' }
#' @keywords internal
#' @noRd
.conj_beta_mv<-function(s1 = NULL, priors = NULL,
                        plot=FALSE, support = NULL, cred.int.level = NULL,
                        calculatingSupport=FALSE){
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list(a=0.5, b=0.5)
  }
  #* `Define dense Support`
  
  if(is.null(support)){
    if(calculatingSupport){return(c(0.0001, 0.9999))}
    support <-seq(0.0001, 0.9999, 0.0001)
  }
  
  out <- list()
  #* `Standardize sample 1 class and names`
  if(is.null(colnames(s1))){
    bins<-(1:ncol(s1))/100
    colnames(s1)<-paste0("b", bins )
    warning(paste0("Assuming unnamed columns represent bins from ", min(bins), " to ", max(bins)))
  }
  if(is.matrix(s1)){ s1<-as.data.frame(s1) }
  
  #* `Reorder columns if they are not in the numeric order`
  histCols<-colnames(s1)
  histCols_bin <- as.numeric(sub("[a-zA-Z_.]+","",colnames(s1)))
  if(any(histCols_bin>1)){
    histCols_bin<-histCols_bin/100
  }
  bins_order<-sort(histCols_bin, index.return=TRUE)$ix
  s1<-s1[,bins_order]
  
  #* `Turn matrix into a vector`
  X1<-rep(histCols_bin[bins_order], as.numeric(round(colSums(s1))) )
  
  #* `get parameters for s1 using method of moments``
  #* y ~ Beta(\alpha, \beta)
  #* \alpha ~ \bar{y}( ( (\bar{y} * (1-\bar{y}))/\bar(var) )-1 )
  #* \beta ~ (1-\bar{y})( ( (\bar{y} * (1-\bar{y}))/\bar(var) )-1 )
  mu1 <- mean(X1) #' \bar{y}
  nu1 <- var(X1)/(nrow(s1)-1) #' \bar{var} the unbiased sample variance
  alpha1 <- mu1*((mu1*(1-mu1))/(nu1) - 1)
  beta1 <- (1-mu1)*((mu1*(1-mu1))/(nu1) - 1)
  
  #* `Add priors`
  a1_prime <- alpha1 + priors$a[1]
  b1_prime <- beta1 + priors$b[1]
  
  #* `calculate density`
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
    hde1 <- (a1_prime-1)/(a1_prime+b1_prime-2)
  }
  
  #* `save summary and parameters`
  out$summary <- data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a1_prime
  out$posterior$b <- b1_prime
  #* `Make Posterior Draws`
  out$posteriorDraws = rbeta(10000, a1_prime, b1_prime)
  out$pdf <- pdf1
  #* `keep data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=support, "prob"=pdf1, "sample"=rep("Sample 1",length(support) ))
  }
  
  return(out)
}



#' @description
#' Internal function for calculating \alpha and \beta of a distribution represented by multi value traits.
#' @param s1 A vector of numerics drawn from a beta distribution.
#' @examples
#' if(FALSE){
#' .conj_beta_sv(s1 = rbeta(100, 5, 10), priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=FALSE)
#' }
#' @keywords internal
#' @noRd
.conj_beta_sv<-function(s1 = NULL, priors = NULL,
                        plot=FALSE, support=NULL, cred.int.level = NULL,
                        calculatingSupport=FALSE){
  
  if(any(c(s1, s2)>1)){stop("Values above 1 cannot be used with the beta distribution")}
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list(a=0.5, b=0.5)
  }
  #* `Define dense Support`
  if(is.null(support)){
    if(calculatingSupport){return(c(0.0001, 0.9999))}
    support <-seq(0.0001, 0.9999, 0.0001)
  }
  out <- list()
  
  #* `get parameters for s1 using method of moments``
  #* y ~ Beta(\alpha, \beta)
  #* \alpha ~ \bar{y}( ( (\bar{y} * (1-\bar{y}))/\bar(var) )-1 )
  #* \beta ~ (1-\bar{y})( ( (\bar{y} * (1-\bar{y}))/\bar(var) )-1 )
  mu1 <- mean(s1) #' \bar{y}
  nu1 <- var(s1)/(length(s1)-1) #' \bar{var} the unbiased sample variance
  alpha1 <- mu1*((mu1*(1-mu1))/(nu1) - 1)
  beta1 <- (1-mu1)*((mu1*(1-mu1))/(nu1) - 1)
  
  #* `Add priors in`
  a1_prime <- priors$a[1] + alpha1
  b1_prime <- priors$b[1] + beta1
  
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
    hde1 <- (a1_prime-1)/(a1_prime+b1_prime-2)
  }
  
  #* `save summary and parameters`
  out$summary <- data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
  out$posterior$a <- a1_prime
  out$posterior$b <- b1_prime
  #* `Make Posterior Draws`
  out$posteriorDraws = rbeta(10000, a1_prime, b1_prime)
  out$pdf <- pdf1
  #* `keep data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=support, "prob"=pdf1, "sample"=rep("Sample 1",length(support) ))
  }
  return(out)
}


