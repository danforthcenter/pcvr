#' @description
#' Internal function for Bayesian T Tests of gaussian data represented by single value traits.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @examples
#' if(FALSE){
#' .conj_gaussian_means_sv(s1=rnorm(100, 50,10), s2= rnorm(100, 60,12),
#'                         priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                         plot=FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'                         cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' }
#' @keywords internal
#' @noRd
.conj_gaussian_means_sv<-function(s1 = NULL,  priors = NULL,
                                  plot=FALSE, support=NULL, cred.int.level=NULL){
  out <- list()
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list( mu=c(0,0),n=c(1,1),s2=c(20,20) )
  }
  #* `Define support if it is missing`
  if(is.null(support)){
    support<-seq(floor(min(c(s1,s2))*0.8), ceiling(max(c(s1,s2))/0.8), length.out=10000)
  }
  #* `Get Mean, Variance, SE, and DF from s1`
  if(length(s1) > 1){
    n1 <- length(s1) # n samples
    m1 = mean(s1) # xbar
    s2_1 = var(s1) # s^2
    
    v1 = priors$n[1] - 1 # prior DF
    n1_n = priors$n[1] + n1 # total N including prior
    m1_n = (n1*m1 + priors$n[1]*priors$mu[1])/n1_n # weighted mean of prior and data
    v1_n = v1 + n1 # degrees of freedom including data
    s2_1_n = ((n1-1)*s2_1 + v1*priors$s2[1] + priors$n[1]*n1*(priors$mu[1] - m1)^2/n1_n)/v1_n # pooled variance
    se1<-sqrt(s2_1_n/n1_n) # standard error of the mean
    
    dens <- extraDistr::dlst(support,v1_n,m1_n, se1)
    pdf1 <- dens/sum(dens)
    hde1_mean <- m1_n
    hdi1_mean <- m1_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v1_n)*se1
    
    out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
    out$posterior$mu <- m1_n
    out$posterior$n <- n1_n
    out$posterior$s2 <- s2_1_n # return variance
    #* `Make Posterior Draws`
    out$posteriorDraws = extraDistr::rlst(10000, v1_n, m1_n, se1)
    out$pdf <- pdf1
    #* `Save data for plotting`
    if(plot){
      out$plot_df <- data.frame("range"=support,
                                "prob"=pdf1,
                                "sample"=rep("Sample 1",length(support))
      )
    } 
  }else{
    stop("s1 must be a numeric of length 2 or greater")
  }
 
  return(out)
}





#' @description
#' Internal function for calculating \alpha and \beta of a distribution represented by multi value traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number representing the "bin".
#' @examples
#' if(FALSE){
#' makeMvGauss<-function(bins=180,mu,sigma){
#'    setNames(data.frame(matrix(hist(rnorm(2000,mu, sigma),
#'    breaks=seq(1,bins,1), plot=FALSE)$counts, nrow=1)),paste0("b",1:(bins-1) ))
#'    }
#' mv_gauss<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=50, sigma=10 )})),
#'                 do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=60, sigma=12 )})))
#' .conj_gaussian_means_mv(s1 = mv_gauss[1:30,], s2= mv_gauss[31:60,],
#'                         priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                         plot=FALSE, rope_range = c(-0.1,0.1), rope_ci = 0.89, 
#'                         cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' }
#' @keywords internal
#' @noRd

.conj_gaussian_means_mv<-function(s1 = NULL, priors = NULL,
                                  plot=FALSE, support=NULL, cred.int.level=NULL){
  out <- list()
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list( mu=c(0,0),n=c(1,1),s2=c(20,20) )
  }
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
  bins_order<-sort(histCols_bin, index.return=TRUE)$ix
  s1<-s1[,bins_order]
  
  #* `Define support if it is missing`
  if(is.null(support)){
    support<-seq(from = min(bins_order), to = max(bins_order), length.out=10000)
  }
  
  #* `Turn s1 matrix into a vector`
  X1<-rep(histCols_bin[bins_order], as.numeric(round(colSums(s1))) )
  
  #* `Get Mean, Variance, SE, and DF from s2`
  n1 <- nrow(s1) # n samples
  m1 = mean(X1) # xbar
  s2_1 = var(X1) # s^2
  
  v1 = priors$n[1] - 1 # prior DF
  n1_n = priors$n[1] + n1 # total N including prior
  m1_n = (n1*m1 + priors$n[1]*priors$mu[1])/n1_n # weighted mean of prior and data
  v1_n = v1 + n1 # degrees of freedom including data
  s2_1_n = ((n1-1)*s2_1 + v1*priors$s2[1] + priors$n[1]*n1*(priors$mu[1] - m1)^2/n1_n)/v1_n # pooled variance
  se1<-sqrt(s2_1_n/n1_n) # standard error of the mean
  
  dens <- extraDistr::dlst(support,v1_n,m1_n, se1)
  pdf1 <- dens/sum(dens)
  hde1_mean <- m1_n
  hdi1_mean <- m1_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v1_n)*se1
  
  out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- m1_n
  out$posterior$n <- n1_n
  out$posterior$s2 <- s2_1_n
  #* `Make Posterior Draws`
  out$posteriorDraws = extraDistr::rlst(10000, v1_n, m1_n, se1)
  out$pdf <- pdf1
  #* `Save data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=support,
                              "prob"=pdf1,
                              "sample"=rep("Sample 1",length(support))
    )
  }
  
  return(out)
  
}

