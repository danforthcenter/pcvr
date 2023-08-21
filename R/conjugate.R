#' Bayesian testing using conjugate priors and method of moments for single or multi value traits.
#' 
#' @description
#' Function to perform bayesian tests and ROPE comparisons using single or multi value traits with several distributions.
#' 
#' 
#' @param s1 A data.frame or matrix of multi value traits or a vector of single value traits.
#' If a multi value trait is used then column names should include a number representing the "bin".
#' @param s2 An optional second sample of the same form as s2.
#' @param method The distribution/method to use.
#' Currently "t", "gaussian", "beta", "lognormal", "poisson", "negbin" (negative binomial), and "dirichlet"
#' are supported. The count distributions (poisson and negative binomial) 
#' are only implemented for single value traits. Dirichlet is only implemented with multi value traits and
#' will return the summary component as a data.table with nested data frames for the ROPE HDI/HDE
#' if rope_range is specified. For brevity the HDE and HDI of the samples are not returned when
#' using the dirichlet method.
#' Note that "t" and "gaussian" both use a T distribution with "t" testing for a difference
#'  of means and "gaussian" testing for a difference in the distributions (similar to a Z test).
#' @param priors Prior distributions described as a list. This varies by method, see details.
#'  By default this is NULL and weak priors (jeffrey's prior where appropriate) are used.
#'  The \code{posterior} part of output can also be recycled as a new prior if Bayesian
#'  updating is appropriate for your use.
#' @param plot Logical, should a ggplot be made and returned.
#' @param rope_range Optional vector specifying a region of practical equivalence.
#' This interval is considered practically equivalent to no effect.
#' Kruschke (2018) suggests c(-0.1, 0.1) as a broadly reasonable ROPE for standardized parameters.
#' That range could also be rescaled by a standard deviation/magnitude for
#' non-standardized parameters, but ultimately this should be informed by your
#' setting and scientific question.
#' See Kruschke (2018) for details on ROPE and other Bayesian methods to aide
#' decision-making \url{https://doi.org/10.3758/s13423-016-1221-4}
#' and \url{https://doi.org/10.1177/251524591877130}.
#' @param rope_ci The credible interval probability to use for ROPE. Defaults to 0.89. 
#' @param cred.int.level The credible interval probability to use 
#' in computing HDI for samples, defaults to 0.89.
#' @param hypothesis Direction of a hypothesis if two samples are provided.
#'  Options are "unequal", "equal", "greater", and "lesser",
#'   read as "sample1 greater than sample2".
#' @param support Optional support vector to include all possible values the random variable
#'  (samples) might take. This defaults to NULL in which case each method will use default
#'  behavior to attempt to calculate a dense support, but it is a good idea to supply this 
#'  with some suitable vector. For example, the Beta method uses \code{seq(0.0001, 0.9999, 0.0001)}
#'  for support. 
#' 
#' @import bayestestR
#' @import ggplot2
#' @import patchwork
#' @import extraDistr
#' @importFrom stats var median rnorm dnorm qbeta rbeta dbeta qbeta rlnorm dlnorm qlnorm qt rgamma qgamma dgamma rpois
#' 
#' @details 
#' 
#' Prior distributions default to be weakly informative and in some cases you may wish to change them. 
#' \itemize{
#'    \item{\strong{"t" and "gaussian":} \code{priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ) },
#'     where mu is the mean, n is the number of prior observations, and s2 is variance}
#'    \item{\strong{"beta":} \code{priors = list( a=c(0.5, 0.5), b=c(0.5, 0.5) )},
#'     where a and b are shape parameters of the beta distribution.}
#'    \item{\strong{"lognormal": } \code{priors = list( mu_log=c(log(10),log(10)),n=c(1,1),sigma_log=c(log(3),log(3)) ) },
#'     where mu_log is the mean on log scale, n is the number of prior observations, and sigma_log is the
#'     standard deviation on log scale }
#'    \item{\strong{"poisson": } \code{priors = list(a=c(0.5,0.5),b=c(0.5,0.5))},
#'     where a and b are shape parameters of the gamma distribution.}
#'    \item{\strong{"negbin": } \code{priors = list(r=c(10,10), a=c(0.5,0.5),b=c(0.5,0.5))},
#'     where r is the r parameter of the negative binomial distribution (representing the number of successes required)
#'      and where a and b are shape parameters of the beta distribution. Note that the r value is not updated.
#'       The conjugate beta prior is only valid when r is fixed and known, which is a limitation for this method.}
#' }
#' 
#' See examples for plots of these prior distributions.
#' 
#' 
#' @examples 
#' makeMvLn<-function(bins=500,mu_log,sigma_log){
#'    setNames(data.frame(matrix(hist(rlnorm(2000,mu_log, sigma_log),
#'       breaks=seq(1,bins,5), plot=FALSE)$counts, nrow=1)),
#'     paste0("b",seq(1,bins,5))[-1] ) }
#' set.seed(123) 
#' mv_ln<-rbind(do.call(rbind,
#'            lapply(1:30, function(i){makeMvLn(mu_log=log(130),
#'              sigma_log=log(1.3) )})),
#'           do.call(rbind,
#'            lapply(1:30, function(i){makeMvLn(mu_log=log(100),
#'              sigma_log=log(1.2) )})))
#'              
#' # lognormal mv
#' ln_mv_ex <- conjugate(s1 = mv_ln[1:30,], s2= mv_ln[31:60,], method = "lognormal",
#'                    priors = list( mu_log=c(log(10),log(10)),n=c(1,1),sigma_log=c(log(3),log(3)) ),
#'                    plot=TRUE, rope_range = c(-40,40), rope_ci = 0.89,
#'                    cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' 
#' # lognormal sv
#' ln_sv_ex <- conjugate(s1=rlnorm(100, log(130), log(1.3)), s2= rlnorm(100, log(100), log(1.6)),
#'  method = "lognormal",
#'                    priors = list( mu_log=c(log(10),log(10)),n=c(1,1),sigma_log=c(log(3),log(3)) ),
#'                    plot=TRUE, rope_range = NULL, rope_ci = 0.89, 
#'                    cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' 
#' # Z test mv example
#'
#' makeMvGauss<-function(bins=180,mu,sigma){
#'    setNames(data.frame(matrix(hist(rnorm(2000,mu, sigma),
#'    breaks=seq(1,bins,1), plot=FALSE)$counts, nrow=1)),
#'    paste0("b",1:(bins-1) ))
#'    }
#' mv_gauss<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=50, sigma=10 )})),
#'                 do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=60, sigma=12 )})))
#'                 
#' gauss_mv_ex <- conjugate(s1=mv_gauss[1:30,], s2= mv_gauss[31:60,], method = "gaussian",
#'                   priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                   plot=TRUE, rope_range = c(-25, 25), rope_ci = 0.89,
#'                   cred.int.level = 0.89, hypothesis="equal", support=NULL)
#'                   
#' # Z test sv example
#' set.seed(123)
#' gauss_sv_ex_bad <- conjugate(s1=rnorm(15, 50,10), s2= rnorm(15, 60,12), method = "gaussian",
#'                   priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                   plot=TRUE, rope_range = c(-10, 10), rope_ci = 0.89, 
#'                   cred.int.level = 0.89, hypothesis="equal", support=NULL)
#'                   
#' # Here the plot clearly shows we have a problem with the default support, so we specify one
#' # naturally the longer the support vector the more time this takes, but supports below 100k length 
#' # tend to be reasonably fast.
#' 
#' gauss_sv_ex <- conjugate(s1=rnorm(15, 50,10), s2= rnorm(15, 60,12), method = "gaussian",
#'                   priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                   plot=TRUE, rope_range = c(-10, 10), rope_ci = 0.89, 
#'                   cred.int.level = 0.89, hypothesis="equal", support=seq(-20,120, 0.01))
#' # Note that the ROPE probability is somewhat unstable here since the distribution of differences
#' # is much wider than the ROPE interval.
#' 
#' # T test mv example 
#' 
#' makeMvGauss<-function(bins=180,mu,sigma){
#'    setNames(data.frame(matrix(hist(rnorm(2000,mu, sigma),
#'     breaks=seq(1,bins,1), plot=FALSE)$counts, nrow=1)),
#'    paste0("b",1:(bins-1) ))
#'    }
#' mv_gauss<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=50, sigma=10 )})),
#'                 do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=60, sigma=12 )})))
#'
#' gaussianMeans_mv_ex <- conjugate(s1=mv_gauss[1:30,], s2= mv_gauss[31:60,], method="t",
#'                         priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                         plot=TRUE, rope_range = c(-5,5), rope_ci = 0.89, 
#'                         cred.int.level = 0.89, hypothesis="equal", support=NULL)
#'                         
#'                         
#' # T test sv example
#' 
#' gaussianMeans_sv_ex <- conjugate(s1=rnorm(10, 50,10), s2= rnorm(10, 60,12), method="t",
#'                         priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                         plot=TRUE, rope_range = c(-5,8), rope_ci = 0.89, 
#'                         cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' 
#' 
#' 
#' # beta mv example
#' 
#' makeMvBeta<-function(n=100,a,b){
#'   setNames(data.frame(matrix(hist(rbeta(2000,a,b),
#'   breaks=seq(0,1,length.out=n), plot=FALSE)$counts, nrow=1)),
#'   paste0("b0.",1:(n-1)))
#' }
#' 
#' mv_beta<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvBeta(n=100, a=5, b=8 )})),
#'                do.call(rbind, lapply(1:30, function(i){makeMvBeta(n=100, a=10, b=3 )})))
#' 
#' beta_mv_ex <- conjugate(s1 = mv_beta[1:30,], s2= mv_beta[31:60,], method="beta",
#'               priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=TRUE, rope_range = c(-0.1, 0.1), rope_ci = 0.89, 
#'               cred.int.level = 0.89, hypothesis="equal")
#' 
#' # beta sv example
#' 
#' beta_sv_ex <- conjugate(s1 = rbeta(20, 5, 5), s2= rbeta(20, 8, 5), method="beta",
#'               priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=TRUE, rope_range = c(-0.1, 0.1), rope_ci = 0.89, 
#'               cred.int.level = 0.89, hypothesis="equal")
#'
#' # poisson sv example
#' 
#' poisson_sv_ex <- conjugate(s1 = rpois(20, 10), s2= rpois(20, 8), method="poisson",
#'               priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=TRUE, rope_range = c(-1, 1), rope_ci = 0.89, 
#'               cred.int.level = 0.89, hypothesis="equal")
#' 
#' # negative binomial sv example
#' # knowing r (required number of successes) is an important caveat for this method.
#' # in the current implementation we suggest using the poisson method for data such as leaf counts
#' 
#' negbin_sv_ex <- conjugate(s1 = rnbinom(20, 10, 0.5), s2= rnbinom(20, 10, 0.25), method="negbin",
#'               priors = list(r = c(10, 10), a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=TRUE, rope_range = c(-1, 1), rope_ci = 0.89, 
#'               cred.int.level = 0.89, hypothesis="equal")
#' 
#' # Dirichlet mv example
#' 
#' makeMvLn<-function(bins=500,mu_log,sigma_log){
#' setNames(data.frame(matrix(hist(rlnorm(2000,mu_log, sigma_log),
#'                                 breaks=seq(1,bins,5), plot=FALSE)$counts, nrow=1)),
#'                                          paste0("b",seq(1,bins,5))[-1] ) }
#' set.seed(123) 
#' mv_ln<-rbind(do.call(rbind,
#'                       lapply(1:30, function(i){makeMvLn(mu_log=log(130),
#'                                           sigma_log=log(1.3) )})),
#'             do.call(rbind, lapply(1:30, function(i){makeMvLn(mu_log=log(100),
#'                                           sigma_log=log(1.2) )})))
#' s1 <- mv_ln[1:30, ]
#' s2 <- mv_ln[31:60, ]
#' diri_ex_1 <- conjugate(s1, s2, method = "dirichlet", priors=NULL, plot=TRUE,
#'       rope_range = c(-0.025, 0.025), rope_ci = 0.89, 
#'       cred.int.level = 0.89, hypothesis="equal")
#'       
#'       
#'       
#' 
#' 
#' # Example usage with plantCV data
#' ## Not run:
#' library(data.table)
#' wide<-read.pcv(
#'  paste0("https://media.githubusercontent.com/media/joshqsumner/",
#'  "pcvrTestData/main/smallPhenotyperRun.csv"),
#'  mode="wide", singleValueOnly =TRUE, multiValPattern = "hist", reader="fread")
#' wide$genotype = substr(wide$barcode, 3,5)
#' wide$genotype = ifelse(wide$genotype == "002", "B73",
#'                        ifelse(wide$genotype == "003", "W605S",
#'                               ifelse(wide$genotype == "004", "MM", "Mo17")))
#' wide$fertilizer = substr(wide$barcode, 8, 8)
#' wide$fertilizer = ifelse(wide$fertilizer == "A", "100",
#'                          ifelse(wide$fertilizer == "B", "50", "0"))
#' wide<-bw.time(wide,timeCol="timestamp", group="barcode")
#' 
#' mo17_sample <- wide[wide$genotype=="Mo17" &
#'                     wide$DAS > 18 &
#'                     wide$fertilizer == 100,
#'                     grepl("hue_freq", colnames(wide))]
#' B73_sample <- wide[wide$genotype=="B73" &
#'                    wide$DAS > 18 &
#'                    wide$fertilizer == 100,
#'                    grepl("hue_freq", colnames(wide))]
#' 
#' hue_res_t <- conjugate(s1 = mo17_sample, s2= B73_sample, method="t",
#'               priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'               plot=TRUE, rope_range = c(-10,10), rope_ci = 0.89, 
#'               cred.int.level = 0.89, hypothesis="equal")
#'               
#' hue_res_ln <- conjugate(s1 = mo17_sample, s2= B73_sample, method="lognormal",
#'               plot=TRUE, rope_range = c(-10,10), rope_ci = 0.89, 
#'               cred.int.level = 0.89, hypothesis="equal")
#' # Picking the right distribution makes a difference,
#' # checking plots of your data (see ?pcv.joyplot for mv traits)
#' # will be useful.
#' 
#' pixels_per_cmsq <- 42.5^2   # pixel per cm^2
#' wide$area_cm2<-wide$area.pixels / pixels_per_cmsq
#' 
#' mo17_area <- wide[wide$genotype=="Mo17" & wide$DAS > 18 & wide$fertilizer == 50, "area_cm2"]
#' B73_area <- wide[wide$genotype=="B73" & wide$DAS > 18 & wide$fertilizer == 50, "area_cm2"]
#' 
#' area_res_t <- conjugate(s1 = mo17_area, s2= B73_area, method="t",
#'               priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'               plot=TRUE, rope_range = c(-5,5), rope_ci = 0.89, 
#'               cred.int.level = 0.89, hypothesis="equal")
#' 
#' ## End(Not run)
#' 
#' # Plots of prior distributions
#' set.seed(123)
#' plot(seq(0,1,0.0001), dbeta(seq(0,1,0.0001), 0.5, 0.5), ylab="Density",
#'  xlab="Support", main="Beta", type="l")
#' 
#' plot(seq(0,10,0.01), dgamma(seq(0,10,0.01), 0.5, 0.5), ylab="Density", xlab="Support",
#'  main="Poisson (gamma prior on Lambda parameter)", type="l")
#'
#' plot(seq(-20,20.001), dnorm(seq(-20,20.001), 0, sqrt(20)), ylab="Density", 
#' xlab="Support", main="T and Gaussian", type="l")
#' 
#' plot(seq(0,70,0.001), dlnorm(seq(0,70,0.001), log(10), log(3)), ylab="Density", 
#' xlab="Support", main="Lognormal", type="l")
#' 
#' @return 
#' 
#' A list with named elements:
#' \itemize{
#'    \item{\strong{summary}: A data frame containing HDI/HDE values for each sample and 
#'    the ROPE as well as posterior probability of the hypothesis.}
#'    \item{\strong{posterior}: A list of updated parameters in the same format as the prior
#'     for the given method. If desired this does allow for Bayesian updating.}
#'    \item{\strong{plot_df}: A data frame of probabilities along the support for each sample.
#'     This is used for making the ggplot.}
#'    \item{\strong{rope_df}: A data frame of draws from the ROPE posterior.}
#'    \item{\strong{plot}: A ggplot showing the distribution of samples and optionally the 
#'    distribution of differences/ROPE}
#' }
#' 
#' @keywords bayesian, conjugate, ROPE
#' @export

conjugate<-function(s1 = NULL, s2= NULL, method = c("t", "gaussian", "beta", "lognormal", "poisson", "negbin", "dirichlet"),
                     priors=NULL, plot=FALSE, rope_range = NULL,
                    rope_ci = 0.89, cred.int.level = 0.89, hypothesis="equal",
                    support=NULL){
  #* `Check sample class`
  if(is.matrix(s1)|is.data.frame(s1)){ vec=FALSE
  } else if(is.vector(s1)) { vec=TRUE 
  } else{ stop("s1 must be a vector, data.frame, or matrix.") }
  
  matched_arg<-match.arg(method, choices = c("t", "gaussian", "beta", "lognormal", "poisson", "negbin", "dirichlet"))
  #* `Pick method`
  
  if(matched_arg == "t" & vec){
    res<-.conj_gaussian_means_sv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis, support)
    
  } else if(matched_arg == "t" & !vec){
    res<-.conj_gaussian_means_mv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis, support)
    
  } else if(matched_arg == "gaussian" & vec){
    res<-.conj_gaussian_sv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis, support)
    
  } else if(matched_arg == "gaussian" & !vec){
    res<-.conj_gaussian_mv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis, support)
    
  } else if(matched_arg == "beta" & vec){
    res<-.conj_beta_sv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis)
    
  } else if(matched_arg == "beta" & !vec){
    res<-.conj_beta_mv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis)
    
  } else if(matched_arg == "lognormal" & vec){
    res<-.conj_lognormal_sv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis, support)
    
  } else if(matched_arg == "lognormal" & !vec){
    res<-.conj_lognormal_mv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis, support)
    
  } else if(matched_arg == "poisson" & vec){
    res<-.conj_poisson_sv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis, support)
    
  } else if(matched_arg == "negbin" & vec){
    res<-.conj_negbin_sv(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis, support)
    
  } else if(matched_arg == "dirichlet" & !vec){
    res<-.conj_diri_mv_1(s1, s2, priors, plot, rope_range, rope_ci, cred.int.level, hypothesis)
  } else {
    stop(paste0("Selected distribution and data type are not supported, review methods and check if you are using SV traits",
                " with an MV only distribution (dirichlet) or MV traits with an SV only distribution (poisson, negbin)"))
  }

  #* `Make plot`
  if(matched_arg != "dirichlet" & plot){
    res <- .conj_general_plot(res, s2, rope_range, support, rope_ci)
  } else if (matched_arg == "dirichlet" & plot){
    res <- .conj_diri_plot(res, s2, rope_range, rope_ci)
  }
  
  return(res)
}

#' ***********************************************************************************************
#' *************** `lognormal distribution (multi value)` ****************************************
#' ***********************************************************************************************
#' log(130), log(1.3)


#' @description
#' Internal function for Bayesian comparison of log-normal data represented by single value traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number representing the "bin".
#' @param s2 An optional second sample of the same form as s2.
#' 
#' @examples 
#' makeMvLn<-function(bins=500,mu_log,sigma_log){
#'    setNames(data.frame(matrix(hist(rlnorm(2000,mu_log, sigma_log), breaks=seq(1,bins,5), plot=FALSE)$counts, nrow=1)), paste0("b",seq(1,bins,5))[-1] ) }
#'    set.seed(123)
#' mv_ln<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvLn(mu_log=log(130), sigma_log=log(1.3) )})),
#'              do.call(rbind, lapply(1:30, function(i){makeMvLn(mu_log=log(100), sigma_log=log(1.2) )})))
#'              
#' .conj_lognormal_mv(s1 = mv_ln[1:30,], s2= mv_ln[31:60,],
#'                    priors = list( mu_log=c(log(10),log(10)),n=c(1,1),sigma_log=c(log(3),log(3)) ),
#'                    plot=FALSE, rope_range = NULL, rope_ci = 0.89,
#'                    cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' @keywords internal
#' @noRd

.conj_lognormal_mv<-function(s1 = NULL, s2= NULL,
                             priors=NULL,
                             plot=FALSE, rope_range = NULL, rope_ci = 0.89,
                             cred.int.level = 0.89, hypothesis="equal", support=NULL){
  out <- list()
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list( mu_log=c(log(10),log(10)),n=c(1,1),sigma_log=c(log(3),log(3)) )
  }
  #* `Standardize sample 1 class and names`
  if(is.null(colnames(s1))){
    bins<-(1:ncol(s1))
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
    support<-seq(min(histCols_bin), max(histCols_bin), length.out=10000)
  }
  
  #* `Turn s1 matrix into a vector`
  X1<-rep(histCols_bin[bins_order], as.numeric(round(colSums(s1))) )
  
  #* `Get mean and variance from s1`
    xbar_1 <- mean(X1)
    s2_1 <- var(X1)/(nrow(s1)-1)
    n1 <-nrow(s1)
    #* `Add prior distribution`
    # n1_n<-n1+priors$n[1]
    # xbar_1_n<-(xbar_1*n1+ priors$mu[1]*priors$n[1])/n1_n
    # s2_1_n <- ( (s2_1*(n1-1)) + (priors$s2[1]*(priors$n[1]-1)) ) / (n1_n-2)
    #* `lognormal method of moments`
    #* Might change prior to use lognormal params and add them in here.
    
    # mu_ls = log(xbar_1_n / sqrt((s2_1_n/xbar_1_n^2)+1))
    # sigma_ls = sqrt(log((s2_1_n/xbar_1_n^2)+1))
    
    mu_ls = log(xbar_1 / sqrt((s2_1/xbar_1^2)+1))
    sigma_ls = sqrt(log((s2_1/xbar_1^2)+1))
    
    n1_n<-n1+priors$n[1]
    mu_ls_n = ((mu_ls * n1) + (priors$mu_log[1] * priors$n[1])) / n1_n
    sigma_ls_n = ((sigma_ls * n1) + (priors$sigma_log[1] * priors$n[1])) / n1_n
    
    #* `posterior`
    dens1 <- dlnorm(support, mu_ls_n, sigma_ls_n)
    pdf1 <- dens1/sum(dens1)
    hde1 <- exp(mu_ls_n)
    hdi1 <- qlnorm(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))), mu_ls_n, sigma_ls_n)
    #* `Store summary`
    out$summary<-data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
    out$posterior$mu_log <- mu_ls_n
    out$posterior$n <- n1_n
    out$posterior$sigma_log <- sigma_ls_n
    #* `save s1 data for plotting`
    if(plot){
      out$plot_df <- data.frame("range"=support,
                                "prob"=pdf1,
                                "sample"=rep("Sample 1",length(support)))
    }
  if(!is.null(s2)){
    #* `Standardize sample 1 class and names`
    if(is.null(colnames(s2))){
      bins<-(1:ncol(s2))
      colnames(s2)<-paste0("b", bins )
      warning(paste0("Assuming unnamed columns represent bins from ", min(bins), " to ", max(bins)))
    }
    if(is.matrix(s2)){ s2<-as.data.frame(s2) }
    
    #* `Reorder columns if they are not in the numeric order`
    histCols<-colnames(s2)
    histCols_bin <- as.numeric(sub("[a-zA-Z_.]+","",colnames(s1)))
    bins_order<-sort(histCols_bin, index.return=TRUE)$ix
    s2<-s2[,bins_order]
    
    #* `Turn s2 matrix into a vector`
    X2<-rep(histCols_bin[bins_order], as.numeric(round(colSums(s2))) )
      
      xbar_2 <- mean(X2)
      s2_2 <- var(X2)/(nrow(s2)-1)
      n2 <-nrow(s2)
      # n2_n<-n2+priors$n[2]
      # xbar_2_n<-(xbar_2*n2+ priors$mu[2]*priors$n[2])/n2_n
      # s2_2_n <- ( (s2_2*(n2-1)) + (priors$s2[2]*(priors$n[2]-1)) ) / (n2_n-2)
      # 
      # mu_ls_2 = log(xbar_2_n / sqrt((s2_2_n/xbar_2_n^2)+1))
      # sigma_ls_2 = sqrt(log((s2_2_n/xbar_2_n^2)+1))
      
      mu_ls_2 = log(xbar_2 / sqrt((s2_2/xbar_2^2)+1))
      sigma_ls_2 = sqrt(log((s2_2/xbar_2^2)+1))
      
      n2_n<-n2+priors$n[2]
      mu_ls_2_n = ((mu_ls_2 * n2) + (priors$mu_log[2] * priors$n[2])) / n2_n
      sigma_ls_2_n = ((sigma_ls_2 * n2) + (priors$sigma_log[2] * priors$n[2])) / n2_n
      
      
      dens2 <- dlnorm(support, mu_ls_2_n, sigma_ls_2_n)
      pdf2 <- dens2/sum(dens2)
      hde2 <- exp(mu_ls_2_n)
      hdi2 <- qlnorm(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))), mu_ls_2_n, sigma_ls_2_n)
      
      out$summary<-data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2],
                              HDE_2=hde2, HDI_2_low = hdi2[1], HDI_2_high = hdi2[2])
      out$posterior$mu_log <- c(mu_ls_n, mu_ls_2_n)
      out$posterior$n <- c(n1_n, n2_n)
      out$posterior$sigma_log <- c(sigma_ls_n, sigma_ls_2_n)
      
      #* `Keep data from s2 for plotting`
      if(plot){
        out$plot_df <- rbind(out$plot_df, data.frame("range"=support,
                                                     "prob"=pdf2,
                                                     "sample"=rep("Sample 2",length(support)))
        )
      }
      #* `Test hypothesis on posterior distributions`
      if(!all(c(as.matrix(s1), as.matrix(s2))==0)){
        post.prob.out <-  .post.prob.from.pdfs(pdf1, pdf2, hypothesis)
        post.prob <-post.prob.out$post.prob
        dirSymbol <-post.prob.out$direction
      }else{
        post.prob <- 1
        dirSymbol = "="
      }
      out$summary$hyp = hypothesis
      out$summary$post.prob = post.prob
      out$dirSymbol = dirSymbol
  }
  #* `ROPE Comparison`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = rlnorm(10000,mu_ls_n, sigma_ls_n)
      if(!is.null(s2)){
        post2 = rlnorm(10000,mu_ls_2_n, sigma_ls_2_n)
        posterior = post1 - post2
      } else {
        posterior=post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  return(out)
}




#' ***********************************************************************************************
#' *************** `lognormal distribution (single value)` ****************************************
#' ***********************************************************************************************

#' @description
#' Internal function for Bayesian comparison of log-normal data represented by single value traits.
#' Lognormal method of moments
#' 
#' \bar{x} \sim e^{\mu + \sigma^2 / 2}
#' \s^2 \sim (e^{\alpha^2} -1) \cdot e^{2 \cdot \mu + \sigma^2}
#' 
#' Calculate Sigma:
#' 
#' \bar{x}^2 \sim e^{2 \cdot \mu + \sigma^2}
#' s^2 \sim (e^\alpha^2 -1) \cdot \bar{x}^2
#' (s^2 / \bar{x}^2) +1 \sim e^\alpha^2
#' \alpha^2 \sim ln((s^2 / \bar{x}^2) +1 )
#' \alpha \sim \sqrt{ln((s^2 / \bar{x}^2) +1 )}
#' 
#' Calculate Mu:
#' 
#' \bar{x}^2 \sim e^{2 \cdot \mu + \sigma^2}
#' 2 \cdot ln(\bar{x}) \sim 2 \cdot \mu + \sigma^2
#' \mu \sim ln(\bar{x} / \sqrt{(s^2 / \bar{x}^2) +1 })
#' 
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' .conj_lognormal_sv(s1=rlnorm(100, log(130), log(1.3)), s2= rlnorm(100, log(100), log(1.6)),
#'                   priors = list( mu_log=c(log(10),log(10)),n=c(1,1),sigma_log=c(log(3),log(3)) ),
#'                   plot=FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89, 
#'                   cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' @keywords internal
#' @noRd

.conj_lognormal_sv<-function(s1 = NULL, s2= NULL, 
                             priors=NULL,
                             plot=FALSE, rope_range = NULL, rope_ci = 0.89,
                             cred.int.level = 0.89, hypothesis="equal", support=NULL){
  out <- list()
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list( mu_log=c(log(10),log(10)),n=c(1,1),sigma_log=c(log(3),log(3)) )
  }
  #* `Define support if it is missing`
  if(is.null(support)){
    support<-seq(floor(min(c(s1,s2))*0.75), ceiling(max(c(s1,s2))/0.75), length.out=10000)
  }
  #* `Get mean and variance from s1`
  if(length(s1)>1){
    xbar_1 <- mean(s1)
    s2_1 <- var(s1)/(length(s1)-1)
    n1 <-length(s1)
    #* `Add prior distribution`
    # n1_n<-n1+priors$n[1]
    # xbar_1_n<-(xbar_1*n1+ priors$mu[1]*priors$n[1])/n1_n
    # s2_1_n <- ( (s2_1*(n1-1)) + (priors$s2[1]*(priors$n[1]-1)) ) / (n1_n-2)
    #* `lognormal method of moments`
    #* Might change prior to use lognormal params and add them in here.
    
    # mu_ls = log(xbar_1_n / sqrt((s2_1_n/xbar_1_n^2)+1))
    # sigma_ls = sqrt(log((s2_1_n/xbar_1_n^2)+1))
    
    mu_ls = log(xbar_1 / sqrt((s2_1/xbar_1^2)+1))
    sigma_ls = sqrt(log((s2_1/xbar_1^2)+1))
    
    n1_n<-n1+priors$n[1]
    mu_ls_n = ((mu_ls * n1) + (priors$mu_log[1] * priors$n[1])) / n1_n
    sigma_ls_n = ((sigma_ls * n1) + (priors$sigma_log[1] * priors$n[1])) / n1_n
    
    #* `posterior`
    dens1 <- dlnorm(support, mu_ls_n, sigma_ls_n)
    pdf1 <- dens1/sum(dens1)
    hde1 <- exp(mu_ls_n)
    hdi1 <- qlnorm(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))), mu_ls_n, sigma_ls_n)
    #* `Store summary`
    out$summary<-data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2])
    out$posterior$mu_log <- mu_ls_n
    out$posterior$n <- n1_n
    out$posterior$sigma_log <- sigma_ls_n
    #* `save s1 data for plotting`
    if(plot){
      out$plot_df <- data.frame("range"=support,
                                "prob"=pdf1,
                                "sample"=rep("Sample 1",length(support)))
    } 
  } else {
    stop("s1 must be a numeric of length 2 or greater")
  }
  if(!is.null(s2)){
    if(length(s2) > 1){
     
      xbar_2 <- mean(s2)
      s2_2 <- var(s2)/(length(s2)-1)
      n2 <-length(s2)
      # n2_n<-n2+priors$n[2]
      # xbar_2_n<-(xbar_2*n2+ priors$mu[2]*priors$n[2])/n2_n
      # s2_2_n <- ( (s2_2*(n2-1)) + (priors$s2[2]*(priors$n[2]-1)) ) / (n2_n-2)
      # 
      # mu_ls_2 = log(xbar_2_n / sqrt((s2_2_n/xbar_2_n^2)+1))
      # sigma_ls_2 = sqrt(log((s2_2_n/xbar_2_n^2)+1))
      
      mu_ls_2 = log(xbar_2 / sqrt((s2_2/xbar_2^2)+1))
      sigma_ls_2 = sqrt(log((s2_2/xbar_2^2)+1))
      
      n2_n<-n2+priors$n[2]
      mu_ls_2_n = ((mu_ls_2 * n2) + (priors$mu_log[2] * priors$n[2])) / n2_n
      sigma_ls_2_n = ((sigma_ls_2 * n2) + (priors$sigma_log[2] * priors$n[2])) / n2_n
      
      
      dens2 <- dlnorm(support, mu_ls_2_n, sigma_ls_2_n)
      pdf2 <- dens2/sum(dens2)
      hde2 <- exp(mu_ls_2_n)
      hdi2 <- qlnorm(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))), mu_ls_2_n, sigma_ls_2_n)
      
      out$summary<-data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2],
                              HDE_2=hde2, HDI_2_low = hdi2[1], HDI_2_high = hdi2[2])
      out$posterior$mu_log <- c(mu_ls_n, mu_ls_2_n)
      out$posterior$n <- c(n1_n, n2_n)
      out$posterior$sigma_log <- c(sigma_ls_n, sigma_ls_2_n)
      
      #* `Keep data from s2 for plotting`
      if(plot){
        out$plot_df <- rbind(out$plot_df, data.frame("range"=support,
                                                     "prob"=pdf2,
                                                     "sample"=rep("Sample 2",length(support)))
        )
      }
      #* `Test hypothesis on posterior distributions`
      if(!all(c(s1, s2)==0)){
        post.prob.out <-  .post.prob.from.pdfs(pdf1, pdf2, hypothesis)
        post.prob <-post.prob.out$post.prob
        dirSymbol <-post.prob.out$direction
      }else{
        post.prob <- 1
        dirSymbol = "="
      }
      out$summary$hyp = hypothesis
      out$summary$post.prob = post.prob
      out$dirSymbol = dirSymbol
    } else{
      stop("s2 must be a numeric of length 2 or greater if provided")
    }
  }
  #* `ROPE Comparison`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = rlnorm(10000,mu_ls, sigma_ls)
      if(!is.null(s2)){
        post2 = rlnorm(10000,mu_ls_2, sigma_ls_2)
        posterior = post1 - post2
      } else {
        posterior=post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  return(out)
}



#' ***********************************************************************************************
#' *************** `Gaussian Dist (single value, Z test)` ****************************************
#' ***********************************************************************************************
#' same as T, just don't use standard error
#' @description
#' Internal function for Bayesian comparisosns of gaussian data represented by single value traits.
#' This version uses the entire posterior distribution instead of the sampling distribution of the mean.
#' In frequentist terms this is analogous to a Z test as opposed to a T test. Generally the T test is desired, 
#' but this is provided for completeness.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' .conj_gaussian_sv(s1=rnorm(100, 50,10), s2= rnorm(100, 60,12),
#'                   priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                   plot=FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'                   cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' @keywords internal
#' @noRd

.conj_gaussian_sv<-function(s1 = NULL, s2= NULL, priors = NULL,
                                  plot=FALSE, rope_range = NULL, rope_ci = 0.89,
                                  cred.int.level = 0.89, hypothesis="equal", support=NULL){
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
    s2_1_n = ((n1-1)*s2_1 + v1*priors$s2[1] + priors$n[1]*n1*(priors$mu[1] - m1)^2 / n1_n)/v1_n # pooled variance
    sigma_1<-sqrt(s2_1_n) # standard error of the mean
    
    post_mean_d1 <- extraDistr::dlst(support, v1_n, m1_n, sigma_1)
    post_mean_pdf1 <- post_mean_d1/sum(post_mean_d1)
    hde1_mean <- m1_n
    hdi1_mean <- m1_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v1_n)*sigma_1
    
    out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
    out$posterior$mu <- m1_n
    out$posterior$n <- n1_n
    out$posterior$s2 <- sigma_1
    #* `Save data for plotting`
    if(plot){
      out$plot_df <- data.frame("range"=support,
                                "prob"=post_mean_pdf1,
                                "sample"=rep("Sample 1",length(support))
      )
    } 
  }else{
    stop("s1 must be a numeric of length 2 or greater")
  }
  #* `Get Mean, Variance, SE, and DF from s2`
  if(!is.null(s2)){
    if(length(s2) > 1){
      n2 <- length(s2) # n samples
      m2 = mean(s2) # xbar
      s2_2 = var(s2) # s^2
      
      v2 = priors$n[2] - 1 # prior DF
      n2_n = priors$n[2] + n2 # total N including prior
      m2_n = (n2*m2 + priors$n[2]*priors$mu[2])/n2_n # weighted mean of prior and data
      v2_n = v2 + n2 # degrees of freedom including data
      s2_2_n = ((n2-1)*s2_2 + v2*priors$s2[2] + priors$n[2]*n2*(priors$mu[2] - m2)^2/n2_n)/v2_n # pooled variance
      sigma_2<-sqrt(s2_2_n) # standard error of the mean
      
      post_mean_d2 <- extraDistr::dlst(support,v2_n,m2_n, sigma_2)
      post_mean_pdf2 <- post_mean_d2/sum(post_mean_d2)
      hde2_mean <- m2_n
      hdi2_mean <- m2_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v2_n)*sigma_2
      
      out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2],
                              HDE_2=hde2_mean, HDI_2_low = hdi2_mean[1], HDI_2_high = hdi2_mean[2])
      out$posterior$mu <- c(m1_n, m2_n)
      out$posterior$n <- c(n1_n, n2_n)
      out$posterior$s2 <- c(sigma_1, sigma_2)
      #* `Keep data from s2 for plotting`
      if(plot){
        out$plot_df <- rbind(out$plot_df, data.frame("range"=support,
                                                     "prob"=post_mean_pdf2,
                                                     "sample"=rep("Sample 2",length(support)))
        )
      } 
      #* `Test hypothesis on posterior distributions`
      if(!all(c(s1, s2)==0)){
        post.prob.out <-  .post.prob.from.pdfs(post_mean_pdf1, post_mean_pdf2, hypothesis)
        post.prob <-post.prob.out$post.prob
        dirSymbol <-post.prob.out$direction
      }else{
        post.prob <- 1
        dirSymbol = "="
      }
      out$summary$hyp = hypothesis
      out$summary$post.prob = post.prob
      out$dirSymbol = dirSymbol
    }else{
      stop("If provided then s2 must be a numeric of length 2 or greater")
    }
  }
  #* `ROPE Comparison`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = extraDistr::rlst(10000,v1_n, m1_n, sigma_1)
      if(!is.null(s2)){
        post2 = extraDistr::rlst(10000,v2_n, m2_n, sigma_2)
        posterior = post1 - post2
      } else {
        posterior=post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  return(out)
}



#' ***********************************************************************************************
#' *************** `Gaussian Dist (multi value, Z test)` ****************************************
#' ***********************************************************************************************

#' @description
#' Internal function for Bayesian comparisosns of gaussian data represented by multi value traits.
#' This version uses the entire posterior distribution instead of the sampling distribution of the mean.
#' In frequentist terms this is analogous to a Z test as opposed to a T test. Generally the T test is desired, 
#' but this is provided for completeness.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' 
#' makeMvGauss<-function(bins=180,mu,sigma){
#'    setNames(data.frame(matrix(hist(rnorm(2000,mu, sigma),
#'    breaks=seq(1,bins,1), plot=FALSE)$counts, nrow=1)),paste0("b",1:(bins-1) ))
#'    }
#' mv_gauss<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=50, sigma=10 )})),
#'                 do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=60, sigma=12 )})))
#'                 
#' .conj_gaussian_mv(s1=mv_gauss[1:30,], s2= mv_gauss[31:60,],
#'                   priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                   plot=FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'                   cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' @keywords internal
#' @noRd

.conj_gaussian_mv<-function(s1 = NULL, s2= NULL, priors = NULL,
                                  plot=FALSE, rope_range = NULL, rope_ci = 0.89,
                                  cred.int.level = 0.89, hypothesis="equal", support=NULL){
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
    support<-seq(min(bins_order), max(bins_order), length.out=10000)
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
  sigma_1<-sqrt(s2_1_n) # standard error of the mean
  
  post_mean_d1 <- extraDistr::dlst(support,v1_n,m1_n, sigma_1)
  post_mean_pdf1 <- post_mean_d1/sum(post_mean_d1)
  hde1_mean <- m1_n
  hdi1_mean <- m1_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v1_n)*sigma_1
  
  out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- m1_n
  out$posterior$n <- n1_n
  out$posterior$s2 <- sigma_1
  #* `Save data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=support,
                              "prob"=post_mean_pdf1,
                              "sample"=rep("Sample 1",length(support))
    )
  }
  
  if(!is.null(s2)){
    #* `Standardize sample 2 class and names`
    if(is.null(colnames(s2))){
      bins<-(1:ncol(s2))/100
      colnames(s2)<-paste0("b", bins )
      warning(paste0("Assuming unnamed columns represent bins from ", min(bins), " to ", max(bins)))
    }
    if(is.matrix(s2)){ s2<-as.data.frame(s2) }
    
    #* `Reorder columns if they are not in the numeric order`
    histCols<-colnames(s2)
    histCols_bin <- as.numeric(sub("[a-zA-Z_.]+","",colnames(s2)))
    bins_order<-sort(histCols_bin, index.return=TRUE)$ix
    s2<-s2[,bins_order]
    
    #* `Turn s2 matrix into a vector`
    X2<-rep(histCols_bin[bins_order], as.numeric(round(colSums(s2))) )
    
    #* `Get Mean, Variance, SE, and DF from s2`
    n2 <- nrow(s2) # n samples
    m2 = mean(X2) # xbar
    s2_2 = var(X2) # s^2
    
    v2 = priors$n[2] - 1 # prior DF
    n2_n = priors$n[2] + n2 # total N including prior
    m2_n = (n2*m2 + priors$n[2]*priors$mu[2])/n2_n # weighted mean of prior and data
    v2_n = v2 + n2 # degrees of freedom including data
    s2_2_n = ((n2-1)*s2_2 + v2*priors$s2[2] + priors$n[2]*n2*(priors$mu[2] - m2)^2/n2_n)/v2_n # pooled variance
    sigma_2<-sqrt(s2_2_n) # standard error of the mean
    
    post_mean_d2 <- extraDistr::dlst(support,v2_n,m2_n, sigma_2)
    post_mean_pdf2 <- post_mean_d2/sum(post_mean_d2)
    hde2_mean <- m2_n
    hdi2_mean <- m2_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v2_n)*sigma_2
    
    out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2],
                            HDE_2=hde2_mean, HDI_2_low = hdi2_mean[1], HDI_2_high = hdi2_mean[2])
    out$posterior$mu <- c(m1_n, m2_n)
    out$posterior$n <- c(n1_n, n2_n)
    out$posterior$s2 <- c(sigma_1, sigma_2)
    #* `Save data for plotting`
    if(plot){
      out$plot_df <- rbind(out$plot_df, data.frame("range"=support,
                                                   "prob"=post_mean_pdf2,
                                                   "sample"=rep("Sample 2",length(support)))
      )
    }
    if(!all(c(as.matrix(s1),as.matrix(s2))==0)){
      post.prob.out <-  .post.prob.from.pdfs(post_mean_pdf1, post_mean_pdf2, hypothesis)
      post.prob <-post.prob.out$post.prob
      dirSymbol <-post.prob.out$direction
    }else{
      post.prob <- 1
    }
    out$summary$hyp = hypothesis
    out$summary$post.prob = post.prob
    out$dirSymbol = dirSymbol
  }
  #* `ROPE Comparison`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = extraDistr::rlst(10000,v1_n, m1_n, sigma_1)
      if(!is.null(s2)){
        post2 = extraDistr::rlst(10000,v2_n, m2_n, sigma_2)
        posterior = post1 - post2
      } else {
        posterior = post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  
  return(out)
  
}



#' ***********************************************************************************************
#' *************** `Gaussian Means (single value)` *********************************************
#' ***********************************************************************************************

#' @description
#' Internal function for Bayesian T Tests of gaussian data represented by single value traits.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' .conj_gaussian_means_sv(s1=rnorm(100, 50,10), s2= rnorm(100, 60,12),
#'                         priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                         plot=FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'                         cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' @keywords internal
#' @noRd
.conj_gaussian_means_sv<-function(s1 = NULL, s2= NULL, priors = NULL,
                            plot=FALSE, rope_range = NULL, rope_ci = 0.89,
                            cred.int.level = 0.89, hypothesis="equal", support=NULL){
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
    
    post_mean_d1 <- extraDistr::dlst(support,v1_n,m1_n, se1)
    post_mean_pdf1 <- post_mean_d1/sum(post_mean_d1)
    hde1_mean <- m1_n
    hdi1_mean <- m1_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v1_n)*se1
    
    out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
    out$posterior$mu <- m1_n
    out$posterior$n <- n1_n
    out$posterior$s2 <- s2_1_n
    #* `Save data for plotting`
    if(plot){
      out$plot_df <- data.frame("range"=support,
                       "prob"=post_mean_pdf1,
                       "sample"=rep("Sample 1",length(support))
      )
    } 
  }else{
    stop("s1 must be a numeric of length 2 or greater")
  }
  #* `Get Mean, Variance, SE, and DF from s2`
  if(!is.null(s2)){
    if(length(s2) > 1){
      n2 <- length(s2) # n samples
      m2 = mean(s2) # xbar
      s2_2 = var(s2) # s^2
      
      v2 = priors$n[2] - 1 # prior DF
      n2_n = priors$n[2] + n1 # total N including prior
      m2_n = (n2*m2 + priors$n[2]*priors$mu[2])/n2_n # weighted mean of prior and data
      v2_n = v2 + n2 # degrees of freedom including data
      s2_2_n = ((n2-1)*s2_2 + v2*priors$s2[2] + priors$n[2]*n2*(priors$mu[2] - m2)^2/n2_n)/v2_n # pooled variance
      se2<-sqrt(s2_2_n/n2_n) # standard error of the mean
      
      post_mean_d2 <- extraDistr::dlst(support,v2_n,m2_n, se2)
      post_mean_pdf2 <- post_mean_d2/sum(post_mean_d2)
      hde2_mean <- m2_n
      hdi2_mean <- m2_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v2_n)*se2
      
      out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2],
                              HDE_2=hde2_mean, HDI_2_low = hdi2_mean[1], HDI_2_high = hdi2_mean[2])
      out$posterior$mu <- c(m1_n, m2_n)
      out$posterior$n <- c(n1_n, n2_n)
      out$posterior$s2 <- c(s2_1_n, s2_2_n)
      #* `Keep data from s2 for plotting`
      if(plot){
        out$plot_df <- rbind(out$plot_df, data.frame("range"=support,
                                  "prob"=post_mean_pdf2,
                                  "sample"=rep("Sample 2",length(support)))
        )
      } 
    #* `Test hypothesis on posterior distributions`
    if(!all(c(s1, s2)==0)){
      post.prob.out <-  .post.prob.from.pdfs(post_mean_pdf1, post_mean_pdf2, hypothesis)
      post.prob <-post.prob.out$post.prob
      dirSymbol <-post.prob.out$direction
    }else{
      post.prob <- 1
      dirSymbol = "="
    }
    out$summary$hyp = hypothesis
    out$summary$post.prob = post.prob
    out$dirSymbol = dirSymbol
    }else{
      stop("If provided then s2 must be a numeric of length 2 or greater")
    }
  }
  #* `ROPE Comparison`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = extraDistr::rlst(10000,v1_n, m1_n, se1)
      if(!is.null(s2)){
        post2 = extraDistr::rlst(10000,v2_n, m2_n, se2)
        posterior = post1 - post2
      } else {
        posterior=post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  return(out)
}





#' ***********************************************************************************************
#' *************** `Gaussian Means (multi value)` *********************************************
#' ***********************************************************************************************

#' @description
#' Internal function for calculating \alpha and \beta of a distribution represented by multi value traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number representing the "bin".
#' @param s2 An optional second sample of the same form as s2.
#' @examples
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
#' @keywords internal
#' @noRd

.conj_gaussian_means_mv<-function(s1 = NULL, s2= NULL, priors = NULL,
                            plot=FALSE, rope_range = NULL, rope_ci = 0.89, 
                            cred.int.level = 0.89, hypothesis="equal", support=NULL){
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
    support<-seq(min(bins_order), max(bins_order), length.out=10000)
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
  
  post_mean_d1 <- extraDistr::dlst(support,v1_n,m1_n, se1)
  post_mean_pdf1 <- post_mean_d1/sum(post_mean_d1)
  hde1_mean <- m1_n
  hdi1_mean <- m1_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v1_n)*se1
  
  out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2])
  out$posterior$mu <- m1_n
  out$posterior$n <- n1_n
  out$posterior$s2 <- s2_1_n
  #* `Save data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=support,
                              "prob"=post_mean_pdf1,
                              "sample"=rep("Sample 1",length(support))
    )
  }
  
  if(!is.null(s2)){
    #* `Standardize sample 2 class and names`
    if(is.null(colnames(s2))){
      bins<-(1:ncol(s2))/100
      colnames(s2)<-paste0("b", bins )
      warning(paste0("Assuming unnamed columns represent bins from ", min(bins), " to ", max(bins)))
    }
    if(is.matrix(s2)){ s2<-as.data.frame(s2) }
    
    #* `Reorder columns if they are not in the numeric order`
    histCols<-colnames(s2)
    histCols_bin <- as.numeric(sub("[a-zA-Z_.]+","",colnames(s2)))
    bins_order<-sort(histCols_bin, index.return=TRUE)$ix
    s2<-s2[,bins_order]
    
    #* `Turn s2 matrix into a vector`
    X2<-rep(histCols_bin[bins_order], as.numeric(round(colSums(s2))) )
    
    #* `Get Mean, Variance, SE, and DF from s2`
    n2 <- nrow(s2) # n samples
    m2 = mean(X2) # xbar
    s2_2 = var(X2) # s^2
    
    v2 = priors$n[2] - 1 # prior DF
    n2_n = priors$n[2] + n2 # total N including prior
    m2_n = (n2*m2 + priors$n[2]*priors$mu[2])/n2_n # weighted mean of prior and data
    v2_n = v2 + n2 # degrees of freedom including data
    s2_2_n = ((n2-1)*s2_2 + v2*priors$s2[2] + priors$n[2]*n2*(priors$mu[2] - m2)^2/n2_n)/v2_n # pooled variance
    se2<-sqrt(s2_2_n/n2_n) # standard error of the mean
    
    post_mean_d2 <- extraDistr::dlst(support,v2_n,m2_n, se2)
    post_mean_pdf2 <- post_mean_d2/sum(post_mean_d2)
    hde2_mean <- m2_n
    hdi2_mean <- m2_n + qt(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),v2_n)*se2
    
    out$summary<-data.frame(HDE_1=hde1_mean, HDI_1_low = hdi1_mean[1], HDI_1_high = hdi1_mean[2],
                            HDE_2=hde2_mean, HDI_2_low = hdi2_mean[1], HDI_2_high = hdi2_mean[2])
    out$posterior$mu <- c(m1_n, m2_n)
    out$posterior$n <- c(n1_n, n2_n)
    out$posterior$s2 <- c(s2_1_n, s2_2_n)
    #* `Save data for plotting`
    if(plot){
      out$plot_df <- rbind(out$plot_df, data.frame("range"=support,
                                "prob"=post_mean_pdf2,
                                "sample"=rep("Sample 2",length(support)))
      )
    }
    if(!all(c(as.matrix(s1),as.matrix(s2))==0)){
      post.prob.out <-  .post.prob.from.pdfs(post_mean_pdf1, post_mean_pdf2, hypothesis)
      post.prob <-post.prob.out$post.prob
      dirSymbol <-post.prob.out$direction
    }else{
      post.prob <- 1
    }
    out$summary$hyp = hypothesis
    out$summary$post.prob = post.prob
    out$dirSymbol = dirSymbol
  }
  #* `ROPE Comparison`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = extraDistr::rlst(10000,v1_n, m1_n, se1)
      if(!is.null(s2)){
      post2 = extraDistr::rlst(10000,v2_n, m2_n, se2)
      posterior = post1 - post2
      } else {
        posterior = post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  
  return(out)

}







#' ***********************************************************************************************
#' *************** `Beta (multi value)` *******************************************************
#' ***********************************************************************************************


#' @description
#' Internal function for calculating \alpha and \beta of a distribution represented by multi value traits.
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number between 0.0001 and 0.9999 representing the "bin".
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' makeMvBeta<-function(n=100,a,b){
#'   setNames(data.frame(matrix(hist(rbeta(2000,a,b),
#'   breaks=seq(0,1,length.out=n), plot=FALSE)$counts, nrow=1)),paste0("b0.",1:(n-1)))
#' }
#' 
#' mv_beta<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvBeta(n=100, a=5, b=8 )})),
#'                do.call(rbind, lapply(1:30, function(i){makeMvBeta(n=100, a=10, b=3 )})))
#' 
#' .conj_beta_mv(s1 = mv_beta[1:30,], s2= mv_beta[31:60,], priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89,
#'               cred.int.level = 0.89, hypothesis="equal")
#' @keywords internal
#' @noRd
.conj_beta_mv<-function(s1 = NULL, s2= NULL, priors = NULL,
                        plot=FALSE, rope_range = NULL, rope_ci = 0.89,
                        cred.int.level = 0.89, hypothesis="equal"){
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list(a=c(0.5,0.5),b=c(0.5,0.5))
  }
  #* `Define dense Support`
  beta_support <-seq(0.0001, 0.9999, 0.0001)
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
  dens1 <- dbeta(beta_support, a1_prime, b1_prime)
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
  
  #* `keep data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=beta_support, "prob"=pdf1, "sample"=rep("Sample 1",length(beta_support) ))
  }
  if(!is.null(s2)){
    #* `Standardize sample 2 class and names`
    if(is.null(colnames(s2))){
      bins<-(1:ncol(s2))/100
      colnames(s2)<-paste0("b", bins )
      warning(paste0("Assuming unnamed columns represent bins from ", min(bins), " to ", max(bins)))
    }
    if(is.matrix(s2)){ s2<-as.data.frame(s2)  }
    
    #* `Reorder columns if they are not in the numeric order`
    histCols<-colnames(s2)
    histCols_bin <- as.numeric(sub("[a-zA-Z_.]+","",colnames(s2)))
    bins_order<-sort(histCols_bin, index.return=TRUE)$ix
    s2<-s2[,bins_order]
    
    #* `Turn matrix into a vector`
    X2<-rep(histCols_bin[bins_order], as.numeric(round(colSums(s2))) )
    
    #* `Method of Moments for Alpha and Beta`
    mu2 <- mean(X2) 
    nu2 <- var(X2)/(nrow(s2)-1)
    alpha2 <- mu2*((mu2*(1-mu2))/(nu2) - 1)
    beta2 <- (1-mu2)*((mu2*(1-mu2))/(nu2) - 1)
    a2_prime <- priors$a[2] + alpha2
    b2_prime <- priors$b[2] + beta2
    dens2 <- dbeta(beta_support, a2_prime, b2_prime)
    pdf2 <- dens2/sum(dens2)
    hdi2 <- qbeta(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),a2_prime,b2_prime)
    if(a2_prime <= 1 & b2_prime > 1){
      hde2 <- 0
    }else if(a2_prime > 1 & b2_prime <= 1){
      hde2 <- 1
    }else{
      hde2 <- (a2_prime-1)/(a2_prime+b2_prime-2)
    }
    out$summary <- data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2], HDE_2=hde2, HDI_2_low = hdi2[1], HDI_2_high = hdi2[2])
    out$posterior$a <- c(a1_prime, a2_prime)
    out$posterior$b <- c(b1_prime, b2_prime)
    
    if(!all(c(as.matrix(s1),as.matrix(s2))==0)){
      post.prob.out <-  .post.prob.from.pdfs(pdf1, pdf2, hypothesis)
      post.prob <-post.prob.out$post.prob
      dirSymbol <-post.prob.out$direction
    }else{
      post.prob <- 1
    }
    out$summary$hyp = hypothesis
    out$summary$post.prob = post.prob
    out$dirSymbol = dirSymbol
  
    #* `Keep data from s2 for plotting`
    if(plot){
      out$plot_df <- rbind(out$plot_df,data.frame("range"=beta_support, "prob"=pdf2, "sample"=rep("Sample 2",length(beta_support)) ))
    }
  }
  #* `Use ROPE method to test if the difference in distributions is meaningful`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 <- rbeta(10000,a1_prime,b1_prime)
      if(!is.null(s2)){
        post2 <- rbeta(10000,a2_prime,b2_prime)
        posterior <- post1 - post2
      } else {
        posterior <- post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  return(out)
}

#' ***********************************************************************************************
#' *************** `Beta (single value)` *******************************************************
#' ***********************************************************************************************


#' @description
#' Internal function for calculating \alpha and \beta of a distribution represented by multi value traits.
#' @param s1 A vector of numerics drawn from a beta distribution.
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' .conj_beta_sv(s1 = rbeta(100, 5, 10), s2= rbeta(100, 10, 5), priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=FALSE, rope_range = c(-0.1, 0.1), rope_ci = 0.89, 
#'               cred.int.level = 0.89, hypothesis="equal")
#' @keywords internal
#' @noRd
.conj_beta_sv<-function(s1 = NULL, s2= NULL, priors = NULL,
                       plot=FALSE, rope_range = NULL, rope_ci = 0.89,
                       cred.int.level = 0.89, hypothesis="equal"){

  if(any(c(s1, s2)>1)){stop("Values above 1 cannot be used with the beta distribution")}
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list(a=c(0.5,0.5),b=c(0.5,0.5))
  }
  #* `Define dense Support`
  beta_support <-seq(0.0001, 0.9999, 0.0001)
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
  dens1 <- dbeta(beta_support, a1_prime, b1_prime)
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
  
  #* `keep data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=beta_support, "prob"=pdf1, "sample"=rep("Sample 1",length(beta_support) ))
  }
  #* `get parameters for s2 if it exists using method of moments`
  if(!is.null(s2)){
    mu2 <- mean(s2) 
    nu2 <- var(s2)/(length(s2)-1)
    alpha2 <- mu2*((mu2*(1-mu2))/(nu2) - 1)
    beta2 <- (1-mu2)*((mu2*(1-mu2))/(nu2) - 1)
    a2_prime <- priors$a[2] + alpha2
    b2_prime <- priors$b[2] + beta2
    dens2 <- dbeta(beta_support, a2_prime, b2_prime)
    pdf2 <- dens2/sum(dens2)
    hdi2 <- qbeta(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),a2_prime,b2_prime)
    if(a2_prime <= 1 & b2_prime > 1){
      hde2 <- 0
    }else if(a2_prime > 1 & b2_prime <= 1){
      hde2 <- 1
    }else{
      hde2 <- (a2_prime-1)/(a2_prime+b2_prime-2)
    }
    out$summary <- data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2], HDE_2=hde2, HDI_2_low = hdi2[1], HDI_2_high = hdi2[2])
    out$posterior$a <- c(a1_prime, a2_prime)
    out$posterior$b <- c(b1_prime, b2_prime)
    
    if(!all(c(s1,s2)==0)){
      post.prob.out <-  .post.prob.from.pdfs(pdf1, pdf2, hypothesis)
      post.prob <-post.prob.out$post.prob
      dirSymbol <-post.prob.out$direction
      } else{
      post.prob <- 1
    }
    out$summary$hyp = hypothesis
    out$summary$post.prob = post.prob
    out$dirSymbol = dirSymbol
    #* `Keep data from s2 for plotting`
    if(plot){
      out$plot_df <- rbind(out$plot_df,data.frame("range"=beta_support, "prob"=pdf2, "sample"=rep("Sample 2",length(beta_support)) ))
    }
  }
  #* `Use ROPE method to test if the difference in distributions is meaningful`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = rbeta(10000,a1_prime,b1_prime)
      if(!is.null(s2)){
        post2 = rbeta(10000,a2_prime,b2_prime)
        posterior = post1 - post2 
      } else {
        posterior = post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  return(out)
}






#' ***********************************************************************************************
#' *************** `Poisson (single value)` *******************************************************
#' ***********************************************************************************************

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
#' .conj_poisson_sv(s1 = rpois(20, 10), s2= rpois(20, 8), priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'               plot=FALSE, rope_range = c(-1,1), rope_ci = 0.89, 
#'               cred.int.level = 0.89, hypothesis="equal")
#' @keywords internal
#' @noRd

.conj_poisson_sv<-function(s1 = NULL, s2= NULL, priors = NULL,
                        plot=FALSE, rope_range = NULL, rope_ci = 0.89, 
                        cred.int.level = 0.89, hypothesis="equal", support=NULL){
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
  
  #* `Use conjugate gamma prior on lambda``
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
  
  #* `keep data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=support, "prob"=pdf1, "sample"=rep("Sample 1",length(support) ))
  }
  #* `get parameters for s2 if it exists using method of moments`
  if(!is.null(s2)){
    a2_prime = priors$a[2]+sum(s2)
    b2_prime = priors$b[2]+length(s2)
    dens2 <- dgamma(support, a2_prime, b2_prime)
    pdf2 <- dens2/sum(dens2)
    hdi2 <- qgamma(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),a2_prime,b2_prime)
    if(a2_prime <= 1 & b2_prime > 1){
      hde2 <- 0
    }else if(a2_prime > 1 & b2_prime <= 1){
      hde2 <- Inf
    }else{
      hde2 <- (a2_prime-1)/b2_prime 
    }
    out$summary <- data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2],
                              HDE_2=hde2, HDI_2_low = hdi2[1], HDI_2_high = hdi2[2])
    out$posterior$a <- c(a1_prime, a2_prime)
    out$posterior$b <- c(b1_prime, b2_prime)
    
    if(!all(c(s1,s2)==0)){
      post.prob.out <-  .post.prob.from.pdfs(pdf1, pdf2, hypothesis)
      post.prob <-post.prob.out$post.prob
      dirSymbol <-post.prob.out$direction
    } else{
      post.prob <- 1
    }
    out$summary$hyp = hypothesis
    out$summary$post.prob = post.prob
    out$dirSymbol = dirSymbol
    #* `Keep data from s2 for plotting`
    if(plot){
      out$plot_df <- rbind(out$plot_df,data.frame("range"=support, "prob"=pdf2, "sample"=rep("Sample 2",length(support)) ))
    }
  }
  #* `Use ROPE method to test if the difference in distributions is meaningful`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = rgamma(10000,a1_prime,b1_prime)
      if(!is.null(s2)){
        post2 = rgamma(10000,a2_prime,b2_prime)
        posterior = post1 - post2 
      } else {
        posterior = post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  return(out)
}




#' ***********************************************************************************************
#' *************** `Neg Binomial (single value)` *******************************************************
#' ***********************************************************************************************

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
#' @keywords internal
#' @noRd

.conj_negbin_sv<-function(s1 = NULL, s2= NULL, priors = NULL,
                           plot=FALSE, rope_range = NULL, rope_ci = 0.89, 
                           cred.int.level = 0.89, hypothesis="equal", support=NULL){

  #* `Check samples`
  if(any(abs(c(s1,s2)-round(c(s1,s2)))>.Machine$double.eps^0.5) | any(c(s1,s2)<0) ){stop("Only positive whole numbers can be used in the Negative Binomial distribution")}
  #* `make default prior if none provided`
  if(is.null(priors)){
    priors <- list(r=c(10, 10), a=c(0.5,0.5),b=c(0.5,0.5)) # beta prior on P
    warning("True value of r for negative binomial distribution has defaulted to 10, you should add a prior including r parameter.")
  }
  #* `Define dense Support if missing`
  if(is.null(support)){
    upper <- max(c(s1,s2))*2
    support <-seq(0, upper, 0.01)
  }
  out <- list()
  
  #* `Use conjugate beta prior on probability`
  #* Note that this is very sensitive to the R value being appropriate
  a1_prime <- priors$a[1] + priors$r[1] * length(s1)
  b1_prime <- priors$b[1] + sum(s1) 
  
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
  
  #* `keep data for plotting`
  if(plot){
    out$plot_df <- data.frame("range"=support, "prob"=pdf1, "sample"=rep("Sample 1",length(support) ))
  }
  #* `get parameters for s2 if it exists using method of moments`
  if(!is.null(s2)){
    #* `Use conjugate beta prior on probability`
    a2_prime <- priors$a[2] + priors$r[2] * length(s2)
    b2_prime <- priors$b[2] + sum(s2) 
    
    dens2 <- dbeta(support, a2_prime, b2_prime)
    pdf2 <- dens2/sum(dens2)
    hdi2 <- qbeta(c((1-cred.int.level)/2, (1-((1-cred.int.level)/2))),a2_prime,b2_prime)
    if(a2_prime <= 1 & b2_prime > 1){
      hde2 <- 0
    }else if(a2_prime > 1 & b2_prime <= 1){
      hde2 <- 1
    }else{
      hde2 <- (a2_prime-1)/(a2_prime+b2_prime-2) # ~ p2 ~ xbar2/var2
    }
    out$summary <- data.frame(HDE_1=hde1, HDI_1_low = hdi1[1], HDI_1_high = hdi1[2],
                              HDE_2=hde2, HDI_2_low = hdi2[1], HDI_2_high = hdi2[2])
    out$posterior$r <- priors$r
    out$posterior$a <- c(a1_prime, a2_prime)
    out$posterior$b <- c(b1_prime, b2_prime)
    
    if(!all(c(s1,s2)==0)){
      post.prob.out <-  .post.prob.from.pdfs(pdf1, pdf2, hypothesis)
      post.prob <-post.prob.out$post.prob
      dirSymbol <-post.prob.out$direction
    } else{
      post.prob <- 1
    }
    out$summary$hyp = hypothesis
    out$summary$post.prob = post.prob
    out$dirSymbol = dirSymbol
    #* `Keep data from s2 for plotting`
    if(plot){
      out$plot_df <- rbind(out$plot_df,data.frame("range"=support, "prob"=pdf2, "sample"=rep("Sample 2",length(support)) ))
    }
  }
  #* `Use ROPE method to test if the difference in distributions is meaningful`
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = rbeta(10000,a1_prime,b1_prime)
      if(!is.null(s2)){
        post2 = rbeta(10000,a2_prime,b2_prime)
        posterior = post1 - post2 
      } else {
        posterior = post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_ci ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi_diff[1], HDI_rope_high = hdi_diff[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        out$rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  return(out)
}


#' ***********************************************************************************************
#' *************** `Dirichlet (multi value)` *******************************************************
#' ***********************************************************************************************

#' @description
#' 
#' The dirichlet distribution is the conjugate prior for a multinomial distribution and can be used to 
#' describe an image histogram without turning the data from discrete to continuous. 
#' This is agnostic to the appearance of the "curve" shown by the histogram at the expense of not providing 
#' single HDE and HDI values. HDE and HDI for the dirichlet are vectors of length n = n_bins.
#' Note that by default the HDE and HDI for each sample are not returned. Those can be estimated as shown in examples.
#' 
#' Note that this returns a data table for convenience in printing the list columns if rope_range is specified.
#' 
#' 
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number representing the "bin".
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' 
#' makeMvLn<-function(bins=500,mu_log,sigma_log){
#' setNames(data.frame(matrix(hist(rlnorm(2000,mu_log, sigma_log),
#'                                 breaks=seq(1,bins,5), plot=FALSE)$counts, nrow=1)),
#'                                          paste0("b",seq(1,bins,5))[-1] ) }
#' set.seed(123) 
#' mv_ln<-rbind(do.call(rbind,
#'                       lapply(1:30, function(i){makeMvLn(mu_log=log(130),
#'                                           sigma_log=log(1.3) )})),
#'             do.call(rbind, lapply(1:30, function(i){makeMvLn(mu_log=log(100),
#'                                           sigma_log=log(1.2) )})))
#' mv_ln$group = rep(c("a", "b"), each = 30)
#' s1 <- mv_ln[mv_ln$group=="a", 1:99]
#' s2 <- mv_ln[mv_ln$group=="b", 1:99]
#' out <- .conj_diri_mv_1(s1, s2, priors=NULL, plot=TRUE, rope_range = c(-0.1, 0.1),
#'       rope_ci = 0.89, cred.int.level = 0.89, hypothesis="equal")
#' dim(out$summary) # summary is still one row dataframe
#' # post.prob 0.542
#' # mean_rope_prob = 1
#' dim(out$summary$HDI_rope[[1]]) # with nested dataframes.
#' 
#' # Calculating HDI and HDE for samples
#' 
#' if(FALSE){
#' 
#'   HDI <- as.data.frame(bayestestR::hdi(as.data.frame(extraDistr::rdirichlet(10000, alpha_vector)), ci = cred.int.level))
#'   HDE_pre <- as.data.frame(bayestestR::hdi(as.data.frame(extraDistr::rdirichlet(10000, alhpa_vector)), ci = 0.001))
#'   HDE <- data.frame("Parameter" = HDE_pre[["Parameter"]],
#'                        "HDE" = rowMeans(HDE_pre[,c("CI_low", "CI_high")]))
#' }
#' 
#' @keywords internal
#' @noRd

.conj_diri_mv_1<-function(s1 = NULL, s2= NULL,
                         priors=NULL,
                         plot=FALSE, rope_range = NULL, rope_ci = 0.89,
                         cred.int.level = 0.89, hypothesis="equal"){
  out <-list()
  
  if(is.null(priors)){
    priors <- list(alpha1 = rep(1, ncol(s1)), # alpha vector priors
                   alpha2 = rep(1, ncol(s1)),
                   prec1_a = 0.5, # gamma priors on precision
                   prec1_b = 0.5,
                   prec2_a = 0.5,
                   prec2_b = 0.5)
  }
  
  if(!is.null(s1)){
    alpha1_prime <- priors$alpha1 + colSums(s1) # updating prior with sum of samples in each bin
    
    precision1 <- sum(alpha1_prime) # reparameterize alpha1 to mean and precision
    mean1 <- alpha1_prime / precision1
    
    prec1_prime_a <- priors$prec1_a + sum(colSums(s1)) # A updates as A' = A + sum(obs)
    prec1_prime_b <- priors$prec1_b + nrow(s1) # B updates as B' = B + n(obs)
    
    out$summary <- data.frame('HDE_1' = NA, 'HDI_1_low' = NA, 'HDI_1_high' = NA)
    
    out$posterior$mean1 = mean1
    out$posterior$precision1 = precision1
    out$posterior$prec1_a <- prec1_prime_a
    out$posterior$prec1_b <- prec1_prime_b
    
    if(plot){
      out$plot_df <- data.frame("bin" = 1:ncol(s1), "prob"=mean1, "sample"=rep("Sample 1",ncol(s1) ))
    }
    
  } else{stop("s1 is required")}
  
  if(!is.null(s2)){
    if(ncol(s1)!=ncol(s2)){stop("s1 and s2 must have the same number of bins (columns)")}
    alpha2_prime <- priors$alpha2 + colSums(s2) # updating prior with sum of samples in each bin
    
    precision2 <- sum(alpha2_prime) # reparameterize alpha1 to mean and precision
    mean2 <- alpha2_prime / precision2
    
    prec2_prime_a <- priors$prec2_a + sum(colSums(s2)) # A updates as A' = A + sum(obs)
    prec2_prime_b <- priors$prec2_b + nrow(s2) # B updates as B' = B + n(obs)
    
    out$summary <- cbind(out$summary, data.frame('HDE_2' = NA, 'HDI_2_low' = NA, 'HDI_2_high' = NA))
    
    out$posterior$mean2 = mean2
    out$posterior$precision2 = precision2
    out$posterior$prec2_a <- prec2_prime_a
    out$posterior$prec2_b <- prec2_prime_b
    
    if(plot){
      out$plot_df <- rbind(out$plot_df, data.frame("bin" = 1:ncol(s2), "prob"=mean2, "sample"=rep("Sample 2",ncol(s2) )))
    }
    
    #* `test hypothesis`
    post.prob.out <-  .post.prob.from.pdfs(mean1, mean2, hypothesis)
    post.prob <-post.prob.out$post.prob
    dirSymbol <-post.prob.out$direction
    
    out$summary <- cbind(out$summary, data.frame('hyp' = hypothesis, 'post.prob' = post.prob))
    out$dirSymbol <- dirSymbol
  }
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      
      #* `difference of mean vectors`
      
      if(!is.null(s2)){
        mean_diff <- mean1 - mean2
        direction <- mean1 >= mean2
      } else {
        mean_diff <- mean1
        direction <- rep(TRUE, length(mean1))
      }
      
      abs_mean_diff <- abs(mean_diff)
      norm_abs_mean_diff <- abs_mean_diff/sum(abs_mean_diff)
      norm_mean_diff <- ifelse(direction, 1, -1) * norm_abs_mean_diff
      
      #* `draw N times from mean vector 1`
      mean1_draws <- extraDistr::rdirichlet(10000, mean1)
      
      if(!is.null(s2)){
        #* `draw N times from mean vector 2`
        mean2_draws <- extraDistr::rdirichlet(10000, mean2)
        posterior <- mean1_draws - mean2_draws
      } else{
        posterior = mean1_draws
      }
      
      rope_probs<-lapply(1:ncol(posterior), function(bin){
        as.numeric(bayestestR::rope(posterior[,bin], range = rope_range, ci_method = "HDI", ci=rope_ci))
      })
      
      mean_rope_prob = sum( unlist(rope_probs)* norm_abs_mean_diff )
      
      HDI_rope <- as.data.frame(bayestestR::hdi(as.data.frame(posterior), ci = rope_ci))
      HDE_rope_pre <- as.data.frame(bayestestR::hdi(as.data.frame(posterior), ci = 0.001))
      HDE_rope <- data.frame("Parameter" = HDE_rope_pre[["Parameter"]],
                             "HDE" = rowMeans(HDE_rope_pre[,c("CI_low", "CI_high")]))
      
      if(plot){
        cis<-c(seq(0.1,rope_ci, 0.2), rope_ci)
        rope_df <- as.data.frame(bayestestR::hdi(as.data.frame(posterior), ci = cis ))
        rope_df$bin <- as.numeric(sub("V","", rope_df$Parameter))
        out$rope_df <- rope_df
      }
      
      out$summary = cbind(out$summary, data.frame('HDE_rope' = list(1), "HDI_rope" = list(1),
                                                  "rope_probs" = list(1), "mean_rope_prob" = mean_rope_prob))
      out$summary <- data.table::as.data.table(out$summary)
      out$summary$HDE_rope <- list(HDE_rope)
      out$summary$HDI_rope <- list(HDI_rope)
      out$summary$rope_probs <- list(rope_probs)
      #out$summary <- as.data.frame(out$summary)
    } else{
      stop("rope must be a vector of length 2")
    }
  }
  return(out)
}


#' ***********************************************************************************************
#' *************** `Calculate Posterior Probability given PDFs` ***********************************
#' ***********************************************************************************************


#' @keywords internal
#' @noRd
.post.prob.from.pdfs<-function(pdf1, pdf2, hypothesis){
  if(hypothesis == "unequal"){
    post.prob <- 1-sum(apply(cbind(pdf1, pdf2),MARGIN=1,function(i) min(i)),na.rm=TRUE)
    dirSymbol="!="
  } else if(hypothesis == "equal"){
    post.prob <- sum(apply(cbind(pdf1, pdf2),MARGIN=1,function(i) min(i)),na.rm=TRUE)
    dirSymbol="="
  }else if(hypothesis == "lesser"){
    direction <- pdf1 <= pdf2
    post.prob <- sum(pdf1 * direction,na.rm=TRUE)
    dirSymbol="<"
  }else if(hypothesis == "greater"){
    direction <- pdf1 >= pdf2
    post.prob <- sum(pdf1 * direction,na.rm=TRUE)
    dirSymbol=">"
  }else{
    stop("hypothesis must be either unequal, equal, lesser, or greater")
  }
  return(list('post.prob'=post.prob, "direction" = dirSymbol))
}

#' ***********************************************************************************************
#' *************** `General Plotting Function` ***********************************
#' ***********************************************************************************************


#' @keywords internal
#' @noRd
.conj_general_plot <- function(res, s2, rope_range, support, rope_ci){
    p <- ggplot2::ggplot(res$plot_df, ggplot2::aes(x=.data$range, y=.data$prob))+
      ggplot2::geom_area(data=res$plot_df[res$plot_df$sample == "Sample 1",],fill="red",alpha=0.5)+
      ggplot2::geom_vline(ggplot2::aes(xintercept=res$summary$HDI_1_low ),color="red",linewidth=1.1)+
      ggplot2::geom_vline(ggplot2::aes(xintercept=res$summary$HDE_1), color="red",linetype="dashed",linewidth=1.1)+
      ggplot2::geom_vline(ggplot2::aes(xintercept=res$summary$HDI_1_high),color="red",linewidth=1.1)+
      ggplot2::labs(x="Posterior Distribution of Random Variable", y="Density", title = "Distribution of Samples",
                    subtitle=paste0("HDE: ",round(res$summary$HDE_1,2),
                                    "\nHDI: [",round(res$summary$HDI_1_low,2),", ",
                                    round(res$summary$HDI_1_high,2),"]"))+
      pcv_theme()
    
    if(!is.null(s2)){
      
      if(res$summary$post.prob<1e-5){post.prob.text = "<1e-5"} else { post.prob.text<-round(res$summary$post.prob,5)}
      
      p <- p +
        ggplot2::geom_area(data=res$plot_df[res$plot_df$sample == "Sample 2",],fill="blue",alpha=0.5)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=res$summary$HDI_2_low),color="blue",linewidth=1.1)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=res$summary$HDE_2),color="blue",linetype="dashed",linewidth=1.1)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=res$summary$HDI_2_high),color="blue",linewidth=1.1)+
        ggplot2::labs(subtitle = paste0(
          "Sample 1:  ",round(res$summary$HDE_1, 2)," [",round(res$summary$HDI_1_low,2),", ",round(res$summary$HDI_1_high,2),"]\n",
          "Sample 2:  ",round(res$summary$HDE_2,2)," [",round(res$summary$HDI_2_low,2),", ",round(res$summary$HDI_2_high,2),"]\n",
          "P[p1",res$dirSymbol,"p2] = ",post.prob.text))
      res<-res[-which(names(res)=="dirSymbol")]
    }
    #* `make x axis range if using default support`
    if(is.null(support)){
      if(!is.null(s2)){
        x_lower<-min( res$summary[,c("HDI_1_low", "HDI_2_low")] )/1.1
        x_upper<-max( res$summary[,c("HDI_1_high", "HDI_2_high")] )*1.1
      } else{
        x_lower<-res$summary$HDI_1_low / 1.1
        x_upper<-res$summary$HDI_1_high * 1.1
      }
      p<-p+ggplot2::coord_cartesian(xlim=c(x_lower, x_upper))
    }
    
    if(length(rope_range) == 2){
      p <- p + ggplot2::ggplot(res$rope_df, ggplot2::aes(x=.data$X))+
        ggplot2::geom_histogram(bins=100,fill="purple",color="purple", alpha=0.7)+
        ggplot2::geom_histogram(data=data.frame("X"=res$rope_df[res$rope_df$X > rope_range[1] & res$rope_df$X < rope_range[2] &
                                                                  res$rope_df$X > res$summary$HDI_rope_low & res$rope_df$X < res$summary$HDI_rope_high,]),
                                bins=100,fill="gray30",color="gray30")+
        ggplot2::geom_segment(ggplot2::aes(x=rope_range[1],xend=rope_range[2],y=0,yend=0),linewidth=2,color="gray70")+
        ggplot2::geom_vline(ggplot2::aes(xintercept=res$summary$HDI_rope_low),linewidth=0.7)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=res$summary$HDE_rope),linetype="dashed",linewidth=0.7)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=res$summary$HDI_rope_high),linewidth=0.7)+
        ggplot2::labs(x="Posterior of Difference", y="Frequency", title = "Distribution of Difference",
                      subtitle =  paste0(
                        "Median Difference of ",round(res$summary$HDE_rope,2),"\n",
                        100*rope_ci,"% CI [",round(res$summary$HDI_rope_low,2),", ",round(res$summary$HDI_rope_high,2),"]\n",
                        rope_ci,"% HDI in [", rope_range[1],", ",rope_range[2], "]: ",round(res$summary$rope_prob,2)))+
        pcv_theme()+
        ggplot2::theme(axis.title.y.left = ggplot2::element_blank(), axis.text.y.left = ggplot2::element_blank())+
        patchwork::plot_layout(widths = c(2,1))
    }
    plot(p)
    res$plot<-p
  return(res)
}

#' ***********************************************************************************************
#' *************** `Dirichlet Plotting Function` ***********************************
#' ***********************************************************************************************


#' @keywords internal
#' @noRd

.conj_diri_plot <- function(res, s2, rope_range, rope_ci){
  
  
  p <- ggplot2::ggplot(res$plot_df, ggplot2::aes(x=.data$bin, y=.data$prob))+
    ggplot2::geom_col(data=res$plot_df[res$plot_df$sample == "Sample 1",], fill="red", alpha=0.5, position="identity")+
    ggplot2::labs(y="Density", title = "Distribution of Samples")+
    pcv_theme()+
    ggplot2::theme(legend.position=c(0.9,0.9),
          legend.title = ggplot2::element_blank(),
          axis.title.x.bottom = ggplot2::element_blank())
  
  if(!is.null(s2)){
    
    if(res$summary$post.prob<1e-5){post.prob.text = "<1e-5"} else { post.prob.text<-round(res$summary$post.prob,5)}
    
    p <- ggplot2::ggplot(res$plot_df, ggplot2::aes(x=.data$bin, y=.data$prob, fill=.data$sample))+
      ggplot2::geom_col(alpha=0.75, position="identity")+
      ggplot2::labs(y="Density", title = "Distribution of Samples",
                    subtitle = paste0("P[p1",res$dirSymbol,"p2] = ",post.prob.text))+
      pcv_theme()+
      ggplot2::theme(legend.position=c(0.9,0.9),
                     legend.title = ggplot2::element_blank(),
                     axis.title.x.bottom = ggplot2::element_blank(), 
                     legend.key.size = ggplot2::unit(0.25, "cm"))
  }
  
  if(length(rope_range) == 2){
    if(!is.null(s2)){
      p <- p +
        ggplot2::labs(subtitle = paste0("P[p1",res$dirSymbol,"p2] = ",post.prob.text, "\n",
                                        rope_ci,"% HDI in [", rope_range[1],", ",rope_range[2], "]: ",
                                        round(res$summary$mean_rope_prob, 2) ))
      res<-res[-which(names(res)=="dirSymbol")]
    }
    
    xLims <- ggplot2::layer_scales(p)$x$range$range
    
    rect_width <- min(diff(unique(res$rope_df$bin)))/2.5
    
    #* this makes the rope probability look incorrect because it does not accurately show the distribution,
    #* it just makes a bar for the HDI, but density within the HDI varies widely for my examples.
    
    rdf <- res$rope_df
    cis<-rev(unique(rdf$CI))
    virPal <- viridis::plasma(length(cis), direction=-1)
    
    rope_plot <- ggplot2::ggplot(rdf, ggplot2::aes(xmin=bin-rect_width, xmax = bin+rect_width,
                                              ymin=CI_low, ymax=CI_high ))+
      lapply(1:length(cis), function(i){
        ggplot2::geom_rect(data = rdf[rdf$CI == cis[[i]], ], ggplot2::aes(fill = as.character(cis[[i]]) ))
      })+
      ggplot2::geom_rect(data=data.frame(x=0,y=0), xmin = -Inf, xmax=Inf, 
                         ymin = min(rope_range), ymax=max(rope_range),
                         fill="gray100", alpha=0.5)+
      lapply(rope_range, function(i){
        ggplot2::geom_hline(yintercept = i, linetype=5, color="gray40", linewidth=0.25)
      })+
      ggplot2::guides(fill=ggplot2::guide_legend(nrow=1, override.aes = list(color=NA)))+
      ggplot2::labs(x="Posterior Distribution of Random Variable")+
      pcv_theme()+
      ggplot2::theme(legend.position=c(0.75, 0.9),
                     legend.title = element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(), 
                     legend.key.size = ggplot2::unit(0.25, "cm"))+
      ggplot2::scale_x_continuous(limits = xLims)+
      ggplot2::scale_fill_manual(values = virPal)
    
  
    layout<-c(patchwork::area(1,1,3,6),
              patchwork::area(4,1,4,6))
    
    p <- p / rope_plot + patchwork::plot_layout(design = layout)
  }
  plot(p)
  
  res$plot<-p
  return(res)
}


