test_that("conjugate single value T works", {
  s1 = c(43.8008289810423, 44.6084228775479, 68.9524219823026, 77.442231894233, 
         45.2302703709121, 53.8005757403944, 33.8292993277826, 59.7018653972819, 
         57.8433869136181, 52.8224097917938) # dput(rnorm(10, 50,10))
  s2 = c(84.7860854952772, 53.38097452501, 52.352235256613, 49.2369049504088, 
         72.7625716991815, 62.6982283802374, 61.2347595388326, 45.298878516913, 
         39.6312400911458, 66.9134811003628) # dput(rnorm(10, 60,12))
  out <- conjugate(s1=s1, s2=s2, method="t",
                          priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
                          plot=FALSE, rope_range = c(-8,8), rope_ci = 0.89,
                          cred.int.level = 0.89, hypothesis="equal", support=NULL)
  
  expect_equal(out$summary$post.prob, 0.7260736, tolerance = 1e-6)
  
  expect_equal(out$summary$rope_prob, 0.6, tolerance = 0.025)
  
  expect_equal(names(out), c("summary", "posterior") )
})

test_that("conjugate multi value T works", {
  makeMvGauss<-function(bins=180,mu,sigma){
     setNames(data.frame(matrix(hist(rnorm(2000,mu, sigma),
      breaks=seq(1,bins,1), plot=FALSE)$counts, nrow=1)),
     paste0("b",1:(bins-1) ))
  }
  set.seed(123)
  mv_gauss<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=50, sigma=10 )})),
                  do.call(rbind, lapply(1:40, function(i){makeMvGauss(bins=180, mu=60, sigma=12 )})))
  
  out <- conjugate(s1=mv_gauss[1:30,], s2= mv_gauss[31:70,], method="t",
                          priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
                          plot=FALSE, rope_range = c(-5,5), rope_ci = 0.89,
                          cred.int.level = 0.89, hypothesis="equal", support=NULL)
  
  expect_equal(out$summary$post.prob, 0.03881883, tolerance = 1e-6)
  
  expect_equal(out$summary$rope_prob, 0.015, tolerance = 0.02)
  
  expect_equal(names(out), c("summary", "posterior") )
})

test_that("conjugate single value gaussian works", {
  s1 = c(43.8008289810423, 44.6084228775479, 68.9524219823026, 77.442231894233, 
         45.2302703709121, 53.8005757403944, 33.8292993277826, 59.7018653972819, 
         57.8433869136181, 52.8224097917938) # dput(rnorm(10, 50,10))
  s2 = c(62.2654459133762, 66.6863571733485, 61.2951438574251, 62.0014980341704, 
         44.0772327229333, 56.169510174076, 71.1378538738675, 55.7547954794673, 
         52.4202653287144, 63.3091644583334, 49.263640809148, 63.2460598779059, 
         60.3804997092304, 25.1210401427447, 42.6563192857856) # dput(rnorm(15, 60,12))
  
  out <- conjugate(s1=s1, s2= s2, method = "gaussian",
                    priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
                    plot=FALSE, rope_range = c(-10, 10), rope_ci = 0.89,
                    cred.int.level = 0.89, hypothesis="equal", support=seq(-20,120, 0.01))
  
  
  expect_equal(out$summary$post.prob, 0.9106556, tolerance = 1e-6)
  
  expect_equal(out$summary$rope_prob, 0.3123245, tolerance = 0.025)
  
  expect_equal(names(out), c("summary", "posterior") )
})

test_that("conjugate multi value T works", {
  makeMvGauss<-function(bins=180,mu,sigma){
    setNames(data.frame(matrix(hist(rnorm(2000,mu, sigma),
                                    breaks=seq(1,bins,1), plot=FALSE)$counts, nrow=1)),
             paste0("b",1:(bins-1) ))
  }
  set.seed(123)
  mv_gauss<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=50, sigma=10 )})),
                  do.call(rbind, lapply(1:40, function(i){makeMvGauss(bins=180, mu=60, sigma=12 )})))
  
  out <- conjugate(s1=mv_gauss[1:30,], s2= mv_gauss[31:70,], method="gaussian",
                   priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
                   plot=FALSE, rope_range = c(-5,5), rope_ci = 0.89,
                   cred.int.level = 0.89, hypothesis="equal", support=NULL)
  
  # summary(unlist(lapply(1:100, function(i){
  #   #set.seed(123)
  #   out <- conjugate(s1=mv_gauss[1:30,], s2= mv_gauss[31:70,], method="gaussian",
  #                    priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
  #                    plot=FALSE, rope_range = c(-5,5), rope_ci = 0.89,
  #                    cred.int.level = 0.89, hypothesis="equal")
  #   out$summary$rope_prob
  # })))
  
  expect_equal(out$summary$post.prob, 0.7179459, tolerance = 1e-6)
  
  expect_equal(out$summary$rope_prob, 0.195, tolerance = 0.015)
  
  expect_equal(names(out), c("summary", "posterior") )
})










test_that("generic plotting works", {
  s1 = c(43.8008289810423, 44.6084228775479, 68.9524219823026, 77.442231894233, 
         45.2302703709121, 53.8005757403944, 33.8292993277826, 59.7018653972819, 
         57.8433869136181, 52.8224097917938) # dput(rnorm(10, 50,10))
  s2 = c(84.7860854952772, 53.38097452501, 52.352235256613, 49.2369049504088, 
         72.7625716991815, 62.6982283802374, 61.2347595388326, 45.298878516913, 
         39.6312400911458, 66.9134811003628) # dput(rnorm(10, 60,12))
  out <- conjugate(s1=s1, s2=s2, method="t",
                   priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
                   plot=TRUE, rope_range = c(-8,8), rope_ci = 0.89,
                   cred.int.level = 0.89, hypothesis="equal", support=NULL)
  
  expect_equal(names(out), c("summary", "posterior", "plot_df", "rope_df", "plot") )
  
  expect_s3_class(out$plot, "ggplot")
  
})

test_that("dirichlet plotting version 1 works and is consistent", {
  makeMvLn<-function(bins=500,mu_log,sigma_log){
  setNames(data.frame(matrix(hist(rlnorm(2000,mu_log, sigma_log),
                                  breaks=seq(1,bins,5), plot=FALSE)$counts, nrow=1)),
                                           paste0("b",seq(1,bins,5))[-1] ) }
  set.seed(123)
  mv_ln<-rbind(do.call(rbind,
                        lapply(1:30, function(i){makeMvLn(mu_log=log(130),
                                            sigma_log=log(1.3) )})),
              do.call(rbind, lapply(1:30, function(i){makeMvLn(mu_log=log(100),
                                            sigma_log=log(1.2) )})))
  s1 <- mv_ln[1:30, ]
  s2 <- mv_ln[31:60, ]
  out <- conjugate(s1, s2, method = "dirichlet", priors=NULL, plot=TRUE,
        rope_range = c(-0.025, 0.025), rope_ci = 0.89,
        cred.int.level = 0.89, hypothesis="equal")
  
  expect_equal(out$summary$post.prob, 0.5420556, tolerance = 1e-6)
  
  expect_equal(out$summary$mean_rope_prob, 0.6407565, tolerance = 1e-6)

  expect_equal(names(out), c("summary", "posterior", "plot_df", "rope_df", "plot") )
  
  expect_s3_class(out$plot, "ggplot")
})


test_that("dirichlet plotting version 2 works and is consistent", {
  makeMvLn<-function(bins=500,mu_log,sigma_log){
    setNames(data.frame(matrix(hist(rlnorm(2000,mu_log, sigma_log),
                                    breaks=seq(1,bins,5), plot=FALSE)$counts, nrow=1)),
             paste0("b",seq(1,bins,5))[-1] ) }
  set.seed(123)
  mv_ln<-rbind(do.call(rbind,
                       lapply(1:30, function(i){makeMvLn(mu_log=log(130),
                                                         sigma_log=log(1.3) )})),
               do.call(rbind, lapply(1:30, function(i){makeMvLn(mu_log=log(100),
                                                                sigma_log=log(1.2) )})))
  s1 <- mv_ln[1:30, ]
  s2 <- mv_ln[31:60, ]
  
  out <- conjugate(s1, s2, method = "dirichlet2", priors=NULL, plot=TRUE,
                         rope_range = c(-0.025, 0.025), rope_ci = 0.89,
                         cred.int.level = 0.89, hypothesis="equal")
  
  expect_equal(out$summary$post.prob, 0.5420556, tolerance = 1e-6)
  
  expect_equal(out$summary$mean_rope_prob, 0.38, tolerance = 0.01)
  
  expect_equal(names(out), c("summary", "posterior", "plot_df", "rope_df", "plot") )
  
  expect_s3_class(out$plot, "ggplot")
})




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
#' diri_ex_2 <- conjugate(s1, s2, method = "dirichlet2", priors=NULL, plot=TRUE,
#'       rope_range = c(-0.025, 0.025), rope_ci = 0.89, 
#'       cred.int.level = 0.89, hypothesis="equal")  
#'       