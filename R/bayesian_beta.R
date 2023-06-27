
library(pcvr)
library(data.table)
sv<-read.pcv("https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv", mode="wide", singleValueOnly = T, reader="fread")
sv$genotype = substr(sv$barcode, 3,5)
sv$genotype = ifelse(sv$genotype == "002", "B73",
                     ifelse(sv$genotype == "003", "W605S",
                            ifelse(sv$genotype == "004", "MM", "Mo17")))
sv$fertilizer = substr(sv$barcode, 8, 8)
sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
                       ifelse(sv$fertilizer == "B", "50", "0"))
sv<-bw.time(sv, plantingDelay = 0, phenotype="area.pixels", cutoff=10, timeCol="timestamp", group=c("barcode", "rotation"), plot=F)
pixels_per_cmsq <- 42.5^2   # pixel per cm^2
sv$area_cm2<-sv$area.pixels / pixels_per_cmsq
sv$height_cm <- sv$height.pixels/42.5

#* generic shared plotting method
#' @export

conjugateBayesian<-function(plot, ...){
  
  #* pick method, use it

  if(plot){
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x=range, y=prob))+
      ggplot2::geom_area(data=plot_df[plot_df$sample == "Sample 1",],fill="red",alpha=0.5)+
      ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[2]] ),color="red",linewidth=1.1)+
      ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[1]]), color="red",linetype="dashed",linewidth=1.1)+
      ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[3]]),color="red",linewidth=1.1)+
      ggplot2::labs(x="Posterior Distribution of Beta Random Variable", y="Density", title = "Distribution of Samples",
                    subtitle=paste0("HDE: ",round(out$summary[[1]],2),"\nHDI: [",round(out$summary[[2]],2),", ",round(out$summary[[3]],2),"]"))+
      pcv_theme()
    
    if(nrow(plot_df[plot_df$sample == "Sample 2",]) == length(beta_support)){
      p <- p +
        ggplot2::geom_area(data=plot_df[plot_df$sample == "Sample 2",],fill="blue",alpha=0.5)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[5]]),color="blue",size=1.1)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[4]]),color="blue",linetype="dashed",size=1.1)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[6]]),color="blue",size=1.1)+
        ggplot2::labs(subtitle = paste0(
          "Sample 1:  ",round(out$summary[[1]],2)," [",round(out$summary[[2]],2),", ",round(out$summary[[3]],2),"]\n",
          "Sample 2:  ",round(out$summary[[4]],2)," [",round(out$summary[[5]],2),", ",round(out$summary[[6]],2),"]\n",
          "P[p1",dirSymbol,"p2] = ",post.prob.text))
      
      if(length(rope_int) == 2){
        p <- p + ggplot2::ggplot(rope_df, ggplot2::aes(x=X))+
          ggplot2::geom_histogram(bins=100,fill="purple",color="purple", alpha=0.7)+
          ggplot2::geom_histogram(data=data.frame("X"=rope_df[rope_df$X > rope_range[1] & rope_df$X < rope_range[2] & rope_df$X > hdi_diff[1] & rope_df$X < hdi_diff[2],]),
                                  bins=100,fill="gray10",color="gray10")+
          ggplot2::geom_segment(ggplot2::aes(x=rope_range[1],xend=rope_range[2],y=0,yend=0),linewidth=2,color="gray70")+
          ggplot2::geom_vline(ggplot2::aes(xintercept=hdi_diff[1]),linewidth=0.7)+
          ggplot2::geom_vline(ggplot2::aes(xintercept=hde_diff),linetype="dashed",linewidth=0.7)+
          ggplot2::geom_vline(ggplot2::aes(xintercept=hdi_diff[2]),linewidth=0.7)+
          ggplot2::labs(x="Posterior of Difference", y="Frequency", title = "Distribution of Difference",
                        subtitle =  paste0(
                          "Median Difference of ",round(hde_diff,2),"\n",100*rope_hdi,"% CI [",round(hdi_diff[1],2),", ",round(hdi_diff[2],2),"]\n",
                          rope_ci,"% HDI in [", rope_range[1],", ",rope_range[2], "]: ",round(rope_prob,2)))+
          pcv_theme()+
          ggplot2::theme(axis.title.y.left = ggplot2::element_blank(), axis.text.y.left = ggplot2::element_blank())+
          patchwork::plot_layout(widths = c(2,1))
      }
    }
    plot(p)
    out$plot<-p
  }
  return(out)
}

#' bayes_beta(dat=df[df$trait=="yii_hist_Fv.over.Fm" & df$id=="tp4",], col="treatment", compare = NULL,
#'            freq="value", bin="label",
#'            plot=T, rope_int=c(-0.05,0.05))
#'            
#' So at first I made this to use long data and with the intention of working mostly with the MV traits
#' Now I'm thinking that it should be good for long or wide data and SV or MV traits.
#' I could ask for arguments as things like longdata[condition, "value"] or wide$col
#' OR I could have some set of args to set like col="area.pixels", etc.
#' 
#' The first is more flexible if that's a concern which I think it is.
#' I also think that encourages better practices with priors.
#' If it can be easily plugged into an lapply then that would make sense to me.
#' something like lapply(fixCompare, function(i){thisFunction(subset1, subset2, etc)})
#' should be accessible/shown in the docs. BUT with long MV trait data that probably does get more complex.
#' For that example I'd have things like:
#' s1<-as.numeric(g1[rep(row.names(g1), g1[[freq]]) & g1$CONDITIONS, bin])
#' for wide MV traits...
#' hm. This might benefit from sv and mv methods for each distribution.
#' 
#' wide MV traits have a vector for each image which should be aggregated like what happens in pcvjoyplot.R
#' The thing is that the inputs are totally different to do that. The s1, s2 option works great for sv traits
#' but for wide mv traits I need:
#' columns to take
#' variable to group by
#' well it could take s1 and s2 as data.frames or matrices. 
#' 
#' 
#' 
#' I also think this has to have a specific function for each distribution that is supported,
#' then the main function can just call those methods.
#' 
#' The distributions that make sense to use include beta, normal, possibly a gamma/lognormal type thing?
#' lognormal is probably super easy to implement just as a logged normal given conjugacy
#' 
#' should I build in a way to do a "ROPE" like comparison when there is one sample?
#' basically just "prob that this posterior is in X interval"? Seems reasonable enough to me to include that.
#' 








#' ***********************************************************************************************
#' *************** `Gaussian Means (single value)` *********************************************
#' ***********************************************************************************************

#' @description
#' Internal function for Bayesian T Tests of gaussian data represented by single value traits.
#' @param s1 A vector of numerics drawn from a gaussian distribution.
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' .conj_gaussian_means_sv(s1=rnorm(100, 50,10), s2= rnorm(100, 60,12), priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'    plot=F, rope_range = c(-0.1, 0.1), rope_ci = 0.89, rope_hdi = 0.89,
#'    cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' @keywords internal
#' @noRd
.conj_gaussian_means_sv<-function(s1 = NULL, s2= NULL, priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
                            plot=F, rope_range = NULL, rope_ci = 0.89, rope_hdi = 0.89,
                            cred.int.level = 0.89, hypothesis="equal", support=NULL){
  out <- list()
  #* `Define support if it is missing`
  if(is.null(support)){
    support<-seq(floor(min(c(s1,s2))*0.8), ceiling(max(c(s1,s2))/0.8), length.out=10000)
  }
  #* `Get Mean, Variance, SE, and DF from s1`
  if(length(s1 > 1)){
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
    
    out$summary<-data.frame(HDE=hde1_mean, HDI_low = hdi1_mean[1], HDI_high = hdi1_mean[2])
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
    if(length(s2 > 1)){
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
      
      out$summary<-data.frame(HDE=hde2_mean, HDI_low = hdi2_mean[1], HDI_high = hdi2_mean[2])
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
      post.prob.out <-  .post.prob.from.pdfs(pdf1, pdf2, hypothesis)
      post.prob <-post.prob.out$post.prob
      dirSymbol <-post.prob.out$direction
    }else{
      post.prob <- 1
      dirSymbol = "="
    }
    if(post.prob<1e-5){post.prob.text = "<1e-5"} else { post.prob.text<-round(post.prob,5)}
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
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_hdi ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi1[1], HDI_rope_high = hdi1[2], rope_prob = rope_prob)
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
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number between 0.0001 and 0.9999 representing the "bin".
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' makeMvGauss<-function(bins=180,mu,sigma){
#'    setNames(data.frame(matrix(hist(rnorm(2000,mu, sigma), breaks=seq(1,bins,1), plot=F)$counts, nrow=1)),paste0("b",1:(bins-1) ))
#'    }
#' mv_gauss<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=50, sigma=10 )})),
#'                 do.call(rbind, lapply(1:30, function(i){makeMvGauss(bins=180, mu=60, sigma=12 )})))
#' .conj_gaussian_means_mv(s1 = mv_gauss[1:30,], s2= mv_gauss[31:60,], priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
#'                         plot=F, rope_range = NULL, rope_ci = 0.89, rope_hdi = 0.89,
#'                         cred.int.level = 0.89, hypothesis="equal", support=NULL)
#' @keywords internal
#' @noRd

.conj_gaussian_means_mv<-function(s1 = NULL, s2= NULL, priors = list( mu=c(0,0),n=c(1,1),s2=c(20,20) ),
                            plot=F, rope_range = NULL, rope_ci = 0.89, rope_hdi = 0.89,
                            cred.int.level = 0.89, hypothesis="equal", support=NULL){
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
  bins_order<-sort(histCols_bin, index.return=T)$ix
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
  
  out$summary<-data.frame(HDE=hde1_mean, HDI_low = hdi1_mean[1], HDI_high = hdi1_mean[2])
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
    bins_order<-sort(histCols_bin, index.return=T)$ix
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
    
    out$summary<-data.frame(HDE=hde2_mean, HDI_low = hdi2_mean[1], HDI_high = hdi2_mean[2])
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
      post.prob.out <-  .post.prob.from.pdfs(pdf1, pdf2, hypothesis)
      post.prob <-post.prob.out$post.prob
      dirSymbol <-post.prob.out$direction
    }else{
      post.prob <- 1
    }
    if(post.prob<1e-5){post.prob.text = "<1e-5"} else { post.prob.text<-round(post.prob,5)}
    out$summary$hyp = hypothesis
    out$summary$post.prob = post.prob
    out$dirSymbol = dirSymbol
    
    #* `ROPE Comparison`
    #* might move this out of the s2 dependent area and add the logic to where posterior is defined.
    #* maybe rope_posterior is an if(!is.null(s2)){current method} else {10000 draws from pdf1}
    #* easier to do that now than later I guess.
    
  }
  
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      post1 = extraDistr::rlst(10000,v1_n, m1_n, se1)
      if(!is.null(s2)){
      post2 = extraDistr::rlst(10000,v2_n, m2_n, se2)
      posterior = post1 - post2
      } else {
        posterior = post1
      }
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_hdi ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi1[1], HDI_rope_high = hdi1[2], rope_prob = rope_prob)
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
#'   setNames(data.frame(matrix(hist(rbeta(2000,a,b), breaks=seq(0,1,length.out=n), plot=F)$counts, nrow=1)),paste0("b0.",1:(n-1)))
#' }
#' 
#' mv_beta<-rbind(do.call(rbind, lapply(1:30, function(i){makeMvBeta(n=100, a=5, b=8 )})),
#'                do.call(rbind, lapply(1:30, function(i){makeMvBeta(n=100, a=10, b=3 )})))
#' 
#' .conj_beta_mv(s1 = mv_beta[1:30,], s2= mv_beta[31:60,], priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
#'                         plot=F, rope_range = c(-0.1, 0.1), rope_ci = 0.89, rope_hdi = 0.89,
#'                         cred.int.level = 0.89, hypothesis="equal")
#' @keywords internal
#' @noRd
.conj_beta_mv<-function(s1 = NULL, s2= NULL, priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
                        plot=F, rope_range = NULL, rope_ci = 0.89, rope_hdi = 0.89,
                        cred.int.level = 0.89, hypothesis="equal"){
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
  bins_order<-sort(histCols_bin, index.return=T)$ix
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
  out$summary <- data.frame(HDE=hde1, HDI_low = hdi1[1], HDI_high = hdi1[2])
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
    bins_order<-sort(histCols_bin, index.return=T)$ix
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
    if(post.prob<1e-5){post.prob.text = "<1e-5"} else { post.prob.text<-round(post.prob,5)}
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
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_hdi ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi1[1], HDI_rope_high = hdi1[2], rope_prob = rope_prob)
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
#'                         plot=F, rope_range = c(-0.1, 0.1), rope_ci = 0.89, rope_hdi = 0.89,
#'                         cred.int.level = 0.89, hypothesis="equal")
#' @keywords internal
#' @noRd
.conj_beta_sv<-function(s1 = NULL, s2= NULL, priors = list(a=c(0.5,0.5),b=c(0.5,0.5)),
                       plot=F, rope_range = NULL, rope_ci = 0.89, rope_hdi = 0.89,
                       cred.int.level = 0.89, hypothesis="equal"){

  if(any(c(s1, s2)>1)){stop("Values above 1 cannot be used with the beta distribution")}
  
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
  out$summary <- data.frame(HDE=hde1, HDI_low = hdi1[1], HDI_high = hdi1[2])
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
    if(post.prob<1e-5){post.prob.text = "<1e-5"} else { post.prob.text<-round(post.prob,5)}
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
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_hdi ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi1[1], HDI_rope_high = hdi1[2], rope_prob = rope_prob)
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
#' *************** `General Functions` *******************************************************
#' ***********************************************************************************************


#' @keywords internal
#' @noRd
.post.prob.from.pdfs<-function(pdf1, pdf2, hypothesis){
  if(hypothesis == "unequal"){
    post.prob <- 1-sum(apply(cbind(pdf1, pdf2),MARGIN=1,function(i) min(i)),na.rm=T)
    dirSymbol="!="
  } else if(hypothesis == "equal"){
    post.prob <- sum(apply(cbind(pdf1, pdf2),MARGIN=1,function(i) min(i)),na.rm=T)
    dirSymbol="="
  }else if(hypothesis == "lesser"){
    direction <- pdf1 <= pdf2
    post.prob <- sum(pdf1 * direction,na.rm=T)
    dirSymbol="<"
  }else if(hypothesis == "greater"){
    direction <- pdf1 >= pdf2
    post.prob <- sum(pdf1 * direction,na.rm=T)
    dirSymbol=">"
  }else{
    stop("hypothesis must be either unequal, equal, lesser, or greater")
  }
  return(list('post.prob'=post.prob, "direction" = dirSymbol))
}



