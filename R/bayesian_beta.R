
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
#' 
#' 
#' I also think this has to have a specific function for each distribution that is supported,
#' then the main function can just call those methods.

conj_beta_sv<-function(){}
#* `arguments`
s1 = sv[sv$DAS==19 & sv$fertilizer==100 & sv$genotype=="B73", "area_cm2"]; s1<-s1/sum(s1) # just for beta purposes
s2 = sv[sv$DAS==19 & sv$fertilizer==100 & sv$genotype=="Mo17", "area_cm2"]; s2<-s2/sum(s2) # just for beta purposes
priors = list(a=c(0.5,0.5),b=c(0.5,0.5)) # jeffrey priors
plot=T
rope_range=c(-0.1,0.1)
rope_ci = 0.89
rope_hdi = 0.89
cred.int.level = 0.89
hypothesis="equal"

#* In function

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
  plot_df <- data.frame("range"=beta_support,
                                      "prob"=pdf1,
                                      "sample"=rep("Sample 1",length(beta_support)
  ))
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
  }else{
    post.prob <- 1
  }
  if(post.prob<1e-5){post.prob = "<1e-5"}
  out$summary$hyp = hypothesis
  out$summary$post.prob = post.prob
  
  #* `Use ROPE method to test if the difference in distributions is meaningful`
  if(!is.logical(rope_range)){
    if(length(rope_range) == 2){
      post1 = rbeta(10000,a1_prime,b1_prime)
      post2 = rbeta(10000,a2_prime,b2_prime)
      posterior = post1 - post2
      hdi_diff <- as.numeric(bayestestR::hdi(posterior, ci = rope_hdi ))[2:3]
      hde_diff <- median(posterior)
      rope_prob <- as.numeric(bayestestR::rope(posterior, range = rope_range, ci_method = "HDI", ci=rope_ci))
      rope_test <-  data.frame(HDE_rope=hde_diff, HDI_rope_low = hdi1[1], HDI_rope_high = hdi1[2], rope_prob = rope_prob)
      out$summary <- cbind(out$summary,rope_test)
      if(plot){
        rope_df <- data.frame("X"=posterior)
      }
    }else{
      stop("rope must be a vector of length 2")
    }
  }
  #* `Keep data from s2 for plotting`
  if(plot){
    plot_df <- rbind(plot_df,data.frame("range"=beta_support,
                              "prob"=pdf2,
                              "sample"=rep("Sample 2",length(beta_support))
    ))
  }
}


if(plot){
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x=range, y=prob))+
    ggplot2::geom_area(data=plot_df[plot_df$sample == "Sample 1",],fill="red",alpha=0.5)+
    ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[2]] ),color="red",size=1.1)+
    ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[1]]), color="red",linetype="dashed",size=1.1)+
    ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[3]]),color="red",size=1.1)+
    ggplot2::labs(x="Posterior Distribution of Beta Random Variable", y="Density", title = "Distribution of Samples",
                  subtitle=paste0("HDE: ",round(out$summary[[1]],2),"\nHDI: [",round(out$summary[[2]],2),", ",round(out$summary[[3]],2),"]"))
  
  if(nrow(plot_df[plot_df$sample == "Sample 2",]) == length(beta_support)){
    p <- p +
      ggplot2::geom_area(data=plot_df[plot_df$sample == "Sample 2",],fill="blue",alpha=0.5)+
      ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[5]]),color="blue",size=1.1)+
      ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[4]]),color="blue",linetype="dashed",size=1.1)+
      ggplot2::geom_vline(ggplot2::aes(xintercept=out$summary[[6]]),color="blue",size=1.1)+
      ggplot2::labs(subtitle = paste0(
        "Sample 1:  ",round(out$summary[[1]],2)," [",round(out$summary[[2]],2),", ",round(out$summary[[3]],2),"]\n",
                     "Sample 2:  ",round(out$summary[[4]],2)," [",round(out$summary[[5]],2),", ",round(out$summary[[6]],2),"]\n",
                     "P[p1",dirSymbol,"p2] = ",round(post.prob,5)))
    
    
    
    if(length(rope_int) == 2){
      p <- p + ggplot2::ggplot(rope_df, ggplot2::aes(x=X))+
        ggplot2::geom_histogram(bins=100,fill="#d4b8e6",color="#d4b8e6")+
        ggplot2::geom_histogram(data=data.frame("X"=rope_df[rope_df$X > rope_range[1] & rope_df$X < rope_range[2] & rope_df$X > hdi_diff[1] & rope_df$X < hdi_diff[2],]),
                                bins=100,fill="gray10",color="gray10")+
        ggplot2::geom_segment(ggplot2::aes(x=rope_range[1],xend=rope_range[2],y=-10,yend=-10),size=2,color="gray70")+
        ggplot2::geom_vline(ggplot2::aes(xintercept=hdi_diff[1]),size=1.1)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=hde_diff),linetype="dashed",size=1.1)+
        ggplot2::geom_vline(ggplot2::aes(xintercept=hdi_diff[2]),size=1.1)+
        ggplot2::labs(x="Posterior of Difference", y="Frequency", title = "Distribution of Difference",
                      subtitle =  paste0(
                        "Difference: ",round(hde_diff,2)," [",round(hdi_diff[1],2),", ",round(hdi_diff[2],2),"]\n",
                        "ROPE Interval: [",rope_range[1],", ",rope_range[2],"]\n",
                        rope_ci,"% HDI inside ROPE:  ",round(rope_prob,2)))+
        plot_layout(widths = c(2,1))
    }
  }
  plot(p)
  out$plot<-p
}

return(out)







