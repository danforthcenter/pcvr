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
#' For brevity the HDE and HDI of the samples are not returned when
#' using the dirichlet method. The "dirichlet" method does not compute an HDI and only compares how much of the posterior
#' distribution (sample 1 - sample 2) is within the ROPE interval. The "dirichlet2" method does make an HDI and run normal 
#' ROPE tests per each bin using the alpha vector from each sample to simulate draws from the dirichlet. Bear in mind
#'  that dirichlet is sensitive to the total number of counts in the data. If the mean vector
#' (alpha_vec ~ mean_vec * precision) were to be used instead of the alpha vector then the distributions
#' in each bin would be much wider.
#' 
#' @details
#' \itemize{
#' \item{\strong{"dirichlet": } \code{priors = list(alpha = rep(1, ncol(s1), prec_a = 0.5, prec_b = 0.5 ))}, 
#'    where alpha is a vector of counts and prec_a/b describe the gamma distribution of precision.
#'    }
#' }
#' 
#' @param s1 A data.frame or matrix of multi value traits. The column names should include a number representing the "bin".
#' @param s2 An optional second sample of the same form as s2.
#' @examples
#' if(FALSE){
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
#' out <- .conj_diri_mv_2(s1, s2, priors=NULL, plot=TRUE, rope_range = c(-0.1, 0.1),
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
#'                        
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
#' diri_ex_1 <- conjugate(s1, s2, method = "dirichlet", priors=NULL, plot=FALSE,
#'       rope_range = c(-0.025, 0.025), rope_ci = 0.89, 
#'       cred.int.level = 0.89, hypothesis="equal")
#'       
#' diri_ex_2 <- conjugate(s1, s2, method = "dirichlet2", priors=NULL, plot=FALSE,
#'       rope_range = c(-0.025, 0.025), rope_ci = 0.89, 
#'       cred.int.level = 0.89, hypothesis="equal")  
#' }
#' }
#' 
#' @keywords internal
#' @noRd

.conj_diri_mv_2<-function(s1 = NULL, priors=NULL,
                          plot=FALSE){
  out <-list()
  
  if(is.null(priors)){
    priors <- list(alpha = rep(1, ncol(s1)), # alpha vector priors
                   prec_a = 0.5, # gamma priors on precision
                   prec_b = 0.5)
  }
  
  if(!is.null(s1)){
    alpha1_prime <- priors$alpha + colSums(s1) # updating prior with sum of samples in each bin
    
    precision1 <- sum(alpha1_prime) # reparameterize alpha1 to mean and precision
    mean1 <- alpha1_prime / precision1
    
    prec1_prime_a <- priors$prec_a + sum(colSums(s1)) # A updates as A' = A + sum(obs)
    prec1_prime_b <- priors$prec_b + nrow(s1) # B updates as B' = B + n(obs)
    
    out$summary <- data.frame('HDE_1' = NA, 'HDI_1_low' = NA, 'HDI_1_high' = NA)
    
    out$pdf = mean1
    out$posterior$precision1 = precision1
    out$posterior$alpha = alpha1_prime
    out$posterior$prec_a <- prec1_prime_a
    out$posterior$prec_b <- prec1_prime_b
    
    out$posteriorDraws = extraDistr::rdirichlet(10000, mean1*ncol(s1))
    #out$pdf <- pdf1
    
    if(plot){
      out$plot_df <- data.frame("bin" = 1:ncol(s1), "prob"=mean1, "sample"=rep("Sample 1",ncol(s1) ))
    }
    
  } else{stop("s1 is required")}
  
  return(out)
}



#' @description
#' Function for dirichlet (option 1) ROPE comparison
#' 
#' @keywords internal
#' @noRd

.diri2_rope <- function(sample_results, rope_range = c(-0.1, 0.1), rope_ci = 0.89, plot){


  if(length(sample_results)==2){
    posterior <- sample_results[[1]]$posteriorDraws - sample_results[[2]]$posteriorDraws
    direction <- sample_results[[1]]$posteriorDraws >= sample_results[[2]]$posteriorDraws
  } else{
    posterior = sample_results[[1]]$posteriorDraws
    direction = rep(TRUE, length(posterior))
  }
  
  abs_mean_diff <- abs(posterior)
  norm_abs_mean_diff <- abs_mean_diff/sum(abs_mean_diff)
  norm_mean_diff <- ifelse(direction, 1, -1) * norm_abs_mean_diff
  
  rope_probs<-lapply(1:ncol(posterior), function(bin){
    as.numeric(bayestestR::rope(posterior[,bin], range = rope_range, ci_method = "HDI", ci=rope_ci))
  })
  
  mean_rope_prob = sum( unlist(rope_probs)* norm_abs_mean_diff )
  
  HDI_rope <- as.data.frame(bayestestR::hdi(as.data.frame(posterior), ci = rope_ci))
  HDE_rope_pre <- as.data.frame(bayestestR::hdi(as.data.frame(posterior), ci = 0.001))
  HDE_rope <- data.frame("Parameter" = HDE_rope_pre[["Parameter"]],
                         "HDE" = rowMeans(HDE_rope_pre[,c("CI_low", "CI_high")]))
  rope_res <- list()
  rope_res$summary <- HDE_rope
  
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
  
}

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
#' @examples
#' if(FALSE){
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
#' }
#' 
#' @keywords internal
#' @noRd

.conj_diri_mv_1<-function(s1 = NULL, priors=NULL, plot=FALSE){
  out <-list()
  
  if(is.null(priors)){
    priors <- list(alpha = rep(1, ncol(s1)), # alpha vector priors
                   prec_a = 0.5, # gamma priors on precision
                   prec_b = 0.5)
  }
  
  if(!is.null(s1)){
    alpha1_prime <- priors$alpha + colSums(s1) # updating prior with sum of samples in each bin
    
    precision1 <- sum(alpha1_prime) # reparameterize alpha1 to mean and precision
    mean1 <- alpha1_prime / precision1
    
    prec1_prime_a <- priors$prec_a + sum(colSums(s1)) # A updates as A' = A + sum(obs)
    prec1_prime_b <- priors$prec_b + nrow(s1) # B updates as B' = B + n(obs)
    
    out$summary <- data.frame('HDE_1' = NA, 'HDI_1_low' = NA, 'HDI_1_high' = NA)
    
    out$pdf = mean1
    out$posterior$precision1 = precision1
    out$posterior$alpha = alpha1_prime
    out$posterior$prec_a <- prec1_prime_a
    out$posterior$prec_b <- prec1_prime_b
    
    if(plot){
      out$plot_df <- data.frame("bin" = 1:ncol(s1), "prob"=mean1, "sample"=rep("Sample 1",ncol(s1) ))
    }
    
  } else{stop("s1 is required")}
  
  return(out)
}

#' @description
#' Function for dirichlet (option 1) ROPE comparison
#' 
#' @keywords internal
#' @noRd

.diri1_rope <- function(sample_results, rope_range = c(-0.1, 0.1), rope_ci = 0.89, plot){
  rope_res <- list()
  
  if(!is.null(rope_range)){
    if(length(rope_range) == 2){
      
      #* `difference of mean vectors`
      
      post1 = sample_results[[1]]$pdf
      if(length(sample_results)==2){
        post2 = sample_results[[2]]$pdf
        posterior = post1 - post2
        direction <- post1 >= post2
      } else {
        posterior = post1
        direction <- rep(TRUE, length(post1))
      }
      
      abs_mean_diff <- abs(posterior)
      norm_abs_mean_diff <- abs_mean_diff/sum(abs_mean_diff)
      norm_mean_diff <- ifelse(direction, 1, -1) * norm_abs_mean_diff
      
      rope.prob <- sum(apply(cbind(norm_abs_mean_diff, max(rope_range)), MARGIN=1, function(i) min(i)), na.rm=T)
      
      if(plot){
        rope_res$rope_df <- data.frame(bin = 1:length(norm_mean_diff), CI = "1")
        rope_res$rope_df$CI_low <- unlist(lapply(1:length(norm_mean_diff), function(i){min(norm_mean_diff[[i]], 0)}))
        rope_res$rope_df$CI_high <- unlist(lapply(1:length(norm_mean_diff), function(i){max(norm_mean_diff[[i]], 0)}))
      }
      
      rope_res$summary = data.frame('HDE_rope' = NA, "HDI_rope" = NA,"rope_probs" = NA, "mean_rope_prob" = rope.prob)
      #rope_res$summary <- data.table::as.data.table(rope_res$summary)
      #out$summary <- as.data.frame(out$summary)
    } else{
      stop("rope must be a vector of length 2")
    }
  }
  return(rope_res)
}



#' @description
#' Plotting function for dirichlet conjugate methods
#' 
#' @keywords internal
#' @noRd

.conj_diri_plot <- function(res, s2, rope_range, rope_ci){
  # If I make two dirichlet options then having two modes for this function would be good.
  
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
    
    xLims <- ggplot2::layer_scales(p)$x$range$range
    
    rect_width <- min(diff(unique(res$rope_df$bin)))/2.5
    
    #* this makes the rope probability look incorrect because it does not accurately show the distribution,
    #* it just makes a bar for the HDI, but density within the HDI varies widely for my examples.
    
    rdf <- res$rope_df
    cis<-rev(unique(rdf$CI))
    virPal <- viridis::plasma(length(cis), direction=-1)
    mode = "2"
    if(length(virPal)==1){
      virPal = "gray40"
      mode = "1"
    }
    
    rope_plot <- ggplot2::ggplot(rdf, ggplot2::aes(xmin=.data$bin-rect_width, xmax = .data$bin+rect_width,
                                                   ymin=.data$CI_low, ymax=.data$CI_high ))+
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
    
    if(mode == "1"){
      rope_plot <- rope_plot +
        ggplot2::theme(legend.position = "none")
    }
    
    if(!is.null(s2)){
      
      if(mode == "1"){
        p <- p +
          ggplot2::labs(subtitle = paste0("P[p1",res$dirSymbol,"p2] = ",post.prob.text, "\n",
                                          "Posterior in [", rope_range[1],", ",rope_range[2], "]: ",
                                          round(res$summary$mean_rope_prob, 2) ))
      } else{
        p <- p +
          ggplot2::labs(subtitle = paste0("P[p1",res$dirSymbol,"p2] = ",post.prob.text, "\n",
                                          rope_ci,"% HDI in [", rope_range[1],", ",rope_range[2], "]: ",
                                          round(res$summary$mean_rope_prob, 2) ))
      }
      res<-res[-which(names(res)=="dirSymbol")]
    }
    
    layout<-c(patchwork::area(1,1,3,6),
              patchwork::area(4,1,4,6))
    
    p <- p / rope_plot + patchwork::plot_layout(design = layout)
  }
  plot(p)
  
  res$plot<-p
  return(res)
}


