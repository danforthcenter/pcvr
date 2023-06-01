#' Plot Histograms for multi value trait plantCV data and optionally compare distributions
#' 
#' @param df Data frame to use. Long format with multi value traits is expected.
#' @param index A multi value trait as a character string. Must be present in `trait`.
#' @param group A length 1 or 2 character vector. This is used for faceting the joyplot and identifying groups for testing.
#' @param method A method to use in comparing distributions/means. Currently "beta", "gaussian", and "ks" are supported. See details for explanations of tests.
#' @param compare Groups to compare. By default this is set to FALSE, which corresponds to no testing. Other values of compare are passed to fixCompare to make t.test comparisons using ggpubr. In short, NULL will run all pairwise T tests, a single value of the X axis variable will compare that level to all other levels of the X variable, alternatively this can be a list as used by ggpubr: list(c("level1", "level2"), c("level1", "level3"))
#' @param priors Parameters for prior distributions if using method = "beta" or "gaussian". If left NULL then wide weak priors are used. This can be set as a single set of parameters or with unique values for each group in `compare`, but it is advised to use the default priors.
#' @param hyp Hypothesis to test for beta and gaussian methods, one of "unequal", "greater", "lesser". Defaults to "unequal" if left NULL.
#' @param bin Column containing histogram (multi value trait) bins. Defaults to "label".
#' @param freq Column containing histogram counts. Defaults to "value"
#' @param trait Column containing phenotype names. Defaults to "trait".
#' @keywords bayesian, ggplot, multi value trait, pcv.joyplot
#' @import ggplot2
#' @import extraDistr
#' 
#' @details 
#' This function is similar to pcv.joyplot but uses different plotting aesthetics.
#' 
#' The method argument is used for statistical testing of groups. There are currently four options, "beta", "gaussian", "ks", and "emd".
#' \itemize{
#'  \item{"beta"}{The beta method is for use with bounded data (0,1) and creates a posterior beta distribution using the histogram counts and bins and priors (which default to beta(0.5,0.5)). The posterior distributions are compared across the groups in `compare` and a posterior probability is returned for each which corresponds to P[hyp]. }
#'  \item{"gaussian"}{The gaussian method is for use with multi value traits that are not bounded 0-1 and are roughly gaussian. Like the beta method this assumes weak priors and uses the histogram data to make a posterior distribution (here a T distribution). Here there will be two posterior probabilities returned for each comparison, one compares the distribution of means only and is roughly analogous to a T test, the other compares the entire posterior distributions and is more similar to a Z test. The posterior probability still takes the form P[hyp]. Note that currently these consider pixels as reps, so it is easy to have very small biological changes at very high posterior probability.}
#'  \item{"ks"}{The ks method performs a ks test on the PDF of each histogram over a set support range. This returns a P value corresponding to a standard KS test that distributions are the same. For this method `hyp` is ignored.}
#'  \item{"emd}{The emd method returns earth mover's distance for compared histograms. This is done via the emdist package and can be very slow. This is included because it is canonical but in the case where EMD is desired it may be worth implementing it directly.}
#' }
#' 
#' @return Returns either a ggplot object or a list containing a ggplot and a dataframe of statistical comparisons (if compare is not FALSE).
#' 
#' @examples 
#' df <- read.pcv("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest2.csv", "long", F)
#' pcv.hists(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"))
#' pcv.hists(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"), method="ks")
#' pcv.hists(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"), method="beta")
#' pcv.hists(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"), method="gaussian")
#' @export


pcv.hists<-function(df = NULL, index = NULL, group = NULL,
                    method=NULL,
                    compare= NULL, priors=NULL,hyp=NULL, # c("unequal", "greater", "lesser")
                    bin="label", freq="value", trait="trait"){
  #* ***** `troubleshooting test values`
  # df=df; index = "index_frequencies_index_ndvi"; group=c("timepoint", "genotype"); method=NULL
  # compare= NULL; priors=NULL; bin="label"; freq="value"; trait="trait"; hyp=NULL
  
  #* ***** `general calculated values`
  if(is.null(index)){sub<-df
  }else{ sub<-df[df[[trait]]==index, ] }
  if(length(unique(sub[[trait]]))>1){warning("More than one trait found, consider an `index` argument")}
  sub$bin = as.numeric(sub[[bin]])
  sub$freq = as.numeric(sub[[freq]])
  
  if(is.null(group)){group = "dummy"; df$dummy = "dummy"; sub$dummy="dummy"}
  if(length(group)==1){sub$fill = sub[[group]]; sub$y = sub[[group]]; facet_layer = ggplot2::facet_wrap(paste0("~",group)) }
  if(length(group)==2){sub$fill = sub[[group[1] ]]; sub$y = sub[[group[1] ]]; facet_layer=ggplot2::facet_grid(as.formula(paste0("~",group[2])))}
  
  sub$grouping<-interaction(sub[,c(group)], drop=T)
  
  # default compare to NULL, but if F then skip all testing 
  if(is.logical(compare) && compare==F){
    doStats=F
  } else if(is.null(compare)){compareTests<-fixCompare(compare,sub,"grouping", T) ; doStats=T
  }else{compareTests=fixCompare(compare,sub,"grouping"); doStats=T}
  if(is.null(hyp)){hyp<-"unequal"}
  
  #* ***** `default joyplot`
  if(is.null(method) | match.arg(method, choices = c("beta", "gaussian", "ks", "mixture", "emd"))=="emd"){
    
    datsp=split(sub, sub$grouping, drop=T)
    bw<-min(diff(sort(as.numeric(unique(sub[[bin]])))))*0.75
    distParams<-lapply(datsp, function(D){
      X1 <- as.numeric(D[rep(rownames(D), D[[freq]]), bin])
      dens<-density(X1, from = min(sub$bin,na.rm=T), to = max(sub$bin,na.rm=T), n = 2^10, bw=bw) 
      return(dens)
    })
    names(distParams)<-names(datsp)
    dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
      dens<-distParams[[nm]]
      out<-data.frame(xdens= dens$x, ydens=dens$y)
      out[,(ncol(out)+1):(ncol(out)+length(group))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
      colnames(out)<-c("xdens", "ydens", group)
      out$y<-out[[group[1] ]]
      out
    }))
    p<-ggplot2::ggplot(data.frame(x=seq(0,1,0.1), y = 0), ggplot2::aes(x=.data$x,y=.data$y))+
      facet_layer+
      lapply(1:length(distParams), function(i){
        dens<-distParams[[i]]
        lapDat<-data.frame(x=dens$x, y = dens$y)
        lapDat[[group[1] ]]=strsplit(names(distParams)[i],"[.]")[[1]][1]
        if(length(group)==2){lapDat[[group[2] ]] = strsplit(names(distParams)[i],"[.]")[[1]][2]}
        lapDat$fill <-lapDat[[group[1] ]]
        ggplot2::geom_area(data=lapDat, ggplot2::aes(x=.data$x, y=.data$y, fill=.data$fill))
      })+
      ggplot2::scale_x_continuous(n.breaks=5, labels = ~round(.,1))+
      ggplot2::scale_fill_viridis_d(option="viridis")+
      ggplot2::labs(x=index, y=freq)
    
    if(match.arg(method, choices = c("beta", "gaussian", "ks", "mixture", "emd"))=="emd" & doStats){
      #* use emd1d/pcv.emd to make an EMD distance matrix/df
      #* this should probably return a long dataframe per normal AND a distance matrix
      
      outStats<-do.call(rbind, lapply(compareTests, function(comp){
        g1<-as.character(comp[1])
        g2<-as.character(comp[2])
        d1<-sub[sub$grouping == g1, ]
        d2<-sub[sub$grouping == g2, ]
        d1s<-aggregate(as.formula(paste0(freq,"~",bin)), d1, sum)[[freq]]
        d2s<-aggregate(as.formula(paste0(freq,"~",bin)), d2, sum)[[freq]]
        emd_res <- emd1d(d1s, d2s)
        data.frame(group1 = g1, group2 = g2, emd = emd_res, method="emd")
        
      }))
      
    }
    
  } else if(match.arg(method, choices = c("beta", "gaussian", "ks", "mixture", "emd"))=="ks"){
    #* ***** `Non parametric joyplot`
    datsp=split(sub, sub$grouping, drop=T)
    bw<-min(diff(sort(as.numeric(unique(sub[[bin]])))))*0.75 # calculating a set bandwidth for all groups
    distParams<-lapply(datsp, function(D){
      X1 <- as.numeric(D[rep(rownames(D), D[[freq]]), bin])
      dens<-density(X1, from = min(sub$bin,na.rm=T), to = max(sub$bin,na.rm=T), n = 2^10, bw=bw)
      return(dens)
    })
    names(distParams)<-names(datsp)
    dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
      dens<-distParams[[nm]]
      out<-data.frame(xdens= dens$x, ydens=dens$y)
      out[,(ncol(out)+1):(ncol(out)+length(group))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
      colnames(out)<-c("xdens", "ydens", group)
      out$y<-out[[group[1] ]]
      out
    }))
    p<-ggplot2::ggplot(data.frame(x=seq(0,1,0.1), y = 0), ggplot2::aes(x=.data$x,y=.data$y))+
      facet_layer+
      lapply(1:length(distParams), function(i){
        dens<-distParams[[i]]
        lapDat<-data.frame(x=dens$x, y = dens$y)
        lapDat[[group[1] ]]=strsplit(names(distParams)[i],"[.]")[[1]][1]
        if(length(group)==2){lapDat[[group[2] ]] = strsplit(names(distParams)[i],"[.]")[[1]][2]}
        lapDat$fill <-lapDat[[group[1] ]]
        ggplot2::geom_area(data=lapDat, ggplot2::aes(x=.data$x, y=.data$y, fill=.data$fill))
      })+
      ggplot2::scale_x_continuous(n.breaks=5, labels = ~round(.,1))+
      ggplot2::scale_fill_viridis_d(option="viridis")+
      ggplot2::labs(x=index, y=freq)
    
    #* ***** `Non Parametric Stats`
    if(doStats){
      outStats<-do.call(rbind, lapply(compareTests, function(comp){
        g1<-as.character(comp[1])
        g2<-as.character(comp[2])
        ks<-suppressWarnings(ks.test(distParams[[g1]]$y, distParams[[g2]]$y))
        ret<-data.frame(group1 = g1, group2 = g2, p = ks$p.value, method="ks", null = "same distribution")
        return(ret)
      }))
    }
  } else if(match.arg(method, choices = c("beta", "gaussian", "ks", "emd"))=="beta"){
    #* ***** `Beta distribution joyplot`
    datsp=split(sub, sub$grouping, drop=T) # split data into panel groups
    
    if(is.null(priors)){ priors = list(a = rep(0.5, times=length(datsp)), b = rep(0.5, times=length(datsp))) } # assume weak prior on everything
    
    distParams<-lapply(1:length(datsp), function(i){ # get beta parameters from histogram data per panel
      D=datsp[[i]]
      X1 <- as.numeric(D[rep(rownames(D), D[[freq]]), bin])
      n1<- length(X1)
      xbar1 = (1/n1) * sum(X1)
      variance1 = (1/(n1))*sum( (X1 - xbar1)^2 )
      alpha1 = (xbar1 * ( (xbar1 * (1-xbar1) )/(variance1)-1 )) + priors$a[i]
      beta1 = ((1-xbar1) * ( (xbar1 * (1-xbar1) )/(variance1)-1 )) + priors$b[i]
      return(c(alpha1, beta1))
    })
    names(distParams)<-names(datsp)
    
    p<-ggplot2::ggplot(data.frame(x=seq(0,1,0.1), y = 0), ggplot2::aes(x=.data$x,y=.data$y))+
      facet_layer+
      lapply(1:length(distParams), function(i){
        pars<-distParams[[i]]
        lapDat<-data.frame(x=seq(0,1,0.1), y = 0)
        lapDat[[group[1] ]]=strsplit(names(distParams)[i],"[.]")[[1]][1]
        if(length(group)==2){lapDat[[group[2] ]] = strsplit(names(distParams)[i],"[.]")[[1]][2]}
        lapDat$fill <-lapDat[[group[1] ]]
        ggplot2::stat_function(data = lapDat, fun=dbeta, args=list( shape1=pars[1], shape2=pars[2]), ggplot2::aes(fill = .data$fill),
                      geom = "polygon", alpha = 0.85)
      })+
      ggplot2::scale_x_continuous(n.breaks=5, labels = ~round(.,1))+
      ggplot2::scale_fill_viridis_d(option="viridis")+
      ggplot2::labs(x=index, y=freq)
    
    
    support<-seq(0.001, 0.999, 0.001)
    
    dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
      pars<-distParams[[nm]]
      dens<-dbeta(support, pars[1], pars[1])
      #* define support, get density across support given these parameters, make dens_df
      out<-data.frame(xdens= support, ydens=dens)
      out[,(ncol(out)+1):(ncol(out)+length(group))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
      colnames(out)<-c("xdens", "ydens", group)
      out$y<-out[[group[1] ]]
      out
    }))
    
    #* ***** `Compare generating beta distributions`
    if(doStats){
      outStats<-do.call(rbind, lapply(compareTests, function(comp){
        #comp=compareTests[[i]]
        g1<-as.character(comp[1])
        g2<-as.character(comp[2])
        
        alpha1<-distParams[[g1]][1]
        beta1<-distParams[[g1]][2]
        alpha2<-distParams[[g2]][1]
        beta2<-distParams[[g2]][2]
        
        my_dense1 <- dbeta(support, alpha1, beta1)
        my_pdf1 <- my_dense1/sum(my_dense1)
        my_dense2 <- dbeta(support, alpha2, beta2)
        my_pdf2 <- my_dense2/sum(my_dense2)
        if(hyp=="unequal"){
          post.prob = 1-sum(apply(cbind(my_pdf1,my_pdf2), MARGIN=1,function(i) min(i)),na.rm=T) # P[pdf1 != pdf2]
        } else if (hyp == "lesser"){
          direction <- my_pdf1 <= my_pdf2
          post.prob <- sum(my_pdf1*direction,na.rm=T)
        } else if(hyp == "greater"){
          direction <- my_pdf1 >= my_pdf2
          post.prob <- sum(my_pdf1*direction,na.rm=T)
        }
        ret<-data.frame(group1 = g1, group2 = g2, post.prob = post.prob, method="beta", null = hyp)
        ret
      }))
    }
    
  } else if(match.arg(method, choices = c("beta", "gaussian", "ks"))=="gaussian"){
    #* ***** `Gaussian distribution joyplot`
    datsp=split(sub, sub$grouping, drop=T)
    if(is.null(priors)){ priors <- list( m=rep(0,length(datsp)), n=rep(1,length(datsp)), s2=rep(20,length(datsp)) ) } # mean, number, and variance
    distParams<-lapply(1:length(datsp), function(i){ 
      D=datsp[[i]]
      X1 <- as.numeric(D[rep(rownames(D), D[[freq]]), bin])
      obs_n = length(X1)
      obs_m = mean(X1)
      obs_s2 = var(X1)
      priorDF = priors$n[i]-1
      n1_n = priors$n[i] + obs_n
      m1_n = (obs_n*obs_m + priors$n[i]*priors$m[1])/n1_n
      newDF = priorDF + obs_n
      s2_n = ((obs_n-1)*obs_s2 + priorDF*priors$s2[1] + priors$n[1]*obs_n*(priors$m[1] - obs_m)^2/n1_n)/newDF # calculate new variance
      list(m = m1_n, s2 = s2_n, n = n1_n, DF = newDF)
    })
    names(distParams)<-names(datsp)
    
    p<-ggplot2::ggplot(data.frame(x=seq(0.1,1,0.1)))+
      facet_layer+
      geom_area(data=sub, aes(x=bin, y = freq), fill="gray40", color="gray40", show.legend=F)+
      lapply(1:length(distParams), function(i){ #* plot distribution of mean on top of total density
        pars<-distParams[[i]]
        lapDat<-data.frame(x=seq(0,1,0.1), y = 0, mu = pars$m)
        lapDat[[group[1] ]]=strsplit(names(distParams)[i],"[.]")[[1]][1]
        if(length(group)==2){lapDat[[group[2] ]] = strsplit(names(distParams)[i],"[.]")[[1]][2]}
        lapDat$fill <-lapDat[[group[1] ]]
        ggplot2::geom_vline(data = lapDat, ggplot2::aes(xintercept = mean(mu), color = .data$fill), linewidth=1)
      })+
      ggplot2::scale_x_continuous(n.breaks=5, labels = ~round(.,1))+
      ggplot2::coord_cartesian(ylim=c(0,15))+
      ggplot2::scale_fill_viridis_d(option="viridis")+
      ggplot2::guides(color="none")+#guide_legend(override.aes = list(linewidth=5)))
      ggplot2::labs(x=index, y=freq)
    
    
    support<-seq(floor(min(sub$bin, na.rm=T)), ceiling(max(sub$bin, na.rm=T)), length.out=1000*round(diff(range(sub$bin, na.rm=T)))  )
    
    dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
      pars<-distParams[[nm]]
      dens<-extraDistr::dlst(support, pars$DF, pars$m, sqrt(pars$s2))
      out<-data.frame(xdens= support, ydens=dens)
      out[,(ncol(out)+1):(ncol(out)+length(group))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
      colnames(out)<-c("xdens", "ydens", group)
      out$y<-out[[group[1] ]]
      out
    }))
    
    #* ***** `Compare generating gaussian distributions`
    if(doStats){
      outStats<-do.call(rbind, lapply(compareTests, function(comp){
        #comp=compareTests[[i]]
        g1<-as.character(comp[1])
        g2<-as.character(comp[2])
        
        m1<-distParams[[g1]]$m
        var1<-distParams[[g1]]$s2
        n1<-distParams[[g1]]$n
        DF1<-distParams[[g1]]$DF
        
        m2<-distParams[[g2]]$m
        var2<-distParams[[g2]]$s2
        n2<-distParams[[g2]]$n
        DF2<-distParams[[g2]]$DF
        
        post_mean_d1 <- extraDistr::dlst(support,DF1, m1, sqrt(var1/n1))
        post_mean_pdf1 <- post_mean_d1/sum(post_mean_d1)
        post_d1 <- extraDistr::dlst(support,DF1, m1, sqrt(var1))
        post_pdf1 <- post_d1/sum(post_d1)
        
        post_mean_d2 <- extraDistr::dlst(support,DF2, m2, sqrt(var2/n2))
        post_mean_pdf2 <- post_mean_d2/sum(post_mean_d2)
        post_d2 <- extraDistr::dlst(support,DF2, m2, sqrt(var2))
        post_pdf2 <- post_d2/sum(post_d2)
        
        if(hyp=="unequal"){
          post.prob.mu = 1-sum(apply(cbind(post_mean_pdf1,post_mean_pdf2),
                                     MARGIN=1,function(i) min(i)),na.rm=T) # P[pdf1 != pdf2]
          post.prob.dist = 1-sum(apply(cbind(post_pdf1,post_pdf2),
                                       MARGIN=1,function(i) min(i)),na.rm=T) # P[pdf1 != pdf2]
        } else if (hyp == "lesser"){
          direction.mu <- post_mean_pdf1 <= post_mean_pdf2
          post.prob.mu <- sum(post_mean_pdf1*direction.mu,na.rm=T)
          direction.dist <-  post_pdf1 <= post_pdf2
          post.prob.dist <- sum(post_pdf1*direction.dist,na.rm=T)
        } else if(hyp == "greater"){
          direction.mu <- post_mean_pdf1 >= post_mean_pdf2
          post.prob.mu <- sum(post_mean_pdf1*direction.mu,na.rm=T)
          direction.dist <-  post_pdf1 >= post_pdf2
          post.prob.dist <- sum(post_pdf1*direction.dist,na.rm=T)
        }
        ret<-data.frame(group1 = g1, group2 = g2, post.prob.mu = post.prob.mu,
                        post.prob.dist = post.prob.dist, method="gaussian", null = hyp)
        ret
      }))
    }
  }
  if(doStats & !is.null(method)){ return(list(p, outStats)) } else { return(p) }
}
