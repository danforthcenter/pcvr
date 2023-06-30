#' Make Joyplots for multi value trait plantCV data and optionally compare distributions
#' 
#' @param df Data frame to use. Long or wide format is accepted.
#' @param index If the data is long then this is a multi value trait as a character string that must be present in `trait`. If the data is wide then this is a string used to find column names to use from the wide data. In the wide case this should include the entire trait name (ie, "hue_frequencies" instead of "hue_freq").
#' @param group A length 1 or 2 character vector. This is used for faceting the joyplot and identifying groups for testing. If this is length 1 then no faceting is done.
#' @param method A method to use in comparing distributions/means. Currently "beta", "gaussian", "emd", and "ks" are supported. See details for explanations of tests.
#' @param compare Groups to compare. By default this is set to FALSE, which corresponds to no testing. Other values of compare are passed to fixCompare to make t.test comparisons using ggpubr. In short, NULL will run all pairwise T tests, a single value of the X axis variable will compare that level to all other levels of the X variable, alternatively this can be a list as used by ggpubr: list(c("level1", "level2"), c("level1", "level3"))
#' @param priors Parameters for prior distributions if using method = "beta" or "gaussian". If left NULL then wide weak priors are used. This can be set as a single set of parameters or with unique values for each group in `compare`, but it is advised to use the default priors.
#' @param hyp Hypothesis to test for beta and gaussian methods, one of "unequal", "greater", "lesser". Defaults to "unequal" if left NULL.
#' @param bin Column containing histogram (multi value trait) bins. Defaults to "label".
#' @param freq Column containing histogram counts. Defaults to "value"
#' @param trait Column containing phenotype names. Defaults to "trait".
#' @param fillx Logical, whether or not to use \code{ggridges::geom_density_ridges_gradient} . Default is T, if F then \code{ggridges::geom_density_ridges} is used instead, with arbitrary fill. Note that \code{ggridges::geom_density_ridges_gradient} may issue a message about deprecated ggplot2 features.
#' @keywords bayesian, ggplot, multi value trait, pcv.hists
#' @import ggplot2
#' @import ggridges
#' 
#' @details 
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
#' wide_beta <- read.pcv("https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest2.csv", "wide", F)
#' pcv.joyplot(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"))
#' pcv.joyplot(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"), method="ks")
#' pcv.joyplot(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"), method="beta")
#' pcv.joyplot(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"), method="gaussian")
#' 
#' wide<-read.pcv("https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv", mode="wide", singleValueOnly = T, multiValPattern = "hist", reader="fread")
#' wide$genotype = substr(wide$barcode, 3,5)
#' wide$genotype = ifelse(wide$genotype == "002", "B73",
#'                        ifelse(wide$genotype == "003", "W605S",
#'                               ifelse(wide$genotype == "004", "MM", "Mo17")))
#' wide$fertilizer = substr(wide$barcode, 8, 8)
#' wide$fertilizer = ifelse(wide$fertilizer == "A", "100",
#'                          ifelse(wide$fertilizer == "B", "50", "0"))
#' p<-pcv.joyplot(wide, index = "hue_frequencies", group=c("genotype", "fertilizer"))
#' # For some color traits it makes sense to show the actual represented color, which can be done easily by adding new fill scales.
#' p+scale_fill_gradientn(colors = scales::hue_pal(l=65)(360))
#' @export


pcv.joyplot<-function(df = NULL, index = NULL, group = NULL,
                      method=NULL,
                      compare= NULL, priors=NULL,hyp=NULL,
                      bin="label", freq="value", trait="trait", fillx=T){
  
  #* ***** `troubleshooting test values`
  # df=df; index = "yii_hist_Fv.over.Fm"; group=c("timepoint", "genotype"); method="emd"
  # compare= NULL; priors=NULL; bin="label"; freq="value"; trait="trait";hyp=NULL; fillx=T    ;index='index_frequencies_index_ari' 
  #* `wide vs long test values`
  #* 
  #* df=hue_wide; index = "hue_frequencies"; group=c("fertilizer", "genotype"); method=NULL
  #* compare= F; priors=NULL; bin="label"; freq="value"; trait="trait";hyp=NULL; fillx=T 
  #* 
  #* df=long; index = "hue_frequencies"; group=c("fertilizer", "genotype"); method="ks"
  #* compare= NULL; priors=NULL; bin="label"; freq="value"; trait="trait";hyp=NULL; fillx=T 
  #* 
  #* ***** `general calculated values`

  if(!is.null(trait) && trait %in% colnames(df)){traitCol = trait ; mode="long" # if there is a trait column then use long options,
    }else{mode="wide"} # else use wide options
  
  #* if long data then subset rows where trait is correct
  if(mode=="long"){
    if(is.null(index)){sub<-df
    }else{ sub<-df[df[[trait]]==index, ] }
    if(length(unique(sub[[trait]]))>1){warning("More than one trait found, consider an `index` argument")}
    sub$bin = as.numeric(sub[[bin]])
    sub$freq = as.numeric(sub[[freq]])
    
  } else if(mode=="wide"){ # if wide then get column names that contain index string
    sub<-df
  }
  
  if(is.null(group)){group = "dummy"; df$dummy = "dummy"; sub$dummy="dummy"}
  if(length(group)==1){sub$fill = sub[[group]]; sub$y = sub[[group]]; facet_layer = list()}#ggplot2::facet_wrap(paste0("~",group)) }
  if(length(group)==2){sub$fill = sub[[group[1] ]]; sub$y = sub[[group[1] ]]; facet_layer=ggplot2::facet_grid(as.formula(paste0("~",group[2])))} # check this change, added group[1] to facet_grid
  
  sub$grouping<-interaction(sub[,c(group)], drop=T)
  
  # default compare to NULL, but if F then skip all testing 
  if(is.logical(compare) && compare==F){
    doStats=F
  } else if(is.null(compare)){compareTests<-fixCompare(compare,sub,"grouping", T) ; doStats=T
  }else{compareTests=fixCompare(compare,sub,"grouping"); doStats=T}
  if(is.null(hyp)){hyp<-"unequal"}
  
  #* ***** `default joyplot`
  if(is.null(method) | match.arg(method, choices = c("beta", "gaussian", "ks", "mixture", "emd"))=="emd"){
    if(mode=="wide"){
      o<-wide.dens.default(d=sub, colPattern = index, group_internal=group)
    } else if(mode=="long"){
      o<-long.dens.default(d=sub, group_internal=group)
    }
    distParams = o[[1]]
    dens_df = o[[2]]
    
    if(match.arg(method, choices = c("beta", "gaussian", "ks", "mixture", "emd"))=="emd" & doStats){
      outStats<-do.call(rbind, lapply(compareTests, function(comp){
        g1<-as.character(comp[1])
        g2<-as.character(comp[2])
        d1<-sub[sub$grouping == g1, ]
        d2<-sub[sub$grouping == g2, ]
        if(mode=="long"){
          d1s<-aggregate(as.formula(paste0(freq,"~",bin)), d1, sum)[[freq]]
          d2s<-aggregate(as.formula(paste0(freq,"~",bin)), d2, sum)[[freq]]
        } else if(mode=="wide"){
          histCols <- colnames(sub)[grepl(index, colnames(sub))]
          histCols_bin <- as.numeric(sub(paste0(index, "[.]?"), "", colnames(sub)[grepl(index, colnames(sub))]))
          bins_order<-sort(histCols_bin, index.return=T)$ix
          histCols <- histCols[bins_order]
          d1s<-colSums(d1[,histCols],na.rm=T)
          d2s<-colSums(d2[,histCols],na.rm=T)
        }
        emd_res <- emd1d(d1s, d2s)
        data.frame(group1 = g1, group2 = g2, emd = emd_res, method="emd")
      }))
    }
  } else if(match.arg(method, choices = c("beta", "gaussian", "ks", "mixture", "emd"))=="ks"){
    #* ***** `Non parametric joyplot`
    
    if(mode=="wide"){
      o<-wide.dens.default(d=sub, colPattern = index, group_internal=group)
    } else if(mode=="long"){
      o<-long.dens.default(d=sub, group_internal=group)
    }
    distParams = o[[1]]
    dens_df = o[[2]]
    
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
    
    if(mode=="wide"){
      o<-wide.dens.beta(d=sub, colPattern = index, group_internal=group)
    } else if(mode=="long"){
      o<-long.dens.beta(d=sub, group_internal=group)
    }
    distParams = o[[1]]
    dens_df = o[[2]]
    
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
    if(mode=="wide"){
      o<-wide.dens.gaussian(d=sub, colPattern = index, group_internal=group)
    } else if(mode=="long"){
      o<-long.dens.gaussian(d=sub, group_internal=group)
    }
    distParams = o[[1]]
    dens_df = o[[2]]
    
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
  if(fillx){dens_df$group = interaction(dens_df[,group])}
  ggridgeLayer<-if(fillx){
    list(suppressMessages(ggridges::geom_density_ridges_gradient(ggplot2::aes(x = xdens, y = y,
                                                        height =ydens,fill=stat(x)),
                                           show.legend=F, stat="identity", rel_min_height = 0.001)),
         ggplot2::scale_fill_viridis_c(option="plasma",
                                        limits = c(min(dens_df$xdens[dens_df$ydens>0.001],na.rm=T),
                                                   max(dens_df$xdens[dens_df$ydens>0.001], na.rm=T)))
         )
  } else{
    list(suppressMessages(ggridges::geom_density_ridges2(ggplot2::aes(x = xdens, y = y, height =ydens, fill = y, color=y),
                                   show.legend=F, stat="identity")),
      ggplot2::scale_color_viridis_d(option="viridis"),
      ggplot2::scale_fill_viridis_d(option="viridis")
      )
  }
  # sub2<-dens_df[dens_df$timepoint=="tp4" & dens_df$genotype=="G6",]
  # ggplot(sub2, aes(x=xdens, y=ydens))+geom_line()
  p<-ggplot2::ggplot(dens_df)+
    facet_layer+
    ggridgeLayer+
    ggplot2::scale_x_continuous(n.breaks=5, labels = ~round(.,1))+
    ggplot2::labs(x=index, y=group[1])+
    pcv_theme()+
    ggplot2::theme(legend.position="none")

  if(doStats & !is.null(method)){ return(list("plot" = p, "stats"=outStats)) } else { return(p) }
}



long.dens.default <- function(d = NULL, group_internal=NULL){
  datsp=split(d, d$grouping, drop=T)
  bw<-min(diff(sort(as.numeric(unique(d$bin )))))*0.75
  distParams<-lapply(datsp, function(D){
    X1 <- as.numeric(D[rep(rownames(D), round(D$freq)), bin])
    dens<-density(X1, from = min(d$bin,na.rm=T), to = max(d$bin,na.rm=T), n = 2^10, bw=bw)
    return(dens)
  })
  names(distParams)<-names(datsp)
  dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
    dens<-distParams[[nm]]
    out<-data.frame(xdens= dens$x, ydens=dens$y)
    out[,(ncol(out)+1):(ncol(out)+length(group_internal))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
    colnames(out)<-c("xdens", "ydens", group_internal)
    out$y<-out[[group_internal[1] ]]
    out
  }))
  return(list(distParams, dens_df))
}

wide.dens.default<-function(d=NULL, colPattern = NULL, group_internal=NULL){
  histCols <- colnames(d)[grepl(colPattern, colnames(d))]
  histCols_bin <- as.numeric(sub(paste0(colPattern, "[.]?"), "", colnames(d)[grepl(colPattern, colnames(d))]))
  bins_order<-sort(histCols_bin, index.return=T)$ix
  histCols <- histCols[bins_order]
  datsp=split(d, d$grouping, drop=T)
  bw<-min(diff(sort(as.numeric(histCols_bin))))*0.75
  
  distParams<-lapply(datsp, function(D){
    D<-D[,histCols]
    X1<-rep(histCols_bin[bins_order], as.numeric(round(colSums(D))) )
    dens<-density(X1, from = min(histCols_bin,na.rm=T), to = max(histCols_bin,na.rm=T), n = 2^10, bw=bw)
    return(dens)
  })
  
  names(distParams)<-names(datsp)
  dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
    dens<-distParams[[nm]]
    out<-data.frame(xdens= dens$x, ydens=dens$y)
    out[,(ncol(out)+1):(ncol(out)+length(group_internal))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
    colnames(out)<-c("xdens", "ydens", group_internal)
    out$y<-out[[group_internal[1] ]]
    out
  }))
  return(list(distParams, dens_df))
}

long.dens.beta<-function(d = NULL, group_internal=NULL){
  datsp=split(d, d$grouping, drop=T) # split data into panel groups
  
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
  support<-seq(0.001, 0.999, 0.001)
  
  dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
    pars<-distParams[[nm]]
    dens<-dbeta(support, pars[1], pars[1])
    out<-data.frame(xdens= support, ydens=dens)
    out[,(ncol(out)+1):(ncol(out)+length(group_internal))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
    colnames(out)<-c("xdens", "ydens", group_internal)
    out$y<-out[[group_internal[1] ]]
    out
  }))
  return(list(distParams, dens_df))
}

wide.dens.beta<-function(d = NULL, colPattern = NULL, group_internal=NULL){
  histCols <- colnames(d)[grepl(colPattern, colnames(d))]
  histCols_bin <- as.numeric(sub(paste0(colPattern, "[.]?"), "", colnames(d)[grepl(colPattern, colnames(d))]))
  bins_order<-sort(histCols_bin, index.return=T)$ix
  histCols <- histCols[bins_order]
  datsp=split(d, d$grouping, drop=T)
  if(is.null(priors)){ priors = list(a = rep(0.5, times=length(datsp)), b = rep(0.5, times=length(datsp))) } # assume weak prior on everything
  bw<-min(diff(sort(as.numeric(histCols_bin))))*0.75
  distParams<-lapply(1:length(datsp), function(i){
    D<-datsp[[i]]
    D<-D[,histCols]
    X1<-rep(histCols_bin[bins_order], as.numeric(round(colSums(D))) )
    n1<- length(X1)
    xbar1 = (1/n1) * sum(X1)
    variance1 = (1/(n1))*sum( (X1 - xbar1)^2 )
    alpha1 = (xbar1 * ( (xbar1 * (1-xbar1) )/(variance1)-1 )) + priors$a[i]
    beta1 = ((1-xbar1) * ( (xbar1 * (1-xbar1) )/(variance1)-1 )) + priors$b[i]
    return(c(alpha1, beta1))
  })
  
  dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
    pars<-distParams[[nm]]
    dens<-dbeta(support, pars[1], pars[1])
    out<-data.frame(xdens= support, ydens=dens)
    out[,(ncol(out)+1):(ncol(out)+length(group_internal))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
    colnames(out)<-c("xdens", "ydens", group_internal)
    out$y<-out[[group_internal[1] ]]
    out
  }))
  return(list(distParams, dens_df))
}


long.dens.gaussian<-function(d=NULL, group_internal=NULL){
  datsp=split(d, d$grouping, drop=T)
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
  support<-seq(floor(min(d$bin, na.rm=T)), ceiling(max(d$bin, na.rm=T)), length.out=1000*round(diff(range(d$bin, na.rm=T)))  )
  
  dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
    pars<-distParams[[nm]]
    dens<-extraDistr::dlst(support, pars$DF, pars$m, sqrt(pars$s2))
    out<-data.frame(xdens= support, ydens=dens)
    out[,(ncol(out)+1):(ncol(out)+length(group_internal))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
    colnames(out)<-c("xdens", "ydens", group_internal)
    out$y<-out[[group_internal[1] ]]
    out
  }))
  return(list(distParams, dens_df))
}


wide.dens.gaussian<-function(d=NULL, colPattern=NULL, group_internal=NULL ){
  histCols <- colnames(d)[grepl(colPattern, colnames(d))]
  histCols_bin <- as.numeric(sub(paste0(colPattern, "[.]?"), "", colnames(d)[grepl(colPattern, colnames(d))]))
  bins_order<-sort(histCols_bin, index.return=T)$ix
  histCols <- histCols[bins_order]
  datsp=split(d, d$grouping, drop=T)
  
  if(is.null(priors)){ priors <- list( m=rep(0,length(datsp)), n=rep(1,length(datsp)), s2=rep(20,length(datsp)) ) } # mean, number, and variance
  distParams<-lapply(1:length(datsp), function(i){ 
    D<-datsp[[i]]
    D<-D[,histCols]
    X1<-rep(histCols_bin[bins_order], as.numeric(round(colSums(D))) )
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
  support<-seq(floor(min(d$bin, na.rm=T)), ceiling(max(d$bin, na.rm=T)), length.out=1000*round(diff(range(d$bin, na.rm=T)))  )
  
  dens_df<-do.call(rbind, lapply(names(distParams), function(nm){
    pars<-distParams[[nm]]
    dens<-extraDistr::dlst(support, pars$DF, pars$m, sqrt(pars$s2))
    out<-data.frame(xdens= support, ydens=dens)
    out[,(ncol(out)+1):(ncol(out)+length(group_internal))]<-lapply(strsplit(nm, "[.]")[[1]], identity)
    colnames(out)<-c("xdens", "ydens", group_internal)
    out$y<-out[[group_internal[1] ]]
    out
  }))
  return(list(distParams, dens_df))
}





