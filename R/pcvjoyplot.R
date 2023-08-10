#' Make Joyplots for multi value trait plantCV data and optionally compare distributions
#' 
#' @param df Data frame to use. Long or wide format is accepted.
#' @param index If the data is long then this is a multi value trait as a
#' character string that must be present in `trait`. 
#' If the data is wide then this is a string used to find column names to use from the wide data.
#'  In the wide case this should include the entire
#'   trait name (ie, "hue_frequencies" instead of "hue_freq").
#' @param group A length 1 or 2 character vector. 
#' This is used for faceting the joyplot and identifying groups for testing. 
#' If this is length 1 then no faceting is done.
#' @param y Optionally a variable to use on the y axis. This is useful when you
#' have three variables to display. This argument will change faceting behavior to
#' add an additional layer of faceting (single length group will be faceted, 
#' length 2 group will be faceted group1 ~ group2).
#' @param method A method to use in comparing distributions/means.
#'  Currently only "ks" is supported, where density is approximated from the 
#'  histogram columns and  KS test is used between those densities. For other
#'  options in comparing multi-value traits see \code{\link{conjugate}} or
#'  \code{\link{pcv.emd}}.
#' @param compare Groups to compare. By default this is set to FALSE, 
#' which corresponds to no testing. Other values of compare are passed to 
#' fixCompare to make t.test comparisons using ggpubr. 
#' In short, NULL will run all pairwise T tests, a single value of the X axis 
#' variable will compare that level to all other levels of the X variable, 
#' alternatively this can be a list as used by ggpubr: list(c("level1", "level2"),
#'  c("level1", "level3"))
#' @param bin Column containing histogram (multi value trait) bins. Defaults to "label".
#' @param freq Column containing histogram counts. Defaults to "value"
#' @param trait Column containing phenotype names. Defaults to "trait".
#' @param fillx Logical, whether or not to use \code{ggridges::geom_density_ridges_gradient}.
#'  Default is T, if F then \code{ggridges::geom_density_ridges} is used instead,
#'   with arbitrary fill. Note that \code{ggridges::geom_density_ridges_gradient} 
#'   may issue a message about deprecated ggplot2 features.
#' @keywords bayesian, ggplot, multi value trait, pcv.hists
#' @import ggplot2
#' @import ggridges
#' @importFrom stats setNames density aggregate as.formula ks.test
#' 
#' @details 
#' The method argument is used for statistical testing of groups.
#' There there is only a "ks" option, with other methods in place for 
#' multi-value trait analysis available in \code{conjugate} and \code{pcv.emd}.
#' \itemize{
#'  \item{"ks": }{The ks method performs a ks test on the PDF of each histogram.
#'  This returns a P value corresponding to a 
#'  standard KS test that distributions are the same.}
#' }

#' 
#' @return Returns either a ggplot object or a list containing a ggplot and a dataframe of statistical comparisons (if compare is not FALSE).
#' 
#' @examples 
#' 
#' ## Not run: 
#' 
#' df <- read.pcv(
#'   "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest2.csv", "long", FALSE)
#' wide_beta <- read.pcv(
#'   "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest2.csv", "wide", FALSE)
#' x <- pcv.joyplot(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"))
#' # pcv.joyplot(df, index = "index_frequencies_index_ndvi", group=c("genotype", "timepoint"),
#' #   method="ks")
#' 
#' wide<-read.pcv(
#' "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv",
#'    mode="wide", singleValueOnly =TRUE, multiValPattern = "hist", reader="fread")
#' wide <- bw.time(wide, mode="DAS",plot=FALSE)
#' wide$genotype = substr(wide$barcode, 3,5)
#' wide$genotype = ifelse(wide$genotype == "002", "B73",
#'                        ifelse(wide$genotype == "003", "W605S",
#'                               ifelse(wide$genotype == "004", "MM", "Mo17")))
#' wide$fertilizer = substr(wide$barcode, 8, 8)
#' wide$fertilizer = ifelse(wide$fertilizer == "A", "100",
#'                          ifelse(wide$fertilizer == "B", "50", "0"))
#' p<-pcv.joyplot(wide[wide$DAS > 15,], index = "hue_frequencies",
#'  group=c("genotype", "fertilizer"), y ="DAS")
#' 
#' # For some color traits it makes sense to show the actual
#' # represented color, which can be done easily by adding new fill scales.
#' p+ggplot2::scale_fill_gradientn(colors = scales::hue_pal(l=65)(360))
#' 
#' ## End(Not run)
#' 
#' @export


pcv.joyplot<-function(df = NULL, index = NULL, group = NULL, y = NULL,
                      method=NULL, compare= NULL,
                      bin="label", freq="value", trait="trait", fillx=TRUE){
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
    if(all(as.matrix(sub[sub[[trait]] == index, freq])<2)){
      sub$freq = sub$freq * 1000/max(sub[sub[[trait]] == index, "freq"])
      if(!is.null(method)){
        warning("Data is being rescaled so that the max count is 1000 for plotting, this can change statistical results")
      }
    }
    
  } else if(mode=="wide"){ # if wide then get column names that contain index string
    sub<-df
    if(all(as.matrix(sub[, grepl(index, colnames(sub))])<2)){
      sub[, grepl(index, colnames(sub))]<-sub[, grepl(index, colnames(sub))] * 1000 / max(rowSums( sub[, grepl(index, colnames(sub))] ))
      if(!is.null(method)){
        warning("Data is being rescaled so that the max count is 1000 for plotting, this can change statistical results")
      }
    }
  }
  
  if(is.null(group)){group = "dummy"; df$dummy = "dummy"; sub$dummy="dummy"}
  if(!is.null(y)){
    if(length(group)==1){sub$fill = sub[[group]]; sub$y = sub[[group]]
    facet_layer=ggplot2::facet_grid(as.formula(paste0("~",group[1]))) }
    if(length(group)==2){sub$fill = sub[[group[1] ]]; sub$y = sub[[group[1] ]]
    facet_layer=ggplot2::facet_grid(as.formula(paste0(group[1], "~",group[2]))) }
  } else { # if y is not provided then one less layer of faceting
    if(length(group)==1){sub$fill = sub[[group]]; sub$y = sub[[group]]
    facet_layer=list() }
    if(length(group)==2){sub$fill = sub[[group[1] ]]; sub$y = sub[[group[1] ]]
    facet_layer=ggplot2::facet_grid(as.formula(paste0("~",group[2]))) }
  }
  
  sub$grouping<-interaction(sub[,c(y,group)], drop=TRUE)
  
  # default compare to NULL, but if F then skip all testing 
  if(is.logical(compare) && compare==FALSE){
    doStats=FALSE
  } else if(is.null(compare)){compareTests<-fixCompare(compare,sub,"grouping", TRUE) ; doStats=TRUE
  }else{compareTests=fixCompare(compare,sub,"grouping"); doStats=TRUE}
  
  #* ***** `default joyplot`
    if(mode=="wide"){
      o<-.wide.dens.default(d=sub, colPattern = index, group_internal=c(y,group))
    } else if(mode=="long"){
      o<-.long.dens.default(d=sub, group_internal=c(y,group), bin_internal = bin, freq_internal= freq)
    }
    distParams = o[[1]]
    dens_df = o[[2]]
    
  if(match.arg(method, choices = c("beta", "gaussian", "ks", "mixture", "emd"))=="ks"){
    #* ***** `default joyplot with KS tests`
    if(mode=="wide"){
      o<-.wide.dens.default(d=sub, colPattern = index, group_internal=c(y,group))
    } else if(mode=="long"){
      o<-.long.dens.default(d=sub, group_internal=c(y,group), bin_internal= bin, freq_internal= freq)
    }
    distParams = o[[1]]
    dens_df = o[[2]]
    
    if(doStats){
      outStats<-do.call(rbind, lapply(compareTests, function(comp){
        g1<-as.character(comp[1])
        g2<-as.character(comp[2])
        ks<-suppressWarnings(ks.test(distParams[[g1]]$y, distParams[[g2]]$y)) # 0.06 with el, 2e-16 with el_r
        ret<-data.frame(group1 = g1, group2 = g2, p = ks$p.value, method="ks", null = "same distribution")
        return(ret)
      }))
    }
  } 
  if(fillx){dens_df$group = interaction(dens_df[,group])}
  ggridgeLayer<-if(fillx){
    x<-NULL # to make R CMD check happy with stat(x)
    list(suppressMessages(ggridges::geom_density_ridges_gradient(ggplot2::aes(x = .data$xdens, y = .data$y,
                                                        height = .data$ydens, fill=stat(x)),
                                           show.legend=FALSE, stat="identity", rel_min_height = 0.001)),
         ggplot2::scale_fill_viridis_c(option="plasma",
                                        limits = c(min(dens_df$xdens[dens_df$ydens>0.001],na.rm=TRUE),
                                                   max(dens_df$xdens[dens_df$ydens>0.001], na.rm=TRUE)))
         )
  } else{
    list(suppressMessages(ggridges::geom_density_ridges2(ggplot2::aes(x = .data$xdens, y = .data$y,
                                                                      height =.data$ydens, fill = .data$y, color=.data$y),
                                   show.legend=FALSE, stat="identity")),
      ggplot2::scale_color_viridis_d(option="viridis"),
      ggplot2::scale_fill_viridis_d(option="viridis")
      )
  }
  p<-ggplot2::ggplot(dens_df)+
    facet_layer+
    ggridgeLayer+
    ggplot2::scale_x_continuous(n.breaks=5, labels = ~round(.,1))+
    ggplot2::labs(x=index, y=c(y,group)[1])+
    pcv_theme()+
    ggplot2::theme(legend.position="none")

  if(doStats & !is.null(method)){ return(list("plot" = p, "stats"=outStats)) } else { return(p) }
}

#' ***********************************************************************************************
#' *************** `Density for Long Data` ****************************************
#' ***********************************************************************************************
#' 
#' @description
#' Internal function for making density of histogram data and returning it in a standard format for 
#' use in KS tests and plotting
#' @param d data with some manipulations to make groups consistent passed from pcv.joyplot
#' @param group_internal grouping passed from pcv.joyplot, c(group, y)
#' @param bin_internal column name to use for bins
#' @param freq_internal column name to use for frequencies (counts) in each bin
#' 
#' @keywords internal
#' @noRd


.long.dens.default <- function(d = NULL, group_internal=NULL, bin_internal= NULL, freq_internal= NULL){
  datsp=split(d, d$grouping, drop=TRUE)
  bw<-min(diff(sort(as.numeric(unique(d$bin )))))*0.75
  distParams<-lapply(datsp, function(D){
    X1 <- as.numeric(D[rep(rownames(D), round(D[[freq_internal]])), bin_internal])
    dens <- density(d$bin, weights = d$freq/sum(d$freq),
                    from = min(d$bin,na.rm=TRUE),
                    to = max(d$bin,na.rm=TRUE),
                    n = 2^10)
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

#' ***********************************************************************************************
#' *************** `Density for Wide Data` ****************************************
#' ***********************************************************************************************
#' 
#' @description
#' Internal function for making density of histogram data and returning it in a standard format for 
#' use in KS tests and plotting
#' @param d data with some manipulations to make groups consistent passed from pcv.joyplot
#' @param colPattern index identifying pattern passed from pcv.joyplot
#' @param group_internal grouping passed from pcv.joyplot, c(group, y)
#' 
#' @keywords internal
#' @noRd


.wide.dens.default<-function(d=NULL, colPattern = NULL, group_internal=NULL){
  histCols <- colnames(d)[grepl(colPattern, colnames(d))]
  histCols_bin <- as.numeric(sub(paste0(colPattern, "[.]?"), "", colnames(d)[grepl(colPattern, colnames(d))]))
  bins_order<-sort(histCols_bin, index.return=TRUE)$ix
  histCols <- histCols[bins_order]
  datsp=split(d, d$grouping, drop=TRUE)
  bw<-min(diff(sort(as.numeric(histCols_bin))))*0.75
  
  distParams<-lapply(datsp, function(D){
    D<-D[,histCols]
    X1<-rep(histCols_bin[bins_order], as.numeric(round(colSums(D))) ) 
    weights <- colSums(D)/sum(colSums(D))
    
    dens<-density(x = histCols_bin[bins_order],
                  weights = weights,
                  from = min(histCols_bin,na.rm=TRUE),
                  to = max(histCols_bin,na.rm=TRUE),
                  n = 2^10)
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




