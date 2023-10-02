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
#' @param id Optionally a variable to show the outline of different replicates.
#' Note that ggridges::geom_density_ridges_gradient does not support transparency,
#' so if fillx is TRUE then only the outer line will show individual IDs.
#' @param method A method to use in comparing distributions/means.
#'  Currently "ks" and "pdf" are supported, where density is approximated from the 
#'  histogram columns and  KS test is used between samples from those densities per each photo.
#'  For the "pdf" method a flat prior is added to the data and the distributions are compared 
#'  depending on the hypothesis provided. For other
#'  options in comparing multi-value traits see \code{\link{conjugate}} or
#'  \code{\link{pcv.emd}}.
#' @param hypothesis A hypothesis for the "pdf" method, must be either unequal, equal, lesser, or greater.
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
#' @import data.table
#' @importFrom stats setNames density aggregate as.formula ks.test
#' 
#' 
#' @return Returns either a ggplot object or a list containing a ggplot and a dataframe of statistical comparisons (if compare is not FALSE).
#' 
#' @examples 
#' 
#' ## Not run: 
#' 
#' df <- read.pcv(
#'   "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest2.csv", "long")
#'   
#' x <- pcv.joyplot(df, index = "index_frequencies_index_ndvi",
#'     group=c("genotype", "timepoint"), method="pdf")
#' 
#' if (FALSE){
#' wide<-read.pcv(
#'   paste0("https://media.githubusercontent.com/media/joshqsumner/",
#'         "pcvrTestData/main/pcv4-multi-value-traits.csv"),
#'    mode="wide")
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
#' }
#' ## End(Not run)
#' 
#' @export


pcv.joyplot<-function(df = NULL, index = NULL, group = NULL, y = NULL, id = NULL,
                      method=NULL, hypothesis = "unequal", compare= NULL,
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
    sub <- sub[, c(group, y, id, "bin", "freq", trait)]
    # if(all(as.matrix(sub[sub[[trait]] == index, freq])<2)){
    #   sub$freq = sub$freq * 1000/max(sub[sub[[trait]] == index, "freq"])
    #   if(!is.null(method)){
    #     warning("Data is being rescaled so that the max count is 1000 for plotting, this can change statistical results")
    #   }
    # }
    
  } else if(mode=="wide"){ # if wide then get column names that contain index string
    #* subset data to only have index columns
    #* turn the data longer
    sub_wide<-data.table::as.data.table(df[, which(colnames(df)%in%c(group, y, id) | grepl(index, colnames(df)) ) ])
    sub <- as.data.frame(data.table::melt(sub_wide, id.vars = c(group, y, id), variable.name = trait, value.name = freq))
    sub[[bin]] <- sub(index,"",sub[[trait]])
    sub$bin <- as.numeric(regmatches(sub[[bin]], regexpr("[0-9].*", sub[[bin]]) ))
    sub[[trait]] <- index
    sub$freq = as.numeric(sub[[freq]])
    sub <- sub[, c(group, y, id, "bin", "freq", trait)]
    
    # if(all(as.matrix(sub[, grepl(index, colnames(sub))])<2)){
    #   sub[, grepl(index, colnames(sub))]<-sub[, grepl(index, colnames(sub))] * 1000 / max(rowSums( sub[, grepl(index, colnames(sub))] ))
    #   if(!is.null(method)){
    #     warning("Data is being rescaled so that the max count is 1000 for plotting, this can change statistical results")
    #   }
    # }
  }
  
  if(is.null(group)){group = "dummy"; df$dummy = "dummy"; sub$dummy="dummy"}
  if(!is.null(y)){
    if(length(group)==1){sub$y = sub[[y]]
    facet_layer=ggplot2::facet_grid(as.formula(paste0("~",group[1]))) }
    if(length(group)==2){sub$y = sub[[y ]]
    facet_layer=ggplot2::facet_grid(as.formula(paste0(group[1], "~",group[2]))) }
  } else { # if y is not provided then one less layer of faceting
    if(length(group)==1){sub$y = sub[[group]]
    facet_layer=list() }
    if(length(group)==2){sub$y = sub[[group[1] ]]
    facet_layer=ggplot2::facet_grid(as.formula(paste0("~",group[2]))) }
  }
  sub$y<-as.character(sub$y)
  sub$grouping<-interaction(sub[,c(y,group)], drop=TRUE)
  
  # default compare to NULL, but if F then skip all testing 
  if(is.logical(compare) && compare==FALSE){
    doStats=FALSE
  } else if(is.null(compare)){compareTests<-fixCompare(compare,sub,"grouping", TRUE) ; doStats=TRUE
  }else{compareTests=fixCompare(compare,sub,"grouping"); doStats=TRUE}
  
  if(!is.null(method) && match.arg(method, choices = c("ks", "pdf"))=="ks" ){
    
    #* ***** `Run KS tests`
    ksVectors <- .makeKSdata(d = sub)

    if(doStats){
      outStats<-do.call(rbind, lapply(compareTests, function(comp){
        g1<-as.character(comp[1])
        g2<-as.character(comp[2])
        ks<-suppressWarnings(ks.test(ksVectors[[g1]], ksVectors[[g2]])) 
        ret<-data.frame(group1 = g1, group2 = g2, p = ks$p.value, method="ks", hypothesis = "same distribution")
        return(ret)
      }))
    }
  } else if(!is.null(method) && match.arg(method, choices = c("ks", "pdf"))=="pdf" ){
  #* ***** `Run PDF comparisons`
  
    pdfs <- .makePDFs(d = sub)
    if(doStats){
      outStats<-do.call(rbind, lapply(compareTests, function(comp){
        g1<-as.character(comp[1])
        g2<-as.character(comp[2])
        prob <- .post.prob.from.pdfs(pdfs[[g1]]$pdf, pdfs[[g2]]$pdf, hypothesis)
        ret<-data.frame(group1 = g1, group2 = g2, p = prob$post.prob, method="pdf", hypothesis = hypothesis)
        return(ret)
      }))
    }
  }
  
  
  #* `if ID is null then aggregate, else draw with ID`
  if(is.null(id)){
    sub <- stats::aggregate(freq ~ . , data=sub, FUN=mean, na.rm=T)
    gg <- ggplot2::ggplot(sub)
  } else{
    sub$id <- sub[[id]]
    gg <-  ggplot2::ggplot(sub, ggplot2::aes(alpha = 0.5, group = interaction(id, y, grouping)))
  }

  #if(fillx){dens_df$group = interaction(dens_df[,group])}
  ggridgeLayer<-if(fillx){
    x<-NULL # to make R CMD check happy with stat(x)
    list(suppressMessages(ggridges::geom_density_ridges_gradient(ggplot2::aes(x = .data$bin, y = .data$y,
                                                        height = .data$freq, fill=stat(x)),
                                           show.legend=FALSE, stat="identity", rel_min_height = 0.001)),
         ggplot2::scale_fill_viridis_c(option="plasma"#,
                                        #limits = c(min(dens_df$xdens[dens_df$ydens>0.001],na.rm=TRUE),
                                        #           max(dens_df$xdens[dens_df$ydens>0.001], na.rm=TRUE))
                                       )
         )
  } else{
    list(suppressMessages(ggridges::geom_density_ridges2(ggplot2::aes(x = .data$bin, y = .data$y,
                                                                      height =.data$freq, fill = .data$bin, color=.data$bin),
                                   show.legend=FALSE, stat="identity")),
      ggplot2::scale_color_viridis_d(option="viridis"),
      ggplot2::scale_fill_viridis_d(option="viridis")
      )
  }
  p<-gg+
    facet_layer+
    ggridgeLayer+
    ggplot2::scale_x_continuous(n.breaks=5, labels = ~round(.,1))+
    ggplot2::labs(x=index, y=c(y,group)[1])+
    pcv_theme()+
    ggplot2::theme(legend.position="none")

  if(doStats & !is.null(method)){ return(list("plot" = p, "stats"=outStats)) } else { return(p) }
}

#' ***********************************************************************************************
#' *************** `KS test vectors` ****************************************
#' ***********************************************************************************************
#' 
#' @description
#' Internal function for making density of histogram data and returning it in a standard format for 
#' use in KS tests. Currently picking how many samples to pull from the distribution is arbitrary.
#' @param d data with some manipulations to make groups consistent passed from pcv.joyplot
#' 
#' @keywords internal
#' @noRd

.makeKSdata<-function(d = NULL){
  datsp=split(d, d$grouping, drop=TRUE)
  ksVectors<-lapply(datsp, function(datsp_iter){
    dens <- density(datsp_iter$bin, weights = datsp_iter$freq/sum(datsp_iter$freq),
                    from = min(d$bin,na.rm=TRUE), # min and max of total data
                    to = max(d$bin,na.rm=TRUE),
                    n = 2^10)
    den_df <- as.data.frame(dens)
    set.seed(123)
    n_bins <- length(unique(datsp_iter$bin))
    n_photos <- nrow(datsp_iter) / n_bins 
    # n <- n_bins * n_photos # should always equal nrow(datsp_iter)
    vec <- sample(den_df$x, size = n_photos, replace=T, prob = den_df$y)
    return(vec)
  })
  return(ksVectors)
}

# d = sub
.makePDFs <- function(d=NULL){
  datsp=split(d, d$grouping, drop=TRUE)
  pdfs<-lapply(datsp, function(datsp_iter){
    ag_df <- stats::aggregate(freq ~ bin+grouping, data=datsp_iter, mean)
    dens <- stats::density(ag_df$bin,
                    weights = (ag_df$freq+ 1/nrow(ag_df))/sum((ag_df$freq+ 1/nrow(ag_df))), # weighting + flat prior
                    from = min(d$bin,na.rm=TRUE), # min and max of total data
                    to = max(d$bin,na.rm=TRUE),
                    n = 2^10)
    dens_df <- data.frame(support = dens$x, dens = dens$y)
    dens_df$pdf <- dens_df$dens/sum(dens_df$dens)
    return(dens_df)
  })
  return(pdfs)
}


