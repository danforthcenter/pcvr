#' Check priors used in ease of use brms functions
#' 
#' @param priors A named list of means for prior distributions.
#' This takes the same input as the prior argument of \code{\link{growthSS}}.
#' @param type Either "density", the default, or a model as would be specified in \code{growthSS}
#' or \code{growthSim} such as "logistic", "gompertz", "monomolecular", "exponential",
#' "linear", "power law", "double logistic", or "double gompertz".
#' If this is a model type then n draws from the prior will be simulated as growth
#' trendlines and densities will be plotted on margins for some distributions.
#' @param n Numeric, if type is a model then how many draws from the prior should be simulated?
#' @param t Numeric, time passed to growthSim. Defaults to 25 (the growthSim default).
#' @keywords Bayesian, brms, prior
#' @return A named list of plots showing prior distributions that \code{growthSS} would use,
#' optionally with a plot of simulated growth curves using draws from those priors.
#' @import ggplot2
#' @import patchwork
#' @importFrom stats rlnorm dlnorm
#' @examples 
#' 
#' ## Not run:
#'
#' set.seed(123)
#' priors = list("A" = c(100, 130), "B" = c(10, 8), "C" = c(0.2, 0.1))
#' plotPrior(priors)
#' 
#' plotPrior(priors, "gompertz")[[1]]
#' 
#' ## End(Not run:)
#'
#' @export

plotPrior<-function(priors, type = "density", n=200, t=25){
  
  densPlots<-lapply(1:length(priors), function(i){
    pri=priors[[i]]
    nm=names(priors)[i]
    
    pri_df <- do.call(rbind, lapply(1:length(pri), function(o){
      prio <- pri[o]
      max=ceiling(max(rlnorm(1000, log(max(pri)), 0.25))*1.1)
      support = seq(0,max,length.out=10000)
      dens = dlnorm(support, meanlog = log(prio), sdlog=0.25)
      pdf = dens/sum(dens)
      data.frame(support = support,
                 dens = pdf,
                 param=nm,
                 item = as.character(o))
    }))
    
    ggplot2::ggplot(pri_df, ggplot2::aes(x = .data$support, y = .data$dens, fill=.data$item, group=.data$item))+
      ggplot2::geom_polygon(alpha = 0.5)+
      ggplot2::theme_minimal()+ggplot2::labs(y="Density", title=nm, fill="Prior")
    
  })
  names(densPlots)<-names(priors)

  if(type=="density"){
    out<-densPlots
  } else{
    
    x_margin_plot = NULL
    y_margin_plot = NULL
    z_margin_plot = NULL
    
    simdf <- do.call(rbind, lapply(1:n,  function(i) {
      iter_params <- .prior_sampler(priors)
      x<-growthSim(model = type, n = 1, t = t, params = iter_params); x$id = paste0("id_",i);x } ))
    
    if(type=="logistic"){
      x_margin_plot = densPlots[["B"]] 
      y_margin_plot = densPlots[["A"]] 
      z_margin_plot = densPlots[["C"]] 
      } else if(type=="gompertz"){
      x_margin_plot = densPlots[["B"]]
      y_margin_plot = densPlots[["A"]]
      z_margin_plot = densPlots[["C"]] 
      } else if(type=="monomolecular"){
      y_margin_plot = densPlots[["A"]] 
      } 
    
    model_plot<-ggplot2::ggplot(simdf,
                    ggplot2::aes(x= .data$time, y=.data$y,group=interaction(.data$id, .data$group), color=.data$group))+
      ggplot2::geom_line(linewidth=0.1)+
      ggplot2::theme_minimal()+
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth=5)))+
      ggplot2::labs(y="Y", title=paste0(n," curves simulated from prior draws"),
                    color="Prior")+
      ggplot2::theme(legend.position=c(0.9,0.9))
    model_plot_solo<-model_plot
    xLims <- ggplot2::layer_scales(model_plot)$x$range$range
    yLims <- ggplot2::layer_scales(model_plot)$y$range$range
    
    if(!is.null(y_margin_plot)){
      y_margin_plot <- y_margin_plot + 
        ggplot2::scale_y_reverse(position = "right")+ 
        ggplot2::scale_x_continuous(position = "top", limits = yLims)+
        ggplot2::labs(x="Asymptote Prior")+
        ggplot2::coord_flip() + 
        ggplot2::theme(plot.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                       axis.text = ggplot2::element_blank(),axis.title.x = ggplot2::element_blank(),
                       legend.position = "none" )
    }
    
    if(!is.null(x_margin_plot)){
      x_margin_plot <- x_margin_plot + 
        ggplot2::labs(x="Inflection Point Prior")+
        ggplot2::theme(plot.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),axis.title.y = ggplot2::element_blank(),
                       legend.position = "none" )+
        ggplot2::coord_cartesian(xlim=xLims)
    }
    
    if(!is.null(z_margin_plot)){
      z_margin_plot <- z_margin_plot + 
        ggplot2::labs(x="Growth Rate Prior")+
        ggplot2::theme(plot.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),axis.title.y = ggplot2::element_blank(),
                       legend.position = "none" )+
        ggplot2::scale_x_continuous(n.breaks=3)
    }
    
    if(!is.null(x_margin_plot) & !is.null(y_margin_plot) & !is.null(z_margin_plot)){
      
      design = c(patchwork::area(1,1,6,6), # model plot
                 patchwork::area(7,1,7,6), # x margin
                 patchwork::area(1,7,6,7), # y margin
                 patchwork::area(7,7,7,7)) # "z" margin
      model_plot <- model_plot + x_margin_plot + y_margin_plot + z_margin_plot + patchwork::plot_layout(design = design)
      
    } else if(!is.null(y_margin_plot) & is.null(x_margin_plot)){

      design = c(patchwork::area(1,1,6,6), # model plot
                 patchwork::area(1,7,6,7)) # y margin
      model_plot <- model_plot + y_margin_plot + patchwork::plot_layout(design = design)
      
    } else if(is.null(y_margin_plot) & !is.null(x_margin_plot)){

      design = c(patchwork::area(1,1,6,6), # model plot
                 patchwork::area(7,1,7,6)) # x margin
      model_plot <- model_plot + x_margin_plot + patchwork::plot_layout(design = design)
      
    }
    if(is(model_plot, "patchwork")){
      densPlots[[length(densPlots)+1]] <- model_plot_solo
    }
    out<-list("simulated"= model_plot, "distributions" = densPlots)
  }

  return(out)
}

#' @description
#' Internal function for drawing from priors
#' @param priors priors as a list
#' @keywords internal
#' @noRd

.prior_sampler<-function(priors){
  lapply(priors, function(pri){ # draw sample from prior
    unlist(lapply(pri, function(mu){
      rlnorm(1,log(mu), 0.25)
    }))
  })
}
