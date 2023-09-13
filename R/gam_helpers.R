#' Helper function for visualizing differences in GAMs fit with \code{mgcv::gam}
#' 
#' @param model A model fit with smooth terms by \code{mgcv::gam}
#' @param support A data.frame of new data to use to make predictions. This should be dense support. See examples.
#' @param g1 A character string for the level of byVar to use as the first group to compare
#' @param g2 The second group to compare (comparison will be g1 - g2)
#' @param byVar Categorical variable name used to separate splines as a string.
#' @param smoothVar The variable that splines were used on. This will often be a time variable.
#' @param cis Confidence interval levels, can be multiple. For example, 0.95 would return Q_0.025 and Q_0.975 columns,
#' and c(0.9, 0.95) would return Q_0.025, Q_0.05, Q_0.95, and Q_0.975 columns. Defaults to \code{seq(0.05, 0.95, 0.05)}
#' @param unconditional Logical, should unconditional variance-covariance be used in calculating standard errors.
#' Defaults to TRUE.
#' @param plot Logical, should a plot of the difference be returned? Defaults to TRUE.
#' 
#' @keywords growth, mgcv, spline, GAM
#' @import mgcv
#' @importFrom stats vcov predict df.residual qt
#' @importFrom viridis mako
#' @import ggplot2
#' @examples 
#' 
#' 
#' ex <- mgcv::gamSim(eg=4,n=400,dist="normal",scale=2,verbose=FALSE)
#' 
#' m <- mgcv::gam(y ~ fac + s(x1, by = fac), data=ex)
#' 
#' support <- expand.grid(x1 = seq(min(ex$y), max(ex$y), length = 400),
#'                        fac = factor(levels(ex$fac)))
#' 
#' out<-spline_diff(m, support, "1", "2", byVar = "fac", smoothVar="x1")
#' dim(out$data)
#' out$plot
#' 

spline_diff <- function(model, support, g1, g2, byVar, smoothVar, cis = seq(0.05, 0.95, 0.05),
                        unconditional = TRUE, plot=TRUE) {
  
  xp <- stats::predict(model, newdata = support, type = 'lpmatrix')
  c1 <- grepl(g1, colnames(xp))
  c2 <- grepl(g2, colnames(xp))
  r1 <- support[[byVar]] == g1
  r2 <- support[[byVar]] == g2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% stats::vcov(model, unconditional = unconditional)) * X))
  df.resid <- stats::df.residual(model)
  
  cis_df <- do.call(cbind, lapply(cis, function(ci){
    crit <- stats::qt( 1-((1-ci)/2) , df.resid, lower.tail = TRUE)
    upr <- dif + (crit * se)
    lwr <- dif - (crit * se)
    setNames(data.frame(upr, lwr), paste0(c("Q_", "Q_"), round(c( 1-((1-ci)/2) , 1-(1-((1-ci)/2))), 3) ))
  }))
  cis_df <- cis_df[,order(colnames(cis_df))]
  #if(length(cis)==1){ colnames(cis_df) <- c("Q_lower", "Q_upper") }
  
  out_df<-cbind(data.frame(g1=g1, g2=g2, mu = dif, se = se), cis_df)
  out_df$df.resid <- df.resid
  
  smoothVarRange <- range(support[[smoothVar]], na.rm=TRUE)
  smoothVarOut <- seq(min(smoothVarRange), max(smoothVarRange), length.out = length(dif))
  out_df[[smoothVar]]<-smoothVarOut
  
  if(plot){
    p <- .plot_spline_diff(out_df)
    out <- list("data"=out_df, "plot" = p)
  } else{
    out <-out_df
  }
  return(out)
}  

#' ***********************************************************************************************
#' *************** `Plot difference in smooths` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for plotting spline_diff output
#' 
#' @keywords internal
#' @noRd

.plot_spline_diff<-function(df){
  
  x = colnames(df)[ncol(df)]
  nms<-colnames(df)
  nms<-as.numeric(sub("Q_", "", nms[grepl("^Q_", nms)]))
  cis<-numeric()
  i=1
  while(length(nms)>1){
    cis[i]<-max(nms)-min(nms)
    nms<-nms[-which.max(nms)]
    nms<-nms[-which.min(nms)]
    i=i+1
  }
  cis<-rev(cis)
  virPal <- viridis::viridis(n=length(cis), option="mako", direction=-1, begin =0.1)
  
  ggplot2::ggplot(df, aes(x=.data[[x]], y = .data[["mu"]]))+
    lapply(length(cis):1, function(i){
      ci = cis[i]
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[paste0("Q_",1-(1-((1-ci)/2)) )]], 
                                        ymax = .data[[paste0("Q_", 1-((1-ci)/2) )]]), fill=virPal[i], alpha=0.5)
    })+
    ggplot2::geom_line()+
    ggplot2::geom_hline(yintercept=0, linetype=5)+
    ggplot2::labs(y = paste0(df[1,"g1"], " - ", df[1, "g2"]))+
    pcv_theme()
}


