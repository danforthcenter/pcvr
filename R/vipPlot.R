#' Plot Variable Influence on Projection
#' 
#' @description This function is used to visualize variable influence on projection (vip) from a plsr model.
#' 
#' @param plsrObject Output from pcv.plsr
#' @param i An index from the plsrObject to use if the plsrObject contains models for several outcomes.
#'    Can be a name or a position. Defaults to 1.
#' @param mean Logical, should the mean be plotted (TRUE)
#'      or should the components be shown individually (FALSE, the default).
#' @param removePattern A pattern to remove to make the wavelength column into a numeric.
#' 
#' @import ggplot2
#' 
#' @keywords PLSR
#' @examples 
#' 
#' ## Not run:
#' 
#' if(FALSE){
#' file = "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest1.csv"
#' df1<-read.pcv(file, "wide", TRUE,
#'    multiValPattern = c("index_frequencies_index_ari", "index_frequencies_index_ci_rededge",
#'    "npq_hist_NPQ", "yii_hist_Fq'/Fm'", "yii_hist_Fv/Fm"))
#' colnames(df1)<-sub("index_frequencies_index_ndvi.", "ndvi_", colnames(df1))
#' # Note you will need pls installed to use pcv.plsr
#' x<-pcv.plsr(df=df1, resps = "area.pixels", spectra = grepl("^ndvi_", colnames(df1)))
#' plotVIP(x)
#' }
#' 
#' ## End(Not run)
#' 
#' @export

plotVIP<-function(plsrObject, i=1, mean=FALSE, removePattern = ".*_"){
  d<-plsrObject[[i]]$vip_df
  d$spectra<-as.numeric(sub(removePattern, "", d$wavelength))
  if(mean){
    d$meanVIP<-rowMeans(as.data.frame(d[, colnames(d)[grepl("VIP_[0-9]+?$", colnames(d))]]))
    p<-ggplot2::ggplot(d,ggplot2::aes(x=.data$spectra, y=.data$meanVIP))+
      ggplot2::geom_line()+
      ggplot2::geom_hline(yintercept=1, linetype=5)+
      ggplot2::labs(title=paste0(plsrObject[[i]]$model_performance$outcome), y="Mean VIP", x="Spectra")
  }else{
    cols<-c(which(grepl("VIP", colnames(d))), which(colnames(d)=="spectra"))
    d2<-as.data.frame(data.table::melt(data.table::as.data.table(d[,cols]), id.vars = "spectra", value.name = "VIP"))
    d2$component<-factor(as.numeric(sub("VIP_", "", d2$variable)))
    p<-ggplot2::ggplot(d2, ggplot2::aes(x=.data$spectra, y=.data$VIP, group=.data$component, color=.data$component))+
      ggplot2::geom_line()+
      ggplot2::geom_hline(yintercept=1, linetype=5)+
      #scale_color_viridis(option="plasma", discrete=TRUE, direction=1, begin=0.1, end=0.9)+
      ggplot2::guides(alpha="none", color = ggplot2::guide_legend(override.aes = list(linewidth=3)))+
      ggplot2::labs(title=paste0(plsrObject[[i]]$model_performance$outcome), y="VIP", x="Spectra", color="Component")+
      pcv_theme()+
      ggplot2::theme(plot.title = ggplot2::element_text(size=16))
  }
  return(p)
}
