#' Boxplots for bellwether data
#' 
#' @param df Data frame to use. Should be in wide format similar to output from read.pcv.bw
#' @param x X axis variable name
#' @param y Y axis variable name
#' @param time Value of timeCol to use in subsetting data. If NULL (default) then timeCol is coerced to numeric and the max is used.
#' @param timeCol Column of time data to use. Defaults to DAS corresponding to Days after Start from read.pcv.bw
#' @param compare Groups to compare. By default this is set to FALSE, which corresponds to no testing. Other values of compare are passed to fixCompare to make t.test comparisons using ggpubr. In short, NULL will run all pairwise T tests, a single value of the X axis variable will compare that level to all other levels of the X variable, alternatively this can be a list as used by ggpubr: list(c("level1", "level2"), c("level1", "level3"))
#' @param fill Optional variable to use for boxplot color. If set to T then the x variable is used. If this is set to a different variable than x then care should be taken with compare.
#' @param ... Additional parameters passed to geom_boxplot
#' @keywords Bellwether, ggplot
#' @import ggplot2
#' @import ggpubr
#' @return Returns a ggplot object
#' @export
#' @examples 
#' ## Not run: 
#' df <- read.pcv.bw("pcvTestData.csv")
#' bellwetherBox(df, x="genotype", y = "area.pixels", time = 19)
#' ## End(Not run:)



bellwetherBox<-function(df=NULL,x=NULL , y=NULL, time=NULL, timeCol="DAS", compare=F, fill = NULL,  ...){
  df<-df[complete.cases(df[,c(x,y,timeCol)]),]
  if(is.null(time)){time=max(as.numeric(df[[timeCol]]), na.rm=T)}
  df<-df[df[[timeCol]] %in% time,]
  df[[x]]<-as.factor(df[[x]])
  appender <- function(string, prefix = timeCol) paste(timeCol, string)
  if(is.logical(fill) && fill){fill=x}
  boxLayer<-if(is.null(fill)){ggplot2::geom_boxplot(show.legend=F, ...)}else{ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[fill]]), ...)}
  
  p<-ggplot2::ggplot(df, ggplot2::aes( x=.data[[x]], y=.data[[y]]))+
    ggplot2::facet_wrap(~.data[[timeCol]],
               nrow = ceiling(sqrt(length(time))),
               labeller = ggplot2::as_labeller(appender) )+
    boxLayer+ 
    pcv_theme()
  if(is.logical(compare) && compare){compare=NULL}
  if(!(is.logical(compare))){
    compare<-fixCompare(compare, df, x)
    p<-p+ggpubr::stat_compare_means(method="t.test", comparisons = compare , label = "p.format")
  }
  return(p)
}
