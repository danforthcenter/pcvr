#' Boxplot for plantCV data in wide format (each phenotype as a column) 
#' 
#' @param df Data frame to use. Can be wide or long format from read.pcv
#' @param x X axis variable(s) name (grouping).
#' If multiple names are provided then the variables are combined with \code{interaction}.
#' @param y Y axis variable name (phenotype)
#' @param fill Optional variable to use for boxplot color.
#' If left NULL then the x variable is used.
#' If this is set to a different variable than x then care should be taken with compare.
#' @param compare Groups to compare. By default this is set to FALSE,
#'  which corresponds to no testing. Other values of compare are passed to
#'  fixCompare to make t.test comparisons using ggpubr.
#'  In short, NULL will run all pairwise T tests, a single value of the X axis 
#'  variable will compare that level to all other levels of the X variable, 
#'  alternatively this can be a list as used by ggpubr: list(c("level1", "level2"),
#'   c("level1", "level3"))
#' @param trait If the data is in long format then `trait` is used to subset data for the phenotype Y.
#'  Defaults to "trait".
#' @param value If the data is in long format then `value` is used for numeric values for the phenotype Y.
#'  Defaults to "value".
#' @param method One of "t.test" or "wilcox.test".
#' @param showPoints Logical, should all points be shown? Defaults to FALSE.
#' @param ... Additional parameters passed to geom_boxplot
#' @keywords Bellwether, ggplot
#' @import ggplot2
#' @import ggpubr
#' @importFrom stats complete.cases
#' @return Returns a ggplot object.
#' @export
#' @examples 
#' 
#' ## Not run: 
#' 
#' sv<-read.pcv(
#' "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv",
#'  reader="fread")
#' sv$genotype = substr(sv$barcode, 3,5)
#' sv$genotype = ifelse(sv$genotype == "002", "B73",
#'               ifelse(sv$genotype == "003", "W605S",
#'               ifelse(sv$genotype == "004", "MM", "Mo17")))
#' sv$fertilizer = substr(sv$barcode, 8, 8)
#' sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
#'               ifelse(sv$fertilizer == "B", "50", "0"))
#'               
#' sv<-bw.time(sv, plantingDelay = 0, phenotype="area", cutoff=10,
#'  timeCol="timestamp", group=c("barcode", "rotation"), plot=TRUE)
#'  
#' pixels_per_cmsq <- 42.5^2   # pixel per cm^2
#' sv$area_cm2<-sv$area / pixels_per_cmsq
#' sv$height_cm <- sv$height/42.5
#' 
#' pcvBox(sv[sv$genotype=='MM' & sv$DAS==15, ],
#'   x="fertilizer", y="area_cm2", compare="0", showPoints = TRUE)
#' 
#' ## End(Not run) 

pcvBox<-function(df=df,x='treatment' , y='area', fill = NULL, compare=FALSE, trait="trait", value="value",method="t.test", showPoints=FALSE, ...){
  if(! (y %in% colnames(df)) && all(c(trait, value) %in% colnames(df))){
    df<-df[df[[trait]]==y,]
    ylab = y
    y = value
  } else {ylab=y}
  if(length(x)>1){
    df[[paste(x,collapse=".")]]<-interaction(df[,x])
    x<-paste(x,collapse=".")
  }else{
    df[[x]]<-as.factor(df[[x]])
  }
  
  df<-df[complete.cases(df[,c(x,y)]),]
  if(is.null(fill)){fill=x}
  boxLayer<-if(!showPoints){ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[fill]]), ...)
    }else{ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[fill]]), outlier.shape=NA, ...)}
  pointLayer<-if(showPoints){ggplot2::geom_jitter(size=0.5,width=0.1)
    }else{ggplot2::theme()}
  p<-ggplot2::ggplot(df, ggplot2::aes( x=.data[[x]], y=.data[[y]]))+ # I miss aes_string
    boxLayer+
    pointLayer+
    ggplot2::labs(x=x, y=ylab)+
    pcv_theme()+
    ggplot2::theme(axis.text.x.bottom = ggplot2::element_text(angle = 0, hjust = 1))
  if(!(is.logical(compare))){
    compare<-fixCompare(compare, df, x)
    p<-p+ggpubr::stat_compare_means(method=method, comparisons = compare , label = "p.format")
  }
  return(p)
}

