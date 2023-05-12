#' Remove outliers from bellwether data using cook's distance
#' 
#' @param df Data frame to use. Should be in the format of output from read.pcv.bw
#' @param phenotype Column to use to classify outliers.
#' @param naTo0 Logical, should NA values to changed to 0.
#' @param group  Grouping variables to find outliers as a character vector. This is typically time  and design variables (DAS, genotype, treatment, etc). These are used as predictors for `phenotype` in a generalized linear model.
#' @param plotGroup Grouping variables for drawing plots if plot=T. Typically this is an identifier for images of a plant over time.
#' @param plot Logical, if TRUE then a list is returned with a ggplot and a dataframe.
#' @param x Optional specification for x axis variable if plot is true. If left NULL (the default) then the first element of `group` is used.
#' @keywords Bellwether, ggplot
#' @import ggplot2
#' @examples 
#' bw<-read.pcv.bw( file="https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/bwTestPhenos.csv", snapshotFile="https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/bwTestSnapshot.csv", designFile="https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/bwTestDesign.csv",metaCol="meta",metaForm="vis_view_angle_zoom_horizontal_gain_exposure_v_new_n_rep",joinSnapshot="id",conversions = list(area=13.2*3.7/46856) )
#' bw<-bw.outliers(df = bw, phenotype="area_adj", naTo0 = F, group = c("DAS", "Genotype"), plotgroup=c('Barcodes',"angle"), plot=T)

bw.outliers<-function(df = NULL,
                      phenotype="area_adj",
                      naTo0 = F,
                      group = c("DAS", "Genotype"),
                      plotgroup=c('Barcodes',"angle"),
                      plot=T,  x=NULL){
  outlierMethod = "cooks"
  if(naTo0){
    df[[phenotype]][is.na(df[[phenotype]])]<-0
  }
  df<-df[complete.cases(df[,c(phenotype,group)]), ]
  outlierForm<-paste("as.numeric(",phenotype,")~", paste(paste0("as.factor(",group,")"),collapse=":"))
  if(outlierMethod=="cooks"){
    cooksd <- cooks.distance(glm(data=df, as.formula(outlierForm)))
    summary(cooksd)
    outlierCutoff<-3*mean(cooksd, na.rm=T)
    cooksd_df<-data.frame("outlier" = cooksd)
    df<-cbind(df, cooksd_df) 
    pctRm<-paste0(100*(1-round(nrow(df[df$outlier < outlierCutoff, ]) / nrow(df), 5)), "% removed as outliers using Cook's Distance")
  }
  if("Genotype" %in% group){ df<-df[df$Genotype != "Empty", ] }
  out<-df[df$outlier < outlierCutoff, -which(colnames(df)=="outlier")]
  rmdf<-df[df$outlier >= outlierCutoff, -which(colnames(df)=="outlier")]
  df$grouping<-interaction(df[,plotgroup])
  out_plotData<-df[df$outlier < outlierCutoff, ]
  rmdf_plotData<-df[df$outlier >= outlierCutoff, ]
  if(plot){
    if(is.null(x)){x = group[1]}
    p<-ggplot2::ggplot()+
      ggplot2::geom_line(data=rmdf_plotData, aes(x=.data[[x]], y=.data[[phenotype]], group=.data[["grouping"]]), linewidth=0.15, color="red")+
      ggplot2::geom_line(data=out_plotData, aes(x=.data[[x]], y=.data[[phenotype]], group=.data[["grouping"]]),linewidth=0.25 )+
      ggplot2::labs(title=pctRm)+
      pcv_theme()
    print(p)
  }
  return(out)
}
