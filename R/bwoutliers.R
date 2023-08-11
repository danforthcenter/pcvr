#' Remove outliers from bellwether data using cook's distance
#' 
#' @param df Data frame to use. Can be in long or wide format.
#' @param phenotype Column to use to classify outliers.
#' @param naTo0 Logical, should NA values to changed to 0.
#' @param group  Grouping variables to find outliers as a character vector.
#' This is typically time  and design variables (DAS, genotype, treatment, etc).
#' These are used as predictors for `phenotype` in a generalized linear model.
#' @param cutoff Cutoff for something being an "outlier" expressed as a multiplier
#'  on the mean of Cooks Distance for this data. This defaults to 3 which tends to be
#'  a good value.
#' @param plotgroup Grouping variables for drawing plots if plot=TRUE.
#' Typically this is an identifier for images of a plant
#' over time and defaults to c('barcode',"rotation").
#' @param plot Logical, if TRUE then a list is returned with a ggplot and a dataframe.
#' @param wide Logical, is the data in wide format? Defaults to TRUE.
#' @param x Optional specification for x axis variable if plot is true.
#' If left NULL (the default) then the first element of `group` is used.
#' @param traitCol Column with phenotype names, defaults to "trait".
#' This should generally not need to be changed from the default.
#' @param valueCol Column with phenotype values, defaults to "value".
#' This should generally not need to be changed from the default.
#' @param idCol Column(s) that identify individuals over time.
#' Defaults to plotGroup.
#' @keywords Bellwether, ggplot
#' @import ggplot2
#' @importFrom stats complete.cases cooks.distance glm as.formula
#' @examples 
#' 
#' ## Not run:
#' 
#' sv<-read.pcv(
#' "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv",
#'  mode="wide", singleValueOnly =TRUE, reader="fread")
#' sv$genotype = substr(sv$barcode, 3,5)
#' sv$genotype = ifelse(sv$genotype == "002", "B73",
#'                      ifelse(sv$genotype == "003", "W605S",
#'                            ifelse(sv$genotype == "004", "MM", "Mo17")))
#' sv$fertilizer = substr(sv$barcode, 8, 8)
#' sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
#'                    ifelse(sv$fertilizer == "B", "50", "0"))
#' sv<-bw.time(sv, plantingDelay = 0, phenotype="area.pixels", cutoff=10, timeCol="timestamp",
#'  group=c("barcode", "rotation"), plot=FALSE)
#' sv<-bw.outliers(df = sv, phenotype="area.pixels", naTo0 =FALSE, 
#'  group = c("DAS", "genotype", "fertilizer"),
#'  plotgroup=c('barcode',"rotation"), plot=TRUE)
#' 
#' 
#' svl<-read.pcv(
#' "https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv",
#'  mode="long", singleValueOnly =TRUE, reader="fread")
#' svl$genotype = substr(svl$barcode, 3,5)
#' svl$genotype = ifelse(svl$genotype == "002", "B73",
#'                      ifelse(svl$genotype == "003", "W605S",
#'                            ifelse(svl$genotype == "004", "MM", "Mo17")))
#' svl$fertilizer = substr(svl$barcode, 8, 8)
#' svl$fertilizer = ifelse(svl$fertilizer == "A", "100",
#'                    ifelse(svl$fertilizer == "B", "50", "0"))
#' svl<-bw.time(svl, plantingDelay = 0, phenotype="area", cutoff=10, timeCol="timestamp",
#'  group=c("barcode", "rotation"), plot=FALSE,wide=FALSE)
#' 
#' svl<-bw.outliers(df = svl, phenotype="area", naTo0 =FALSE,
#'   group = c("DAS", "genotype", "fertilizer"),
#'   plotgroup=c('barcode',"rotation"), plot=TRUE, wide=FALSE)
#' 
#' ## End(Not run)
#' 
#' @return The input dataframe with outliers removed.
#' @export


bw.outliers<-function(df = NULL,
                      phenotype,
                      naTo0 =FALSE,
                      group = c(),
                      cutoff = 3,
                      plotgroup=c('barcode',"rotation"),
                      plot=TRUE,  wide =TRUE, x=NULL, traitCol="trait", valueCol="value", idCol=NULL){
  # df = svl; phenotype="area"; naTo0 = F; group = c("DAS", "genotype", "fertilizer"); plotgroup=c("barcode","rotation"); plot=FALSE
  # wide = F ; traitCol="trait"; valueCol="value"; idCol=c("barcode","rotation")
  if(is.null(idCol)){idCol = plotgroup}
  outlierMethod = "cooks"
  if(wide){
    if(naTo0){
      df[[phenotype]][is.na(df[[phenotype]])]<-0
    }
    df<-df[complete.cases(df[,c(phenotype,group)]), ]
    outlierForm<-paste("as.numeric(",phenotype,")~", paste(paste0("as.factor(",group,")"),collapse=":"))
    if(outlierMethod=="cooks"){
      cooksd <- cooks.distance(glm(data=df, as.formula(outlierForm)))
      # summary(cooksd)
      outlierCutoff<-cutoff*mean(cooksd, na.rm=TRUE)
      cooksd[is.na(cooksd)] <- outlierCutoff - 0.1 # keeping NAs by assigning a value below cutoff.
      cooksd_df<-data.frame("outlier" = cooksd)
      df<-cbind(df, cooksd_df) 
      pctRm<-paste0(100*(1-round(nrow(df[df$outlier < outlierCutoff, ]) / nrow(df), 5)), "% removed as outliers using Cook's Distance")
    }
  } else{ # long data version
    if(naTo0){
      df[df[[traitCol]]==phenotype, valueCol][is.na(df[df[[traitCol]]==phenotype, valueCol])]<-0
    }
    subdf<-df[complete.cases( df[df[[traitCol]]==phenotype, c(valueCol, traitCol, group)] ) & df[[traitCol]]==phenotype , ]
    outlierForm<-paste("as.numeric(",valueCol,")~", paste(paste0("as.factor(",group,")"),collapse=":"))
    if(outlierMethod=="cooks"){
      cooksd <- cooks.distance(glm(data=subdf, as.formula(outlierForm)))
      outlierCutoff<-3*mean(cooksd, na.rm=TRUE)
      cooksd[is.na(cooksd)] <- outlierCutoff - 0.1 # keeping NAs by assigning a value below cutoff.
      cooksd_df<-data.frame("outlier" = cooksd)
      subdf<-cbind(subdf, cooksd_df)
      pctRm<-paste0(100*(1-round(nrow(subdf[subdf$outlier < outlierCutoff, ]) / nrow(subdf), 5)),
                    "% removed as outliers using Cook's Distance")
      ids<-unique(interaction(subdf[,c(idCol)]))
      df<-merge(df, subdf, all.x = TRUE)
    }
  }
  out<-df[df$outlier < outlierCutoff, -which(colnames(df)=="outlier")]
  if(plot & wide){
    df$grouping<-interaction(df[,plotgroup])
    out_plotData<-df[df$outlier < outlierCutoff, ]
    rmdf_plotData<-df[df$outlier >= outlierCutoff, ]
    if(is.null(x)){x = group[1]}
    p<-ggplot2::ggplot()+
      ggplot2::geom_line(data=rmdf_plotData, ggplot2::aes(x=.data[[x]], y=.data[[phenotype]], group=.data[["grouping"]]), linewidth=0.15, color="red")+
      ggplot2::geom_line(data=out_plotData, ggplot2::aes(x=.data[[x]], y=.data[[phenotype]], group=.data[["grouping"]]),linewidth=0.25 )+
      ggplot2::labs(title=pctRm)+
      pcv_theme()
    print(p)
  } else if(plot & !wide){
    plotdf<-df[df[[traitCol]]==phenotype,]
    plotdf$grouping<-interaction(plotdf[,plotgroup])
    out_plotData<-plotdf[plotdf$outlier < outlierCutoff, ]
    rmdf_plotData<-plotdf[plotdf$outlier >= outlierCutoff, ]
    if(is.null(x)){x = group[1]}
    p<-ggplot2::ggplot()+
      ggplot2::geom_line(data=rmdf_plotData, ggplot2::aes(x=.data[[x]], y=.data[[valueCol]], group=.data[["grouping"]]), linewidth=0.15, color="red")+
      ggplot2::geom_line(data=out_plotData, ggplot2::aes(x=.data[[x]], y=.data[[valueCol]], group=.data[["grouping"]]),linewidth=0.25 )+
      ggplot2::labs(title=pctRm)+
      pcv_theme()
    print(p)
    
  }
  return(out)
}
