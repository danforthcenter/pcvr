#' Remove outliers from bellwether data using cook's distance
#' 
#' @param df Data frame to use. Can be in long or wide format.
#' @param phenotype Column to use to classify outliers. If this is length > 1 then
#' it is taken as the multi-value traits to use. See examples.
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
#' @param x Optional specification for x axis variable if plot is true.
#' If left NULL (the default) then the first element of `group` is used.
#' @param traitCol Column with phenotype names, defaults to "trait".
#' This should generally not need to be changed from the default.
#'    If this and valueCol are present in colnames(df) then the data
#'    is assumed to be in long format.
#' @param valueCol Column with phenotype values, defaults to "value".
#' This should generally not need to be changed from the default.
#' @param labelCol Column with phenotype labels for long data, defaults to "label".
#' This should generally not need to be changed from the default.
#' @param idCol Column(s) that identify individuals over time.
#' Defaults to plotGroup.
#' @param ncp Optionally specify the number of principle components to be used for MV data outlier detection with cooks distance.
#' If left NULL (the default) then 3 will be used.
#' @keywords Bellwether, ggplot, outliers
#' @import ggplot2
#' @import data.table
#' @importFrom stats complete.cases cooks.distance glm as.formula lm
#' @examples 
#' 
#' ## Not run:
#' 
#' sv<-read.pcv(
#' "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv",
#'  reader="fread")
#' sv$genotype = substr(sv$barcode, 3,5)
#' sv$genotype = ifelse(sv$genotype == "002", "B73",
#'                      ifelse(sv$genotype == "003", "W605S",
#'                            ifelse(sv$genotype == "004", "MM", "Mo17")))
#' sv$fertilizer = substr(sv$barcode, 8, 8)
#' sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
#'                    ifelse(sv$fertilizer == "B", "50", "0"))
#' sv<-bw.time(sv, plantingDelay = 0, phenotype="area_pixels", cutoff=10, timeCol="timestamp",
#'  group=c("barcode", "rotation"), plot=FALSE)
#' sv<-bw.outliers(df = sv, phenotype="area_pixels", naTo0 =FALSE, 
#'  group = c("DAS", "genotype", "fertilizer"),
#'  plotgroup=c('barcode',"rotation"), plot=TRUE)
#' 
#' if(FALSE){
#' svl<-read.pcv(
#' "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcv4-single-value-traits.csv",
#'  mode="long", reader="fread")
#' svl$genotype = substr(svl$barcode, 3,5)
#' svl$genotype = ifelse(svl$genotype == "002", "B73",
#'                      ifelse(svl$genotype == "003", "W605S",
#'                            ifelse(svl$genotype == "004", "MM", "Mo17")))
#' svl$fertilizer = substr(svl$barcode, 8, 8)
#' svl$fertilizer = ifelse(svl$fertilizer == "A", "100",
#'                    ifelse(svl$fertilizer == "B", "50", "0"))
#' svl<-bw.time(svl, plantingDelay = 0, phenotype="area_pixels", cutoff=10, timeCol="timestamp",
#'  group=c("barcode", "rotation"), plot=FALSE)
#' 
#' svl<-bw.outliers(df = svl, phenotype="area_pixels", naTo0 =FALSE,
#'   group = c("DAS", "genotype", "fertilizer"),
#'   plotgroup=c('barcode',"rotation"), plot=TRUE)
#'   
#' mvw<-read.pcv(paste0("https://media.githubusercontent.com/media/joshqsumner/",
#'   "pcvrTestData/main/pcv4-multi-value-traits.csv"), mode="wide")
#'  mvw$genotype = substr(mvw$barcode, 3,5)
#' mvw$genotype = ifelse(mvw$genotype == "002", "B73",
#'                      ifelse(mvw$genotype == "003", "W605S",
#'                             ifelse(mvw$genotype == "004", "MM", "Mo17")))
#' mvw$fertilizer = substr(mvw$barcode, 8, 8)
#' mvw$fertilizer = ifelse(mvw$fertilizer == "A", "100",
#'                        ifelse(mvw$fertilizer == "B", "50", "0"))
#' mvw<-bw.time(mvw,timeCol="timestamp", group="barcode", plot = FALSE)
#' 
#' phenotypes = which(grepl("hue_freq", colnames(mvw)))
#' 
#' mvw2 <- bw.outliers(df = mvw, phenotype = phenotypes, naTo0 = FALSE,
#'     group = c("DAS", "genotype", "fertilizer"), cutoff = 3, plotgroup=c("barcode", "rotation"))
#' 
#' 
#' mvl<-read.pcv(paste0("https://media.githubusercontent.com/media/joshqsumner/",
#'   "pcvrTestData/main/pcv4-multi-value-traits.csv"), mode="long")
#'  mvl$genotype = substr(mvl$barcode, 3,5)
#' mvl$genotype = ifelse(mvl$genotype == "002", "B73",
#'                      ifelse(mvl$genotype == "003", "W605S",
#'                             ifelse(mvl$genotype == "004", "MM", "Mo17")))
#' mvl$fertilizer = substr(mvl$barcode, 8, 8)
#' mvl$fertilizer = ifelse(mvl$fertilizer == "A", "100",
#'                        ifelse(mvl$fertilizer == "B", "50", "0"))
#' mvl<-bw.time(mvl,timeCol="timestamp", group="barcode", plot = FALSE)
#' 
#' mvl2 <- bw.outliers(df = mvl, phenotype = "hue_frequencies", naTo0 = FALSE,
#'     group = c("DAS", "genotype", "fertilizer"), cutoff = 3, plotgroup=c("barcode", "rotation"))
#' }
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
                      plot=TRUE, x=NULL, traitCol="trait", valueCol="value", labelCol="label", idCol=NULL, ncp = NULL){
  # df = svl; phenotype="area_pixels"; naTo0 = F; group = c("DAS", "genotype", "fertilizer"); plotgroup=c("barcode","rotation"); plot=FALSE
  # traitCol="trait"; valueCol="value"; idCol=c("barcode","rotation")
  if(all(c(traitCol, valueCol) %in% colnames(df))){
    wide=FALSE
  } else{ wide = TRUE}
  if(is.null(phenotype)){stop("A phenotype must be provided")}
  
  if( (wide & length(phenotype)>1) | (!wide && length(unique(interaction(df[df[[traitCol]]==phenotype, 
                                                                            colnames(df) %in% c(traitCol, labelCol)])))>1 ) ){
    mv = TRUE
  } else { mv = FALSE }
  
  
  if(is.null(idCol)){idCol = plotgroup}
  outlierMethod = "cooks"
  
  
  if(wide & !mv){ # wide data single value
    res<-.wide_sv_cooks_bw.outliers(df, naTo0, phenotype, group, cutoff)
    df <- res[["data"]]
    pctRm <- res[["pctRm"]]
    
  } else if (!wide & !mv){ # long data single value
    res<-.long_sv_cooks_bw.outliers(df, naTo0, phenotype, group, cutoff, traitCol, valueCol, labelCol, idCol)
    df <- res[["data"]]
    pctRm <- res[["pctRm"]]
    
  } else if(wide & mv){ # wide multi value data
    res<-.wide_mv_cooks_bw.outliers(df, naTo0, phenotype, group, cutoff, ncp)
    df <- res[["data"]]
    pctRm <- res[["pctRm"]]
    
  } else if(!wide & mv){ # long multi value data
    res<-.long_mv_cooks_bw.outliers(df, naTo0, phenotype, group, cutoff, ncp, traitCol, valueCol, labelCol, idCol)
    df <- res[["data"]]
    pctRm <- res[["pctRm"]]
    
  }
  
  
  out<-df[which(!df$outlier), -which(grepl("outlier", colnames(df)))]
  
  if(plot & wide & !mv){
    df$grouping<-interaction(df[,plotgroup])
    out_plotData<-df[!df$outlier, ]
    rmdf_plotData<-df[df$outlier, ]
    if(is.null(x)){x = group[1]}
    p<-ggplot2::ggplot()+
      ggplot2::geom_line(data=rmdf_plotData, ggplot2::aes(x=.data[[x]], y=.data[[phenotype]], group=.data[["grouping"]]), linewidth=0.15, color="red")+
      ggplot2::geom_line(data=out_plotData, ggplot2::aes(x=.data[[x]], y=.data[[phenotype]], group=.data[["grouping"]]),linewidth=0.25 )+
      ggplot2::labs(title=pctRm)+
      pcv_theme()
    print(p)
  } else if(plot & !wide & !mv){
    plotdf<-df[df[[traitCol]]==phenotype,]
    plotdf$grouping<-interaction(plotdf[,plotgroup])
    out_plotData<-plotdf[!plotdf$outlier, ]
    rmdf_plotData<-plotdf[plotdf$outlier, ]
    if(is.null(x)){x = group[1]}
    p<-ggplot2::ggplot()+
      ggplot2::geom_line(data=rmdf_plotData, ggplot2::aes(x=.data[[x]], y=.data[[valueCol]], group=.data[["grouping"]]), linewidth=0.15, color="red")+
      ggplot2::geom_line(data=out_plotData, ggplot2::aes(x=.data[[x]], y=.data[[valueCol]], group=.data[["grouping"]]),linewidth=0.25 )+
      ggplot2::labs(title=pctRm)+
      pcv_theme()
    print(p)
    
  } else if(plot & wide & mv){
    #* lengthen data, make bar plot?
    
    plotdf <- suppressWarnings(as.data.frame(data.table::melt(data.table::as.data.table(df), measure.vars = phenotype,
                                                           variable.name = traitCol, value.name = valueCol)))
    
    plotdf$bin <- as.numeric(regmatches(plotdf$trait,regexpr("[0-9]+",plotdf$trait)))
    
    plotdf$grouping<-interaction(plotdf[,plotgroup])
    out_plotData<-plotdf[!plotdf$outlier, ]
    rmdf_plotData<-plotdf[plotdf$outlier, ]
    
    if(is.null(x)){x = group[1]}
    
    p<-ggplot2::ggplot()+
      ggplot2::geom_col(data = rmdf_plotData, ggplot2::aes(x = .data[['bin']], y=.data[[valueCol]] ), position="identity",
                        fill="red", alpha=0.25)+
      ggplot2::geom_col(data = out_plotData, ggplot2::aes(x = .data[['bin']], y=.data[[valueCol]] ), position="identity",
                        alpha=0.25)+
      ggplot2::labs(title=pctRm)+
      pcv_theme()
    print(p)
    
  } else if(plot & !wide & mv){
    
    plotdf <- df
    plotdf$grouping<-interaction(plotdf[,plotgroup])
    out_plotData<-plotdf[!plotdf$outlier, ]
    rmdf_plotData<-plotdf[plotdf$outlier, ]
    
    p<-ggplot2::ggplot()+
      ggplot2::geom_col(data = rmdf_plotData, ggplot2::aes(x = .data[[traitCol]], y=.data[[valueCol]] ), position="identity",
                        fill="red", alpha=0.25)+
      ggplot2::geom_col(data = out_plotData, ggplot2::aes(x = .data[[traitCol]], y=.data[[valueCol]] ), position="identity",
                        alpha=0.25)+
      ggplot2::labs(title=pctRm)+
      pcv_theme()
    print(p)
  }
  
  
  
  return(out)
}


#' ***********************************************************************************************
#' *************** `wide SV cooks distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in wide SV data.
#' 
#' @keywords internal
#' @noRd

.wide_sv_cooks_bw.outliers <- function(df, naTo0, phenotype, group, cutoff){
  if(naTo0){
    df[[phenotype]][is.na(df[[phenotype]])]<-0
  }
  df<-df[complete.cases(df[,c(phenotype,group)]), ]
  outlierForm<-paste("as.numeric(",phenotype,")~", paste(paste0("as.factor(",group,")"),collapse=":"))
    cooksd <- cooks.distance(glm(data=df, as.formula(outlierForm)))
    # summary(cooksd)
    outlierCutoff<-cutoff*mean(cooksd, na.rm=TRUE)
    cooksd[is.na(cooksd)] <- outlierCutoff - 0.1 # keeping NAs by assigning a value below cutoff.
    cooksd_df<-data.frame("outlier" = cooksd)
    df<-cbind(df, cooksd_df)
    df$outlier <- df$outlier > outlierCutoff
    pctRm<-paste0(100*(round(nrow(df[df$outlier, ]) / nrow(df), 5)), "% removed as outliers using Cook's Distance")
    
    return(list('data'=df, 'pctRm'=pctRm))
}

#' ***********************************************************************************************
#' *************** `long SV cooks distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in wide SV data.
#' 
#' @keywords internal
#' @noRd

.long_sv_cooks_bw.outliers <- function(df, naTo0, phenotype, group, cutoff, traitCol, valueCol, labelCol, idCol){
  if(naTo0){
    df[df[[traitCol]]==phenotype, valueCol][is.na(df[df[[traitCol]]==phenotype, valueCol])]<-0
  }
  subdf<-df[complete.cases( df[df[[traitCol]]==phenotype, c(valueCol, traitCol, group)] ) & df[[traitCol]]==phenotype , ]
  outlierForm<-paste("as.numeric(",valueCol,")~", paste(paste0("as.factor(",group,")"),collapse=":"))
    cooksd <- cooks.distance(glm(data=subdf, as.formula(outlierForm)))
    outlierCutoff<-cutoff*mean(cooksd, na.rm=TRUE)
    cooksd[is.na(cooksd)] <- outlierCutoff - 0.1 # keeping NAs by assigning a value below cutoff.
    cooksd_df<-data.frame("outlier" = cooksd)
    subdf<-cbind(subdf, cooksd_df)
    subdf<-subdf[, c(group, idCol, "outlier")]
    subdf<-subdf[!duplicated(subdf[,c(group, idCol)]),] # if there are multiple images per day this will change data.
    subdf$outlier <- subdf$outlier > outlierCutoff
    pctRm<-paste0(100*(round(nrow(subdf[subdf$outlier, ]) / nrow(subdf), 5)),
                  "% removed as outliers using Cook's Distance")
    #* take IDs using plotgroup and label all phenotype rows
    df<-merge(df, subdf, all.x = TRUE)
  return(list('data'=df, 'pctRm'=pctRm))
}

#' ***********************************************************************************************
#' *************** `wide MV cooks distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in wide MV data.
#' 
#' @keywords internal
#' @noRd

.wide_mv_cooks_bw.outliers <- function(df, naTo0, phenotype, group, cutoff, ncp){
  
  if(naTo0){
    df[,phenotype][is.na(df[,phenotype])]<-0
  }
  
  phenos_df <- df[, phenotype]
  if(is.null(ncp)){
    use_ncp <- min(min(dim(phenos_df))-1, 3)
  } else{
    use_ncp <- ncp
  }
  pca<-FactoMineR::PCA(phenos_df, ncp=use_ncp, graph=FALSE)
  pc1Var<-round(pca$eig[1,2], 3)
  pc2Var<-round(pca$eig[2,2], 3)
  coords<-as.data.frame(pca$ind)
  coords<-coords[,grepl("coord", colnames(coords))]
  colnames(coords)<-gsub("coord.Dim.", "pc", colnames(coords))
  pca_cols <- colnames(coords)
  df <- cbind(df, coords)
  
  if(is.null(ncp)){
    message(paste0("Using ", use_ncp, " PCs comprising ", round(pca$eig[use_ncp, 3], 3),
                   "% of variation"))
  }
  
  df<-df[complete.cases(df[,c(pca_cols,group)]), ]
  
  outlierForm<-paste("cbind(",paste0("pc",1:use_ncp, collapse=","),")~",
                     paste(paste0("as.factor(",group,")"),collapse=":"))
  cooksd <- cooks.distance(lm(data=df, as.formula(outlierForm)))
  
  df <- df[, -which(colnames(df) %in% c(paste0("pc",1:use_ncp)))]
  
  if(length(cutoff)==1){
    cutoff <- rep(cutoff, use_ncp)
  }
  outlierCutoffs<-cutoff*colMeans(cooksd, na.rm=TRUE)
  
  outlierMatrix <- do.call(cbind, lapply(1:ncol(cooksd), function(i){
    cooks_vec <- cooksd[,i]
    cooks_vec[is.na(cooks_vec)] <- outlierCutoffs[i] - 0.1
    setNames(data.frame(cooks_vec > outlierCutoffs[i]), paste0("outlier_",i))
  }))
  outlierMatrix$outlier <- unlist(lapply(1:nrow(outlierMatrix), function(i){
    any(outlierMatrix[i,]) # could be a more nuanced rule
  }))
  
  df<-cbind(df, outlierMatrix)
  pctRm<-paste0(100*(round(nrow(df[df$outlier, ]) / nrow(df), 5)), "% removed as outliers using Cook's Distance\non first ", use_ncp, " Principle Components.")
  return(list('data'=df, 'pctRm'=pctRm))
}


#' ***********************************************************************************************
#' *************** `long MV cooks distance` ****************************************
#' ***********************************************************************************************
#' @description
#' Internal function for outlier detection in long MV data.
#' 
#' @keywords internal
#' @noRd

.long_mv_cooks_bw.outliers <- function(df, naTo0, phenotype, group, cutoff, ncp, traitCol, valueCol, labelCol, idCol){
  #* widen data
  dcast_form <- as.formula(paste0("... ~ ", traitCol, "+", labelCol))
  dfw <- as.data.frame(data.table::dcast(data.table::as.data.table(df[df[[traitCol]]==phenotype,]),
                                         dcast_form,  value.var = valueCol, sep="."))
  phenotypew <- which(grepl(phenotype, colnames(dfw)))
  #* call .wide method on dfw
  wide_res <- .wide_mv_cooks_bw.outliers(dfw, naTo0, phenotypew, group, cutoff, ncp)
  pctRm <- wide_res[["wide_res"]]
  sub_df <- wide_res[["data"]]
  #* label long data based on .wide output
  kept <- unique(as.character(interaction(sub_df[!sub_df$outlier, c(group, idCol)])))
  df_ids <- as.character(interaction(df[,c(group, idCol)]))
  df$outlier <- !df_ids %in% kept
  
  return(list('data'=df, 'pctRm'=pctRm))
}






