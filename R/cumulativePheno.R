#' Reduce phenotypes in longitudinal data to cumulative sums of phenotypes.
#' 
#' 
#' 
#' @description Often in bellwether experiments we are curious about the effect of some treatment vs control. For certain routes in analysing the data this requires considering phenotypes as relative differences compared to a control.
#' 
#' @param df Dataframe to use, this is expected to be in wide format although in the future long may also be supported.
#' @param phenotypes A character vector of column names for the phenotypes that should be compared against control.
#' @param grouping A character vector of column names that identify groups in the data. Defaults to "barcode". These groups will be calibrated separately, with the exception of the group that identifies a control within the greater hierarchy.
#' @param control A column name for the variable to be used to select the control observations. If left NULL (the default) then this will be taken as the first string in the group argument.
#' @param controlGroup The level of the control variable to compare groups against.
#' @param method The method or methods to use, any of "proportion", "difference", or "zscore". These methods will be appended to the added column names ('phenotype_method').
#' @param traitCol Column with phenotype names, defaults to "trait". This should generally not need to be changed from the default.
#' @param valueCol Column with phenotype values, defaults to "value". This should generally not need to be changed from the default.
#' @return A dataframe with cumulative sum columns added for specified phenotypes
#' @keywords single-value-traits
#' @examples 
#' 
#' sv<-read.pcv("https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv", mode="wide", singleValueOnly = T, reader="fread")
#' sv$genotype = substr(sv$barcode, 3,5)
#' sv$genotype = ifelse(sv$genotype == "002", "B73",
#'               ifelse(sv$genotype == "003", "W605S",
#'               ifelse(sv$genotype == "004", "MM", "Mo17")))
#' sv$fertilizer = substr(sv$barcode, 8, 8)
#' sv$fertilizer = ifelse(sv$fertilizer == "A", "100",
#'               ifelse(sv$fertilizer == "B", "50", "0"))
#'               
#' sv<-bw.time(sv, plantingDelay = 0, phenotype="area.pixels", cutoff=10, timeCol="timestamp", group=c("barcode", "rotation"), plot=T)
#' sv<-bw.outliers(sv, phenotype="area.pixels", group = c("DAS", "genotype", "fertilizer"), plotgroup = c("barcode", "rotation"))
#' phenotypes <- c('area.pixels', 'convex_hull_area.pixels', 'convex_hull_vertices', 'ellipse_angle.degrees', 'ellipse_eccentricity', 'ellipse_major_axis.pixels', 'ellipse_minor_axis.pixels', 'height.pixels', 'hue_circular_mean.degrees', 'hue_circular_std.degrees', 'hue_median.degrees', 'longest_path.pixels', 'perimeter.pixels', 'solidity', 'width.pixels')
#' phenoForm<-paste0("cbind(", paste0(phenotypes, collapse=", "), ")")
#' groupForm<-"DAS+DAP+barcode+genotype+fertilizer"
#' form<-as.formula(paste0(phenoForm, "~", groupForm))
#' sv<-aggregate(form, data=sv, mean, na.rm=T)
#' pixels_per_cmsq <- 42.5^2   # pixel per cm^2
#' sv$area_cm2<-sv$area.pixels / pixels_per_cmsq
#' sv$height_cm <- sv$height.pixels/42.5
#' df = sv
#' phenotypes = c("area_cm2", "height_cm")
#' group = c("barcode")
#' timeCol = "DAS"
#' df<-cumulativePheno(df, phenotypes, group, timeCol)
#' 
#' 
#' sv_l<-read.pcv("https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv", mode="long", singleValueOnly = T, reader="fread")
#' sv_l$genotype = substr(sv$barcode, 3,5)
#' sv_l$genotype = ifelse(sv_l$genotype == "002", "B73",
#'               ifelse(sv_l$genotype == "003", "W605S",
#'               ifelse(sv_l$genotype == "004", "MM", "Mo17")))
#' sv_l$fertilizer = substr(sv_l$barcode, 8, 8)
#' sv_l$fertilizer = ifelse(sv_l$fertilizer == "A", "100",
#'               ifelse(sv_l$fertilizer == "B", "50", "0"))
#' sv_l<-bw.time(sv_l, plantingDelay = 0, phenotype="area", cutoff=10, timeCol="timestamp", group=c("barcode", "rotation"), wide=F)
#' sv_l<-cumulativePheno(sv_l, phenotypes=c("area", "height"), group=c("barcode", "rotation"), timeCol="DAS", wide=F)
#' 
#' 
#' @export
#' 
cumulativePheno<-function(df, phenotypes=NULL, group="barcode", timeCol="DAS", wide=T, traitCol="trait", valueCol="value"){
  
  if(length(group)>1){
    df$GROUP = as.character(interaction(df[, group]))
    group="GROUP"}
  
  if(!wide){ # note this is untested so far
    dat_sp<-split(df, df[[group]])
    out<-do.call(rbind, lapply(split(df, df[[group]]), function(d){
      newRows<-do.call(rbind, lapply(phenotypes, function(pheno){
        di<-d[d[[traitCol]] == pheno, ]
        di[[valueCol]]<-cumsum(di[[valueCol]])
        di[[traitCol]]<-paste0(pheno, "_csum")
        di
      }))
      rbind(d, newRows)
    }))
  } else{
    dat_sp<-split(df, df[[group]])
    out<-do.call(rbind, lapply(dat_sp, function(d){
      d<-d[sort(d[[timeCol]], index.return=T)$ix,]
      d2<-setNames(as.data.frame(do.call(cbind, lapply(phenotypes, function(pheno){cumsum( d[[pheno]] )}))), paste0(phenotypes, "_csum"))
      cbind(d,d2)
    }))
  }
  return(out)
}
