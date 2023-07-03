#' Calculate relative tolerance of some phenotype(s) relative to control
#' 
#' @description Often in bellwether experiments we are curious about the effect of some treatment vs control. For certain routes in analysing the data this requires considering phenotypes as relative differences compared to a control.
#' 
#' @param df Dataframe to use, this is expected to be in wide format although in the future long may also be supported.
#' @param phenotypes A character vector of column names for the phenotypes that should be compared against control.
#' @param grouping A character vector of column names that identify groups in the data. These groups will be calibrated separately, with the exception of the group that identifies a control within the greater hierarchy. Note that for levels of grouping where the control group does not exist the output will be NA.
#' @param control A column name for the variable to be used to select the control observations. If left NULL (the default) then this will be taken as the first string in the group argument.
#' @param controlGroup The level of the control variable to compare groups against.
#' @param method The method or methods to use, any of "proportion", "difference", or "zscore". These methods will be appended to the added column names ('phenotype_method').
#' @param naTo0 Logical, should NA and Inf values be replaced with 0? This is useful if output are going to be used in a cumulative step, but otherwise should be left False
#' @param wide Logical, is the input data in wide format? Defaults to TRUE.
#' @param traitCol Column with phenotype names, defaults to "trait". This should generally not need to be changed from the default.
#' @param valueCol Column with phenotype values, defaults to "value". This should generally not need to be changed from the default.
#' @return A dataframe with relative tolerance columns added.
#' @keywords single-value-trait
#' @examples 
#' 
#' ## Not run:
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
#' df=sv
#' phenotypes=phenotypes = c("area_cm2", "height_cm")
#' grouping = c("fertilizer", "genotype", "DAS") # need to include genotype so that the interaction effect is included
#' controlGroup = "100" 
#' control = "fertilizer"
#' method = c("proportion", "difference", "zscore")
#' naTo0=T
#' rt <-relativeTolerance(df, phenotypes, grouping, control, controlGroup, method, naTo0)
#' sum(is.na(rt$))
#' 
#' ## End(Not run)
#' 
#' @export
#' 

relativeTolerance<-function(df, phenotypes=NULL, grouping=NULL, control=NULL,
                            controlGroup = NULL, method = c("proportion", "difference", "zscore"), naTo0=F,
                            wide=T, traitCol="trait", valueCol="value"){

if(is.null(grouping)){grouping=control}
if(is.null(control)){control=grouping[1]}
if(is.null(controlGroup)){controlGroup = unique(df[[control]])[1]}
method = match.arg(method, c("proportion", "difference", "zscore"), several.ok = T)

if(control %in% grouping){
  group_no_control<-grouping[grouping!=control]
} else{
  group_no_control=grouping
}

group_no_control_factor<-interaction(df[,group_no_control])
datsp<-split(x=df, f=group_no_control_factor)

if(wide){
df2<-do.call(rbind, lapply(1:length(datsp), function(i){
  d = datsp[[i]]
  d2<-do.call(cbind, lapply(phenotypes, function(pheno){
    control_mean_value = mean( d[ d[[control]]==controlGroup, pheno ], na.rm=T)
   if(is.na(control_mean_value)){ values=NA
    }else{
      values<-c()
      values<-unlist(lapply(method, function(mthd){
        if(mthd=="proportion"){
          values = c(values, d[, pheno] / control_mean_value)
        }
        if(mthd == "difference"){
          values = c(values, d[, pheno] - control_mean_value)
        }
        if(mthd =="zscore"){
          values = c(values, (d[, pheno] - control_mean_value) / sd(d[ d[[control]]==controlGroup, pheno ], na.rm=T))
        }
        return(values)
      }))
    }
    if(naTo0){
      values[is.na(values)|is.infinite(values)]<-0
    }
    return( setNames(as.data.frame(matrix(values, ncol=length(method), nrow = nrow(d))), paste0(pheno, "_",method)) )
  }))
  d2<-cbind(d,d2)
  return(d2)
}))
} else { # long version
  
  df2<-do.call(rbind, lapply(1:length(datsp), function(i){
    d = datsp[[i]]
    d2<-do.call(cbind, lapply(phenotypes, function(pheno){
      control_mean_value = mean( d[ d[[control]]==controlGroup & d[[traitCol]]==pheno, valueCol ], na.rm=T)
      inputData <- d[d[[traitCol]]==pheno, ]
      if(is.na(control_mean_value)){ return(inputData)
      }else{
        newRows<-do.call(rbind, lapply(method, function(mthd){
          if(mthd=="proportion"){
            inputData[[valueCol]] = d[d[[traitCol]]==pheno, valueCol] / control_mean_value
            inputData[[traitCol]] = paste0(pheno,"_proportion")
          }
          if(mthd == "difference"){
            inputData[[valueCol]] =  d[d[[traitCol]]==pheno, valueCol] - control_mean_value
            inputData[[traitCol]] = paste0(pheno,"_difference")
          }
          if(mthd =="zscore"){
            inputData[[valueCol]] =   (d[d[[traitCol]]==pheno, valueCol] - control_mean_value) / sd(d[ d[[control]]==controlGroup & d[[traitCol]]==pheno, valueCol ], na.rm=T)
            inputData[[traitCol]] = paste0(pheno,"_zscore")
          }
          return(values)
        }))
        if(naTo0){
          newRows[[valueCol]][is.na(newRows[[valueCol]])]<-0
        }
        return(rbind(inputData, newRows))
      }
    }))
    return(d2)
  }))
}
return(df2)
}





