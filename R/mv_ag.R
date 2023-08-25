#' Multi Value Trait Aggregation function
#' 
#' @description EMD can get very heavy with large datasets. For an example
#' lemnatech dataset filtering for images from every 5th day there are
#' 6332^2 = 40,094,224 pairwise EMD values. In long format that's a 40 million row dataframe,
#' which is unwieldy. This function is to help reduce the size of datasets before
#' comparing histograms and moving on with matrix methods or network analysis.
#' 
#' @param df A dataframe with multi value traits. This can be in wide or long format,
#' data is assumed to be long if traitCol, valueCol, and labelCol are present.
#' @param group Vector of column names for variables which uniquely identify groups
#' in the data to summarize data over. Typically this would be the design variables
#' and a time variable.
#' @param mvCols Either a vector of column names/positions representing multi value
#' traits or a character string that identifies the multi value trait columns as a
#' regex pattern. Defaults to "frequencies".
#' @param n_per_group Number of rows to return for each group.
#' @param outRows Optionally this is a different way to specify how many rows to return.
#' This will often not be exact so that groups have the same number of observations each.
#' @param keep A vector of single value traits to also average over groups, if there are
#' a mix of single and multi value traits in your data.
#' @param parallel Optionally the groups can be run in parallel with this number of cores,
#' defaults to 1 if the "mc.cores" option is not set globally.
#' @param traitCol Column with phenotype names, defaults to "trait".
#' @param labelCol Column with phenotype labels (units), defaults to "label".
#' @param valueCol Column with phenotype values, defaults to "value".
#' @param id Column that uniquely identifies images if the data is in long format. 
#' This is ignored when data is in wide format.
#' @keywords emd, earth mover's distance, multi-value trait, network
#' @import parallel
#' @importFrom stats setNames aggregate as.formula
#' @examples 
#' 
#' ## Not run:
#' 
#' hue_wide<-read.pcv(paste0(
#' "https://media.githubusercontent.com/media/joshqsumner/",
#' "pcvrTestData/main/pcv4-multi-value-traits.csv"),
#'  mode="wide", reader="fread")
#'  
#' hue_wide$genotype = substr(hue_wide$barcode, 3,5)
#' hue_wide$genotype = ifelse(hue_wide$genotype == "002", "B73",
#'                            ifelse(hue_wide$genotype == "003", "W605S",
#'                                   ifelse(hue_wide$genotype == "004", "MM", "Mo17")))
#' hue_wide$fertilizer = substr(hue_wide$barcode, 8, 8)
#' hue_wide$fertilizer = ifelse(hue_wide$fertilizer == "A", "100",
#'                              ifelse(hue_wide$fertilizer == "B", "50", "0"))
#' hue_wide<-bw.time(hue_wide,timeCol="timestamp", group="barcode")
#' phenotypes <- colnames(hue_wide)[grepl("hue_frequencies", colnames(hue_wide))]
#' phenoForm<-paste0("cbind(", paste0(phenotypes, collapse=", "), ")")
#' groupForm<-"DAS+barcode+genotype+fertilizer"
#' form<-as.formula(paste0(phenoForm, "~", groupForm))
#' hue_wide<-aggregate(form, data=hue_wide, mean, na.rm=TRUE)
#' dim(hue_wide)
#' hue_ag1<-mv_ag(df=hue_wide, group = c("DAS", "genotype", "fertilizer"),
#'  n_per_group=2)
#' dim(hue_ag1)
#' hue_ag2<-mv_ag(hue_wide, group = c("DAS", "genotype", "fertilizer"),
#'  n_per_group=1)
#' dim(hue_ag2)
#' 
#' hue_long<-read.pcv(paste0(
#' "https://media.githubusercontent.com/media/joshqsumner/",
#' "pcvrTestData/main/pcv4-multi-value-traits.csv"),
#'  reader="fread")
#' hue_long$genotype = substr(hue_long$barcode, 3,5)
#' hue_long$genotype = ifelse(hue_long$genotype == "002", "B73",
#'                            ifelse(hue_long$genotype == "003", "W605S",
#'                                   ifelse(hue_long$genotype == "004", "MM", "Mo17")))
#' hue_long$fertilizer = substr(hue_long$barcode, 8, 8)
#' hue_long$fertilizer = ifelse(hue_long$fertilizer == "A", "100",
#'                              ifelse(hue_long$fertilizer == "B", "50", "0"))
#' hue_long<-bw.time(hue_long,timeCol="timestamp", group="barcode")
#' dim(hue_long)
#' hue_ag3 <- mv_ag(hue_long, group = c("DAS", "genotype", "fertilizer"),
#'   mvCols = "hue", n_per_group = 2)
#' dim(hue_ag3)
#' head(hue_ag3)
#' 
#' ## End(Not run)
#' 
#' @return Returns a dataframe summarized by the specified groups over the multi-value traits.
#' @export

mv_ag<-function(df, group, mvCols="frequencies", n_per_group=1, outRows=NULL, keep=NULL, parallel=getOption("mc.cores",1),
                traitCol="trait", labelCol="label", valueCol="value", id = "image"){
  
  #* ***** [decide if data is long or wide]
  if(all(c(traitCol, valueCol, labelCol) %in% colnames(df))){
    long = TRUE
  } else if(!any(c(traitCol, valueCol, labelCol) %in% colnames(df))){
    long = FALSE } else{
      found <- c('traitCol', 'valueCol', 'labelCol')[which(c(traitCol, valueCol, labelCol) %in% colnames(df1))]
      stop(paste0( paste(found, collapse = ", "), " found in column names of data but either all or none of traitCol, valueCol, and labelCol are expected."))
    }
  
  
  #* ***** [calculated values]
  multi_group=FALSE
  if(length(group)>1){
    original_group=group
    df$INTERNAL_MULTI_GROUP = as.character(interaction(df[,group]))
    group="INTERNAL_MULTI_GROUP"
    multi_group=TRUE
  }
  
  
  #* ***** [wide column selection]
  if(!long){
    if(length(mvCols)==1 && is.character(mvCols)){ mvCols = colnames(df)[grepl(mvCols, colnames(df))] }
    if(is.numeric(mvCols)){mvCols<-colnames(df)[mvCols]}
  } else{
    #* ***** [long trait identification]
    
    df <- df[grepl(mvCols, df[[traitCol]]), ]
    if(length(unique(df[[traitCol]])) > 1){
      stop(paste0("In long format mvCols should only match one trait, ", 
                  mvCols, " matches ", paste0(unique(df[[traitCol]]), collapse=", ") ))
    }
  }
  
  #* ***** [split data by group]
  
  dat_sp <-split(x=df, f=df[[group]])
  if(!is.null(outRows)){n_per_group = round(length(dat_sp)/outRows) }
  
  #* ***** [aggregate wide format data]
  if(!long){
    out<-do.call(rbind, parallel::mclapply(dat_sp, function(d){
      mv<-as.matrix(d[,mvCols], rownames.force=TRUE)
      mv<-mv/rowSums(mv) # rescale everything to sum to 1
      if(nrow(mv) < n_per_group){ iter_n = nrow(mv) } else{ iter_n = n_per_group }
      if(is.null(rownames(mv))){rownames(mv)<-1:nrow(mv)} # should be redundant
      nms<-sample(rownames(mv), nrow(mv), replace=FALSE)
      if(nrow(mv)>1 & iter_n>1){ 
        index<-cut(1:nrow(mv), iter_n)
        nms_split<-split(nms, index)
      } else if(nrow(mv)>1 & iter_n==1 ){
        nms_split = list(rownames(mv))
      } else{
          nms_split = list(rownames(mv))
        }
    
      mv_ag<-data.frame(do.call(rbind, lapply(nms_split, function(rwnms){
        mvi<-matrix(mv[rwnms,], nrow=length(rwnms))
        if(nrow(mvi)>1){
            matrix(colMeans(mvi), nrow=1)
        } else{
            mvi
          }
      })))
      colnames(mv_ag)<-mvCols
      
      if(!is.null(keep)){
        kept<-data.frame(do.call(rbind, lapply(nms_split, function(rwnms){
          kp<-as.matrix(d[rwnms,keep])
          if(nrow(kp)>1){
            matrix(colMeans(kp), nrow=1)
          } else{
            kp
          }
        })))
        colnames(kept)<-keep
        mv_ag<-cbind(kept, mv_ag)
      }
      
      mv_ag<-cbind(setNames(data.frame(rep(d[1,group], nrow(mv_ag))), group), mv_ag)
      return(mv_ag)
    }, mc.cores=parallel) )
  } else{
    #* ***** [aggregate long data]
    
    out<-do.call(rbind, parallel::mclapply(dat_sp, function(d){
      #* get unique images
      IDS <- unique(d[[id]])
      #* define number of possible groups
      if(length(IDS) < n_per_group){ iter_n = length(IDS) } else{ iter_n = n_per_group }
      #* rescale values to sum to 1
      d <- do.call(rbind, lapply(IDS, function(i){
        id_d <- d[d[[id]]==IDS,]
        id_d[[valueCol]] <- id_d[[valueCol]]/sum(id_d[[valueCol]], na.rm=T)
        return(id_d)
      }))
      #* separate IDS into groups
      if(length(IDS)>1 & iter_n>1){
        index<-cut(1:length(IDS), iter_n)
        ids_split<-split(IDS, index)
      } else if(length(IDS)>1 & iter_n==1 ){
        ids_split = list(IDS)
      } else{
        ids_split = list(IDS)
      }
      
      #* mean of bin per groups
      mv_ag<-data.frame(do.call(rbind, lapply(ids_split, function(ids){
        d_group <- d[d[[id]] %in% ids,]
        aggregate(as.formula(paste0(valueCol, "~",paste0(c(labelCol, traitCol, group),collapse="+"))), d_group, mean, na.rm=T)
      })))
      #* add back in kept traits
      if(!is.null(keep)){
        kept<-d[d[[traitCol]]%in% keep, ]
        mv_ag<-cbind(mv_ag, kept)
      }
      rownames(mv_ag)<-NULL
      return(mv_ag)
    }, mc.cores=parallel) )
  }

  if(multi_group){
    group_df = setNames(data.frame(out[[group]]), "group")
    group_df<-setNames(as.data.frame(do.call(rbind,lapply(group_df$group, function(s) matrix(strsplit(s,split="[.]")[[1]],nrow=1)))),original_group)
    out<-cbind(group_df, out[,-which(colnames(out)=="INTERNAL_MULTI_GROUP")])
  }
  
  return(out)
}





