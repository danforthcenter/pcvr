#' Summarizing function for multi-value traits
#' 
#' goal: 
#' EMD is very heavy with large datasets. For SINC2 I picked images from every 5th day and still had 6332 rows in wide format
#'    which translates to 6332^2 = 40,094,224 pairwise EMD values. In long format that's a 40 million row dataframe. 
#'    Making a network out of that is a tall order. I am currently trying it but I don't really expect it to work well.
#'    The solution here is probably at least twofold. 
#'       1: the pcv.net function should have a filter argument like net.plot that immediately filters the data before passing
#'           it into igraph
#'       2: There should be a standardized way to summarize images. I wish I could talk to Jorge about this since we both are
#'           skeptical of combining histograms in general but I can brainstorm about this and come back to Katie/Ella/Malia/Noah
#'           with ideas. Summarizing color histograms across images will be the focus of this function. Output should be a new
#'           long/wide dataframe with fewer rows (observations). 
#'           But the method for summarization depends on the hypothesis. Immediately there are three that come to mind:
#'                2A: Comparing over time within treatment group (looking for change over time of a group)
#'                2B: Comparing across treatment groups at a given time (looking for differences in design)
#'                2C: Comparing over time across treatment groups (looking for differences in design and over time)
#'           I could just ask people to subset data [arg] then give a group [arg] c(design, time)/design/time
#'             and split the data by that group
#'             summarize each split chunk into X summarized observations [arg]. Summarization could be summing, 
#'             normalizing then summing, density fitting and draw generating,
#'             normalizing and averaging, could back convert to an original scale?
#'           I think that works. The number of pieces for each group part seems clunky though. Maybe an outlength would work best.
#'           Then I could give time estimates for reasonable outlengths (I'd have to iterate and benchmark but that's probably helpful)
#'           
#'           Logistically how does the n_per_group work? If that is 1 then it's an easy sum/mean/whatever.
#'           If it is 2+ then I could randomly split the data into that many groups then do a sum/mean/whatever.
#'               If the n_per_group is higher than the actual nrow(splitData) then warn and use nrow(splitData) (calibration still matters for scale)
#' 
#' @description Cut from above
#' @param df
#' @keywords emd, earth mover's distance, multi-value trait, network
#' @examples 
#' library(pcvr)
#' hue_wide<-read.pcv("https://media.githubusercontent.com/media/joshqsumner/pcvrTestData/main/smallPhenotyperRun.csv", mode="wide", singleValueOnly = T, multiValPattern = "hist", reader="fread")
#' hue_wide$genotype = substr(hue_wide$barcode, 3,5)
#' hue_wide$genotype = ifelse(hue_wide$genotype == "002", "B73",
#'                            ifelse(hue_wide$genotype == "003", "W605S",
#'                                   ifelse(hue_wide$genotype == "004", "MM", "Mo17")))
#' hue_wide$fertilizer = substr(hue_wide$barcode, 8, 8)
#' hue_wide$fertilizer = ifelse(hue_wide$fertilizer == "A", "100",
#'                              ifelse(hue_wide$fertilizer == "B", "50", "0"))
#' hue_wide<-bw.time(hue_wide,timeCol="timestamp", group="barcode")
#' phenotypes <- colnames(hue_wide)[19:225]
#' phenoForm<-paste0("cbind(", paste0(phenotypes, collapse=", "), ")")
#' groupForm<-"DAS+barcode+genotype+fertilizer"
#' form<-as.formula(paste0(phenoForm, "~", groupForm))
#' hue_wide<-aggregate(form, data=hue_wide, mean, na.rm=T)
#' @details Cut from above
#' @return Returns a dataframe summarized by the specified groups over the multi-value traits.
#' @export


#* [args]

df=hue_wide # data
group = c("DAS", "genotype", "fertilizer") # grouping
mvCols = "hue_frequencies" # columns to take
n_per_group = 1 # first option is to say how many rows should be preserved, default to 1 but plan to generally use more?
outRows = 100 # if n_per_group is NULL OR this is not NULL then use this?
keep = NULL # column names of other phenotypes to keep, these will also be mean aggregated
#* could parallelize I guess?
#* I'm wondering if these need to distinguish between mv and sv traits. I think they do so that the sum of the MV trait can be constant.
#* [function]
mv_ag<-function(df, group, mvCols="frequencies", n_per_group=1, outRows=NULL, keep=NULL){
  #* ***** [calculated values]
  multi_group=F
  if(length(group)>1){
    original_group=group
    df$GROUP = as.character(interaction(df[,group]))
    group="INTERNAL_MULTI_GROUP"
    multi_group=T
  }
  dat_sp <-split(x=df, f=df[[group]])
  if(!is.null(outRows)){n_per_group = round(length(dat_sp)/outRows) }
  if(length(mvCols)==1 && is.character(mvCols)){ mvCols = colnames(df)[grepl(mvCols, colnames(df))] }
  if(is.numeric(mvCols)){mvCols<-colnames(df)[mvCols]}
  #* ***** [do stuff]
  out<-do.call(rbind, lapply(dat_sp, function(d){
    mv<-as.matrix(d[,mvCols])
    mv<-mv/rowSums(mv) # rescale everything to sum to 1
    if(nrow(mv) < n_per_group){ iter_n = nrow(mv) } else{ iter_n = n_per_group }
    
    nms<-sample(rownames(mv), nrow(mv), replace=F)
    index<-cut(1:nrow(mv), iter_n)
    nms_split<-split(nms, index)
  
    mv_ag<-data.frame(do.call(rbind, lapply(nms_split, function(rwnms){
      mvi<-mv[rwnms,]
      matrix(colMeans(mvi), nrow=1)
    })))
    colnames(mv_ag)<-mvCols
    
    if(!is.null(keep)){
      kept<-data.frame(do.call(rbind, lapply(nms_split, function(rwnms){
        kp<-as.matrix(d[rwnms,keep])
        matrix(colMeans(kp), nrow=1)
      })))
      colnames(kept)<-keep
      mv_ag<-cbind(kept, mv_ag)
    }
    
    mv_ag<-cbind(setNames(data.frame(rep(d[1,group], nrow(mv_ag))), group), mv_ag)
    return(mv_ag)
  }))
  if(multi_group){
    group_df = setNames(data.frame(out[[group]]), "group")
    group_df<-setNames(as.data.frame(do.call(rbind,lapply(group_df$group, function(s) matrix(strsplit(s,split="[.]")[[1]],nrow=1)))),original_group)
    out<-cbind(group_df, out)
  }
  return(out)
}





