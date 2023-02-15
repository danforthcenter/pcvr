#' Growth curve plot for bellwether data
#' 
#' @param df Data frame to use. Should be in wide format similar to output from read.pcv.bw
#' @param x X axis variable name (time variable)
#' @param y Y axis variable name (phenotype)
#' @param group Set of variables to use in grouping observations. These taken together should identify a unique plant (or unique plant at a unique angle) across time.
#' @param warn Logical, should a warning be printed if there are duplicate values for a timepoint (`x`) using the `group` variables. If TRUE then the warning will suggest how to find duplicates to amend `group`.
#' @param color Optional variable to use for line colors.
#' @param ... Additional parameters passed to geom_line
#' @keywords Bellwether, ggplot
#' @import ggplot2
#' @export
#' @examples 
#' 
#' 

 bellwetherGrowth<-function(df=NULL,x='DAS' , y='area_adj', group=c('Barcodes',"angle"),warn=T, color=NULL, ...){
   tab<-table(interaction(df[,c(x,group)]))
   if(any(tab>1) & warn){
     dataname<-deparse(substitute(df))
     nms<-names(tab)[which(as.numeric(tab)>1)]
     dupString<-paste0(dataname,"[duplicated(interaction(", paste( paste0(dataname,"$",c(x,group)), collapse=", "), ")),]")
     firstDup<-paste0(dataname,"[interaction(", paste( paste0(dataname,"$",c(x,group)), collapse=", "), ")=='",nms[1],"',]")
     eval(parse(text=dupString))
     w<-paste0("There are ", length(nms), " observations that are not uniquely identified on",
               " the `x` variable by specified `groups`.\nThe max number of duplicates is ",
               max(tab,na.rm=T),".\nRun `", dupString, "` to see the duplicated rows,\n",
               " or ",firstDup, " to see the first duplicated instance.",
               "\n\n Plot output will group these duplicates and show the mean of ", y)
     warning(w)
   }
   df$bgGrowth = interaction(bw[,group])
   
   lineLayer<-if(!is.null(color)){ ggplot2::geom_line(ggplot2::aes(x=.data[[x]], y=.data[[y]], group=.data$bgGrowth, color = .data[[color]]),linewidth=0.25, ... )
   } else{ ggplot2::geom_line(ggplot2::aes(x=.data[[x]], y=.data[[y]], group=.data$bgGrowth),linewidth=0.25, ... ) }
   
   ggplot2::ggplot(df)+
     lineLayer+
     pcv_theme()
 }
 