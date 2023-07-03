#' Variance partitioning using Fully Random Effects Models
#' 
#' df, des, phenotypes, timeCol, cor, returnData, combine, markSingular, time
#' 
#' @param df Dataframe containing phenotypes and design variables, optionally over time.
#' @param des Design variables to partition variance for as a character vector.
#' @param phenotypes Phenotype column names (data is assumed to be in wide format) as a character vector.
#' @param timeCol A column of the data that denotes time for longitudinal experiments. If left NULL (the default) then all data is assumed to be from one timepoint.
#' @param cor Logical, should a correlation plot be made? Defaults to TRUE.
#' @param returnData Logical, should the used to make plots be returned? Defaults to FALSE.
#' @param combine Logical, should plots be combined with patchwork? Defaults to T, which works well when there is a single timepoint being used.
#' @param markSingular Logical, should singular fits be marked in the variance explained plot? This is FALSE by default but it is good practice to check with TRUE in some situations. If TRUE this will add white markings to the plot where models had singular fits, which is the most common problem with this type of model.
#' @param time If the data contains multiple timepoints then which should be used? This can be left NULL which will use the maximum time if \code{timeCol} is specified. If a single number is provided then that time value will be used. Multiple numbers will include those timepoints. The string "all" will include all timepoints.
#' 
#' 
#' @import lme4
#' @import ggplot2
#' @import scales
#' @import patchwork
#' 
#' @keywords read.csv, pcv, wide, long
#' 
#' @return Returns either a plot (if returnData=F) or a list with a plot and data/a list of dataframes (depending on returnData and cor).
#' 
#' @examples 
#' 
#' ## Not run:
#' 
#' x<-read.pcv("/home/jsumner/Desktop/shares/mgehan_share/llima/Maize_Project_2022/nir_maize_first_exp_results.csv",
#'    filters = list("trait in area, convex_hull_area, solidity, perimeter, width, height, longest_path, convex_hull_vertices, ellipse_major_axis, ellipse_minor_axis, ellipse_angle, ellipse_eccentricity, hue_circular_mean, hue_circular_std, hue_median", "sample in default"))
#'    x$image<-sub("/shares/mgehan_share/raw_data/raw_image/maize_project_2022/MG001_E_060722_corrected/", "",x$image)
#'    write.csv(x, "/home/jsumner/Desktop/stargate/fahlgren_lab/pcvrTestData/bw_example.csv", row.names=F)
#'    
#' 
#' bw<-read.pcv("/home/jsumner/Desktop/stargate/fahlgren_lab/pcvrTestData/bw_example.csv", "long")
#' bw<-bw[nchar(bw$barcode)==13,]
#' bars<-bw$barcode
#' bars<-substr(bars,3,13)
#' genotype<-substr(bars, 1,3)
#' bars<-substr(bars,4,11)
#' meta2<-substr(bars,1,2)
#' id<-substr(bars,3,8)
#' bw$genotype = genotype
#' bw$meta = meta2
#' bw$id = id
#' bw$date<-lubridate::ymd(substr(bw$timestamp,1,10))
#' begin <- min(bw$date)
#' bw$day <- as.numeric(bw$date - begin)+1
#' 
#' des=c("genotype", "meta")
#' phenotypes = colnames(bw)[19:33]
#' timeCol = "day"
#' 
#' frem(bw, des, phenotypes, timeCol, T, F, F, T, "all")
#' frem(bw, des, phenotypes, timeCol, T, F, F, T, "all")
#' frem(bw, des, phenotypes, timeCol, T, F, F, F, NULL)
#' frem(bw, des, phenotypes, timeCol, T, F, F, T, 1)
#' 
#' ## End(Not run)
#' 
#' @export

frem<-function(df, des, phenotypes, timeCol=NULL, cor=T, returnData=F, combine=T, markSingular = F,time=NULL){
  #* `check for values`
  if(any(missing(df), missing(des), missing(phenotypes))){
    stop("df, des, phenotypes, and timeCol arguments need to be specified.")
  }
  if(is.null(timeCol)){
    timeCol = "dummy_x_axis"
    df[[timeCol]]=1
    dummyX=T
  } else{dummyX=F}
  #* `Make formulas`
  ext <- FALSE
  if(length(des)==2){
    ind_fmla <- paste0("(1|",des[1],")+(1|",des[2],")+(1|",des[1],":",des[2],")")
  }else{
    ind_fmla <- paste(paste0("(1|",des,")"),collapse = "+")
    ext <- TRUE
  }
  #* `Find time and subset data`
  if(is.null(time)){
    dat <- na.omit(df[df[[timeCol]]==max(df[[timeCol]]),c(des,phenotypes,timeCol)])
    LONGITUDINAL=F
  } else if(is.numeric(time) & length(time)==1){
    dat <- na.omit(df[df[[timeCol]]==time, c(des,phenotypes,timeCol)])
    LONGITUDINAL=F
  } else if(is.numeric(time) & length(time)>1){
    LONGITUDINAL=T
    dat <- na.omit(df[df[[timeCol]] %in% time, c(des,phenotypes,timeCol)])
  } else if(time=="all"){
    LONGITUDINAL=T
    dat <-na.omit(df[,c(des,phenotypes,timeCol)])
  }
  
  #* `Partition Variance`

  H2<-data.frame(do.call(rbind, lapply(sort(unique(dat[[timeCol]])) , function(tm){
    sub<-dat[dat[[timeCol]]==tm, ]
    do.call(rbind, lapply(phenotypes, function(e){
      
      fmla <- as.formula(paste0("as.numeric(",e,") ~ ",ind_fmla))
      model <- suppressMessages(lme4::lmer(fmla,data = sub))
      if(length(model@optinfo$conv$lme4)>=1){
        singular<-any(grepl("isSingular",model@optinfo$conv$lme4$messages))
      } else {singular=F}
      re<- lme4::VarCorr(model)
      res<-attr(lme4::VarCorr(model), "sc")^2
      
      if(!ext){
        interaction.var <- as.numeric(attr(re[[which( grepl(":", names(re))) ]],"stddev"))^2
        des1.var <- as.numeric(attr(re[[des[1]]],"stddev"))^2
        des2.var <- as.numeric(attr(re[[des[2]]],"stddev"))^2
        
        tot.var<-sum(as.numeric(re),res)
        unexp <- 1-sum(as.numeric(re))/sum(as.numeric(re),res)
        
        h2 <- c((des1.var/tot.var),
                (des2.var/tot.var),
                (interaction.var/tot.var),
                unexp,
                tm,
                singular)
      }else{
        var <- lapply(des,function(i){as.numeric(attr(re[[i]],"stddev"))^2})
        
        tot.var <- sum(as.numeric(re),res)
        unexp <- 1-sum(as.numeric(re))/sum(as.numeric(re),res)
        
        h2 <- c(unlist(var)/tot.var,unexp, tm, singular)
      }
      return(h2)
    }))
  })))

  H2$Phenotypes <- rep(phenotypes, length.out=nrow(H2))
  if(!ext){
    colnames(H2) <- c(des[1],des[2],"Interaction","Unexplained",timeCol, "singular", "Phenotypes")
  }else{
    colnames(H2) <- c(des,"Unexplained", timeCol, "singular", "Phenotypes")
  }
  ordering<-H2[H2[[timeCol]]==max(H2[[timeCol]]), ]
  H2$Phenotypes <-  ordered(H2$Phenotypes,levels= ordering$Phenotypes[order(ordering$Unexplained)] )
  H2_melt <- data.frame(data.table::melt(as.data.table(H2),id=c("Phenotypes", "singular", timeCol)))
  
  if(!ext){
    H2_melt$variable <- ordered(H2_melt$variable,levels=c("Unexplained",des[1],des[2],"Interaction"))
  }else{
    H2_melt$variable <- ordered(H2_melt$variable,levels=c("Unexplained",des))
  }
  anova_dat <- H2_melt
  
  #* `Plot Variance`
  
  if(!LONGITUDINAL){
    p <- ggplot2::ggplot(data=anova_dat)+
      ggplot2::geom_col(ggplot2::aes(y =Phenotypes, x = value, fill=variable))+
      ggplot2::xlab("Variance Explained")+
      ggplot2::guides(fill=ggplot2::guide_legend(title="", reverse=T))+
      ggplot2::theme_minimal()+
      ggplot2::scale_x_continuous(expand=c(0,0,0,0), labels = scales::percent_format())+
      ggplot2::theme(axis.text = ggplot2::element_text(size = 14),
            axis.title.y= ggplot2::element_blank())+
      ggplot2::theme(axis.ticks.length=ggplot2::unit(0.2,"cm"))+
      ggplot2::theme(legend.position = "bottom")
    
    if(markSingular){
      p<-p+ggplot2::geom_point(data=anova_dat[as.logical(anova_dat$singular) & anova_dat$variable=="Unexplained",],
                               aes(x=0.99, y=Phenotypes), color="white", shape=0)
    }
    if(dummyX){
      p<-p+ggplot2::theme(axis.title.x.bottom=ggplot2::element_blank())
    }
    
  } else{
    p <- ggplot(data = anova_dat, aes(x=.data[[timeCol]], y=value, fill=variable))+
      ggplot2::geom_area()+
      ggplot2::facet_wrap(~Phenotypes)+
      ggplot2::ylab("Variance Explained")+
      ggplot2::guides(fill=ggplot2::guide_legend(title="", reverse=T))+
      ggplot2::scale_y_continuous(expand=c(0,0,0,0), labels = scales::percent_format())+
      ggplot2::scale_x_continuous(expand=c(0,0,0,0), labels=~round(.) )+
      ggplot2::theme_minimal()+
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10),
                     axis.title.y= ggplot2::element_blank(),
                     legend.position = "bottom")
    if(markSingular){
      p<-p+ggplot2::geom_vline(data = anova_dat[anova_dat$variable=="Unexplained" & as.logical(anova_dat$singular), ],
                   ggplot2::aes(xintercept=.data[[timeCol]]), color="white", linetype=5, linewidth=0.1)
    }
    if(dummyX){
      p<-p+ggplot2::theme(axis.title.x.bottom=ggplot2::elemnt_blank())
    }
  }
  
  #* `Get Correlations`
  
  if(cor){
    #if(!LONGITUDINAL){
      corr <- cor(apply(dat[,(colnames(dat) %in% phenotypes)],2,as.numeric),use="complete.obs",method="spearman")
      unexp <- anova_dat[anova_dat$variable == "Unexplained" & anova_dat[[timeCol]]==max(anova_dat[[timeCol]]) ,]
      corr <- corr[as.character(unexp$Phenotypes[order(unexp$value)]), as.character(unexp$Phenotypes[order(unexp$value)])]

      x<-na.omit(as.data.frame(suppressWarnings(data.table::melt(as.data.table(corr), variable.name = "Var2"))))
      x$Var1<-rep(rownames(corr),length.out=nrow(x))
      
      x$Var1 <- ordered(x$Var1,levels=unexp$Phenotypes[order(unexp$value)])
      x$Var2 <- ordered(x$Var2,levels=unexp$Phenotypes[order(unexp$value)])
      
      p2 <- ggplot2::ggplot(x,ggplot2::aes(Var1,Var2))+
        ggplot2::geom_point(ggplot2::aes(color=value),size=4)+
        ggplot2::scale_color_gradient2(limits=c(-1,1),
                                       midpoint = 0)+
        ggplot2::theme_minimal()+
        ggplot2::guides(color = ggplot2::guide_colourbar(barwidth = 15,title = NULL))+
        ggplot2::labs(x="Correlations", y="")+
        ggplot2::theme(legend.position="bottom",
                       axis.text.x.bottom = ggplot2::element_text(angle=90,hjust=1))
      if(combine){
        p2<-p2+ggplot2::theme(axis.text.y.left = ggplot2::element_blank())
      }
    # } else{ # Longitudinal option
    #   x_cor<-do.call(rbind, lapply(sort(unique(na.omit(df[[timeCol]]))), function(tm){
    #     corr <- cor(apply(df[df[[timeCol]]==tm, (colnames(df) %in% phenotypes)],2,as.numeric),
    #                 use="complete.obs",method="spearman")
    #     unexp <- anova_dat[anova_dat[[timeCol]]==tm & anova_dat$variable == "Unexplained",]
    #     corr <- corr[as.character(unexp$Phenotypes[order(unexp$value)]), as.character(unexp$Phenotypes[order(unexp$value)])]
    #     x<-na.omit(as.data.frame(suppressWarnings(data.table::melt(as.data.table(corr), variable.name = "Var2"))))
    #     x$Var1<-rep(rownames(corr),length.out=nrow(x))
    #     x$Var1 <- ordered(x$Var1,levels=unexp$Phenotypes[order(unexp$value)])
    #     x$Var2 <- ordered(x$Var2,levels=unexp$Phenotypes[order(unexp$value)])
    #     x[[timeCol]]<-tm
    #     return(x)
    #   } ))
    #   ggplot(x_cor[x_cor$Var2 != x_cor$Var1, ], aes(x=.data[[timeCol]], y=value,group=interaction(Var1, Var2), color = Var1))+
    #     geom_line()+
    #     facet_wrap(~Var2)+
    #     theme_minimal()
    # }
  }

  if(cor){
    if(combine){
      plot<-p+p2
    } else{
      plot<-list(p,p2)
    }
  }else{
    plot<-p
  }
  if(returnData){
    if(cor){
      out_data<-list(H2, x)
    } else{
      out_data<-H2
    }
    out<-list(plot, out_data)
  } else{
    out<-plot
  }
  return(out)
}







