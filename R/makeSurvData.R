#' Function to reshape phenotype data to survival data based on some pcvrForm
#' 
#' @param df a dataframe to use
#' @param form a formula describing the survival analysis model (see growthSS)
#' @param model The distribution to use (model from .survModelParser )
#' 
#' @return A list including a dataframe and elements of the formula parsed.
#' 
#' @examples 
#' 
#' df <- growthSim("logistic", n=20, t=25,
#'   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#'   form <- y > 100 ~ time|id/group
#'   model = "weibull" # or binomial
#' .makeSurvData(df, form, model)
#' 
#' @keywords internal
#' @noRd

.makeSurvData <- function(df, form , model = "weibull"){
  #* `general pcvr formula parsing`
  parsed_form <- .parsePcvrForm(form, df)
  y <- parsed_form$y
  x <- parsed_form$x
  individual <- parsed_form$individual
  group <- parsed_form$group
  # USEGROUP <- parsed_form$USEG
  # USEINDIVIDUAL <- parsed_form$USEID
  df <- parsed_form$data
  #* `further survival formula steps`
  y_var <- trimws(strsplit(y, ">|[|]")[[1]][1])
  y_condition <- as.numeric(trimws(strsplit(y, ">|[|]")[[1]][2]))
  df[[y_var]] <- as.numeric(df[[y_var]] >= y_condition)
  
  df$remove_interaction <- interaction(df[[individual]], df[[group]])
  
  if(model=="weibull"){
    out_df <- do.call(rbind, lapply(unique(df$remove_interaction), function(i){
      sub <- df[df$remove_interaction == i, ]
      sub$censor <- ifelse(sub[[x]] == max(df[[x]]) & sub[[y_var]]==0, 1, 0 )
      sub[sub$censor == 1 | sub[[x]]==min(c(sub[sub[[y_var]]==1, x], Inf)), ]
    }))
    colnames(out_df)[which(colnames(out_df)==y_var)]<-"event"
    
  } else if(model =="binomial"){
    out_df <- do.call(rbind, lapply(unique(df[[x]]), function(time){
      sub <- df[df[[x]]==time, ]
      if(time != unique(df[[x]])[1] ){
        prev <- df[df[[x]]==time-1, ] 
      } else{
        prev <- sub
      }
      lt <- stats::setNames(aggregate(as.formula(paste0(y_var, " ~ ", group)), sub, sum ), c("group", "n_events"))
      lt$n_no_event <- aggregate(as.formula(paste0(y_var, " ~ ", group)), sub, function(x){sum(x==0)} )[,2]
      lt$n_eligible <- aggregate(as.formula(paste0(y_var, " ~ ", group)), prev, length)[,2]
      lt$pct_event <- lt$n_events / lt$n_eligible
      lt[[x]] <- time
      lt
    }))
    out_df[[x]] <- factor(out_df[[x]])
  }
  out_df <- out_df[,!grepl("remove_interaction", colnames(out_df))]
  ret <- list("data"= out_df, "y_var" = y_var, "y_cutoff"= y_condition,
              "x"=x, "group"=group,"individual"=individual)
  return(ret)
}





