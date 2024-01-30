#' Ease of use nlrq starter function for 6 growth model parameterizations
#' 
#' Internal to growthSS
#' 
#' @examples 
#' 
#' simdf<-growthSim("logistic", n=20, t=25,
#'   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' 
#' ss<-.mgcvSS(model = "gam", form=y~time|id/group, df=simdf)
#' names(ss) # formula, df, pcvrForm
#' 
#' @keywords internal
#' @noRd

.mgcvSS<-function(model="gam", form, df){
  #* `parse form argument`
  y=as.character(form)[2]
  x<-as.character(form)[3]
  group = NULL
  if(grepl("\\|", x) & grepl("\\/",x)){
    x3<-trimws(strsplit(x, "[|]|[/]")[[1]])
    x<-x3[1]
    individual = x3[2]
    group = x3[length(x3)] 
    message(paste0("Individual is not used with type = 'gam'."))
    if(length(unique(df[[group]]))==1){USEGROUP=FALSE}else{USEGROUP=TRUE} # if there is only one group then ignore grouping for parameter and variance formulas
  } else if (grepl("\\|", x)){
    x2<-trimws(strsplit(x, "[|]")[[1]])
    x<-x2[1]
    group = x2[2] # individual will not be used here 
    USEGROUP=TRUE # leave x as is, in this parameterization this means no term[group] syntax
    message(paste0("Individual is not used with type = 'gam', interpreting ", group, ", as a group."))
  } else {
    # leave x as is, in this parameterization this means no term[group] syntax
    USEGROUP=FALSE
  }
  if(USEGROUP){
    df[[group]]<-factor(df[[group]])
    df[[paste0(group,"_numericLabel")]]<-unclass(df[[group]])
  }
  #* `assemble gam formula`
  if(USEGROUP){
    gam_form <- stats::as.formula(paste0(y, "~0+",group,"+s(", x, ", by=", group, ")"))
  } else{
    gam_form <- stats::as.formula(paste0(y, "~0+s(", x, ")"))
  }
  #* `return list`
  out<-list()
  out[["formula"]] <- gam_form
  out[["df"]]<-df
  out[["pcvrForm"]]<-form
  return(out)
}
