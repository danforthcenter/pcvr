#' Hypothesis testing for frequentist \code{fitGrowth} models.
#' 
#' @param ss A list output from \link{growthSS}. This is not required if \code{test} is given
#' in \code{brms::hypothesis} style as a written statement.
#' @param fit A non-brms model (or list of nlrq models) output from \link{fitGrowth}.
#' @param test A description of the hypothesis to test. This can take two main forms, either the parameter names
#' to vary before comparing a nested model ("A", "B", "C") using an anova or a hypothesis test/list of hypothesis tests 
#' written as character strings. The latter method is not implemented for \code{nlrq} models.
#' If this is a vector of parameters to test in the model then they should be parameters which
#' vary by group in your original model and that you want to test against a null model where they
#' do not vary by group. Alternatively for nlrq models this can be a comparison of model terms
#' written as \code{"group_X|tau|par > group_Y|tau|par"}, which uses a fat tailed T distribution to make 
#' comparisons on the means of each quantile estimate. For GAMs these tests compare the model with
#' splines either by group or interacting with group to a model that ignores the grouping in the data.
#' If this is a list of hypothesis tests then they should describe tests similar to "A.group1 - A.group2*1.1" 
#' and can be thought of as contrasts.
#' @keywords hypothesis, nlme, nls, nlrq
#' @importFrom stats getCall logLik pchisq anova as.formula setNames vcov coef pt
#' @importFrom nlme pdIdent corAR1 fixef
#' @importFrom extraDistr qlst dlst
#' @importFrom methods is
#' @importFrom car deltaMethod
#' 
#' @details
#' For nls and nlme models an anova is run and returned as part of a list along with the null model.
#'  For nlrq models several assumptions are made and a likelihood ratio test for each tau 
#'  is run and returned as a list.
#' 
#' 
#' @return A list containing an anova object comparing non-linear growth models and the null model.
#' 
#' @examples 
#' 
#' ## Not run:
#' set.seed(123)
#' simdf<-growthSim("logistic", n=20, t=25,
#'    params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-suppressMessages(growthSS(model = "logistic", form=y~time|id/group,
#'   df=simdf, type="nlrq"))
#' fit <- fitGrowth(ss)
#' testGrowth(ss, fit, "A")
#' 
#' ## End(Not run)
#' @export

testGrowth<-function(ss=NULL, fit, test = "A"){
  
  if(all(unlist(lapply(test, nchar)) <= 2)){
    method = "anova"
  } else{method = "contrast"}
  
  
  if(method=="contrast"){
    
    if(!is.null(ss)){
      if(ss$type %in% c("nlrq", "mgcv")){
        stop("These hypotheses are only implemented for nls and nlme models.")
      } else if(ss$type=="brms"){
        stop("For brms model tests use brms::hypothesis")
      } 
    }
    
    res <- .nlhypothesis(fit, test)
    
  } else if(method=="anova"){
    if(ss$model=="gam"){
      #* do gam things
      if(ss$type == "nls"){
        res <- .lmGamAnova(ss, fit)
      } else if(ss$type == "nlrq"){
        res <- .rqGamAnova(ss, fit)
      } else if(ss$type == "nlme"){
        res <- .lmeGamAnova(ss, fit)
      } else if(ss$type == "mgcv"){
        res <- .mgcvGamAnova(ss, fit)
      } else if(ss$type=="brms"){
        stop("For brms model tests use brms::hypothesis")
      } 
    } else if(ss$type %in% c("nls", "nlme")){
      res <- .nlsAnova(ss, fit, test_pars = test)
    } else if(ss$type == "nlrq"){
      if(any(test %in% c("A", "B", "C")) ){
        res <- .nlrqTest(ss, fit, test_pars = test) 
      } else{
        res<-.nlrqTest2(ss, fit, test_pars = test)
      }
    } else if(ss$type=="brms"){
      stop("For brms model tests use brms::hypothesis")
    } 
  }
  
  return(res)
}

#' (non)linear hypothesis function for nls and nlrq models
#' @examples
#' if(FALSE){
#' set.seed(123)
#' logistic_df<-growthSim("logistic", n=20, t=25,
#'                   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "logistic", form=y~time|id/group, 
#'           df=logistic_df, type = "nls")
#' fit <- fitGrowth(ss)
#' .nlhypothesis(fit, test = list("A1 - A2", "B1-B2"))
#' .nlhypothesis(fit, test = "A1 + 10 - A2 * 1.1")
#' .nlhypothesis(fit, test = "A1/B1 - A2/B2")
#' }
#' @keywords internal
#' @noRd

.nlhypothesis <- function(fit, test){
  if(methods::is(fit, "nlme")){
    coefs <- nlme::fixef(fit)
  } else if(methods::is(fit, "nls")){
    coefs <- stats::coef(fit)
  } else{
    stop("These hypotheses are only implemented for nls and nlme models.")
  }
  dfresid <- summary(fit)$df[2]   
  vcMat <- stats::vcov(fit)
  hypotheses <- data.frame(form = unlist(test))
  colnames(hypotheses) <- c("Form") 
  val <- do.call(rbind, lapply(1:nrow(hypotheses), function(i){
    car::deltaMethod(object=coefs, g=as.character(hypotheses$Form[i]), 
                     vcov.=vcMat, level = 0.95 )
  }))
  row.names(val) <- 1:nrow(hypotheses)
  val <- cbind(hypotheses, val)
  lenVal <- ncol(val)
  out <- val[, 1:(lenVal-2)]
  out$"t-value" <- abs(out$Estimate/out$SE)
  residual_DF <- ifelse(is.null(dfresid), Inf, dfresid)
  out$"p-value" <- 2 * pt(out$"t-value", residual_DF, lower.tail = F)
  if(residual_DF == Inf) colnames(out)[length(colnames(out)) - 1] <- "Z-value"
  return(out)
}



#' mgcv gam testing function
#' @examples
#' if(FALSE){
#' set.seed(123)
#' logistic_df<-growthSim("logistic", n=20, t=25,
#'                   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "gam", form=y~time|id/group, 
#'           df=logistic_df, type = "mgcv")
#' fit <- fitGrowth(ss)
#' .mgcvGamAnova(ss, fit)$anova
#' }
#' @keywords internal
#' @noRd

.mgcvGamAnova <- function(ss, fit){
  #* `Get x variable`
  RHS <- as.character(ss$formula)[3]
  x <- sub(",", "", sub("s\\(", "", regmatches(RHS, regexpr("s\\(.*,", RHS))))
  ssNew <- ss
  #* ***** `Make Null formula`
  ssNew$formula <- stats::as.formula(paste0("y ~ s(",x,")"))
  #* `rerun fitGrowth with new formula`
  nullMod <- fitGrowth(ssNew)
  #* `compare models and return values`
  anv <- stats::anova(nullMod,fit, test="F")
  out <- list("anova" = anv, "nullMod"=nullMod)
  return(out)
}



#' lme gam testing function
#' @examples
#' if(FALSE){
#' set.seed(123)
#' logistic_df<-growthSim("logistic", n=20, t=25,
#'                   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "gam", form=y~time|id/group, 
#'           df=logistic_df, type = "nlme", sigma="power")
#' fit <- fitGrowth(ss)
#' .lmeGamAnova(ss, fit)$anova
#' }
#' @keywords internal
#' @noRd

.lmeGamAnova <- function(ss, fit){
  #* `Get x variable`
  x<- trimws(sub("\\*.*","",as.character(ss$formula$model)[3]))
  ssNew <- ss
  #* ***** `Make Null formulas`
  newdf <- ss$df
  newdf$dummy <- "A"
  model_form <- ss$formula$model # fixed effects should be constant for likelihood testing
  # stats::as.formula(paste0(sub("\\*", "+", ss$formula$model)[c(2,1,3)], collapse=""))
  #* random effects formula
  random_form <- stats::setNames(list(nlme::pdIdent(~splines - 1, data = newdf) ), "dummy")
  #* variance formula
  weights_form <- .nlme_sigma_form(as.list(ss$call)$sigma, x, "dummy")
  #* correlation formula
  correlation_form <- nlme::corAR1(0.8, form = stats::as.formula(paste0("~ 1 | dummy")))
  #* add all to newSS as list
  ssNew$formula <- list("model"=model_form, "random"=random_form,
                   "weights" = weights_form, "cor_form" = correlation_form)
  ssNew$df <- newdf
  #* `rerun fitGrowth with new formula`
  nullMod <- fitGrowth(ssNew)
  #* `compare models and return values`
  anv <- stats::anova(nullMod,fit)
  out <- list("anova" = anv, "nullMod"=nullMod)
  return(out)
}


#' rq gam testing function for multiple taus
#' @examples
#' if(FALSE){
#' set.seed(123)
#' logistic_df<-growthSim("logistic", n=20, t=25,
#'                   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "gam", form=y~time|id/group, 
#'           df=logistic_df, type = "nlrq", tau=c(0.25, 0.5, 0.75))
#' fit <- fitGrowth(ss, cores=3)
#' .rqGamAnova(ss, fit)$'0.25'$anova
#' .rqGamAnova(ss, fit[[1]])$anova
#' }
#' @keywords internal
#' @noRd

.rqGamAnova <- function(ss, fit){
  if(methods::is(fit, "rq")){
    tau <- as.character(fit$tau)
    fit <- stats::setNames(list(fit), tau)
  }
  taus <- names(fit)
  res <- lapply(taus, function(tau){
    f <- fit[[tau]]
    #* `remove grouping from ss$formula`
    rhs <- as.character(ss$formula)[3]
    newRhs <- trimws(sub("\\*.*", "", rhs))
    newFormula <- stats::as.formula( paste0(as.character(ss$formula)[2], "~", newRhs) )
    #* `Make new SS object`
    ssNew <- ss
    ssNew$formula <- newFormula
    ssNew$taus <- as.numeric(tau)
    #* `rerun fitGrowth with new formula`
    nullMod <- fitGrowth(ssNew)
    #* `compare models and return values`
    anv <- stats::anova(nullMod,f)
    out <- list("anova" = anv, "nullMod"=nullMod)
    return(out)
  })
  names(res) <- taus
  if(length(res)==1){ res <- res[[1]] }
  return(res)
}



#' lm gam testing function
#' @examples
#' if(FALSE){
#' set.seed(123)
#' logistic_df<-growthSim("logistic", n=20, t=25,
#'                   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "gam", form=y~time|id/group, 
#'           df=logistic_df, type = "nls")
#' fit <- fitGrowth(ss)
#' .lmGamAnova(ss, fit)$anova
#' }
#' @keywords internal
#' @noRd

.lmGamAnova <- function(ss, fit){
  #* `remove grouping from ss$formula`
  rhs <- as.character(ss$formula)[3]
  newRhs <- trimws(sub("\\*.*", "", rhs))
  newFormula <- stats::as.formula( paste0(as.character(ss$formula)[2], "~", newRhs) )
  #* `Make new SS object`
  ssNew <- ss
  ssNew$formula <- newFormula
  #* `rerun fitGrowth with new formula`
  nullMod <- fitGrowth(ssNew)
  #* `compare models and return values`
  anv <- stats::anova(nullMod,fit)
  out <- list("anova" = anv, "nullMod"=nullMod)
  return(out)
}


#' nls and nlme testing function
#' @examples
#' if(FALSE){
#' set.seed(123)
#' logistic_df<-growthSim("logistic", n=20, t=25,
#'                   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "logistic", form=y~time|id/group, 
#'           df=logistic_df, type = "nls")
#' fit <- fitGrowth(ss)
#' .nlsAnova(ss, fit, test_pars = "A")$anova
#' }
#' @keywords internal
#' @noRd

.nlsAnova <- function(ss, fit, test_pars = "A"){
  #* `Get parameters to vary in comparison model`
  xForm <- as.character(ss$formula)[3]
  rMatches <- gregexpr(".2?\\[", xForm)
  original_grouped_pars <- sub("\\[", "", regmatches(xForm, rMatches)[[1]])
  null_pars = original_grouped_pars[!original_grouped_pars %in% test_pars]
  #* `match call for previous model, updating pars`
  lcall <- as.list(ss$call)
  lcall$pars = null_pars
  lcall$df = ss$df
  new_call <- as.call(lcall)
  #* `rerun growthSS and fitGrowth with updated parameters`
  nullSS <- suppressMessages(eval(new_call))
  nullMod <- fitGrowth(nullSS)
  #* `compare models and return values`
  anv <- stats::anova(nullMod,fit)
  out <- list("anova" = anv, "nullMod"=nullMod)
  return(out)
}

#' pseudo LRT function for nlrq
#' @examples
#' if(FALSE){
#' set.seed(123)
#' logistic_df<-growthSim("logistic", n=20, t=25,
#'                   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "logistic", form=y~time|id/group, 
#'           df=logistic_df, type = "nlrq", tau=c(0.25, 0.5, 0.75))
#' fit <- fitGrowth(ss)
#' 
#' .nlrqTest(ss, fit, test_pars = "A")$anova
#' }
#' @keywords internal
#' @noRd

.nlrqTest <- function(ss, fit, test_pars = "A"){
#* `Get parameters to vary in comparison model`
xForm <- as.character(ss$formula)[3]
rMatches <- gregexpr(".2?\\[", xForm)
original_grouped_pars <- sub("\\[", "", regmatches(xForm, rMatches)[[1]])
null_pars = original_grouped_pars[!original_grouped_pars %in% test_pars]
#* `match call for previous model, updating pars`
lcall <- as.list(ss$call)
lcall$pars = null_pars
lcall$df = ss$df
new_call <- as.call(lcall)
#* `rerun growthSS and fitGrowth with updated parameters`
nullSS <- suppressMessages(eval(new_call))
nullMods <- fitGrowth(nullSS)
if(methods::is(fit, "nlrq")){
  fit <- list(fit)
  names(fit)<-ss$taus
  nullMods <- list(nullMods)
  names(nullMods)<-ss$taus
}

#* `arrange models for comparisons`

modsList <- lapply(ss$taus, function(tau) list(fit[[paste0(tau)]], nullMods[[paste0(tau)]]))

res <- lapply(modsList, function(modsPair){
  .nlrq_pseudoLRT(modsPair)
})
names(res)<-ss$taus
return(res)
}


#' pseudo LRT function for nlrq, aiming to use empirical likelihood given more time
#' @param nested_models list of models with the same tau generated by .nlrqTest
#' @keywords internal
#' @noRd

.nlrq_pseudoLRT <- function(nested_models){
  nModels = length(nested_models) # currently always 2, but might change in the future.
  rval <- matrix(rep(NA, 5 * 2), ncol = 5)
  colnames(rval) <- c("#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)")
  rownames(rval) <- 1:2
  logL <- lapply(nested_models, stats::logLik)
  rval[,1] <- as.numeric(sapply(logL, function(x) attr(x, "df")))  
  rval[,2] <- sapply(logL, as.numeric)
  rval[2:nModels, 3] <- rval[2:nModels, 1] - rval[1:(nModels-1), 1]
  rval[2:nModels, 4] <- 2 * abs(rval[2:nModels, 2] - rval[1:(nModels-1), 2])
  rval[,5] <- stats::pchisq(rval[,4], round(abs(rval[,3])), lower.tail = FALSE)
  
  variables <- lapply(nested_models, function(x){deparse(as.list(stats::getCall(x))$formula)})
  header <- paste("Model ", format(1:nModels),": ", variables, sep="", collapse="\n")
  out<-structure(as.data.frame(rval), heading = header,
                 class = c("anova", "data.frame"))
  return(out)
}


#' Parameterized test for nlrq models using SE
#' @examples
#' if(FALSE){
#' set.seed(123)
#' logistic_df<-growthSim("logistic", n=20, t=25,
#'                   params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' ss<-growthSS(model = "logistic", form=y~time|id/group, 
#'           df=logistic_df, type = "nlrq", tau=c(0.25, 0.5, 0.75))
#' fit <- fitGrowth(ss)
#' 
#' .nlrqTest2(ss, fit, test_pars = "a|0.5|A > b|0.5|A")
#' }
#' @keywords internal
#' @noRd

.nlrqTest2 <- function(ss, fit, test_pars = "a|0.5|A > b|0.5|A"){
  #* `Get parameters to vary in comparison model`
  xForm <- as.character(ss$formula)[3]
  rMatches <- gregexpr(".2?\\[", xForm)
  original_grouped_pars <- sub("\\[", "", regmatches(xForm, rMatches)[[1]])
  #* `Get name of grouping variable in formula`
  rMatches2 <- regexpr("\\[[a-zA-Z0-9.]*\\]", xForm)
  groupVar <- sub("\\]", "", sub("\\[", "", regmatches(xForm, rMatches2)))
  #* `parse test_pars formula`
  split_form <- lapply(strsplit(test_pars, ">|<")[[1]], function(s){
      stats::setNames(strsplit(trimws(s), "\\|")[[1]], c("group", "tau", "par"))
      })
  direction <- ifelse(grepl(">", test_pars), "greater", "lesser")
  #* `select models based on tau`
  if(!methods::is(fit, "nlrq")){
    fit1 <- fit[[split_form[[1]]["tau"] ]]
    fit2 <- fit[[split_form[[2]]["tau"] ]]
    fits <- list(fit1, fit2)
  } else{ fits <- list(fit, fit)}
  #* `get model data for each group in formula`
  mdf <- do.call(rbind, lapply(1:length(split_form), function(i){
    sf <- split_form[[i]]
    #* `get model parameters`
    mdf <- as.data.frame(summary(fits[[i]])$coefficients)
    #* `replace 1-nGroups with sorted group names`
    mdf$par <- substr(rownames(mdf),1,1)
    mdf$numericGroup<-sub("^[A-C]", "", rownames(mdf))
    
    gNames <- stats::setNames(data.frame(sort(unique(ss$df[[groupVar]])), 1:length(unique(ss$df[[groupVar]]))),
                       c(groupVar, "numericGroup"))
    mdf <- merge(mdf, gNames, by="numericGroup")
    return( mdf[mdf[[groupVar]]==sf[["group"]] & mdf[["par"]]==sf[["par"]],
                c("Value", "Std. Error", "t value", "Pr(>|t|)", "par", "group")] )
  }))
  #* `make T distributions for comparisons`
  support1 <- extraDistr::qlst(c(0.0001, 0.9999), df = 5, mu = mdf[1,"Value"], sigma = mdf[1, "Std. Error"])
  support2 <- extraDistr::qlst(c(0.0001, 0.9999), df = 5, mu = mdf[2,"Value"], sigma = mdf[2, "Std. Error"])
  mn <- 0.99 * min(c(support1, support2))
  mx <- 1.01 * max(c(support1, support2))
  support <- seq(mn, mx, length.out=10000)
  t1<-extraDistr::dlst(support, df = 5, mu = mdf[1,"Value"], sigma = mdf[1, "Std. Error"])
  t2<-extraDistr::dlst(support, df = 5, mu = mdf[2,"Value"], sigma = mdf[2, "Std. Error"])
  pdf1 <- t1/sum(t1)
  pdf2 <- t2/sum(t2)
  res <- .post.prob.from.pdfs(pdf1, pdf2, direction)$post.prob
  mdf[[paste0("prob.", direction)]]<-c(res, NA)
  return(mdf)
}










