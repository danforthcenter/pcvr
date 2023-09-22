#' Ease of use nlrq starter function for 6 growth model parameterizations
#' 
#' Internal to growthSS
#' 
#' @examples 
#' 
#' simdf<-growthSim("logistic", n=20, t=25,
#' params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
#' 
#' ss<-.nlrqSS(model = "logistic", form=y~time|id/group,
#'   tau=0.5, df=simdf, start=NULL)
#' dim(ss$df)
#' ss[c("formula", "taus", "start", "pcvrForm")]
#' 
#' @keywords internal
#' @noRd
if(FALSE){
  
library(nlme)
ex<-growthSim("logistic", n=20, t=25,
                 params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))
x="time"
y="y"
as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C))"))

ex$group<-factor(ex$group) # group MUST be a factor for nlme (typical)
# sigma ~ none
nl1<-nlme(model = y ~ A/(1 + exp((B - time)/C)),
     data = ex,
     fixed = list(
       A ~ 0 + group,
       B ~ 0 + group,
       C ~ 0 + group),
     random = A ~ 1,
     groups = ~group,
     weights = varIdent(form = ~1|group),
     start = c(150, 150,10, 10,3, 3) )
# linear sigma using power
nl2<-nlme(model = y ~ A/(1 + exp((B - time)/C)),
     data = ex,
     fixed = list(
       A ~ 0 + group,
       B ~ 0 + group,
       C ~ 0 + group),
     random = A ~ 1,
     groups = ~group,
     weights = varPower(form = ~time|group),
     start = c(150, 150,10, 10,3, 3) )
# linear sigma using exp
nl3<-nlme(model = y ~ A/(1 + exp((B - time)/C)),
     data = ex,
     fixed = list(
       A ~ 0 + group,
       B ~ 0 + group,
       C ~ 0 + group),
     random = A ~ 1,
     groups = ~group,
     weights = varExp(form = ~time|group),
     start = c(150, 150,10, 10,3, 3) )

anova(nl1, nl2, nl3)

}




