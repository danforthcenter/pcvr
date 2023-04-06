#' Ease of use brms starter functions
#' 
#' maybe I should make a data simulating function for each supported distribution as well? 
#' something that takes the parameters, the number of samples per group, and returns a dataframe?
#' then the self starter functions should make the bf() call for each supported distribution
#' the self starter function could also return a list of inits for the model. I think that makes sense
#' as a component of the self-starter function because it is a self starting component.
#' This could optionally make a list of priors for you taking similar arguments as the data simulating function
#' and just centering a lognormal distribution on those means. The default would be not to do that though. I guess.
#' It could also take an argument for how many cores/chains to use? Not sure about that yet.
#' 
#' Then the actual brm call would be something like
#' brm(ss$bf, family = student, prior = prior, data = df, iter = 2000, 
#' cores = 4, chains = 4, backend = "cmdstanr",  
#' control = list(adapt_delta = 0.999, max_treedepth = 20),
#' init =ss$init)
#' 
#' @param df Data frame to use. Can be wide or long format from read.pcv
#' @keywords Bayesian, brms
#' @import brms
#' @import ggplot2
#' @examples 
#' 
#' @export
#' 
#' library(brms)
#' 
#' 

simdf<-growthSim("logistic", n=20, t=25, params = list("A"=c(200,160), "B"=c(13, 11), "C"=c(3, 3.5)))

prior_spline <- prior(lognormal(log(130), .25),nlpar = "phi1") +
  prior(lognormal(log(12), .25), nlpar = "phi2") +
  prior(lognormal(log(3), .25), nlpar = "phi3") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit_spline <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                     sigma~s(time,by=group), 
                     phi1 + phi2 + phi3 ~ 0+group,
                     autocor = ~arma(~time|id:group,1,1), 
                     nl = TRUE),
                  family = student, prior = prior_spline, data = simdf, iter = 2000, 
                  cores = 4, chains = 4, backend = "cmdstanr",  
                  control = list(adapt_delta = 0.999, max_treedepth = 20),
                  init = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))}) 


#* `args`
model = c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law")
form = y~time
#* maybe if form has a pipe then it should parse that for group and individual rather than using the argument options. It would 
#* make sense to me to be able to just give a formula and a name for a kind of heteroskedastic sub model

form = y~time|id/group

sigma = "spline" #c("homo", "linear", "spline")]
group = "group"
autocor=T
individual = c("id", "group")
data=simdf

#* bf components:
#* growth formula
#* sigma forula
#* parameter group formula: phi1 + phi2 + phi3 ~ 0+treatment,
#* autocorrelation formula: NAMED autocor = ~arma(~time|sample:treatment,1,1), 
#* nl statement

sigmas<-c("homo", "linear", "spline")
models<-c("logistic", "gompertz", "monomolecular", "exponential", "linear", "power law")
y=as.character(form)[2]
x<-as.character(form)[3]
if(grepl("\\|", x) | grepl("\\/",x)){
  x3<-trimws(str_split(x, "[|]|[/]")[[1]])
  x<-x3[1]
  individual = x3[2]
  group = x3[3]
}

sigmaForm<-if(match.arg(sigma, choices=sigmas)=="homo"){
  as.formula(paste0("sigma ~ 0+", group))
} else if (match.arg(sigma, choices=sigmas)=="linear"){
  as.formula(paste0("sigma ~ ", x, "+", x, ":",group))
} else if(match.arg(sigma, choices=sigmas)=="spline"){
  as.formula(paste0("sigma ~ s(",x,", by=", group, ")"))
}
if(match.arg(model, models)=="logistic"){
  form_fun<-form_logistic
} else if (match.arg(model, models)=="gompertz"){
  form_fun<-form_gompertz
} else if (match.arg(model, models)=="monomolecular"){
  form_fun<-form_monomolecular
} else if (match.arg(model, models)=="exponential"){
  form_fun<-form_exponential
} else if (match.arg(model, models)=="linear"){
  form_fun<-form_linear
} else if (match.arg(model, models)=="power law"){
  form_fun<-form_powerlaw
}
growthForm = form_fun(x,y)





form_logistic<-function(x, y){
  return(as.formula(paste0(y," ~ A/(1+exp((B-",x,")/C))")))
}
form_gompertz<-function(x, y){
  return(as.formula(paste0(y,"A*exp(-B*exp(-C*",x,"))")))
}
form_monomolecular<-function(x, y){
  return(as.formula(paste0(y,"~A-A*exp(-B*",x,")")))
}
form_exponential<-function(x, y){
  return(as.formula(paste0(y," ~ A*exp(B*",x, ")")))
}
form_linear<-function(x, y){
  return(as.formula(paste0(y," ~ A*",x)))
}
form_powerlaw<-function(x, y){
  return(as.formula(paste0(y," ~ A*",x,"^B")))
}










