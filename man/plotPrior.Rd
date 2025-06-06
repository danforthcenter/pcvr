% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPrior.R
\name{plotPrior}
\alias{plotPrior}
\title{Check priors used in ease of use brms functions}
\usage{
plotPrior(priors, type = "density", n = 200, t = 25)
}
\arguments{
\item{priors}{A named list of means for prior distributions.
This takes the same input as the prior argument of \code{\link{growthSS}}.
Alternatively, if given the output of growthSS this will preform a prior predictive check
and return a plot from \code{\link{growthPlot}} of that check ignoring all other arguments.
Note that all priors must be
proper in that case (non-flat) and the fit is likely to be strange looking due to how thick
tailed the default priors from \code{\link{growthSS}} are.}

\item{type}{Either "density", the default, or a model as would be specified in \code{growthSS}
or \code{growthSim} such as "logistic", "gompertz", "monomolecular", "exponential",
"linear", "power law", "double logistic", or "double gompertz".
If this is a model type then n draws from the prior will be simulated as growth
trendlines and densities will be plotted on margins for some distributions.}

\item{n}{Numeric, if type is a model then how many draws from the prior should be simulated?}

\item{t}{Numeric, time passed to growthSim. Defaults to 25 (the growthSim default).}
}
\value{
A named list of plots showing prior distributions that \code{growthSS} would use,
optionally with a plot of simulated growth curves using draws from those priors.
}
\description{
Check priors used in ease of use brms functions
}
\examples{

set.seed(123)
priors <- list("A" = c(100, 130), "B" = c(10, 8), "C" = c(0.2, 0.1))
plotPrior(priors)

plotPrior(priors, "gompertz")[[1]]

}
\seealso{
\link{barg} for Bayesian model reporting metrics, \link{growthSim} for simulating data using
similar specification.
}
\keyword{Bayesian}
\keyword{brms}
\keyword{priors}
