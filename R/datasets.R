#' Gompertz growth model fit with brms to a small bellwether dataset.
#'
#' A \code{brmsfit} object for use with examples in vignettes, 
#'    documentation, and familiarization with models fit using \code{fitGrowth}.
#'    The model uses the gompertz growth formula and has 3 main parameters with
#'    several parameters to account for heteroskedasticity. The parameters starting
#'    with (b_) A, B, or C are from the gompertz formula and correspond to the
#'    asymptote, inflection point, and growth rate respectively.
#'    The model uses groups corresponding to fertilizer level (0, 50, 100) and
#'    genotype (MM, Mo17, B73, W605S). See brms documentation for full details 
#'    on \code{brmsfit} model objects. 
#'    \code{summary(fit)} will show a detailed summary of the model and some 
#'    fit diagnostics. The smooth term and population level effects are 
#'    explained briefly below.
#'
#' \itemize{
#'   \item Estimate: The mean of draws for this parameter. Use \code{as.data.frame(fit)} or
#'    \code{brms::as_draws(fit)} to access all draws.
#'   \item Est.Error: The error of the estimate.
#'   \item l-95\% CI: The lower bound on a 95 percent credible interval for this parameter.
#'   \item u-95\%: The upper bound on a 95 percent credible interval for this parameter.
#'   \item Rhat: A measure of convergence, this should be very close to 1 if a model fit well.
#'   Specifically, Rhat measures the between and within chain estimates. If your chains
#'   failed to mix then Rhat may be above 1.05, which is considered a viable cutoff for
#'   calling something a poor fit.
#'   \item Bulk_ESS: The bulk estimated sample size (ESS). That is, the effective number of samples
#'   that are being used for mean/median estimates. Broadly, ESS is a number of 
#'   post-warmup draws adjusted for the correlation between subsequent draws due to 
#'   the MCMC algorithm. These will depend on the length of your chains, the number
#'   of chains, and the quality of your model fit. There are several ideas about 
#'   how large of an ESS is large enough with low cutoffs being around 100. Broadly,
#'   larger is better.
#'   \item Tail_ESS: The ESS from tails of the distribution. Estimating tails of a
#'   posterior distribution is potentially difficult depending on the posterior, so
#'   this will be a lower number than Bulk_ESS but the same guidelines are used in 
#'   deciding if it is "large enough".
#' }
#'
#' @docType data
#' @keywords datasets
#' @name bw_vignette_fit
#' @usage data(bw_vignette_fit)
#' @format A brmsfit with 4000 draws across 257 total parameters
NULL


#' pcv.emd output from a subset of bellwether data.
#'
#' A list containing data and plot elements.
#' 
#' The "data" component is a data.frame with 363,609 rows and 9 columns.
#'
#' \itemize{
#'   \item i: The i-th row being compared
#'   \item j: The j-th row being compared
#'   \item emd: The Earth Mover's Distance between hue histograms between rows i and j.
#'   \item fertilizer_i: Fertilizer level for row i
#'   (these metadata terms were controlled by the \code{reorder} argument in \code{pcv.emd}).
#'   \item genotype_i: Genotype for row i
#'   \item  DAS_i: DAS for row i
#'   \item fertilizer_j : Fertilizer level for row j
#'   \item genotype_j : Genotype for row j
#'   \item DAS_j: DAS for row j
#' }
#'
#' @docType data
#' @keywords datasets
#' @name bw_vignette_small_emd
#' @usage data(bw_vignette_small_emd)
#' @format A list containing data and plot elements.
NULL